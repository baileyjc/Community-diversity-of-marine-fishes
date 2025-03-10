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

    ##   |                  |          |   0%  |                  |          |   3%                                                                          |                  |.         |   7% [Bringing everything together load modifying files packages]             |                  |.         |  10%                                                                          |                  |.         |  13% [Bringing everything together load in modifying files]                   |                  |..        |  17%                                                                          |                  |..        |  20% [Check species names across files]                                       |                  |..        |  23%                                                                          |                  |...       |  27% [Modify environment data]                                                |                  |...       |  30%                                                                          |                  |...       |  33% [Modify incidence matrices]                                              |                  |....      |  37%                                                                          |                  |....      |  40% [Modify phylogeny]                                                       |                  |....      |  43%                                                                          |                  |.....     |  47% [Modify trait data]                                                      |                  |.....     |  50%                                                                          |                  |.....     |  53% [Modify Site_type data frames]                                           |                  |......    |  57%                                                                          |                  |......    |  60% [Modify Site_type trait data]                                            |                  |......    |  63%                                                                          |                  |.......   |  67% [Site_type trait data tests]                                             |                  |.......   |  70%                                                                          |                  |.......   |  73% [Modify site trait data frames]                                          |                  |........  |  77%                                                                          |                  |........  |  80% [Site trait data tests]                                                  |                  |........  |  83%                                                                          |                  |......... |  87% [unnamed-chunk-3]                                                        |                  |......... |  90%                                                                          |                  |......... |  93% [unnamed-chunk-4]                                                        |                  |..........|  97%                                                                          |                  |..........| 100% [Bringing everything together load out modified files and session info]

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
```

### Edit input files

``` r
strat_fish <- stratified_fish[,-c(1,10)]
mix_fish <- mixed_fish[,-c(1,10)]
oc_fish <- ocean_fish[,-c(1,8)]

env <- env[order(env$X),]
presabs_lake <- presabs_lake[order(row.names(presabs_lake)),]
surveyed_sites_lake <- surveyed_sites_lake[order(row.names(surveyed_sites_lake)),]

Site_type_group_ref <- env[,19]
Site_type_group <- env[surveyed_sites,19]
```

# Species Diversity

## SD alpha diversity

### SD alpha bar plot

``` r
#Plot species richness with a bar plot
SR_env$Site_type <- factor(SR_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

(SR_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = reorder(X, row_sum, decreasing = T), y = row_sum, color = Site_type, fill = Site_type, )) + 
  geom_bar(stat = 'identity',
    alpha = 1) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 14), 
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Site", y="Species Richness", colour = "Site type:", fill = "Site type:"))
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

### SD alpha outliers for each site type

``` r
outlier_SD_alpha <- SR_env[surveyed_sites,] %>%
  group_by(Site_type) %>%
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
(outlier_SD_alpha_plot <- ggplot(outlier_SD_alpha, aes(x = Site_type, y = row_sum, fill = Site_type)) +
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
  labs(y = "Species Richness", x = "Site type", color = "Outlier:", tag = "a"))
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20outliers-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/outlier_SD_alpha.jpg", outlier_SD_alpha_plot, width = 6.26, height = 6, units = "in")
```

### Identify VIF of environmental and geographical variables

``` r
# Determine which env variables have variance
environment <- c("temperature_median", "salinity_median", "oxygen_median", "pH_median")
nzv <- nearZeroVar(env[,environment])
nzv
```

    ## integer(0)

``` r
# All env variables have variance

mod <-  lm(row_sum ~ salinity_median + oxygen_median + temperature_median + pH_median, SR_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ salinity_median + oxygen_median + temperature_median + 
    ##     pH_median, data = SR_env)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -25.9542 -10.3708  -0.4769   8.9325  29.5684 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -541.813    201.886  -2.684   0.0178 *
    ## salinity_median       0.201      1.973   0.102   0.9203  
    ## oxygen_median         6.058      5.343   1.134   0.2758  
    ## temperature_median   -4.950      3.874  -1.278   0.2221  
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
    ## salinity_median     1 9172.6  9172.6 34.6016 3.989e-05 ***
    ## oxygen_median       1 3837.2  3837.2 14.4750  0.001933 ** 
    ## temperature_median  1   13.2    13.2  0.0499  0.826414    
    ## pH_median           1 1789.6  1789.6  6.7509  0.021049 *  
    ## Residuals          14 3711.3   265.1                      
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
mod <-  lm(row_sum ~ salinity_median + oxygen_median + temperature_median, SR_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = SR_env)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -31.549  -9.598  -4.173   9.619  47.290 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)        -131.1738   147.7497  -0.888  0.38866   
    ## salinity_median       4.1712     1.4680   2.841  0.01238 * 
    ## oxygen_median        15.2171     4.7218   3.223  0.00569 **
    ## temperature_median   -0.7881     4.1484  -0.190  0.85187   
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

### SD alpha & site type

``` r
### Site_type
# Run ANOVA
anova_logSRic_result <- aov(log(row_sum) ~ Site_type, data = SR_env[surveyed_sites,])
summary(anova_logSRic_result)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Site_type    2 22.346  11.173   27.41 2.51e-06 ***
    ## Residuals   19  7.745   0.408                     
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
    ## Fit: aov(formula = log(row_sum) ~ Site_type, data = SR_env[surveyed_sites, ])
    ## 
    ## $Site_type
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
Site_type <- tukey_logSRic_result$Site_type

# q-values
(OM_q <- (abs(Site_type[1,1]))/USE)
```

    ## [1] 0.1259392

``` r
(MS_q <- (abs(Site_type[3,1]))/BSE)
```

    ## [1] 9.222747

``` r
(SO_q <- (abs(Site_type[2,1]))/USE)
```

    ## [1] 8.664544

``` r
# Run ANOVA
anova_SRic_result <- aov(row_sum ~ Site_type, data = SR_env[surveyed_sites,])
summary(anova_SRic_result)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Site_type    2  10743    5371   11.01 0.000667 ***
    ## Residuals   19   9268     488                     
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
    ## Fit: aov(formula = row_sum ~ Site_type, data = SR_env[surveyed_sites, ])
    ## 
    ## $Site_type
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
Site_type <- tukey_SRic_result$Site_type

# q-values
(OM_q <- (abs(Site_type[1,1]))/USE)
```

    ## [1] 0.1235068

``` r
(MS_q <- (abs(Site_type[3,1]))/BSE)
```

    ## [1] 5.939085

``` r
(SO_q <- (abs(Site_type[2,1]))/USE)
```

    ## [1] 5.375017

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx
```

### SD alpha with env and geo linear models & ANOVAs

``` r
### Environmental
# Surveyed sites
SD_logSRic_env_am <- aov(log(row_sum) ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = SR_env[surveyed_sites_env,])
summary(SD_logSRic_env_am)
```

    ##                              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median               1 22.260  22.260  67.896 3.53e-05 ***
    ## oxygen_median                 1  2.548   2.548   7.771   0.0236 *  
    ## temperature_median            1  0.023   0.023   0.069   0.7992    
    ## Site_type                     2  0.144   0.072   0.219   0.8079    
    ## salinity_median:Site_type     2  0.437   0.219   0.667   0.5397    
    ## oxygen_median:Site_type       2  0.416   0.208   0.634   0.5553    
    ## temperature_median:Site_type  1  0.041   0.041   0.125   0.7332    
    ## Residuals                     8  2.623   0.328                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_env_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.0002469618 0.1655141764 1.0000000000 1.0000000000 1.0000000000
    ## [6] 1.0000000000 1.0000000000           NA

``` r
SD_logSRic_env_lm <- lm(log(row_sum) ~ salinity_median + oxygen_median + temperature_median, SR_env[surveyed_sites_env,])
summary(SD_logSRic_env_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = SR_env[surveyed_sites_env, ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.73012 -0.20864 -0.06934  0.27642  0.88318 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -7.63922    3.81124  -2.004  0.06343 .  
    ## salinity_median     0.26035    0.03787   6.876 5.27e-06 ***
    ## oxygen_median       0.39484    0.12180   3.242  0.00548 ** 
    ## temperature_median  0.03262    0.10701   0.305  0.76470    
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

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##       2.537061e-01       2.109736e-05       2.190101e-02       1.000000e+00

``` r
SD_SRic_env_am <- aov(row_sum ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = SR_env[surveyed_sites_env,])
summary(SD_SRic_env_am)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)   
    ## salinity_median               1   9173    9173  23.348 0.0013 **
    ## oxygen_median                 1   3837    3837   9.767 0.0141 * 
    ## temperature_median            1     13      13   0.034 0.8589   
    ## Site_type                     2     71      36   0.091 0.9144   
    ## salinity_median:Site_type     2    690     345   0.878 0.4521   
    ## oxygen_median:Site_type       2   1151     575   1.465 0.2870   
    ## temperature_median:Site_type  1    446     446   1.136 0.3177   
    ## Residuals                     8   3143     393                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_env_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.009108864 0.098807167 1.000000000 1.000000000 1.000000000 1.000000000
    ## [7] 1.000000000          NA

``` r
SD_SRic_env_lm <- lm(row_sum ~ salinity_median + oxygen_median + temperature_median, SR_env[surveyed_sites_env,])
summary(SD_SRic_env_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = SR_env[surveyed_sites_env, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -31.549  -9.598  -4.173   9.619  47.290 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)        -131.1738   147.7497  -0.888  0.38866   
    ## salinity_median       4.1712     1.4680   2.841  0.01238 * 
    ## oxygen_median        15.2171     4.7218   3.223  0.00569 **
    ## temperature_median   -0.7881     4.1484  -0.190  0.85187   
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

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##         1.00000000         0.04950985         0.02277051         1.00000000

``` r
# Mixed and stratified lakes
SD_logSRic_env_MS_am <- aov(log(row_sum) ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_env_MS_am)
```

    ##                              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median               1 18.016  18.016  54.951 7.53e-05 ***
    ## oxygen_median                 1  2.168   2.168   6.613    0.033 *  
    ## temperature_median            1  0.019   0.019   0.057    0.818    
    ## Site_type                     1  0.119   0.119   0.364    0.563    
    ## salinity_median:Site_type     1  0.461   0.461   1.407    0.270    
    ## oxygen_median:Site_type       1  0.338   0.338   1.031    0.340    
    ## temperature_median:Site_type  1  0.041   0.041   0.125    0.733    
    ## Residuals                     8  2.623   0.328                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_env_MS_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.0005268214 0.2313246604 1.0000000000 1.0000000000 1.0000000000
    ## [6] 1.0000000000 1.0000000000           NA

``` r
SD_logSRic_env_MS_lm <- lm(log(row_sum) ~ salinity_median + oxygen_median + temperature_median, SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_env_MS_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = SR_env[mixed_stratified_lakes, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.71621 -0.27163 -0.06681  0.37624  0.91222 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -7.73929    4.51110  -1.716   0.1119    
    ## salinity_median     0.26303    0.04265   6.168 4.81e-05 ***
    ## oxygen_median       0.41650    0.15412   2.702   0.0192 *  
    ## temperature_median  0.03109    0.12473   0.249   0.8074    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5464 on 12 degrees of freedom
    ## Multiple R-squared:  0.8494, Adjusted R-squared:  0.8117 
    ## F-statistic: 22.56 on 3 and 12 DF,  p-value: 3.194e-05

``` r
p_values <- summary(SD_logSRic_env_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##       0.4476551535       0.0001925785       0.0768885411       1.0000000000

``` r
SD_SRic_env_MS_am <- aov(row_sum ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_SRic_env_MS_am)
```

    ##                              Df Sum Sq Mean Sq F value  Pr(>F)   
    ## salinity_median               1   6790    6790  17.283 0.00318 **
    ## oxygen_median                 1   3039    3039   7.735 0.02388 * 
    ## temperature_median            1     72      72   0.182 0.68086   
    ## Site_type                     1      1       1   0.002 0.96894   
    ## salinity_median:Site_type     1    772     772   1.966 0.19845   
    ## oxygen_median:Site_type       1    629     629   1.602 0.24125   
    ## temperature_median:Site_type  1    446     446   1.136 0.31766   
    ## Residuals                     8   3143     393                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_env_MS_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.02223932 0.16716232 1.00000000 1.00000000 1.00000000 1.00000000 1.00000000
    ## [8]         NA

``` r
SD_SRic_env_MS_lm <- lm(row_sum ~ salinity_median + oxygen_median + temperature_median, SR_env[mixed_stratified_lakes,])
summary(SD_SRic_env_MS_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = SR_env[mixed_stratified_lakes, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -31.897  -9.129  -2.596   7.932  45.355 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)         -94.826    168.390  -0.563   0.5837  
    ## salinity_median       4.176      1.592   2.623   0.0223 *
    ## oxygen_median        15.034      5.753   2.613   0.0227 *
    ## temperature_median   -1.931      4.656  -0.415   0.6857  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 20.39 on 12 degrees of freedom
    ## Multiple R-squared:  0.6648, Adjusted R-squared:  0.581 
    ## F-statistic: 7.934 on 3 and 12 DF,  p-value: 0.003509

``` r
p_values <- summary(SD_SRic_env_MS_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##         1.00000000         0.08902474         0.09067730         1.00000000

``` r
# Ocean sites and mixed lakes
SD_logSRic_env_OM_am <- aov(log(row_sum) ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_env_OM_am)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)
    ## salinity_median            1 1.1092  1.1092   3.328  0.142
    ## oxygen_median              1 0.2667  0.2667   0.800  0.422
    ## temperature_median         1 0.0326  0.0326   0.098  0.770
    ## Site_type                  1 0.2353  0.2353   0.706  0.448
    ## salinity_median:Site_type  1 0.0058  0.0058   0.017  0.901
    ## oxygen_median:Site_type    1 0.0398  0.0398   0.119  0.747
    ## Residuals                  4 1.3333  0.3333

``` r
p_values <- summary(SD_logSRic_env_OM_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.853145 1.000000 1.000000 1.000000 1.000000 1.000000       NA

``` r
SD_logSRic_env_OM_lm <- lm(log(row_sum) ~ salinity_median + oxygen_median + temperature_median, SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_env_OM_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = SR_env[ocean_mixed_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.52504 -0.28816 -0.07468  0.20202  0.64145 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -8.29859   18.14153  -0.457    0.661
    ## salinity_median     0.39271    0.46891   0.838    0.430
    ## oxygen_median       0.43003    0.37824   1.137    0.293
    ## temperature_median -0.09415    0.25046  -0.376    0.718
    ## 
    ## Residual standard error: 0.4802 on 7 degrees of freedom
    ## Multiple R-squared:  0.466,  Adjusted R-squared:  0.2371 
    ## F-statistic: 2.036 on 3 and 7 DF,  p-value: 0.1975

``` r
p_values <- summary(SD_logSRic_env_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
SD_SRic_env_OM_am <- aov(row_sum ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_SRic_env_OM_am)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)
    ## salinity_median            1 1693.9  1693.9   2.216  0.211
    ## oxygen_median              1  839.7   839.7   1.098  0.354
    ## temperature_median         1  284.5   284.5   0.372  0.575
    ## Site_type                  1  569.3   569.3   0.745  0.437
    ## salinity_median:Site_type  1    3.6     3.6   0.005  0.948
    ## oxygen_median:Site_type    1  574.8   574.8   0.752  0.435
    ## Residuals                  4 3058.3   764.6

``` r
p_values <- summary(SD_SRic_env_OM_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1]  1  1  1  1  1  1 NA

``` r
SD_SRic_env_OM_lm <- lm(row_sum ~ salinity_median + oxygen_median + temperature_median, SR_env[ocean_mixed_sites_env,])
summary(SD_SRic_env_OM_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = SR_env[ocean_mixed_sites_env, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -24.059 -18.215  -3.151  11.182  38.375 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         -37.380    926.047  -0.040    0.969
    ## salinity_median       7.098     23.936   0.297    0.775
    ## oxygen_median        26.197     19.308   1.357    0.217
    ## temperature_median   -8.798     12.785  -0.688    0.514
    ## 
    ## Residual standard error: 24.51 on 7 degrees of freedom
    ## Multiple R-squared:  0.4012, Adjusted R-squared:  0.1446 
    ## F-statistic: 1.563 on 3 and 7 DF,  p-value: 0.2814

``` r
p_values <- summary(SD_SRic_env_OM_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          1.0000000          0.8678638          1.0000000

``` r
# Stratified lakes and ocean sites
SD_logSRic_env_SO_am <- aov(log(row_sum) ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = SR_env[ocean_stratified_sites_env,])
summary(SD_logSRic_env_SO_am)
```

    ##                           Df Sum Sq Mean Sq F value  Pr(>F)   
    ## salinity_median            1 12.655  12.655  39.255 0.00331 **
    ## oxygen_median              1  1.731   1.731   5.370 0.08138 . 
    ## temperature_median         1  0.104   0.104   0.324 0.59987   
    ## Site_type                  1  0.398   0.398   1.235 0.32870   
    ## salinity_median:Site_type  1  0.000   0.000   0.002 0.97088   
    ## oxygen_median:Site_type    1  0.149   0.149   0.463 0.53352   
    ## Residuals                  4  1.290   0.322                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_env_SO_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.01986745 0.48826519 1.00000000 1.00000000 1.00000000 1.00000000         NA

``` r
SD_logSRic_env_SO_lm <- lm(log(row_sum) ~ salinity_median + oxygen_median + temperature_median, SR_env[ocean_stratified_sites_env,])
summary(SD_logSRic_env_SO_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = SR_env[ocean_stratified_sites_env, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.74751 -0.14210 -0.02821  0.14682  0.84158 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -8.58582    4.18049  -2.054 0.079087 .  
    ## salinity_median     0.25128    0.04292   5.854 0.000628 ***
    ## oxygen_median       0.34943    0.13583   2.573 0.036866 *  
    ## temperature_median  0.07539    0.11959   0.630 0.548449    
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

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##        0.316349339        0.002511461        0.147465030        1.000000000

``` r
SD_SRic_env_SO_am <- aov(row_sum ~ (salinity_median + oxygen_median + temperature_median) * Site_type, data = SR_env[ocean_stratified_sites_env,])
summary(SD_SRic_env_SO_am)
```

    ##                           Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median            1   4627    4627 218.879 0.000122 ***
    ## oxygen_median              1   2476    2476 117.119 0.000414 ***
    ## temperature_median         1     60      60   2.848 0.166784    
    ## Site_type                  1    171     171   8.088 0.046684 *  
    ## salinity_median:Site_type  1      3       3   0.155 0.714093    
    ## oxygen_median:Site_type    1    674     674  31.880 0.004845 ** 
    ## Residuals                  4     85      21                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_env_SO_am)[[1]][, "Pr(>F)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.0007290897 0.0024815321 1.0000000000 0.2801063330 1.0000000000
    ## [6] 0.0290727615           NA

``` r
SD_SRic_env_SO_lm <- lm(row_sum ~ salinity_median + oxygen_median + temperature_median, SR_env[ocean_stratified_sites_env,])
summary(SD_SRic_env_SO_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = SR_env[ocean_stratified_sites_env, ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -12.023  -6.847  -2.633   3.538  17.340 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)        -206.9563    94.1879  -2.197  0.06399 . 
    ## salinity_median       4.2382     0.9671   4.383  0.00322 **
    ## oxygen_median        13.2060     3.0602   4.315  0.00350 **
    ## temperature_median    1.8111     2.6945   0.672  0.52306   
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

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##         0.25595351         0.01289859         0.01399934         1.00000000

``` r
# Ocean sites
SD_logSRic_env_O_lm <- lm(log(row_sum) ~ salinity_median + oxygen_median + temperature_median, SR_env[ocean_sites,])
summary(SD_logSRic_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = SR_env[ocean_sites, ])
    ## 
    ## Residuals:
    ## ALL 3 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)         -2.2547        NaN     NaN      NaN
    ## salinity_median      0.0550        NaN     NaN      NaN
    ## oxygen_median        0.8351        NaN     NaN      NaN
    ## temperature_median       NA         NA      NA       NA
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
    ##                 Df   Sum Sq  Mean Sq F value Pr(>F)
    ## salinity_median  1 0.001745 0.001745     NaN    NaN
    ## oxygen_median    1 0.162841 0.162841     NaN    NaN
    ## Residuals        0 0.000000      NaN

``` r
p_values <- summary(SD_logSRic_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##     (Intercept) salinity_median   oxygen_median 
    ##             NaN             NaN             NaN

``` r
SD_SRic_env_O_lm <- lm(row_sum ~ salinity_median + oxygen_median + temperature_median, SR_env[ocean_sites,])
summary(SD_SRic_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = SR_env[ocean_sites, ])
    ## 
    ## Residuals:
    ## ALL 3 residuals are 0: no residual degrees of freedom!
    ## 
    ## Coefficients: (1 not defined because of singularities)
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -447.624        NaN     NaN      NaN
    ## salinity_median       5.997        NaN     NaN      NaN
    ## oxygen_median        57.123        NaN     NaN      NaN
    ## temperature_median       NA         NA      NA       NA
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
    ##                 Df Sum Sq Mean Sq F value Pr(>F)
    ## salinity_median  1  10.73   10.73     NaN    NaN
    ## oxygen_median    1 761.94  761.94     NaN    NaN
    ## Residuals        0   0.00     NaN

``` r
p_values <- summary(SD_SRic_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##     (Intercept) salinity_median   oxygen_median 
    ##             NaN             NaN             NaN

``` r
# Mixed lakes
SD_logSRic_env_M_lm <- lm(log(row_sum) ~ salinity_median + oxygen_median + temperature_median, SR_env[mixed_lakes,])
summary(SD_logSRic_env_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = SR_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##     FLK     HLO     LLN     MLN     NLN     NLU     OLO     ULN 
    ##  0.1488 -0.2163  0.5206 -0.3148  0.6425 -0.4922 -0.4585  0.1700 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -14.6394    27.2353  -0.538    0.619
    ## salinity_median      0.5768     0.6332   0.911    0.414
    ## oxygen_median        0.5189     0.5122   1.013    0.368
    ## temperature_median  -0.0965     0.3930  -0.246    0.818
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
    ## salinity_median     1 0.97758 0.97758  2.9327 0.1620
    ## oxygen_median       1 0.33089 0.33089  0.9927 0.3755
    ## temperature_median  1 0.02010 0.02010  0.0603 0.8181
    ## Residuals           4 1.33334 0.33334

``` r
p_values <- summary(SD_logSRic_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
SD_SRic_env_M_lm <- lm(row_sum ~ salinity_median + oxygen_median + temperature_median, SR_env[mixed_lakes,])
summary(SD_SRic_env_M_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = SR_env[mixed_lakes, ])
    ## 
    ## Residuals:
    ##     FLK     HLO     LLN     MLN     NLN     NLU     OLO     ULN 
    ##   7.440 -13.112  28.023 -11.289  30.672 -23.921 -20.011   2.198 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)          -32.78    1304.37  -0.025    0.981
    ## salinity_median       12.84      30.32   0.424    0.694
    ## oxygen_median         27.71      24.53   1.130    0.322
    ## temperature_median   -15.30      18.82  -0.813    0.462
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
    ## salinity_median     1 1580.33 1580.33  2.0669 0.2239
    ## oxygen_median       1  868.84  868.84  1.1364 0.3465
    ## temperature_median  1  505.41  505.41  0.6610 0.4618
    ## Residuals           4 3058.30  764.57

``` r
p_values <- summary(SD_SRic_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
# Stratified lakes
SD_logSRic_env_S_lm <- lm(log(row_sum) ~ salinity_median + oxygen_median + temperature_median, SR_env[stratified_lakes,])
summary(SD_logSRic_env_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = SR_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##      BCM      CLM      GLK      HLM      NLK      OTM      SLN      TLN 
    ##  0.37953  0.04308  0.17921  0.24651 -0.28497 -0.87230 -0.12917  0.43811 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -3.19314    6.27810  -0.509    0.638
    ## salinity_median     0.14298    0.09750   1.466    0.216
    ## oxygen_median      -0.16727    0.42987  -0.389    0.717
    ## temperature_median  0.04973    0.13905   0.358    0.739
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
    ##                    Df  Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median     1 2.39771 2.39771  7.4375 0.0526 .
    ## oxygen_median       1 0.05939 0.05939  0.1842 0.6899  
    ## temperature_median  1 0.04124 0.04124  0.1279 0.7387  
    ## Residuals           4 1.28953 0.32238                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_logSRic_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          0.8657002          1.0000000          1.0000000

``` r
SD_SRic_env_S_lm <- lm(row_sum ~ salinity_median + oxygen_median + temperature_median, SR_env[stratified_lakes,])
summary(SD_SRic_env_S_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = SR_env[stratified_lakes, ])
    ## 
    ## Residuals:
    ##     BCM     CLM     GLK     HLM     NLK     OTM     SLN     TLN 
    ##  2.3022  1.9724  0.9048  2.9002 -3.2143 -6.6194 -1.4092  3.1633 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -12.67663   50.83889  -0.249    0.815
    ## salinity_median      1.04089    0.78954   1.318    0.258
    ## oxygen_median       -2.41766    3.48100  -0.695    0.526
    ## temperature_median  -0.02215    1.12601  -0.020    0.985
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
    ## salinity_median     1 181.181 181.181  8.5705 0.04292 *
    ## oxygen_median       1  10.250  10.250  0.4849 0.52456  
    ## temperature_median  1   0.008   0.008  0.0004 0.98525  
    ## Residuals           4  84.560  21.140                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
### Geographical
# Surveyed sites
SD_logSRic_geo_am <- aov(log(row_sum) ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = SR_env[surveyed_sites,])
summary(SD_logSRic_geo_am)
```

    ##                                   Df Sum Sq Mean Sq F value  Pr(>F)    
    ## distance_to_ocean_min_m            1 17.419  17.419  46.808 4.5e-05 ***
    ## max_depth                          1  0.021   0.021   0.055 0.81908    
    ## logArea                            1  0.136   0.136   0.366 0.55882    
    ## Site_type                          2  7.445   3.723  10.003 0.00411 ** 
    ## distance_to_ocean_min_m:Site_type  2  0.218   0.109   0.292 0.75273    
    ## max_depth:Site_type                2  0.505   0.253   0.679 0.52907    
    ## logArea:Site_type                  2  0.625   0.313   0.840 0.45998    
    ## Residuals                         10  3.721   0.372                    
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
SD_SRic_geo_am <- aov(row_sum ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = SR_env[surveyed_sites,])
summary(SD_SRic_geo_am)
```

    ##                                   Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m            1   6762    6762  20.745 0.00105 **
    ## max_depth                          1    428     428   1.314 0.27838   
    ## logArea                            1    699     699   2.146 0.17367   
    ## Site_type                          2   5631    2816   8.638 0.00662 **
    ## distance_to_ocean_min_m:Site_type  2    149      74   0.229 0.79972   
    ## max_depth:Site_type                2   1554     777   2.383 0.14244   
    ## logArea:Site_type                  2   1527     764   2.343 0.14638   
    ## Residuals                         10   3260     326                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_geo_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
SD_SRic_geo_lm <- lm(row_sum ~  max_depth + logArea * Site_type, SR_env[surveyed_sites,])
summary(SD_SRic_geo_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ max_depth + logArea * Site_type, data = SR_env[surveyed_sites, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -27.6617  -7.6425   0.1206   9.7574  23.2474 
    ## 
    ## Coefficients:
    ##                              Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)                   24.5607    33.3088   0.737   0.4723  
    ## max_depth                      0.8677     0.4735   1.833   0.0868 .
    ## logArea                        1.3683     2.8789   0.475   0.6414  
    ## Site_typeMixed              -111.7171    52.8812  -2.113   0.0518 .
    ## Site_typeStratified           46.0307    69.7863   0.660   0.5195  
    ## logArea:Site_typeMixed        12.0791     5.0741   2.381   0.0310 *
    ## logArea:Site_typeStratified   -9.5152     6.8430  -1.391   0.1847  
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

    ##                 (Intercept)                   max_depth 
    ##                   1.0000000                   0.6075203 
    ##                     logArea              Site_typeMixed 
    ##                   1.0000000                   0.3626745 
    ##         Site_typeStratified      logArea:Site_typeMixed 
    ##                   1.0000000                   0.2168482 
    ## logArea:Site_typeStratified 
    ##                   1.0000000

``` r
# Mixed and stratified lakes
SD_logSRic_geo_MS_am <- aov(log(row_sum) ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_geo_MS_am)
```

    ##                                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m            1 12.814  12.814  34.386 0.000377 ***
    ## max_depth                          1  0.075   0.075   0.202 0.664772    
    ## logArea                            1  2.999   2.999   8.048 0.021915 *  
    ## Site_type                          1  4.624   4.624  12.410 0.007814 ** 
    ## distance_to_ocean_min_m:Site_type  1  0.000   0.000   0.000 0.991222    
    ## max_depth:Site_type                1  0.246   0.246   0.659 0.440430    
    ## logArea:Site_type                  1  0.046   0.046   0.123 0.734712    
    ## Residuals                          8  2.981   0.373                     
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
SD_SRic_geo_MS_am <- aov(row_sum ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_SRic_geo_MS_am)
```

    ##                                   Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m            1   4715    4715  21.003 0.00180 **
    ## max_depth                          1      1       1   0.005 0.94705   
    ## logArea                            1   4426    4426  19.718 0.00217 **
    ## Site_type                          1   2800    2800  12.474 0.00771 **
    ## distance_to_ocean_min_m:Site_type  1    110     110   0.491 0.50343   
    ## max_depth:Site_type                1    703     703   3.131 0.11476   
    ## logArea:Site_type                  1    340     340   1.514 0.25342   
    ## Residuals                          8   1796     224                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
p_values <- summary(SD_SRic_geo_MS_am)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## numeric(0)

``` r
SD_SRic_geo_MS_lm <- lm(row_sum ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, SR_env[mixed_stratified_lakes,])
summary(SD_SRic_geo_MS_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ (distance_to_ocean_min_m + max_depth + 
    ##     logArea) * Site_type, data = SR_env[mixed_stratified_lakes, 
    ##     ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -28.0861  -5.3118   0.2086   4.9038  22.0988 
    ## 
    ## Coefficients:
    ##                                              Estimate Std. Error t value
    ## (Intercept)                                 -99.47259   41.77624  -2.381
    ## distance_to_ocean_min_m                      -0.08445    0.24707  -0.342
    ## max_depth                                     0.09009    1.07945   0.083
    ## logArea                                      16.35700    5.54093   2.952
    ## Site_typeStratified                         108.13429  101.28289   1.068
    ## distance_to_ocean_min_m:Site_typeStratified   0.02938    0.25922   0.113
    ## max_depth:Site_typeStratified                -0.10032    1.55089  -0.065
    ## logArea:Site_typeStratified                 -15.61396   12.68790  -1.231
    ##                                             Pr(>|t|)  
    ## (Intercept)                                   0.0445 *
    ## distance_to_ocean_min_m                       0.7413  
    ## max_depth                                     0.9355  
    ## logArea                                       0.0184 *
    ## Site_typeStratified                           0.3168  
    ## distance_to_ocean_min_m:Site_typeStratified   0.9126  
    ## max_depth:Site_typeStratified                 0.9500  
    ## logArea:Site_typeStratified                   0.2534  
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

    ##                                 (Intercept) 
    ##                                   0.3557648 
    ##                     distance_to_ocean_min_m 
    ##                                   1.0000000 
    ##                                   max_depth 
    ##                                   1.0000000 
    ##                                     logArea 
    ##                                   0.1469491 
    ##                         Site_typeStratified 
    ##                                   1.0000000 
    ## distance_to_ocean_min_m:Site_typeStratified 
    ##                                   1.0000000 
    ##               max_depth:Site_typeStratified 
    ##                                   1.0000000 
    ##                 logArea:Site_typeStratified 
    ##                                   1.0000000

``` r
# Ocean sites and mixed lakes
SD_logSRic_geo_OM_am <- aov(log(row_sum) ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = SR_env[ocean_mixed_sites,])
summary(SD_logSRic_geo_OM_am)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m            1 0.0712  0.0712   0.309 0.5981  
    ## max_depth                          1 1.6259  1.6259   7.063 0.0376 *
    ## logArea                            1 0.1182  0.1182   0.514 0.5005  
    ## Site_type                          1 0.0121  0.0121   0.052 0.8266  
    ## distance_to_ocean_min_m:Site_type  1 0.0101  0.0101   0.044 0.8410  
    ## max_depth:Site_type                1 0.1230  0.1230   0.534 0.4923  
    ## logArea:Site_type                  1 0.6182  0.6182   2.686 0.1524  
    ## Residuals                          6 1.3812  0.2302                 
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
SD_SRic_geo_OM_am <- aov(row_sum ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = SR_env[ocean_mixed_sites,])
summary(SD_SRic_geo_OM_am)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m            1     14      14   0.027 0.8751  
    ## max_depth                          1   3536    3536   6.853 0.0397 *
    ## logArea                            1    702     702   1.361 0.2876  
    ## Site_type                          1     16      16   0.032 0.8641  
    ## distance_to_ocean_min_m:Site_type  1      2       2   0.005 0.9476  
    ## max_depth:Site_type                1    120     120   0.232 0.6468  
    ## logArea:Site_type                  1   1509    1509   2.925 0.1381  
    ## Residuals                          6   3096     516                 
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
SD_logSRic_geo_SO_am <- aov(log(row_sum) ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = SR_env[ocean_stratified_sites,])
summary(SD_logSRic_geo_SO_am)
```

    ##                                   Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m            1 14.817  14.817  28.858 0.00171 **
    ## max_depth                          1  0.265   0.265   0.515 0.49978   
    ## logArea                            1  0.161   0.161   0.313 0.59604   
    ## Site_type                          1  1.849   1.849   3.601 0.10653   
    ## distance_to_ocean_min_m:Site_type  1  0.009   0.009   0.017 0.90111   
    ## max_depth:Site_type                1  0.162   0.162   0.316 0.59413   
    ## logArea:Site_type                  1  0.041   0.041   0.079 0.78809   
    ## Residuals                          6  3.081   0.513                   
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
SD_SRic_geo_SO_am <- aov(row_sum ~ (distance_to_ocean_min_m + max_depth + logArea) * Site_type, data = SR_env[ocean_stratified_sites,])
summary(SD_SRic_geo_SO_am)
```

    ##                                   Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m            1   5520    5520  20.349 0.00406 **
    ## max_depth                          1    432     432   1.594 0.25355   
    ## logArea                            1    293     293   1.080 0.33883   
    ## Site_type                          1   1572    1572   5.794 0.05279 . 
    ## distance_to_ocean_min_m:Site_type  1      1       1   0.005 0.94809   
    ## max_depth:Site_type                1    853     853   3.146 0.12647   
    ## logArea:Site_type                  1      3       3   0.009 0.92647   
    ## Residuals                          6   1627     271                   
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

### SD alpha with env and geo correlated variables

``` r
# Log Area
SD_alpha_LA_plot <- ggplot(data = SR_env[surveyed_sites,], mapping = aes(y = log(row_sum), x = logArea, color = Site_type, fill = Site_type)) + 
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
  labs(x="Log Area (m"^"2"~")", y="Log SRic", colour = "Site type:", fill = "Site type:", tag = "c")
(SD_alpha_LA_plot <- SD_alpha_LA_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_alpha_LA_plot.jpg", plot = SD_alpha_LA_plot, width = 6.26, height = 6, units = "in")

# Distance from the ocean mean
SD_alpha_D_plot <- ggplot(data = SR_env[surveyed_sites,], mapping = aes(y = row_sum, x = distance_to_ocean_min_m, color = Site_type, fill = Site_type)) + 
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
  labs(x="Isolation (m)", y="SRic", colour = "Site type:", fill = "Site type:", tag = "b")
(SD_alpha_D_plot <- SD_alpha_D_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_alpha_D_plot.jpg", plot = SD_alpha_D_plot, width = 6.26, height = 6, units = "in")

# Max depth
SD_alpha_MD_plot <- ggplot(data = SR_env[surveyed_sites,], mapping = aes(y = row_sum, x = max_depth, color = Site_type, fill = Site_type)) + 
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
  labs(x="Age (m)", y="SRic", colour = "Site type:", fill = "Site type:", tag = "a")
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
(SD_beta_ref_BD <- betadisper(SD_beta_ref_dist$Btotal, Site_type_group_ref))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_ref_dist$Btotal, group =
    ## Site_type_group_ref)
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
    ## 3 Stratified vs Reference  1 0.6261283 1.909357 0.2143091   0.119      0.714
    ## 4          Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.042      0.252
    ## 5      Mixed vs Reference  1 0.5979203 1.966332 0.2193017   0.120      0.720
    ## 6      Ocean vs Reference  1 0.5599217 1.628213 0.2456489   0.150      0.900
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
(SD_beta_BD <- betadisper(SD_beta_dist$Btotal, Site_type_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_dist$Btotal, group = Site_type_group)
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
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.046      0.138

``` r
# Without LCN
SD_beta_wo_LCN_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_LCN,], abund = F)
(SD_beta_wo_LCN_BD <- betadisper(SD_beta_wo_LCN_dist$Btotal, Site_type_group[-8]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_wo_LCN_dist$Btotal, group =
    ## Site_type_group[-8])
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
    ## 1 Stratified vs Mixed  1 1.3139866 4.158150 0.2289964   0.002      0.006   *
    ## 2 Stratified vs Ocean  1 1.3841204 4.397984 0.2856207   0.002      0.006   *
    ## 3      Mixed vs Ocean  1 0.6396174 2.135322 0.1625633   0.007      0.021   .

``` r
# Without TLN and HLM
SD_beta_wo_TLN_HLM_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_TLN_HLM,], abund = F)
(SD_beta_wo_TLN_HLM_BD <- betadisper(SD_beta_wo_TLN_HLM_dist$Btotal, Site_type_group[-c(5,21)]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_wo_TLN_HLM_dist$Btotal, group =
    ## Site_type_group[-c(5, 21)])
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
    ## 2 Stratified vs Ocean  1 1.3260066 4.147587 0.2931657   0.009      0.027   .
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.057      0.171

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

### SD beta dispersions

- The dissimilarity between sites of the same site type

``` r
# Dispersion within-groups
SD_beta_BD_dist <- SD_beta_BD$distances
SD_beta_BD_dist <- as.data.frame(SD_beta_BD_dist)
SD_beta_BD_dist$X <- row.names(SD_beta_BD_dist)
SD_beta_BD_dist_env <- merge(SD_beta_BD_dist, env[surveyed_sites,], by = "X", sort = F)

SD_beta_BD_dist_env$Site_type <- factor(SD_beta_BD_dist_env$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

(SD_beta_BD_dist_env_plot <- ggplot(SD_beta_BD_dist_env, aes(x = Site_type, y = SD_beta_BD_dist, color = Site_type, fill = Site_type)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Site_type, fill = Site_type)) +
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
  labs(color = "Site_type", tag = "a"))
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
    ## Run 1 stress 0.1082722 
    ## Run 2 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307073  max resid 0.07844779 
    ## Run 3 stress 0.1138394 
    ## Run 4 stress 0.1063128 
    ## Run 5 stress 0.1051402 
    ## ... Procrustes: rmse 0.02306853  max resid 0.0784679 
    ## Run 6 stress 0.113594 
    ## Run 7 stress 0.1063128 
    ## Run 8 stress 0.1051402 
    ## ... Procrustes: rmse 0.02306826  max resid 0.07846811 
    ## Run 9 stress 0.1063128 
    ## Run 10 stress 0.1135599 
    ## Run 11 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307294  max resid 0.07848206 
    ## Run 12 stress 0.1088547 
    ## Run 13 stress 0.113594 
    ## Run 14 stress 0.1088548 
    ## Run 15 stress 0.113594 
    ## Run 16 stress 0.1097989 
    ## Run 17 stress 0.1135597 
    ## Run 18 stress 0.1084964 
    ## Run 19 stress 0.1063128 
    ## Run 20 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214365  max resid 0.0789568 
    ## Run 21 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214785  max resid 0.07898547 
    ## Run 22 stress 0.1051402 
    ## ... Procrustes: rmse 0.02306877  max resid 0.07846957 
    ## Run 23 stress 0.1083502 
    ## Run 24 stress 0.1050338 
    ## ... Procrustes: rmse 0.008873992  max resid 0.0248306 
    ## Run 25 stress 0.1063129 
    ## Run 26 stress 0.1063128 
    ## Run 27 stress 0.1135597 
    ## Run 28 stress 0.1051402 
    ## ... Procrustes: rmse 0.02306971  max resid 0.07847087 
    ## Run 29 stress 0.1138394 
    ## Run 30 stress 0.1079381 
    ## Run 31 stress 0.1050338 
    ## ... Procrustes: rmse 0.008878153  max resid 0.02485679 
    ## Run 32 stress 0.1135598 
    ## Run 33 stress 0.1049903 
    ## ... New best solution
    ## ... Procrustes: rmse 4.713729e-05  max resid 0.0001338538 
    ## ... Similar to previous best
    ## Run 34 stress 0.1050144 
    ## ... Procrustes: rmse 0.006660525  max resid 0.02376584 
    ## Run 35 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307792  max resid 0.07849244 
    ## Run 36 stress 0.1083988 
    ## Run 37 stress 0.1138984 
    ## Run 38 stress 0.1174544 
    ## Run 39 stress 0.1063128 
    ## Run 40 stress 0.1097988 
    ## Run 41 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307808  max resid 0.07848025 
    ## Run 42 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307566  max resid 0.07848004 
    ## Run 43 stress 0.1138397 
    ## Run 44 stress 0.1051923 
    ## ... Procrustes: rmse 0.0221515  max resid 0.07903013 
    ## Run 45 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214741  max resid 0.07899934 
    ## Run 46 stress 0.1051923 
    ## ... Procrustes: rmse 0.022148  max resid 0.07900622 
    ## Run 47 stress 0.1050338 
    ## ... Procrustes: rmse 0.008888813  max resid 0.02495194 
    ## Run 48 stress 0.1051402 
    ## ... Procrustes: rmse 0.02308499  max resid 0.07850867 
    ## Run 49 stress 0.1135598 
    ## Run 50 stress 0.1050144 
    ## ... Procrustes: rmse 0.006657973  max resid 0.02374481 
    ## Run 51 stress 0.1088548 
    ## Run 52 stress 0.1082722 
    ## Run 53 stress 0.113594 
    ## Run 54 stress 0.1063128 
    ## Run 55 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307682  max resid 0.07848005 
    ## Run 56 stress 0.1082722 
    ## Run 57 stress 0.1050338 
    ## ... Procrustes: rmse 0.008888546  max resid 0.02496593 
    ## Run 58 stress 0.1084964 
    ## Run 59 stress 0.1063128 
    ## Run 60 stress 0.1135596 
    ## Run 61 stress 0.1051924 
    ## ... Procrustes: rmse 0.02215756  max resid 0.0790719 
    ## Run 62 stress 0.1049903 
    ## ... Procrustes: rmse 5.15851e-05  max resid 0.0001294736 
    ## ... Similar to previous best
    ## Run 63 stress 0.1136647 
    ## Run 64 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307464  max resid 0.07850398 
    ## Run 65 stress 0.1063128 
    ## Run 66 stress 0.1136646 
    ## Run 67 stress 0.113899 
    ## Run 68 stress 0.1051923 
    ## ... Procrustes: rmse 0.0221485  max resid 0.07900731 
    ## Run 69 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214408  max resid 0.07897705 
    ## Run 70 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214517  max resid 0.07898539 
    ## Run 71 stress 0.116054 
    ## Run 72 stress 0.1097987 
    ## Run 73 stress 0.1050338 
    ## ... Procrustes: rmse 0.008889176  max resid 0.02496781 
    ## Run 74 stress 0.1138398 
    ## Run 75 stress 0.1051402 
    ## ... Procrustes: rmse 0.02308043  max resid 0.07850338 
    ## Run 76 stress 0.1049903 
    ## ... Procrustes: rmse 3.670096e-05  max resid 8.933527e-05 
    ## ... Similar to previous best
    ## Run 77 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214803  max resid 0.07900458 
    ## Run 78 stress 0.1051402 
    ## ... Procrustes: rmse 0.02308526  max resid 0.07852018 
    ## Run 79 stress 0.1084965 
    ## Run 80 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307619  max resid 0.07848637 
    ## Run 81 stress 0.1063128 
    ## Run 82 stress 0.1050144 
    ## ... Procrustes: rmse 0.006658092  max resid 0.02374901 
    ## Run 83 stress 0.1136646 
    ## Run 84 stress 0.1063128 
    ## Run 85 stress 0.113594 
    ## Run 86 stress 0.1050144 
    ## ... Procrustes: rmse 0.006662066  max resid 0.02377381 
    ## Run 87 stress 0.1138393 
    ## Run 88 stress 0.1050338 
    ## ... Procrustes: rmse 0.008890542  max resid 0.02496576 
    ## Run 89 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 0.00671726  max resid 0.02572039 
    ## Run 90 stress 0.10835 
    ## Run 91 stress 0.1050338 
    ## ... Procrustes: rmse 0.006252567  max resid 0.02184525 
    ## Run 92 stress 0.1135597 
    ## Run 93 stress 0.1084965 
    ## Run 94 stress 0.1050144 
    ## ... Procrustes: rmse 0.009651259  max resid 0.02603771 
    ## Run 95 stress 0.1160547 
    ## Run 96 stress 0.1063128 
    ## Run 97 stress 0.1135597 
    ## Run 98 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218072  max resid 0.0768805 
    ## Run 99 stress 0.1063128 
    ## Run 100 stress 0.1082722 
    ## Run 101 stress 0.1136647 
    ## Run 102 stress 0.1063128 
    ## Run 103 stress 0.1174532 
    ## Run 104 stress 0.1088547 
    ## Run 105 stress 0.1063128 
    ## Run 106 stress 0.1138394 
    ## Run 107 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284609  max resid 0.07734067 
    ## Run 108 stress 0.1136647 
    ## Run 109 stress 0.1136647 
    ## Run 110 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181211  max resid 0.07690355 
    ## Run 111 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284441  max resid 0.07733335 
    ## Run 112 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284434  max resid 0.0773212 
    ## Run 113 stress 0.1135598 
    ## Run 114 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284846  max resid 0.07731794 
    ## Run 115 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181411  max resid 0.07690347 
    ## Run 116 stress 0.113594 
    ## Run 117 stress 0.10835 
    ## Run 118 stress 0.1138393 
    ## Run 119 stress 0.1063128 
    ## Run 120 stress 0.1050144 
    ## ... Procrustes: rmse 0.009652152  max resid 0.02604522 
    ## Run 121 stress 0.1063128 
    ## Run 122 stress 0.1050338 
    ## ... Procrustes: rmse 0.006258646  max resid 0.02187759 
    ## Run 123 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284577  max resid 0.07727136 
    ## Run 124 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284583  max resid 0.0773459 
    ## Run 125 stress 0.1082722 
    ## Run 126 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181302  max resid 0.07689951 
    ## Run 127 stress 0.1138394 
    ## Run 128 stress 0.1088547 
    ## Run 129 stress 0.10835 
    ## Run 130 stress 0.110272 
    ## Run 131 stress 0.1083503 
    ## Run 132 stress 0.1050144 
    ## ... Procrustes: rmse 0.009655953  max resid 0.02607541 
    ## Run 133 stress 0.1084964 
    ## Run 134 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181241  max resid 0.07689851 
    ## Run 135 stress 0.1160539 
    ## Run 136 stress 0.1058423 
    ## Run 137 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181878  max resid 0.07692078 
    ## Run 138 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228448  max resid 0.07733361 
    ## Run 139 stress 0.1102725 
    ## Run 140 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284716  max resid 0.07731794 
    ## Run 141 stress 0.1136646 
    ## Run 142 stress 0.1079381 
    ## Run 143 stress 0.1079382 
    ## Run 144 stress 0.1083501 
    ## Run 145 stress 0.1136646 
    ## Run 146 stress 0.1063128 
    ## Run 147 stress 0.1174536 
    ## Run 148 stress 0.1088547 
    ## Run 149 stress 0.1083501 
    ## Run 150 stress 0.1165302 
    ## Run 151 stress 0.1063128 
    ## Run 152 stress 0.1063128 
    ## Run 153 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284529  max resid 0.07733344 
    ## Run 154 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284731  max resid 0.07737654 
    ## Run 155 stress 0.1160548 
    ## Run 156 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180825  max resid 0.07688145 
    ## Run 157 stress 0.1063128 
    ## Run 158 stress 0.113594 
    ## Run 159 stress 0.1088548 
    ## Run 160 stress 0.1097987 
    ## Run 161 stress 0.1097987 
    ## Run 162 stress 0.1079383 
    ## Run 163 stress 0.113594 
    ## Run 164 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180847  max resid 0.07688406 
    ## Run 165 stress 0.1088547 
    ## Run 166 stress 0.1097988 
    ## Run 167 stress 0.1083501 
    ## Run 168 stress 0.1135597 
    ## Run 169 stress 0.1136646 
    ## Run 170 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284464  max resid 0.07733437 
    ## Run 171 stress 0.1063128 
    ## Run 172 stress 0.1082722 
    ## Run 173 stress 0.1079382 
    ## Run 174 stress 0.1082722 
    ## Run 175 stress 0.105871 
    ## Run 176 stress 0.1084964 
    ## Run 177 stress 0.1083995 
    ## Run 178 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284553  max resid 0.07733491 
    ## Run 179 stress 0.1084965 
    ## Run 180 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284462  max resid 0.07731947 
    ## Run 181 stress 0.1136646 
    ## Run 182 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180936  max resid 0.07688751 
    ## Run 183 stress 0.1063128 
    ## Run 184 stress 0.1135598 
    ## Run 185 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284702  max resid 0.0772884 
    ## Run 186 stress 0.1136646 
    ## Run 187 stress 0.1136646 
    ## Run 188 stress 0.1063128 
    ## Run 189 stress 0.1135597 
    ## Run 190 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180845  max resid 0.07688323 
    ## Run 191 stress 0.1083502 
    ## Run 192 stress 0.1050144 
    ## ... Procrustes: rmse 0.009666233  max resid 0.02613468 
    ## Run 193 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218087  max resid 0.07686421 
    ## Run 194 stress 0.113594 
    ## Run 195 stress 0.1136646 
    ## Run 196 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284159  max resid 0.07734073 
    ## Run 197 stress 0.1160555 
    ## Run 198 stress 0.1135597 
    ## Run 199 stress 0.1058423 
    ## Run 200 stress 0.1083503 
    ## Run 201 stress 0.1083992 
    ## Run 202 stress 0.113594 
    ## Run 203 stress 0.1082723 
    ## Run 204 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284892  max resid 0.07737131 
    ## Run 205 stress 0.1097987 
    ## Run 206 stress 0.1050338 
    ## ... Procrustes: rmse 0.006255554  max resid 0.02186335 
    ## Run 207 stress 0.1079381 
    ## Run 208 stress 0.1088548 
    ## Run 209 stress 0.1135597 
    ## Run 210 stress 0.1063128 
    ## Run 211 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181519  max resid 0.0769097 
    ## Run 212 stress 0.1083499 
    ## Run 213 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181082  max resid 0.07689143 
    ## Run 214 stress 0.1136647 
    ## Run 215 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 3.657369e-05  max resid 0.0001041828 
    ## ... Similar to previous best
    ## Run 216 stress 0.1051403 
    ## ... Procrustes: rmse 0.0218012  max resid 0.07684724 
    ## Run 217 stress 0.1079381 
    ## Run 218 stress 0.113594 
    ## Run 219 stress 0.113594 
    ## Run 220 stress 0.1084965 
    ## Run 221 stress 0.1135599 
    ## Run 222 stress 0.1082722 
    ## Run 223 stress 0.1135596 
    ## Run 224 stress 0.113594 
    ## Run 225 stress 0.1136646 
    ## Run 226 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228374  max resid 0.07727802 
    ## Run 227 stress 0.1174536 
    ## Run 228 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283735  max resid 0.07727991 
    ## Run 229 stress 0.1063128 
    ## Run 230 stress 0.1138396 
    ## Run 231 stress 0.1135597 
    ## Run 232 stress 0.1050338 
    ## ... Procrustes: rmse 0.0062568  max resid 0.0218524 
    ## Run 233 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180349  max resid 0.07685677 
    ## Run 234 stress 0.1135941 
    ## Run 235 stress 0.1097987 
    ## Run 236 stress 0.1138984 
    ## Run 237 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283537  max resid 0.07725666 
    ## Run 238 stress 0.1160548 
    ## Run 239 stress 0.1050338 
    ## ... Procrustes: rmse 0.006254569  max resid 0.02183821 
    ## Run 240 stress 0.1097988 
    ## Run 241 stress 0.1160538 
    ## Run 242 stress 0.1097988 
    ## Run 243 stress 0.1136647 
    ## Run 244 stress 0.1063128 
    ## Run 245 stress 0.1049903 
    ## ... Procrustes: rmse 0.006737538  max resid 0.02581911 
    ## Run 246 stress 0.1079381 
    ## Run 247 stress 0.1088547 
    ## Run 248 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179869  max resid 0.07684177 
    ## Run 249 stress 0.1082722 
    ## Run 250 stress 0.1160546 
    ## Run 251 stress 0.1136648 
    ## Run 252 stress 0.1082722 
    ## Run 253 stress 0.1049903 
    ## ... Procrustes: rmse 0.006734299  max resid 0.02580668 
    ## Run 254 stress 0.1138396 
    ## Run 255 stress 0.1160547 
    ## Run 256 stress 0.1049903 
    ## ... Procrustes: rmse 0.006735814  max resid 0.02581462 
    ## Run 257 stress 0.1050144 
    ## ... Procrustes: rmse 0.009675316  max resid 0.02615485 
    ## Run 258 stress 0.1063128 
    ## Run 259 stress 0.1050144 
    ## ... Procrustes: rmse 0.009654943  max resid 0.02606708 
    ## Run 260 stress 0.1084965 
    ## Run 261 stress 0.1050144 
    ## ... Procrustes: rmse 0.009669126  max resid 0.02610388 
    ## Run 262 stress 0.1088547 
    ## Run 263 stress 0.1135598 
    ## Run 264 stress 0.1174531 
    ## Run 265 stress 0.1050144 
    ## ... Procrustes: rmse 0.009675005  max resid 0.02614146 
    ## Run 266 stress 0.1160546 
    ## Run 267 stress 0.1136646 
    ## Run 268 stress 0.1097988 
    ## Run 269 stress 0.10835 
    ## Run 270 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283757  max resid 0.07729638 
    ## Run 271 stress 0.1135596 
    ## Run 272 stress 0.1058423 
    ## Run 273 stress 0.1058423 
    ## Run 274 stress 0.1135597 
    ## Run 275 stress 0.1097987 
    ## Run 276 stress 0.1138392 
    ## Run 277 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179731  max resid 0.07683615 
    ## Run 278 stress 0.1058423 
    ## Run 279 stress 0.1084965 
    ## Run 280 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283659  max resid 0.07727798 
    ## Run 281 stress 0.1049903 
    ## ... Procrustes: rmse 0.006723531  max resid 0.02575193 
    ## Run 282 stress 0.1051923 
    ## ... Procrustes: rmse 0.022833  max resid 0.07728824 
    ## Run 283 stress 0.1063128 
    ## Run 284 stress 0.1309382 
    ## Run 285 stress 0.1051924 
    ## ... Procrustes: rmse 0.0228353  max resid 0.07732443 
    ## Run 286 stress 0.1050144 
    ## ... Procrustes: rmse 0.009674109  max resid 0.02614029 
    ## Run 287 stress 0.1049791 
    ## ... Procrustes: rmse 4.381908e-05  max resid 0.0001202389 
    ## ... Similar to previous best
    ## Run 288 stress 0.1051924 
    ## ... Procrustes: rmse 0.02283981  max resid 0.07729281 
    ## Run 289 stress 0.1138395 
    ## Run 290 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180081  max resid 0.0768455 
    ## Run 291 stress 0.1050338 
    ## ... Procrustes: rmse 0.006255824  max resid 0.02184866 
    ## Run 292 stress 0.1135599 
    ## Run 293 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283653  max resid 0.07727693 
    ## Run 294 stress 0.1135597 
    ## Run 295 stress 0.1084964 
    ## Run 296 stress 0.1135599 
    ## Run 297 stress 0.1135597 
    ## Run 298 stress 0.1063128 
    ## Run 299 stress 0.1084965 
    ## Run 300 stress 0.1079381 
    ## Run 301 stress 0.1049903 
    ## ... Procrustes: rmse 0.00676278  max resid 0.02593681 
    ## Run 302 stress 0.1135598 
    ## Run 303 stress 0.1136647 
    ## Run 304 stress 0.1088547 
    ## Run 305 stress 0.117453 
    ## Run 306 stress 0.1084965 
    ## Run 307 stress 0.1063128 
    ## Run 308 stress 0.1088547 
    ## Run 309 stress 0.1097987 
    ## Run 310 stress 0.1138985 
    ## Run 311 stress 0.1050144 
    ## ... Procrustes: rmse 0.009683737  max resid 0.02618364 
    ## Run 312 stress 0.113594 
    ## Run 313 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179154  max resid 0.07681606 
    ## Run 314 stress 0.113594 
    ## Run 315 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179803  max resid 0.07683922 
    ## Run 316 stress 0.1309391 
    ## Run 317 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179827  max resid 0.07683578 
    ## Run 318 stress 0.1083501 
    ## Run 319 stress 0.1058423 
    ## Run 320 stress 0.107938 
    ## Run 321 stress 0.1135598 
    ## Run 322 stress 0.1049791 
    ## ... Procrustes: rmse 3.160085e-05  max resid 8.80804e-05 
    ## ... Similar to previous best
    ## Run 323 stress 0.1174531 
    ## Run 324 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179513  max resid 0.07682796 
    ## Run 325 stress 0.1050144 
    ## ... Procrustes: rmse 0.00967734  max resid 0.02616184 
    ## Run 326 stress 0.1051402 
    ## ... Procrustes: rmse 0.0217968  max resid 0.07683469 
    ## Run 327 stress 0.1088547 
    ## Run 328 stress 0.1097988 
    ## Run 329 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180178  max resid 0.07684985 
    ## Run 330 stress 0.1063128 
    ## Run 331 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179885  max resid 0.07684066 
    ## Run 332 stress 0.11356 
    ## Run 333 stress 0.1135599 
    ## Run 334 stress 0.1097987 
    ## Run 335 stress 0.1050144 
    ## ... Procrustes: rmse 0.009667975  max resid 0.02610066 
    ## Run 336 stress 0.1050144 
    ## ... Procrustes: rmse 0.009679595  max resid 0.02616893 
    ## Run 337 stress 0.1083995 
    ## Run 338 stress 0.1079382 
    ## Run 339 stress 0.113594 
    ## Run 340 stress 0.1063129 
    ## Run 341 stress 0.1050338 
    ## ... Procrustes: rmse 0.006257546  max resid 0.02185142 
    ## Run 342 stress 0.1050144 
    ## ... Procrustes: rmse 0.009662451  max resid 0.02606441 
    ## Run 343 stress 0.1082722 
    ## Run 344 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179773  max resid 0.07683193 
    ## Run 345 stress 0.1079381 
    ## Run 346 stress 0.1138395 
    ## Run 347 stress 0.1102723 
    ## Run 348 stress 0.1050338 
    ## ... Procrustes: rmse 0.006263189  max resid 0.0218969 
    ## Run 349 stress 0.1136646 
    ## Run 350 stress 0.1160552 
    ## Run 351 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179875  max resid 0.07684078 
    ## Run 352 stress 0.116055 
    ## Run 353 stress 0.1135596 
    ## Run 354 stress 0.1138395 
    ## Run 355 stress 0.1138392 
    ## Run 356 stress 0.1063128 
    ## Run 357 stress 0.1160545 
    ## Run 358 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179772  max resid 0.0768306 
    ## Run 359 stress 0.1079382 
    ## Run 360 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180038  max resid 0.07684535 
    ## Run 361 stress 0.1050338 
    ## ... Procrustes: rmse 0.006257734  max resid 0.02185704 
    ## Run 362 stress 0.1084964 
    ## Run 363 stress 0.1138395 
    ## Run 364 stress 0.1102722 
    ## Run 365 stress 0.1097988 
    ## Run 366 stress 0.1049903 
    ## ... Procrustes: rmse 0.006778948  max resid 0.02601263 
    ## Run 367 stress 0.1063128 
    ## Run 368 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180105  max resid 0.07685415 
    ## Run 369 stress 0.1084964 
    ## Run 370 stress 0.1097989 
    ## Run 371 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283975  max resid 0.07729066 
    ## Run 372 stress 0.1084965 
    ## Run 373 stress 0.1138392 
    ## Run 374 stress 0.1136646 
    ## Run 375 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180455  max resid 0.07686149 
    ## Run 376 stress 0.1063128 
    ## Run 377 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.742397e-05  max resid 5.451117e-05 
    ## ... Similar to previous best
    ## Run 378 stress 0.1063128 
    ## Run 379 stress 0.1102723 
    ## Run 380 stress 0.1102722 
    ## Run 381 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179903  max resid 0.07684718 
    ## Run 382 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283721  max resid 0.07728484 
    ## Run 383 stress 0.1049903 
    ## ... Procrustes: rmse 0.006702598  max resid 0.02568468 
    ## Run 384 stress 0.1097988 
    ## Run 385 stress 0.1082722 
    ## Run 386 stress 0.1049791 
    ## ... Procrustes: rmse 2.034587e-05  max resid 6.265083e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.1049791 
    ## ... Procrustes: rmse 3.644875e-05  max resid 9.889179e-05 
    ## ... Similar to previous best
    ## Run 388 stress 0.11356 
    ## Run 389 stress 0.1136647 
    ## Run 390 stress 0.1138393 
    ## Run 391 stress 0.1082722 
    ## Run 392 stress 0.1063128 
    ## Run 393 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283771  max resid 0.0773028 
    ## Run 394 stress 0.1083994 
    ## Run 395 stress 0.1136646 
    ## Run 396 stress 0.1138394 
    ## Run 397 stress 0.1082722 
    ## Run 398 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179978  max resid 0.07683313 
    ## Run 399 stress 0.1135598 
    ## Run 400 stress 0.1084965 
    ## Run 401 stress 0.10835 
    ## Run 402 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180199  max resid 0.0768579 
    ## Run 403 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283724  max resid 0.07729912 
    ## Run 404 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283624  max resid 0.0772823 
    ## Run 405 stress 0.1082722 
    ## Run 406 stress 0.1160539 
    ## Run 407 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180114  max resid 0.07685708 
    ## Run 408 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179873  max resid 0.07682096 
    ## Run 409 stress 0.1049791 
    ## ... Procrustes: rmse 2.550647e-05  max resid 7.788268e-05 
    ## ... Similar to previous best
    ## Run 410 stress 0.1102724 
    ## Run 411 stress 0.1084965 
    ## Run 412 stress 0.1050338 
    ## ... Procrustes: rmse 0.006254946  max resid 0.02185643 
    ## Run 413 stress 0.1309391 
    ## Run 414 stress 0.1135597 
    ## Run 415 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283657  max resid 0.07727944 
    ## Run 416 stress 0.1050338 
    ## ... Procrustes: rmse 0.006253128  max resid 0.02185053 
    ## Run 417 stress 0.1136646 
    ## Run 418 stress 0.1082722 
    ## Run 419 stress 0.1082723 
    ## Run 420 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283981  max resid 0.07730693 
    ## Run 421 stress 0.1049903 
    ## ... Procrustes: rmse 0.006778568  max resid 0.0260361 
    ## Run 422 stress 0.1063128 
    ## Run 423 stress 0.1063129 
    ## Run 424 stress 0.1084965 
    ## Run 425 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180162  max resid 0.0768521 
    ## Run 426 stress 0.1174543 
    ## Run 427 stress 0.10835 
    ## Run 428 stress 0.1136647 
    ## Run 429 stress 0.1050338 
    ## ... Procrustes: rmse 0.006259085  max resid 0.02187375 
    ## Run 430 stress 0.1049903 
    ## ... Procrustes: rmse 0.006719347  max resid 0.02576571 
    ## Run 431 stress 0.1309396 
    ## Run 432 stress 0.1063128 
    ## Run 433 stress 0.1050144 
    ## ... Procrustes: rmse 0.009656792  max resid 0.0260572 
    ## Run 434 stress 0.1084964 
    ## Run 435 stress 0.1050338 
    ## ... Procrustes: rmse 0.00626016  max resid 0.02189254 
    ## Run 436 stress 0.1160542 
    ## Run 437 stress 0.1097989 
    ## Run 438 stress 0.1088547 
    ## Run 439 stress 0.113594 
    ## Run 440 stress 0.1084965 
    ## Run 441 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283776  max resid 0.07728978 
    ## Run 442 stress 0.1083994 
    ## Run 443 stress 0.1088547 
    ## Run 444 stress 0.10835 
    ## Run 445 stress 0.1063128 
    ## Run 446 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283701  max resid 0.07729114 
    ## Run 447 stress 0.1084965 
    ## Run 448 stress 0.1136646 
    ## Run 449 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283417  max resid 0.07725988 
    ## Run 450 stress 0.1063128 
    ## Run 451 stress 0.1082722 
    ## Run 452 stress 0.1063128 
    ## Run 453 stress 0.1160551 
    ## Run 454 stress 0.1084964 
    ## Run 455 stress 0.1097987 
    ## Run 456 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283581  max resid 0.07727242 
    ## Run 457 stress 0.1063128 
    ## Run 458 stress 0.1050144 
    ## ... Procrustes: rmse 0.009659912  max resid 0.02607277 
    ## Run 459 stress 0.1084965 
    ## Run 460 stress 0.1102721 
    ## Run 461 stress 0.1135597 
    ## Run 462 stress 0.1079384 
    ## Run 463 stress 0.1138984 
    ## Run 464 stress 0.1063128 
    ## Run 465 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179805  max resid 0.07681745 
    ## Run 466 stress 0.1135598 
    ## Run 467 stress 0.1097988 
    ## Run 468 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180176  max resid 0.07685659 
    ## Run 469 stress 0.1082722 
    ## Run 470 stress 0.1135597 
    ## Run 471 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284106  max resid 0.07729394 
    ## Run 472 stress 0.1135598 
    ## Run 473 stress 0.1058423 
    ## Run 474 stress 0.1063128 
    ## Run 475 stress 0.1138984 
    ## Run 476 stress 0.1136647 
    ## Run 477 stress 0.1049791 
    ## ... Procrustes: rmse 2.873227e-05  max resid 8.557167e-05 
    ## ... Similar to previous best
    ## Run 478 stress 0.1050338 
    ## ... Procrustes: rmse 0.006256009  max resid 0.02186282 
    ## Run 479 stress 0.1084964 
    ## Run 480 stress 0.1083501 
    ## Run 481 stress 0.1079384 
    ## Run 482 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283735  max resid 0.07728869 
    ## Run 483 stress 0.1084966 
    ## Run 484 stress 0.1097987 
    ## Run 485 stress 0.1102721 
    ## Run 486 stress 0.1138392 
    ## Run 487 stress 0.1063128 
    ## Run 488 stress 0.107938 
    ## Run 489 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180323  max resid 0.07686177 
    ## Run 490 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218002  max resid 0.07685156 
    ## Run 491 stress 0.1084965 
    ## Run 492 stress 0.1049903 
    ## ... Procrustes: rmse 0.006699933  max resid 0.02567299 
    ## Run 493 stress 0.1079382 
    ## Run 494 stress 0.1102722 
    ## Run 495 stress 0.1102722 
    ## Run 496 stress 0.1051403 
    ## ... Procrustes: rmse 0.02179187  max resid 0.07682679 
    ## Run 497 stress 0.1102722 
    ## Run 498 stress 0.1088547 
    ## Run 499 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228369  max resid 0.07728468 
    ## Run 500 stress 0.1084964 
    ## *** Best solution repeated 5 times

``` r
# Surveyed sites 
SD_beta_NMDS <- metaMDS(SD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09021173 
    ## Run 1 stress 0.09087341 
    ## Run 2 stress 0.09099536 
    ## Run 3 stress 0.08926064 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04364507  max resid 0.1169943 
    ## Run 4 stress 0.1092127 
    ## Run 5 stress 0.09503437 
    ## Run 6 stress 0.08946655 
    ## ... Procrustes: rmse 0.033208  max resid 0.1163238 
    ## Run 7 stress 0.1076304 
    ## Run 8 stress 0.08926074 
    ## ... Procrustes: rmse 0.0001473169  max resid 0.0004564591 
    ## ... Similar to previous best
    ## Run 9 stress 0.1111386 
    ## Run 10 stress 0.1111388 
    ## Run 11 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183077  max resid 0.03971665 
    ## Run 12 stress 0.08926068 
    ## ... Procrustes: rmse 5.64398e-05  max resid 0.0001846404 
    ## ... Similar to previous best
    ## Run 13 stress 0.1060091 
    ## Run 14 stress 0.1064429 
    ## Run 15 stress 0.08938968 
    ## ... Procrustes: rmse 0.03594803  max resid 0.1182765 
    ## Run 16 stress 0.1060431 
    ## Run 17 stress 0.09503471 
    ## Run 18 stress 0.08938962 
    ## ... Procrustes: rmse 0.03594994  max resid 0.1182844 
    ## Run 19 stress 0.0901871 
    ## Run 20 stress 0.108839 
    ## Run 21 stress 0.09594312 
    ## Run 22 stress 0.09503432 
    ## Run 23 stress 0.1061301 
    ## Run 24 stress 0.1093318 
    ## Run 25 stress 0.08938554 
    ## ... Procrustes: rmse 0.01176518  max resid 0.03974109 
    ## Run 26 stress 0.1092364 
    ## Run 27 stress 0.1063191 
    ## Run 28 stress 0.109188 
    ## Run 29 stress 0.1093478 
    ## Run 30 stress 0.09044608 
    ## Run 31 stress 0.1052647 
    ## Run 32 stress 0.1062565 
    ## Run 33 stress 0.1076303 
    ## Run 34 stress 0.08946346 
    ## ... Procrustes: rmse 0.03748085  max resid 0.1177197 
    ## Run 35 stress 0.0903913 
    ## Run 36 stress 0.0903914 
    ## Run 37 stress 0.08926082 
    ## ... Procrustes: rmse 0.0002062901  max resid 0.0006456972 
    ## ... Similar to previous best
    ## Run 38 stress 0.09039138 
    ## Run 39 stress 0.08926083 
    ## ... Procrustes: rmse 0.0002127131  max resid 0.0006913684 
    ## ... Similar to previous best
    ## Run 40 stress 0.1111381 
    ## Run 41 stress 0.09109108 
    ## Run 42 stress 0.08938541 
    ## ... Procrustes: rmse 0.01184378  max resid 0.03968657 
    ## Run 43 stress 0.09018717 
    ## Run 44 stress 0.1052647 
    ## Run 45 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595448  max resid 0.1182909 
    ## Run 46 stress 0.08926085 
    ## ... Procrustes: rmse 0.000241355  max resid 0.0007828007 
    ## ... Similar to previous best
    ## Run 47 stress 0.08926105 
    ## ... Procrustes: rmse 0.0003590552  max resid 0.001187653 
    ## ... Similar to previous best
    ## Run 48 stress 0.08926099 
    ## ... Procrustes: rmse 0.0003292286  max resid 0.001083765 
    ## ... Similar to previous best
    ## Run 49 stress 0.09130088 
    ## Run 50 stress 0.08938967 
    ## ... Procrustes: rmse 0.03593697  max resid 0.1182644 
    ## Run 51 stress 0.08926096 
    ## ... Procrustes: rmse 0.000307548  max resid 0.001012665 
    ## ... Similar to previous best
    ## Run 52 stress 0.1075422 
    ## Run 53 stress 0.09039135 
    ## Run 54 stress 0.1061303 
    ## Run 55 stress 0.08938544 
    ## ... Procrustes: rmse 0.01179271  max resid 0.03973113 
    ## Run 56 stress 0.09503429 
    ## Run 57 stress 0.08938963 
    ## ... Procrustes: rmse 0.03595006  max resid 0.1182857 
    ## Run 58 stress 0.09130087 
    ## Run 59 stress 0.08951722 
    ## ... Procrustes: rmse 0.03504514  max resid 0.1156255 
    ## Run 60 stress 0.0892611 
    ## ... Procrustes: rmse 0.0003834738  max resid 0.001204989 
    ## ... Similar to previous best
    ## Run 61 stress 0.09099543 
    ## Run 62 stress 0.09021171 
    ## Run 63 stress 0.09178294 
    ## Run 64 stress 0.1052648 
    ## Run 65 stress 0.1056895 
    ## Run 66 stress 0.08938971 
    ## ... Procrustes: rmse 0.03592941  max resid 0.118252 
    ## Run 67 stress 0.1103717 
    ## Run 68 stress 0.1060429 
    ## Run 69 stress 0.09087352 
    ## Run 70 stress 0.1108324 
    ## Run 71 stress 0.1052648 
    ## Run 72 stress 0.1101446 
    ## Run 73 stress 0.09594324 
    ## Run 74 stress 0.09039131 
    ## Run 75 stress 0.08926065 
    ## ... Procrustes: rmse 2.993389e-05  max resid 9.226518e-05 
    ## ... Similar to previous best
    ## Run 76 stress 0.08946329 
    ## ... Procrustes: rmse 0.03752437  max resid 0.1177854 
    ## Run 77 stress 0.1075421 
    ## Run 78 stress 0.1111369 
    ## Run 79 stress 0.1089905 
    ## Run 80 stress 0.08938541 
    ## ... Procrustes: rmse 0.01185614  max resid 0.03967376 
    ## Run 81 stress 0.1092363 
    ## Run 82 stress 0.09021165 
    ## Run 83 stress 0.1052647 
    ## Run 84 stress 0.08938547 
    ## ... Procrustes: rmse 0.01186232  max resid 0.03966994 
    ## Run 85 stress 0.09503453 
    ## Run 86 stress 0.08946329 
    ## ... Procrustes: rmse 0.03751441  max resid 0.1177684 
    ## Run 87 stress 0.1062569 
    ## Run 88 stress 0.08926105 
    ## ... Procrustes: rmse 0.0003559391  max resid 0.001179469 
    ## ... Similar to previous best
    ## Run 89 stress 0.09099604 
    ## Run 90 stress 0.09088615 
    ## Run 91 stress 0.09503429 
    ## Run 92 stress 0.08926066 
    ## ... Procrustes: rmse 4.338229e-05  max resid 0.0001408174 
    ## ... Similar to previous best
    ## Run 93 stress 0.1056907 
    ## Run 94 stress 0.1056894 
    ## Run 95 stress 0.1074431 
    ## Run 96 stress 0.08946338 
    ## ... Procrustes: rmse 0.03748909  max resid 0.1177326 
    ## Run 97 stress 0.08938537 
    ## ... Procrustes: rmse 0.01182952  max resid 0.03970359 
    ## Run 98 stress 0.1066199 
    ## Run 99 stress 0.09044608 
    ## Run 100 stress 0.09039133 
    ## Run 101 stress 0.1087567 
    ## Run 102 stress 0.1092362 
    ## Run 103 stress 0.105265 
    ## Run 104 stress 0.09178282 
    ## Run 105 stress 0.08938544 
    ## ... Procrustes: rmse 0.01177092  max resid 0.03969462 
    ## Run 106 stress 0.1087568 
    ## Run 107 stress 0.1067503 
    ## Run 108 stress 0.0950344 
    ## Run 109 stress 0.1087864 
    ## Run 110 stress 0.09503463 
    ## Run 111 stress 0.1092128 
    ## Run 112 stress 0.08946342 
    ## ... Procrustes: rmse 0.03748258  max resid 0.1177199 
    ## Run 113 stress 0.08938963 
    ## ... Procrustes: rmse 0.03595525  max resid 0.1182932 
    ## Run 114 stress 0.1091237 
    ## Run 115 stress 0.09088608 
    ## Run 116 stress 0.1105346 
    ## Run 117 stress 0.1092092 
    ## Run 118 stress 0.10569 
    ## Run 119 stress 0.08946343 
    ## ... Procrustes: rmse 0.03751078  max resid 0.1177586 
    ## Run 120 stress 0.0917828 
    ## Run 121 stress 0.1079014 
    ## Run 122 stress 0.09503442 
    ## Run 123 stress 0.08938963 
    ## ... Procrustes: rmse 0.03597661  max resid 0.1183254 
    ## Run 124 stress 0.08951727 
    ## ... Procrustes: rmse 0.03504313  max resid 0.115645 
    ## Run 125 stress 0.09592174 
    ## Run 126 stress 0.09099533 
    ## Run 127 stress 0.08938961 
    ## ... Procrustes: rmse 0.03596037  max resid 0.1183006 
    ## Run 128 stress 0.08938545 
    ## ... Procrustes: rmse 0.0117886  max resid 0.03973292 
    ## Run 129 stress 0.09087341 
    ## Run 130 stress 0.09592142 
    ## Run 131 stress 0.1076304 
    ## Run 132 stress 0.09592134 
    ## Run 133 stress 0.1060433 
    ## Run 134 stress 0.08926066 
    ## ... Procrustes: rmse 0.0001851134  max resid 0.0004891975 
    ## ... Similar to previous best
    ## Run 135 stress 0.08926067 
    ## ... Procrustes: rmse 0.0001904567  max resid 0.0004699197 
    ## ... Similar to previous best
    ## Run 136 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595782  max resid 0.1182971 
    ## Run 137 stress 0.09018708 
    ## Run 138 stress 0.1052648 
    ## Run 139 stress 0.08946329 
    ## ... Procrustes: rmse 0.03751452  max resid 0.1177713 
    ## Run 140 stress 0.08938544 
    ## ... Procrustes: rmse 0.01188616  max resid 0.03972581 
    ## Run 141 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183351  max resid 0.03971173 
    ## Run 142 stress 0.1092096 
    ## Run 143 stress 0.08926074 
    ## ... Procrustes: rmse 0.0001463764  max resid 0.0004607561 
    ## ... Similar to previous best
    ## Run 144 stress 0.08938538 
    ## ... Procrustes: rmse 0.01184318  max resid 0.03969617 
    ## Run 145 stress 0.09775567 
    ## Run 146 stress 0.1074432 
    ## Run 147 stress 0.1101446 
    ## Run 148 stress 0.09044617 
    ## Run 149 stress 0.09503427 
    ## Run 150 stress 0.1071321 
    ## Run 151 stress 0.09109119 
    ## Run 152 stress 0.08946652 
    ## ... Procrustes: rmse 0.03324045  max resid 0.1163791 
    ## Run 153 stress 0.10613 
    ## Run 154 stress 0.08938569 
    ## ... Procrustes: rmse 0.01171163  max resid 0.03960879 
    ## Run 155 stress 0.1074431 
    ## Run 156 stress 0.08938976 
    ## ... Procrustes: rmse 0.03592085  max resid 0.1182375 
    ## Run 157 stress 0.0892608 
    ## ... Procrustes: rmse 0.0001982318  max resid 0.0006242368 
    ## ... Similar to previous best
    ## Run 158 stress 0.1056897 
    ## Run 159 stress 0.1061297 
    ## Run 160 stress 0.08938965 
    ## ... Procrustes: rmse 0.03594234  max resid 0.118273 
    ## Run 161 stress 0.09044628 
    ## Run 162 stress 0.09503442 
    ## Run 163 stress 0.08946331 
    ## ... Procrustes: rmse 0.0375056  max resid 0.1177558 
    ## Run 164 stress 0.09039135 
    ## Run 165 stress 0.1061296 
    ## Run 166 stress 0.08938537 
    ## ... Procrustes: rmse 0.01183264  max resid 0.03968759 
    ## Run 167 stress 0.09087343 
    ## Run 168 stress 0.1052653 
    ## Run 169 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595531  max resid 0.118292 
    ## Run 170 stress 0.1076305 
    ## Run 171 stress 0.1096136 
    ## Run 172 stress 0.08926086 
    ## ... Procrustes: rmse 0.0003691233  max resid 0.001095087 
    ## ... Similar to previous best
    ## Run 173 stress 0.08926077 
    ## ... Procrustes: rmse 0.0001818837  max resid 0.0005735416 
    ## ... Similar to previous best
    ## Run 174 stress 0.1071322 
    ## Run 175 stress 0.1071321 
    ## Run 176 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183331  max resid 0.03971077 
    ## Run 177 stress 0.09044627 
    ## Run 178 stress 0.09503437 
    ## Run 179 stress 0.1075422 
    ## Run 180 stress 0.1104182 
    ## Run 181 stress 0.1075422 
    ## Run 182 stress 0.08938966 
    ## ... Procrustes: rmse 0.0359391  max resid 0.1182667 
    ## Run 183 stress 0.08938543 
    ## ... Procrustes: rmse 0.01181156  max resid 0.03972015 
    ## Run 184 stress 0.1104181 
    ## Run 185 stress 0.09503441 
    ## Run 186 stress 0.1067515 
    ## Run 187 stress 0.1056893 
    ## Run 188 stress 0.09503449 
    ## Run 189 stress 0.08951716 
    ## ... Procrustes: rmse 0.03507342  max resid 0.1156514 
    ## Run 190 stress 0.1052649 
    ## Run 191 stress 0.08926081 
    ## ... Procrustes: rmse 0.0002097584  max resid 0.0006834921 
    ## ... Similar to previous best
    ## Run 192 stress 0.08938555 
    ## ... Procrustes: rmse 0.0117753  max resid 0.03974472 
    ## Run 193 stress 0.1075421 
    ## Run 194 stress 0.1087862 
    ## Run 195 stress 0.0894665 
    ## ... Procrustes: rmse 0.03323393  max resid 0.1163695 
    ## Run 196 stress 0.1052648 
    ## Run 197 stress 0.08951716 
    ## ... Procrustes: rmse 0.03507067  max resid 0.1156582 
    ## Run 198 stress 0.1074435 
    ## Run 199 stress 0.09109109 
    ## Run 200 stress 0.0893855 
    ## ... Procrustes: rmse 0.01176224  max resid 0.03962432 
    ## Run 201 stress 0.1067392 
    ## Run 202 stress 0.1101451 
    ## Run 203 stress 0.1075422 
    ## Run 204 stress 0.0892609 
    ## ... Procrustes: rmse 0.0002697972  max resid 0.0008682722 
    ## ... Similar to previous best
    ## Run 205 stress 0.1079019 
    ## Run 206 stress 0.09130087 
    ## Run 207 stress 0.110418 
    ## Run 208 stress 0.1087864 
    ## Run 209 stress 0.08946652 
    ## ... Procrustes: rmse 0.03324558  max resid 0.1164001 
    ## Run 210 stress 0.109213 
    ## Run 211 stress 0.1108638 
    ## Run 212 stress 0.09503433 
    ## Run 213 stress 0.08938544 
    ## ... Procrustes: rmse 0.01185758  max resid 0.03975643 
    ## Run 214 stress 0.08938546 
    ## ... Procrustes: rmse 0.01183305  max resid 0.0397348 
    ## Run 215 stress 0.08926103 
    ## ... Procrustes: rmse 0.0004837125  max resid 0.001433052 
    ## ... Similar to previous best
    ## Run 216 stress 0.1076303 
    ## Run 217 stress 0.09088608 
    ## Run 218 stress 0.08946654 
    ## ... Procrustes: rmse 0.03321082  max resid 0.116329 
    ## Run 219 stress 0.08951722 
    ## ... Procrustes: rmse 0.03505752  max resid 0.1156447 
    ## Run 220 stress 0.09130087 
    ## Run 221 stress 0.1052647 
    ## Run 222 stress 0.09039131 
    ## Run 223 stress 0.1086567 
    ## Run 224 stress 0.1092091 
    ## Run 225 stress 0.08926076 
    ## ... Procrustes: rmse 0.0001692911  max resid 0.0005484776 
    ## ... Similar to previous best
    ## Run 226 stress 0.106839 
    ## Run 227 stress 0.08938967 
    ## ... Procrustes: rmse 0.0359356  max resid 0.1182627 
    ## Run 228 stress 0.08938961 
    ## ... Procrustes: rmse 0.0359717  max resid 0.1183146 
    ## Run 229 stress 0.08926106 
    ## ... Procrustes: rmse 0.0003590029  max resid 0.001116157 
    ## ... Similar to previous best
    ## Run 230 stress 0.08926072 
    ## ... Procrustes: rmse 0.0002747431  max resid 0.0008280343 
    ## ... Similar to previous best
    ## Run 231 stress 0.09503442 
    ## Run 232 stress 0.09180894 
    ## Run 233 stress 0.1065382 
    ## Run 234 stress 0.08926109 
    ## ... Procrustes: rmse 0.000299739  max resid 0.0009806552 
    ## ... Similar to previous best
    ## Run 235 stress 0.09130099 
    ## Run 236 stress 0.08938963 
    ## ... Procrustes: rmse 0.03597948  max resid 0.118331 
    ## Run 237 stress 0.08946329 
    ## ... Procrustes: rmse 0.03753922  max resid 0.1178084 
    ## Run 238 stress 0.1075421 
    ## Run 239 stress 0.1092947 
    ## Run 240 stress 0.1056908 
    ## Run 241 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183877  max resid 0.03971924 
    ## Run 242 stress 0.08938976 
    ## ... Procrustes: rmse 0.0359215  max resid 0.1182388 
    ## Run 243 stress 0.08938545 
    ## ... Procrustes: rmse 0.01179453  max resid 0.03972664 
    ## Run 244 stress 0.0908736 
    ## Run 245 stress 0.09594341 
    ## Run 246 stress 0.09503427 
    ## Run 247 stress 0.1093479 
    ## Run 248 stress 0.1056901 
    ## Run 249 stress 0.08938543 
    ## ... Procrustes: rmse 0.0118471  max resid 0.03974428 
    ## Run 250 stress 0.0950345 
    ## Run 251 stress 0.1074431 
    ## Run 252 stress 0.1071321 
    ## Run 253 stress 0.1071321 
    ## Run 254 stress 0.08938543 
    ## ... Procrustes: rmse 0.01178525  max resid 0.03965522 
    ## Run 255 stress 0.1071322 
    ## Run 256 stress 0.1075422 
    ## Run 257 stress 0.1075421 
    ## Run 258 stress 0.08938539 
    ## ... Procrustes: rmse 0.01181442  max resid 0.03970838 
    ## Run 259 stress 0.08926086 
    ## ... Procrustes: rmse 0.0002469578  max resid 0.0007954819 
    ## ... Similar to previous best
    ## Run 260 stress 0.08926067 
    ## ... Procrustes: rmse 0.0001700496  max resid 0.0004493965 
    ## ... Similar to previous best
    ## Run 261 stress 0.08926073 
    ## ... Procrustes: rmse 0.0001926491  max resid 0.0005356542 
    ## ... Similar to previous best
    ## Run 262 stress 0.1052647 
    ## Run 263 stress 0.09130093 
    ## Run 264 stress 0.1068388 
    ## Run 265 stress 0.1071026 
    ## Run 266 stress 0.08926091 
    ## ... Procrustes: rmse 0.0002777051  max resid 0.0008762371 
    ## ... Similar to previous best
    ## Run 267 stress 0.09503474 
    ## Run 268 stress 0.08951724 
    ## ... Procrustes: rmse 0.03509518  max resid 0.1156868 
    ## Run 269 stress 0.1052649 
    ## Run 270 stress 0.09775555 
    ## Run 271 stress 0.09088623 
    ## Run 272 stress 0.1091883 
    ## Run 273 stress 0.1105345 
    ## Run 274 stress 0.1074431 
    ## Run 275 stress 0.09503453 
    ## Run 276 stress 0.1092213 
    ## Run 277 stress 0.08926121 
    ## ... Procrustes: rmse 0.0003879087  max resid 0.001267218 
    ## ... Similar to previous best
    ## Run 278 stress 0.1052647 
    ## Run 279 stress 0.08946658 
    ## ... Procrustes: rmse 0.03320726  max resid 0.1163259 
    ## Run 280 stress 0.09130088 
    ## Run 281 stress 0.1075421 
    ## Run 282 stress 0.09099576 
    ## Run 283 stress 0.1075421 
    ## Run 284 stress 0.1071321 
    ## Run 285 stress 0.09039129 
    ## Run 286 stress 0.1079019 
    ## Run 287 stress 0.08938538 
    ## ... Procrustes: rmse 0.01181651  max resid 0.03968033 
    ## Run 288 stress 0.105265 
    ## Run 289 stress 0.09592144 
    ## Run 290 stress 0.1092129 
    ## Run 291 stress 0.1074432 
    ## Run 292 stress 0.0893855 
    ## ... Procrustes: rmse 0.0117878  max resid 0.03974878 
    ## Run 293 stress 0.08926065 
    ## ... Procrustes: rmse 0.0001413928  max resid 0.0003300033 
    ## ... Similar to previous best
    ## Run 294 stress 0.0892608 
    ## ... Procrustes: rmse 0.0002008244  max resid 0.0006361294 
    ## ... Similar to previous best
    ## Run 295 stress 0.1104475 
    ## Run 296 stress 0.090886 
    ## Run 297 stress 0.09044608 
    ## Run 298 stress 0.08946665 
    ## ... Procrustes: rmse 0.03319196  max resid 0.1163012 
    ## Run 299 stress 0.08946658 
    ## ... Procrustes: rmse 0.03320164  max resid 0.1163098 
    ## Run 300 stress 0.1086561 
    ## Run 301 stress 0.1067602 
    ## Run 302 stress 0.10569 
    ## Run 303 stress 0.1092096 
    ## Run 304 stress 0.09775604 
    ## Run 305 stress 0.1074434 
    ## Run 306 stress 0.09021166 
    ## Run 307 stress 0.09503417 
    ## Run 308 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183778  max resid 0.03970982 
    ## Run 309 stress 0.1052649 
    ## Run 310 stress 0.08938967 
    ## ... Procrustes: rmse 0.03598559  max resid 0.1183495 
    ## Run 311 stress 0.08946653 
    ## ... Procrustes: rmse 0.03321555  max resid 0.1163414 
    ## Run 312 stress 0.08946338 
    ## ... Procrustes: rmse 0.03749304  max resid 0.117735 
    ## Run 313 stress 0.08951717 
    ## ... Procrustes: rmse 0.03506128  max resid 0.1156457 
    ## Run 314 stress 0.1067893 
    ## Run 315 stress 0.0950344 
    ## Run 316 stress 0.1063194 
    ## Run 317 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 3.226395e-05  max resid 0.0001131122 
    ## ... Similar to previous best
    ## Run 318 stress 0.1052649 
    ## Run 319 stress 0.08938969 
    ## ... Procrustes: rmse 0.03592627  max resid 0.1182498 
    ## Run 320 stress 0.09109108 
    ## Run 321 stress 0.08938543 
    ## ... Procrustes: rmse 0.01187424  max resid 0.03968699 
    ## Run 322 stress 0.09130094 
    ## Run 323 stress 0.08926094 
    ## ... Procrustes: rmse 0.0003205836  max resid 0.001073881 
    ## ... Similar to previous best
    ## Run 324 stress 0.09018708 
    ## Run 325 stress 0.09109108 
    ## Run 326 stress 0.08938561 
    ## ... Procrustes: rmse 0.01176364  max resid 0.03975969 
    ## Run 327 stress 0.09180883 
    ## Run 328 stress 0.09594352 
    ## Run 329 stress 0.08938539 
    ## ... Procrustes: rmse 0.01186558  max resid 0.03974559 
    ## Run 330 stress 0.08946332 
    ## ... Procrustes: rmse 0.03750524  max resid 0.1177573 
    ## Run 331 stress 0.09109111 
    ## Run 332 stress 0.1071321 
    ## Run 333 stress 0.1092092 
    ## Run 334 stress 0.1085132 
    ## Run 335 stress 0.09228483 
    ## Run 336 stress 0.0893854 
    ## ... Procrustes: rmse 0.01185637  max resid 0.03975149 
    ## Run 337 stress 0.08938539 
    ## ... Procrustes: rmse 0.01184723  max resid 0.03973015 
    ## Run 338 stress 0.08926073 
    ## ... Procrustes: rmse 0.0001716282  max resid 0.0005574065 
    ## ... Similar to previous best
    ## Run 339 stress 0.0893898 
    ## ... Procrustes: rmse 0.03591109  max resid 0.1182276 
    ## Run 340 stress 0.08938544 
    ## ... Procrustes: rmse 0.01188102  max resid 0.03968237 
    ## Run 341 stress 0.0918088 
    ## Run 342 stress 0.0893897 
    ## ... Procrustes: rmse 0.03592362  max resid 0.1182457 
    ## Run 343 stress 0.111895 
    ## Run 344 stress 0.08946653 
    ## ... Procrustes: rmse 0.0332403  max resid 0.1163939 
    ## Run 345 stress 0.1074434 
    ## Run 346 stress 0.1067406 
    ## Run 347 stress 0.105265 
    ## Run 348 stress 0.09130085 
    ## Run 349 stress 0.08938964 
    ## ... Procrustes: rmse 0.03597488  max resid 0.1183317 
    ## Run 350 stress 0.1065371 
    ## Run 351 stress 0.09503463 
    ## Run 352 stress 0.1060089 
    ## Run 353 stress 0.09099527 
    ## Run 354 stress 0.1067367 
    ## Run 355 stress 0.08946331 
    ## ... Procrustes: rmse 0.03750484  max resid 0.1177655 
    ## Run 356 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183077  max resid 0.03970992 
    ## Run 357 stress 0.09039137 
    ## Run 358 stress 0.1076303 
    ## Run 359 stress 0.08951722 
    ## ... Procrustes: rmse 0.03504644  max resid 0.1156398 
    ## Run 360 stress 0.09018707 
    ## Run 361 stress 0.1076304 
    ## Run 362 stress 0.09018707 
    ## Run 363 stress 0.08926082 
    ## ... Procrustes: rmse 0.0002499771  max resid 0.0008199991 
    ## ... Similar to previous best
    ## Run 364 stress 0.1087574 
    ## Run 365 stress 0.09180886 
    ## Run 366 stress 0.1063198 
    ## Run 367 stress 0.1092564 
    ## Run 368 stress 0.10892 
    ## Run 369 stress 0.08926105 
    ## ... Procrustes: rmse 0.0003792319  max resid 0.001272941 
    ## ... Similar to previous best
    ## Run 370 stress 0.08938543 
    ## ... Procrustes: rmse 0.01187667  max resid 0.0396933 
    ## Run 371 stress 0.09503417 
    ## Run 372 stress 0.1079017 
    ## Run 373 stress 0.1071321 
    ## Run 374 stress 0.1071322 
    ## Run 375 stress 0.09087341 
    ## Run 376 stress 0.09039129 
    ## Run 377 stress 0.09021166 
    ## Run 378 stress 0.09039136 
    ## Run 379 stress 0.1052647 
    ## Run 380 stress 0.09130091 
    ## Run 381 stress 0.1091235 
    ## Run 382 stress 0.09503415 
    ## Run 383 stress 0.0892611 
    ## ... Procrustes: rmse 0.0003698346  max resid 0.001240425 
    ## ... Similar to previous best
    ## Run 384 stress 0.1080638 
    ## Run 385 stress 0.09099529 
    ## Run 386 stress 0.1060096 
    ## Run 387 stress 0.08926074 
    ## ... Procrustes: rmse 0.000183966  max resid 0.0006179666 
    ## ... Similar to previous best
    ## Run 388 stress 0.1086563 
    ## Run 389 stress 0.09039129 
    ## Run 390 stress 0.08938544 
    ## ... Procrustes: rmse 0.01186711  max resid 0.03977073 
    ## Run 391 stress 0.1052649 
    ## Run 392 stress 0.09180888 
    ## Run 393 stress 0.09503425 
    ## Run 394 stress 0.09044613 
    ## Run 395 stress 0.1086564 
    ## Run 396 stress 0.0950342 
    ## Run 397 stress 0.08951716 
    ## ... Procrustes: rmse 0.03506358  max resid 0.1156384 
    ## Run 398 stress 0.08946666 
    ## ... Procrustes: rmse 0.03318925  max resid 0.1162984 
    ## Run 399 stress 0.1085125 
    ## Run 400 stress 0.1052647 
    ## Run 401 stress 0.106043 
    ## Run 402 stress 0.1067883 
    ## Run 403 stress 0.08938973 
    ## ... Procrustes: rmse 0.03591971  max resid 0.1182362 
    ## Run 404 stress 0.08938539 
    ## ... Procrustes: rmse 0.0118541  max resid 0.03974351 
    ## Run 405 stress 0.1056903 
    ## Run 406 stress 0.08938544 
    ## ... Procrustes: rmse 0.01184019  max resid 0.03976696 
    ## Run 407 stress 0.0893855 
    ## ... Procrustes: rmse 0.01189501  max resid 0.03967178 
    ## Run 408 stress 0.09039138 
    ## Run 409 stress 0.1060436 
    ## Run 410 stress 0.1052653 
    ## Run 411 stress 0.1092092 
    ## Run 412 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595951  max resid 0.118301 
    ## Run 413 stress 0.0922849 
    ## Run 414 stress 0.1074433 
    ## Run 415 stress 0.1087575 
    ## Run 416 stress 0.08946651 
    ## ... Procrustes: rmse 0.03323307  max resid 0.1163701 
    ## Run 417 stress 0.1052647 
    ## Run 418 stress 0.1056901 
    ## Run 419 stress 0.1092129 
    ## Run 420 stress 0.1076305 
    ## Run 421 stress 0.09130094 
    ## Run 422 stress 0.08946662 
    ## ... Procrustes: rmse 0.03319082  max resid 0.1162987 
    ## Run 423 stress 0.1074433 
    ## Run 424 stress 0.09044608 
    ## Run 425 stress 0.09039132 
    ## Run 426 stress 0.08946328 
    ## ... Procrustes: rmse 0.03752729  max resid 0.117785 
    ## Run 427 stress 0.1052646 
    ## Run 428 stress 0.08938551 
    ## ... Procrustes: rmse 0.01186319  max resid 0.0396293 
    ## Run 429 stress 0.08926069 
    ## ... Procrustes: rmse 0.0001158184  max resid 0.000392319 
    ## ... Similar to previous best
    ## Run 430 stress 0.1067586 
    ## Run 431 stress 0.1064189 
    ## Run 432 stress 0.08946665 
    ## ... Procrustes: rmse 0.03318586  max resid 0.1162921 
    ## Run 433 stress 0.09088618 
    ## Run 434 stress 0.09087342 
    ## Run 435 stress 0.1067513 
    ## Run 436 stress 0.09021165 
    ## Run 437 stress 0.1079014 
    ## Run 438 stress 0.106043 
    ## Run 439 stress 0.08926089 
    ## ... Procrustes: rmse 0.0002941692  max resid 0.0009825075 
    ## ... Similar to previous best
    ## Run 440 stress 0.09180898 
    ## Run 441 stress 0.1071321 
    ## Run 442 stress 0.1061297 
    ## Run 443 stress 0.1056897 
    ## Run 444 stress 0.1060094 
    ## Run 445 stress 0.1076304 
    ## Run 446 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183105  max resid 0.03971604 
    ## Run 447 stress 0.1091828 
    ## Run 448 stress 0.09018709 
    ## Run 449 stress 0.1091825 
    ## Run 450 stress 0.08946666 
    ## ... Procrustes: rmse 0.03318689  max resid 0.1162983 
    ## Run 451 stress 0.09503415 
    ## Run 452 stress 0.1079017 
    ## Run 453 stress 0.1056901 
    ## Run 454 stress 0.09088609 
    ## Run 455 stress 0.08946652 
    ## ... Procrustes: rmse 0.03321563  max resid 0.1163375 
    ## Run 456 stress 0.1085131 
    ## Run 457 stress 0.1075421 
    ## Run 458 stress 0.08926092 
    ## ... Procrustes: rmse 0.0003910801  max resid 0.001084866 
    ## ... Similar to previous best
    ## Run 459 stress 0.089261 
    ## ... Procrustes: rmse 0.0003359521  max resid 0.001105437 
    ## ... Similar to previous best
    ## Run 460 stress 0.1089201 
    ## Run 461 stress 0.1067511 
    ## Run 462 stress 0.08926089 
    ## ... Procrustes: rmse 0.0002760491  max resid 0.0009196601 
    ## ... Similar to previous best
    ## Run 463 stress 0.0959214 
    ## Run 464 stress 0.1087578 
    ## Run 465 stress 0.110832 
    ## Run 466 stress 0.0910911 
    ## Run 467 stress 0.1091825 
    ## Run 468 stress 0.09088601 
    ## Run 469 stress 0.1091884 
    ## Run 470 stress 0.1111369 
    ## Run 471 stress 0.1089907 
    ## Run 472 stress 0.09088607 
    ## Run 473 stress 0.1067499 
    ## Run 474 stress 0.09178305 
    ## Run 475 stress 0.1060437 
    ## Run 476 stress 0.09099593 
    ## Run 477 stress 0.1074431 
    ## Run 478 stress 0.08938544 
    ## ... Procrustes: rmse 0.01189389  max resid 0.03970488 
    ## Run 479 stress 0.09099575 
    ## Run 480 stress 0.09099608 
    ## Run 481 stress 0.08946348 
    ## ... Procrustes: rmse 0.03746947  max resid 0.1176981 
    ## Run 482 stress 0.09262353 
    ## Run 483 stress 0.0901871 
    ## Run 484 stress 0.1092362 
    ## Run 485 stress 0.08926073 
    ## ... Procrustes: rmse 0.0001589412  max resid 0.0005292934 
    ## ... Similar to previous best
    ## Run 486 stress 0.1060437 
    ## Run 487 stress 0.08926111 
    ## ... Procrustes: rmse 0.0004906171  max resid 0.001409013 
    ## ... Similar to previous best
    ## Run 488 stress 0.1086562 
    ## Run 489 stress 0.0893855 
    ## ... Procrustes: rmse 0.01185228  max resid 0.03962609 
    ## Run 490 stress 0.08946652 
    ## ... Procrustes: rmse 0.0332367  max resid 0.1163813 
    ## Run 491 stress 0.08926096 
    ## ... Procrustes: rmse 0.0003054081  max resid 0.001017122 
    ## ... Similar to previous best
    ## Run 492 stress 0.1075422 
    ## Run 493 stress 0.08926073 
    ## ... Procrustes: rmse 0.0001653625  max resid 0.0005569083 
    ## ... Similar to previous best
    ## Run 494 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 2.061029e-05  max resid 5.965258e-05 
    ## ... Similar to previous best
    ## Run 495 stress 0.1111371 
    ## Run 496 stress 0.1064193 
    ## Run 497 stress 0.09503425 
    ## Run 498 stress 0.08938549 
    ## ... Procrustes: rmse 0.01178817  max resid 0.03973354 
    ## Run 499 stress 0.1071323 
    ## Run 500 stress 0.08926066 
    ## ... Procrustes: rmse 6.947916e-05  max resid 0.0002172911 
    ## ... Similar to previous best
    ## *** Best solution repeated 2 times

``` r
round(SD_beta_NMDS$stress, digits = 2)
```

    ## [1] 0.09

``` r
SD_beta_rep_NMDS <- metaMDS(SD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.3323292 
    ## Run 1 stress 0.3358504 
    ## Run 2 stress 0.334617 
    ## Run 3 stress 0.3372972 
    ## Run 4 stress 0.3297162 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1583427  max resid 0.3225316 
    ## Run 5 stress 0.339168 
    ## Run 6 stress 0.3398654 
    ## Run 7 stress 0.3434048 
    ## Run 8 stress 0.3301673 
    ## ... Procrustes: rmse 0.1865879  max resid 0.327492 
    ## Run 9 stress 0.3462408 
    ## Run 10 stress 0.3356718 
    ## Run 11 stress 0.3388978 
    ## Run 12 stress 0.3269391 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1886852  max resid 0.3531058 
    ## Run 13 stress 0.3411712 
    ## Run 14 stress 0.333416 
    ## Run 15 stress 0.3345968 
    ## Run 16 stress 0.3319122 
    ## Run 17 stress 0.34015 
    ## Run 18 stress 0.3310525 
    ## Run 19 stress 0.342847 
    ## Run 20 stress 0.3393844 
    ## Run 21 stress 0.3290954 
    ## Run 22 stress 0.3312386 
    ## Run 23 stress 0.3302667 
    ## Run 24 stress 0.3310577 
    ## Run 25 stress 0.3311131 
    ## Run 26 stress 0.3468516 
    ## Run 27 stress 0.3285562 
    ## Run 28 stress 0.3281589 
    ## Run 29 stress 0.3414469 
    ## Run 30 stress 0.3397019 
    ## Run 31 stress 0.3401059 
    ## Run 32 stress 0.3513006 
    ## Run 33 stress 0.3335814 
    ## Run 34 stress 0.3435156 
    ## Run 35 stress 0.3308837 
    ## Run 36 stress 0.3341112 
    ## Run 37 stress 0.3277483 
    ## Run 38 stress 0.332953 
    ## Run 39 stress 0.3379651 
    ## Run 40 stress 0.340508 
    ## Run 41 stress 0.3442455 
    ## Run 42 stress 0.3365778 
    ## Run 43 stress 0.333842 
    ## Run 44 stress 0.3493546 
    ## Run 45 stress 0.3331838 
    ## Run 46 stress 0.3458905 
    ## Run 47 stress 0.3290636 
    ## Run 48 stress 0.3375798 
    ## Run 49 stress 0.3356252 
    ## Run 50 stress 0.3358014 
    ## Run 51 stress 0.3342495 
    ## Run 52 stress 0.3283209 
    ## Run 53 stress 0.3454152 
    ## Run 54 stress 0.3391481 
    ## Run 55 stress 0.3401843 
    ## Run 56 stress 0.3395001 
    ## Run 57 stress 0.3464923 
    ## Run 58 stress 0.3342014 
    ## Run 59 stress 0.3457104 
    ## Run 60 stress 0.3448915 
    ## Run 61 stress 0.3301168 
    ## Run 62 stress 0.3374971 
    ## Run 63 stress 0.3320426 
    ## Run 64 stress 0.3314206 
    ## Run 65 stress 0.3421621 
    ## Run 66 stress 0.341967 
    ## Run 67 stress 0.3310718 
    ## Run 68 stress 0.3377787 
    ## Run 69 stress 0.3298809 
    ## Run 70 stress 0.3264142 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1574647  max resid 0.4061205 
    ## Run 71 stress 0.333308 
    ## Run 72 stress 0.3332965 
    ## Run 73 stress 0.334892 
    ## Run 74 stress 0.3451706 
    ## Run 75 stress 0.3358413 
    ## Run 76 stress 0.3408192 
    ## Run 77 stress 0.3428354 
    ## Run 78 stress 0.3361546 
    ## Run 79 stress 0.3319461 
    ## Run 80 stress 0.3345171 
    ## Run 81 stress 0.3381458 
    ## Run 82 stress 0.330873 
    ## Run 83 stress 0.3428902 
    ## Run 84 stress 0.3305386 
    ## Run 85 stress 0.3807946 
    ## Run 86 stress 0.3392537 
    ## Run 87 stress 0.3312102 
    ## Run 88 stress 0.3305544 
    ## Run 89 stress 0.3350279 
    ## Run 90 stress 0.3337854 
    ## Run 91 stress 0.3330022 
    ## Run 92 stress 0.3474676 
    ## Run 93 stress 0.3288212 
    ## Run 94 stress 0.3250528 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1563775  max resid 0.3821772 
    ## Run 95 stress 0.3310961 
    ## Run 96 stress 0.3436568 
    ## Run 97 stress 0.3391791 
    ## Run 98 stress 0.3337695 
    ## Run 99 stress 0.3270848 
    ## Run 100 stress 0.3390876 
    ## Run 101 stress 0.3491185 
    ## Run 102 stress 0.3325213 
    ## Run 103 stress 0.3392048 
    ## Run 104 stress 0.3470347 
    ## Run 105 stress 0.3300529 
    ## Run 106 stress 0.343172 
    ## Run 107 stress 0.3465162 
    ## Run 108 stress 0.3285386 
    ## Run 109 stress 0.3439424 
    ## Run 110 stress 0.344262 
    ## Run 111 stress 0.3287655 
    ## Run 112 stress 0.3411487 
    ## Run 113 stress 0.3342564 
    ## Run 114 stress 0.3371386 
    ## Run 115 stress 0.3507644 
    ## Run 116 stress 0.3365472 
    ## Run 117 stress 0.3304125 
    ## Run 118 stress 0.3361783 
    ## Run 119 stress 0.3315104 
    ## Run 120 stress 0.3432836 
    ## Run 121 stress 0.3368883 
    ## Run 122 stress 0.3321134 
    ## Run 123 stress 0.3270107 
    ## Run 124 stress 0.346469 
    ## Run 125 stress 0.3279003 
    ## Run 126 stress 0.3298303 
    ## Run 127 stress 0.3298538 
    ## Run 128 stress 0.3362836 
    ## Run 129 stress 0.3406055 
    ## Run 130 stress 0.325933 
    ## Run 131 stress 0.3272677 
    ## Run 132 stress 0.331174 
    ## Run 133 stress 0.3386109 
    ## Run 134 stress 0.328569 
    ## Run 135 stress 0.3324728 
    ## Run 136 stress 0.3427511 
    ## Run 137 stress 0.3348819 
    ## Run 138 stress 0.3309907 
    ## Run 139 stress 0.3389601 
    ## Run 140 stress 0.3479368 
    ## Run 141 stress 0.3360003 
    ## Run 142 stress 0.34894 
    ## Run 143 stress 0.3296388 
    ## Run 144 stress 0.3298098 
    ## Run 145 stress 0.332135 
    ## Run 146 stress 0.3342611 
    ## Run 147 stress 0.3383529 
    ## Run 148 stress 0.3506216 
    ## Run 149 stress 0.3416183 
    ## Run 150 stress 0.3337716 
    ## Run 151 stress 0.3348515 
    ## Run 152 stress 0.3434094 
    ## Run 153 stress 0.3321025 
    ## Run 154 stress 0.3393886 
    ## Run 155 stress 0.3309038 
    ## Run 156 stress 0.345828 
    ## Run 157 stress 0.330778 
    ## Run 158 stress 0.3480134 
    ## Run 159 stress 0.3395556 
    ## Run 160 stress 0.3324158 
    ## Run 161 stress 0.3434432 
    ## Run 162 stress 0.3310425 
    ## Run 163 stress 0.3316008 
    ## Run 164 stress 0.3557017 
    ## Run 165 stress 0.3326938 
    ## Run 166 stress 0.3311382 
    ## Run 167 stress 0.3346165 
    ## Run 168 stress 0.3430385 
    ## Run 169 stress 0.332206 
    ## Run 170 stress 0.3369884 
    ## Run 171 stress 0.333254 
    ## Run 172 stress 0.3406809 
    ## Run 173 stress 0.3326185 
    ## Run 174 stress 0.346851 
    ## Run 175 stress 0.3452951 
    ## Run 176 stress 0.3366304 
    ## Run 177 stress 0.3346523 
    ## Run 178 stress 0.3409457 
    ## Run 179 stress 0.345634 
    ## Run 180 stress 0.3452937 
    ## Run 181 stress 0.336646 
    ## Run 182 stress 0.3351095 
    ## Run 183 stress 0.3424005 
    ## Run 184 stress 0.3421509 
    ## Run 185 stress 0.3340734 
    ## Run 186 stress 0.3305334 
    ## Run 187 stress 0.3326252 
    ## Run 188 stress 0.3329656 
    ## Run 189 stress 0.3317837 
    ## Run 190 stress 0.3433156 
    ## Run 191 stress 0.3393199 
    ## Run 192 stress 0.3410703 
    ## Run 193 stress 0.3325706 
    ## Run 194 stress 0.3367263 
    ## Run 195 stress 0.3313894 
    ## Run 196 stress 0.3406977 
    ## Run 197 stress 0.3309986 
    ## Run 198 stress 0.3358409 
    ## Run 199 stress 0.3366335 
    ## Run 200 stress 0.3321991 
    ## Run 201 stress 0.3439732 
    ## Run 202 stress 0.3333549 
    ## Run 203 stress 0.3373625 
    ## Run 204 stress 0.3363764 
    ## Run 205 stress 0.3308381 
    ## Run 206 stress 0.3346246 
    ## Run 207 stress 0.3336925 
    ## Run 208 stress 0.3431478 
    ## Run 209 stress 0.3271571 
    ## Run 210 stress 0.3528613 
    ## Run 211 stress 0.3327477 
    ## Run 212 stress 0.3332444 
    ## Run 213 stress 0.3368167 
    ## Run 214 stress 0.3346299 
    ## Run 215 stress 0.3391297 
    ## Run 216 stress 0.3284 
    ## Run 217 stress 0.3308337 
    ## Run 218 stress 0.3342033 
    ## Run 219 stress 0.3293337 
    ## Run 220 stress 0.3339814 
    ## Run 221 stress 0.332494 
    ## Run 222 stress 0.3291793 
    ## Run 223 stress 0.3480134 
    ## Run 224 stress 0.3330526 
    ## Run 225 stress 0.3419669 
    ## Run 226 stress 0.3433727 
    ## Run 227 stress 0.329409 
    ## Run 228 stress 0.3358357 
    ## Run 229 stress 0.3475088 
    ## Run 230 stress 0.3377128 
    ## Run 231 stress 0.3300873 
    ## Run 232 stress 0.3467259 
    ## Run 233 stress 0.3338645 
    ## Run 234 stress 0.3427193 
    ## Run 235 stress 0.3479322 
    ## Run 236 stress 0.3432571 
    ## Run 237 stress 0.3348671 
    ## Run 238 stress 0.3345102 
    ## Run 239 stress 0.3483238 
    ## Run 240 stress 0.3372378 
    ## Run 241 stress 0.3403856 
    ## Run 242 stress 0.3434974 
    ## Run 243 stress 0.3330391 
    ## Run 244 stress 0.3365704 
    ## Run 245 stress 0.3313614 
    ## Run 246 stress 0.343685 
    ## Run 247 stress 0.3434895 
    ## Run 248 stress 0.3333507 
    ## Run 249 stress 0.3416026 
    ## Run 250 stress 0.3349724 
    ## Run 251 stress 0.3462695 
    ## Run 252 stress 0.3342658 
    ## Run 253 stress 0.3487501 
    ## Run 254 stress 0.3497071 
    ## Run 255 stress 0.3420315 
    ## Run 256 stress 0.3275701 
    ## Run 257 stress 0.3502014 
    ## Run 258 stress 0.3361576 
    ## Run 259 stress 0.338336 
    ## Run 260 stress 0.3304614 
    ## Run 261 stress 0.3299169 
    ## Run 262 stress 0.3332459 
    ## Run 263 stress 0.3430494 
    ## Run 264 stress 0.336565 
    ## Run 265 stress 0.3335466 
    ## Run 266 stress 0.3412228 
    ## Run 267 stress 0.3303562 
    ## Run 268 stress 0.3327043 
    ## Run 269 stress 0.3323947 
    ## Run 270 stress 0.340907 
    ## Run 271 stress 0.336685 
    ## Run 272 stress 0.3426605 
    ## Run 273 stress 0.351669 
    ## Run 274 stress 0.3402737 
    ## Run 275 stress 0.3322247 
    ## Run 276 stress 0.3353011 
    ## Run 277 stress 0.3454144 
    ## Run 278 stress 0.3300076 
    ## Run 279 stress 0.3290871 
    ## Run 280 stress 0.3389328 
    ## Run 281 stress 0.348539 
    ## Run 282 stress 0.3347496 
    ## Run 283 stress 0.3411516 
    ## Run 284 stress 0.334578 
    ## Run 285 stress 0.3319639 
    ## Run 286 stress 0.3315082 
    ## Run 287 stress 0.3299424 
    ## Run 288 stress 0.3298275 
    ## Run 289 stress 0.3289783 
    ## Run 290 stress 0.3338068 
    ## Run 291 stress 0.3384556 
    ## Run 292 stress 0.3298803 
    ## Run 293 stress 0.3437866 
    ## Run 294 stress 0.3335977 
    ## Run 295 stress 0.3302526 
    ## Run 296 stress 0.340583 
    ## Run 297 stress 0.3349667 
    ## Run 298 stress 0.3297059 
    ## Run 299 stress 0.3308752 
    ## Run 300 stress 0.3365177 
    ## Run 301 stress 0.3371792 
    ## Run 302 stress 0.3294159 
    ## Run 303 stress 0.3360299 
    ## Run 304 stress 0.3355156 
    ## Run 305 stress 0.3357494 
    ## Run 306 stress 0.330751 
    ## Run 307 stress 0.336442 
    ## Run 308 stress 0.3320805 
    ## Run 309 stress 0.3418644 
    ## Run 310 stress 0.3313753 
    ## Run 311 stress 0.3332189 
    ## Run 312 stress 0.3507199 
    ## Run 313 stress 0.3396736 
    ## Run 314 stress 0.3324478 
    ## Run 315 stress 0.3389701 
    ## Run 316 stress 0.3316881 
    ## Run 317 stress 0.3407663 
    ## Run 318 stress 0.3362008 
    ## Run 319 stress 0.3359864 
    ## Run 320 stress 0.3318653 
    ## Run 321 stress 0.3300982 
    ## Run 322 stress 0.3410184 
    ## Run 323 stress 0.3341829 
    ## Run 324 stress 0.3412052 
    ## Run 325 stress 0.3734771 
    ## Run 326 stress 0.3386645 
    ## Run 327 stress 0.3329592 
    ## Run 328 stress 0.3287456 
    ## Run 329 stress 0.3323984 
    ## Run 330 stress 0.3321325 
    ## Run 331 stress 0.3349179 
    ## Run 332 stress 0.3332329 
    ## Run 333 stress 0.3322616 
    ## Run 334 stress 0.3390547 
    ## Run 335 stress 0.3318828 
    ## Run 336 stress 0.3391891 
    ## Run 337 stress 0.3370543 
    ## Run 338 stress 0.3308032 
    ## Run 339 stress 0.3374225 
    ## Run 340 stress 0.3387826 
    ## Run 341 stress 0.3335314 
    ## Run 342 stress 0.3355799 
    ## Run 343 stress 0.3328538 
    ## Run 344 stress 0.3348115 
    ## Run 345 stress 0.3447847 
    ## Run 346 stress 0.3396927 
    ## Run 347 stress 0.3323351 
    ## Run 348 stress 0.326519 
    ## Run 349 stress 0.3297562 
    ## Run 350 stress 0.3312818 
    ## Run 351 stress 0.3470229 
    ## Run 352 stress 0.3407229 
    ## Run 353 stress 0.3270871 
    ## Run 354 stress 0.3444493 
    ## Run 355 stress 0.3335096 
    ## Run 356 stress 0.3244065 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1511168  max resid 0.3823134 
    ## Run 357 stress 0.3344932 
    ## Run 358 stress 0.3412998 
    ## Run 359 stress 0.3388882 
    ## Run 360 stress 0.3368526 
    ## Run 361 stress 0.3466823 
    ## Run 362 stress 0.3372717 
    ## Run 363 stress 0.3302815 
    ## Run 364 stress 0.3343458 
    ## Run 365 stress 0.3480692 
    ## Run 366 stress 0.3410526 
    ## Run 367 stress 0.3346609 
    ## Run 368 stress 0.3418101 
    ## Run 369 stress 0.3412869 
    ## Run 370 stress 0.3358158 
    ## Run 371 stress 0.3444249 
    ## Run 372 stress 0.3380906 
    ## Run 373 stress 0.3281479 
    ## Run 374 stress 0.3503097 
    ## Run 375 stress 0.3422047 
    ## Run 376 stress 0.3334126 
    ## Run 377 stress 0.3320979 
    ## Run 378 stress 0.3355468 
    ## Run 379 stress 0.3331058 
    ## Run 380 stress 0.3386694 
    ## Run 381 stress 0.3428824 
    ## Run 382 stress 0.3441593 
    ## Run 383 stress 0.3286202 
    ## Run 384 stress 0.3299238 
    ## Run 385 stress 0.3479175 
    ## Run 386 stress 0.3457708 
    ## Run 387 stress 0.3275729 
    ## Run 388 stress 0.340352 
    ## Run 389 stress 0.329005 
    ## Run 390 stress 0.3269794 
    ## Run 391 stress 0.329873 
    ## Run 392 stress 0.3393307 
    ## Run 393 stress 0.3304198 
    ## Run 394 stress 0.3264628 
    ## Run 395 stress 0.3410569 
    ## Run 396 stress 0.3414719 
    ## Run 397 stress 0.3417419 
    ## Run 398 stress 0.3454756 
    ## Run 399 stress 0.3322601 
    ## Run 400 stress 0.3397785 
    ## Run 401 stress 0.3358288 
    ## Run 402 stress 0.3437615 
    ## Run 403 stress 0.334615 
    ## Run 404 stress 0.3455822 
    ## Run 405 stress 0.3342834 
    ## Run 406 stress 0.3345191 
    ## Run 407 stress 0.3327682 
    ## Run 408 stress 0.3307486 
    ## Run 409 stress 0.3437864 
    ## Run 410 stress 0.3291193 
    ## Run 411 stress 0.3329984 
    ## Run 412 stress 0.3381624 
    ## Run 413 stress 0.3356707 
    ## Run 414 stress 0.3450129 
    ## Run 415 stress 0.3363885 
    ## Run 416 stress 0.331 
    ## Run 417 stress 0.3417594 
    ## Run 418 stress 0.3312903 
    ## Run 419 stress 0.3376557 
    ## Run 420 stress 0.3363467 
    ## Run 421 stress 0.349862 
    ## Run 422 stress 0.3344471 
    ## Run 423 stress 0.3330992 
    ## Run 424 stress 0.3453897 
    ## Run 425 stress 0.343171 
    ## Run 426 stress 0.3428825 
    ## Run 427 stress 0.3479538 
    ## Run 428 stress 0.334774 
    ## Run 429 stress 0.3285926 
    ## Run 430 stress 0.3366752 
    ## Run 431 stress 0.3340261 
    ## Run 432 stress 0.350242 
    ## Run 433 stress 0.3336259 
    ## Run 434 stress 0.3293173 
    ## Run 435 stress 0.3330004 
    ## Run 436 stress 0.3505966 
    ## Run 437 stress 0.3505096 
    ## Run 438 stress 0.3422803 
    ## Run 439 stress 0.3361661 
    ## Run 440 stress 0.3376359 
    ## Run 441 stress 0.3375998 
    ## Run 442 stress 0.3402389 
    ## Run 443 stress 0.3378554 
    ## Run 444 stress 0.3302534 
    ## Run 445 stress 0.3421584 
    ## Run 446 stress 0.3424331 
    ## Run 447 stress 0.3366894 
    ## Run 448 stress 0.3433502 
    ## Run 449 stress 0.3502531 
    ## Run 450 stress 0.3320904 
    ## Run 451 stress 0.3441349 
    ## Run 452 stress 0.3282055 
    ## Run 453 stress 0.3381721 
    ## Run 454 stress 0.3483069 
    ## Run 455 stress 0.3506088 
    ## Run 456 stress 0.330685 
    ## Run 457 stress 0.3300425 
    ## Run 458 stress 0.3357012 
    ## Run 459 stress 0.339142 
    ## Run 460 stress 0.3349757 
    ## Run 461 stress 0.3418123 
    ## Run 462 stress 0.3423049 
    ## Run 463 stress 0.336879 
    ## Run 464 stress 0.3459171 
    ## Run 465 stress 0.3318513 
    ## Run 466 stress 0.3308013 
    ## Run 467 stress 0.34758 
    ## Run 468 stress 0.3356644 
    ## Run 469 stress 0.3324663 
    ## Run 470 stress 0.3364948 
    ## Run 471 stress 0.3296851 
    ## Run 472 stress 0.3360881 
    ## Run 473 stress 0.3350501 
    ## Run 474 stress 0.3354068 
    ## Run 475 stress 0.3423356 
    ## Run 476 stress 0.3304567 
    ## Run 477 stress 0.3373464 
    ## Run 478 stress 0.337875 
    ## Run 479 stress 0.3308343 
    ## Run 480 stress 0.3327952 
    ## Run 481 stress 0.3438055 
    ## Run 482 stress 0.3451825 
    ## Run 483 stress 0.3442202 
    ## Run 484 stress 0.3371389 
    ## Run 485 stress 0.3380025 
    ## Run 486 stress 0.3421212 
    ## Run 487 stress 0.3341492 
    ## Run 488 stress 0.3392754 
    ## Run 489 stress 0.3412558 
    ## Run 490 stress 0.338277 
    ## Run 491 stress 0.3362402 
    ## Run 492 stress 0.3305353 
    ## Run 493 stress 0.3277674 
    ## Run 494 stress 0.3294955 
    ## Run 495 stress 0.3343785 
    ## Run 496 stress 0.3361428 
    ## Run 497 stress 0.3316326 
    ## Run 498 stress 0.332034 
    ## Run 499 stress 0.3289339 
    ## Run 500 stress 0.3383598 
    ## *** Best solution was not repeated -- monoMDS stopping criteria:
    ##      1: no. of iterations >= maxit
    ##    499: stress ratio > sratmax

``` r
SD_beta_ric_NMDS <- metaMDS(SD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.01879552 
    ## Run 1 stress 0.02509475 
    ## Run 2 stress 0.01883689 
    ## ... Procrustes: rmse 0.004972737  max resid 0.01023381 
    ## Run 3 stress 0.02486131 
    ## Run 4 stress 0.01884275 
    ## ... Procrustes: rmse 0.004452013  max resid 0.008825664 
    ## ... Similar to previous best
    ## Run 5 stress 0.02492212 
    ## Run 6 stress 0.01880774 
    ## ... Procrustes: rmse 0.003332019  max resid 0.006859562 
    ## ... Similar to previous best
    ## Run 7 stress 0.02509473 
    ## Run 8 stress 0.01879596 
    ## ... Procrustes: rmse 0.0001988156  max resid 0.0004080858 
    ## ... Similar to previous best
    ## Run 9 stress 0.018796 
    ## ... Procrustes: rmse 0.001326854  max resid 0.002732491 
    ## ... Similar to previous best
    ## Run 10 stress 0.02492191 
    ## Run 11 stress 0.01880301 
    ## ... Procrustes: rmse 0.002549602  max resid 0.005247403 
    ## ... Similar to previous best
    ## Run 12 stress 0.01879571 
    ## ... Procrustes: rmse 0.001215344  max resid 0.002503009 
    ## ... Similar to previous best
    ## Run 13 stress 0.02492181 
    ## Run 14 stress 0.01879557 
    ## ... Procrustes: rmse 0.001159556  max resid 0.002388218 
    ## ... Similar to previous best
    ## Run 15 stress 0.01879555 
    ## ... Procrustes: rmse 0.001147448  max resid 0.002363295 
    ## ... Similar to previous best
    ## Run 16 stress 0.02509436 
    ## Run 17 stress 0.01879584 
    ## ... Procrustes: rmse 0.001283161  max resid 0.002642586 
    ## ... Similar to previous best
    ## Run 18 stress 0.01880106 
    ## ... Procrustes: rmse 0.00238916  max resid 0.004918456 
    ## ... Similar to previous best
    ## Run 19 stress 0.01881072 
    ## ... Procrustes: rmse 0.003591027  max resid 0.007390888 
    ## ... Similar to previous best
    ## Run 20 stress 0.01888645 
    ## ... Procrustes: rmse 0.007299814  max resid 0.0147733 
    ## Run 21 stress 0.02520883 
    ## Run 22 stress 0.02471498 
    ## Run 23 stress 0.01883616 
    ## ... Procrustes: rmse 0.00373835  max resid 0.007315802 
    ## ... Similar to previous best
    ## Run 24 stress 0.02509475 
    ## Run 25 stress 0.01879561 
    ## ... Procrustes: rmse 4.281735e-05  max resid 8.750467e-05 
    ## ... Similar to previous best
    ## Run 26 stress 0.01879588 
    ## ... Procrustes: rmse 0.001286942  max resid 0.002650473 
    ## ... Similar to previous best
    ## Run 27 stress 0.02509442 
    ## Run 28 stress 0.02492187 
    ## Run 29 stress 0.01879568 
    ## ... Procrustes: rmse 7.792599e-05  max resid 0.0001599019 
    ## ... Similar to previous best
    ## Run 30 stress 0.01883595 
    ## ... Procrustes: rmse 0.003720022  max resid 0.007279428 
    ## ... Similar to previous best
    ## Run 31 stress 0.02509442 
    ## Run 32 stress 0.02520903 
    ## Run 33 stress 0.0187955 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001124371  max resid 0.002315796 
    ## ... Similar to previous best
    ## Run 34 stress 0.02520921 
    ## Run 35 stress 0.01879572 
    ## ... Procrustes: rmse 0.00122371  max resid 0.002522374 
    ## ... Similar to previous best
    ## Run 36 stress 0.01879555 
    ## ... Procrustes: rmse 2.281094e-05  max resid 4.681935e-05 
    ## ... Similar to previous best
    ## Run 37 stress 0.02509455 
    ## Run 38 stress 0.02486118 
    ## Run 39 stress 0.01879556 
    ## ... Procrustes: rmse 2.934814e-05  max resid 6.046031e-05 
    ## ... Similar to previous best
    ## Run 40 stress 0.01930224 
    ## Run 41 stress 0.01879808 
    ## ... Procrustes: rmse 0.000785562  max resid 0.001619287 
    ## ... Similar to previous best
    ## Run 42 stress 0.01879565 
    ## ... Procrustes: rmse 6.967987e-05  max resid 0.0001435133 
    ## ... Similar to previous best
    ## Run 43 stress 0.01879588 
    ## ... Procrustes: rmse 0.0001738076  max resid 0.0003580235 
    ## ... Similar to previous best
    ## Run 44 stress 0.0250946 
    ## Run 45 stress 0.02486128 
    ## Run 46 stress 0.02492324 
    ## Run 47 stress 0.0250944 
    ## Run 48 stress 0.02492175 
    ## Run 49 stress 0.01883479 
    ## ... Procrustes: rmse 0.004363984  max resid 0.008980075 
    ## ... Similar to previous best
    ## Run 50 stress 0.01879571 
    ## ... Procrustes: rmse 0.0001031773  max resid 0.0002125785 
    ## ... Similar to previous best
    ## Run 51 stress 0.02486133 
    ## Run 52 stress 0.0250946 
    ## Run 53 stress 0.01879559 
    ## ... Procrustes: rmse 0.001156866  max resid 0.002385039 
    ## ... Similar to previous best
    ## Run 54 stress 0.01879528 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001367705  max resid 0.000281823 
    ## ... Similar to previous best
    ## Run 55 stress 0.0187954 
    ## ... Procrustes: rmse 7.895217e-05  max resid 0.0001625905 
    ## ... Similar to previous best
    ## Run 56 stress 0.01880468 
    ## ... Procrustes: rmse 0.001963053  max resid 0.004045885 
    ## ... Similar to previous best
    ## Run 57 stress 0.02509473 
    ## Run 58 stress 0.02486148 
    ## Run 59 stress 0.01879574 
    ## ... Procrustes: rmse 0.0002527985  max resid 0.0005207856 
    ## ... Similar to previous best
    ## Run 60 stress 0.01879568 
    ## ... Procrustes: rmse 0.0002243902  max resid 0.0004623373 
    ## ... Similar to previous best
    ## Run 61 stress 0.0188122 
    ## ... Procrustes: rmse 0.002760003  max resid 0.005685176 
    ## ... Similar to previous best
    ## Run 62 stress 0.01883671 
    ## ... Procrustes: rmse 0.004557649  max resid 0.009379295 
    ## ... Similar to previous best
    ## Run 63 stress 0.0188679 
    ## ... Procrustes: rmse 0.00532311  max resid 0.01064472 
    ## Run 64 stress 0.02486108 
    ## Run 65 stress 0.01879563 
    ## ... Procrustes: rmse 0.0001963815  max resid 0.000404379 
    ## ... Similar to previous best
    ## Run 66 stress 0.01879577 
    ## ... Procrustes: rmse 0.0002645695  max resid 0.0005450793 
    ## ... Similar to previous best
    ## Run 67 stress 0.01892974 
    ## ... Procrustes: rmse 0.00869453  max resid 0.01785245 
    ## Run 68 stress 0.01879943 
    ## ... Procrustes: rmse 0.0005609279  max resid 0.001149077 
    ## ... Similar to previous best
    ## Run 69 stress 0.01887827 
    ## ... Procrustes: rmse 0.006718812  max resid 0.01380898 
    ## Run 70 stress 0.02486112 
    ## Run 71 stress 0.01879579 
    ## ... Procrustes: rmse 0.0002619801  max resid 0.0005393537 
    ## ... Similar to previous best
    ## Run 72 stress 0.2710607 
    ## Run 73 stress 0.02486115 
    ## Run 74 stress 0.3690982 
    ## Run 75 stress 0.01881765 
    ## ... Procrustes: rmse 0.003009742  max resid 0.006195974 
    ## ... Similar to previous best
    ## Run 76 stress 0.01879577 
    ## ... Procrustes: rmse 0.0002631681  max resid 0.0005421653 
    ## ... Similar to previous best
    ## Run 77 stress 0.02492152 
    ## Run 78 stress 0.01879542 
    ## ... Procrustes: rmse 0.0009352572  max resid 0.001928641 
    ## ... Similar to previous best
    ## Run 79 stress 0.02493315 
    ## Run 80 stress 0.01879613 
    ## ... Procrustes: rmse 0.0004046839  max resid 0.0008335911 
    ## ... Similar to previous best
    ## Run 81 stress 0.02509451 
    ## Run 82 stress 0.02492197 
    ## Run 83 stress 0.01880161 
    ## ... Procrustes: rmse 0.001569541  max resid 0.003234263 
    ## ... Similar to previous best
    ## Run 84 stress 0.02492203 
    ## Run 85 stress 0.02509478 
    ## Run 86 stress 0.01879584 
    ## ... Procrustes: rmse 0.0002785214  max resid 0.0005737867 
    ## ... Similar to previous best
    ## Run 87 stress 0.02486109 
    ## Run 88 stress 0.01881103 
    ## ... Procrustes: rmse 0.002693708  max resid 0.00554922 
    ## ... Similar to previous best
    ## Run 89 stress 0.01884288 
    ## ... Procrustes: rmse 0.003535053  max resid 0.006857333 
    ## ... Similar to previous best
    ## Run 90 stress 0.01883876 
    ## ... Procrustes: rmse 0.003079937  max resid 0.005872284 
    ## ... Similar to previous best
    ## Run 91 stress 0.01879578 
    ## ... Procrustes: rmse 0.001100554  max resid 0.002268167 
    ## ... Similar to previous best
    ## Run 92 stress 0.01885964 
    ## ... Procrustes: rmse 0.004794022  max resid 0.009534548 
    ## ... Similar to previous best
    ## Run 93 stress 0.01879601 
    ## ... Procrustes: rmse 0.0003469607  max resid 0.0007143309 
    ## ... Similar to previous best
    ## Run 94 stress 0.01882993 
    ## ... Procrustes: rmse 0.002049011  max resid 0.004770628 
    ## ... Similar to previous best
    ## Run 95 stress 0.01883607 
    ## ... Procrustes: rmse 0.004553674  max resid 0.009366955 
    ## ... Similar to previous best
    ## Run 96 stress 0.0250943 
    ## Run 97 stress 0.02486113 
    ## Run 98 stress 0.02509464 
    ## Run 99 stress 0.02509492 
    ## Run 100 stress 0.01879539 
    ## ... Procrustes: rmse 7.024478e-05  max resid 0.0001445252 
    ## ... Similar to previous best
    ## Run 101 stress 0.0250947 
    ## Run 102 stress 0.0187958 
    ## ... Procrustes: rmse 0.001118163  max resid 0.002304429 
    ## ... Similar to previous best
    ## Run 103 stress 0.01887385 
    ## ... Procrustes: rmse 0.005654392  max resid 0.0113404 
    ## Run 104 stress 0.02486151 
    ## Run 105 stress 0.01879586 
    ## ... Procrustes: rmse 0.0002879635  max resid 0.0005928785 
    ## ... Similar to previous best
    ## Run 106 stress 0.01879597 
    ## ... Procrustes: rmse 0.0003335294  max resid 0.0006867482 
    ## ... Similar to previous best
    ## Run 107 stress 0.01879579 
    ## ... Procrustes: rmse 0.0002740529  max resid 0.0005646113 
    ## ... Similar to previous best
    ## Run 108 stress 0.02509434 
    ## Run 109 stress 0.02509446 
    ## Run 110 stress 0.01884717 
    ## ... Procrustes: rmse 0.00299633  max resid 0.005694735 
    ## ... Similar to previous best
    ## Run 111 stress 0.02492202 
    ## Run 112 stress 0.01928682 
    ## ... Procrustes: rmse 0.01720581  max resid 0.03542107 
    ## Run 113 stress 0.01879565 
    ## ... Procrustes: rmse 0.0002033591  max resid 0.0004190123 
    ## ... Similar to previous best
    ## Run 114 stress 0.02509466 
    ## Run 115 stress 0.02520901 
    ## Run 116 stress 0.01879581 
    ## ... Procrustes: rmse 0.000268091  max resid 0.0005518714 
    ## ... Similar to previous best
    ## Run 117 stress 0.01880166 
    ## ... Procrustes: rmse 0.001406538  max resid 0.002863173 
    ## ... Similar to previous best
    ## Run 118 stress 0.02520895 
    ## Run 119 stress 0.0188186 
    ## ... Procrustes: rmse 0.003276045  max resid 0.006743754 
    ## ... Similar to previous best
    ## Run 120 stress 0.02025749 
    ## Run 121 stress 0.02486111 
    ## Run 122 stress 0.01879568 
    ## ... Procrustes: rmse 0.0002225009  max resid 0.0004583201 
    ## ... Similar to previous best
    ## Run 123 stress 0.02486122 
    ## Run 124 stress 0.01895859 
    ## ... Procrustes: rmse 0.009611121  max resid 0.01972538 
    ## Run 125 stress 0.02509456 
    ## Run 126 stress 0.01879594 
    ## ... Procrustes: rmse 0.001170965  max resid 0.002412419 
    ## ... Similar to previous best
    ## Run 127 stress 0.01888648 
    ## ... Procrustes: rmse 0.006985191  max resid 0.01435965 
    ## Run 128 stress 0.02509452 
    ## Run 129 stress 0.02509459 
    ## Run 130 stress 0.02509473 
    ## Run 131 stress 0.02509478 
    ## Run 132 stress 0.02492202 
    ## Run 133 stress 0.01939748 
    ## Run 134 stress 0.02492188 
    ## Run 135 stress 0.02486134 
    ## Run 136 stress 0.01879564 
    ## ... Procrustes: rmse 0.001037491  max resid 0.002138681 
    ## ... Similar to previous best
    ## Run 137 stress 0.01882993 
    ## ... Procrustes: rmse 0.002040483  max resid 0.00474549 
    ## ... Similar to previous best
    ## Run 138 stress 0.01886313 
    ## ... Procrustes: rmse 0.005031038  max resid 0.01003151 
    ## Run 139 stress 0.02014456 
    ## Run 140 stress 0.01883437 
    ## ... Procrustes: rmse 0.002655434  max resid 0.005335801 
    ## ... Similar to previous best
    ## Run 141 stress 0.02509441 
    ## Run 142 stress 0.01880304 
    ## ... Procrustes: rmse 0.001780636  max resid 0.003669427 
    ## ... Similar to previous best
    ## Run 143 stress 0.02509439 
    ## Run 144 stress 0.01879581 
    ## ... Procrustes: rmse 0.0002796518  max resid 0.0005761837 
    ## ... Similar to previous best
    ## Run 145 stress 0.02486126 
    ## Run 146 stress 0.0250945 
    ## Run 147 stress 0.02492165 
    ## Run 148 stress 0.02492167 
    ## Run 149 stress 0.018796 
    ## ... Procrustes: rmse 0.0003604818  max resid 0.0007425979 
    ## ... Similar to previous best
    ## Run 150 stress 0.02509448 
    ## Run 151 stress 0.0248613 
    ## Run 152 stress 0.0249219 
    ## Run 153 stress 0.01883332 
    ## ... Procrustes: rmse 0.00439616  max resid 0.009046033 
    ## ... Similar to previous best
    ## Run 154 stress 0.01879562 
    ## ... Procrustes: rmse 0.0001936693  max resid 0.000399036 
    ## ... Similar to previous best
    ## Run 155 stress 0.01879564 
    ## ... Procrustes: rmse 0.0001998502  max resid 0.0004114791 
    ## ... Similar to previous best
    ## Run 156 stress 0.01883838 
    ## ... Procrustes: rmse 0.002725964  max resid 0.00537445 
    ## ... Similar to previous best
    ## Run 157 stress 0.02509428 
    ## Run 158 stress 0.01880195 
    ## ... Procrustes: rmse 0.001592857  max resid 0.003283066 
    ## ... Similar to previous best
    ## Run 159 stress 0.02492146 
    ## Run 160 stress 0.01879581 
    ## ... Procrustes: rmse 0.0002832731  max resid 0.0005835982 
    ## ... Similar to previous best
    ## Run 161 stress 0.02509485 
    ## Run 162 stress 0.02509421 
    ## Run 163 stress 0.0187957 
    ## ... Procrustes: rmse 0.001073719  max resid 0.002213091 
    ## ... Similar to previous best
    ## Run 164 stress 0.01879572 
    ## ... Procrustes: rmse 0.0002414414  max resid 0.0004974166 
    ## ... Similar to previous best
    ## Run 165 stress 0.02509465 
    ## Run 166 stress 0.02509476 
    ## Run 167 stress 0.01946828 
    ## Run 168 stress 0.02492175 
    ## Run 169 stress 0.02509441 
    ## Run 170 stress 0.02492489 
    ## Run 171 stress 0.01879549 
    ## ... Procrustes: rmse 0.0001294531  max resid 0.0002667326 
    ## ... Similar to previous best
    ## Run 172 stress 0.01879577 
    ## ... Procrustes: rmse 0.001107445  max resid 0.002282396 
    ## ... Similar to previous best
    ## Run 173 stress 0.02509448 
    ## Run 174 stress 0.02492173 
    ## Run 175 stress 0.02492221 
    ## Run 176 stress 0.01882989 
    ## ... Procrustes: rmse 0.002036923  max resid 0.004737384 
    ## ... Similar to previous best
    ## Run 177 stress 0.02509463 
    ## Run 178 stress 0.01879583 
    ## ... Procrustes: rmse 0.0002918061  max resid 0.0006011697 
    ## ... Similar to previous best
    ## Run 179 stress 0.02520886 
    ## Run 180 stress 0.01879591 
    ## ... Procrustes: rmse 0.001154665  max resid 0.002379417 
    ## ... Similar to previous best
    ## Run 181 stress 0.01879586 
    ## ... Procrustes: rmse 0.0002984957  max resid 0.0006147555 
    ## ... Similar to previous best
    ## Run 182 stress 0.018834 
    ## ... Procrustes: rmse 0.002608517  max resid 0.005299822 
    ## ... Similar to previous best
    ## Run 183 stress 0.018796 
    ## ... Procrustes: rmse 0.001186675  max resid 0.002445338 
    ## ... Similar to previous best
    ## Run 184 stress 0.0250947 
    ## Run 185 stress 0.01884247 
    ## ... Procrustes: rmse 0.003489221  max resid 0.006758499 
    ## ... Similar to previous best
    ## Run 186 stress 0.01884611 
    ## ... Procrustes: rmse 0.003803037  max resid 0.007432298 
    ## ... Similar to previous best
    ## Run 187 stress 0.01885167 
    ## ... Procrustes: rmse 0.004123226  max resid 0.008118326 
    ## ... Similar to previous best
    ## Run 188 stress 0.0252088 
    ## Run 189 stress 0.02509472 
    ## Run 190 stress 0.02509473 
    ## Run 191 stress 0.0187959 
    ## ... Procrustes: rmse 0.0003147864  max resid 0.0006482772 
    ## ... Similar to previous best
    ## Run 192 stress 0.02492196 
    ## Run 193 stress 0.2039098 
    ## Run 194 stress 0.01887735 
    ## ... Procrustes: rmse 0.005814869  max resid 0.01167272 
    ## Run 195 stress 0.02486139 
    ## Run 196 stress 0.01885412 
    ## ... Procrustes: rmse 0.005588944  max resid 0.01149363 
    ## Run 197 stress 0.02486145 
    ## Run 198 stress 0.02000047 
    ## Run 199 stress 0.02492158 
    ## Run 200 stress 0.0250947 
    ## Run 201 stress 0.01879566 
    ## ... Procrustes: rmse 0.000214655  max resid 0.0004421986 
    ## ... Similar to previous best
    ## Run 202 stress 0.01883144 
    ## ... Procrustes: rmse 0.002274576  max resid 0.005036464 
    ## ... Similar to previous best
    ## Run 203 stress 0.02509446 
    ## Run 204 stress 0.01879563 
    ## ... Procrustes: rmse 0.0001985871  max resid 0.000409116 
    ## ... Similar to previous best
    ## Run 205 stress 0.01880348 
    ## ... Procrustes: rmse 0.001841569  max resid 0.003795099 
    ## ... Similar to previous best
    ## Run 206 stress 0.01879599 
    ## ... Procrustes: rmse 0.0003396687  max resid 0.0006997745 
    ## ... Similar to previous best
    ## Run 207 stress 0.02492185 
    ## Run 208 stress 0.01881131 
    ## ... Procrustes: rmse 0.002714007  max resid 0.005591373 
    ## ... Similar to previous best
    ## Run 209 stress 0.01879929 
    ## ... Procrustes: rmse 0.00118604  max resid 0.002443954 
    ## ... Similar to previous best
    ## Run 210 stress 0.01879553 
    ## ... Procrustes: rmse 0.00014936  max resid 0.0003077328 
    ## ... Similar to previous best
    ## Run 211 stress 0.01879552 
    ## ... Procrustes: rmse 0.0001430538  max resid 0.0002946681 
    ## ... Similar to previous best
    ## Run 212 stress 0.02486136 
    ## Run 213 stress 0.01879559 
    ## ... Procrustes: rmse 0.0001809486  max resid 0.0003728419 
    ## ... Similar to previous best
    ## Run 214 stress 0.02492199 
    ## Run 215 stress 0.01879885 
    ## ... Procrustes: rmse 0.001100497  max resid 0.002267007 
    ## ... Similar to previous best
    ## Run 216 stress 0.01879601 
    ## ... Procrustes: rmse 0.0003457546  max resid 0.0007122924 
    ## ... Similar to previous best
    ## Run 217 stress 0.02486098 
    ## Run 218 stress 0.02509474 
    ## Run 219 stress 0.02520892 
    ## Run 220 stress 0.02509484 
    ## Run 221 stress 0.01879567 
    ## ... Procrustes: rmse 0.00106127  max resid 0.002187541 
    ## ... Similar to previous best
    ## Run 222 stress 0.01879587 
    ## ... Procrustes: rmse 0.0003052201  max resid 0.0006287797 
    ## ... Similar to previous best
    ## Run 223 stress 0.02486133 
    ## Run 224 stress 0.01879559 
    ## ... Procrustes: rmse 0.0001828392  max resid 0.0003767065 
    ## ... Similar to previous best
    ## Run 225 stress 0.2038604 
    ## Run 226 stress 0.02509439 
    ## Run 227 stress 0.02494974 
    ## Run 228 stress 0.02509456 
    ## Run 229 stress 0.01881671 
    ## ... Procrustes: rmse 0.003148722  max resid 0.006481208 
    ## ... Similar to previous best
    ## Run 230 stress 0.01879577 
    ## ... Procrustes: rmse 0.0002460953  max resid 0.0005065276 
    ## ... Similar to previous best
    ## Run 231 stress 0.02520898 
    ## Run 232 stress 0.02492169 
    ## Run 233 stress 0.02520906 
    ## Run 234 stress 0.01879595 
    ## ... Procrustes: rmse 0.0003269629  max resid 0.0006736058 
    ## ... Similar to previous best
    ## Run 235 stress 0.02045748 
    ## Run 236 stress 0.01880882 
    ## ... Procrustes: rmse 0.002461273  max resid 0.005070401 
    ## ... Similar to previous best
    ## Run 237 stress 0.01879576 
    ## ... Procrustes: rmse 0.001103081  max resid 0.002273415 
    ## ... Similar to previous best
    ## Run 238 stress 0.02492196 
    ## Run 239 stress 0.01879539 
    ## ... Procrustes: rmse 7.1401e-05  max resid 0.0001471335 
    ## ... Similar to previous best
    ## Run 240 stress 0.01880047 
    ## ... Procrustes: rmse 0.00124452  max resid 0.002576434 
    ## ... Similar to previous best
    ## Run 241 stress 0.02509431 
    ## Run 242 stress 0.01879767 
    ## ... Procrustes: rmse 0.0008459954  max resid 0.001742887 
    ## ... Similar to previous best
    ## Run 243 stress 0.02492202 
    ## Run 244 stress 0.02492178 
    ## Run 245 stress 0.01885299 
    ## ... Procrustes: rmse 0.004302012  max resid 0.008496267 
    ## ... Similar to previous best
    ## Run 246 stress 0.0249219 
    ## Run 247 stress 0.01879582 
    ## ... Procrustes: rmse 0.0002864025  max resid 0.0005900377 
    ## ... Similar to previous best
    ## Run 248 stress 0.02493667 
    ## Run 249 stress 0.02492174 
    ## Run 250 stress 0.01879755 
    ## ... Procrustes: rmse 0.0008202613  max resid 0.001690067 
    ## ... Similar to previous best
    ## Run 251 stress 0.0187983 
    ## ... Procrustes: rmse 0.0009883638  max resid 0.00203787 
    ## ... Similar to previous best
    ## Run 252 stress 0.01879706 
    ## ... Procrustes: rmse 0.0006919169  max resid 0.001425267 
    ## ... Similar to previous best
    ## Run 253 stress 0.0250947 
    ## Run 254 stress 0.01879603 
    ## ... Procrustes: rmse 0.0003597233  max resid 0.0007411134 
    ## ... Similar to previous best
    ## Run 255 stress 0.02493691 
    ## Run 256 stress 0.025209 
    ## Run 257 stress 0.02486119 
    ## Run 258 stress 0.02486115 
    ## Run 259 stress 0.3726997 
    ## Run 260 stress 0.02486131 
    ## Run 261 stress 0.02520906 
    ## Run 262 stress 0.02509459 
    ## Run 263 stress 0.01879583 
    ## ... Procrustes: rmse 0.0002902701  max resid 0.0005980285 
    ## ... Similar to previous best
    ## Run 264 stress 0.02492186 
    ## Run 265 stress 0.02509466 
    ## Run 266 stress 0.01879587 
    ## ... Procrustes: rmse 0.0003020393  max resid 0.0006220766 
    ## ... Similar to previous best
    ## Run 267 stress 0.02509473 
    ## Run 268 stress 0.02492173 
    ## Run 269 stress 0.02520903 
    ## Run 270 stress 0.01888861 
    ## ... Procrustes: rmse 0.007159065  max resid 0.0147118 
    ## Run 271 stress 0.01879598 
    ## ... Procrustes: rmse 0.001192302  max resid 0.002456768 
    ## ... Similar to previous best
    ## Run 272 stress 0.01879645 
    ## ... Procrustes: rmse 0.0005168323  max resid 0.001064723 
    ## ... Similar to previous best
    ## Run 273 stress 0.02486108 
    ## Run 274 stress 0.02520891 
    ## Run 275 stress 0.01886995 
    ## ... Procrustes: rmse 0.005436166  max resid 0.01088259 
    ## Run 276 stress 0.02509462 
    ## Run 277 stress 0.02520917 
    ## Run 278 stress 0.01879568 
    ## ... Procrustes: rmse 0.0002229732  max resid 0.0004594215 
    ## ... Similar to previous best
    ## Run 279 stress 0.02509446 
    ## Run 280 stress 0.01900566 
    ## ... Procrustes: rmse 0.01101581  max resid 0.02258892 
    ## Run 281 stress 0.0187953 
    ## ... Procrustes: rmse 1.483289e-05  max resid 3.002932e-05 
    ## ... Similar to previous best
    ## Run 282 stress 0.01887105 
    ## ... Procrustes: rmse 0.005502362  max resid 0.01101988 
    ## Run 283 stress 0.0249736 
    ## Run 284 stress 0.02492175 
    ## Run 285 stress 0.025209 
    ## Run 286 stress 0.02509461 
    ## Run 287 stress 0.01884602 
    ## ... Procrustes: rmse 0.005140835  max resid 0.01057247 
    ## Run 288 stress 0.0249216 
    ## Run 289 stress 0.02509462 
    ## Run 290 stress 0.02509465 
    ## Run 291 stress 0.0193158 
    ## Run 292 stress 0.02492174 
    ## Run 293 stress 0.02509476 
    ## Run 294 stress 0.02509468 
    ## Run 295 stress 0.01883493 
    ## ... Procrustes: rmse 0.0026868  max resid 0.005362395 
    ## ... Similar to previous best
    ## Run 296 stress 0.01882856 
    ## ... Procrustes: rmse 0.001856171  max resid 0.004361044 
    ## ... Similar to previous best
    ## Run 297 stress 0.01879586 
    ## ... Procrustes: rmse 0.0003011161  max resid 0.0006203892 
    ## ... Similar to previous best
    ## Run 298 stress 0.02510789 
    ## Run 299 stress 0.02509462 
    ## Run 300 stress 0.01885975 
    ## ... Procrustes: rmse 0.004668919  max resid 0.009272059 
    ## ... Similar to previous best
    ## Run 301 stress 0.02492186 
    ## Run 302 stress 0.01879584 
    ## ... Procrustes: rmse 0.000281719  max resid 0.0005804425 
    ## ... Similar to previous best
    ## Run 303 stress 0.02509453 
    ## Run 304 stress 0.02509462 
    ## Run 305 stress 0.024922 
    ## Run 306 stress 0.02509483 
    ## Run 307 stress 0.02486148 
    ## Run 308 stress 0.01883442 
    ## ... Procrustes: rmse 0.00266178  max resid 0.005339793 
    ## ... Similar to previous best
    ## Run 309 stress 0.01879552 
    ## ... Procrustes: rmse 0.0001445729  max resid 0.0002978137 
    ## ... Similar to previous best
    ## Run 310 stress 0.02509464 
    ## Run 311 stress 0.02520898 
    ## Run 312 stress 0.02520897 
    ## Run 313 stress 0.0250946 
    ## Run 314 stress 0.01879588 
    ## ... Procrustes: rmse 0.000312775  max resid 0.0006443407 
    ## ... Similar to previous best
    ## Run 315 stress 0.01883597 
    ## ... Procrustes: rmse 0.002841677  max resid 0.005466998 
    ## ... Similar to previous best
    ## Run 316 stress 0.01879591 
    ## ... Procrustes: rmse 0.0003219357  max resid 0.0006632168 
    ## ... Similar to previous best
    ## Run 317 stress 0.01879541 
    ## ... Procrustes: rmse 8.184698e-05  max resid 0.000168612 
    ## ... Similar to previous best
    ## Run 318 stress 0.02509472 
    ## Run 319 stress 0.02492183 
    ## Run 320 stress 0.01879586 
    ## ... Procrustes: rmse 0.0003009098  max resid 0.0006197691 
    ## ... Similar to previous best
    ## Run 321 stress 0.01879585 
    ## ... Procrustes: rmse 0.0003004641  max resid 0.0006189716 
    ## ... Similar to previous best
    ## Run 322 stress 0.02509469 
    ## Run 323 stress 0.01879567 
    ## ... Procrustes: rmse 0.0002176615  max resid 0.0004484222 
    ## ... Similar to previous best
    ## Run 324 stress 0.01880926 
    ## ... Procrustes: rmse 0.00246372  max resid 0.005076935 
    ## ... Similar to previous best
    ## Run 325 stress 0.01890129 
    ## ... Procrustes: rmse 0.007647936  max resid 0.01571403 
    ## Run 326 stress 0.0187957 
    ## ... Procrustes: rmse 0.0002310639  max resid 0.0004760912 
    ## ... Similar to previous best
    ## Run 327 stress 0.02509454 
    ## Run 328 stress 0.01887677 
    ## ... Procrustes: rmse 0.005726791  max resid 0.0114875 
    ## Run 329 stress 0.0188504 
    ## ... Procrustes: rmse 0.00538829  max resid 0.011083 
    ## Run 330 stress 0.02509436 
    ## Run 331 stress 0.02509461 
    ## Run 332 stress 0.018862 
    ## ... Procrustes: rmse 0.004880315  max resid 0.009711822 
    ## ... Similar to previous best
    ## Run 333 stress 0.01879565 
    ## ... Procrustes: rmse 0.0002021389  max resid 0.0004160872 
    ## ... Similar to previous best
    ## Run 334 stress 0.01885547 
    ## ... Procrustes: rmse 0.004517253  max resid 0.008951081 
    ## ... Similar to previous best
    ## Run 335 stress 0.01882569 
    ## ... Procrustes: rmse 0.003896993  max resid 0.008020514 
    ## ... Similar to previous best
    ## Run 336 stress 0.01879581 
    ## ... Procrustes: rmse 0.0002809092  max resid 0.000578739 
    ## ... Similar to previous best
    ## Run 337 stress 0.02486115 
    ## Run 338 stress 0.02492176 
    ## Run 339 stress 0.01879597 
    ## ... Procrustes: rmse 0.0003473312  max resid 0.0007154764 
    ## ... Similar to previous best
    ## Run 340 stress 0.02486117 
    ## Run 341 stress 0.01879607 
    ## ... Procrustes: rmse 0.001210355  max resid 0.002493013 
    ## ... Similar to previous best
    ## Run 342 stress 0.02509441 
    ## Run 343 stress 0.02486137 
    ## Run 344 stress 0.01900889 
    ## ... Procrustes: rmse 0.01078464  max resid 0.0219699 
    ## Run 345 stress 0.01879591 
    ## ... Procrustes: rmse 0.0003233748  max resid 0.0006661578 
    ## ... Similar to previous best
    ## Run 346 stress 0.01879558 
    ## ... Procrustes: rmse 0.001018903  max resid 0.002100458 
    ## ... Similar to previous best
    ## Run 347 stress 0.0249219 
    ## Run 348 stress 0.01882755 
    ## ... Procrustes: rmse 0.00400591  max resid 0.008242753 
    ## ... Similar to previous best
    ## Run 349 stress 0.02486114 
    ## Run 350 stress 0.01883115 
    ## ... Procrustes: rmse 0.002137715  max resid 0.004888417 
    ## ... Similar to previous best
    ## Run 351 stress 0.0249217 
    ## Run 352 stress 0.01890285 
    ## ... Procrustes: rmse 0.007722043  max resid 0.0158631 
    ## Run 353 stress 0.02492173 
    ## Run 354 stress 0.3654897 
    ## Run 355 stress 0.02520913 
    ## Run 356 stress 0.02486104 
    ## Run 357 stress 0.01879589 
    ## ... Procrustes: rmse 0.000303837  max resid 0.0006255738 
    ## ... Similar to previous best
    ## Run 358 stress 0.02509446 
    ## Run 359 stress 0.01879561 
    ## ... Procrustes: rmse 0.0001908913  max resid 0.0003932449 
    ## ... Similar to previous best
    ## Run 360 stress 0.01879571 
    ## ... Procrustes: rmse 0.001082456  max resid 0.002231051 
    ## ... Similar to previous best
    ## Run 361 stress 0.01894789 
    ## ... Procrustes: rmse 0.009232758  max resid 0.0189556 
    ## Run 362 stress 0.02486114 
    ## Run 363 stress 0.01882952 
    ## ... Procrustes: rmse 0.001979942  max resid 0.004618973 
    ## ... Similar to previous best
    ## Run 364 stress 0.01883687 
    ## ... Procrustes: rmse 0.004615339  max resid 0.009496556 
    ## ... Similar to previous best
    ## Run 365 stress 0.01883917 
    ## ... Procrustes: rmse 0.003182245  max resid 0.006096006 
    ## ... Similar to previous best
    ## Run 366 stress 0.01884638 
    ## ... Procrustes: rmse 0.005176137  max resid 0.01064639 
    ## Run 367 stress 0.01879593 
    ## ... Procrustes: rmse 0.001171106  max resid 0.0024132 
    ## ... Similar to previous best
    ## Run 368 stress 0.01879583 
    ## ... Procrustes: rmse 0.001130682  max resid 0.00233015 
    ## ... Similar to previous best
    ## Run 369 stress 0.01880476 
    ## ... Procrustes: rmse 0.001965873  max resid 0.004051685 
    ## ... Similar to previous best
    ## Run 370 stress 0.01879568 
    ## ... Procrustes: rmse 0.0002216401  max resid 0.0004566698 
    ## ... Similar to previous best
    ## Run 371 stress 0.01879797 
    ## ... Procrustes: rmse 0.0009173876  max resid 0.00189065 
    ## ... Similar to previous best
    ## Run 372 stress 0.01879598 
    ## ... Procrustes: rmse 0.0003436942  max resid 0.0007081 
    ## ... Similar to previous best
    ## Run 373 stress 0.02492196 
    ## Run 374 stress 0.02509474 
    ## Run 375 stress 0.01882719 
    ## ... Procrustes: rmse 0.003990251  max resid 0.008212708 
    ## ... Similar to previous best
    ## Run 376 stress 0.02492173 
    ## Run 377 stress 0.01879556 
    ## ... Procrustes: rmse 0.001002675  max resid 0.0020671 
    ## ... Similar to previous best
    ## Run 378 stress 0.02509457 
    ## Run 379 stress 0.0187957 
    ## ... Procrustes: rmse 0.001077142  max resid 0.002220058 
    ## ... Similar to previous best
    ## Run 380 stress 0.02492172 
    ## Run 381 stress 0.01879554 
    ## ... Procrustes: rmse 0.0001480952  max resid 0.000304762 
    ## ... Similar to previous best
    ## Run 382 stress 0.01880778 
    ## ... Procrustes: rmse 0.001959343  max resid 0.004037141 
    ## ... Similar to previous best
    ## Run 383 stress 0.0250947 
    ## Run 384 stress 0.01879586 
    ## ... Procrustes: rmse 0.0002899581  max resid 0.0005969441 
    ## ... Similar to previous best
    ## Run 385 stress 0.0188344 
    ## ... Procrustes: rmse 0.0024498  max resid 0.005181813 
    ## ... Similar to previous best
    ## Run 386 stress 0.01881956 
    ## ... Procrustes: rmse 0.001624534  max resid 0.003073929 
    ## ... Similar to previous best
    ## Run 387 stress 0.02486115 
    ## Run 388 stress 0.01879863 
    ## ... Procrustes: rmse 0.001035802  max resid 0.00213597 
    ## ... Similar to previous best
    ## Run 389 stress 0.01879603 
    ## ... Procrustes: rmse 0.000338152  max resid 0.000695911 
    ## ... Similar to previous best
    ## Run 390 stress 0.01885326 
    ## ... Procrustes: rmse 0.00403353  max resid 0.008273578 
    ## ... Similar to previous best
    ## Run 391 stress 0.02509461 
    ## Run 392 stress 0.02492544 
    ## Run 393 stress 0.02492163 
    ## Run 394 stress 0.01887907 
    ## ... Procrustes: rmse 0.005947935  max resid 0.01195402 
    ## Run 395 stress 0.01879587 
    ## ... Procrustes: rmse 0.0002951233  max resid 0.0006076647 
    ## ... Similar to previous best
    ## Run 396 stress 0.0249218 
    ## Run 397 stress 0.3156979 
    ## Run 398 stress 0.02492184 
    ## Run 399 stress 0.0250946 
    ## Run 400 stress 0.02486102 
    ## Run 401 stress 0.02492192 
    ## Run 402 stress 0.01879587 
    ## ... Procrustes: rmse 0.000303844  max resid 0.0006260052 
    ## ... Similar to previous best
    ## Run 403 stress 0.02509459 
    ## Run 404 stress 0.02492197 
    ## Run 405 stress 0.01880004 
    ## ... Procrustes: rmse 0.001254891  max resid 0.002584644 
    ## ... Similar to previous best
    ## Run 406 stress 0.01880593 
    ## ... Procrustes: rmse 0.002144378  max resid 0.004419994 
    ## ... Similar to previous best
    ## Run 407 stress 0.02520882 
    ## Run 408 stress 0.01879549 
    ## ... Procrustes: rmse 0.0009695532  max resid 0.001999174 
    ## ... Similar to previous best
    ## Run 409 stress 0.01879596 
    ## ... Procrustes: rmse 0.0003348875  max resid 0.0006899542 
    ## ... Similar to previous best
    ## Run 410 stress 0.02509478 
    ## Run 411 stress 0.02492199 
    ## Run 412 stress 0.0191461 
    ## ... Procrustes: rmse 0.01412517  max resid 0.02900314 
    ## Run 413 stress 0.01879611 
    ## ... Procrustes: rmse 0.0003834289  max resid 0.000789885 
    ## ... Similar to previous best
    ## Run 414 stress 0.01879608 
    ## ... Procrustes: rmse 0.0003722707  max resid 0.0007669349 
    ## ... Similar to previous best
    ## Run 415 stress 0.01883175 
    ## ... Procrustes: rmse 0.002317074  max resid 0.005072084 
    ## ... Similar to previous best
    ## Run 416 stress 0.02486116 
    ## Run 417 stress 0.02492191 
    ## Run 418 stress 0.02486126 
    ## Run 419 stress 0.01879562 
    ## ... Procrustes: rmse 0.0001911136  max resid 0.000393683 
    ## ... Similar to previous best
    ## Run 420 stress 0.02492178 
    ## Run 421 stress 0.01879572 
    ## ... Procrustes: rmse 0.0002342563  max resid 0.0004826643 
    ## ... Similar to previous best
    ## Run 422 stress 0.01887894 
    ## ... Procrustes: rmse 0.005943084  max resid 0.01194353 
    ## Run 423 stress 0.0187958 
    ## ... Procrustes: rmse 0.001104404  max resid 0.002276059 
    ## ... Similar to previous best
    ## Run 424 stress 0.01885816 
    ## ... Procrustes: rmse 0.004706668  max resid 0.009349692 
    ## ... Similar to previous best
    ## Run 425 stress 0.01879579 
    ## ... Procrustes: rmse 0.0002718161  max resid 0.0005598697 
    ## ... Similar to previous best
    ## Run 426 stress 0.0188917 
    ## ... Procrustes: rmse 0.00658657  max resid 0.01328829 
    ## Run 427 stress 0.01879585 
    ## ... Procrustes: rmse 0.0002933697  max resid 0.0006044386 
    ## ... Similar to previous best
    ## Run 428 stress 0.01879604 
    ## ... Procrustes: rmse 0.0003590011  max resid 0.0007396046 
    ## ... Similar to previous best
    ## Run 429 stress 0.0188015 
    ## ... Procrustes: rmse 0.001555786  max resid 0.003205807 
    ## ... Similar to previous best
    ## Run 430 stress 0.01888179 
    ## ... Procrustes: rmse 0.005929476  max resid 0.01191727 
    ## Run 431 stress 0.02239632 
    ## Run 432 stress 0.2240558 
    ## Run 433 stress 0.0188604 
    ## ... Procrustes: rmse 0.005902306  max resid 0.0121361 
    ## Run 434 stress 0.01879599 
    ## ... Procrustes: rmse 0.001184103  max resid 0.002439906 
    ## ... Similar to previous best
    ## Run 435 stress 0.02492184 
    ## Run 436 stress 0.02509429 
    ## Run 437 stress 0.0250946 
    ## Run 438 stress 0.0250944 
    ## Run 439 stress 0.01905865 
    ## ... Procrustes: rmse 0.01189389  max resid 0.02427092 
    ## Run 440 stress 0.02509442 
    ## Run 441 stress 0.02520901 
    ## Run 442 stress 0.01879556 
    ## ... Procrustes: rmse 0.001008882  max resid 0.002079912 
    ## ... Similar to previous best
    ## Run 443 stress 0.01879566 
    ## ... Procrustes: rmse 0.001056841  max resid 0.002178422 
    ## ... Similar to previous best
    ## Run 444 stress 0.01883662 
    ## ... Procrustes: rmse 0.002915436  max resid 0.005517741 
    ## ... Similar to previous best
    ## Run 445 stress 0.02486137 
    ## Run 446 stress 0.01887513 
    ## ... Procrustes: rmse 0.006579062  max resid 0.01352543 
    ## Run 447 stress 0.02509469 
    ## Run 448 stress 0.01879586 
    ## ... Procrustes: rmse 0.0003030983  max resid 0.0006244393 
    ## ... Similar to previous best
    ## Run 449 stress 0.01879581 
    ## ... Procrustes: rmse 0.0002829811  max resid 0.0005829781 
    ## ... Similar to previous best
    ## Run 450 stress 0.01879591 
    ## ... Procrustes: rmse 0.001166824  max resid 0.00240441 
    ## ... Similar to previous best
    ## Run 451 stress 0.01885338 
    ## ... Procrustes: rmse 0.004349822  max resid 0.008597364 
    ## ... Similar to previous best
    ## Run 452 stress 0.01890385 
    ## ... Procrustes: rmse 0.007146638  max resid 0.01445531 
    ## Run 453 stress 0.02509468 
    ## Run 454 stress 0.02520871 
    ## Run 455 stress 0.02509442 
    ## Run 456 stress 0.01879564 
    ## ... Procrustes: rmse 0.001046763  max resid 0.002157739 
    ## ... Similar to previous best
    ## Run 457 stress 0.02509431 
    ## Run 458 stress 0.01883375 
    ## ... Procrustes: rmse 0.004434959  max resid 0.009125335 
    ## ... Similar to previous best
    ## Run 459 stress 0.01880411 
    ## ... Procrustes: rmse 0.001919635  max resid 0.003956074 
    ## ... Similar to previous best
    ## Run 460 stress 0.02492175 
    ## Run 461 stress 0.02492196 
    ## Run 462 stress 0.01884596 
    ## ... Procrustes: rmse 0.003793631  max resid 0.007411257 
    ## ... Similar to previous best
    ## Run 463 stress 0.0187959 
    ## ... Procrustes: rmse 0.0003076068  max resid 0.0006333532 
    ## ... Similar to previous best
    ## Run 464 stress 0.380826 
    ## Run 465 stress 0.02509461 
    ## Run 466 stress 0.01880529 
    ## ... Procrustes: rmse 0.002064073  max resid 0.004253684 
    ## ... Similar to previous best
    ## Run 467 stress 0.01894154 
    ## ... Procrustes: rmse 0.009080669  max resid 0.01863953 
    ## Run 468 stress 0.01879584 
    ## ... Procrustes: rmse 0.0002959402  max resid 0.0006096663 
    ## ... Similar to previous best
    ## Run 469 stress 0.01879556 
    ## ... Procrustes: rmse 0.0001651331  max resid 0.0003402162 
    ## ... Similar to previous best
    ## Run 470 stress 0.3726972 
    ## Run 471 stress 0.01880241 
    ## ... Procrustes: rmse 0.001462109  max resid 0.002954093 
    ## ... Similar to previous best
    ## Run 472 stress 0.01889316 
    ## ... Procrustes: rmse 0.00618423  max resid 0.01245228 
    ## Run 473 stress 0.01888159 
    ## ... Procrustes: rmse 0.006080002  max resid 0.01223011 
    ## Run 474 stress 0.02492181 
    ## Run 475 stress 0.01879572 
    ## ... Procrustes: rmse 0.0002355936  max resid 0.0004856767 
    ## ... Similar to previous best
    ## Run 476 stress 0.02509461 
    ## Run 477 stress 0.02520895 
    ## Run 478 stress 0.01879574 
    ## ... Procrustes: rmse 0.001069829  max resid 0.002205083 
    ## ... Similar to previous best
    ## Run 479 stress 0.02486132 
    ## Run 480 stress 0.01885542 
    ## ... Procrustes: rmse 0.005654929  max resid 0.01162923 
    ## Run 481 stress 0.02486144 
    ## Run 482 stress 0.01879575 
    ## ... Procrustes: rmse 0.0002562932  max resid 0.0005280198 
    ## ... Similar to previous best
    ## Run 483 stress 0.02486151 
    ## Run 484 stress 0.01879555 
    ## ... Procrustes: rmse 0.000153825  max resid 0.0003169015 
    ## ... Similar to previous best
    ## Run 485 stress 0.01897202 
    ## ... Procrustes: rmse 0.009792871  max resid 0.0200827 
    ## Run 486 stress 0.02492165 
    ## Run 487 stress 0.01879567 
    ## ... Procrustes: rmse 0.0002209064  max resid 0.0004551468 
    ## ... Similar to previous best
    ## Run 488 stress 0.01879583 
    ## ... Procrustes: rmse 0.0002863776  max resid 0.0005898612 
    ## ... Similar to previous best
    ## Run 489 stress 0.01879594 
    ## ... Procrustes: rmse 0.0003215177  max resid 0.0006620163 
    ## ... Similar to previous best
    ## Run 490 stress 0.01885909 
    ## ... Procrustes: rmse 0.004772935  max resid 0.009489441 
    ## ... Similar to previous best
    ## Run 491 stress 0.02509444 
    ## Run 492 stress 0.01879569 
    ## ... Procrustes: rmse 0.0002294385  max resid 0.0004727356 
    ## ... Similar to previous best
    ## Run 493 stress 0.01891087 
    ## ... Procrustes: rmse 0.007457736  max resid 0.01510148 
    ## Run 494 stress 0.01881032 
    ## ... Procrustes: rmse 0.002578839  max resid 0.005313767 
    ## ... Similar to previous best
    ## Run 495 stress 0.01879835 
    ## ... Procrustes: rmse 0.0009870718  max resid 0.002033847 
    ## ... Similar to previous best
    ## Run 496 stress 0.01879578 
    ## ... Procrustes: rmse 0.0002573552  max resid 0.0005302197 
    ## ... Similar to previous best
    ## Run 497 stress 0.0249217 
    ## Run 498 stress 0.3727005 
    ## Run 499 stress 0.01880013 
    ## ... Procrustes: rmse 0.001334215  max resid 0.002749349 
    ## ... Similar to previous best
    ## Run 500 stress 0.0187973 
    ## ... Procrustes: rmse 0.0007584677  max resid 0.001562464 
    ## ... Similar to previous best
    ## *** Best solution repeated 196 times

``` r
# Mixed and stratified lakes
SD_beta_MS_NMDS <- metaMDS(SD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09308959 
    ## Run 2 stress 0.08440255 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01442919  max resid 0.04295692 
    ## Run 3 stress 0.1038044 
    ## Run 4 stress 0.08503537 
    ## Run 5 stress 0.09286083 
    ## Run 6 stress 0.09539192 
    ## Run 7 stress 0.08773487 
    ## Run 8 stress 0.08773471 
    ## Run 9 stress 0.0944731 
    ## Run 10 stress 0.09416153 
    ## Run 11 stress 0.09030394 
    ## Run 12 stress 0.1052851 
    ## Run 13 stress 0.08503462 
    ## Run 14 stress 0.09969967 
    ## Run 15 stress 0.08773488 
    ## Run 16 stress 0.09268341 
    ## Run 17 stress 0.09407976 
    ## Run 18 stress 0.09407982 
    ## Run 19 stress 0.09145325 
    ## Run 20 stress 0.09159096 
    ## Run 21 stress 0.08503505 
    ## Run 22 stress 0.09464485 
    ## Run 23 stress 0.08973863 
    ## Run 24 stress 0.0930897 
    ## Run 25 stress 0.09408005 
    ## Run 26 stress 0.08440261 
    ## ... Procrustes: rmse 8.026312e-05  max resid 0.0001476555 
    ## ... Similar to previous best
    ## Run 27 stress 0.09030396 
    ## Run 28 stress 0.08973872 
    ## Run 29 stress 0.09977535 
    ## Run 30 stress 0.09407981 
    ## Run 31 stress 0.09416143 
    ## Run 32 stress 0.09268429 
    ## Run 33 stress 0.08503455 
    ## Run 34 stress 0.09268361 
    ## Run 35 stress 0.097212 
    ## Run 36 stress 0.09159094 
    ## Run 37 stress 0.1038048 
    ## Run 38 stress 0.09159109 
    ## Run 39 stress 0.08973892 
    ## Run 40 stress 0.09268352 
    ## Run 41 stress 0.09145321 
    ## Run 42 stress 0.09268365 
    ## Run 43 stress 0.08503496 
    ## Run 44 stress 0.09159087 
    ## Run 45 stress 0.09374331 
    ## Run 46 stress 0.09407991 
    ## Run 47 stress 0.09286104 
    ## Run 48 stress 0.09416135 
    ## Run 49 stress 0.09465913 
    ## Run 50 stress 0.09030394 
    ## Run 51 stress 0.09721209 
    ## Run 52 stress 0.09465915 
    ## Run 53 stress 0.09465906 
    ## Run 54 stress 0.09337236 
    ## Run 55 stress 0.09417791 
    ## Run 56 stress 0.08773489 
    ## Run 57 stress 0.09337283 
    ## Run 58 stress 0.09030407 
    ## Run 59 stress 0.09760824 
    ## Run 60 stress 0.09145331 
    ## Run 61 stress 0.08440264 
    ## ... Procrustes: rmse 0.0003312032  max resid 0.0006500353 
    ## ... Similar to previous best
    ## Run 62 stress 0.09407996 
    ## Run 63 stress 0.09669853 
    ## Run 64 stress 0.08440261 
    ## ... Procrustes: rmse 7.067098e-05  max resid 0.0001359706 
    ## ... Similar to previous best
    ## Run 65 stress 0.08440269 
    ## ... Procrustes: rmse 0.0001432178  max resid 0.0002988624 
    ## ... Similar to previous best
    ## Run 66 stress 0.08973869 
    ## Run 67 stress 0.09030399 
    ## Run 68 stress 0.08440259 
    ## ... Procrustes: rmse 5.285584e-05  max resid 0.0001067367 
    ## ... Similar to previous best
    ## Run 69 stress 0.09268387 
    ## Run 70 stress 0.08440261 
    ## ... Procrustes: rmse 8.349081e-05  max resid 0.000207294 
    ## ... Similar to previous best
    ## Run 71 stress 0.08503537 
    ## Run 72 stress 0.1052849 
    ## Run 73 stress 0.09286083 
    ## Run 74 stress 0.08503491 
    ## Run 75 stress 0.0877347 
    ## Run 76 stress 0.08973883 
    ## Run 77 stress 0.09168929 
    ## Run 78 stress 0.09539195 
    ## Run 79 stress 0.09030399 
    ## Run 80 stress 0.3148536 
    ## Run 81 stress 0.09407969 
    ## Run 82 stress 0.105455 
    ## Run 83 stress 0.09465911 
    ## Run 84 stress 0.08503477 
    ## Run 85 stress 0.09407997 
    ## Run 86 stress 0.09030408 
    ## Run 87 stress 0.08973874 
    ## Run 88 stress 0.09977532 
    ## Run 89 stress 0.09286083 
    ## Run 90 stress 0.09030411 
    ## Run 91 stress 0.09374305 
    ## Run 92 stress 0.08503457 
    ## Run 93 stress 0.08503476 
    ## Run 94 stress 0.09465905 
    ## Run 95 stress 0.0959054 
    ## Run 96 stress 0.09030405 
    ## Run 97 stress 0.08773485 
    ## Run 98 stress 0.09407975 
    ## Run 99 stress 0.09308946 
    ## Run 100 stress 0.09590887 
    ## Run 101 stress 0.09407969 
    ## Run 102 stress 0.09539209 
    ## Run 103 stress 0.09465905 
    ## Run 104 stress 0.08440253 
    ## ... New best solution
    ## ... Procrustes: rmse 3.776146e-05  max resid 6.474153e-05 
    ## ... Similar to previous best
    ## Run 105 stress 0.09407126 
    ## Run 106 stress 0.09407984 
    ## Run 107 stress 0.1038045 
    ## Run 108 stress 0.1053425 
    ## Run 109 stress 0.09380584 
    ## Run 110 stress 0.08773471 
    ## Run 111 stress 0.09168942 
    ## Run 112 stress 0.0916893 
    ## Run 113 stress 0.09407118 
    ## Run 114 stress 0.08973893 
    ## Run 115 stress 0.09030405 
    ## Run 116 stress 0.09535463 
    ## Run 117 stress 0.09159086 
    ## Run 118 stress 0.09030398 
    ## Run 119 stress 0.0937431 
    ## Run 120 stress 0.08773472 
    ## Run 121 stress 0.09308963 
    ## Run 122 stress 0.09535403 
    ## Run 123 stress 0.09416146 
    ## Run 124 stress 0.09669911 
    ## Run 125 stress 0.09408007 
    ## Run 126 stress 0.09760829 
    ## Run 127 stress 0.09535428 
    ## Run 128 stress 0.0940342 
    ## Run 129 stress 0.09168932 
    ## Run 130 stress 0.09308955 
    ## Run 131 stress 0.08503605 
    ## Run 132 stress 0.09417786 
    ## Run 133 stress 0.1026215 
    ## Run 134 stress 0.09268333 
    ## Run 135 stress 0.09464446 
    ## Run 136 stress 0.0953549 
    ## Run 137 stress 0.09417791 
    ## Run 138 stress 0.09286085 
    ## Run 139 stress 0.09168945 
    ## Run 140 stress 0.0903041 
    ## Run 141 stress 0.09447319 
    ## Run 142 stress 0.09337282 
    ## Run 143 stress 0.09159087 
    ## Run 144 stress 0.09030408 
    ## Run 145 stress 0.0938185 
    ## Run 146 stress 0.09030401 
    ## Run 147 stress 0.09416146 
    ## Run 148 stress 0.09969971 
    ## Run 149 stress 0.09030397 
    ## Run 150 stress 0.08440265 
    ## ... Procrustes: rmse 0.0001469626  max resid 0.0002985776 
    ## ... Similar to previous best
    ## Run 151 stress 0.09145324 
    ## Run 152 stress 0.0850362 
    ## Run 153 stress 0.09970004 
    ## Run 154 stress 0.08773476 
    ## Run 155 stress 0.09268341 
    ## Run 156 stress 0.09539207 
    ## Run 157 stress 0.0940798 
    ## Run 158 stress 0.09030398 
    ## Run 159 stress 0.09286084 
    ## Run 160 stress 0.08503662 
    ## Run 161 stress 0.0877347 
    ## Run 162 stress 0.09535463 
    ## Run 163 stress 0.09159083 
    ## Run 164 stress 0.1095796 
    ## Run 165 stress 0.09374346 
    ## Run 166 stress 0.085035 
    ## Run 167 stress 0.09539205 
    ## Run 168 stress 0.0930896 
    ## Run 169 stress 0.09407111 
    ## Run 170 stress 0.09030396 
    ## Run 171 stress 0.09168929 
    ## Run 172 stress 0.09168928 
    ## Run 173 stress 0.09145325 
    ## Run 174 stress 0.09030413 
    ## Run 175 stress 0.09407121 
    ## Run 176 stress 0.09407976 
    ## Run 177 stress 0.09969981 
    ## Run 178 stress 0.09308964 
    ## Run 179 stress 0.09407124 
    ## Run 180 stress 0.09030408 
    ## Run 181 stress 0.09712976 
    ## Run 182 stress 0.08503685 
    ## Run 183 stress 0.09286085 
    ## Run 184 stress 0.09268388 
    ## Run 185 stress 0.09407976 
    ## Run 186 stress 0.09445464 
    ## Run 187 stress 0.09465904 
    ## Run 188 stress 0.09407108 
    ## Run 189 stress 0.0850348 
    ## Run 190 stress 0.09286084 
    ## Run 191 stress 0.09030406 
    ## Run 192 stress 0.09030395 
    ## Run 193 stress 0.09535474 
    ## Run 194 stress 0.09030401 
    ## Run 195 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 5.130053e-05  max resid 0.0001067403 
    ## ... Similar to previous best
    ## Run 196 stress 0.09337276 
    ## Run 197 stress 0.09408 
    ## Run 198 stress 0.08503562 
    ## Run 199 stress 0.08773466 
    ## Run 200 stress 0.0976082 
    ## Run 201 stress 0.08503496 
    ## Run 202 stress 0.09030409 
    ## Run 203 stress 0.09030406 
    ## Run 204 stress 0.09969978 
    ## Run 205 stress 0.09030393 
    ## Run 206 stress 0.08440263 
    ## ... Procrustes: rmse 0.0001693317  max resid 0.0003620045 
    ## ... Similar to previous best
    ## Run 207 stress 0.09407991 
    ## Run 208 stress 0.09535502 
    ## Run 209 stress 0.09407985 
    ## Run 210 stress 0.09030403 
    ## Run 211 stress 0.08440265 
    ## ... Procrustes: rmse 0.0001908776  max resid 0.0004031084 
    ## ... Similar to previous best
    ## Run 212 stress 0.09465909 
    ## Run 213 stress 0.09268363 
    ## Run 214 stress 0.1038086 
    ## Run 215 stress 0.09416135 
    ## Run 216 stress 0.09407966 
    ## Run 217 stress 0.09286083 
    ## Run 218 stress 0.09308957 
    ## Run 219 stress 0.08773465 
    ## Run 220 stress 0.09969985 
    ## Run 221 stress 0.0946592 
    ## Run 222 stress 0.08440257 
    ## ... Procrustes: rmse 0.0001184359  max resid 0.0002224409 
    ## ... Similar to previous best
    ## Run 223 stress 0.09374197 
    ## Run 224 stress 0.09308948 
    ## Run 225 stress 0.0916893 
    ## Run 226 stress 0.09590812 
    ## Run 227 stress 0.09268393 
    ## Run 228 stress 0.09286088 
    ## Run 229 stress 0.09030413 
    ## Run 230 stress 0.09407977 
    ## Run 231 stress 0.1004611 
    ## Run 232 stress 0.09308954 
    ## Run 233 stress 0.0944732 
    ## Run 234 stress 0.09407986 
    ## Run 235 stress 0.09145323 
    ## Run 236 stress 0.3217457 
    ## Run 237 stress 0.08973864 
    ## Run 238 stress 0.09268387 
    ## Run 239 stress 0.09145324 
    ## Run 240 stress 0.08503471 
    ## Run 241 stress 0.09286084 
    ## Run 242 stress 0.09374204 
    ## Run 243 stress 0.0877348 
    ## Run 244 stress 0.1038048 
    ## Run 245 stress 0.09030408 
    ## Run 246 stress 0.09168935 
    ## Run 247 stress 0.09417774 
    ## Run 248 stress 0.09168937 
    ## Run 249 stress 0.09407988 
    ## Run 250 stress 0.08440262 
    ## ... Procrustes: rmse 0.0001510835  max resid 0.0003178509 
    ## ... Similar to previous best
    ## Run 251 stress 0.09407121 
    ## Run 252 stress 0.09268391 
    ## Run 253 stress 0.09408009 
    ## Run 254 stress 0.09464447 
    ## Run 255 stress 0.09445463 
    ## Run 256 stress 0.09337308 
    ## Run 257 stress 0.08503488 
    ## Run 258 stress 0.2955827 
    ## Run 259 stress 0.09145335 
    ## Run 260 stress 0.2908825 
    ## Run 261 stress 0.1038084 
    ## Run 262 stress 0.09168952 
    ## Run 263 stress 0.08503579 
    ## Run 264 stress 0.09407119 
    ## Run 265 stress 0.09408005 
    ## Run 266 stress 0.09168938 
    ## Run 267 stress 0.09030399 
    ## Run 268 stress 0.09145322 
    ## Run 269 stress 0.09407987 
    ## Run 270 stress 0.09159084 
    ## Run 271 stress 0.09337292 
    ## Run 272 stress 0.09337238 
    ## Run 273 stress 0.09721197 
    ## Run 274 stress 0.09908183 
    ## Run 275 stress 0.08440254 
    ## ... Procrustes: rmse 4.843278e-05  max resid 0.0001000137 
    ## ... Similar to previous best
    ## Run 276 stress 0.08773473 
    ## Run 277 stress 0.0940802 
    ## Run 278 stress 0.09030393 
    ## Run 279 stress 0.09590812 
    ## Run 280 stress 0.08503592 
    ## Run 281 stress 0.08973862 
    ## Run 282 stress 0.1038054 
    ## Run 283 stress 0.1026215 
    ## Run 284 stress 0.1019601 
    ## Run 285 stress 0.08773468 
    ## Run 286 stress 0.09464425 
    ## Run 287 stress 0.09030403 
    ## Run 288 stress 0.09145325 
    ## Run 289 stress 0.09374332 
    ## Run 290 stress 0.09407993 
    ## Run 291 stress 0.09308985 
    ## Run 292 stress 0.09030397 
    ## Run 293 stress 0.0940054 
    ## Run 294 stress 0.09465907 
    ## Run 295 stress 0.09407969 
    ## Run 296 stress 0.09416141 
    ## Run 297 stress 0.09407115 
    ## Run 298 stress 0.09268407 
    ## Run 299 stress 0.09145326 
    ## Run 300 stress 0.09159087 
    ## Run 301 stress 0.3040735 
    ## Run 302 stress 0.09465907 
    ## Run 303 stress 0.09416145 
    ## Run 304 stress 0.08973885 
    ## Run 305 stress 0.09145334 
    ## Run 306 stress 0.09145322 
    ## Run 307 stress 0.09407965 
    ## Run 308 stress 0.0926834 
    ## Run 309 stress 0.09268329 
    ## Run 310 stress 0.09308973 
    ## Run 311 stress 0.08503506 
    ## Run 312 stress 0.09464465 
    ## Run 313 stress 0.09030401 
    ## Run 314 stress 0.09308971 
    ## Run 315 stress 0.08440257 
    ## ... Procrustes: rmse 0.0001065687  max resid 0.0001775721 
    ## ... Similar to previous best
    ## Run 316 stress 0.08973876 
    ## Run 317 stress 0.09535448 
    ## Run 318 stress 0.09407971 
    ## Run 319 stress 0.09308967 
    ## Run 320 stress 0.1084714 
    ## Run 321 stress 0.09159083 
    ## Run 322 stress 0.09407981 
    ## Run 323 stress 0.1038046 
    ## Run 324 stress 0.08440254 
    ## ... Procrustes: rmse 0.0001343433  max resid 0.0002581945 
    ## ... Similar to previous best
    ## Run 325 stress 0.09286097 
    ## Run 326 stress 0.09408004 
    ## Run 327 stress 0.09969969 
    ## Run 328 stress 0.09465927 
    ## Run 329 stress 0.0937427 
    ## Run 330 stress 0.09168942 
    ## Run 331 stress 0.08773484 
    ## Run 332 stress 0.09286087 
    ## Run 333 stress 0.09030403 
    ## Run 334 stress 0.09416146 
    ## Run 335 stress 0.2492117 
    ## Run 336 stress 0.09712951 
    ## Run 337 stress 0.08503499 
    ## Run 338 stress 0.09030409 
    ## Run 339 stress 0.08973883 
    ## Run 340 stress 0.08503473 
    ## Run 341 stress 0.09145326 
    ## Run 342 stress 0.103804 
    ## Run 343 stress 0.08440257 
    ## ... Procrustes: rmse 0.0001201226  max resid 0.0002238631 
    ## ... Similar to previous best
    ## Run 344 stress 0.08440257 
    ## ... Procrustes: rmse 0.0001105648  max resid 0.0002240815 
    ## ... Similar to previous best
    ## Run 345 stress 0.0850346 
    ## Run 346 stress 0.0915909 
    ## Run 347 stress 0.09030408 
    ## Run 348 stress 0.09030415 
    ## Run 349 stress 0.09030398 
    ## Run 350 stress 0.09321792 
    ## Run 351 stress 0.09416138 
    ## Run 352 stress 0.08503498 
    ## Run 353 stress 0.09908188 
    ## Run 354 stress 0.1052855 
    ## Run 355 stress 0.09030403 
    ## Run 356 stress 0.09407976 
    ## Run 357 stress 0.09308956 
    ## Run 358 stress 0.09145321 
    ## Run 359 stress 0.09286084 
    ## Run 360 stress 0.09308949 
    ## Run 361 stress 0.09030394 
    ## Run 362 stress 0.09030409 
    ## Run 363 stress 0.0930895 
    ## Run 364 stress 0.09145323 
    ## Run 365 stress 0.09408005 
    ## Run 366 stress 0.09969963 
    ## Run 367 stress 0.08440269 
    ## ... Procrustes: rmse 0.0002194361  max resid 0.0003766358 
    ## ... Similar to previous best
    ## Run 368 stress 0.09321829 
    ## Run 369 stress 0.08503519 
    ## Run 370 stress 0.09969983 
    ## Run 371 stress 0.09407995 
    ## Run 372 stress 0.09268365 
    ## Run 373 stress 0.1026216 
    ## Run 374 stress 0.09407999 
    ## Run 375 stress 0.09374224 
    ## Run 376 stress 0.08503637 
    ## Run 377 stress 0.09286084 
    ## Run 378 stress 0.08773484 
    ## Run 379 stress 0.09465905 
    ## Run 380 stress 0.09030393 
    ## Run 381 stress 0.09610759 
    ## Run 382 stress 0.08503642 
    ## Run 383 stress 0.09464435 
    ## Run 384 stress 0.09030416 
    ## Run 385 stress 0.09286086 
    ## Run 386 stress 0.0850348 
    ## Run 387 stress 0.08503517 
    ## Run 388 stress 0.09969967 
    ## Run 389 stress 0.09381844 
    ## Run 390 stress 0.09407967 
    ## Run 391 stress 0.09416135 
    ## Run 392 stress 0.0944547 
    ## Run 393 stress 0.09030394 
    ## Run 394 stress 0.09268393 
    ## Run 395 stress 0.09464409 
    ## Run 396 stress 0.09030399 
    ## Run 397 stress 0.08440271 
    ## ... Procrustes: rmse 0.0002737278  max resid 0.0005520723 
    ## ... Similar to previous best
    ## Run 398 stress 0.0915909 
    ## Run 399 stress 0.09417776 
    ## Run 400 stress 0.0937427 
    ## Run 401 stress 0.09308963 
    ## Run 402 stress 0.09400534 
    ## Run 403 stress 0.09407989 
    ## Run 404 stress 0.09465907 
    ## Run 405 stress 0.09417791 
    ## Run 406 stress 0.09168942 
    ## Run 407 stress 0.09145337 
    ## Run 408 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001035969  max resid 0.0002041154 
    ## ... Similar to previous best
    ## Run 409 stress 0.09337236 
    ## Run 410 stress 0.09407301 
    ## Run 411 stress 0.09407971 
    ## Run 412 stress 0.09168952 
    ## Run 413 stress 0.08773469 
    ## Run 414 stress 0.08973864 
    ## Run 415 stress 0.08503519 
    ## Run 416 stress 0.090304 
    ## Run 417 stress 0.09268376 
    ## Run 418 stress 0.09669864 
    ## Run 419 stress 0.09159083 
    ## Run 420 stress 0.1026217 
    ## Run 421 stress 0.1038086 
    ## Run 422 stress 0.08973866 
    ## Run 423 stress 0.09268362 
    ## Run 424 stress 0.08973866 
    ## Run 425 stress 0.08503523 
    ## Run 426 stress 0.0933727 
    ## Run 427 stress 0.09721198 
    ## Run 428 stress 0.08503466 
    ## Run 429 stress 0.08440253 
    ## ... Procrustes: rmse 6.34253e-05  max resid 0.0001207106 
    ## ... Similar to previous best
    ## Run 430 stress 0.09539195 
    ## Run 431 stress 0.09268412 
    ## Run 432 stress 0.09286084 
    ## Run 433 stress 0.09159088 
    ## Run 434 stress 0.09374277 
    ## Run 435 stress 0.09308947 
    ## Run 436 stress 0.09030406 
    ## Run 437 stress 0.09407985 
    ## Run 438 stress 0.09374248 
    ## Run 439 stress 0.09030395 
    ## Run 440 stress 0.2919882 
    ## Run 441 stress 0.09145324 
    ## Run 442 stress 0.08503593 
    ## Run 443 stress 0.0926837 
    ## Run 444 stress 0.09407974 
    ## Run 445 stress 0.09407975 
    ## Run 446 stress 0.09337251 
    ## Run 447 stress 0.09337236 
    ## Run 448 stress 0.09464458 
    ## Run 449 stress 0.08503494 
    ## Run 450 stress 0.08503486 
    ## Run 451 stress 0.08503579 
    ## Run 452 stress 0.09407966 
    ## Run 453 stress 0.09268347 
    ## Run 454 stress 0.08503585 
    ## Run 455 stress 0.09308995 
    ## Run 456 stress 0.09407979 
    ## Run 457 stress 0.09610773 
    ## Run 458 stress 0.09407978 
    ## Run 459 stress 0.09268375 
    ## Run 460 stress 0.0946591 
    ## Run 461 stress 0.09760814 
    ## Run 462 stress 0.09308968 
    ## Run 463 stress 0.09969968 
    ## Run 464 stress 0.09590651 
    ## Run 465 stress 0.09374278 
    ## Run 466 stress 0.08773479 
    ## Run 467 stress 0.09374253 
    ## Run 468 stress 0.08440254 
    ## ... Procrustes: rmse 8.088473e-05  max resid 0.0001720876 
    ## ... Similar to previous best
    ## Run 469 stress 0.09268393 
    ## Run 470 stress 0.08503462 
    ## Run 471 stress 0.09590686 
    ## Run 472 stress 0.09969985 
    ## Run 473 stress 0.09416131 
    ## Run 474 stress 0.09168931 
    ## Run 475 stress 0.09145323 
    ## Run 476 stress 0.09030398 
    ## Run 477 stress 0.09969967 
    ## Run 478 stress 0.08440263 
    ## ... Procrustes: rmse 0.0001641326  max resid 0.000301285 
    ## ... Similar to previous best
    ## Run 479 stress 0.09286104 
    ## Run 480 stress 0.09159086 
    ## Run 481 stress 0.0914533 
    ## Run 482 stress 0.09969975 
    ## Run 483 stress 0.08773469 
    ## Run 484 stress 0.09308954 
    ## Run 485 stress 0.09030416 
    ## Run 486 stress 0.08773466 
    ## Run 487 stress 0.0850347 
    ## Run 488 stress 0.08440267 
    ## ... Procrustes: rmse 0.0001968271  max resid 0.0003491117 
    ## ... Similar to previous best
    ## Run 489 stress 0.0971295 
    ## Run 490 stress 0.09145339 
    ## Run 491 stress 0.09465905 
    ## Run 492 stress 0.09465919 
    ## Run 493 stress 0.08503515 
    ## Run 494 stress 0.09030411 
    ## Run 495 stress 0.09145324 
    ## Run 496 stress 0.08440251 
    ## ... Procrustes: rmse 1.270018e-05  max resid 2.675763e-05 
    ## ... Similar to previous best
    ## Run 497 stress 0.09969985 
    ## Run 498 stress 0.09159086 
    ## Run 499 stress 0.09321603 
    ## Run 500 stress 0.09308953 
    ## *** Best solution repeated 18 times

``` r
# Ocean sites and mixed lakes
SD_beta_OM_NMDS <- metaMDS(SD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07365783 
    ## Run 1 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001059661  max resid 0.0002516394 
    ## ... Similar to previous best
    ## Run 2 stress 0.07629244 
    ## Run 3 stress 0.08000469 
    ## Run 4 stress 0.08000469 
    ## Run 5 stress 0.08233769 
    ## Run 6 stress 0.08000472 
    ## Run 7 stress 0.07365783 
    ## ... Procrustes: rmse 2.51408e-05  max resid 5.839224e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.07629233 
    ## Run 9 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743986  max resid 0.05102104 
    ## Run 10 stress 0.0823377 
    ## Run 11 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001265569  max resid 0.0002954006 
    ## ... Similar to previous best
    ## Run 12 stress 0.08000477 
    ## Run 13 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744031  max resid 0.05104675 
    ## Run 14 stress 0.07629238 
    ## Run 15 stress 0.0762924 
    ## Run 16 stress 0.08000472 
    ## Run 17 stress 0.07378232 
    ## ... Procrustes: rmse 0.01745763  max resid 0.05112909 
    ## Run 18 stress 0.08000475 
    ## Run 19 stress 0.07365783 
    ## ... Procrustes: rmse 2.171844e-05  max resid 5.073375e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.07629236 
    ## Run 21 stress 0.07629239 
    ## Run 22 stress 0.07365783 
    ## ... Procrustes: rmse 3.909051e-06  max resid 7.492674e-06 
    ## ... Similar to previous best
    ## Run 23 stress 0.07629239 
    ## Run 24 stress 0.08233771 
    ## Run 25 stress 0.08000468 
    ## Run 26 stress 0.0737823 
    ## ... Procrustes: rmse 0.01743943  max resid 0.05104865 
    ## Run 27 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744168  max resid 0.05103065 
    ## Run 28 stress 0.07365783 
    ## ... Procrustes: rmse 5.258978e-06  max resid 8.876999e-06 
    ## ... Similar to previous best
    ## Run 29 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744538  max resid 0.05105732 
    ## Run 30 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744373  max resid 0.05103896 
    ## Run 31 stress 0.07365785 
    ## ... Procrustes: rmse 8.246127e-05  max resid 0.0001950008 
    ## ... Similar to previous best
    ## Run 32 stress 0.08000467 
    ## Run 33 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001292608  max resid 0.0003084612 
    ## ... Similar to previous best
    ## Run 34 stress 0.08000468 
    ## Run 35 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 3.012879e-06  max resid 5.695315e-06 
    ## ... Similar to previous best
    ## Run 36 stress 0.08000471 
    ## Run 37 stress 0.07629248 
    ## Run 38 stress 0.07629236 
    ## Run 39 stress 0.08000468 
    ## Run 40 stress 0.07629233 
    ## Run 41 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744682  max resid 0.05105016 
    ## Run 42 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744181  max resid 0.05105065 
    ## Run 43 stress 0.07365784 
    ## ... Procrustes: rmse 4.081675e-05  max resid 9.478353e-05 
    ## ... Similar to previous best
    ## Run 44 stress 0.07378232 
    ## ... Procrustes: rmse 0.01746361  max resid 0.05115032 
    ## Run 45 stress 0.07365783 
    ## ... Procrustes: rmse 3.531744e-06  max resid 6.947045e-06 
    ## ... Similar to previous best
    ## Run 46 stress 0.07629233 
    ## Run 47 stress 0.08000472 
    ## Run 48 stress 0.07629244 
    ## Run 49 stress 0.08288059 
    ## Run 50 stress 0.08000474 
    ## Run 51 stress 0.07365783 
    ## ... Procrustes: rmse 3.047629e-05  max resid 6.983578e-05 
    ## ... Similar to previous best
    ## Run 52 stress 0.07365783 
    ## ... Procrustes: rmse 3.576516e-06  max resid 8.464733e-06 
    ## ... Similar to previous best
    ## Run 53 stress 0.07629233 
    ## Run 54 stress 0.07365787 
    ## ... Procrustes: rmse 6.743714e-05  max resid 0.0001609392 
    ## ... Similar to previous best
    ## Run 55 stress 0.07629234 
    ## Run 56 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744523  max resid 0.05104002 
    ## Run 57 stress 0.07378231 
    ## ... Procrustes: rmse 0.01744117  max resid 0.05105261 
    ## Run 58 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001144399  max resid 0.0002706026 
    ## ... Similar to previous best
    ## Run 59 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001663416  max resid 0.0003900514 
    ## ... Similar to previous best
    ## Run 60 stress 0.07629237 
    ## Run 61 stress 0.08288039 
    ## Run 62 stress 0.07629241 
    ## Run 63 stress 0.0800047 
    ## Run 64 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744442  max resid 0.05103809 
    ## Run 65 stress 0.07732937 
    ## Run 66 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744944  max resid 0.05106247 
    ## Run 67 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744619  max resid 0.05105345 
    ## Run 68 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744381  max resid 0.05103624 
    ## Run 69 stress 0.07365787 
    ## ... Procrustes: rmse 9.757522e-05  max resid 0.0002316221 
    ## ... Similar to previous best
    ## Run 70 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744291  max resid 0.05105475 
    ## Run 71 stress 0.07629235 
    ## Run 72 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001352273  max resid 0.0003183981 
    ## ... Similar to previous best
    ## Run 73 stress 0.0773293 
    ## Run 74 stress 0.07629242 
    ## Run 75 stress 0.07629248 
    ## Run 76 stress 0.07629241 
    ## Run 77 stress 0.0800047 
    ## Run 78 stress 0.0800047 
    ## Run 79 stress 0.07629245 
    ## Run 80 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001357864  max resid 0.0003135066 
    ## ... Similar to previous best
    ## Run 81 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744628  max resid 0.05104571 
    ## Run 82 stress 0.07365785 
    ## ... Procrustes: rmse 6.146236e-05  max resid 0.0001467427 
    ## ... Similar to previous best
    ## Run 83 stress 0.08233758 
    ## Run 84 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745032  max resid 0.05107729 
    ## Run 85 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744682  max resid 0.05106143 
    ## Run 86 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744832  max resid 0.05106283 
    ## Run 87 stress 0.07732942 
    ## Run 88 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001068975  max resid 0.0002525051 
    ## ... Similar to previous best
    ## Run 89 stress 0.07629237 
    ## Run 90 stress 0.08000476 
    ## Run 91 stress 0.07732934 
    ## Run 92 stress 0.07732928 
    ## Run 93 stress 0.08233753 
    ## Run 94 stress 0.07365784 
    ## ... Procrustes: rmse 6.004324e-05  max resid 0.0001378156 
    ## ... Similar to previous best
    ## Run 95 stress 0.0800047 
    ## Run 96 stress 0.07365784 
    ## ... Procrustes: rmse 4.713426e-05  max resid 0.0001063322 
    ## ... Similar to previous best
    ## Run 97 stress 0.07732944 
    ## Run 98 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744395  max resid 0.05102232 
    ## Run 99 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744407  max resid 0.05103243 
    ## Run 100 stress 0.07732938 
    ## Run 101 stress 0.07629235 
    ## Run 102 stress 0.07732928 
    ## Run 103 stress 0.07365784 
    ## ... Procrustes: rmse 3.974732e-05  max resid 9.219659e-05 
    ## ... Similar to previous best
    ## Run 104 stress 0.07629248 
    ## Run 105 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001035597  max resid 0.0002422819 
    ## ... Similar to previous best
    ## Run 106 stress 0.08000468 
    ## Run 107 stress 0.07365784 
    ## ... Procrustes: rmse 5.874444e-05  max resid 0.0001346346 
    ## ... Similar to previous best
    ## Run 108 stress 0.07629237 
    ## Run 109 stress 0.07629244 
    ## Run 110 stress 0.0762924 
    ## Run 111 stress 0.08000469 
    ## Run 112 stress 0.07732937 
    ## Run 113 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174442  max resid 0.05103901 
    ## Run 114 stress 0.07629235 
    ## Run 115 stress 0.08000468 
    ## Run 116 stress 0.07629233 
    ## Run 117 stress 0.07365784 
    ## ... Procrustes: rmse 4.411031e-05  max resid 0.0001055511 
    ## ... Similar to previous best
    ## Run 118 stress 0.07365784 
    ## ... Procrustes: rmse 4.368538e-05  max resid 0.0001021487 
    ## ... Similar to previous best
    ## Run 119 stress 0.2747113 
    ## Run 120 stress 0.07365784 
    ## ... Procrustes: rmse 4.553064e-05  max resid 0.0001077683 
    ## ... Similar to previous best
    ## Run 121 stress 0.07629233 
    ## Run 122 stress 0.07732946 
    ## Run 123 stress 0.08233757 
    ## Run 124 stress 0.07365783 
    ## ... Procrustes: rmse 2.047067e-05  max resid 4.401463e-05 
    ## ... Similar to previous best
    ## Run 125 stress 0.07629238 
    ## Run 126 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744243  max resid 0.05102395 
    ## Run 127 stress 0.07629248 
    ## Run 128 stress 0.07365794 
    ## ... Procrustes: rmse 0.0001708755  max resid 0.0004070231 
    ## ... Similar to previous best
    ## Run 129 stress 0.08000471 
    ## Run 130 stress 0.0773294 
    ## Run 131 stress 0.07365784 
    ## ... Procrustes: rmse 4.126089e-05  max resid 9.774655e-05 
    ## ... Similar to previous best
    ## Run 132 stress 0.07732932 
    ## Run 133 stress 0.07365784 
    ## ... Procrustes: rmse 2.165325e-05  max resid 4.38256e-05 
    ## ... Similar to previous best
    ## Run 134 stress 0.07365787 
    ## ... Procrustes: rmse 0.000106009  max resid 0.000250888 
    ## ... Similar to previous best
    ## Run 135 stress 0.07365786 
    ## ... Procrustes: rmse 0.0001007647  max resid 0.0002353725 
    ## ... Similar to previous best
    ## Run 136 stress 0.07365784 
    ## ... Procrustes: rmse 3.094235e-05  max resid 6.911448e-05 
    ## ... Similar to previous best
    ## Run 137 stress 0.07365786 
    ## ... Procrustes: rmse 7.128308e-05  max resid 0.0001696827 
    ## ... Similar to previous best
    ## Run 138 stress 0.07365784 
    ## ... Procrustes: rmse 4.838671e-05  max resid 0.0001112411 
    ## ... Similar to previous best
    ## Run 139 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001578403  max resid 0.0003740526 
    ## ... Similar to previous best
    ## Run 140 stress 0.07732939 
    ## Run 141 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001383071  max resid 0.0003285375 
    ## ... Similar to previous best
    ## Run 142 stress 0.2656487 
    ## Run 143 stress 0.07629235 
    ## Run 144 stress 0.0828804 
    ## Run 145 stress 0.08233753 
    ## Run 146 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744486  max resid 0.05105081 
    ## Run 147 stress 0.08000468 
    ## Run 148 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001184039  max resid 0.0002794743 
    ## ... Similar to previous best
    ## Run 149 stress 0.08288043 
    ## Run 150 stress 0.07629234 
    ## Run 151 stress 0.07378229 
    ## ... Procrustes: rmse 0.01746436  max resid 0.05111002 
    ## Run 152 stress 0.07365783 
    ## ... Procrustes: rmse 7.131225e-06  max resid 1.574766e-05 
    ## ... Similar to previous best
    ## Run 153 stress 0.08288062 
    ## Run 154 stress 0.0800047 
    ## Run 155 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744573  max resid 0.05104345 
    ## Run 156 stress 0.08000471 
    ## Run 157 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744641  max resid 0.05105341 
    ## Run 158 stress 0.07365783 
    ## ... Procrustes: rmse 1.970623e-05  max resid 4.644643e-05 
    ## ... Similar to previous best
    ## Run 159 stress 0.07365785 
    ## ... Procrustes: rmse 6.683989e-05  max resid 0.0001585648 
    ## ... Similar to previous best
    ## Run 160 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001434028  max resid 0.0003377842 
    ## ... Similar to previous best
    ## Run 161 stress 0.07378228 
    ## ... Procrustes: rmse 0.01743389  max resid 0.05097028 
    ## Run 162 stress 0.07629234 
    ## Run 163 stress 0.08000467 
    ## Run 164 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001308219  max resid 0.0003071015 
    ## ... Similar to previous best
    ## Run 165 stress 0.07365785 
    ## ... Procrustes: rmse 6.844312e-05  max resid 0.0001633324 
    ## ... Similar to previous best
    ## Run 166 stress 0.08000481 
    ## Run 167 stress 0.07629233 
    ## Run 168 stress 0.07629237 
    ## Run 169 stress 0.08000474 
    ## Run 170 stress 0.07378231 
    ## ... Procrustes: rmse 0.01743947  max resid 0.05105052 
    ## Run 171 stress 0.07365788 
    ## ... Procrustes: rmse 9.474237e-05  max resid 0.0002148153 
    ## ... Similar to previous best
    ## Run 172 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744526  max resid 0.05106033 
    ## Run 173 stress 0.07365784 
    ## ... Procrustes: rmse 2.311302e-05  max resid 5.284843e-05 
    ## ... Similar to previous best
    ## Run 174 stress 0.07629238 
    ## Run 175 stress 0.08233766 
    ## Run 176 stress 0.07629247 
    ## Run 177 stress 0.07629238 
    ## Run 178 stress 0.07732943 
    ## Run 179 stress 0.07378231 
    ## ... Procrustes: rmse 0.01743927  max resid 0.05105005 
    ## Run 180 stress 0.07629243 
    ## Run 181 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744678  max resid 0.05104805 
    ## Run 182 stress 0.0737823 
    ## ... Procrustes: rmse 0.0174387  max resid 0.05104763 
    ## Run 183 stress 0.07629239 
    ## Run 184 stress 0.07629237 
    ## Run 185 stress 0.07629237 
    ## Run 186 stress 0.07378231 
    ## ... Procrustes: rmse 0.01745863  max resid 0.05107009 
    ## Run 187 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744507  max resid 0.05105013 
    ## Run 188 stress 0.08000468 
    ## Run 189 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744031  max resid 0.05105132 
    ## Run 190 stress 0.07629233 
    ## Run 191 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744224  max resid 0.0510529 
    ## Run 192 stress 0.07732937 
    ## Run 193 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744917  max resid 0.05107316 
    ## Run 194 stress 0.08000468 
    ## Run 195 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744918  max resid 0.05106524 
    ## Run 196 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744241  max resid 0.05102235 
    ## Run 197 stress 0.07629238 
    ## Run 198 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 3.248758e-06  max resid 7.260008e-06 
    ## ... Similar to previous best
    ## Run 199 stress 0.07365786 
    ## ... Procrustes: rmse 9.432613e-05  max resid 0.0002197807 
    ## ... Similar to previous best
    ## Run 200 stress 0.08000468 
    ## Run 201 stress 0.2409736 
    ## Run 202 stress 0.07365788 
    ## ... Procrustes: rmse 7.216385e-05  max resid 0.0001714135 
    ## ... Similar to previous best
    ## Run 203 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744977  max resid 0.0510452 
    ## Run 204 stress 0.07629234 
    ## Run 205 stress 0.07365797 
    ## ... Procrustes: rmse 4.569918e-05  max resid 9.329546e-05 
    ## ... Similar to previous best
    ## Run 206 stress 0.07629233 
    ## Run 207 stress 0.07365786 
    ## ... Procrustes: rmse 7.518425e-05  max resid 0.00017939 
    ## ... Similar to previous best
    ## Run 208 stress 0.08233761 
    ## Run 209 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744327  max resid 0.05103225 
    ## Run 210 stress 0.07365786 
    ## ... Procrustes: rmse 8.875133e-05  max resid 0.0002099081 
    ## ... Similar to previous best
    ## Run 211 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744701  max resid 0.05105064 
    ## Run 212 stress 0.08233755 
    ## Run 213 stress 0.08233753 
    ## Run 214 stress 0.07629233 
    ## Run 215 stress 0.07629243 
    ## Run 216 stress 0.08000468 
    ## Run 217 stress 0.08000479 
    ## Run 218 stress 0.07365784 
    ## ... Procrustes: rmse 4.709221e-05  max resid 0.0001111433 
    ## ... Similar to previous best
    ## Run 219 stress 0.07365783 
    ## ... Procrustes: rmse 4.248185e-06  max resid 9.16473e-06 
    ## ... Similar to previous best
    ## Run 220 stress 0.07629233 
    ## Run 221 stress 0.07629234 
    ## Run 222 stress 0.07629234 
    ## Run 223 stress 0.0800047 
    ## Run 224 stress 0.07732934 
    ## Run 225 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001144529  max resid 0.0002722063 
    ## ... Similar to previous best
    ## Run 226 stress 0.07365783 
    ## ... Procrustes: rmse 2.569094e-05  max resid 6.018069e-05 
    ## ... Similar to previous best
    ## Run 227 stress 0.08233771 
    ## Run 228 stress 0.3578484 
    ## Run 229 stress 0.0737823 
    ## ... Procrustes: rmse 0.017444  max resid 0.05102802 
    ## Run 230 stress 0.07629233 
    ## Run 231 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745346  max resid 0.05107445 
    ## Run 232 stress 0.07378234 
    ## ... Procrustes: rmse 0.01746495  max resid 0.05115568 
    ## Run 233 stress 0.07629235 
    ## Run 234 stress 0.07629234 
    ## Run 235 stress 0.08000476 
    ## Run 236 stress 0.0800047 
    ## Run 237 stress 0.08000477 
    ## Run 238 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744522  max resid 0.05106712 
    ## Run 239 stress 0.07365784 
    ## ... Procrustes: rmse 6.027343e-05  max resid 0.0001403561 
    ## ... Similar to previous best
    ## Run 240 stress 0.07365794 
    ## ... Procrustes: rmse 0.000168752  max resid 0.0004019637 
    ## ... Similar to previous best
    ## Run 241 stress 0.07365783 
    ## ... Procrustes: rmse 1.522701e-05  max resid 3.438008e-05 
    ## ... Similar to previous best
    ## Run 242 stress 0.07365785 
    ## ... Procrustes: rmse 7.672323e-05  max resid 0.0001816701 
    ## ... Similar to previous best
    ## Run 243 stress 0.08000484 
    ## Run 244 stress 0.08000468 
    ## Run 245 stress 0.07629234 
    ## Run 246 stress 0.07629241 
    ## Run 247 stress 0.07365786 
    ## ... Procrustes: rmse 8.513253e-05  max resid 0.0002014984 
    ## ... Similar to previous best
    ## Run 248 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001406368  max resid 0.0003276416 
    ## ... Similar to previous best
    ## Run 249 stress 0.3578481 
    ## Run 250 stress 0.07629236 
    ## Run 251 stress 0.08000479 
    ## Run 252 stress 0.07629234 
    ## Run 253 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744563  max resid 0.05103845 
    ## Run 254 stress 0.07365783 
    ## ... Procrustes: rmse 3.189497e-05  max resid 7.56538e-05 
    ## ... Similar to previous best
    ## Run 255 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174432  max resid 0.05105458 
    ## Run 256 stress 0.07629233 
    ## Run 257 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744771  max resid 0.05105015 
    ## Run 258 stress 0.07365784 
    ## ... Procrustes: rmse 5.023475e-05  max resid 0.00011302 
    ## ... Similar to previous best
    ## Run 259 stress 0.07365785 
    ## ... Procrustes: rmse 8.105781e-05  max resid 0.000192331 
    ## ... Similar to previous best
    ## Run 260 stress 0.07732932 
    ## Run 261 stress 0.07629238 
    ## Run 262 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744399  max resid 0.05103358 
    ## Run 263 stress 0.07365783 
    ## ... Procrustes: rmse 1.424021e-05  max resid 3.099906e-05 
    ## ... Similar to previous best
    ## Run 264 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174439  max resid 0.05102519 
    ## Run 265 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744674  max resid 0.05105552 
    ## Run 266 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174466  max resid 0.05104998 
    ## Run 267 stress 0.07629233 
    ## Run 268 stress 0.07629248 
    ## Run 269 stress 0.07378231 
    ## ... Procrustes: rmse 0.01746218  max resid 0.05108179 
    ## Run 270 stress 0.08000474 
    ## Run 271 stress 0.07732932 
    ## Run 272 stress 0.07365784 
    ## ... Procrustes: rmse 5.525304e-05  max resid 0.0001301306 
    ## ... Similar to previous best
    ## Run 273 stress 0.07365786 
    ## ... Procrustes: rmse 8.4203e-05  max resid 0.0001996609 
    ## ... Similar to previous best
    ## Run 274 stress 0.07732932 
    ## Run 275 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744205  max resid 0.05102888 
    ## Run 276 stress 0.0737823 
    ## ... Procrustes: rmse 0.0174417  max resid 0.05105478 
    ## Run 277 stress 0.07365785 
    ## ... Procrustes: rmse 4.750152e-05  max resid 0.0001042159 
    ## ... Similar to previous best
    ## Run 278 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744592  max resid 0.05105768 
    ## Run 279 stress 0.07365784 
    ## ... Procrustes: rmse 5.336549e-05  max resid 0.0001265127 
    ## ... Similar to previous best
    ## Run 280 stress 0.07629233 
    ## Run 281 stress 0.0762924 
    ## Run 282 stress 0.07629245 
    ## Run 283 stress 0.07629234 
    ## Run 284 stress 0.07629238 
    ## Run 285 stress 0.07629234 
    ## Run 286 stress 0.07365784 
    ## ... Procrustes: rmse 5.630892e-05  max resid 0.0001323105 
    ## ... Similar to previous best
    ## Run 287 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001010376  max resid 0.0002416088 
    ## ... Similar to previous best
    ## Run 288 stress 0.08000471 
    ## Run 289 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001019104  max resid 0.000239711 
    ## ... Similar to previous best
    ## Run 290 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744371  max resid 0.05102387 
    ## Run 291 stress 0.07629233 
    ## Run 292 stress 0.07365783 
    ## ... Procrustes: rmse 1.181794e-05  max resid 2.763435e-05 
    ## ... Similar to previous best
    ## Run 293 stress 0.0762924 
    ## Run 294 stress 0.08000475 
    ## Run 295 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745479  max resid 0.0510643 
    ## Run 296 stress 0.357848 
    ## Run 297 stress 0.07732936 
    ## Run 298 stress 0.07629233 
    ## Run 299 stress 0.07629245 
    ## Run 300 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745259  max resid 0.05105177 
    ## Run 301 stress 0.08000471 
    ## Run 302 stress 0.07629237 
    ## Run 303 stress 0.08000479 
    ## Run 304 stress 0.0737823 
    ## ... Procrustes: rmse 0.01745975  max resid 0.05107589 
    ## Run 305 stress 0.08000476 
    ## Run 306 stress 0.0773294 
    ## Run 307 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744735  max resid 0.05105389 
    ## Run 308 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744522  max resid 0.05104327 
    ## Run 309 stress 0.07732926 
    ## Run 310 stress 0.07629233 
    ## Run 311 stress 0.07629233 
    ## Run 312 stress 0.2997109 
    ## Run 313 stress 0.08000469 
    ## Run 314 stress 0.0800047 
    ## Run 315 stress 0.07365785 
    ## ... Procrustes: rmse 6.86198e-05  max resid 0.0001643655 
    ## ... Similar to previous best
    ## Run 316 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744589  max resid 0.05104397 
    ## Run 317 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744883  max resid 0.05107091 
    ## Run 318 stress 0.07732934 
    ## Run 319 stress 0.08000472 
    ## Run 320 stress 0.07732933 
    ## Run 321 stress 0.07629245 
    ## Run 322 stress 0.08233764 
    ## Run 323 stress 0.08000468 
    ## Run 324 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001340399  max resid 0.0003144255 
    ## ... Similar to previous best
    ## Run 325 stress 0.2409737 
    ## Run 326 stress 0.08233768 
    ## Run 327 stress 0.07732929 
    ## Run 328 stress 0.08288045 
    ## Run 329 stress 0.0737823 
    ## ... Procrustes: rmse 0.0174411  max resid 0.05105402 
    ## Run 330 stress 0.07365783 
    ## ... Procrustes: rmse 1.58553e-05  max resid 3.059609e-05 
    ## ... Similar to previous best
    ## Run 331 stress 0.07629233 
    ## Run 332 stress 0.0773293 
    ## Run 333 stress 0.08233774 
    ## Run 334 stress 0.07629238 
    ## Run 335 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744345  max resid 0.05103489 
    ## Run 336 stress 0.07629238 
    ## Run 337 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745  max resid 0.05107385 
    ## Run 338 stress 0.08000467 
    ## Run 339 stress 0.07365784 
    ## ... Procrustes: rmse 4.185136e-05  max resid 9.883515e-05 
    ## ... Similar to previous best
    ## Run 340 stress 0.08000468 
    ## Run 341 stress 0.07378231 
    ## ... Procrustes: rmse 0.01745253  max resid 0.05110153 
    ## Run 342 stress 0.08000468 
    ## Run 343 stress 0.08000467 
    ## Run 344 stress 0.07365783 
    ## ... Procrustes: rmse 1.693678e-06  max resid 2.790264e-06 
    ## ... Similar to previous best
    ## Run 345 stress 0.08000481 
    ## Run 346 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744576  max resid 0.05104695 
    ## Run 347 stress 0.07365783 
    ## ... Procrustes: rmse 2.795026e-05  max resid 6.663771e-05 
    ## ... Similar to previous best
    ## Run 348 stress 0.07629237 
    ## Run 349 stress 0.07629236 
    ## Run 350 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744806  max resid 0.05105835 
    ## Run 351 stress 0.08000475 
    ## Run 352 stress 0.08233754 
    ## Run 353 stress 0.08000473 
    ## Run 354 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744155  max resid 0.05100353 
    ## Run 355 stress 0.07629244 
    ## Run 356 stress 0.08000478 
    ## Run 357 stress 0.07365784 
    ## ... Procrustes: rmse 4.661719e-05  max resid 0.0001104032 
    ## ... Similar to previous best
    ## Run 358 stress 0.07629246 
    ## Run 359 stress 0.07629248 
    ## Run 360 stress 0.07629234 
    ## Run 361 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174427  max resid 0.05103004 
    ## Run 362 stress 0.07629235 
    ## Run 363 stress 0.07365789 
    ## ... Procrustes: rmse 0.000131782  max resid 0.0003112133 
    ## ... Similar to previous best
    ## Run 364 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744702  max resid 0.05105083 
    ## Run 365 stress 0.07365783 
    ## ... Procrustes: rmse 3.035678e-06  max resid 7.257261e-06 
    ## ... Similar to previous best
    ## Run 366 stress 0.07629237 
    ## Run 367 stress 0.08288049 
    ## Run 368 stress 0.08000475 
    ## Run 369 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744687  max resid 0.05105198 
    ## Run 370 stress 0.07629246 
    ## Run 371 stress 0.07629233 
    ## Run 372 stress 0.07732925 
    ## Run 373 stress 0.07365783 
    ## ... Procrustes: rmse 4.763408e-06  max resid 9.546079e-06 
    ## ... Similar to previous best
    ## Run 374 stress 0.07629243 
    ## Run 375 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744491  max resid 0.051033 
    ## Run 376 stress 0.07629233 
    ## Run 377 stress 0.08288047 
    ## Run 378 stress 0.08000467 
    ## Run 379 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744411  max resid 0.05105895 
    ## Run 380 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744618  max resid 0.0510577 
    ## Run 381 stress 0.07629242 
    ## Run 382 stress 0.07732926 
    ## Run 383 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744593  max resid 0.05104258 
    ## Run 384 stress 0.07365783 
    ## ... Procrustes: rmse 1.732933e-05  max resid 4.0508e-05 
    ## ... Similar to previous best
    ## Run 385 stress 0.07732928 
    ## Run 386 stress 0.07365783 
    ## ... Procrustes: rmse 1.690224e-05  max resid 4.071393e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.07629244 
    ## Run 388 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001382407  max resid 0.0003283424 
    ## ... Similar to previous best
    ## Run 389 stress 0.07629235 
    ## Run 390 stress 0.07365785 
    ## ... Procrustes: rmse 7.975273e-05  max resid 0.0001878984 
    ## ... Similar to previous best
    ## Run 391 stress 0.08000468 
    ## Run 392 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743889  max resid 0.05101536 
    ## Run 393 stress 0.07365783 
    ## ... Procrustes: rmse 3.146602e-05  max resid 7.388857e-05 
    ## ... Similar to previous best
    ## Run 394 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001051249  max resid 0.0002515 
    ## ... Similar to previous best
    ## Run 395 stress 0.07365783 
    ## ... Procrustes: rmse 3.482512e-05  max resid 8.249101e-05 
    ## ... Similar to previous best
    ## Run 396 stress 0.07365784 
    ## ... Procrustes: rmse 1.578104e-05  max resid 2.800895e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.07365783 
    ## ... Procrustes: rmse 1.761324e-05  max resid 4.136605e-05 
    ## ... Similar to previous best
    ## Run 398 stress 0.08000474 
    ## Run 399 stress 0.08000469 
    ## Run 400 stress 0.07629235 
    ## Run 401 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744003  max resid 0.05101898 
    ## Run 402 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174467  max resid 0.05104878 
    ## Run 403 stress 0.08000469 
    ## Run 404 stress 0.07365784 
    ## ... Procrustes: rmse 4.545723e-05  max resid 0.0001080907 
    ## ... Similar to previous best
    ## Run 405 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744988  max resid 0.05106899 
    ## Run 406 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744477  max resid 0.05105659 
    ## Run 407 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174452  max resid 0.0510397 
    ## Run 408 stress 0.07732929 
    ## Run 409 stress 0.07732939 
    ## Run 410 stress 0.08000468 
    ## Run 411 stress 0.07365783 
    ## ... Procrustes: rmse 2.297229e-05  max resid 5.452769e-05 
    ## ... Similar to previous best
    ## Run 412 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744653  max resid 0.05105959 
    ## Run 413 stress 0.07365784 
    ## ... Procrustes: rmse 5.479311e-05  max resid 0.0001262082 
    ## ... Similar to previous best
    ## Run 414 stress 0.0737823 
    ## ... Procrustes: rmse 0.01743443  max resid 0.05096179 
    ## Run 415 stress 0.07365784 
    ## ... Procrustes: rmse 2.189832e-05  max resid 4.392861e-05 
    ## ... Similar to previous best
    ## Run 416 stress 0.08288043 
    ## Run 417 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001046546  max resid 0.0002483135 
    ## ... Similar to previous best
    ## Run 418 stress 0.07365785 
    ## ... Procrustes: rmse 7.284112e-05  max resid 0.0001700771 
    ## ... Similar to previous best
    ## Run 419 stress 0.0762924 
    ## Run 420 stress 0.07365784 
    ## ... Procrustes: rmse 1.146803e-05  max resid 2.069645e-05 
    ## ... Similar to previous best
    ## Run 421 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744229  max resid 0.05103122 
    ## Run 422 stress 0.07629244 
    ## Run 423 stress 0.3578481 
    ## Run 424 stress 0.07365783 
    ## ... Procrustes: rmse 3.130043e-05  max resid 7.331678e-05 
    ## ... Similar to previous best
    ## Run 425 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744213  max resid 0.0510303 
    ## Run 426 stress 0.07629242 
    ## Run 427 stress 0.0800047 
    ## Run 428 stress 0.07365784 
    ## ... Procrustes: rmse 2.370323e-05  max resid 5.479085e-05 
    ## ... Similar to previous best
    ## Run 429 stress 0.07378232 
    ## ... Procrustes: rmse 0.01744804  max resid 0.05101938 
    ## Run 430 stress 0.07365784 
    ## ... Procrustes: rmse 5.589694e-05  max resid 0.0001327091 
    ## ... Similar to previous best
    ## Run 431 stress 0.07629237 
    ## Run 432 stress 0.07629234 
    ## Run 433 stress 0.07629234 
    ## Run 434 stress 0.0800047 
    ## Run 435 stress 0.07629236 
    ## Run 436 stress 0.07629239 
    ## Run 437 stress 0.0737823 
    ## ... Procrustes: rmse 0.01745367  max resid 0.0511002 
    ## Run 438 stress 0.3265737 
    ## Run 439 stress 0.2574761 
    ## Run 440 stress 0.07365784 
    ## ... Procrustes: rmse 5.518081e-05  max resid 0.0001300684 
    ## ... Similar to previous best
    ## Run 441 stress 0.08288038 
    ## Run 442 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744849  max resid 0.05106915 
    ## Run 443 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744609  max resid 0.05105644 
    ## Run 444 stress 0.08000475 
    ## Run 445 stress 0.07629236 
    ## Run 446 stress 0.08233751 
    ## Run 447 stress 0.08000469 
    ## Run 448 stress 0.07629243 
    ## Run 449 stress 0.07365785 
    ## ... Procrustes: rmse 7.625687e-05  max resid 0.0001862768 
    ## ... Similar to previous best
    ## Run 450 stress 0.08000467 
    ## Run 451 stress 0.07629233 
    ## Run 452 stress 0.07365783 
    ## ... Procrustes: rmse 8.145791e-06  max resid 1.917351e-05 
    ## ... Similar to previous best
    ## Run 453 stress 0.07378231 
    ## ... Procrustes: rmse 0.01744013  max resid 0.05105221 
    ## Run 454 stress 0.07365784 
    ## ... Procrustes: rmse 3.481613e-05  max resid 8.200455e-05 
    ## ... Similar to previous best
    ## Run 455 stress 0.08000471 
    ## Run 456 stress 0.08233756 
    ## Run 457 stress 0.07365784 
    ## ... Procrustes: rmse 5.352251e-05  max resid 0.000122097 
    ## ... Similar to previous best
    ## Run 458 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001223058  max resid 0.0002833633 
    ## ... Similar to previous best
    ## Run 459 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744492  max resid 0.05103429 
    ## Run 460 stress 0.07365785 
    ## ... Procrustes: rmse 6.549915e-05  max resid 0.0001555939 
    ## ... Similar to previous best
    ## Run 461 stress 0.07365783 
    ## ... Procrustes: rmse 5.90826e-06  max resid 1.404803e-05 
    ## ... Similar to previous best
    ## Run 462 stress 0.08000467 
    ## Run 463 stress 0.07629235 
    ## Run 464 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001035253  max resid 0.0002450952 
    ## ... Similar to previous best
    ## Run 465 stress 0.08000472 
    ## Run 466 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744951  max resid 0.05106104 
    ## Run 467 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744244  max resid 0.05101724 
    ## Run 468 stress 0.07629233 
    ## Run 469 stress 0.08000468 
    ## Run 470 stress 0.07365784 
    ## ... Procrustes: rmse 6.206804e-05  max resid 0.000145955 
    ## ... Similar to previous best
    ## Run 471 stress 0.07629234 
    ## Run 472 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001337156  max resid 0.0003164747 
    ## ... Similar to previous best
    ## Run 473 stress 0.07378231 
    ## ... Procrustes: rmse 0.0174432  max resid 0.0510022 
    ## Run 474 stress 0.08000474 
    ## Run 475 stress 0.07629234 
    ## Run 476 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001538482  max resid 0.000361276 
    ## ... Similar to previous best
    ## Run 477 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001526241  max resid 0.0003589469 
    ## ... Similar to previous best
    ## Run 478 stress 0.07365784 
    ## ... Procrustes: rmse 4.266837e-05  max resid 9.808326e-05 
    ## ... Similar to previous best
    ## Run 479 stress 0.07629244 
    ## Run 480 stress 0.07365783 
    ## ... Procrustes: rmse 1.398765e-05  max resid 3.502599e-05 
    ## ... Similar to previous best
    ## Run 481 stress 0.0762924 
    ## Run 482 stress 0.07365783 
    ## ... Procrustes: rmse 2.4157e-05  max resid 5.681168e-05 
    ## ... Similar to previous best
    ## Run 483 stress 0.07365784 
    ## ... Procrustes: rmse 3.790323e-05  max resid 9.076281e-05 
    ## ... Similar to previous best
    ## Run 484 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744541  max resid 0.05104043 
    ## Run 485 stress 0.07365783 
    ## ... Procrustes: rmse 2.977533e-06  max resid 6.712657e-06 
    ## ... Similar to previous best
    ## Run 486 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744327  max resid 0.05101643 
    ## Run 487 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001424895  max resid 0.0003407673 
    ## ... Similar to previous best
    ## Run 488 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174464  max resid 0.05105556 
    ## Run 489 stress 0.07629238 
    ## Run 490 stress 0.07378226 
    ## ... Procrustes: rmse 0.01746015  max resid 0.05110358 
    ## Run 491 stress 0.08288046 
    ## Run 492 stress 0.07629238 
    ## Run 493 stress 0.07629236 
    ## Run 494 stress 0.07365785 
    ## ... Procrustes: rmse 6.861694e-05  max resid 0.000159411 
    ## ... Similar to previous best
    ## Run 495 stress 0.08000472 
    ## Run 496 stress 0.07732936 
    ## Run 497 stress 0.07629249 
    ## Run 498 stress 0.07365784 
    ## ... Procrustes: rmse 3.702696e-05  max resid 8.841523e-05 
    ## ... Similar to previous best
    ## Run 499 stress 0.07378231 
    ## ... Procrustes: rmse 0.01744105  max resid 0.05105114 
    ## Run 500 stress 0.07732931 
    ## *** Best solution repeated 78 times

``` r
# Stratified lakes and ocean sites
SD_beta_SO_NMDS <- metaMDS(SD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07970526 
    ## Run 1 stress 0.07250814 
    ## ... New best solution
    ## ... Procrustes: rmse 0.08998187  max resid 0.2576221 
    ## Run 2 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04031666  max resid 0.1229608 
    ## Run 3 stress 0.08340287 
    ## Run 4 stress 0.07250813 
    ## Run 5 stress 0.08340291 
    ## Run 6 stress 0.06942776 
    ## ... Procrustes: rmse 3.27561e-06  max resid 8.487142e-06 
    ## ... Similar to previous best
    ## Run 7 stress 0.08340288 
    ## Run 8 stress 0.06942776 
    ## ... Procrustes: rmse 1.265807e-05  max resid 3.277398e-05 
    ## ... Similar to previous best
    ## Run 9 stress 0.07250812 
    ## Run 10 stress 0.07970525 
    ## Run 11 stress 0.07970525 
    ## Run 12 stress 0.07250812 
    ## Run 13 stress 0.07970525 
    ## Run 14 stress 0.07250815 
    ## Run 15 stress 0.07970526 
    ## Run 16 stress 0.07250812 
    ## Run 17 stress 0.07428313 
    ## Run 18 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321345  max resid 0.03324083 
    ## Run 19 stress 0.07250812 
    ## Run 20 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323408  max resid 0.0332834 
    ## Run 21 stress 0.08340288 
    ## Run 22 stress 0.08340291 
    ## Run 23 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325636  max resid 0.03333081 
    ## Run 24 stress 0.08448446 
    ## Run 25 stress 0.06978198 
    ## ... Procrustes: rmse 0.01326698  max resid 0.03335665 
    ## Run 26 stress 0.08340289 
    ## Run 27 stress 0.08448444 
    ## Run 28 stress 0.07970525 
    ## Run 29 stress 0.06942776 
    ## ... Procrustes: rmse 2.559045e-05  max resid 6.614581e-05 
    ## ... Similar to previous best
    ## Run 30 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325668  max resid 0.03333083 
    ## Run 31 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319561  max resid 0.03320256 
    ## Run 32 stress 0.06942776 
    ## ... Procrustes: rmse 3.031678e-05  max resid 7.764079e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.08340292 
    ## Run 34 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323391  max resid 0.03328534 
    ## Run 35 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323451  max resid 0.03328406 
    ## Run 36 stress 0.06942776 
    ## ... Procrustes: rmse 2.932048e-05  max resid 7.593032e-05 
    ## ... Similar to previous best
    ## Run 37 stress 0.07250815 
    ## Run 38 stress 0.07970526 
    ## Run 39 stress 0.069782 
    ## ... Procrustes: rmse 0.01328098  max resid 0.03337799 
    ## Run 40 stress 0.07428312 
    ## Run 41 stress 0.07428315 
    ## Run 42 stress 0.08340287 
    ## Run 43 stress 0.06978192 
    ## ... Procrustes: rmse 0.01318974  max resid 0.03318527 
    ## Run 44 stress 0.06942777 
    ## ... Procrustes: rmse 4.888412e-05  max resid 0.0001260572 
    ## ... Similar to previous best
    ## Run 45 stress 0.06978191 
    ## ... Procrustes: rmse 0.01322507  max resid 0.0332691 
    ## Run 46 stress 0.06942778 
    ## ... Procrustes: rmse 6.009185e-05  max resid 0.0001546539 
    ## ... Similar to previous best
    ## Run 47 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324455  max resid 0.0333063 
    ## Run 48 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319993  max resid 0.03321095 
    ## Run 49 stress 0.07250812 
    ## Run 50 stress 0.07428312 
    ## Run 51 stress 0.0834029 
    ## Run 52 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326608  max resid 0.03334839 
    ## Run 53 stress 0.06942777 
    ## ... Procrustes: rmse 3.076555e-05  max resid 7.960939e-05 
    ## ... Similar to previous best
    ## Run 54 stress 0.08340287 
    ## Run 55 stress 0.07250812 
    ## Run 56 stress 0.07250812 
    ## Run 57 stress 0.0742832 
    ## Run 58 stress 0.08340287 
    ## Run 59 stress 0.06942777 
    ## ... Procrustes: rmse 3.652495e-05  max resid 9.372678e-05 
    ## ... Similar to previous best
    ## Run 60 stress 0.08340287 
    ## Run 61 stress 0.08448461 
    ## Run 62 stress 0.07844932 
    ## Run 63 stress 0.06942778 
    ## ... Procrustes: rmse 4.951366e-05  max resid 0.0001300636 
    ## ... Similar to previous best
    ## Run 64 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322022  max resid 0.03325831 
    ## Run 65 stress 0.08340296 
    ## Run 66 stress 0.07970525 
    ## Run 67 stress 0.06978197 
    ## ... Procrustes: rmse 0.01321264  max resid 0.0332456 
    ## Run 68 stress 0.07970525 
    ## Run 69 stress 0.06978196 
    ## ... Procrustes: rmse 0.0131753  max resid 0.03315896 
    ## Run 70 stress 0.06942778 
    ## ... Procrustes: rmse 6.892876e-05  max resid 0.0001773157 
    ## ... Similar to previous best
    ## Run 71 stress 0.07970525 
    ## Run 72 stress 0.06942776 
    ## ... Procrustes: rmse 3.402417e-06  max resid 7.422997e-06 
    ## ... Similar to previous best
    ## Run 73 stress 0.07250816 
    ## Run 74 stress 0.07250812 
    ## Run 75 stress 0.06942777 
    ## ... Procrustes: rmse 4.700001e-05  max resid 0.0001211153 
    ## ... Similar to previous best
    ## Run 76 stress 0.07970527 
    ## Run 77 stress 0.06942776 
    ## ... Procrustes: rmse 1.087197e-05  max resid 2.813199e-05 
    ## ... Similar to previous best
    ## Run 78 stress 0.08340287 
    ## Run 79 stress 0.06942776 
    ## ... Procrustes: rmse 1.481071e-05  max resid 3.850883e-05 
    ## ... Similar to previous best
    ## Run 80 stress 0.06978203 
    ## ... Procrustes: rmse 0.01315642  max resid 0.0331178 
    ## Run 81 stress 0.07428313 
    ## Run 82 stress 0.08340286 
    ## Run 83 stress 0.06978199 
    ## ... Procrustes: rmse 0.01316494  max resid 0.03313406 
    ## Run 84 stress 0.06942776 
    ## ... Procrustes: rmse 1.229588e-06  max resid 2.634591e-06 
    ## ... Similar to previous best
    ## Run 85 stress 0.06978193 
    ## ... Procrustes: rmse 0.01324821  max resid 0.03331696 
    ## Run 86 stress 0.06942776 
    ## ... Procrustes: rmse 6.739207e-06  max resid 1.777156e-05 
    ## ... Similar to previous best
    ## Run 87 stress 0.08340287 
    ## Run 88 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325973  max resid 0.0333378 
    ## Run 89 stress 0.08340301 
    ## Run 90 stress 0.07250813 
    ## Run 91 stress 0.07970525 
    ## Run 92 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319738  max resid 0.03320844 
    ## Run 93 stress 0.07970526 
    ## Run 94 stress 0.08340287 
    ## Run 95 stress 0.08340287 
    ## Run 96 stress 0.06942776 
    ## ... Procrustes: rmse 2.689377e-05  max resid 6.946483e-05 
    ## ... Similar to previous best
    ## Run 97 stress 0.07428316 
    ## Run 98 stress 0.07250812 
    ## Run 99 stress 0.08340287 
    ## Run 100 stress 0.08448448 
    ## Run 101 stress 0.07250816 
    ## Run 102 stress 0.07428314 
    ## Run 103 stress 0.0834029 
    ## Run 104 stress 0.08340286 
    ## Run 105 stress 0.07428317 
    ## Run 106 stress 0.08340287 
    ## Run 107 stress 0.07970527 
    ## Run 108 stress 0.08340287 
    ## Run 109 stress 0.06942776 
    ## ... Procrustes: rmse 1.857086e-05  max resid 4.799599e-05 
    ## ... Similar to previous best
    ## Run 110 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325725  max resid 0.03333181 
    ## Run 111 stress 0.07970525 
    ## Run 112 stress 0.07250812 
    ## Run 113 stress 0.08340287 
    ## Run 114 stress 0.08448437 
    ## Run 115 stress 0.06942776 
    ## ... Procrustes: rmse 2.098811e-05  max resid 5.361038e-05 
    ## ... Similar to previous best
    ## Run 116 stress 0.07970526 
    ## Run 117 stress 0.07428312 
    ## Run 118 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324937  max resid 0.03331981 
    ## Run 119 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.282509e-06  max resid 4.514857e-06 
    ## ... Similar to previous best
    ## Run 120 stress 0.06942776 
    ## ... Procrustes: rmse 4.723168e-06  max resid 1.013719e-05 
    ## ... Similar to previous best
    ## Run 121 stress 0.0844844 
    ## Run 122 stress 0.07428313 
    ## Run 123 stress 0.08340287 
    ## Run 124 stress 0.07970526 
    ## Run 125 stress 0.07428315 
    ## Run 126 stress 0.08340288 
    ## Run 127 stress 0.08448437 
    ## Run 128 stress 0.07428314 
    ## Run 129 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326975  max resid 0.03336089 
    ## Run 130 stress 0.08340288 
    ## Run 131 stress 0.07970525 
    ## Run 132 stress 0.07428314 
    ## Run 133 stress 0.08340292 
    ## Run 134 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323705  max resid 0.03329418 
    ## Run 135 stress 0.069782 
    ## ... Procrustes: rmse 0.01327893  max resid 0.03337692 
    ## Run 136 stress 0.07250812 
    ## Run 137 stress 0.06942776 
    ## ... Procrustes: rmse 4.308693e-06  max resid 1.127209e-05 
    ## ... Similar to previous best
    ## Run 138 stress 0.07844945 
    ## Run 139 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323725  max resid 0.03329189 
    ## Run 140 stress 0.07970525 
    ## Run 141 stress 0.07250812 
    ## Run 142 stress 0.07250812 
    ## Run 143 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322928  max resid 0.03327702 
    ## Run 144 stress 0.07970525 
    ## Run 145 stress 0.07250812 
    ## Run 146 stress 0.0834029 
    ## Run 147 stress 0.08340289 
    ## Run 148 stress 0.06942776 
    ## ... Procrustes: rmse 4.772165e-06  max resid 1.181103e-05 
    ## ... Similar to previous best
    ## Run 149 stress 0.07428313 
    ## Run 150 stress 0.07970526 
    ## Run 151 stress 0.07970525 
    ## Run 152 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326091  max resid 0.03334117 
    ## Run 153 stress 0.07250812 
    ## Run 154 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327742  max resid 0.03337639 
    ## Run 155 stress 0.07970525 
    ## Run 156 stress 0.0834029 
    ## Run 157 stress 0.07970525 
    ## Run 158 stress 0.06942776 
    ## ... Procrustes: rmse 3.844933e-06  max resid 9.942929e-06 
    ## ... Similar to previous best
    ## Run 159 stress 0.07250812 
    ## Run 160 stress 0.08448437 
    ## Run 161 stress 0.08340288 
    ## Run 162 stress 0.06942777 
    ## ... Procrustes: rmse 4.017424e-05  max resid 0.0001034603 
    ## ... Similar to previous best
    ## Run 163 stress 0.06978197 
    ## ... Procrustes: rmse 0.01317172  max resid 0.03315256 
    ## Run 164 stress 0.07428313 
    ## Run 165 stress 0.06978191 
    ## ... Procrustes: rmse 0.01320168  max resid 0.03321873 
    ## Run 166 stress 0.07250812 
    ## Run 167 stress 0.07844967 
    ## Run 168 stress 0.07970527 
    ## Run 169 stress 0.08448441 
    ## Run 170 stress 0.0834029 
    ## Run 171 stress 0.07250812 
    ## Run 172 stress 0.07428313 
    ## Run 173 stress 0.06942776 
    ## ... Procrustes: rmse 2.38894e-05  max resid 6.129233e-05 
    ## ... Similar to previous best
    ## Run 174 stress 0.07428313 
    ## Run 175 stress 0.07970525 
    ## Run 176 stress 0.07428317 
    ## Run 177 stress 0.08340289 
    ## Run 178 stress 0.07250812 
    ## Run 179 stress 0.06942776 
    ## ... Procrustes: rmse 8.558445e-06  max resid 2.240321e-05 
    ## ... Similar to previous best
    ## Run 180 stress 0.07250812 
    ## Run 181 stress 0.06942776 
    ## ... Procrustes: rmse 4.635815e-06  max resid 1.199922e-05 
    ## ... Similar to previous best
    ## Run 182 stress 0.07970526 
    ## Run 183 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327616  max resid 0.03337083 
    ## Run 184 stress 0.07250814 
    ## Run 185 stress 0.06942776 
    ## ... Procrustes: rmse 3.387883e-05  max resid 8.825182e-05 
    ## ... Similar to previous best
    ## Run 186 stress 0.08340287 
    ## Run 187 stress 0.07428312 
    ## Run 188 stress 0.07250813 
    ## Run 189 stress 0.07428314 
    ## Run 190 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327367  max resid 0.03336744 
    ## Run 191 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325102  max resid 0.03332662 
    ## Run 192 stress 0.07970526 
    ## Run 193 stress 0.07250814 
    ## Run 194 stress 0.08340287 
    ## Run 195 stress 0.07250813 
    ## Run 196 stress 0.07428313 
    ## Run 197 stress 0.08340289 
    ## Run 198 stress 0.07428316 
    ## Run 199 stress 0.06942776 
    ## ... Procrustes: rmse 2.55239e-06  max resid 3.795441e-06 
    ## ... Similar to previous best
    ## Run 200 stress 0.07428313 
    ## Run 201 stress 0.07970525 
    ## Run 202 stress 0.08448435 
    ## Run 203 stress 0.07970525 
    ## Run 204 stress 0.07250812 
    ## Run 205 stress 0.06978193 
    ## ... Procrustes: rmse 0.01319038  max resid 0.03319301 
    ## Run 206 stress 0.07250815 
    ## Run 207 stress 0.07428312 
    ## Run 208 stress 0.06942776 
    ## ... Procrustes: rmse 1.025539e-05  max resid 2.844827e-05 
    ## ... Similar to previous best
    ## Run 209 stress 0.08340286 
    ## Run 210 stress 0.07970527 
    ## Run 211 stress 0.06978205 
    ## ... Procrustes: rmse 0.01329531  max resid 0.03341244 
    ## Run 212 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325427  max resid 0.03333195 
    ## Run 213 stress 0.08340287 
    ## Run 214 stress 0.07428318 
    ## Run 215 stress 0.07970525 
    ## Run 216 stress 0.07844944 
    ## Run 217 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324449  max resid 0.03330909 
    ## Run 218 stress 0.07428315 
    ## Run 219 stress 0.07428314 
    ## Run 220 stress 0.07250815 
    ## Run 221 stress 0.07428312 
    ## Run 222 stress 0.06942777 
    ## ... Procrustes: rmse 3.095284e-05  max resid 8.133516e-05 
    ## ... Similar to previous best
    ## Run 223 stress 0.07428315 
    ## Run 224 stress 0.07970525 
    ## Run 225 stress 0.07250813 
    ## Run 226 stress 0.06942778 
    ## ... Procrustes: rmse 5.998547e-05  max resid 0.0001539717 
    ## ... Similar to previous best
    ## Run 227 stress 0.07250812 
    ## Run 228 stress 0.08340289 
    ## Run 229 stress 0.07970525 
    ## Run 230 stress 0.08340296 
    ## Run 231 stress 0.07428316 
    ## Run 232 stress 0.07428314 
    ## Run 233 stress 0.07428313 
    ## Run 234 stress 0.07970525 
    ## Run 235 stress 0.0742832 
    ## Run 236 stress 0.07428313 
    ## Run 237 stress 0.07428315 
    ## Run 238 stress 0.07428313 
    ## Run 239 stress 0.06978195 
    ## ... Procrustes: rmse 0.01317754  max resid 0.03316662 
    ## Run 240 stress 0.07250813 
    ## Run 241 stress 0.07250814 
    ## Run 242 stress 0.07428313 
    ## Run 243 stress 0.08448445 
    ## Run 244 stress 0.07250812 
    ## Run 245 stress 0.08340287 
    ## Run 246 stress 0.07250812 
    ## Run 247 stress 0.07844949 
    ## Run 248 stress 0.0697819 
    ## ... Procrustes: rmse 0.013233  max resid 0.03328506 
    ## Run 249 stress 0.06942776 
    ## ... Procrustes: rmse 7.871147e-06  max resid 1.990777e-05 
    ## ... Similar to previous best
    ## Run 250 stress 0.06942779 
    ## ... Procrustes: rmse 7.44684e-05  max resid 0.0001910862 
    ## ... Similar to previous best
    ## Run 251 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 1.018898e-06  max resid 2.491899e-06 
    ## ... Similar to previous best
    ## Run 252 stress 0.06978203 
    ## ... Procrustes: rmse 0.01328116  max resid 0.03338712 
    ## Run 253 stress 0.07250814 
    ## Run 254 stress 0.07970526 
    ## Run 255 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324607  max resid 0.03331158 
    ## Run 256 stress 0.07428318 
    ## Run 257 stress 0.06942785 
    ## ... Procrustes: rmse 5.11398e-05  max resid 0.0001448549 
    ## ... Similar to previous best
    ## Run 258 stress 0.06942776 
    ## ... Procrustes: rmse 1.086199e-05  max resid 2.801065e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.07428313 
    ## Run 260 stress 0.07428312 
    ## Run 261 stress 0.08340297 
    ## Run 262 stress 0.08340289 
    ## Run 263 stress 0.07428313 
    ## Run 264 stress 0.08340287 
    ## Run 265 stress 0.07250812 
    ## Run 266 stress 0.07970526 
    ## Run 267 stress 0.07844919 
    ## Run 268 stress 0.06978201 
    ## ... Procrustes: rmse 0.01328544  max resid 0.03339059 
    ## Run 269 stress 0.08340287 
    ## Run 270 stress 0.06942776 
    ## ... Procrustes: rmse 2.04476e-05  max resid 5.281085e-05 
    ## ... Similar to previous best
    ## Run 271 stress 0.08340286 
    ## Run 272 stress 0.07970525 
    ## Run 273 stress 0.07428313 
    ## Run 274 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325509  max resid 0.03332693 
    ## Run 275 stress 0.07428317 
    ## Run 276 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320769  max resid 0.03323005 
    ## Run 277 stress 0.07428313 
    ## Run 278 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319663  max resid 0.0332115 
    ## Run 279 stress 0.07970525 
    ## Run 280 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324303  max resid 0.03330318 
    ## Run 281 stress 0.06942778 
    ## ... Procrustes: rmse 6.470266e-05  max resid 0.0001666741 
    ## ... Similar to previous best
    ## Run 282 stress 0.07428312 
    ## Run 283 stress 0.08448446 
    ## Run 284 stress 0.06942776 
    ## ... Procrustes: rmse 3.857763e-06  max resid 6.209192e-06 
    ## ... Similar to previous best
    ## Run 285 stress 0.06942776 
    ## ... Procrustes: rmse 7.263687e-06  max resid 1.873963e-05 
    ## ... Similar to previous best
    ## Run 286 stress 0.06942776 
    ## ... Procrustes: rmse 6.298098e-06  max resid 1.63101e-05 
    ## ... Similar to previous best
    ## Run 287 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322199  max resid 0.03326181 
    ## Run 288 stress 0.07970525 
    ## Run 289 stress 0.08448455 
    ## Run 290 stress 0.06942776 
    ## ... Procrustes: rmse 1.120905e-05  max resid 2.892541e-05 
    ## ... Similar to previous best
    ## Run 291 stress 0.06942776 
    ## ... Procrustes: rmse 1.410957e-05  max resid 3.703329e-05 
    ## ... Similar to previous best
    ## Run 292 stress 0.08340289 
    ## Run 293 stress 0.06942777 
    ## ... Procrustes: rmse 2.046615e-05  max resid 5.480816e-05 
    ## ... Similar to previous best
    ## Run 294 stress 0.07428314 
    ## Run 295 stress 0.07970526 
    ## Run 296 stress 0.07428313 
    ## Run 297 stress 0.07970526 
    ## Run 298 stress 0.07428317 
    ## Run 299 stress 0.07428313 
    ## Run 300 stress 0.06978193 
    ## ... Procrustes: rmse 0.01324079  max resid 0.03329666 
    ## Run 301 stress 0.08340288 
    ## Run 302 stress 0.07428335 
    ## Run 303 stress 0.06978195 
    ## ... Procrustes: rmse 0.01317971  max resid 0.03316895 
    ## Run 304 stress 0.07970526 
    ## Run 305 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325391  max resid 0.0333265 
    ## Run 306 stress 0.08340289 
    ## Run 307 stress 0.07844954 
    ## Run 308 stress 0.06978199 
    ## ... Procrustes: rmse 0.01327871  max resid 0.03337782 
    ## Run 309 stress 0.06942778 
    ## ... Procrustes: rmse 6.427767e-05  max resid 0.000165296 
    ## ... Similar to previous best
    ## Run 310 stress 0.07844933 
    ## Run 311 stress 0.07970525 
    ## Run 312 stress 0.07428313 
    ## Run 313 stress 0.07970525 
    ## Run 314 stress 0.07250815 
    ## Run 315 stress 0.08340288 
    ## Run 316 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321474  max resid 0.03324375 
    ## Run 317 stress 0.06978191 
    ## ... Procrustes: rmse 0.01320296  max resid 0.03322103 
    ## Run 318 stress 0.0834029 
    ## Run 319 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323779  max resid 0.03329453 
    ## Run 320 stress 0.08340292 
    ## Run 321 stress 0.08448456 
    ## Run 322 stress 0.08340288 
    ## Run 323 stress 0.07428313 
    ## Run 324 stress 0.08340286 
    ## Run 325 stress 0.06942776 
    ## ... Procrustes: rmse 1.278897e-05  max resid 3.322879e-05 
    ## ... Similar to previous best
    ## Run 326 stress 0.06942777 
    ## ... Procrustes: rmse 4.891574e-05  max resid 0.0001258831 
    ## ... Similar to previous best
    ## Run 327 stress 0.07428316 
    ## Run 328 stress 0.07250812 
    ## Run 329 stress 0.08340287 
    ## Run 330 stress 0.08340287 
    ## Run 331 stress 0.08340287 
    ## Run 332 stress 0.07250813 
    ## Run 333 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319717  max resid 0.03320688 
    ## Run 334 stress 0.07970525 
    ## Run 335 stress 0.06942776 
    ## ... Procrustes: rmse 6.863783e-07  max resid 1.274581e-06 
    ## ... Similar to previous best
    ## Run 336 stress 0.07428313 
    ## Run 337 stress 0.07970525 
    ## Run 338 stress 0.07844958 
    ## Run 339 stress 0.06978193 
    ## ... Procrustes: rmse 0.01318905  max resid 0.03319157 
    ## Run 340 stress 0.07970525 
    ## Run 341 stress 0.07250812 
    ## Run 342 stress 0.07250813 
    ## Run 343 stress 0.06978195 
    ## ... Procrustes: rmse 0.01318335  max resid 0.03317737 
    ## Run 344 stress 0.06942776 
    ## ... Procrustes: rmse 7.020495e-06  max resid 1.815357e-05 
    ## ... Similar to previous best
    ## Run 345 stress 0.08340287 
    ## Run 346 stress 0.08340289 
    ## Run 347 stress 0.06942777 
    ## ... Procrustes: rmse 3.415272e-05  max resid 8.945637e-05 
    ## ... Similar to previous best
    ## Run 348 stress 0.08448445 
    ## Run 349 stress 0.07970525 
    ## Run 350 stress 0.06942776 
    ## ... Procrustes: rmse 2.113145e-05  max resid 5.204384e-05 
    ## ... Similar to previous best
    ## Run 351 stress 0.06978201 
    ## ... Procrustes: rmse 0.01328406  max resid 0.03338728 
    ## Run 352 stress 0.08340287 
    ## Run 353 stress 0.07970525 
    ## Run 354 stress 0.08340296 
    ## Run 355 stress 0.06942776 
    ## ... Procrustes: rmse 6.457972e-06  max resid 1.669474e-05 
    ## ... Similar to previous best
    ## Run 356 stress 0.08340286 
    ## Run 357 stress 0.07970525 
    ## Run 358 stress 0.07250812 
    ## Run 359 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325213  max resid 0.0333245 
    ## Run 360 stress 0.07428321 
    ## Run 361 stress 0.07428314 
    ## Run 362 stress 0.06942776 
    ## ... Procrustes: rmse 8.570527e-07  max resid 2.041015e-06 
    ## ... Similar to previous best
    ## Run 363 stress 0.07428314 
    ## Run 364 stress 0.07970526 
    ## Run 365 stress 0.07428312 
    ## Run 366 stress 0.08340288 
    ## Run 367 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326033  max resid 0.03334293 
    ## Run 368 stress 0.07428316 
    ## Run 369 stress 0.06942776 
    ## ... Procrustes: rmse 2.207426e-05  max resid 5.699311e-05 
    ## ... Similar to previous best
    ## Run 370 stress 0.07970525 
    ## Run 371 stress 0.08448444 
    ## Run 372 stress 0.07428316 
    ## Run 373 stress 0.07428314 
    ## Run 374 stress 0.08340288 
    ## Run 375 stress 0.06942776 
    ## ... Procrustes: rmse 3.14174e-06  max resid 8.109758e-06 
    ## ... Similar to previous best
    ## Run 376 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327581  max resid 0.03337427 
    ## Run 377 stress 0.07970526 
    ## Run 378 stress 0.07970526 
    ## Run 379 stress 0.08340288 
    ## Run 380 stress 0.07250812 
    ## Run 381 stress 0.08340288 
    ## Run 382 stress 0.07250814 
    ## Run 383 stress 0.07428313 
    ## Run 384 stress 0.07428313 
    ## Run 385 stress 0.07428313 
    ## Run 386 stress 0.07970525 
    ## Run 387 stress 0.08340287 
    ## Run 388 stress 0.07428315 
    ## Run 389 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323003  max resid 0.0332778 
    ## Run 390 stress 0.06942777 
    ## ... Procrustes: rmse 4.643636e-05  max resid 0.0001202895 
    ## ... Similar to previous best
    ## Run 391 stress 0.08448437 
    ## Run 392 stress 0.06978193 
    ## ... Procrustes: rmse 0.01318686  max resid 0.03318647 
    ## Run 393 stress 0.07250812 
    ## Run 394 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327392  max resid 0.03336445 
    ## Run 395 stress 0.07250813 
    ## Run 396 stress 0.07970526 
    ## Run 397 stress 0.07428314 
    ## Run 398 stress 0.08340292 
    ## Run 399 stress 0.07970525 
    ## Run 400 stress 0.07250813 
    ## Run 401 stress 0.07428315 
    ## Run 402 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324562  max resid 0.03330925 
    ## Run 403 stress 0.07428313 
    ## Run 404 stress 0.07970526 
    ## Run 405 stress 0.07428313 
    ## Run 406 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132369  max resid 0.03329254 
    ## Run 407 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326219  max resid 0.0333435 
    ## Run 408 stress 0.07970526 
    ## Run 409 stress 0.07250812 
    ## Run 410 stress 0.08340293 
    ## Run 411 stress 0.07250812 
    ## Run 412 stress 0.0834029 
    ## Run 413 stress 0.07428314 
    ## Run 414 stress 0.07970525 
    ## Run 415 stress 0.07970525 
    ## Run 416 stress 0.07250812 
    ## Run 417 stress 0.0834029 
    ## Run 418 stress 0.07428313 
    ## Run 419 stress 0.0697819 
    ## ... Procrustes: rmse 0.013236  max resid 0.03329095 
    ## Run 420 stress 0.08340287 
    ## Run 421 stress 0.07250813 
    ## Run 422 stress 0.07250812 
    ## Run 423 stress 0.07250814 
    ## Run 424 stress 0.07250812 
    ## Run 425 stress 0.07428316 
    ## Run 426 stress 0.08340292 
    ## Run 427 stress 0.07250812 
    ## Run 428 stress 0.07428313 
    ## Run 429 stress 0.0784495 
    ## Run 430 stress 0.06942776 
    ## ... Procrustes: rmse 5.426909e-06  max resid 1.403654e-05 
    ## ... Similar to previous best
    ## Run 431 stress 0.06942777 
    ## ... Procrustes: rmse 4.764197e-05  max resid 0.0001225315 
    ## ... Similar to previous best
    ## Run 432 stress 0.07970526 
    ## Run 433 stress 0.06942778 
    ## ... Procrustes: rmse 5.815371e-05  max resid 0.0001498642 
    ## ... Similar to previous best
    ## Run 434 stress 0.069782 
    ## ... Procrustes: rmse 0.01316233  max resid 0.03313193 
    ## Run 435 stress 0.07428313 
    ## Run 436 stress 0.0834029 
    ## Run 437 stress 0.07250813 
    ## Run 438 stress 0.08340289 
    ## Run 439 stress 0.07428313 
    ## Run 440 stress 0.08340288 
    ## Run 441 stress 0.08340286 
    ## Run 442 stress 0.07428315 
    ## Run 443 stress 0.08340287 
    ## Run 444 stress 0.07250812 
    ## Run 445 stress 0.06942776 
    ## ... Procrustes: rmse 2.842132e-05  max resid 7.318546e-05 
    ## ... Similar to previous best
    ## Run 446 stress 0.06978194 
    ## ... Procrustes: rmse 0.0132599  max resid 0.03333918 
    ## Run 447 stress 0.06942778 
    ## ... Procrustes: rmse 6.613184e-05  max resid 0.0001708073 
    ## ... Similar to previous best
    ## Run 448 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324576  max resid 0.03331138 
    ## Run 449 stress 0.06942776 
    ## ... Procrustes: rmse 3.034756e-06  max resid 8.242044e-06 
    ## ... Similar to previous best
    ## Run 450 stress 0.07428313 
    ## Run 451 stress 0.07428313 
    ## Run 452 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323265  max resid 0.03328384 
    ## Run 453 stress 0.07428313 
    ## Run 454 stress 0.08340287 
    ## Run 455 stress 0.06942777 
    ## ... Procrustes: rmse 4.107533e-05  max resid 0.0001057698 
    ## ... Similar to previous best
    ## Run 456 stress 0.07428313 
    ## Run 457 stress 0.07428313 
    ## Run 458 stress 0.07428322 
    ## Run 459 stress 0.07250817 
    ## Run 460 stress 0.07844926 
    ## Run 461 stress 0.08340286 
    ## Run 462 stress 0.07428313 
    ## Run 463 stress 0.06978197 
    ## ... Procrustes: rmse 0.01317134  max resid 0.03315091 
    ## Run 464 stress 0.08340289 
    ## Run 465 stress 0.07250813 
    ## Run 466 stress 0.08340291 
    ## Run 467 stress 0.07428315 
    ## Run 468 stress 0.07428313 
    ## Run 469 stress 0.06978195 
    ## ... Procrustes: rmse 0.01318541  max resid 0.03317911 
    ## Run 470 stress 0.07250812 
    ## Run 471 stress 0.07250812 
    ## Run 472 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325895  max resid 0.03333822 
    ## Run 473 stress 0.08340291 
    ## Run 474 stress 0.07428317 
    ## Run 475 stress 0.06942776 
    ## ... Procrustes: rmse 1.334644e-05  max resid 3.502983e-05 
    ## ... Similar to previous best
    ## Run 476 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324331  max resid 0.03330802 
    ## Run 477 stress 0.08340287 
    ## Run 478 stress 0.07428317 
    ## Run 479 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323085  max resid 0.03327811 
    ## Run 480 stress 0.06942776 
    ## ... Procrustes: rmse 1.738314e-05  max resid 4.509031e-05 
    ## ... Similar to previous best
    ## Run 481 stress 0.06978193 
    ## ... Procrustes: rmse 0.01318912  max resid 0.03318932 
    ## Run 482 stress 0.08340288 
    ## Run 483 stress 0.08340288 
    ## Run 484 stress 0.07970525 
    ## Run 485 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323664  max resid 0.03329069 
    ## Run 486 stress 0.06942776 
    ## ... Procrustes: rmse 1.835045e-05  max resid 4.691927e-05 
    ## ... Similar to previous best
    ## Run 487 stress 0.06942778 
    ## ... Procrustes: rmse 6.272993e-05  max resid 0.0001625099 
    ## ... Similar to previous best
    ## Run 488 stress 0.07250812 
    ## Run 489 stress 0.06978201 
    ## ... Procrustes: rmse 0.0132861  max resid 0.03339431 
    ## Run 490 stress 0.06978193 
    ## ... Procrustes: rmse 0.01318611  max resid 0.03318269 
    ## Run 491 stress 0.07428313 
    ## Run 492 stress 0.07428313 
    ## Run 493 stress 0.06942776 
    ## ... Procrustes: rmse 5.694677e-06  max resid 1.462935e-05 
    ## ... Similar to previous best
    ## Run 494 stress 0.08448437 
    ## Run 495 stress 0.06978197 
    ## ... Procrustes: rmse 0.01317622  max resid 0.03315894 
    ## Run 496 stress 0.08340289 
    ## Run 497 stress 0.08340289 
    ## Run 498 stress 0.06942777 
    ## ... Procrustes: rmse 3.636689e-05  max resid 9.381364e-05 
    ## ... Similar to previous best
    ## Run 499 stress 0.06978202 
    ## ... Procrustes: rmse 0.01327869  max resid 0.03338107 
    ## Run 500 stress 0.07970525 
    ## *** Best solution repeated 36 times

``` r
# Stratified lakes
SD_beta_S_NMDS <- metaMDS(SD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1410449 
    ## Run 1 stress 0.1445811 
    ## Run 2 stress 0.1551863 
    ## Run 3 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 0.2673131  max resid 0.5425958 
    ## Run 4 stress 0.1407298 
    ## Run 5 stress 0.1407298 
    ## Run 6 stress 0.1771969 
    ## Run 7 stress 0.1415299 
    ## Run 8 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 9.021358e-07  max resid 1.77068e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.1320258 
    ## ... Procrustes: rmse 1.99024e-06  max resid 4.164428e-06 
    ## ... Similar to previous best
    ## Run 10 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 7.685121e-07  max resid 1.277279e-06 
    ## ... Similar to previous best
    ## Run 11 stress 0.1383682 
    ## Run 12 stress 0.1970303 
    ## Run 13 stress 0.1407208 
    ## Run 14 stress 0.1572668 
    ## Run 15 stress 0.1551863 
    ## Run 16 stress 0.2885554 
    ## Run 17 stress 0.1693717 
    ## Run 18 stress 0.1551863 
    ## Run 19 stress 0.1693717 
    ## Run 20 stress 0.1415299 
    ## Run 21 stress 0.1407298 
    ## Run 22 stress 0.1410448 
    ## Run 23 stress 0.1434197 
    ## Run 24 stress 0.2385029 
    ## Run 25 stress 0.1320258 
    ## ... Procrustes: rmse 1.390287e-06  max resid 2.355699e-06 
    ## ... Similar to previous best
    ## Run 26 stress 0.1320258 
    ## ... Procrustes: rmse 6.442922e-07  max resid 1.327677e-06 
    ## ... Similar to previous best
    ## Run 27 stress 0.3071885 
    ## Run 28 stress 0.1320258 
    ## ... Procrustes: rmse 1.387189e-06  max resid 2.651134e-06 
    ## ... Similar to previous best
    ## Run 29 stress 0.1434197 
    ## Run 30 stress 0.1407298 
    ## Run 31 stress 0.1407207 
    ## Run 32 stress 0.1415299 
    ## Run 33 stress 0.2416997 
    ## Run 34 stress 0.1407298 
    ## Run 35 stress 0.2578639 
    ## Run 36 stress 0.1407298 
    ## Run 37 stress 0.1693717 
    ## Run 38 stress 0.1407298 
    ## Run 39 stress 0.1663544 
    ## Run 40 stress 0.1407207 
    ## Run 41 stress 0.2224299 
    ## Run 42 stress 0.141045 
    ## Run 43 stress 0.1407298 
    ## Run 44 stress 0.1383682 
    ## Run 45 stress 0.141045 
    ## Run 46 stress 0.1572668 
    ## Run 47 stress 0.1551863 
    ## Run 48 stress 0.1410453 
    ## Run 49 stress 0.1383681 
    ## Run 50 stress 0.1407298 
    ## Run 51 stress 0.1415299 
    ## Run 52 stress 0.1663544 
    ## Run 53 stress 0.2885554 
    ## Run 54 stress 0.1415299 
    ## Run 55 stress 0.1771969 
    ## Run 56 stress 0.1383681 
    ## Run 57 stress 0.1415299 
    ## Run 58 stress 0.1383681 
    ## Run 59 stress 0.1587623 
    ## Run 60 stress 0.1320258 
    ## ... Procrustes: rmse 2.242204e-06  max resid 4.592218e-06 
    ## ... Similar to previous best
    ## Run 61 stress 0.1320258 
    ## ... Procrustes: rmse 1.408002e-06  max resid 2.890859e-06 
    ## ... Similar to previous best
    ## Run 62 stress 0.1383682 
    ## Run 63 stress 0.1410451 
    ## Run 64 stress 0.1572668 
    ## Run 65 stress 0.1583016 
    ## Run 66 stress 0.1929519 
    ## Run 67 stress 0.1410448 
    ## Run 68 stress 0.1572668 
    ## Run 69 stress 0.1970303 
    ## Run 70 stress 0.1415299 
    ## Run 71 stress 0.1410449 
    ## Run 72 stress 0.1830332 
    ## Run 73 stress 0.1410451 
    ## Run 74 stress 0.1572668 
    ## Run 75 stress 0.1587623 
    ## Run 76 stress 0.1583016 
    ## Run 77 stress 0.1693717 
    ## Run 78 stress 0.1407208 
    ## Run 79 stress 0.1962476 
    ## Run 80 stress 0.1383681 
    ## Run 81 stress 0.1415299 
    ## Run 82 stress 0.1572668 
    ## Run 83 stress 0.1320258 
    ## ... Procrustes: rmse 3.260853e-07  max resid 5.717287e-07 
    ## ... Similar to previous best
    ## Run 84 stress 0.1410448 
    ## Run 85 stress 0.1407298 
    ## Run 86 stress 0.1551863 
    ## Run 87 stress 0.2385029 
    ## Run 88 stress 0.1598154 
    ## Run 89 stress 0.3083098 
    ## Run 90 stress 0.1415299 
    ## Run 91 stress 0.1407298 
    ## Run 92 stress 0.1628606 
    ## Run 93 stress 0.1415299 
    ## Run 94 stress 0.1320258 
    ## ... Procrustes: rmse 7.233589e-07  max resid 1.513782e-06 
    ## ... Similar to previous best
    ## Run 95 stress 0.1693717 
    ## Run 96 stress 0.1771966 
    ## Run 97 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 6.827554e-07  max resid 1.388187e-06 
    ## ... Similar to previous best
    ## Run 98 stress 0.1551863 
    ## Run 99 stress 0.1572668 
    ## Run 100 stress 0.2085297 
    ## Run 101 stress 0.1551863 
    ## Run 102 stress 0.1830332 
    ## Run 103 stress 0.1383681 
    ## Run 104 stress 0.1581642 
    ## Run 105 stress 0.1383682 
    ## Run 106 stress 0.1802752 
    ## Run 107 stress 0.1415299 
    ## Run 108 stress 0.2600627 
    ## Run 109 stress 0.1572668 
    ## Run 110 stress 0.1320258 
    ## ... Procrustes: rmse 1.968427e-06  max resid 3.923186e-06 
    ## ... Similar to previous best
    ## Run 111 stress 0.1693717 
    ## Run 112 stress 0.1415299 
    ## Run 113 stress 0.1320258 
    ## ... Procrustes: rmse 2.278258e-06  max resid 4.599195e-06 
    ## ... Similar to previous best
    ## Run 114 stress 0.2015531 
    ## Run 115 stress 0.1415299 
    ## Run 116 stress 0.1771966 
    ## Run 117 stress 0.2629507 
    ## Run 118 stress 0.1415299 
    ## Run 119 stress 0.1415299 
    ## Run 120 stress 0.1693717 
    ## Run 121 stress 0.1383681 
    ## Run 122 stress 0.1415299 
    ## Run 123 stress 0.1320258 
    ## ... Procrustes: rmse 3.683331e-07  max resid 7.522022e-07 
    ## ... Similar to previous best
    ## Run 124 stress 0.3047186 
    ## Run 125 stress 0.1693717 
    ## Run 126 stress 0.1572668 
    ## Run 127 stress 0.1320258 
    ## ... Procrustes: rmse 1.285576e-06  max resid 2.865381e-06 
    ## ... Similar to previous best
    ## Run 128 stress 0.1320258 
    ## ... Procrustes: rmse 1.678788e-06  max resid 3.502693e-06 
    ## ... Similar to previous best
    ## Run 129 stress 0.1383681 
    ## Run 130 stress 0.1410453 
    ## Run 131 stress 0.1407298 
    ## Run 132 stress 0.1410454 
    ## Run 133 stress 0.1410449 
    ## Run 134 stress 0.1410448 
    ## Run 135 stress 0.1970303 
    ## Run 136 stress 0.1551863 
    ## Run 137 stress 0.1962476 
    ## Run 138 stress 0.2227526 
    ## Run 139 stress 0.1383681 
    ## Run 140 stress 0.1551863 
    ## Run 141 stress 0.1320258 
    ## ... Procrustes: rmse 2.689494e-06  max resid 5.521289e-06 
    ## ... Similar to previous best
    ## Run 142 stress 0.2550043 
    ## Run 143 stress 0.1383682 
    ## Run 144 stress 0.1628606 
    ## Run 145 stress 0.2227527 
    ## Run 146 stress 0.1415299 
    ## Run 147 stress 0.2224299 
    ## Run 148 stress 0.1407298 
    ## Run 149 stress 0.1551863 
    ## Run 150 stress 0.1693717 
    ## Run 151 stress 0.1415299 
    ## Run 152 stress 0.1320258 
    ## ... Procrustes: rmse 2.296804e-06  max resid 4.677541e-06 
    ## ... Similar to previous best
    ## Run 153 stress 0.1776793 
    ## Run 154 stress 0.1415299 
    ## Run 155 stress 0.1663544 
    ## Run 156 stress 0.1551863 
    ## Run 157 stress 0.1410448 
    ## Run 158 stress 0.1583016 
    ## Run 159 stress 0.1320258 
    ## ... Procrustes: rmse 2.302404e-06  max resid 4.745642e-06 
    ## ... Similar to previous best
    ## Run 160 stress 0.1320258 
    ## ... Procrustes: rmse 3.218875e-07  max resid 5.292467e-07 
    ## ... Similar to previous best
    ## Run 161 stress 0.1415299 
    ## Run 162 stress 0.1320258 
    ## ... Procrustes: rmse 1.663241e-06  max resid 3.215367e-06 
    ## ... Similar to previous best
    ## Run 163 stress 0.3083098 
    ## Run 164 stress 0.1445811 
    ## Run 165 stress 0.2008163 
    ## Run 166 stress 0.1407298 
    ## Run 167 stress 0.1410448 
    ## Run 168 stress 0.1320258 
    ## ... Procrustes: rmse 1.524186e-07  max resid 2.890247e-07 
    ## ... Similar to previous best
    ## Run 169 stress 0.2224298 
    ## Run 170 stress 0.1383681 
    ## Run 171 stress 0.1551863 
    ## Run 172 stress 0.1320258 
    ## ... Procrustes: rmse 2.394681e-06  max resid 4.921279e-06 
    ## ... Similar to previous best
    ## Run 173 stress 0.1583016 
    ## Run 174 stress 0.1415299 
    ## Run 175 stress 0.1572668 
    ## Run 176 stress 0.1410451 
    ## Run 177 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 2.428592e-07  max resid 5.213326e-07 
    ## ... Similar to previous best
    ## Run 178 stress 0.1320258 
    ## ... Procrustes: rmse 1.059241e-06  max resid 2.148531e-06 
    ## ... Similar to previous best
    ## Run 179 stress 0.1407298 
    ## Run 180 stress 0.2520602 
    ## Run 181 stress 0.1434197 
    ## Run 182 stress 0.1415299 
    ## Run 183 stress 0.1383682 
    ## Run 184 stress 0.3080565 
    ## Run 185 stress 0.1320258 
    ## ... Procrustes: rmse 1.609739e-06  max resid 3.291235e-06 
    ## ... Similar to previous best
    ## Run 186 stress 0.1320258 
    ## ... Procrustes: rmse 1.661613e-06  max resid 3.36885e-06 
    ## ... Similar to previous best
    ## Run 187 stress 0.1410448 
    ## Run 188 stress 0.1383681 
    ## Run 189 stress 0.1383682 
    ## Run 190 stress 0.1320258 
    ## ... Procrustes: rmse 6.365561e-07  max resid 1.054417e-06 
    ## ... Similar to previous best
    ## Run 191 stress 0.2227527 
    ## Run 192 stress 0.1407298 
    ## Run 193 stress 0.2514078 
    ## Run 194 stress 0.1572668 
    ## Run 195 stress 0.1415299 
    ## Run 196 stress 0.1320258 
    ## ... Procrustes: rmse 3.457884e-06  max resid 7.355297e-06 
    ## ... Similar to previous best
    ## Run 197 stress 0.1572668 
    ## Run 198 stress 0.1410448 
    ## Run 199 stress 0.1320258 
    ## ... Procrustes: rmse 1.995309e-06  max resid 4.059865e-06 
    ## ... Similar to previous best
    ## Run 200 stress 0.2687489 
    ## Run 201 stress 0.1383681 
    ## Run 202 stress 0.1410448 
    ## Run 203 stress 0.1407298 
    ## Run 204 stress 0.2501066 
    ## Run 205 stress 0.1407298 
    ## Run 206 stress 0.1407208 
    ## Run 207 stress 0.1320258 
    ## ... Procrustes: rmse 5.638372e-07  max resid 8.019729e-07 
    ## ... Similar to previous best
    ## Run 208 stress 0.1415299 
    ## Run 209 stress 0.1445811 
    ## Run 210 stress 0.1572668 
    ## Run 211 stress 0.1320258 
    ## ... Procrustes: rmse 1.139914e-06  max resid 2.340755e-06 
    ## ... Similar to previous best
    ## Run 212 stress 0.1587623 
    ## Run 213 stress 0.1628606 
    ## Run 214 stress 0.1407207 
    ## Run 215 stress 0.1407298 
    ## Run 216 stress 0.1572668 
    ## Run 217 stress 0.2384387 
    ## Run 218 stress 0.2578639 
    ## Run 219 stress 0.1383682 
    ## Run 220 stress 0.1445811 
    ## Run 221 stress 0.1693717 
    ## Run 222 stress 0.1551863 
    ## Run 223 stress 0.1320258 
    ## ... Procrustes: rmse 2.561061e-06  max resid 5.236729e-06 
    ## ... Similar to previous best
    ## Run 224 stress 0.1410448 
    ## Run 225 stress 0.1383681 
    ## Run 226 stress 0.1551863 
    ## Run 227 stress 0.1383681 
    ## Run 228 stress 0.1693717 
    ## Run 229 stress 0.1445811 
    ## Run 230 stress 0.2085297 
    ## Run 231 stress 0.1407298 
    ## Run 232 stress 0.1320258 
    ## ... Procrustes: rmse 1.47124e-06  max resid 3.376248e-06 
    ## ... Similar to previous best
    ## Run 233 stress 0.1572668 
    ## Run 234 stress 0.1551863 
    ## Run 235 stress 0.1415299 
    ## Run 236 stress 0.1771969 
    ## Run 237 stress 0.2612469 
    ## Run 238 stress 0.1320258 
    ## ... Procrustes: rmse 5.042232e-07  max resid 8.164139e-07 
    ## ... Similar to previous best
    ## Run 239 stress 0.1410449 
    ## Run 240 stress 0.1320258 
    ## ... Procrustes: rmse 1.51144e-06  max resid 3.029452e-06 
    ## ... Similar to previous best
    ## Run 241 stress 0.2510079 
    ## Run 242 stress 0.2085297 
    ## Run 243 stress 0.2385029 
    ## Run 244 stress 0.1407298 
    ## Run 245 stress 0.2008159 
    ## Run 246 stress 0.2015531 
    ## Run 247 stress 0.1796864 
    ## Run 248 stress 0.1410449 
    ## Run 249 stress 0.1445811 
    ## Run 250 stress 0.1410451 
    ## Run 251 stress 0.1320258 
    ## ... Procrustes: rmse 9.479235e-07  max resid 1.374061e-06 
    ## ... Similar to previous best
    ## Run 252 stress 0.1663544 
    ## Run 253 stress 0.1383682 
    ## Run 254 stress 0.1410448 
    ## Run 255 stress 0.1962476 
    ## Run 256 stress 0.1410448 
    ## Run 257 stress 0.1383681 
    ## Run 258 stress 0.1445811 
    ## Run 259 stress 0.1410448 
    ## Run 260 stress 0.1320258 
    ## ... Procrustes: rmse 1.430298e-06  max resid 2.918066e-06 
    ## ... Similar to previous best
    ## Run 261 stress 0.1830332 
    ## Run 262 stress 0.1320258 
    ## ... Procrustes: rmse 1.167794e-06  max resid 2.318241e-06 
    ## ... Similar to previous best
    ## Run 263 stress 0.1407207 
    ## Run 264 stress 0.1320258 
    ## ... Procrustes: rmse 1.89032e-06  max resid 3.814594e-06 
    ## ... Similar to previous best
    ## Run 265 stress 0.1628606 
    ## Run 266 stress 0.2385029 
    ## Run 267 stress 0.1410451 
    ## Run 268 stress 0.1383682 
    ## Run 269 stress 0.2416997 
    ## Run 270 stress 0.1320258 
    ## ... Procrustes: rmse 2.588045e-06  max resid 5.196281e-06 
    ## ... Similar to previous best
    ## Run 271 stress 0.1998955 
    ## Run 272 stress 0.1407298 
    ## Run 273 stress 0.1663544 
    ## Run 274 stress 0.1410448 
    ## Run 275 stress 0.1962476 
    ## Run 276 stress 0.1830332 
    ## Run 277 stress 0.1407207 
    ## Run 278 stress 0.1693717 
    ## Run 279 stress 0.1771966 
    ## Run 280 stress 0.1320258 
    ## ... Procrustes: rmse 4.562834e-07  max resid 9.176401e-07 
    ## ... Similar to previous best
    ## Run 281 stress 0.1320258 
    ## ... Procrustes: rmse 3.139869e-07  max resid 3.949518e-07 
    ## ... Similar to previous best
    ## Run 282 stress 0.1407298 
    ## Run 283 stress 0.1383681 
    ## Run 284 stress 0.1320258 
    ## ... Procrustes: rmse 1.998478e-06  max resid 4.0127e-06 
    ## ... Similar to previous best
    ## Run 285 stress 0.1415299 
    ## Run 286 stress 0.1407298 
    ## Run 287 stress 0.1583016 
    ## Run 288 stress 0.1410448 
    ## Run 289 stress 0.1383681 
    ## Run 290 stress 0.1970303 
    ## Run 291 stress 0.1415299 
    ## Run 292 stress 0.2227526 
    ## Run 293 stress 0.1410449 
    ## Run 294 stress 0.1802752 
    ## Run 295 stress 0.1320258 
    ## ... Procrustes: rmse 2.520158e-06  max resid 5.449145e-06 
    ## ... Similar to previous best
    ## Run 296 stress 0.1410448 
    ## Run 297 stress 0.1771969 
    ## Run 298 stress 0.1551863 
    ## Run 299 stress 0.2842805 
    ## Run 300 stress 0.1320258 
    ## ... Procrustes: rmse 1.489185e-06  max resid 3.037062e-06 
    ## ... Similar to previous best
    ## Run 301 stress 0.1415299 
    ## Run 302 stress 0.1415299 
    ## Run 303 stress 0.2466649 
    ## Run 304 stress 0.1415299 
    ## Run 305 stress 0.1407298 
    ## Run 306 stress 0.1771966 
    ## Run 307 stress 0.1410448 
    ## Run 308 stress 0.1410452 
    ## Run 309 stress 0.2385029 
    ## Run 310 stress 0.1572668 
    ## Run 311 stress 0.2085297 
    ## Run 312 stress 0.1407208 
    ## Run 313 stress 0.1583016 
    ## Run 314 stress 0.1445811 
    ## Run 315 stress 0.1802587 
    ## Run 316 stress 0.1551863 
    ## Run 317 stress 0.1663544 
    ## Run 318 stress 0.1796864 
    ## Run 319 stress 0.1572668 
    ## Run 320 stress 0.1410449 
    ## Run 321 stress 0.1320258 
    ## ... Procrustes: rmse 6.336499e-07  max resid 1.114025e-06 
    ## ... Similar to previous best
    ## Run 322 stress 0.1581642 
    ## Run 323 stress 0.1415299 
    ## Run 324 stress 0.1445811 
    ## Run 325 stress 0.2422742 
    ## Run 326 stress 0.1320258 
    ## ... Procrustes: rmse 1.878965e-06  max resid 3.753382e-06 
    ## ... Similar to previous best
    ## Run 327 stress 0.141045 
    ## Run 328 stress 0.1771966 
    ## Run 329 stress 0.1776793 
    ## Run 330 stress 0.2385029 
    ## Run 331 stress 0.2227527 
    ## Run 332 stress 0.1320258 
    ## ... Procrustes: rmse 4.014335e-07  max resid 6.878396e-07 
    ## ... Similar to previous best
    ## Run 333 stress 0.1320258 
    ## ... Procrustes: rmse 1.59943e-06  max resid 3.247053e-06 
    ## ... Similar to previous best
    ## Run 334 stress 0.1415299 
    ## Run 335 stress 0.1407298 
    ## Run 336 stress 0.1415299 
    ## Run 337 stress 0.2008162 
    ## Run 338 stress 0.2224298 
    ## Run 339 stress 0.1407298 
    ## Run 340 stress 0.1410448 
    ## Run 341 stress 0.2384387 
    ## Run 342 stress 0.1830332 
    ## Run 343 stress 0.1320258 
    ## ... Procrustes: rmse 1.300398e-06  max resid 1.867571e-06 
    ## ... Similar to previous best
    ## Run 344 stress 0.2224298 
    ## Run 345 stress 0.1583016 
    ## Run 346 stress 0.1410452 
    ## Run 347 stress 0.1320258 
    ## ... Procrustes: rmse 3.543141e-07  max resid 6.004692e-07 
    ## ... Similar to previous best
    ## Run 348 stress 0.1383681 
    ## Run 349 stress 0.2227526 
    ## Run 350 stress 0.1415299 
    ## Run 351 stress 0.1970303 
    ## Run 352 stress 0.1434197 
    ## Run 353 stress 0.1407298 
    ## Run 354 stress 0.2550045 
    ## Run 355 stress 0.2496382 
    ## Run 356 stress 0.1383682 
    ## Run 357 stress 0.1415299 
    ## Run 358 stress 0.2385029 
    ## Run 359 stress 0.1320258 
    ## ... Procrustes: rmse 2.348984e-06  max resid 4.774118e-06 
    ## ... Similar to previous best
    ## Run 360 stress 0.1320258 
    ## ... Procrustes: rmse 1.90985e-06  max resid 3.86944e-06 
    ## ... Similar to previous best
    ## Run 361 stress 0.1407298 
    ## Run 362 stress 0.1320258 
    ## ... Procrustes: rmse 1.328963e-06  max resid 2.112227e-06 
    ## ... Similar to previous best
    ## Run 363 stress 0.2384387 
    ## Run 364 stress 0.1572668 
    ## Run 365 stress 0.1802752 
    ## Run 366 stress 0.1410449 
    ## Run 367 stress 0.1693717 
    ## Run 368 stress 0.1383682 
    ## Run 369 stress 0.2550044 
    ## Run 370 stress 0.1383681 
    ## Run 371 stress 0.1628606 
    ## Run 372 stress 0.1572668 
    ## Run 373 stress 0.1583016 
    ## Run 374 stress 0.2422742 
    ## Run 375 stress 0.1551863 
    ## Run 376 stress 0.1551863 
    ## Run 377 stress 0.1598154 
    ## Run 378 stress 0.1320258 
    ## ... Procrustes: rmse 1.684507e-06  max resid 3.263744e-06 
    ## ... Similar to previous best
    ## Run 379 stress 0.1407298 
    ## Run 380 stress 0.1572668 
    ## Run 381 stress 0.1434197 
    ## Run 382 stress 0.1320258 
    ## ... Procrustes: rmse 8.708897e-07  max resid 1.756366e-06 
    ## ... Similar to previous best
    ## Run 383 stress 0.2224298 
    ## Run 384 stress 0.1320258 
    ## ... Procrustes: rmse 1.265919e-06  max resid 2.822299e-06 
    ## ... Similar to previous best
    ## Run 385 stress 0.1572668 
    ## Run 386 stress 0.1830332 
    ## Run 387 stress 0.3083098 
    ## Run 388 stress 0.1415299 
    ## Run 389 stress 0.1551863 
    ## Run 390 stress 0.1929519 
    ## Run 391 stress 0.2384387 
    ## Run 392 stress 0.1415299 
    ## Run 393 stress 0.1693717 
    ## Run 394 stress 0.1572668 
    ## Run 395 stress 0.1320258 
    ## ... Procrustes: rmse 2.849679e-06  max resid 5.830202e-06 
    ## ... Similar to previous best
    ## Run 396 stress 0.1445811 
    ## Run 397 stress 0.1415299 
    ## Run 398 stress 0.1410453 
    ## Run 399 stress 0.141045 
    ## Run 400 stress 0.1410452 
    ## Run 401 stress 0.1572668 
    ## Run 402 stress 0.1445811 
    ## Run 403 stress 0.1407206 
    ## Run 404 stress 0.1663544 
    ## Run 405 stress 0.1572668 
    ## Run 406 stress 0.1445811 
    ## Run 407 stress 0.1962476 
    ## Run 408 stress 0.1587623 
    ## Run 409 stress 0.1572668 
    ## Run 410 stress 0.1551863 
    ## Run 411 stress 0.1410448 
    ## Run 412 stress 0.1962476 
    ## Run 413 stress 0.285216 
    ## Run 414 stress 0.1410451 
    ## Run 415 stress 0.1407298 
    ## Run 416 stress 0.1383681 
    ## Run 417 stress 0.1320258 
    ## ... Procrustes: rmse 1.082182e-06  max resid 2.238723e-06 
    ## ... Similar to previous best
    ## Run 418 stress 0.1407298 
    ## Run 419 stress 0.1320258 
    ## ... Procrustes: rmse 5.473008e-07  max resid 9.044049e-07 
    ## ... Similar to previous best
    ## Run 420 stress 0.1551863 
    ## Run 421 stress 0.1583016 
    ## Run 422 stress 0.1551863 
    ## Run 423 stress 0.2008162 
    ## Run 424 stress 0.1771969 
    ## Run 425 stress 0.1320258 
    ## ... Procrustes: rmse 2.86706e-06  max resid 5.901271e-06 
    ## ... Similar to previous best
    ## Run 426 stress 0.1383682 
    ## Run 427 stress 0.1407298 
    ## Run 428 stress 0.1551863 
    ## Run 429 stress 0.3083098 
    ## Run 430 stress 0.1320258 
    ## ... Procrustes: rmse 5.574508e-07  max resid 1.080083e-06 
    ## ... Similar to previous best
    ## Run 431 stress 0.1776793 
    ## Run 432 stress 0.1320258 
    ## ... Procrustes: rmse 2.010042e-06  max resid 4.125946e-06 
    ## ... Similar to previous best
    ## Run 433 stress 0.1410448 
    ## Run 434 stress 0.1929519 
    ## Run 435 stress 0.1320258 
    ## ... Procrustes: rmse 9.083446e-07  max resid 1.7622e-06 
    ## ... Similar to previous best
    ## Run 436 stress 0.1693717 
    ## Run 437 stress 0.1320258 
    ## ... Procrustes: rmse 1.343648e-06  max resid 2.03506e-06 
    ## ... Similar to previous best
    ## Run 438 stress 0.1434197 
    ## Run 439 stress 0.1663544 
    ## Run 440 stress 0.2842805 
    ## Run 441 stress 0.2456775 
    ## Run 442 stress 0.1415299 
    ## Run 443 stress 0.1771969 
    ## Run 444 stress 0.1583016 
    ## Run 445 stress 0.2085297 
    ## Run 446 stress 0.1581642 
    ## Run 447 stress 0.1320258 
    ## ... Procrustes: rmse 8.400344e-07  max resid 1.701941e-06 
    ## ... Similar to previous best
    ## Run 448 stress 0.1581642 
    ## Run 449 stress 0.1410452 
    ## Run 450 stress 0.1407298 
    ## Run 451 stress 0.1415299 
    ## Run 452 stress 0.1320258 
    ## ... Procrustes: rmse 8.331002e-07  max resid 1.604264e-06 
    ## ... Similar to previous best
    ## Run 453 stress 0.1551863 
    ## Run 454 stress 0.1572668 
    ## Run 455 stress 0.1383682 
    ## Run 456 stress 0.1320258 
    ## ... Procrustes: rmse 1.812406e-06  max resid 3.671578e-06 
    ## ... Similar to previous best
    ## Run 457 stress 0.2224298 
    ## Run 458 stress 0.1998955 
    ## Run 459 stress 0.1693717 
    ## Run 460 stress 0.1407298 
    ## Run 461 stress 0.1383681 
    ## Run 462 stress 0.1970303 
    ## Run 463 stress 0.1383682 
    ## Run 464 stress 0.1320258 
    ## ... Procrustes: rmse 3.965291e-07  max resid 6.528904e-07 
    ## ... Similar to previous best
    ## Run 465 stress 0.1407298 
    ## Run 466 stress 0.1572668 
    ## Run 467 stress 0.1320258 
    ## ... Procrustes: rmse 1.068009e-06  max resid 2.10221e-06 
    ## ... Similar to previous best
    ## Run 468 stress 0.1320258 
    ## ... Procrustes: rmse 3.869284e-07  max resid 7.080727e-07 
    ## ... Similar to previous best
    ## Run 469 stress 0.1693717 
    ## Run 470 stress 0.1407298 
    ## Run 471 stress 0.1572668 
    ## Run 472 stress 0.1583016 
    ## Run 473 stress 0.1551863 
    ## Run 474 stress 0.1693717 
    ## Run 475 stress 0.1551863 
    ## Run 476 stress 0.1410449 
    ## Run 477 stress 0.1410449 
    ## Run 478 stress 0.1407298 
    ## Run 479 stress 0.1830332 
    ## Run 480 stress 0.2384387 
    ## Run 481 stress 0.1320258 
    ## ... Procrustes: rmse 1.207063e-06  max resid 2.742441e-06 
    ## ... Similar to previous best
    ## Run 482 stress 0.2015531 
    ## Run 483 stress 0.1407298 
    ## Run 484 stress 0.1970303 
    ## Run 485 stress 0.2008159 
    ## Run 486 stress 0.1383682 
    ## Run 487 stress 0.1771969 
    ## Run 488 stress 0.1830332 
    ## Run 489 stress 0.1445811 
    ## Run 490 stress 0.1628606 
    ## Run 491 stress 0.1320258 
    ## ... Procrustes: rmse 5.011312e-07  max resid 6.699313e-07 
    ## ... Similar to previous best
    ## Run 492 stress 0.1587623 
    ## Run 493 stress 0.1663544 
    ## Run 494 stress 0.1572668 
    ## Run 495 stress 0.1572668 
    ## Run 496 stress 0.1663545 
    ## Run 497 stress 0.1929519 
    ## Run 498 stress 0.1407207 
    ## Run 499 stress 0.1771966 
    ## Run 500 stress 0.1415299 
    ## *** Best solution repeated 51 times

``` r
### Environmental
# Surveyed sites 
SD_beta_env_NMDS <- metaMDS(SD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07669178 
    ## Run 1 stress 0.08314399 
    ## Run 2 stress 0.08097769 
    ## Run 3 stress 0.07951482 
    ## Run 4 stress 0.0809777 
    ## Run 5 stress 0.08252129 
    ## Run 6 stress 0.08717125 
    ## Run 7 stress 0.08182115 
    ## Run 8 stress 0.08175702 
    ## Run 9 stress 0.08175694 
    ## Run 10 stress 0.08175691 
    ## Run 11 stress 0.08471586 
    ## Run 12 stress 0.0895756 
    ## Run 13 stress 0.08175691 
    ## Run 14 stress 0.0826545 
    ## Run 15 stress 0.08364079 
    ## Run 16 stress 0.08265415 
    ## Run 17 stress 0.08421026 
    ## Run 18 stress 0.084337 
    ## Run 19 stress 0.08175705 
    ## Run 20 stress 0.08499304 
    ## Run 21 stress 0.08314374 
    ## Run 22 stress 0.0818938 
    ## Run 23 stress 0.0821734 
    ## Run 24 stress 0.08175699 
    ## Run 25 stress 0.08337164 
    ## Run 26 stress 0.08360265 
    ## Run 27 stress 0.08252147 
    ## Run 28 stress 0.08470949 
    ## Run 29 stress 0.07669164 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0005407989  max resid 0.0009354444 
    ## ... Similar to previous best
    ## Run 30 stress 0.07948997 
    ## Run 31 stress 0.08405262 
    ## Run 32 stress 0.08243375 
    ## Run 33 stress 0.08835173 
    ## Run 34 stress 0.08003174 
    ## Run 35 stress 0.07960691 
    ## Run 36 stress 0.0821852 
    ## Run 37 stress 0.07943329 
    ## Run 38 stress 0.08175698 
    ## Run 39 stress 0.08249602 
    ## Run 40 stress 0.08305449 
    ## Run 41 stress 0.08665277 
    ## Run 42 stress 0.08175691 
    ## Run 43 stress 0.0794901 
    ## Run 44 stress 0.08175697 
    ## Run 45 stress 0.08175708 
    ## Run 46 stress 0.08283568 
    ## Run 47 stress 0.08217346 
    ## Run 48 stress 0.08252118 
    ## Run 49 stress 0.07868355 
    ## Run 50 stress 0.07951477 
    ## Run 51 stress 0.08175712 
    ## Run 52 stress 0.0842103 
    ## Run 53 stress 0.08530163 
    ## Run 54 stress 0.08250126 
    ## Run 55 stress 0.08674793 
    ## Run 56 stress 0.08175708 
    ## Run 57 stress 0.08189398 
    ## Run 58 stress 0.0836027 
    ## Run 59 stress 0.08182125 
    ## Run 60 stress 0.07960667 
    ## Run 61 stress 0.07871729 
    ## Run 62 stress 0.08104907 
    ## Run 63 stress 0.08887949 
    ## Run 64 stress 0.08355837 
    ## Run 65 stress 0.08317315 
    ## Run 66 stress 0.0837577 
    ## Run 67 stress 0.08175701 
    ## Run 68 stress 0.08097762 
    ## Run 69 stress 0.08477409 
    ## Run 70 stress 0.07871749 
    ## Run 71 stress 0.08305432 
    ## Run 72 stress 0.08684085 
    ## Run 73 stress 0.0794334 
    ## Run 74 stress 0.0800319 
    ## Run 75 stress 0.08527051 
    ## Run 76 stress 0.08433685 
    ## Run 77 stress 0.08003201 
    ## Run 78 stress 0.08507001 
    ## Run 79 stress 0.08375724 
    ## Run 80 stress 0.08217345 
    ## Run 81 stress 0.07960657 
    ## Run 82 stress 0.08265429 
    ## Run 83 stress 0.08175695 
    ## Run 84 stress 0.08317312 
    ## Run 85 stress 0.08527024 
    ## Run 86 stress 0.081757 
    ## Run 87 stress 0.08337156 
    ## Run 88 stress 0.08355859 
    ## Run 89 stress 0.08317353 
    ## Run 90 stress 0.08175709 
    ## Run 91 stress 0.08252123 
    ## Run 92 stress 0.08175699 
    ## Run 93 stress 0.08433696 
    ## Run 94 stress 0.08104916 
    ## Run 95 stress 0.08957549 
    ## Run 96 stress 0.07948995 
    ## Run 97 stress 0.07951488 
    ## Run 98 stress 0.08375774 
    ## Run 99 stress 0.08317371 
    ## Run 100 stress 0.08577475 
    ## Run 101 stress 0.08233204 
    ## Run 102 stress 0.07960658 
    ## Run 103 stress 0.08175699 
    ## Run 104 stress 0.07731515 
    ## Run 105 stress 0.08097762 
    ## Run 106 stress 0.08375756 
    ## Run 107 stress 0.08104929 
    ## Run 108 stress 0.08547268 
    ## Run 109 stress 0.08835139 
    ## Run 110 stress 0.08217355 
    ## Run 111 stress 0.08182122 
    ## Run 112 stress 0.08189406 
    ## Run 113 stress 0.07669158 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001826432  max resid 0.0003995064 
    ## ... Similar to previous best
    ## Run 114 stress 0.08012958 
    ## Run 115 stress 0.07669152 
    ## ... New best solution
    ## ... Procrustes: rmse 7.20453e-05  max resid 0.0001247798 
    ## ... Similar to previous best
    ## Run 116 stress 0.08337197 
    ## Run 117 stress 0.07948999 
    ## Run 118 stress 0.08305189 
    ## Run 119 stress 0.08584987 
    ## Run 120 stress 0.07871759 
    ## Run 121 stress 0.08104927 
    ## Run 122 stress 0.07949022 
    ## Run 123 stress 0.08013384 
    ## Run 124 stress 0.07951479 
    ## Run 125 stress 0.08464766 
    ## Run 126 stress 0.08144925 
    ## Run 127 stress 0.08433687 
    ## Run 128 stress 0.08182138 
    ## Run 129 stress 0.07871751 
    ## Run 130 stress 0.08433693 
    ## Run 131 stress 0.08175698 
    ## Run 132 stress 0.08314412 
    ## Run 133 stress 0.08012983 
    ## Run 134 stress 0.08243374 
    ## Run 135 stress 0.07731517 
    ## Run 136 stress 0.08012963 
    ## Run 137 stress 0.08012978 
    ## Run 138 stress 0.08104886 
    ## Run 139 stress 0.08835135 
    ## Run 140 stress 0.08189381 
    ## Run 141 stress 0.07871742 
    ## Run 142 stress 0.08360267 
    ## Run 143 stress 0.07960696 
    ## Run 144 stress 0.08249619 
    ## Run 145 stress 0.08433684 
    ## Run 146 stress 0.08218537 
    ## Run 147 stress 0.0817569 
    ## Run 148 stress 0.09012098 
    ## Run 149 stress 0.08459681 
    ## Run 150 stress 0.08274934 
    ## Run 151 stress 0.08243365 
    ## Run 152 stress 0.09070243 
    ## Run 153 stress 0.08684145 
    ## Run 154 stress 0.08265443 
    ## Run 155 stress 0.08433685 
    ## Run 156 stress 0.08433686 
    ## Run 157 stress 0.08012994 
    ## Run 158 stress 0.07951496 
    ## Run 159 stress 0.08527024 
    ## Run 160 stress 0.08241486 
    ## Run 161 stress 0.07669166 
    ## ... Procrustes: rmse 0.0002182184  max resid 0.0004369731 
    ## ... Similar to previous best
    ## Run 162 stress 0.08097759 
    ## Run 163 stress 0.08370302 
    ## Run 164 stress 0.08265418 
    ## Run 165 stress 0.0825212 
    ## Run 166 stress 0.08464777 
    ## Run 167 stress 0.07669154 
    ## ... Procrustes: rmse 0.0001610044  max resid 0.0002818048 
    ## ... Similar to previous best
    ## Run 168 stress 0.08252144 
    ## Run 169 stress 0.08283586 
    ## Run 170 stress 0.07669157 
    ## ... Procrustes: rmse 0.0001894537  max resid 0.0003209339 
    ## ... Similar to previous best
    ## Run 171 stress 0.07731517 
    ## Run 172 stress 0.07960693 
    ## Run 173 stress 0.07731531 
    ## Run 174 stress 0.08421009 
    ## Run 175 stress 0.07669169 
    ## ... Procrustes: rmse 0.0002583926  max resid 0.0004412167 
    ## ... Similar to previous best
    ## Run 176 stress 0.07943368 
    ## Run 177 stress 0.08510826 
    ## Run 178 stress 0.08217355 
    ## Run 179 stress 0.08283578 
    ## Run 180 stress 0.07731519 
    ## Run 181 stress 0.08250142 
    ## Run 182 stress 0.08337205 
    ## Run 183 stress 0.08684191 
    ## Run 184 stress 0.08175713 
    ## Run 185 stress 0.08252119 
    ## Run 186 stress 0.08339526 
    ## Run 187 stress 0.07871749 
    ## Run 188 stress 0.08104905 
    ## Run 189 stress 0.08684099 
    ## Run 190 stress 0.08233237 
    ## Run 191 stress 0.08252123 
    ## Run 192 stress 0.08252147 
    ## Run 193 stress 0.08012974 
    ## Run 194 stress 0.08104883 
    ## Run 195 stress 0.0830546 
    ## Run 196 stress 0.08249618 
    ## Run 197 stress 0.08500356 
    ## Run 198 stress 0.08305336 
    ## Run 199 stress 0.08012975 
    ## Run 200 stress 0.08329768 
    ## Run 201 stress 0.08835166 
    ## Run 202 stress 0.08314416 
    ## Run 203 stress 0.08477411 
    ## Run 204 stress 0.08012963 
    ## Run 205 stress 0.08003217 
    ## Run 206 stress 0.08283545 
    ## Run 207 stress 0.08012946 
    ## Run 208 stress 0.08360278 
    ## Run 209 stress 0.08252118 
    ## Run 210 stress 0.08305535 
    ## Run 211 stress 0.08362752 
    ## Run 212 stress 0.07960681 
    ## Run 213 stress 0.0800323 
    ## Run 214 stress 0.08464781 
    ## Run 215 stress 0.08012944 
    ## Run 216 stress 0.08012952 
    ## Run 217 stress 0.08929141 
    ## Run 218 stress 0.07669166 
    ## ... Procrustes: rmse 0.0002304768  max resid 0.0003990043 
    ## ... Similar to previous best
    ## Run 219 stress 0.08175691 
    ## Run 220 stress 0.07949027 
    ## Run 221 stress 0.07868352 
    ## Run 222 stress 0.08433688 
    ## Run 223 stress 0.08252119 
    ## Run 224 stress 0.07871766 
    ## Run 225 stress 0.08182129 
    ## Run 226 stress 0.08684106 
    ## Run 227 stress 0.08175717 
    ## Run 228 stress 0.07871734 
    ## Run 229 stress 0.08265431 
    ## Run 230 stress 0.08265465 
    ## Run 231 stress 0.0895751 
    ## Run 232 stress 0.08218539 
    ## Run 233 stress 0.07948982 
    ## Run 234 stress 0.08283594 
    ## Run 235 stress 0.08175691 
    ## Run 236 stress 0.07871745 
    ## Run 237 stress 0.07871733 
    ## Run 238 stress 0.08500343 
    ## Run 239 stress 0.08252116 
    ## Run 240 stress 0.08182107 
    ## Run 241 stress 0.07951486 
    ## Run 242 stress 0.08471602 
    ## Run 243 stress 0.08443287 
    ## Run 244 stress 0.08283565 
    ## Run 245 stress 0.08233226 
    ## Run 246 stress 0.08243366 
    ## Run 247 stress 0.08957508 
    ## Run 248 stress 0.08360265 
    ## Run 249 stress 0.08249605 
    ## Run 250 stress 0.08283554 
    ## Run 251 stress 0.08485158 
    ## Run 252 stress 0.08305473 
    ## Run 253 stress 0.08493069 
    ## Run 254 stress 0.0794342 
    ## Run 255 stress 0.08175694 
    ## Run 256 stress 0.07871743 
    ## Run 257 stress 0.08013003 
    ## Run 258 stress 0.07669155 
    ## ... Procrustes: rmse 0.00014539  max resid 0.0003188872 
    ## ... Similar to previous best
    ## Run 259 stress 0.07871744 
    ## Run 260 stress 0.08218523 
    ## Run 261 stress 0.08144986 
    ## Run 262 stress 0.07960686 
    ## Run 263 stress 0.08421012 
    ## Run 264 stress 0.08175709 
    ## Run 265 stress 0.07731544 
    ## Run 266 stress 0.08835132 
    ## Run 267 stress 0.08012968 
    ## Run 268 stress 0.08104895 
    ## Run 269 stress 0.07868373 
    ## Run 270 stress 0.08003191 
    ## Run 271 stress 0.08507063 
    ## Run 272 stress 0.08673735 
    ## Run 273 stress 0.08104893 
    ## Run 274 stress 0.08265419 
    ## Run 275 stress 0.08329765 
    ## Run 276 stress 0.07731515 
    ## Run 277 stress 0.08249604 
    ## Run 278 stress 0.08104903 
    ## Run 279 stress 0.08898825 
    ## Run 280 stress 0.08433686 
    ## Run 281 stress 0.0838298 
    ## Run 282 stress 0.07669166 
    ## ... Procrustes: rmse 0.0002547076  max resid 0.00056954 
    ## ... Similar to previous best
    ## Run 283 stress 0.0801298 
    ## Run 284 stress 0.0809776 
    ## Run 285 stress 0.08252134 
    ## Run 286 stress 0.07669175 
    ## ... Procrustes: rmse 0.0003122367  max resid 0.0005420729 
    ## ... Similar to previous best
    ## Run 287 stress 0.08249599 
    ## Run 288 stress 0.0821736 
    ## Run 289 stress 0.08337161 
    ## Run 290 stress 0.07868384 
    ## Run 291 stress 0.07731524 
    ## Run 292 stress 0.07943369 
    ## Run 293 stress 0.08684114 
    ## Run 294 stress 0.08337145 
    ## Run 295 stress 0.07731513 
    ## Run 296 stress 0.07669182 
    ## ... Procrustes: rmse 0.0003588478  max resid 0.0006243537 
    ## ... Similar to previous best
    ## Run 297 stress 0.08243364 
    ## Run 298 stress 0.07669158 
    ## ... Procrustes: rmse 0.0001989717  max resid 0.0003394898 
    ## ... Similar to previous best
    ## Run 299 stress 0.08337202 
    ## Run 300 stress 0.08265411 
    ## Run 301 stress 0.07960664 
    ## Run 302 stress 0.08144874 
    ## Run 303 stress 0.0869213 
    ## Run 304 stress 0.08013691 
    ## Run 305 stress 0.08097748 
    ## Run 306 stress 0.07960669 
    ## Run 307 stress 0.07871747 
    ## Run 308 stress 0.08476969 
    ## Run 309 stress 0.08421022 
    ## Run 310 stress 0.07669157 
    ## ... Procrustes: rmse 0.0002040527  max resid 0.0003761106 
    ## ... Similar to previous best
    ## Run 311 stress 0.08175704 
    ## Run 312 stress 0.08283552 
    ## Run 313 stress 0.08305374 
    ## Run 314 stress 0.07669176 
    ## ... Procrustes: rmse 0.0002506005  max resid 0.0004980303 
    ## ... Similar to previous best
    ## Run 315 stress 0.08217339 
    ## Run 316 stress 0.07871731 
    ## Run 317 stress 0.07669169 
    ## ... Procrustes: rmse 0.0002867006  max resid 0.0004978888 
    ## ... Similar to previous best
    ## Run 318 stress 0.08305392 
    ## Run 319 stress 0.08527014 
    ## Run 320 stress 0.08379894 
    ## Run 321 stress 0.08144946 
    ## Run 322 stress 0.08305464 
    ## Run 323 stress 0.07948951 
    ## Run 324 stress 0.08265463 
    ## Run 325 stress 0.07731516 
    ## Run 326 stress 0.09103948 
    ## Run 327 stress 0.08012954 
    ## Run 328 stress 0.08317297 
    ## Run 329 stress 0.08003223 
    ## Run 330 stress 0.08553068 
    ## Run 331 stress 0.08317303 
    ## Run 332 stress 0.08477006 
    ## Run 333 stress 0.08175712 
    ## Run 334 stress 0.08243374 
    ## Run 335 stress 0.07868358 
    ## Run 336 stress 0.08144851 
    ## Run 337 stress 0.07871761 
    ## Run 338 stress 0.08250139 
    ## Run 339 stress 0.07669153 
    ## ... Procrustes: rmse 7.726245e-05  max resid 0.0001202171 
    ## ... Similar to previous best
    ## Run 340 stress 0.08003204 
    ## Run 341 stress 0.08097774 
    ## Run 342 stress 0.08097755 
    ## Run 343 stress 0.08835145 
    ## Run 344 stress 0.08360277 
    ## Run 345 stress 0.08355834 
    ## Run 346 stress 0.08471616 
    ## Run 347 stress 0.08433697 
    ## Run 348 stress 0.07949011 
    ## Run 349 stress 0.08144909 
    ## Run 350 stress 0.08360259 
    ## Run 351 stress 0.07669154 
    ## ... Procrustes: rmse 0.000134079  max resid 0.0002309275 
    ## ... Similar to previous best
    ## Run 352 stress 0.08375722 
    ## Run 353 stress 0.07960667 
    ## Run 354 stress 0.08684134 
    ## Run 355 stress 0.08499309 
    ## Run 356 stress 0.08210327 
    ## Run 357 stress 0.08686177 
    ## Run 358 stress 0.07669171 
    ## ... Procrustes: rmse 0.0002734482  max resid 0.0005473244 
    ## ... Similar to previous best
    ## Run 359 stress 0.0924631 
    ## Run 360 stress 0.08433706 
    ## Run 361 stress 0.07669159 
    ## ... Procrustes: rmse 0.0002238189  max resid 0.0003925121 
    ## ... Similar to previous best
    ## Run 362 stress 0.08682686 
    ## Run 363 stress 0.08175695 
    ## Run 364 stress 0.08012936 
    ## Run 365 stress 0.08433701 
    ## Run 366 stress 0.08898802 
    ## Run 367 stress 0.07669151 
    ## ... New best solution
    ## ... Procrustes: rmse 4.142013e-05  max resid 9.290179e-05 
    ## ... Similar to previous best
    ## Run 368 stress 0.07669177 
    ## ... Procrustes: rmse 0.0002486721  max resid 0.0004501762 
    ## ... Similar to previous best
    ## Run 369 stress 0.07960674 
    ## Run 370 stress 0.08012984 
    ## Run 371 stress 0.08097764 
    ## Run 372 stress 0.07731531 
    ## Run 373 stress 0.08175691 
    ## Run 374 stress 0.08968587 
    ## Run 375 stress 0.07871725 
    ## Run 376 stress 0.08682738 
    ## Run 377 stress 0.08175709 
    ## Run 378 stress 0.09024594 
    ## Run 379 stress 0.08470923 
    ## Run 380 stress 0.07669158 
    ## ... Procrustes: rmse 0.0001588195  max resid 0.0003502745 
    ## ... Similar to previous best
    ## Run 381 stress 0.08252118 
    ## Run 382 stress 0.07871735 
    ## Run 383 stress 0.08317298 
    ## Run 384 stress 0.08433685 
    ## Run 385 stress 0.08003212 
    ## Run 386 stress 0.08305274 
    ## Run 387 stress 0.07731532 
    ## Run 388 stress 0.08835187 
    ## Run 389 stress 0.09034816 
    ## Run 390 stress 0.08500377 
    ## Run 391 stress 0.08314428 
    ## Run 392 stress 0.07871726 
    ## Run 393 stress 0.07669174 
    ## ... Procrustes: rmse 0.0002613859  max resid 0.0004776688 
    ## ... Similar to previous best
    ## Run 394 stress 0.07669181 
    ## ... Procrustes: rmse 0.0002678303  max resid 0.0004628252 
    ## ... Similar to previous best
    ## Run 395 stress 0.07871726 
    ## Run 396 stress 0.08175694 
    ## Run 397 stress 0.08241518 
    ## Run 398 stress 0.08527045 
    ## Run 399 stress 0.09054952 
    ## Run 400 stress 0.08305325 
    ## Run 401 stress 0.08012947 
    ## Run 402 stress 0.0850146 
    ## Run 403 stress 0.08370311 
    ## Run 404 stress 0.08182122 
    ## Run 405 stress 0.08097773 
    ## Run 406 stress 0.08433697 
    ## Run 407 stress 0.0787173 
    ## Run 408 stress 0.0837572 
    ## Run 409 stress 0.07951485 
    ## Run 410 stress 0.08317348 
    ## Run 411 stress 0.08265471 
    ## Run 412 stress 0.08433685 
    ## Run 413 stress 0.08144791 
    ## Run 414 stress 0.07951475 
    ## Run 415 stress 0.08305396 
    ## Run 416 stress 0.08243393 
    ## Run 417 stress 0.08459632 
    ## Run 418 stress 0.07960666 
    ## Run 419 stress 0.08175708 
    ## Run 420 stress 0.08249615 
    ## Run 421 stress 0.08433694 
    ## Run 422 stress 0.08684093 
    ## Run 423 stress 0.08433684 
    ## Run 424 stress 0.07669151 
    ## ... Procrustes: rmse 1.043482e-05  max resid 1.766679e-05 
    ## ... Similar to previous best
    ## Run 425 stress 0.0796067 
    ## Run 426 stress 0.08305323 
    ## Run 427 stress 0.08252116 
    ## Run 428 stress 0.08013589 
    ## Run 429 stress 0.08684124 
    ## Run 430 stress 0.08241492 
    ## Run 431 stress 0.08335413 
    ## Run 432 stress 0.07948972 
    ## Run 433 stress 0.08144951 
    ## Run 434 stress 0.08968591 
    ## Run 435 stress 0.0831732 
    ## Run 436 stress 0.086921 
    ## Run 437 stress 0.08265488 
    ## Run 438 stress 0.08957577 
    ## Run 439 stress 0.08012966 
    ## Run 440 stress 0.0809774 
    ## Run 441 stress 0.07731518 
    ## Run 442 stress 0.08249606 
    ## Run 443 stress 0.079515 
    ## Run 444 stress 0.08582918 
    ## Run 445 stress 0.08317344 
    ## Run 446 stress 0.08182124 
    ## Run 447 stress 0.07943366 
    ## Run 448 stress 0.07948973 
    ## Run 449 stress 0.08405291 
    ## Run 450 stress 0.08097768 
    ## Run 451 stress 0.08175698 
    ## Run 452 stress 0.07868341 
    ## Run 453 stress 0.0892749 
    ## Run 454 stress 0.08189368 
    ## Run 455 stress 0.07669149 
    ## ... New best solution
    ## ... Procrustes: rmse 5.056272e-05  max resid 9.203004e-05 
    ## ... Similar to previous best
    ## Run 456 stress 0.08370295 
    ## Run 457 stress 0.08360267 
    ## Run 458 stress 0.08362752 
    ## Run 459 stress 0.08265378 
    ## Run 460 stress 0.0766917 
    ## ... Procrustes: rmse 0.0002745403  max resid 0.0004746661 
    ## ... Similar to previous best
    ## Run 461 stress 0.08175709 
    ## Run 462 stress 0.08548254 
    ## Run 463 stress 0.07951488 
    ## Run 464 stress 0.07669156 
    ## ... Procrustes: rmse 0.0001445529  max resid 0.0002553677 
    ## ... Similar to previous best
    ## Run 465 stress 0.08175688 
    ## Run 466 stress 0.08175695 
    ## Run 467 stress 0.08233294 
    ## Run 468 stress 0.08360262 
    ## Run 469 stress 0.08063923 
    ## Run 470 stress 0.07669184 
    ## ... Procrustes: rmse 0.0003491472  max resid 0.0006032514 
    ## ... Similar to previous best
    ## Run 471 stress 0.0858418 
    ## Run 472 stress 0.0826538 
    ## Run 473 stress 0.08012968 
    ## Run 474 stress 0.08283538 
    ## Run 475 stress 0.08584019 
    ## Run 476 stress 0.08305495 
    ## Run 477 stress 0.08314381 
    ## Run 478 stress 0.08250096 
    ## Run 479 stress 0.08243375 
    ## Run 480 stress 0.08175701 
    ## Run 481 stress 0.08317299 
    ## Run 482 stress 0.0817571 
    ## Run 483 stress 0.08104899 
    ## Run 484 stress 0.07731518 
    ## Run 485 stress 0.08421046 
    ## Run 486 stress 0.08527031 
    ## Run 487 stress 0.07669157 
    ## ... Procrustes: rmse 0.0001650719  max resid 0.0003113109 
    ## ... Similar to previous best
    ## Run 488 stress 0.080136 
    ## Run 489 stress 0.07871751 
    ## Run 490 stress 0.08097774 
    ## Run 491 stress 0.07868334 
    ## Run 492 stress 0.08317535 
    ## Run 493 stress 0.08684133 
    ## Run 494 stress 0.08217341 
    ## Run 495 stress 0.07669163 
    ## ... Procrustes: rmse 0.000226245  max resid 0.0004302992 
    ## ... Similar to previous best
    ## Run 496 stress 0.08314397 
    ## Run 497 stress 0.0801296 
    ## Run 498 stress 0.08013373 
    ## Run 499 stress 0.07871741 
    ## Run 500 stress 0.08553058 
    ## *** Best solution repeated 6 times

``` r
# Mixed and stratified lakes
SD_beta_env_MS_NMDS <- metaMDS(SD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09145321 
    ## Run 2 stress 0.08440261 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01444142  max resid 0.04302193 
    ## Run 3 stress 0.08503472 
    ## Run 4 stress 0.08503465 
    ## Run 5 stress 0.08973869 
    ## Run 6 stress 0.09669831 
    ## Run 7 stress 0.0944547 
    ## Run 8 stress 0.09610759 
    ## Run 9 stress 0.09030404 
    ## Run 10 stress 0.0897389 
    ## Run 11 stress 0.09969978 
    ## Run 12 stress 0.08503495 
    ## Run 13 stress 0.09030393 
    ## Run 14 stress 0.08503457 
    ## Run 15 stress 0.09030396 
    ## Run 16 stress 0.09407972 
    ## Run 17 stress 0.09308948 
    ## Run 18 stress 0.08773481 
    ## Run 19 stress 0.09721198 
    ## Run 20 stress 0.09539193 
    ## Run 21 stress 0.09145322 
    ## Run 22 stress 0.09168935 
    ## Run 23 stress 0.1072862 
    ## Run 24 stress 0.08973863 
    ## Run 25 stress 0.0940798 
    ## Run 26 stress 0.09447296 
    ## Run 27 stress 0.09416136 
    ## Run 28 stress 0.08973868 
    ## Run 29 stress 0.1026216 
    ## Run 30 stress 0.09321767 
    ## Run 31 stress 0.09030413 
    ## Run 32 stress 0.09308952 
    ## Run 33 stress 0.2465515 
    ## Run 34 stress 0.09381847 
    ## Run 35 stress 0.09610771 
    ## Run 36 stress 0.09145325 
    ## Run 37 stress 0.3042514 
    ## Run 38 stress 0.08503654 
    ## Run 39 stress 0.08503463 
    ## Run 40 stress 0.08973867 
    ## Run 41 stress 0.09268372 
    ## Run 42 stress 0.08503533 
    ## Run 43 stress 0.09407991 
    ## Run 44 stress 0.09286083 
    ## Run 45 stress 0.09416154 
    ## Run 46 stress 0.09407987 
    ## Run 47 stress 0.09145322 
    ## Run 48 stress 0.09407978 
    ## Run 49 stress 0.0877347 
    ## Run 50 stress 0.09286088 
    ## Run 51 stress 0.09159087 
    ## Run 52 stress 0.09535516 
    ## Run 53 stress 0.09159083 
    ## Run 54 stress 0.09465914 
    ## Run 55 stress 0.0916893 
    ## Run 56 stress 0.09539194 
    ## Run 57 stress 0.09308965 
    ## Run 58 stress 0.08503489 
    ## Run 59 stress 0.1038033 
    ## Run 60 stress 0.2951789 
    ## Run 61 stress 0.09286089 
    ## Run 62 stress 0.09374207 
    ## Run 63 stress 0.09030394 
    ## Run 64 stress 0.08440276 
    ## ... Procrustes: rmse 0.0002347623  max resid 0.0006046007 
    ## ... Similar to previous best
    ## Run 65 stress 0.09721199 
    ## Run 66 stress 0.09403413 
    ## Run 67 stress 0.09286086 
    ## Run 68 stress 0.08503484 
    ## Run 69 stress 0.09374309 
    ## Run 70 stress 0.0933723 
    ## Run 71 stress 0.08503602 
    ## Run 72 stress 0.09030395 
    ## Run 73 stress 0.08503579 
    ## Run 74 stress 0.08773468 
    ## Run 75 stress 0.09159085 
    ## Run 76 stress 0.09407375 
    ## Run 77 stress 0.08440254 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002509904  max resid 0.0005095744 
    ## ... Similar to previous best
    ## Run 78 stress 0.08503621 
    ## Run 79 stress 0.09408004 
    ## Run 80 stress 0.0903041 
    ## Run 81 stress 0.09337311 
    ## Run 82 stress 0.08503475 
    ## Run 83 stress 0.09712941 
    ## Run 84 stress 0.08503474 
    ## Run 85 stress 0.09030408 
    ## Run 86 stress 0.08503465 
    ## Run 87 stress 0.09030397 
    ## Run 88 stress 0.0997 
    ## Run 89 stress 0.1072864 
    ## Run 90 stress 0.09286087 
    ## Run 91 stress 0.08773478 
    ## Run 92 stress 0.0937425 
    ## Run 93 stress 0.09145326 
    ## Run 94 stress 0.0903041 
    ## Run 95 stress 0.09168935 
    ## Run 96 stress 0.09030413 
    ## Run 97 stress 0.09030397 
    ## Run 98 stress 0.09286086 
    ## Run 99 stress 0.09539196 
    ## Run 100 stress 0.08973876 
    ## Run 101 stress 0.08773469 
    ## Run 102 stress 0.09268361 
    ## Run 103 stress 0.0933727 
    ## Run 104 stress 0.09286087 
    ## Run 105 stress 0.09445468 
    ## Run 106 stress 0.08773479 
    ## Run 107 stress 0.09590942 
    ## Run 108 stress 0.09168926 
    ## Run 109 stress 0.08973886 
    ## Run 110 stress 0.09590983 
    ## Run 111 stress 0.08440254 
    ## ... Procrustes: rmse 5.962099e-05  max resid 0.0001620909 
    ## ... Similar to previous best
    ## Run 112 stress 0.09400539 
    ## Run 113 stress 0.09030408 
    ## Run 114 stress 0.09030394 
    ## Run 115 stress 0.08503675 
    ## Run 116 stress 0.08503642 
    ## Run 117 stress 0.09669838 
    ## Run 118 stress 0.08503629 
    ## Run 119 stress 0.08503631 
    ## Run 120 stress 0.09407975 
    ## Run 121 stress 0.09337254 
    ## Run 122 stress 0.09380582 
    ## Run 123 stress 0.09374187 
    ## Run 124 stress 0.08440254 
    ## ... New best solution
    ## ... Procrustes: rmse 5.773209e-05  max resid 0.000171308 
    ## ... Similar to previous best
    ## Run 125 stress 0.09465914 
    ## Run 126 stress 0.09721206 
    ## Run 127 stress 0.08503614 
    ## Run 128 stress 0.1038042 
    ## Run 129 stress 0.09159098 
    ## Run 130 stress 0.09445468 
    ## Run 131 stress 0.08973885 
    ## Run 132 stress 0.08440252 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001768126  max resid 0.0003583636 
    ## ... Similar to previous best
    ## Run 133 stress 0.09286086 
    ## Run 134 stress 0.0926836 
    ## Run 135 stress 0.09268351 
    ## Run 136 stress 0.091591 
    ## Run 137 stress 0.09337232 
    ## Run 138 stress 0.09380595 
    ## Run 139 stress 0.08773468 
    ## Run 140 stress 0.1095799 
    ## Run 141 stress 0.2510691 
    ## Run 142 stress 0.08503484 
    ## Run 143 stress 0.09535438 
    ## Run 144 stress 0.09168927 
    ## Run 145 stress 0.08973866 
    ## Run 146 stress 0.08503587 
    ## Run 147 stress 0.09445484 
    ## Run 148 stress 0.09286087 
    ## Run 149 stress 0.09416139 
    ## Run 150 stress 0.09030395 
    ## Run 151 stress 0.09030393 
    ## Run 152 stress 0.09145321 
    ## Run 153 stress 0.09159084 
    ## Run 154 stress 0.1026216 
    ## Run 155 stress 0.1005995 
    ## Run 156 stress 0.09030398 
    ## Run 157 stress 0.09286083 
    ## Run 158 stress 0.090304 
    ## Run 159 stress 0.09407969 
    ## Run 160 stress 0.09159086 
    ## Run 161 stress 0.08503565 
    ## Run 162 stress 0.08440274 
    ## ... Procrustes: rmse 0.0002122574  max resid 0.0003609541 
    ## ... Similar to previous best
    ## Run 163 stress 0.09374335 
    ## Run 164 stress 0.0941614 
    ## Run 165 stress 0.08503697 
    ## Run 166 stress 0.08503456 
    ## Run 167 stress 0.08503457 
    ## Run 168 stress 0.0850346 
    ## Run 169 stress 0.0940798 
    ## Run 170 stress 0.09337244 
    ## Run 171 stress 0.08503711 
    ## Run 172 stress 0.09465905 
    ## Run 173 stress 0.09268345 
    ## Run 174 stress 0.09286083 
    ## Run 175 stress 0.09721199 
    ## Run 176 stress 0.1038081 
    ## Run 177 stress 0.0933726 
    ## Run 178 stress 0.09030413 
    ## Run 179 stress 0.09969971 
    ## Run 180 stress 0.08773482 
    ## Run 181 stress 0.08503463 
    ## Run 182 stress 0.08440255 
    ## ... Procrustes: rmse 6.254763e-05  max resid 0.0001820849 
    ## ... Similar to previous best
    ## Run 183 stress 0.08503609 
    ## Run 184 stress 0.08773475 
    ## Run 185 stress 0.09465909 
    ## Run 186 stress 0.09286083 
    ## Run 187 stress 0.09145322 
    ## Run 188 stress 0.09447316 
    ## Run 189 stress 0.08503524 
    ## Run 190 stress 0.08440258 
    ## ... Procrustes: rmse 0.0002371711  max resid 0.0004752251 
    ## ... Similar to previous best
    ## Run 191 stress 0.08440268 
    ## ... Procrustes: rmse 0.0003068209  max resid 0.0006262713 
    ## ... Similar to previous best
    ## Run 192 stress 0.09535472 
    ## Run 193 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 8.722316e-05  max resid 0.0001851239 
    ## ... Similar to previous best
    ## Run 194 stress 0.09374298 
    ## Run 195 stress 0.08440258 
    ## ... Procrustes: rmse 0.0001311404  max resid 0.0002604449 
    ## ... Similar to previous best
    ## Run 196 stress 0.09308981 
    ## Run 197 stress 0.09465908 
    ## Run 198 stress 0.09268345 
    ## Run 199 stress 0.09590594 
    ## Run 200 stress 0.2992112 
    ## Run 201 stress 0.09030409 
    ## Run 202 stress 0.09407977 
    ## Run 203 stress 0.09268356 
    ## Run 204 stress 0.09145321 
    ## Run 205 stress 0.09407984 
    ## Run 206 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 1.797271e-05  max resid 3.427039e-05 
    ## ... Similar to previous best
    ## Run 207 stress 0.09286083 
    ## Run 208 stress 0.09969972 
    ## Run 209 stress 0.09400547 
    ## Run 210 stress 0.09374241 
    ## Run 211 stress 0.08503604 
    ## Run 212 stress 0.09268351 
    ## Run 213 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001636533  max resid 0.0003049351 
    ## ... Similar to previous best
    ## Run 214 stress 0.09535448 
    ## Run 215 stress 0.09380579 
    ## Run 216 stress 0.0915909 
    ## Run 217 stress 0.1104221 
    ## Run 218 stress 0.0926836 
    ## Run 219 stress 0.0903041 
    ## Run 220 stress 0.08440253 
    ## ... Procrustes: rmse 9.791063e-05  max resid 0.0001748291 
    ## ... Similar to previous best
    ## Run 221 stress 0.09030399 
    ## Run 222 stress 0.09030401 
    ## Run 223 stress 0.09464461 
    ## Run 224 stress 0.09374252 
    ## Run 225 stress 0.09407972 
    ## Run 226 stress 0.09286094 
    ## Run 227 stress 0.08773473 
    ## Run 228 stress 0.09407969 
    ## Run 229 stress 0.08773472 
    ## Run 230 stress 0.0972121 
    ## Run 231 stress 0.09400521 
    ## Run 232 stress 0.09407989 
    ## Run 233 stress 0.09159086 
    ## Run 234 stress 0.085035 
    ## Run 235 stress 0.09168955 
    ## Run 236 stress 0.09408002 
    ## Run 237 stress 0.09168947 
    ## Run 238 stress 0.09465909 
    ## Run 239 stress 0.09539192 
    ## Run 240 stress 0.3047131 
    ## Run 241 stress 0.09030406 
    ## Run 242 stress 0.100461 
    ## Run 243 stress 0.09030394 
    ## Run 244 stress 0.09268377 
    ## Run 245 stress 0.09416148 
    ## Run 246 stress 0.0946443 
    ## Run 247 stress 0.08503549 
    ## Run 248 stress 0.0926835 
    ## Run 249 stress 0.1038043 
    ## Run 250 stress 0.08503608 
    ## Run 251 stress 0.2483807 
    ## Run 252 stress 0.09145324 
    ## Run 253 stress 0.09268325 
    ## Run 254 stress 0.08773475 
    ## Run 255 stress 0.08973864 
    ## Run 256 stress 0.09145322 
    ## Run 257 stress 0.08503522 
    ## Run 258 stress 0.09168946 
    ## Run 259 stress 0.0940799 
    ## Run 260 stress 0.0926835 
    ## Run 261 stress 0.09908189 
    ## Run 262 stress 0.09030406 
    ## Run 263 stress 0.1038034 
    ## Run 264 stress 0.0959097 
    ## Run 265 stress 0.09416138 
    ## Run 266 stress 0.09416133 
    ## Run 267 stress 0.1038047 
    ## Run 268 stress 0.09145326 
    ## Run 269 stress 0.0940797 
    ## Run 270 stress 0.09337272 
    ## Run 271 stress 0.08973863 
    ## Run 272 stress 0.09400517 
    ## Run 273 stress 0.09268337 
    ## Run 274 stress 0.08773478 
    ## Run 275 stress 0.097212 
    ## Run 276 stress 0.09145323 
    ## Run 277 stress 0.09447326 
    ## Run 278 stress 0.09145324 
    ## Run 279 stress 0.09030409 
    ## Run 280 stress 0.09535535 
    ## Run 281 stress 0.09407128 
    ## Run 282 stress 0.09590956 
    ## Run 283 stress 0.09159085 
    ## Run 284 stress 0.09337258 
    ## Run 285 stress 0.3309718 
    ## Run 286 stress 0.08503561 
    ## Run 287 stress 0.08503602 
    ## Run 288 stress 0.103805 
    ## Run 289 stress 0.09030399 
    ## Run 290 stress 0.08973897 
    ## Run 291 stress 0.09268388 
    ## Run 292 stress 0.08440265 
    ## ... Procrustes: rmse 0.0002373881  max resid 0.0004383618 
    ## ... Similar to previous best
    ## Run 293 stress 0.1038044 
    ## Run 294 stress 0.09030402 
    ## Run 295 stress 0.09308956 
    ## Run 296 stress 0.09539193 
    ## Run 297 stress 0.08503472 
    ## Run 298 stress 0.09969968 
    ## Run 299 stress 0.09030394 
    ## Run 300 stress 0.09337309 
    ## Run 301 stress 0.09969978 
    ## Run 302 stress 0.09308946 
    ## Run 303 stress 0.09030397 
    ## Run 304 stress 0.09145321 
    ## Run 305 stress 0.09539194 
    ## Run 306 stress 0.08503504 
    ## Run 307 stress 0.09407976 
    ## Run 308 stress 0.09030406 
    ## Run 309 stress 0.09268371 
    ## Run 310 stress 0.09374267 
    ## Run 311 stress 0.09337253 
    ## Run 312 stress 0.08773479 
    ## Run 313 stress 0.09308952 
    ## Run 314 stress 0.09374363 
    ## Run 315 stress 0.09416152 
    ## Run 316 stress 0.09159086 
    ## Run 317 stress 0.09539194 
    ## Run 318 stress 0.09407978 
    ## Run 319 stress 0.09417791 
    ## Run 320 stress 0.08440258 
    ## ... Procrustes: rmse 0.0001752367  max resid 0.0003244839 
    ## ... Similar to previous best
    ## Run 321 stress 0.08973886 
    ## Run 322 stress 0.08503674 
    ## Run 323 stress 0.09407986 
    ## Run 324 stress 0.09416137 
    ## Run 325 stress 0.09145322 
    ## Run 326 stress 0.1052849 
    ## Run 327 stress 0.09159086 
    ## Run 328 stress 0.08503461 
    ## Run 329 stress 0.09374263 
    ## Run 330 stress 0.085035 
    ## Run 331 stress 0.09168929 
    ## Run 332 stress 0.09535513 
    ## Run 333 stress 0.08973864 
    ## Run 334 stress 0.09337238 
    ## Run 335 stress 0.09030402 
    ## Run 336 stress 0.1038043 
    ## Run 337 stress 0.09286085 
    ## Run 338 stress 0.09145321 
    ## Run 339 stress 0.08973873 
    ## Run 340 stress 0.09374402 
    ## Run 341 stress 0.09030398 
    ## Run 342 stress 0.09308963 
    ## Run 343 stress 0.09539208 
    ## Run 344 stress 0.1038038 
    ## Run 345 stress 0.08973863 
    ## Run 346 stress 0.09407966 
    ## Run 347 stress 0.09286089 
    ## Run 348 stress 0.09286087 
    ## Run 349 stress 0.09159088 
    ## Run 350 stress 0.09308966 
    ## Run 351 stress 0.1038044 
    ## Run 352 stress 0.09447309 
    ## Run 353 stress 0.08503509 
    ## Run 354 stress 0.100461 
    ## Run 355 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001249456  max resid 0.0002524271 
    ## ... Similar to previous best
    ## Run 356 stress 0.09407985 
    ## Run 357 stress 0.08440252 
    ## ... Procrustes: rmse 6.829473e-05  max resid 0.0001293812 
    ## ... Similar to previous best
    ## Run 358 stress 0.09030397 
    ## Run 359 stress 0.08973878 
    ## Run 360 stress 0.09030409 
    ## Run 361 stress 0.09286088 
    ## Run 362 stress 0.09465904 
    ## Run 363 stress 0.09465907 
    ## Run 364 stress 0.09337278 
    ## Run 365 stress 0.09465907 
    ## Run 366 stress 0.09464502 
    ## Run 367 stress 0.09445482 
    ## Run 368 stress 0.09535485 
    ## Run 369 stress 0.09969991 
    ## Run 370 stress 0.08503497 
    ## Run 371 stress 0.09535467 
    ## Run 372 stress 0.09321488 
    ## Run 373 stress 0.08503462 
    ## Run 374 stress 0.08973875 
    ## Run 375 stress 0.09416132 
    ## Run 376 stress 0.1038037 
    ## Run 377 stress 0.09465909 
    ## Run 378 stress 0.08973875 
    ## Run 379 stress 0.09465905 
    ## Run 380 stress 0.09407979 
    ## Run 381 stress 0.3411806 
    ## Run 382 stress 0.09159084 
    ## Run 383 stress 0.08503513 
    ## Run 384 stress 0.0926835 
    ## Run 385 stress 0.09168939 
    ## Run 386 stress 0.09286093 
    ## Run 387 stress 0.09268334 
    ## Run 388 stress 0.08973864 
    ## Run 389 stress 0.09407972 
    ## Run 390 stress 0.08440267 
    ## ... Procrustes: rmse 0.0002483704  max resid 0.0004543267 
    ## ... Similar to previous best
    ## Run 391 stress 0.1117251 
    ## Run 392 stress 0.09030395 
    ## Run 393 stress 0.08503462 
    ## Run 394 stress 0.1038047 
    ## Run 395 stress 0.09168929 
    ## Run 396 stress 0.08503456 
    ## Run 397 stress 0.0996999 
    ## Run 398 stress 0.08503504 
    ## Run 399 stress 0.09407466 
    ## Run 400 stress 0.09535475 
    ## Run 401 stress 0.09168954 
    ## Run 402 stress 0.09168927 
    ## Run 403 stress 0.09030394 
    ## Run 404 stress 0.09969998 
    ## Run 405 stress 0.09286086 
    ## Run 406 stress 0.09268359 
    ## Run 407 stress 0.09145321 
    ## Run 408 stress 0.09374281 
    ## Run 409 stress 0.09407977 
    ## Run 410 stress 0.09969974 
    ## Run 411 stress 0.2738609 
    ## Run 412 stress 0.09159084 
    ## Run 413 stress 0.09416149 
    ## Run 414 stress 0.0940797 
    ## Run 415 stress 0.09159091 
    ## Run 416 stress 0.08773483 
    ## Run 417 stress 0.09268349 
    ## Run 418 stress 0.09400518 
    ## Run 419 stress 0.09407115 
    ## Run 420 stress 0.09407117 
    ## Run 421 stress 0.09969989 
    ## Run 422 stress 0.09308969 
    ## Run 423 stress 0.09030395 
    ## Run 424 stress 0.09168966 
    ## Run 425 stress 0.08503561 
    ## Run 426 stress 0.09030395 
    ## Run 427 stress 0.09168945 
    ## Run 428 stress 0.09721197 
    ## Run 429 stress 0.09416136 
    ## Run 430 stress 0.08440262 
    ## ... Procrustes: rmse 0.0002107481  max resid 0.0003866485 
    ## ... Similar to previous best
    ## Run 431 stress 0.09403397 
    ## Run 432 stress 0.09308972 
    ## Run 433 stress 0.09337241 
    ## Run 434 stress 0.09168952 
    ## Run 435 stress 0.08503573 
    ## Run 436 stress 0.09030394 
    ## Run 437 stress 0.09145333 
    ## Run 438 stress 0.09286084 
    ## Run 439 stress 0.09030396 
    ## Run 440 stress 0.1026217 
    ## Run 441 stress 0.08773468 
    ## Run 442 stress 0.09407965 
    ## Run 443 stress 0.09539196 
    ## Run 444 stress 0.09590483 
    ## Run 445 stress 0.09374319 
    ## Run 446 stress 0.09321712 
    ## Run 447 stress 0.09268346 
    ## Run 448 stress 0.08773472 
    ## Run 449 stress 0.08503536 
    ## Run 450 stress 0.09539192 
    ## Run 451 stress 0.09159084 
    ## Run 452 stress 0.09337278 
    ## Run 453 stress 0.0903041 
    ## Run 454 stress 0.09407976 
    ## Run 455 stress 0.08773465 
    ## Run 456 stress 0.09539197 
    ## Run 457 stress 0.08440259 
    ## ... Procrustes: rmse 0.000181761  max resid 0.0003462501 
    ## ... Similar to previous best
    ## Run 458 stress 0.09268363 
    ## Run 459 stress 0.09374312 
    ## Run 460 stress 0.08773471 
    ## Run 461 stress 0.09030393 
    ## Run 462 stress 0.08440252 
    ## ... Procrustes: rmse 5.389927e-05  max resid 0.0001301378 
    ## ... Similar to previous best
    ## Run 463 stress 0.09539197 
    ## Run 464 stress 0.09381844 
    ## Run 465 stress 0.08503704 
    ## Run 466 stress 0.09030414 
    ## Run 467 stress 0.08503463 
    ## Run 468 stress 0.3271789 
    ## Run 469 stress 0.09159087 
    ## Run 470 stress 0.09407979 
    ## Run 471 stress 0.09465905 
    ## Run 472 stress 0.09381845 
    ## Run 473 stress 0.09407977 
    ## Run 474 stress 0.09760836 
    ## Run 475 stress 0.0897387 
    ## Run 476 stress 0.09168944 
    ## Run 477 stress 0.09590925 
    ## Run 478 stress 0.09168928 
    ## Run 479 stress 0.2465476 
    ## Run 480 stress 0.09408004 
    ## Run 481 stress 0.08773468 
    ## Run 482 stress 0.08503473 
    ## Run 483 stress 0.09969965 
    ## Run 484 stress 0.08503607 
    ## Run 485 stress 0.08503584 
    ## Run 486 stress 0.09417763 
    ## Run 487 stress 0.11245 
    ## Run 488 stress 0.08503541 
    ## Run 489 stress 0.09308981 
    ## Run 490 stress 0.09407994 
    ## Run 491 stress 0.09337229 
    ## Run 492 stress 0.09030394 
    ## Run 493 stress 0.0933725 
    ## Run 494 stress 0.09030409 
    ## Run 495 stress 0.09721201 
    ## Run 496 stress 0.08503539 
    ## Run 497 stress 0.09337297 
    ## Run 498 stress 0.08773466 
    ## Run 499 stress 0.09286093 
    ## Run 500 stress 0.09407967 
    ## *** Best solution repeated 11 times

``` r
# Ocean sites and mixed lakes
SD_beta_env_OM_NMDS <- metaMDS(SD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.06623458 
    ## Run 1 stress 0.06623458 
    ## ... Procrustes: rmse 0.0001000538  max resid 0.0001543287 
    ## ... Similar to previous best
    ## Run 2 stress 0.0647787 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1519232  max resid 0.264469 
    ## Run 3 stress 0.06477868 
    ## ... New best solution
    ## ... Procrustes: rmse 3.455486e-05  max resid 5.761473e-05 
    ## ... Similar to previous best
    ## Run 4 stress 0.06623459 
    ## Run 5 stress 0.06123361 
    ## ... New best solution
    ## ... Procrustes: rmse 0.08990843  max resid 0.2169337 
    ## Run 6 stress 0.06477964 
    ## Run 7 stress 0.07411951 
    ## Run 8 stress 0.06124496 
    ## ... Procrustes: rmse 0.01505088  max resid 0.04168424 
    ## Run 9 stress 0.06124516 
    ## ... Procrustes: rmse 0.01522659  max resid 0.042154 
    ## Run 10 stress 0.07166979 
    ## Run 11 stress 0.07411943 
    ## Run 12 stress 0.07411945 
    ## Run 13 stress 0.06477871 
    ## Run 14 stress 0.06477998 
    ## Run 15 stress 0.2381016 
    ## Run 16 stress 0.07166965 
    ## Run 17 stress 0.0741195 
    ## Run 18 stress 0.07166975 
    ## Run 19 stress 0.0647796 
    ## Run 20 stress 0.07411943 
    ## Run 21 stress 0.2381016 
    ## Run 22 stress 0.06123362 
    ## ... Procrustes: rmse 9.38309e-05  max resid 0.0001273189 
    ## ... Similar to previous best
    ## Run 23 stress 0.09623062 
    ## Run 24 stress 0.06623464 
    ## Run 25 stress 0.06477984 
    ## Run 26 stress 0.08226158 
    ## Run 27 stress 0.06477838 
    ## Run 28 stress 0.3270116 
    ## Run 29 stress 0.07411944 
    ## Run 30 stress 0.06623457 
    ## Run 31 stress 0.07411939 
    ## Run 32 stress 0.06477882 
    ## Run 33 stress 0.07411939 
    ## Run 34 stress 0.0612453 
    ## ... Procrustes: rmse 0.01424642  max resid 0.03952452 
    ## Run 35 stress 0.06623463 
    ## Run 36 stress 0.2345379 
    ## Run 37 stress 0.07166971 
    ## Run 38 stress 0.07411945 
    ## Run 39 stress 0.2970708 
    ## Run 40 stress 0.06623458 
    ## Run 41 stress 0.07166977 
    ## Run 42 stress 0.06477887 
    ## Run 43 stress 0.06623457 
    ## Run 44 stress 0.2361811 
    ## Run 45 stress 0.07411939 
    ## Run 46 stress 0.06477943 
    ## Run 47 stress 0.07166968 
    ## Run 48 stress 0.06623459 
    ## Run 49 stress 0.06477936 
    ## Run 50 stress 0.07166966 
    ## Run 51 stress 0.09093467 
    ## Run 52 stress 0.07166988 
    ## Run 53 stress 0.06623457 
    ## Run 54 stress 0.0612452 
    ## ... Procrustes: rmse 0.01431033  max resid 0.03969847 
    ## Run 55 stress 0.07411954 
    ## Run 56 stress 0.07166965 
    ## Run 57 stress 0.08226159 
    ## Run 58 stress 0.07166965 
    ## Run 59 stress 0.07166978 
    ## Run 60 stress 0.06124497 
    ## ... Procrustes: rmse 0.01505885  max resid 0.04170545 
    ## Run 61 stress 0.07411939 
    ## Run 62 stress 0.08226159 
    ## Run 63 stress 0.06623458 
    ## Run 64 stress 0.07411941 
    ## Run 65 stress 0.06623465 
    ## Run 66 stress 0.06124512 
    ## ... Procrustes: rmse 0.01519085  max resid 0.04205986 
    ## Run 67 stress 0.06623457 
    ## Run 68 stress 0.06477838 
    ## Run 69 stress 0.0741194 
    ## Run 70 stress 0.07411939 
    ## Run 71 stress 0.07411938 
    ## Run 72 stress 0.0662346 
    ## Run 73 stress 0.06124505 
    ## ... Procrustes: rmse 0.01513942  max resid 0.04192032 
    ## Run 74 stress 0.06123364 
    ## ... Procrustes: rmse 9.022462e-05  max resid 0.0001964346 
    ## ... Similar to previous best
    ## Run 75 stress 0.06623465 
    ## Run 76 stress 0.08226155 
    ## Run 77 stress 0.2381016 
    ## Run 78 stress 0.07411938 
    ## Run 79 stress 0.06123361 
    ## ... Procrustes: rmse 4.913917e-05  max resid 0.000117213 
    ## ... Similar to previous best
    ## Run 80 stress 0.07166975 
    ## Run 81 stress 0.07411938 
    ## Run 82 stress 0.08226172 
    ## Run 83 stress 0.06124489 
    ## ... Procrustes: rmse 0.0146305  max resid 0.04055742 
    ## Run 84 stress 0.06477961 
    ## Run 85 stress 0.07411945 
    ## Run 86 stress 0.06124518 
    ## ... Procrustes: rmse 0.01435422  max resid 0.03981281 
    ## Run 87 stress 0.06477853 
    ## Run 88 stress 0.06623464 
    ## Run 89 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 1.823652e-05  max resid 2.512356e-05 
    ## ... Similar to previous best
    ## Run 90 stress 0.07411944 
    ## Run 91 stress 0.08226157 
    ## Run 92 stress 0.07166978 
    ## Run 93 stress 0.07411944 
    ## Run 94 stress 0.06477931 
    ## Run 95 stress 0.07411949 
    ## Run 96 stress 0.06477938 
    ## Run 97 stress 0.06623458 
    ## Run 98 stress 0.06123363 
    ## ... Procrustes: rmse 4.917842e-05  max resid 6.609707e-05 
    ## ... Similar to previous best
    ## Run 99 stress 0.06124513 
    ## ... Procrustes: rmse 0.01518819  max resid 0.04205029 
    ## Run 100 stress 0.07166982 
    ## Run 101 stress 0.06623461 
    ## Run 102 stress 0.06124518 
    ## ... Procrustes: rmse 0.01432219  max resid 0.03972888 
    ## Run 103 stress 0.06477844 
    ## Run 104 stress 0.0612452 
    ## ... Procrustes: rmse 0.01523746  max resid 0.04218483 
    ## Run 105 stress 0.06477946 
    ## Run 106 stress 0.06623458 
    ## Run 107 stress 0.0612336 
    ## ... Procrustes: rmse 1.670118e-05  max resid 3.486815e-05 
    ## ... Similar to previous best
    ## Run 108 stress 0.08226158 
    ## Run 109 stress 0.07166979 
    ## Run 110 stress 0.08226162 
    ## Run 111 stress 0.08226161 
    ## Run 112 stress 0.06477874 
    ## Run 113 stress 0.07166966 
    ## Run 114 stress 0.06623457 
    ## Run 115 stress 0.07166966 
    ## Run 116 stress 0.06623458 
    ## Run 117 stress 0.06623458 
    ## Run 118 stress 0.07411942 
    ## Run 119 stress 0.07411948 
    ## Run 120 stress 0.07411942 
    ## Run 121 stress 0.08226175 
    ## Run 122 stress 0.06124501 
    ## ... Procrustes: rmse 0.01509411  max resid 0.04179904 
    ## Run 123 stress 0.06477965 
    ## Run 124 stress 0.0612451 
    ## ... Procrustes: rmse 0.0151577  max resid 0.04196819 
    ## Run 125 stress 0.06623457 
    ## Run 126 stress 0.07411952 
    ## Run 127 stress 0.06623457 
    ## Run 128 stress 0.07411938 
    ## Run 129 stress 0.06477942 
    ## Run 130 stress 0.06124499 
    ## ... Procrustes: rmse 0.01507567  max resid 0.04174977 
    ## Run 131 stress 0.0647791 
    ## Run 132 stress 0.07166987 
    ## Run 133 stress 0.06623458 
    ## Run 134 stress 0.06124498 
    ## ... Procrustes: rmse 0.01506243  max resid 0.04171421 
    ## Run 135 stress 0.06124513 
    ## ... Procrustes: rmse 0.0151978  max resid 0.04207678 
    ## Run 136 stress 0.07411938 
    ## Run 137 stress 0.2362366 
    ## Run 138 stress 0.06623457 
    ## Run 139 stress 0.06124528 
    ## ... Procrustes: rmse 0.01529681  max resid 0.04234144 
    ## Run 140 stress 0.06124501 
    ## ... Procrustes: rmse 0.01509247  max resid 0.04179456 
    ## Run 141 stress 0.06623459 
    ## Run 142 stress 0.07411941 
    ## Run 143 stress 0.07411948 
    ## Run 144 stress 0.07411945 
    ## Run 145 stress 0.2548998 
    ## Run 146 stress 0.07166974 
    ## Run 147 stress 0.0741194 
    ## Run 148 stress 0.06477931 
    ## Run 149 stress 0.2362366 
    ## Run 150 stress 0.08226157 
    ## Run 151 stress 0.07166967 
    ## Run 152 stress 0.09469048 
    ## Run 153 stress 0.06623465 
    ## Run 154 stress 0.07411948 
    ## Run 155 stress 0.07411947 
    ## Run 156 stress 0.06477881 
    ## Run 157 stress 0.06477828 
    ## Run 158 stress 0.06477939 
    ## Run 159 stress 0.08226164 
    ## Run 160 stress 0.06623458 
    ## Run 161 stress 0.06623458 
    ## Run 162 stress 0.06623459 
    ## Run 163 stress 0.07411947 
    ## Run 164 stress 0.06623458 
    ## Run 165 stress 0.06477957 
    ## Run 166 stress 0.06477947 
    ## Run 167 stress 0.07411945 
    ## Run 168 stress 0.06124515 
    ## ... Procrustes: rmse 0.01433982  max resid 0.03977642 
    ## Run 169 stress 0.07166968 
    ## Run 170 stress 0.07166982 
    ## Run 171 stress 0.07411939 
    ## Run 172 stress 0.06477874 
    ## Run 173 stress 0.0741195 
    ## Run 174 stress 0.3148812 
    ## Run 175 stress 0.06623459 
    ## Run 176 stress 0.06124489 
    ## ... Procrustes: rmse 0.01492653  max resid 0.04135048 
    ## Run 177 stress 0.3407343 
    ## Run 178 stress 0.07166971 
    ## Run 179 stress 0.0741194 
    ## Run 180 stress 0.06478007 
    ## Run 181 stress 0.08226157 
    ## Run 182 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 3.422688e-06  max resid 6.813747e-06 
    ## ... Similar to previous best
    ## Run 183 stress 0.06477895 
    ## Run 184 stress 0.06623457 
    ## Run 185 stress 0.3266727 
    ## Run 186 stress 0.06123361 
    ## ... Procrustes: rmse 2.153033e-05  max resid 5.161639e-05 
    ## ... Similar to previous best
    ## Run 187 stress 0.0662346 
    ## Run 188 stress 0.07166974 
    ## Run 189 stress 0.06477865 
    ## Run 190 stress 0.0716697 
    ## Run 191 stress 0.07411942 
    ## Run 192 stress 0.09468991 
    ## Run 193 stress 0.07411943 
    ## Run 194 stress 0.06623463 
    ## Run 195 stress 0.07411946 
    ## Run 196 stress 0.06623465 
    ## Run 197 stress 0.06477898 
    ## Run 198 stress 0.07411938 
    ## Run 199 stress 0.07166979 
    ## Run 200 stress 0.07411941 
    ## Run 201 stress 0.0822616 
    ## Run 202 stress 0.07166964 
    ## Run 203 stress 0.07411939 
    ## Run 204 stress 0.07411939 
    ## Run 205 stress 0.06477931 
    ## Run 206 stress 0.06623461 
    ## Run 207 stress 0.06478008 
    ## Run 208 stress 0.06623457 
    ## Run 209 stress 0.0647788 
    ## Run 210 stress 0.06124493 
    ## ... Procrustes: rmse 0.01498614  max resid 0.04150934 
    ## Run 211 stress 0.06623466 
    ## Run 212 stress 0.08226158 
    ## Run 213 stress 0.0647794 
    ## Run 214 stress 0.06477852 
    ## Run 215 stress 0.07411939 
    ## Run 216 stress 0.07166964 
    ## Run 217 stress 0.06477842 
    ## Run 218 stress 0.06477972 
    ## Run 219 stress 0.07411945 
    ## Run 220 stress 0.06477868 
    ## Run 221 stress 0.06123363 
    ## ... Procrustes: rmse 7.569281e-05  max resid 0.000162429 
    ## ... Similar to previous best
    ## Run 222 stress 0.08226162 
    ## Run 223 stress 0.07411942 
    ## Run 224 stress 0.06477846 
    ## Run 225 stress 0.06623458 
    ## Run 226 stress 0.06477937 
    ## Run 227 stress 0.07166969 
    ## Run 228 stress 0.06477882 
    ## Run 229 stress 0.06623464 
    ## Run 230 stress 0.06477905 
    ## Run 231 stress 0.0662346 
    ## Run 232 stress 0.06477887 
    ## Run 233 stress 0.06623462 
    ## Run 234 stress 0.2362366 
    ## Run 235 stress 0.2591203 
    ## Run 236 stress 0.06477842 
    ## Run 237 stress 0.06477954 
    ## Run 238 stress 0.06124508 
    ## ... Procrustes: rmse 0.01515339  max resid 0.04195846 
    ## Run 239 stress 0.08226156 
    ## Run 240 stress 0.06623457 
    ## Run 241 stress 0.07411951 
    ## Run 242 stress 0.06124487 
    ## ... Procrustes: rmse 0.01488341  max resid 0.04123477 
    ## Run 243 stress 0.06124502 
    ## ... Procrustes: rmse 0.01510591  max resid 0.04183055 
    ## Run 244 stress 0.07411947 
    ## Run 245 stress 0.3407342 
    ## Run 246 stress 0.3112505 
    ## Run 247 stress 0.06623461 
    ## Run 248 stress 0.07411943 
    ## Run 249 stress 0.06477828 
    ## Run 250 stress 0.06124519 
    ## ... Procrustes: rmse 0.01523725  max resid 0.04218189 
    ## Run 251 stress 0.0612336 
    ## ... Procrustes: rmse 1.728416e-05  max resid 2.167457e-05 
    ## ... Similar to previous best
    ## Run 252 stress 0.06623457 
    ## Run 253 stress 0.06477887 
    ## Run 254 stress 0.06623459 
    ## Run 255 stress 0.06623459 
    ## Run 256 stress 0.0647793 
    ## Run 257 stress 0.08226177 
    ## Run 258 stress 0.06477901 
    ## Run 259 stress 0.06623457 
    ## Run 260 stress 0.06123361 
    ## ... Procrustes: rmse 3.125693e-05  max resid 6.574892e-05 
    ## ... Similar to previous best
    ## Run 261 stress 0.07166971 
    ## Run 262 stress 0.07411939 
    ## Run 263 stress 0.0662346 
    ## Run 264 stress 0.08226176 
    ## Run 265 stress 0.07411944 
    ## Run 266 stress 0.0741195 
    ## Run 267 stress 0.07411949 
    ## Run 268 stress 0.07166967 
    ## Run 269 stress 0.07411938 
    ## Run 270 stress 0.07166974 
    ## Run 271 stress 0.06124499 
    ## ... Procrustes: rmse 0.01507256  max resid 0.04174182 
    ## Run 272 stress 0.06124501 
    ## ... Procrustes: rmse 0.0150972  max resid 0.04180749 
    ## Run 273 stress 0.07411941 
    ## Run 274 stress 0.06623457 
    ## Run 275 stress 0.06477914 
    ## Run 276 stress 0.06623459 
    ## Run 277 stress 0.07411949 
    ## Run 278 stress 0.06623458 
    ## Run 279 stress 0.06124517 
    ## ... Procrustes: rmse 0.01522585  max resid 0.0421519 
    ## Run 280 stress 0.06477898 
    ## Run 281 stress 0.0741195 
    ## Run 282 stress 0.07411952 
    ## Run 283 stress 0.0662346 
    ## Run 284 stress 0.06623457 
    ## Run 285 stress 0.06623457 
    ## Run 286 stress 0.0647783 
    ## Run 287 stress 0.07166975 
    ## Run 288 stress 0.06623462 
    ## Run 289 stress 0.07411947 
    ## Run 290 stress 0.06124511 
    ## ... Procrustes: rmse 0.01517544  max resid 0.04201829 
    ## Run 291 stress 0.07166974 
    ## Run 292 stress 0.0741195 
    ## Run 293 stress 0.06124514 
    ## ... Procrustes: rmse 0.01434512  max resid 0.03979064 
    ## Run 294 stress 0.07411951 
    ## Run 295 stress 0.0741194 
    ## Run 296 stress 0.06124504 
    ## ... Procrustes: rmse 0.01512258  max resid 0.0418749 
    ## Run 297 stress 0.08226157 
    ## Run 298 stress 0.07411938 
    ## Run 299 stress 0.06623458 
    ## Run 300 stress 0.06124495 
    ## ... Procrustes: rmse 0.01502663  max resid 0.04161894 
    ## Run 301 stress 0.07166975 
    ## Run 302 stress 0.06623457 
    ## Run 303 stress 0.08226166 
    ## Run 304 stress 0.0662346 
    ## Run 305 stress 0.08226167 
    ## Run 306 stress 0.06477846 
    ## Run 307 stress 0.06623458 
    ## Run 308 stress 0.06623463 
    ## Run 309 stress 0.07166984 
    ## Run 310 stress 0.06124516 
    ## ... Procrustes: rmse 0.01513371  max resid 0.04190138 
    ## Run 311 stress 0.07166965 
    ## Run 312 stress 0.07411952 
    ## Run 313 stress 0.0662346 
    ## Run 314 stress 0.06623462 
    ## Run 315 stress 0.07411938 
    ## Run 316 stress 0.3390013 
    ## Run 317 stress 0.07411942 
    ## Run 318 stress 0.0662346 
    ## Run 319 stress 0.06123361 
    ## ... Procrustes: rmse 1.407946e-05  max resid 2.132109e-05 
    ## ... Similar to previous best
    ## Run 320 stress 0.0662346 
    ## Run 321 stress 0.3121118 
    ## Run 322 stress 0.2548998 
    ## Run 323 stress 0.06477848 
    ## Run 324 stress 0.06623458 
    ## Run 325 stress 0.06124521 
    ## ... Procrustes: rmse 0.0152174  max resid 0.04212625 
    ## Run 326 stress 0.06124491 
    ## ... Procrustes: rmse 0.01497487  max resid 0.04148015 
    ## Run 327 stress 0.2361811 
    ## Run 328 stress 0.0822616 
    ## Run 329 stress 0.07166995 
    ## Run 330 stress 0.07166966 
    ## Run 331 stress 0.07411938 
    ## Run 332 stress 0.08226166 
    ## Run 333 stress 0.07411942 
    ## Run 334 stress 0.06124505 
    ## ... Procrustes: rmse 0.01513539  max resid 0.04190976 
    ## Run 335 stress 0.08226164 
    ## Run 336 stress 0.07411938 
    ## Run 337 stress 0.06623457 
    ## Run 338 stress 0.06623457 
    ## Run 339 stress 0.09469047 
    ## Run 340 stress 0.0612453 
    ## ... Procrustes: rmse 0.01423519  max resid 0.03949505 
    ## Run 341 stress 0.08226156 
    ## Run 342 stress 0.06124494 
    ## ... Procrustes: rmse 0.01499761  max resid 0.04153986 
    ## Run 343 stress 0.07166973 
    ## Run 344 stress 0.2585395 
    ## Run 345 stress 0.07411939 
    ## Run 346 stress 0.06124515 
    ## ... Procrustes: rmse 0.01520654  max resid 0.04209941 
    ## Run 347 stress 0.07166965 
    ## Run 348 stress 0.07166964 
    ## Run 349 stress 0.07411945 
    ## Run 350 stress 0.08226161 
    ## Run 351 stress 0.06477975 
    ## Run 352 stress 0.07411938 
    ## Run 353 stress 0.06623462 
    ## Run 354 stress 0.06124513 
    ## ... Procrustes: rmse 0.01519744  max resid 0.04207554 
    ## Run 355 stress 0.07411939 
    ## Run 356 stress 0.07411952 
    ## Run 357 stress 0.06623458 
    ## Run 358 stress 0.07411947 
    ## Run 359 stress 0.09300801 
    ## Run 360 stress 0.06477887 
    ## Run 361 stress 0.07411948 
    ## Run 362 stress 0.07411948 
    ## Run 363 stress 0.06477924 
    ## Run 364 stress 0.07166965 
    ## Run 365 stress 0.0612336 
    ## ... Procrustes: rmse 5.717892e-06  max resid 8.908335e-06 
    ## ... Similar to previous best
    ## Run 366 stress 0.06477852 
    ## Run 367 stress 0.07166965 
    ## Run 368 stress 0.0716697 
    ## Run 369 stress 0.06124486 
    ## ... Procrustes: rmse 0.01483685  max resid 0.04111008 
    ## Run 370 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 1.552563e-05  max resid 2.014615e-05 
    ## ... Similar to previous best
    ## Run 371 stress 0.08226163 
    ## Run 372 stress 0.06124492 
    ## ... Procrustes: rmse 0.01476065  max resid 0.04090977 
    ## Run 373 stress 0.07411948 
    ## Run 374 stress 0.0647795 
    ## Run 375 stress 0.06623459 
    ## Run 376 stress 0.07411938 
    ## Run 377 stress 0.07411943 
    ## Run 378 stress 0.06477829 
    ## Run 379 stress 0.06623458 
    ## Run 380 stress 0.06123361 
    ## ... Procrustes: rmse 5.086799e-05  max resid 7.233145e-05 
    ## ... Similar to previous best
    ## Run 381 stress 0.06477839 
    ## Run 382 stress 0.061245 
    ## ... Procrustes: rmse 0.01508017  max resid 0.04176179 
    ## Run 383 stress 0.06124524 
    ## ... Procrustes: rmse 0.015256  max resid 0.04223056 
    ## Run 384 stress 0.0741194 
    ## Run 385 stress 0.0741195 
    ## Run 386 stress 0.06623463 
    ## Run 387 stress 0.312836 
    ## Run 388 stress 0.06623464 
    ## Run 389 stress 0.06623459 
    ## Run 390 stress 0.06477877 
    ## Run 391 stress 0.07411941 
    ## Run 392 stress 0.07411945 
    ## Run 393 stress 0.0741195 
    ## Run 394 stress 0.06123361 
    ## ... Procrustes: rmse 2.187255e-05  max resid 3.236361e-05 
    ## ... Similar to previous best
    ## Run 395 stress 0.06123362 
    ## ... Procrustes: rmse 4.532784e-05  max resid 6.325796e-05 
    ## ... Similar to previous best
    ## Run 396 stress 0.07411944 
    ## Run 397 stress 0.0647789 
    ## Run 398 stress 0.09468853 
    ## Run 399 stress 0.06477827 
    ## Run 400 stress 0.07411939 
    ## Run 401 stress 0.0647797 
    ## Run 402 stress 0.06477984 
    ## Run 403 stress 0.06623457 
    ## Run 404 stress 0.06623464 
    ## Run 405 stress 0.06623458 
    ## Run 406 stress 0.07166967 
    ## Run 407 stress 0.0612336 
    ## ... Procrustes: rmse 5.974506e-06  max resid 9.442119e-06 
    ## ... Similar to previous best
    ## Run 408 stress 0.0612336 
    ## ... Procrustes: rmse 1.870038e-05  max resid 2.37048e-05 
    ## ... Similar to previous best
    ## Run 409 stress 0.07166975 
    ## Run 410 stress 0.07166969 
    ## Run 411 stress 0.3035107 
    ## Run 412 stress 0.06124506 
    ## ... Procrustes: rmse 0.01513595  max resid 0.0419111 
    ## Run 413 stress 0.06623458 
    ## Run 414 stress 0.07166981 
    ## Run 415 stress 0.06124498 
    ## ... Procrustes: rmse 0.01506134  max resid 0.04171147 
    ## Run 416 stress 0.061245 
    ## ... Procrustes: rmse 0.01507303  max resid 0.04174333 
    ## Run 417 stress 0.06623457 
    ## Run 418 stress 0.06477993 
    ## Run 419 stress 0.06123361 
    ## ... Procrustes: rmse 4.220738e-05  max resid 9.509805e-05 
    ## ... Similar to previous best
    ## Run 420 stress 0.08226158 
    ## Run 421 stress 0.06124514 
    ## ... Procrustes: rmse 0.01519728  max resid 0.04207523 
    ## Run 422 stress 0.06623462 
    ## Run 423 stress 0.06124496 
    ## ... Procrustes: rmse 0.01502375  max resid 0.04161043 
    ## Run 424 stress 0.07166981 
    ## Run 425 stress 0.07166984 
    ## Run 426 stress 0.06623463 
    ## Run 427 stress 0.07166965 
    ## Run 428 stress 0.08226166 
    ## Run 429 stress 0.06623461 
    ## Run 430 stress 0.0741194 
    ## Run 431 stress 0.07411942 
    ## Run 432 stress 0.06124487 
    ## ... Procrustes: rmse 0.0148478  max resid 0.04113986 
    ## Run 433 stress 0.06623459 
    ## Run 434 stress 0.06124493 
    ## ... Procrustes: rmse 0.01499466  max resid 0.04153326 
    ## Run 435 stress 0.07166971 
    ## Run 436 stress 0.06623458 
    ## Run 437 stress 0.07166968 
    ## Run 438 stress 0.07411946 
    ## Run 439 stress 0.06477913 
    ## Run 440 stress 0.06477854 
    ## Run 441 stress 0.06124503 
    ## ... Procrustes: rmse 0.01510785  max resid 0.04183628 
    ## Run 442 stress 0.06123361 
    ## ... Procrustes: rmse 2.461862e-05  max resid 3.20406e-05 
    ## ... Similar to previous best
    ## Run 443 stress 0.09093869 
    ## Run 444 stress 0.08226159 
    ## Run 445 stress 0.07166968 
    ## Run 446 stress 0.07166968 
    ## Run 447 stress 0.06477977 
    ## Run 448 stress 0.07166973 
    ## Run 449 stress 0.2961559 
    ## Run 450 stress 0.06623459 
    ## Run 451 stress 0.07411942 
    ## Run 452 stress 0.0612336 
    ## ... Procrustes: rmse 1.003539e-05  max resid 1.235738e-05 
    ## ... Similar to previous best
    ## Run 453 stress 0.08226158 
    ## Run 454 stress 0.0662346 
    ## Run 455 stress 0.06124493 
    ## ... Procrustes: rmse 0.0149919  max resid 0.04152586 
    ## Run 456 stress 0.0662346 
    ## Run 457 stress 0.07411941 
    ## Run 458 stress 0.06124533 
    ## ... Procrustes: rmse 0.01421156  max resid 0.03943051 
    ## Run 459 stress 0.06477955 
    ## Run 460 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 9.009891e-06  max resid 1.816934e-05 
    ## ... Similar to previous best
    ## Run 461 stress 0.06477841 
    ## Run 462 stress 0.07411947 
    ## Run 463 stress 0.07166964 
    ## Run 464 stress 0.06124513 
    ## ... Procrustes: rmse 0.01519016  max resid 0.04205684 
    ## Run 465 stress 0.07411948 
    ## Run 466 stress 0.06124491 
    ## ... Procrustes: rmse 0.01493997  max resid 0.04138522 
    ## Run 467 stress 0.0612451 
    ## ... Procrustes: rmse 0.01516762  max resid 0.04199596 
    ## Run 468 stress 0.06124493 
    ## ... Procrustes: rmse 0.01499841  max resid 0.04154344 
    ## Run 469 stress 0.0662346 
    ## Run 470 stress 0.07411942 
    ## Run 471 stress 0.07411938 
    ## Run 472 stress 0.06623457 
    ## Run 473 stress 0.07166978 
    ## Run 474 stress 0.07411938 
    ## Run 475 stress 0.07166975 
    ## Run 476 stress 0.07411947 
    ## Run 477 stress 0.07411948 
    ## Run 478 stress 0.06623459 
    ## Run 479 stress 0.08226164 
    ## Run 480 stress 0.06124505 
    ## ... Procrustes: rmse 0.01512101  max resid 0.04187232 
    ## Run 481 stress 0.06123361 
    ## ... Procrustes: rmse 2.035017e-05  max resid 2.731013e-05 
    ## ... Similar to previous best
    ## Run 482 stress 0.07166986 
    ## Run 483 stress 0.0822617 
    ## Run 484 stress 0.06623463 
    ## Run 485 stress 0.0741194 
    ## Run 486 stress 0.07411955 
    ## Run 487 stress 0.0612336 
    ## ... Procrustes: rmse 5.842561e-06  max resid 1.299736e-05 
    ## ... Similar to previous best
    ## Run 488 stress 0.3032918 
    ## Run 489 stress 0.07411942 
    ## Run 490 stress 0.06477952 
    ## Run 491 stress 0.06623459 
    ## Run 492 stress 0.2592302 
    ## Run 493 stress 0.2345379 
    ## Run 494 stress 0.06124487 
    ## ... Procrustes: rmse 0.01487064  max resid 0.04120124 
    ## Run 495 stress 0.0741194 
    ## Run 496 stress 0.08226156 
    ## Run 497 stress 0.0741194 
    ## Run 498 stress 0.06477977 
    ## Run 499 stress 0.07411945 
    ## Run 500 stress 0.07166979 
    ## *** Best solution repeated 3 times

``` r
# Stratified lakes and ocean sites
SD_beta_env_SO_NMDS <- metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.59483e-05 
    ## Run 1 stress 8.39807e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001580213  max resid 0.0002448616 
    ## ... Similar to previous best
    ## Run 2 stress 9.72512e-05 
    ## ... Procrustes: rmse 0.0001563809  max resid 0.0003123254 
    ## ... Similar to previous best
    ## Run 3 stress 9.781484e-05 
    ## ... Procrustes: rmse 0.0001699595  max resid 0.0002267186 
    ## ... Similar to previous best
    ## Run 4 stress 9.508724e-05 
    ## ... Procrustes: rmse 0.0001256798  max resid 0.0002026536 
    ## ... Similar to previous best
    ## Run 5 stress 9.213618e-05 
    ## ... Procrustes: rmse 0.0001446337  max resid 0.000260323 
    ## ... Similar to previous best
    ## Run 6 stress 8.926567e-05 
    ## ... Procrustes: rmse 0.0001486353  max resid 0.0002557376 
    ## ... Similar to previous best
    ## Run 7 stress 9.926632e-05 
    ## ... Procrustes: rmse 9.717279e-05  max resid 0.0002324338 
    ## ... Similar to previous best
    ## Run 8 stress 6.697224e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000137534  max resid 0.0002659578 
    ## ... Similar to previous best
    ## Run 9 stress 9.446237e-05 
    ## ... Procrustes: rmse 0.0001763296  max resid 0.0003524381 
    ## ... Similar to previous best
    ## Run 10 stress 9.507435e-05 
    ## ... Procrustes: rmse 0.0001195811  max resid 0.0002648321 
    ## ... Similar to previous best
    ## Run 11 stress 9.468093e-05 
    ## ... Procrustes: rmse 0.0001756143  max resid 0.0003600844 
    ## ... Similar to previous best
    ## Run 12 stress 9.027327e-05 
    ## ... Procrustes: rmse 0.0001433288  max resid 0.0002705288 
    ## ... Similar to previous best
    ## Run 13 stress 8.291617e-05 
    ## ... Procrustes: rmse 0.0001223632  max resid 0.0001939931 
    ## ... Similar to previous best
    ## Run 14 stress 9.275766e-05 
    ## ... Procrustes: rmse 0.0001283874  max resid 0.0002629264 
    ## ... Similar to previous best
    ## Run 15 stress 9.264313e-05 
    ## ... Procrustes: rmse 0.0001433415  max resid 0.0002505287 
    ## ... Similar to previous best
    ## Run 16 stress 9.386158e-05 
    ## ... Procrustes: rmse 0.0001630505  max resid 0.000249093 
    ## ... Similar to previous best
    ## Run 17 stress 9.794724e-05 
    ## ... Procrustes: rmse 0.0001851683  max resid 0.0003442381 
    ## ... Similar to previous best
    ## Run 18 stress 9.437709e-05 
    ## ... Procrustes: rmse 0.0001547865  max resid 0.0002327061 
    ## ... Similar to previous best
    ## Run 19 stress 8.957711e-05 
    ## ... Procrustes: rmse 0.0001378305  max resid 0.000264335 
    ## ... Similar to previous best
    ## Run 20 stress 8.9764e-05 
    ## ... Procrustes: rmse 0.0001679142  max resid 0.0003521498 
    ## ... Similar to previous best
    ## Run 21 stress 9.464312e-05 
    ## ... Procrustes: rmse 0.0001246134  max resid 0.0002313659 
    ## ... Similar to previous best
    ## Run 22 stress 9.74816e-05 
    ## ... Procrustes: rmse 0.0001899819  max resid 0.0002857032 
    ## ... Similar to previous best
    ## Run 23 stress 9.954837e-05 
    ## ... Procrustes: rmse 0.0002005291  max resid 0.0003815849 
    ## ... Similar to previous best
    ## Run 24 stress 8.69582e-05 
    ## ... Procrustes: rmse 0.0001753874  max resid 0.0003471899 
    ## ... Similar to previous best
    ## Run 25 stress 9.235217e-05 
    ## ... Procrustes: rmse 0.0001192142  max resid 0.0002213947 
    ## ... Similar to previous best
    ## Run 26 stress 9.564972e-05 
    ## ... Procrustes: rmse 0.0001318383  max resid 0.0002109215 
    ## ... Similar to previous best
    ## Run 27 stress 9.36907e-05 
    ## ... Procrustes: rmse 0.000121143  max resid 0.0002322797 
    ## ... Similar to previous best
    ## Run 28 stress 9.956361e-05 
    ## ... Procrustes: rmse 0.0001138331  max resid 0.0001683886 
    ## ... Similar to previous best
    ## Run 29 stress 9.644151e-05 
    ## ... Procrustes: rmse 0.0001278516  max resid 0.0002454197 
    ## ... Similar to previous best
    ## Run 30 stress 9.918597e-05 
    ## ... Procrustes: rmse 0.0001596898  max resid 0.0002530417 
    ## ... Similar to previous best
    ## Run 31 stress 9.629769e-05 
    ## ... Procrustes: rmse 0.000110443  max resid 0.0002110523 
    ## ... Similar to previous best
    ## Run 32 stress 9.551382e-05 
    ## ... Procrustes: rmse 0.0001950433  max resid 0.0002882687 
    ## ... Similar to previous best
    ## Run 33 stress 9.820579e-05 
    ## ... Procrustes: rmse 0.0001659707  max resid 0.0003015795 
    ## ... Similar to previous best
    ## Run 34 stress 9.185498e-05 
    ## ... Procrustes: rmse 0.0001507821  max resid 0.0002641901 
    ## ... Similar to previous best
    ## Run 35 stress 9.669139e-05 
    ## ... Procrustes: rmse 8.98177e-05  max resid 0.0001262367 
    ## ... Similar to previous best
    ## Run 36 stress 9.356881e-05 
    ## ... Procrustes: rmse 0.0001206037  max resid 0.0002139417 
    ## ... Similar to previous best
    ## Run 37 stress 9.737003e-05 
    ## ... Procrustes: rmse 0.0001811169  max resid 0.0003796594 
    ## ... Similar to previous best
    ## Run 38 stress 9.745177e-05 
    ## ... Procrustes: rmse 0.0001502619  max resid 0.0002739916 
    ## ... Similar to previous best
    ## Run 39 stress 8.846252e-05 
    ## ... Procrustes: rmse 0.0001524636  max resid 0.0002342131 
    ## ... Similar to previous best
    ## Run 40 stress 8.78246e-05 
    ## ... Procrustes: rmse 0.0001521474  max resid 0.0003091756 
    ## ... Similar to previous best
    ## Run 41 stress 8.265501e-05 
    ## ... Procrustes: rmse 9.696672e-05  max resid 0.0002180257 
    ## ... Similar to previous best
    ## Run 42 stress 9.945541e-05 
    ## ... Procrustes: rmse 9.726668e-05  max resid 0.0001561614 
    ## ... Similar to previous best
    ## Run 43 stress 9.738277e-05 
    ## ... Procrustes: rmse 0.0001900003  max resid 0.0003501947 
    ## ... Similar to previous best
    ## Run 44 stress 9.771435e-05 
    ## ... Procrustes: rmse 8.65562e-05  max resid 0.0001143503 
    ## ... Similar to previous best
    ## Run 45 stress 9.35977e-05 
    ## ... Procrustes: rmse 9.839223e-05  max resid 0.0001826626 
    ## ... Similar to previous best
    ## Run 46 stress 9.866379e-05 
    ## ... Procrustes: rmse 0.0001592434  max resid 0.0003376617 
    ## ... Similar to previous best
    ## Run 47 stress 9.310223e-05 
    ## ... Procrustes: rmse 0.0001513128  max resid 0.000286139 
    ## ... Similar to previous best
    ## Run 48 stress 9.733648e-05 
    ## ... Procrustes: rmse 0.0001800925  max resid 0.0003640427 
    ## ... Similar to previous best
    ## Run 49 stress 0.2469684 
    ## Run 50 stress 9.937972e-05 
    ## ... Procrustes: rmse 0.0001545428  max resid 0.0002617149 
    ## ... Similar to previous best
    ## Run 51 stress 9.019405e-05 
    ## ... Procrustes: rmse 8.376953e-05  max resid 0.0001285772 
    ## ... Similar to previous best
    ## Run 52 stress 8.832531e-05 
    ## ... Procrustes: rmse 0.0001135886  max resid 0.0001812048 
    ## ... Similar to previous best
    ## Run 53 stress 9.989238e-05 
    ## ... Procrustes: rmse 0.0001139165  max resid 0.0001723389 
    ## ... Similar to previous best
    ## Run 54 stress 9.892463e-05 
    ## ... Procrustes: rmse 0.000199991  max resid 0.0002967283 
    ## ... Similar to previous best
    ## Run 55 stress 9.552697e-05 
    ## ... Procrustes: rmse 0.0001857993  max resid 0.0003451591 
    ## ... Similar to previous best
    ## Run 56 stress 8.92545e-05 
    ## ... Procrustes: rmse 0.0001628586  max resid 0.0003468539 
    ## ... Similar to previous best
    ## Run 57 stress 9.076984e-05 
    ## ... Procrustes: rmse 0.0001390618  max resid 0.0002387064 
    ## ... Similar to previous best
    ## Run 58 stress 9.356987e-05 
    ## ... Procrustes: rmse 0.0001596336  max resid 0.0003350049 
    ## ... Similar to previous best
    ## Run 59 stress 9.203608e-05 
    ## ... Procrustes: rmse 9.576757e-05  max resid 0.0001424027 
    ## ... Similar to previous best
    ## Run 60 stress 9.934826e-05 
    ## ... Procrustes: rmse 0.0001083844  max resid 0.0001587407 
    ## ... Similar to previous best
    ## Run 61 stress 9.506404e-05 
    ## ... Procrustes: rmse 0.0001937593  max resid 0.0002883958 
    ## ... Similar to previous best
    ## Run 62 stress 8.806341e-05 
    ## ... Procrustes: rmse 0.0001042321  max resid 0.0001749699 
    ## ... Similar to previous best
    ## Run 63 stress 8.884375e-05 
    ## ... Procrustes: rmse 0.0001139151  max resid 0.0002175949 
    ## ... Similar to previous best
    ## Run 64 stress 9.332823e-05 
    ## ... Procrustes: rmse 0.0001000777  max resid 0.0001684711 
    ## ... Similar to previous best
    ## Run 65 stress 9.258206e-05 
    ## ... Procrustes: rmse 5.879255e-05  max resid 9.943322e-05 
    ## ... Similar to previous best
    ## Run 66 stress 9.95249e-05 
    ## ... Procrustes: rmse 0.0001659577  max resid 0.0002930704 
    ## ... Similar to previous best
    ## Run 67 stress 9.012035e-05 
    ## ... Procrustes: rmse 0.0001388095  max resid 0.0002646826 
    ## ... Similar to previous best
    ## Run 68 stress 9.370628e-05 
    ## ... Procrustes: rmse 9.708032e-05  max resid 0.0001446411 
    ## ... Similar to previous best
    ## Run 69 stress 9.722077e-05 
    ## ... Procrustes: rmse 0.0002047855  max resid 0.0003003576 
    ## ... Similar to previous best
    ## Run 70 stress 9.756179e-05 
    ## ... Procrustes: rmse 0.0001526897  max resid 0.0002768995 
    ## ... Similar to previous best
    ## Run 71 stress 7.929623e-05 
    ## ... Procrustes: rmse 0.0001193477  max resid 0.0002274761 
    ## ... Similar to previous best
    ## Run 72 stress 9.657799e-05 
    ## ... Procrustes: rmse 0.0001676626  max resid 0.0002577333 
    ## ... Similar to previous best
    ## Run 73 stress 9.759854e-05 
    ## ... Procrustes: rmse 0.0001933823  max resid 0.0003666246 
    ## ... Similar to previous best
    ## Run 74 stress 9.806077e-05 
    ## ... Procrustes: rmse 0.0001731338  max resid 0.0002652071 
    ## ... Similar to previous best
    ## Run 75 stress 9.598594e-05 
    ## ... Procrustes: rmse 0.0001559689  max resid 0.0002765043 
    ## ... Similar to previous best
    ## Run 76 stress 9.627526e-05 
    ## ... Procrustes: rmse 0.0001250406  max resid 0.0002303888 
    ## ... Similar to previous best
    ## Run 77 stress 9.714087e-05 
    ## ... Procrustes: rmse 0.0001292864  max resid 0.0002361067 
    ## ... Similar to previous best
    ## Run 78 stress 9.024793e-05 
    ## ... Procrustes: rmse 0.0001499258  max resid 0.000322238 
    ## ... Similar to previous best
    ## Run 79 stress 9.331058e-05 
    ## ... Procrustes: rmse 0.0001750324  max resid 0.0003597352 
    ## ... Similar to previous best
    ## Run 80 stress 9.197527e-05 
    ## ... Procrustes: rmse 0.0001278123  max resid 0.000208449 
    ## ... Similar to previous best
    ## Run 81 stress 8.721763e-05 
    ## ... Procrustes: rmse 7.318905e-05  max resid 0.000103959 
    ## ... Similar to previous best
    ## Run 82 stress 0.2282323 
    ## Run 83 stress 8.917902e-05 
    ## ... Procrustes: rmse 9.841528e-05  max resid 0.0001826927 
    ## ... Similar to previous best
    ## Run 84 stress 9.423103e-05 
    ## ... Procrustes: rmse 0.0001701129  max resid 0.0003408802 
    ## ... Similar to previous best
    ## Run 85 stress 9.421792e-05 
    ## ... Procrustes: rmse 8.697799e-05  max resid 0.0001283344 
    ## ... Similar to previous best
    ## Run 86 stress 9.437455e-05 
    ## ... Procrustes: rmse 0.0001568881  max resid 0.0002658701 
    ## ... Similar to previous best
    ## Run 87 stress 9.915663e-05 
    ## ... Procrustes: rmse 0.000187218  max resid 0.0002826971 
    ## ... Similar to previous best
    ## Run 88 stress 9.58779e-05 
    ## ... Procrustes: rmse 0.0001730768  max resid 0.0003548954 
    ## ... Similar to previous best
    ## Run 89 stress 8.06687e-05 
    ## ... Procrustes: rmse 0.0001496857  max resid 0.0003367115 
    ## ... Similar to previous best
    ## Run 90 stress 9.610367e-05 
    ## ... Procrustes: rmse 0.0001920972  max resid 0.0002861584 
    ## ... Similar to previous best
    ## Run 91 stress 8.995841e-05 
    ## ... Procrustes: rmse 0.0001456686  max resid 0.000273896 
    ## ... Similar to previous best
    ## Run 92 stress 9.7934e-05 
    ## ... Procrustes: rmse 0.0002069716  max resid 0.0003032546 
    ## ... Similar to previous best
    ## Run 93 stress 0.2361264 
    ## Run 94 stress 9.947409e-05 
    ## ... Procrustes: rmse 0.0001874468  max resid 0.0003745548 
    ## ... Similar to previous best
    ## Run 95 stress 9.495284e-05 
    ## ... Procrustes: rmse 0.0001577756  max resid 0.0003324563 
    ## ... Similar to previous best
    ## Run 96 stress 9.688871e-05 
    ## ... Procrustes: rmse 0.0001066535  max resid 0.000156615 
    ## ... Similar to previous best
    ## Run 97 stress 9.971838e-05 
    ## ... Procrustes: rmse 0.0002048813  max resid 0.0003017397 
    ## ... Similar to previous best
    ## Run 98 stress 8.903161e-05 
    ## ... Procrustes: rmse 0.0001628781  max resid 0.0002604284 
    ## ... Similar to previous best
    ## Run 99 stress 9.799353e-05 
    ## ... Procrustes: rmse 0.0001453957  max resid 0.0002364119 
    ## ... Similar to previous best
    ## Run 100 stress 9.856996e-05 
    ## ... Procrustes: rmse 0.000111012  max resid 0.0001606263 
    ## ... Similar to previous best
    ## Run 101 stress 9.642057e-05 
    ## ... Procrustes: rmse 0.0001686967  max resid 0.0002594741 
    ## ... Similar to previous best
    ## Run 102 stress 9.518609e-05 
    ## ... Procrustes: rmse 0.0001698022  max resid 0.0002991899 
    ## ... Similar to previous best
    ## Run 103 stress 9.603687e-05 
    ## ... Procrustes: rmse 0.0001494147  max resid 0.0002536409 
    ## ... Similar to previous best
    ## Run 104 stress 7.701173e-05 
    ## ... Procrustes: rmse 0.0001537236  max resid 0.0002518367 
    ## ... Similar to previous best
    ## Run 105 stress 9.448394e-05 
    ## ... Procrustes: rmse 0.0001359007  max resid 0.0002536768 
    ## ... Similar to previous best
    ## Run 106 stress 9.788877e-05 
    ## ... Procrustes: rmse 0.0001903875  max resid 0.0003515612 
    ## ... Similar to previous best
    ## Run 107 stress 9.340275e-05 
    ## ... Procrustes: rmse 0.0001871443  max resid 0.0003636998 
    ## ... Similar to previous best
    ## Run 108 stress 9.046638e-05 
    ## ... Procrustes: rmse 0.0001725583  max resid 0.0002695951 
    ## ... Similar to previous best
    ## Run 109 stress 8.958292e-05 
    ## ... Procrustes: rmse 0.0001351658  max resid 0.0002607963 
    ## ... Similar to previous best
    ## Run 110 stress 9.447041e-05 
    ## ... Procrustes: rmse 0.0001003688  max resid 0.0001529441 
    ## ... Similar to previous best
    ## Run 111 stress 9.533203e-05 
    ## ... Procrustes: rmse 0.0001697751  max resid 0.0002979103 
    ## ... Similar to previous best
    ## Run 112 stress 9.175996e-05 
    ## ... Procrustes: rmse 0.0001668968  max resid 0.0003173041 
    ## ... Similar to previous best
    ## Run 113 stress 9.70541e-05 
    ## ... Procrustes: rmse 0.0001984601  max resid 0.0003773823 
    ## ... Similar to previous best
    ## Run 114 stress 9.540124e-05 
    ## ... Procrustes: rmse 0.0001905789  max resid 0.0003658283 
    ## ... Similar to previous best
    ## Run 115 stress 9.829286e-05 
    ## ... Procrustes: rmse 0.0002073978  max resid 0.0003036281 
    ## ... Similar to previous best
    ## Run 116 stress 9.59924e-05 
    ## ... Procrustes: rmse 0.0001900016  max resid 0.0002879175 
    ## ... Similar to previous best
    ## Run 117 stress 8.449311e-05 
    ## ... Procrustes: rmse 0.000173124  max resid 0.0003502796 
    ## ... Similar to previous best
    ## Run 118 stress 9.62549e-05 
    ## ... Procrustes: rmse 0.0001474954  max resid 0.0002625699 
    ## ... Similar to previous best
    ## Run 119 stress 8.726894e-05 
    ## ... Procrustes: rmse 8.647016e-05  max resid 0.0001304644 
    ## ... Similar to previous best
    ## Run 120 stress 9.584197e-05 
    ## ... Procrustes: rmse 0.0001651226  max resid 0.00034757 
    ## ... Similar to previous best
    ## Run 121 stress 9.166332e-05 
    ## ... Procrustes: rmse 8.772889e-05  max resid 0.0001658056 
    ## ... Similar to previous best
    ## Run 122 stress 9.911299e-05 
    ## ... Procrustes: rmse 0.0001807114  max resid 0.0003303114 
    ## ... Similar to previous best
    ## Run 123 stress 9.724823e-05 
    ## ... Procrustes: rmse 0.000179494  max resid 0.0003662952 
    ## ... Similar to previous best
    ## Run 124 stress 9.993618e-05 
    ## ... Procrustes: rmse 0.0001143217  max resid 0.0001689106 
    ## ... Similar to previous best
    ## Run 125 stress 9.034827e-05 
    ## ... Procrustes: rmse 8.657111e-05  max resid 0.0001337225 
    ## ... Similar to previous best
    ## Run 126 stress 8.537323e-05 
    ## ... Procrustes: rmse 0.0001137637  max resid 0.0002167837 
    ## ... Similar to previous best
    ## Run 127 stress 9.492053e-05 
    ## ... Procrustes: rmse 0.0001259466  max resid 0.0002330641 
    ## ... Similar to previous best
    ## Run 128 stress 9.38978e-05 
    ## ... Procrustes: rmse 0.0001456502  max resid 0.0002694563 
    ## ... Similar to previous best
    ## Run 129 stress 9.595778e-05 
    ## ... Procrustes: rmse 0.0001797002  max resid 0.0003638207 
    ## ... Similar to previous best
    ## Run 130 stress 9.27166e-05 
    ## ... Procrustes: rmse 0.0001323506  max resid 0.0002487655 
    ## ... Similar to previous best
    ## Run 131 stress 9.281208e-05 
    ## ... Procrustes: rmse 0.0001267303  max resid 0.0002201397 
    ## ... Similar to previous best
    ## Run 132 stress 9.626027e-05 
    ## ... Procrustes: rmse 0.0001603824  max resid 0.0003485893 
    ## ... Similar to previous best
    ## Run 133 stress 9.256073e-05 
    ## ... Procrustes: rmse 0.0001432356  max resid 0.0002679868 
    ## ... Similar to previous best
    ## Run 134 stress 9.274378e-05 
    ## ... Procrustes: rmse 0.0001024556  max resid 0.0002024305 
    ## ... Similar to previous best
    ## Run 135 stress 9.254565e-05 
    ## ... Procrustes: rmse 0.0001844343  max resid 0.0002820793 
    ## ... Similar to previous best
    ## Run 136 stress 9.772411e-05 
    ## ... Procrustes: rmse 0.0001090496  max resid 0.0001656527 
    ## ... Similar to previous best
    ## Run 137 stress 9.939639e-05 
    ## ... Procrustes: rmse 0.0001149252  max resid 0.0002147169 
    ## ... Similar to previous best
    ## Run 138 stress 9.320692e-05 
    ## ... Procrustes: rmse 0.000156126  max resid 0.0002396654 
    ## ... Similar to previous best
    ## Run 139 stress 9.537613e-05 
    ## ... Procrustes: rmse 0.0001048864  max resid 0.0001626713 
    ## ... Similar to previous best
    ## Run 140 stress 9.407552e-05 
    ## ... Procrustes: rmse 0.0001822191  max resid 0.0003489249 
    ## ... Similar to previous best
    ## Run 141 stress 9.802756e-05 
    ## ... Procrustes: rmse 0.0001834379  max resid 0.0003831972 
    ## ... Similar to previous best
    ## Run 142 stress 9.674933e-05 
    ## ... Procrustes: rmse 0.0001500793  max resid 0.0002813525 
    ## ... Similar to previous best
    ## Run 143 stress 9.790585e-05 
    ## ... Procrustes: rmse 0.0001569883  max resid 0.0002789069 
    ## ... Similar to previous best
    ## Run 144 stress 9.604368e-05 
    ## ... Procrustes: rmse 0.0001216408  max resid 0.0002301783 
    ## ... Similar to previous best
    ## Run 145 stress 7.430988e-05 
    ## ... Procrustes: rmse 0.0001442072  max resid 0.0002527234 
    ## ... Similar to previous best
    ## Run 146 stress 9.819566e-05 
    ## ... Procrustes: rmse 0.0001766171  max resid 0.0002543284 
    ## ... Similar to previous best
    ## Run 147 stress 9.792727e-05 
    ## ... Procrustes: rmse 0.0001932132  max resid 0.0003306684 
    ## ... Similar to previous best
    ## Run 148 stress 9.621267e-05 
    ## ... Procrustes: rmse 0.0001339817  max resid 0.0002202002 
    ## ... Similar to previous best
    ## Run 149 stress 9.62798e-05 
    ## ... Procrustes: rmse 0.0001032995  max resid 0.0001537536 
    ## ... Similar to previous best
    ## Run 150 stress 9.92205e-05 
    ## ... Procrustes: rmse 0.0001821562  max resid 0.0003338985 
    ## ... Similar to previous best
    ## Run 151 stress 8.818062e-05 
    ## ... Procrustes: rmse 0.0001455754  max resid 0.000225666 
    ## ... Similar to previous best
    ## Run 152 stress 8.942802e-05 
    ## ... Procrustes: rmse 0.0001203835  max resid 0.000221337 
    ## ... Similar to previous best
    ## Run 153 stress 9.628572e-05 
    ## ... Procrustes: rmse 0.0001064164  max resid 0.0001613345 
    ## ... Similar to previous best
    ## Run 154 stress 9.324527e-05 
    ## ... Procrustes: rmse 0.0001670061  max resid 0.0002982037 
    ## ... Similar to previous best
    ## Run 155 stress 9.935794e-05 
    ## ... Procrustes: rmse 0.000197327  max resid 0.000292648 
    ## ... Similar to previous best
    ## Run 156 stress 0.2721502 
    ## Run 157 stress 9.577105e-05 
    ## ... Procrustes: rmse 0.0001972486  max resid 0.0003746231 
    ## ... Similar to previous best
    ## Run 158 stress 9.389887e-05 
    ## ... Procrustes: rmse 0.0001879593  max resid 0.0003653189 
    ## ... Similar to previous best
    ## Run 159 stress 9.876102e-05 
    ## ... Procrustes: rmse 0.000176315  max resid 0.0003058032 
    ## ... Similar to previous best
    ## Run 160 stress 9.166764e-05 
    ## ... Procrustes: rmse 0.0001466562  max resid 0.0002938444 
    ## ... Similar to previous best
    ## Run 161 stress 9.796953e-05 
    ## ... Procrustes: rmse 0.0001972658  max resid 0.0002936615 
    ## ... Similar to previous best
    ## Run 162 stress 9.726066e-05 
    ## ... Procrustes: rmse 0.0001610984  max resid 0.0002495547 
    ## ... Similar to previous best
    ## Run 163 stress 9.419613e-05 
    ## ... Procrustes: rmse 7.523786e-05  max resid 0.0001222696 
    ## ... Similar to previous best
    ## Run 164 stress 9.583114e-05 
    ## ... Procrustes: rmse 0.0001496705  max resid 0.0002530115 
    ## ... Similar to previous best
    ## Run 165 stress 9.861665e-05 
    ## ... Procrustes: rmse 0.000208711  max resid 0.0003049449 
    ## ... Similar to previous best
    ## Run 166 stress 9.665747e-05 
    ## ... Procrustes: rmse 0.000194986  max resid 0.0002922828 
    ## ... Similar to previous best
    ## Run 167 stress 9.746522e-05 
    ## ... Procrustes: rmse 0.000197683  max resid 0.0002948347 
    ## ... Similar to previous best
    ## Run 168 stress 8.400108e-05 
    ## ... Procrustes: rmse 0.0001050915  max resid 0.0002228898 
    ## ... Similar to previous best
    ## Run 169 stress 9.512245e-05 
    ## ... Procrustes: rmse 0.0001409362  max resid 0.0002645785 
    ## ... Similar to previous best
    ## Run 170 stress 9.517101e-05 
    ## ... Procrustes: rmse 0.0001926082  max resid 0.0002904508 
    ## ... Similar to previous best
    ## Run 171 stress 9.900612e-05 
    ## ... Procrustes: rmse 0.0001251131  max resid 0.000240072 
    ## ... Similar to previous best
    ## Run 172 stress 9.754255e-05 
    ## ... Procrustes: rmse 0.0001994042  max resid 0.0003773808 
    ## ... Similar to previous best
    ## Run 173 stress 9.621092e-05 
    ## ... Procrustes: rmse 0.0001588716  max resid 0.000281236 
    ## ... Similar to previous best
    ## Run 174 stress 9.719438e-05 
    ## ... Procrustes: rmse 0.0001717858  max resid 0.0002625909 
    ## ... Similar to previous best
    ## Run 175 stress 9.768396e-05 
    ## ... Procrustes: rmse 0.0001093169  max resid 0.0001645271 
    ## ... Similar to previous best
    ## Run 176 stress 9.429186e-05 
    ## ... Procrustes: rmse 0.0002000556  max resid 0.0002958728 
    ## ... Similar to previous best
    ## Run 177 stress 9.676127e-05 
    ## ... Procrustes: rmse 0.000188636  max resid 0.0003379078 
    ## ... Similar to previous best
    ## Run 178 stress 9.907428e-05 
    ## ... Procrustes: rmse 0.0001997213  max resid 0.0003843261 
    ## ... Similar to previous best
    ## Run 179 stress 9.969584e-05 
    ## ... Procrustes: rmse 0.0001946151  max resid 0.0002923673 
    ## ... Similar to previous best
    ## Run 180 stress 7.751714e-05 
    ## ... Procrustes: rmse 0.0001410174  max resid 0.0002651078 
    ## ... Similar to previous best
    ## Run 181 stress 9.772629e-05 
    ## ... Procrustes: rmse 0.0001567278  max resid 0.0002797227 
    ## ... Similar to previous best
    ## Run 182 stress 9.877051e-05 
    ## ... Procrustes: rmse 0.0001763287  max resid 0.0003308988 
    ## ... Similar to previous best
    ## Run 183 stress 8.871066e-05 
    ## ... Procrustes: rmse 0.0001086852  max resid 0.0002096292 
    ## ... Similar to previous best
    ## Run 184 stress 9.976542e-05 
    ## ... Procrustes: rmse 0.0001876546  max resid 0.00031797 
    ## ... Similar to previous best
    ## Run 185 stress 9.250353e-05 
    ## ... Procrustes: rmse 0.0001508112  max resid 0.000328463 
    ## ... Similar to previous best
    ## Run 186 stress 9.572965e-05 
    ## ... Procrustes: rmse 0.0001953447  max resid 0.0002899135 
    ## ... Similar to previous best
    ## Run 187 stress 9.709723e-05 
    ## ... Procrustes: rmse 0.000159524  max resid 0.0003006091 
    ## ... Similar to previous best
    ## Run 188 stress 9.322763e-05 
    ## ... Procrustes: rmse 0.0001711673  max resid 0.0003579056 
    ## ... Similar to previous best
    ## Run 189 stress 9.912407e-05 
    ## ... Procrustes: rmse 0.0001600048  max resid 0.0002488083 
    ## ... Similar to previous best
    ## Run 190 stress 9.717248e-05 
    ## ... Procrustes: rmse 0.0001982623  max resid 0.0002920875 
    ## ... Similar to previous best
    ## Run 191 stress 9.100679e-05 
    ## ... Procrustes: rmse 0.0001191672  max resid 0.0002434541 
    ## ... Similar to previous best
    ## Run 192 stress 9.557555e-05 
    ## ... Procrustes: rmse 0.0001774952  max resid 0.0003605665 
    ## ... Similar to previous best
    ## Run 193 stress 8.764298e-05 
    ## ... Procrustes: rmse 0.0001284401  max resid 0.0002016868 
    ## ... Similar to previous best
    ## Run 194 stress 9.715428e-05 
    ## ... Procrustes: rmse 0.0001825405  max resid 0.0003683525 
    ## ... Similar to previous best
    ## Run 195 stress 8.715919e-05 
    ## ... Procrustes: rmse 0.0001599549  max resid 0.000343364 
    ## ... Similar to previous best
    ## Run 196 stress 9.899537e-05 
    ## ... Procrustes: rmse 7.465416e-05  max resid 0.0001155665 
    ## ... Similar to previous best
    ## Run 197 stress 9.645284e-05 
    ## ... Procrustes: rmse 0.0001808982  max resid 0.0003661741 
    ## ... Similar to previous best
    ## Run 198 stress 9.597639e-05 
    ## ... Procrustes: rmse 0.0001468273  max resid 0.0002239739 
    ## ... Similar to previous best
    ## Run 199 stress 9.843717e-05 
    ## ... Procrustes: rmse 0.0001621288  max resid 0.0002861162 
    ## ... Similar to previous best
    ## Run 200 stress 9.62191e-05 
    ## ... Procrustes: rmse 9.966264e-05  max resid 0.0001503735 
    ## ... Similar to previous best
    ## Run 201 stress 9.03514e-05 
    ## ... Procrustes: rmse 0.0001655033  max resid 0.0002771638 
    ## ... Similar to previous best
    ## Run 202 stress 9.510953e-05 
    ## ... Procrustes: rmse 0.0001480547  max resid 0.0002733342 
    ## ... Similar to previous best
    ## Run 203 stress 9.402173e-05 
    ## ... Procrustes: rmse 0.0001569152  max resid 0.0003029811 
    ## ... Similar to previous best
    ## Run 204 stress 9.935269e-05 
    ## ... Procrustes: rmse 0.0001807433  max resid 0.0003783992 
    ## ... Similar to previous best
    ## Run 205 stress 9.903294e-05 
    ## ... Procrustes: rmse 0.0001511963  max resid 0.0002759841 
    ## ... Similar to previous best
    ## Run 206 stress 9.766566e-05 
    ## ... Procrustes: rmse 0.0001268564  max resid 0.0002346642 
    ## ... Similar to previous best
    ## Run 207 stress 8.732064e-05 
    ## ... Procrustes: rmse 0.0001358328  max resid 0.0002812115 
    ## ... Similar to previous best
    ## Run 208 stress 9.536143e-05 
    ## ... Procrustes: rmse 0.0001655653  max resid 0.000293801 
    ## ... Similar to previous best
    ## Run 209 stress 7.871629e-05 
    ## ... Procrustes: rmse 9.33724e-05  max resid 0.0001903243 
    ## ... Similar to previous best
    ## Run 210 stress 9.92406e-05 
    ## ... Procrustes: rmse 0.0001272414  max resid 0.0001861962 
    ## ... Similar to previous best
    ## Run 211 stress 9.760879e-05 
    ## ... Procrustes: rmse 0.0001931608  max resid 0.0002876759 
    ## ... Similar to previous best
    ## Run 212 stress 9.614023e-05 
    ## ... Procrustes: rmse 0.0002036494  max resid 0.000300133 
    ## ... Similar to previous best
    ## Run 213 stress 9.78823e-05 
    ## ... Procrustes: rmse 0.0001831164  max resid 0.0003828296 
    ## ... Similar to previous best
    ## Run 214 stress 9.760887e-05 
    ## ... Procrustes: rmse 0.0001077767  max resid 0.0001653811 
    ## ... Similar to previous best
    ## Run 215 stress 9.622999e-05 
    ## ... Procrustes: rmse 0.0001804238  max resid 0.0003651975 
    ## ... Similar to previous best
    ## Run 216 stress 9.791983e-05 
    ## ... Procrustes: rmse 0.0001939518  max resid 0.0002888784 
    ## ... Similar to previous best
    ## Run 217 stress 9.148048e-05 
    ## ... Procrustes: rmse 0.0001635073  max resid 0.0003039639 
    ## ... Similar to previous best
    ## Run 218 stress 9.314582e-05 
    ## ... Procrustes: rmse 9.919049e-05  max resid 0.0001475318 
    ## ... Similar to previous best
    ## Run 219 stress 9.412329e-05 
    ## ... Procrustes: rmse 0.0001940166  max resid 0.0003715756 
    ## ... Similar to previous best
    ## Run 220 stress 9.59141e-05 
    ## ... Procrustes: rmse 0.0001779806  max resid 0.0003618759 
    ## ... Similar to previous best
    ## Run 221 stress 9.115183e-05 
    ## ... Procrustes: rmse 0.0001195701  max resid 0.0002259416 
    ## ... Similar to previous best
    ## Run 222 stress 9.975319e-05 
    ## ... Procrustes: rmse 0.0001165744  max resid 0.0001876476 
    ## ... Similar to previous best
    ## Run 223 stress 9.504633e-05 
    ## ... Procrustes: rmse 0.0001719318  max resid 0.0003443581 
    ## ... Similar to previous best
    ## Run 224 stress 9.421955e-05 
    ## ... Procrustes: rmse 0.0001835013  max resid 0.0003445485 
    ## ... Similar to previous best
    ## Run 225 stress 9.157692e-05 
    ## ... Procrustes: rmse 9.008581e-05  max resid 0.0001424991 
    ## ... Similar to previous best
    ## Run 226 stress 9.939708e-05 
    ## ... Procrustes: rmse 0.0001061128  max resid 0.0001638741 
    ## ... Similar to previous best
    ## Run 227 stress 8.795376e-05 
    ## ... Procrustes: rmse 0.0001599068  max resid 0.0003406024 
    ## ... Similar to previous best
    ## Run 228 stress 9.847306e-05 
    ## ... Procrustes: rmse 0.0001945765  max resid 0.0002922805 
    ## ... Similar to previous best
    ## Run 229 stress 9.920328e-05 
    ## ... Procrustes: rmse 0.0001340049  max resid 0.0002431947 
    ## ... Similar to previous best
    ## Run 230 stress 9.57182e-05 
    ## ... Procrustes: rmse 0.0001789176  max resid 0.0003646895 
    ## ... Similar to previous best
    ## Run 231 stress 8.168471e-05 
    ## ... Procrustes: rmse 7.39778e-05  max resid 0.0001159107 
    ## ... Similar to previous best
    ## Run 232 stress 9.625077e-05 
    ## ... Procrustes: rmse 0.0001676928  max resid 0.0002553527 
    ## ... Similar to previous best
    ## Run 233 stress 9.074238e-05 
    ## ... Procrustes: rmse 0.0001722096  max resid 0.0003249323 
    ## ... Similar to previous best
    ## Run 234 stress 9.943684e-05 
    ## ... Procrustes: rmse 0.0001615765  max resid 0.0002940496 
    ## ... Similar to previous best
    ## Run 235 stress 9.463714e-05 
    ## ... Procrustes: rmse 0.0001920577  max resid 0.0002892544 
    ## ... Similar to previous best
    ## Run 236 stress 8.582873e-05 
    ## ... Procrustes: rmse 0.0001534223  max resid 0.0003085898 
    ## ... Similar to previous best
    ## Run 237 stress 9.76754e-05 
    ## ... Procrustes: rmse 0.0001842048  max resid 0.0003712447 
    ## ... Similar to previous best
    ## Run 238 stress 9.083194e-05 
    ## ... Procrustes: rmse 0.0001158129  max resid 0.000221224 
    ## ... Similar to previous best
    ## Run 239 stress 9.872522e-05 
    ## ... Procrustes: rmse 0.0001845046  max resid 0.000344862 
    ## ... Similar to previous best
    ## Run 240 stress 9.79882e-05 
    ## ... Procrustes: rmse 0.0001724482  max resid 0.000346914 
    ## ... Similar to previous best
    ## Run 241 stress 9.056656e-05 
    ## ... Procrustes: rmse 0.0001385704  max resid 0.0002610027 
    ## ... Similar to previous best
    ## Run 242 stress 9.684546e-05 
    ## ... Procrustes: rmse 0.0001561889  max resid 0.0002805745 
    ## ... Similar to previous best
    ## Run 243 stress 9.963063e-05 
    ## ... Procrustes: rmse 0.0001077665  max resid 0.0001673431 
    ## ... Similar to previous best
    ## Run 244 stress 9.663999e-05 
    ## ... Procrustes: rmse 0.0001173973  max resid 0.0002325147 
    ## ... Similar to previous best
    ## Run 245 stress 9.288136e-05 
    ## ... Procrustes: rmse 0.0001734743  max resid 0.0002915135 
    ## ... Similar to previous best
    ## Run 246 stress 8.443009e-05 
    ## ... Procrustes: rmse 0.0001738965  max resid 0.0003481653 
    ## ... Similar to previous best
    ## Run 247 stress 9.418921e-05 
    ## ... Procrustes: rmse 0.0001023081  max resid 0.0001513301 
    ## ... Similar to previous best
    ## Run 248 stress 9.159825e-05 
    ## ... Procrustes: rmse 0.0001621621  max resid 0.0003437653 
    ## ... Similar to previous best
    ## Run 249 stress 9.555002e-05 
    ## ... Procrustes: rmse 0.0002022169  max resid 0.0002982519 
    ## ... Similar to previous best
    ## Run 250 stress 8.853265e-05 
    ## ... Procrustes: rmse 8.481924e-05  max resid 0.0001228043 
    ## ... Similar to previous best
    ## Run 251 stress 9.836988e-05 
    ## ... Procrustes: rmse 0.0001717557  max resid 0.0002792217 
    ## ... Similar to previous best
    ## Run 252 stress 9.756607e-05 
    ## ... Procrustes: rmse 0.0001784214  max resid 0.0003260201 
    ## ... Similar to previous best
    ## Run 253 stress 9.916073e-05 
    ## ... Procrustes: rmse 0.0001539871  max resid 0.0002452172 
    ## ... Similar to previous best
    ## Run 254 stress 9.6267e-05 
    ## ... Procrustes: rmse 0.0001679238  max resid 0.0002762282 
    ## ... Similar to previous best
    ## Run 255 stress 9.598118e-05 
    ## ... Procrustes: rmse 0.0001481718  max resid 0.0002852357 
    ## ... Similar to previous best
    ## Run 256 stress 9.697264e-05 
    ## ... Procrustes: rmse 0.000173339  max resid 0.0003587862 
    ## ... Similar to previous best
    ## Run 257 stress 9.798407e-05 
    ## ... Procrustes: rmse 0.0001792452  max resid 0.0003266541 
    ## ... Similar to previous best
    ## Run 258 stress 9.773726e-05 
    ## ... Procrustes: rmse 0.0001559779  max resid 0.0002693041 
    ## ... Similar to previous best
    ## Run 259 stress 9.261756e-05 
    ## ... Procrustes: rmse 0.0001934824  max resid 0.0003721361 
    ## ... Similar to previous best
    ## Run 260 stress 9.367158e-05 
    ## ... Procrustes: rmse 0.0001531352  max resid 0.0002711539 
    ## ... Similar to previous best
    ## Run 261 stress 9.468058e-05 
    ## ... Procrustes: rmse 0.0001933971  max resid 0.0002875066 
    ## ... Similar to previous best
    ## Run 262 stress 9.723244e-05 
    ## ... Procrustes: rmse 0.0001521305  max resid 0.000276847 
    ## ... Similar to previous best
    ## Run 263 stress 8.113194e-05 
    ## ... Procrustes: rmse 0.0001499675  max resid 0.0003318457 
    ## ... Similar to previous best
    ## Run 264 stress 9.905315e-05 
    ## ... Procrustes: rmse 0.0001729616  max resid 0.0003144693 
    ## ... Similar to previous best
    ## Run 265 stress 9.963333e-05 
    ## ... Procrustes: rmse 0.0002070367  max resid 0.0003025458 
    ## ... Similar to previous best
    ## Run 266 stress 9.161191e-05 
    ## ... Procrustes: rmse 8.409945e-05  max resid 0.0001124441 
    ## ... Similar to previous best
    ## Run 267 stress 9.135648e-05 
    ## ... Procrustes: rmse 0.0001336082  max resid 0.0002588675 
    ## ... Similar to previous best
    ## Run 268 stress 9.552655e-05 
    ## ... Procrustes: rmse 0.000157261  max resid 0.0002780681 
    ## ... Similar to previous best
    ## Run 269 stress 9.512649e-05 
    ## ... Procrustes: rmse 0.0001693878  max resid 0.0003131743 
    ## ... Similar to previous best
    ## Run 270 stress 9.942781e-05 
    ## ... Procrustes: rmse 0.0001491394  max resid 0.0002731444 
    ## ... Similar to previous best
    ## Run 271 stress 9.479799e-05 
    ## ... Procrustes: rmse 0.0001125515  max resid 0.000228691 
    ## ... Similar to previous best
    ## Run 272 stress 7.842852e-05 
    ## ... Procrustes: rmse 0.0001010935  max resid 0.0002206534 
    ## ... Similar to previous best
    ## Run 273 stress 0.2413829 
    ## Run 274 stress 9.561056e-05 
    ## ... Procrustes: rmse 0.0001744888  max resid 0.0003459756 
    ## ... Similar to previous best
    ## Run 275 stress 9.795214e-05 
    ## ... Procrustes: rmse 0.0001109641  max resid 0.0001666553 
    ## ... Similar to previous best
    ## Run 276 stress 9.544253e-05 
    ## ... Procrustes: rmse 0.0001607398  max resid 0.0002469604 
    ## ... Similar to previous best
    ## Run 277 stress 0.206342 
    ## Run 278 stress 8.980303e-05 
    ## ... Procrustes: rmse 0.000114415  max resid 0.0002199771 
    ## ... Similar to previous best
    ## Run 279 stress 9.139779e-05 
    ## ... Procrustes: rmse 0.0001595142  max resid 0.0002673327 
    ## ... Similar to previous best
    ## Run 280 stress 7.402704e-05 
    ## ... Procrustes: rmse 0.0001214752  max resid 0.000294902 
    ## ... Similar to previous best
    ## Run 281 stress 9.780586e-05 
    ## ... Procrustes: rmse 0.0001831309  max resid 0.0002769147 
    ## ... Similar to previous best
    ## Run 282 stress 9.589689e-05 
    ## ... Procrustes: rmse 0.0001791525  max resid 0.0003624936 
    ## ... Similar to previous best
    ## Run 283 stress 7.335807e-05 
    ## ... Procrustes: rmse 0.0001315697  max resid 0.0002946782 
    ## ... Similar to previous best
    ## Run 284 stress 9.679798e-05 
    ## ... Procrustes: rmse 0.0001336812  max resid 0.0002822884 
    ## ... Similar to previous best
    ## Run 285 stress 9.890402e-05 
    ## ... Procrustes: rmse 0.0001650271  max resid 0.0003349662 
    ## ... Similar to previous best
    ## Run 286 stress 9.42941e-05 
    ## ... Procrustes: rmse 0.0001732172  max resid 0.0003173087 
    ## ... Similar to previous best
    ## Run 287 stress 8.755753e-05 
    ## ... Procrustes: rmse 0.000162243  max resid 0.0003471288 
    ## ... Similar to previous best
    ## Run 288 stress 9.811436e-05 
    ## ... Procrustes: rmse 0.0001628032  max resid 0.0002475008 
    ## ... Similar to previous best
    ## Run 289 stress 9.787747e-05 
    ## ... Procrustes: rmse 0.0001869227  max resid 0.0003418668 
    ## ... Similar to previous best
    ## Run 290 stress 9.732332e-05 
    ## ... Procrustes: rmse 9.507968e-05  max resid 0.0001282913 
    ## ... Similar to previous best
    ## Run 291 stress 9.705787e-05 
    ## ... Procrustes: rmse 0.000129265  max resid 0.0001939822 
    ## ... Similar to previous best
    ## Run 292 stress 9.61994e-05 
    ## ... Procrustes: rmse 0.0001946649  max resid 0.0002920141 
    ## ... Similar to previous best
    ## Run 293 stress 9.912029e-05 
    ## ... Procrustes: rmse 0.000184189  max resid 0.0003719155 
    ## ... Similar to previous best
    ## Run 294 stress 9.888584e-05 
    ## ... Procrustes: rmse 0.0001652987  max resid 0.0002941506 
    ## ... Similar to previous best
    ## Run 295 stress 9.887008e-05 
    ## ... Procrustes: rmse 0.0001832774  max resid 0.0003815916 
    ## ... Similar to previous best
    ## Run 296 stress 8.780853e-05 
    ## ... Procrustes: rmse 0.0001253679  max resid 0.00027724 
    ## ... Similar to previous best
    ## Run 297 stress 9.741326e-05 
    ## ... Procrustes: rmse 0.0001658807  max resid 0.0002972633 
    ## ... Similar to previous best
    ## Run 298 stress 9.637035e-05 
    ## ... Procrustes: rmse 0.0001736848  max resid 0.0003022285 
    ## ... Similar to previous best
    ## Run 299 stress 9.425441e-05 
    ## ... Procrustes: rmse 0.0001340488  max resid 0.0002028751 
    ## ... Similar to previous best
    ## Run 300 stress 8.415457e-05 
    ## ... Procrustes: rmse 0.0001224356  max resid 0.0002479397 
    ## ... Similar to previous best
    ## Run 301 stress 9.386893e-05 
    ## ... Procrustes: rmse 0.0001840083  max resid 0.000316737 
    ## ... Similar to previous best
    ## Run 302 stress 8.698672e-05 
    ## ... Procrustes: rmse 0.000110376  max resid 0.0002125816 
    ## ... Similar to previous best
    ## Run 303 stress 9.534016e-05 
    ## ... Procrustes: rmse 0.0001486998  max resid 0.0002718338 
    ## ... Similar to previous best
    ## Run 304 stress 9.527921e-05 
    ## ... Procrustes: rmse 0.0001034688  max resid 0.0001568878 
    ## ... Similar to previous best
    ## Run 305 stress 9.330105e-05 
    ## ... Procrustes: rmse 0.0001399483  max resid 0.0002589703 
    ## ... Similar to previous best
    ## Run 306 stress 9.728845e-05 
    ## ... Procrustes: rmse 0.0001505953  max resid 0.0002636392 
    ## ... Similar to previous best
    ## Run 307 stress 9.811467e-05 
    ## ... Procrustes: rmse 0.0001399637  max resid 0.0002945283 
    ## ... Similar to previous best
    ## Run 308 stress 9.418474e-05 
    ## ... Procrustes: rmse 0.000188312  max resid 0.0002854199 
    ## ... Similar to previous best
    ## Run 309 stress 9.500666e-05 
    ## ... Procrustes: rmse 0.0001904721  max resid 0.0002873827 
    ## ... Similar to previous best
    ## Run 310 stress 8.485223e-05 
    ## ... Procrustes: rmse 0.0001475592  max resid 0.0002310782 
    ## ... Similar to previous best
    ## Run 311 stress 9.780982e-05 
    ## ... Procrustes: rmse 0.00019917  max resid 0.0002958152 
    ## ... Similar to previous best
    ## Run 312 stress 9.730652e-05 
    ## ... Procrustes: rmse 0.0001428762  max resid 0.0002620175 
    ## ... Similar to previous best
    ## Run 313 stress 9.980899e-05 
    ## ... Procrustes: rmse 0.0001342412  max resid 0.0002440623 
    ## ... Similar to previous best
    ## Run 314 stress 9.69404e-05 
    ## ... Procrustes: rmse 0.0001492867  max resid 0.0002733527 
    ## ... Similar to previous best
    ## Run 315 stress 9.770802e-05 
    ## ... Procrustes: rmse 0.0001522435  max resid 0.0002765746 
    ## ... Similar to previous best
    ## Run 316 stress 9.553714e-05 
    ## ... Procrustes: rmse 0.0001862559  max resid 0.0003470719 
    ## ... Similar to previous best
    ## Run 317 stress 9.018647e-05 
    ## ... Procrustes: rmse 0.0001315379  max resid 0.0002543982 
    ## ... Similar to previous best
    ## Run 318 stress 9.315816e-05 
    ## ... Procrustes: rmse 0.0001699806  max resid 0.0003654024 
    ## ... Similar to previous best
    ## Run 319 stress 9.651934e-05 
    ## ... Procrustes: rmse 0.0001514258  max resid 0.0002763367 
    ## ... Similar to previous best
    ## Run 320 stress 9.471337e-05 
    ## ... Procrustes: rmse 0.0001499228  max resid 0.0003107113 
    ## ... Similar to previous best
    ## Run 321 stress 9.613266e-05 
    ## ... Procrustes: rmse 0.0001020563  max resid 0.0001545811 
    ## ... Similar to previous best
    ## Run 322 stress 9.112494e-05 
    ## ... Procrustes: rmse 0.0001669026  max resid 0.0003494941 
    ## ... Similar to previous best
    ## Run 323 stress 9.045292e-05 
    ## ... Procrustes: rmse 0.0001137307  max resid 0.0002314649 
    ## ... Similar to previous best
    ## Run 324 stress 9.206299e-05 
    ## ... Procrustes: rmse 0.0001792716  max resid 0.0003368964 
    ## ... Similar to previous best
    ## Run 325 stress 8.908511e-05 
    ## ... Procrustes: rmse 9.232594e-05  max resid 0.0001626901 
    ## ... Similar to previous best
    ## Run 326 stress 8.992121e-05 
    ## ... Procrustes: rmse 9.305375e-05  max resid 0.0001395326 
    ## ... Similar to previous best
    ## Run 327 stress 9.624318e-05 
    ## ... Procrustes: rmse 0.0001523691  max resid 0.0002581045 
    ## ... Similar to previous best
    ## Run 328 stress 9.511153e-05 
    ## ... Procrustes: rmse 0.0001745797  max resid 0.0003586484 
    ## ... Similar to previous best
    ## Run 329 stress 9.803675e-05 
    ## ... Procrustes: rmse 0.0001888099  max resid 0.0003406281 
    ## ... Similar to previous best
    ## Run 330 stress 9.363568e-05 
    ## ... Procrustes: rmse 0.0001841326  max resid 0.0002821311 
    ## ... Similar to previous best
    ## Run 331 stress 7.48943e-05 
    ## ... Procrustes: rmse 0.0001163395  max resid 0.0002414991 
    ## ... Similar to previous best
    ## Run 332 stress 8.950124e-05 
    ## ... Procrustes: rmse 0.0001633023  max resid 0.0003339503 
    ## ... Similar to previous best
    ## Run 333 stress 9.843931e-05 
    ## ... Procrustes: rmse 0.0001318138  max resid 0.0002405074 
    ## ... Similar to previous best
    ## Run 334 stress 9.76966e-05 
    ## ... Procrustes: rmse 0.0001500482  max resid 0.0002740573 
    ## ... Similar to previous best
    ## Run 335 stress 8.969384e-05 
    ## ... Procrustes: rmse 0.0001661566  max resid 0.0003486245 
    ## ... Similar to previous best
    ## Run 336 stress 9.422487e-05 
    ## ... Procrustes: rmse 8.464039e-05  max resid 0.0001295658 
    ## ... Similar to previous best
    ## Run 337 stress 9.793995e-05 
    ## ... Procrustes: rmse 0.0001531329  max resid 0.0002779835 
    ## ... Similar to previous best
    ## Run 338 stress 9.825253e-05 
    ## ... Procrustes: rmse 0.0002028994  max resid 0.0003811389 
    ## ... Similar to previous best
    ## Run 339 stress 9.439782e-05 
    ## ... Procrustes: rmse 9.688803e-05  max resid 0.0001571139 
    ## ... Similar to previous best
    ## Run 340 stress 9.501166e-05 
    ## ... Procrustes: rmse 0.0001776366  max resid 0.0003640869 
    ## ... Similar to previous best
    ## Run 341 stress 8.502802e-05 
    ## ... Procrustes: rmse 0.0001446818  max resid 0.0002523361 
    ## ... Similar to previous best
    ## Run 342 stress 9.01098e-05 
    ## ... Procrustes: rmse 0.0001463841  max resid 0.0002902603 
    ## ... Similar to previous best
    ## Run 343 stress 9.830416e-05 
    ## ... Procrustes: rmse 0.0001352612  max resid 0.0002567091 
    ## ... Similar to previous best
    ## Run 344 stress 9.413192e-05 
    ## ... Procrustes: rmse 0.0001654532  max resid 0.0003403918 
    ## ... Similar to previous best
    ## Run 345 stress 9.819182e-05 
    ## ... Procrustes: rmse 0.0001322838  max resid 0.0002404243 
    ## ... Similar to previous best
    ## Run 346 stress 9.905287e-05 
    ## ... Procrustes: rmse 0.0001525059  max resid 0.0002766767 
    ## ... Similar to previous best
    ## Run 347 stress 9.813204e-05 
    ## ... Procrustes: rmse 0.0001650282  max resid 0.0002513975 
    ## ... Similar to previous best
    ## Run 348 stress 9.612309e-05 
    ## ... Procrustes: rmse 0.0001918891  max resid 0.0002891265 
    ## ... Similar to previous best
    ## Run 349 stress 9.624937e-05 
    ## ... Procrustes: rmse 0.0001491213  max resid 0.0002731487 
    ## ... Similar to previous best
    ## Run 350 stress 9.805411e-05 
    ## ... Procrustes: rmse 0.00016649  max resid 0.0003149332 
    ## ... Similar to previous best
    ## Run 351 stress 9.686383e-05 
    ## ... Procrustes: rmse 0.0001822441  max resid 0.0003677485 
    ## ... Similar to previous best
    ## Run 352 stress 9.905272e-05 
    ## ... Procrustes: rmse 0.0001786113  max resid 0.0003090179 
    ## ... Similar to previous best
    ## Run 353 stress 9.776212e-05 
    ## ... Procrustes: rmse 0.0001691598  max resid 0.0002760888 
    ## ... Similar to previous best
    ## Run 354 stress 9.743726e-05 
    ## ... Procrustes: rmse 0.0001705069  max resid 0.0002592244 
    ## ... Similar to previous best
    ## Run 355 stress 9.883466e-05 
    ## ... Procrustes: rmse 0.000154527  max resid 0.0002776809 
    ## ... Similar to previous best
    ## Run 356 stress 9.868668e-05 
    ## ... Procrustes: rmse 0.0001606214  max resid 0.000267253 
    ## ... Similar to previous best
    ## Run 357 stress 9.664754e-05 
    ## ... Procrustes: rmse 0.0001950252  max resid 0.0002923363 
    ## ... Similar to previous best
    ## Run 358 stress 9.441982e-05 
    ## ... Procrustes: rmse 0.0001909579  max resid 0.0002884686 
    ## ... Similar to previous best
    ## Run 359 stress 9.164781e-05 
    ## ... Procrustes: rmse 0.0001688178  max resid 0.000364436 
    ## ... Similar to previous best
    ## Run 360 stress 9.697221e-05 
    ## ... Procrustes: rmse 0.0001710679  max resid 0.0002618052 
    ## ... Similar to previous best
    ## Run 361 stress 8.400449e-05 
    ## ... Procrustes: rmse 0.0001484782  max resid 0.0003326022 
    ## ... Similar to previous best
    ## Run 362 stress 9.384684e-05 
    ## ... Procrustes: rmse 0.0001447705  max resid 0.0002693059 
    ## ... Similar to previous best
    ## Run 363 stress 8.798286e-05 
    ## ... Procrustes: rmse 0.0001328343  max resid 0.0002042046 
    ## ... Similar to previous best
    ## Run 364 stress 9.503895e-05 
    ## ... Procrustes: rmse 0.0002005152  max resid 0.0002964694 
    ## ... Similar to previous best
    ## Run 365 stress 9.822324e-05 
    ## ... Procrustes: rmse 7.466165e-05  max resid 0.000106183 
    ## ... Similar to previous best
    ## Run 366 stress 9.624805e-05 
    ## ... Procrustes: rmse 0.0001627948  max resid 0.0003525959 
    ## ... Similar to previous best
    ## Run 367 stress 9.925185e-05 
    ## ... Procrustes: rmse 0.0001960982  max resid 0.0003373718 
    ## ... Similar to previous best
    ## Run 368 stress 9.986875e-05 
    ## ... Procrustes: rmse 0.0001143282  max resid 0.0001733327 
    ## ... Similar to previous best
    ## Run 369 stress 9.917757e-05 
    ## ... Procrustes: rmse 0.0001565563  max resid 0.0002801989 
    ## ... Similar to previous best
    ## Run 370 stress 9.826736e-05 
    ## ... Procrustes: rmse 0.0001847694  max resid 0.0003699358 
    ## ... Similar to previous best
    ## Run 371 stress 9.694043e-05 
    ## ... Procrustes: rmse 0.0001913618  max resid 0.0003291247 
    ## ... Similar to previous best
    ## Run 372 stress 9.504807e-05 
    ## ... Procrustes: rmse 0.0001955749  max resid 0.000371941 
    ## ... Similar to previous best
    ## Run 373 stress 8.673078e-05 
    ## ... Procrustes: rmse 0.0001023088  max resid 0.0001873137 
    ## ... Similar to previous best
    ## Run 374 stress 9.378212e-05 
    ## ... Procrustes: rmse 9.755555e-05  max resid 0.0001417876 
    ## ... Similar to previous best
    ## Run 375 stress 8.943003e-05 
    ## ... Procrustes: rmse 0.0001321618  max resid 0.000250379 
    ## ... Similar to previous best
    ## Run 376 stress 9.809191e-05 
    ## ... Procrustes: rmse 0.000168495  max resid 0.0002752881 
    ## ... Similar to previous best
    ## Run 377 stress 9.717957e-05 
    ## ... Procrustes: rmse 0.0001753169  max resid 0.0003074787 
    ## ... Similar to previous best
    ## Run 378 stress 9.864838e-05 
    ## ... Procrustes: rmse 0.0001802011  max resid 0.0003594407 
    ## ... Similar to previous best
    ## Run 379 stress 9.567999e-05 
    ## ... Procrustes: rmse 0.000186777  max resid 0.000342219 
    ## ... Similar to previous best
    ## Run 380 stress 9.42293e-05 
    ## ... Procrustes: rmse 0.0001645172  max resid 0.0002721674 
    ## ... Similar to previous best
    ## Run 381 stress 9.700042e-05 
    ## ... Procrustes: rmse 0.0001033029  max resid 0.0001525693 
    ## ... Similar to previous best
    ## Run 382 stress 9.911862e-05 
    ## ... Procrustes: rmse 0.0001608782  max resid 0.0002479728 
    ## ... Similar to previous best
    ## Run 383 stress 9.796912e-05 
    ## ... Procrustes: rmse 0.0001910392  max resid 0.0003538538 
    ## ... Similar to previous best
    ## Run 384 stress 9.570156e-05 
    ## ... Procrustes: rmse 9.654293e-05  max resid 0.0001420668 
    ## ... Similar to previous best
    ## Run 385 stress 9.587324e-05 
    ## ... Procrustes: rmse 0.0001649562  max resid 0.0002762871 
    ## ... Similar to previous best
    ## Run 386 stress 9.733765e-05 
    ## ... Procrustes: rmse 0.000195255  max resid 0.0002928635 
    ## ... Similar to previous best
    ## Run 387 stress 9.882164e-05 
    ## ... Procrustes: rmse 0.0001100489  max resid 0.0001593241 
    ## ... Similar to previous best
    ## Run 388 stress 9.407659e-05 
    ## ... Procrustes: rmse 0.000198133  max resid 0.0002958203 
    ## ... Similar to previous best
    ## Run 389 stress 9.11316e-05 
    ## ... Procrustes: rmse 8.970256e-05  max resid 0.0001366107 
    ## ... Similar to previous best
    ## Run 390 stress 8.580264e-05 
    ## ... Procrustes: rmse 0.0001509522  max resid 0.0002807623 
    ## ... Similar to previous best
    ## Run 391 stress 9.870301e-05 
    ## ... Procrustes: rmse 0.0001251478  max resid 0.0002410156 
    ## ... Similar to previous best
    ## Run 392 stress 9.419878e-05 
    ## ... Procrustes: rmse 0.0001997545  max resid 0.0002960489 
    ## ... Similar to previous best
    ## Run 393 stress 9.887862e-05 
    ## ... Procrustes: rmse 0.0001038516  max resid 0.0001446036 
    ## ... Similar to previous best
    ## Run 394 stress 8.542116e-05 
    ## ... Procrustes: rmse 8.351683e-05  max resid 0.0001281375 
    ## ... Similar to previous best
    ## Run 395 stress 8.717834e-05 
    ## ... Procrustes: rmse 0.0001760778  max resid 0.0002616472 
    ## ... Similar to previous best
    ## Run 396 stress 8.688265e-05 
    ## ... Procrustes: rmse 0.0001347893  max resid 0.0002543381 
    ## ... Similar to previous best
    ## Run 397 stress 9.708103e-05 
    ## ... Procrustes: rmse 0.0001597475  max resid 0.0003507642 
    ## ... Similar to previous best
    ## Run 398 stress 9.358958e-05 
    ## ... Procrustes: rmse 0.0001000043  max resid 0.0002071104 
    ## ... Similar to previous best
    ## Run 399 stress 9.452254e-05 
    ## ... Procrustes: rmse 0.0001184027  max resid 0.0002285525 
    ## ... Similar to previous best
    ## Run 400 stress 8.793423e-05 
    ## ... Procrustes: rmse 0.0001513596  max resid 0.0003367972 
    ## ... Similar to previous best
    ## Run 401 stress 9.472569e-05 
    ## ... Procrustes: rmse 0.0001726924  max resid 0.0003474126 
    ## ... Similar to previous best
    ## Run 402 stress 9.515508e-05 
    ## ... Procrustes: rmse 0.0001170163  max resid 0.0002366519 
    ## ... Similar to previous best
    ## Run 403 stress 9.865916e-05 
    ## ... Procrustes: rmse 0.0001848423  max resid 0.000369275 
    ## ... Similar to previous best
    ## Run 404 stress 7.10305e-05 
    ## ... Procrustes: rmse 8.777525e-05  max resid 0.0001546481 
    ## ... Similar to previous best
    ## Run 405 stress 9.822223e-05 
    ## ... Procrustes: rmse 0.0001494399  max resid 0.0002555683 
    ## ... Similar to previous best
    ## Run 406 stress 9.800177e-05 
    ## ... Procrustes: rmse 0.000120273  max resid 0.000231241 
    ## ... Similar to previous best
    ## Run 407 stress 9.506264e-05 
    ## ... Procrustes: rmse 0.0001594309  max resid 0.0003347139 
    ## ... Similar to previous best
    ## Run 408 stress 9.910902e-05 
    ## ... Procrustes: rmse 0.0001850974  max resid 0.0003690414 
    ## ... Similar to previous best
    ## Run 409 stress 9.280748e-05 
    ## ... Procrustes: rmse 0.0001498595  max resid 0.00026692 
    ## ... Similar to previous best
    ## Run 410 stress 9.766282e-05 
    ## ... Procrustes: rmse 0.0001961845  max resid 0.0002920987 
    ## ... Similar to previous best
    ## Run 411 stress 9.939858e-05 
    ## ... Procrustes: rmse 0.0001813632  max resid 0.0003161465 
    ## ... Similar to previous best
    ## Run 412 stress 9.581516e-05 
    ## ... Procrustes: rmse 0.0001855735  max resid 0.0002817908 
    ## ... Similar to previous best
    ## Run 413 stress 0.2362213 
    ## Run 414 stress 9.645795e-05 
    ## ... Procrustes: rmse 0.0001585129  max resid 0.0002806322 
    ## ... Similar to previous best
    ## Run 415 stress 8.480178e-05 
    ## ... Procrustes: rmse 8.156427e-05  max resid 0.0001247451 
    ## ... Similar to previous best
    ## Run 416 stress 9.799322e-05 
    ## ... Procrustes: rmse 0.0001093786  max resid 0.0001693052 
    ## ... Similar to previous best
    ## Run 417 stress 8.796355e-05 
    ## ... Procrustes: rmse 0.000134322  max resid 0.0002376343 
    ## ... Similar to previous best
    ## Run 418 stress 8.631917e-05 
    ## ... Procrustes: rmse 8.374087e-05  max resid 0.0001267861 
    ## ... Similar to previous best
    ## Run 419 stress 9.834562e-05 
    ## ... Procrustes: rmse 0.0001061134  max resid 0.0001601638 
    ## ... Similar to previous best
    ## Run 420 stress 9.373389e-05 
    ## ... Procrustes: rmse 0.0001895593  max resid 0.0002872753 
    ## ... Similar to previous best
    ## Run 421 stress 9.642699e-05 
    ## ... Procrustes: rmse 0.0001155508  max resid 0.000182891 
    ## ... Similar to previous best
    ## Run 422 stress 9.925786e-05 
    ## ... Procrustes: rmse 0.0001510396  max resid 0.0002756375 
    ## ... Similar to previous best
    ## Run 423 stress 8.947795e-05 
    ## ... Procrustes: rmse 0.0001590322  max resid 0.0002903088 
    ## ... Similar to previous best
    ## Run 424 stress 9.946246e-05 
    ## ... Procrustes: rmse 0.0001939131  max resid 0.0003563844 
    ## ... Similar to previous best
    ## Run 425 stress 9.740464e-05 
    ## ... Procrustes: rmse 0.0001709905  max resid 0.0002622478 
    ## ... Similar to previous best
    ## Run 426 stress 9.696859e-05 
    ## ... Procrustes: rmse 0.0001626299  max resid 0.0002534797 
    ## ... Similar to previous best
    ## Run 427 stress 9.568223e-05 
    ## ... Procrustes: rmse 0.0001964929  max resid 0.0003728964 
    ## ... Similar to previous best
    ## Run 428 stress 9.298575e-05 
    ## ... Procrustes: rmse 0.0001334903  max resid 0.0002572887 
    ## ... Similar to previous best
    ## Run 429 stress 9.25725e-05 
    ## ... Procrustes: rmse 0.0001732927  max resid 0.0003212942 
    ## ... Similar to previous best
    ## Run 430 stress 9.198359e-05 
    ## ... Procrustes: rmse 6.405824e-05  max resid 9.037074e-05 
    ## ... Similar to previous best
    ## Run 431 stress 9.013452e-05 
    ## ... Procrustes: rmse 9.32145e-05  max resid 0.0001415553 
    ## ... Similar to previous best
    ## Run 432 stress 9.849796e-05 
    ## ... Procrustes: rmse 0.0001532739  max resid 0.0002793845 
    ## ... Similar to previous best
    ## Run 433 stress 9.476742e-05 
    ## ... Procrustes: rmse 0.0001732635  max resid 0.0003667347 
    ## ... Similar to previous best
    ## Run 434 stress 9.773464e-05 
    ## ... Procrustes: rmse 0.0001426454  max resid 0.0002253098 
    ## ... Similar to previous best
    ## Run 435 stress 9.94249e-05 
    ## ... Procrustes: rmse 0.0001659035  max resid 0.0002923389 
    ## ... Similar to previous best
    ## Run 436 stress 9.567826e-05 
    ## ... Procrustes: rmse 0.0001055062  max resid 0.000159933 
    ## ... Similar to previous best
    ## Run 437 stress 9.373097e-05 
    ## ... Procrustes: rmse 0.0001632222  max resid 0.0002698633 
    ## ... Similar to previous best
    ## Run 438 stress 9.743117e-05 
    ## ... Procrustes: rmse 9.750401e-05  max resid 0.0001495648 
    ## ... Similar to previous best
    ## Run 439 stress 8.790243e-05 
    ## ... Procrustes: rmse 0.0001753627  max resid 0.0002742735 
    ## ... Similar to previous best
    ## Run 440 stress 9.93641e-05 
    ## ... Procrustes: rmse 0.000106405  max resid 0.0001726391 
    ## ... Similar to previous best
    ## Run 441 stress 9.659795e-05 
    ## ... Procrustes: rmse 0.0001346718  max resid 0.0002596702 
    ## ... Similar to previous best
    ## Run 442 stress 9.668488e-05 
    ## ... Procrustes: rmse 0.0001078822  max resid 0.0001622052 
    ## ... Similar to previous best
    ## Run 443 stress 9.794645e-05 
    ## ... Procrustes: rmse 0.0001995404  max resid 0.000296673 
    ## ... Similar to previous best
    ## Run 444 stress 8.741382e-05 
    ## ... Procrustes: rmse 8.707626e-05  max resid 0.0001278199 
    ## ... Similar to previous best
    ## Run 445 stress 9.165165e-05 
    ## ... Procrustes: rmse 0.0001717336  max resid 0.0003201709 
    ## ... Similar to previous best
    ## Run 446 stress 9.24761e-05 
    ## ... Procrustes: rmse 8.181406e-05  max resid 0.0001379232 
    ## ... Similar to previous best
    ## Run 447 stress 8.783011e-05 
    ## ... Procrustes: rmse 0.0001659653  max resid 0.0003437061 
    ## ... Similar to previous best
    ## Run 448 stress 9.151201e-05 
    ## ... Procrustes: rmse 0.0001233303  max resid 0.0002064652 
    ## ... Similar to previous best
    ## Run 449 stress 9.71895e-05 
    ## ... Procrustes: rmse 0.0002063233  max resid 0.0003029699 
    ## ... Similar to previous best
    ## Run 450 stress 9.760102e-05 
    ## ... Procrustes: rmse 0.0001846457  max resid 0.0002824626 
    ## ... Similar to previous best
    ## Run 451 stress 9.540613e-05 
    ## ... Procrustes: rmse 0.0001871608  max resid 0.0003226364 
    ## ... Similar to previous best
    ## Run 452 stress 9.858982e-05 
    ## ... Procrustes: rmse 0.0002012945  max resid 0.0002981086 
    ## ... Similar to previous best
    ## Run 453 stress 9.696783e-05 
    ## ... Procrustes: rmse 0.0001646958  max resid 0.0002494039 
    ## ... Similar to previous best
    ## Run 454 stress 9.989652e-05 
    ## ... Procrustes: rmse 0.0001881928  max resid 0.0003741127 
    ## ... Similar to previous best
    ## Run 455 stress 9.805507e-05 
    ## ... Procrustes: rmse 0.0001946213  max resid 0.0002919503 
    ## ... Similar to previous best
    ## Run 456 stress 9.536315e-05 
    ## ... Procrustes: rmse 0.0001025274  max resid 0.0001575407 
    ## ... Similar to previous best
    ## Run 457 stress 9.435759e-05 
    ## ... Procrustes: rmse 0.0001621932  max resid 0.0003090691 
    ## ... Similar to previous best
    ## Run 458 stress 9.350465e-05 
    ## ... Procrustes: rmse 9.059501e-05  max resid 0.0001514276 
    ## ... Similar to previous best
    ## Run 459 stress 8.939812e-05 
    ## ... Procrustes: rmse 0.0001632888  max resid 0.0003487809 
    ## ... Similar to previous best
    ## Run 460 stress 9.409956e-05 
    ## ... Procrustes: rmse 0.0001733467  max resid 0.0003576799 
    ## ... Similar to previous best
    ## Run 461 stress 8.167978e-05 
    ## ... Procrustes: rmse 0.0001125672  max resid 0.0002340905 
    ## ... Similar to previous best
    ## Run 462 stress 9.694361e-05 
    ## ... Procrustes: rmse 0.000190256  max resid 0.0002830261 
    ## ... Similar to previous best
    ## Run 463 stress 9.86394e-05 
    ## ... Procrustes: rmse 0.0002026061  max resid 0.0002957111 
    ## ... Similar to previous best
    ## Run 464 stress 9.975182e-05 
    ## ... Procrustes: rmse 0.0001917343  max resid 0.0003796609 
    ## ... Similar to previous best
    ## Run 465 stress 9.385668e-05 
    ## ... Procrustes: rmse 0.0001835248  max resid 0.0002821252 
    ## ... Similar to previous best
    ## Run 466 stress 9.654889e-05 
    ## ... Procrustes: rmse 0.0001575947  max resid 0.0002893537 
    ## ... Similar to previous best
    ## Run 467 stress 8.982122e-05 
    ## ... Procrustes: rmse 7.778477e-05  max resid 0.0001086201 
    ## ... Similar to previous best
    ## Run 468 stress 9.591744e-05 
    ## ... Procrustes: rmse 0.0001651333  max resid 0.0003110334 
    ## ... Similar to previous best
    ## Run 469 stress 0.2307363 
    ## Run 470 stress 9.793546e-05 
    ## ... Procrustes: rmse 0.0001727516  max resid 0.0002466243 
    ## ... Similar to previous best
    ## Run 471 stress 9.333342e-05 
    ## ... Procrustes: rmse 0.0001853189  max resid 0.0002808493 
    ## ... Similar to previous best
    ## Run 472 stress 9.726136e-05 
    ## ... Procrustes: rmse 0.00017205  max resid 0.0002632257 
    ## ... Similar to previous best
    ## Run 473 stress 9.896106e-05 
    ## ... Procrustes: rmse 0.0001837822  max resid 0.0003719591 
    ## ... Similar to previous best
    ## Run 474 stress 8.841073e-05 
    ## ... Procrustes: rmse 8.893957e-05  max resid 0.0001332594 
    ## ... Similar to previous best
    ## Run 475 stress 9.283567e-05 
    ## ... Procrustes: rmse 0.0001439473  max resid 0.0002477813 
    ## ... Similar to previous best
    ## Run 476 stress 9.067093e-05 
    ## ... Procrustes: rmse 0.0001820375  max resid 0.0003595737 
    ## ... Similar to previous best
    ## Run 477 stress 9.477857e-05 
    ## ... Procrustes: rmse 0.0001702094  max resid 0.0003185229 
    ## ... Similar to previous best
    ## Run 478 stress 9.466484e-05 
    ## ... Procrustes: rmse 0.0001617087  max resid 0.0002691271 
    ## ... Similar to previous best
    ## Run 479 stress 9.542628e-05 
    ## ... Procrustes: rmse 0.0001374097  max resid 0.0002554112 
    ## ... Similar to previous best
    ## Run 480 stress 8.840111e-05 
    ## ... Procrustes: rmse 0.0001289241  max resid 0.0002263666 
    ## ... Similar to previous best
    ## Run 481 stress 9.999802e-05 
    ## ... Procrustes: rmse 9.666655e-05  max resid 0.000128653 
    ## ... Similar to previous best
    ## Run 482 stress 9.455056e-05 
    ## ... Procrustes: rmse 0.0001414197  max resid 0.0002661302 
    ## ... Similar to previous best
    ## Run 483 stress 8.646646e-05 
    ## ... Procrustes: rmse 7.318138e-05  max resid 0.000112135 
    ## ... Similar to previous best
    ## Run 484 stress 9.376097e-05 
    ## ... Procrustes: rmse 8.782473e-05  max resid 0.0001214203 
    ## ... Similar to previous best
    ## Run 485 stress 9.648024e-05 
    ## ... Procrustes: rmse 0.0001791515  max resid 0.0003645108 
    ## ... Similar to previous best
    ## Run 486 stress 9.561019e-05 
    ## ... Procrustes: rmse 0.0001668056  max resid 0.0002733832 
    ## ... Similar to previous best
    ## Run 487 stress 9.374264e-05 
    ## ... Procrustes: rmse 0.0001114524  max resid 0.0002132286 
    ## ... Similar to previous best
    ## Run 488 stress 9.340169e-05 
    ## ... Procrustes: rmse 0.0001749595  max resid 0.0003220801 
    ## ... Similar to previous best
    ## Run 489 stress 9.507068e-05 
    ## ... Procrustes: rmse 0.000201092  max resid 0.0002981301 
    ## ... Similar to previous best
    ## Run 490 stress 9.827768e-05 
    ## ... Procrustes: rmse 0.000180532  max resid 0.0003530872 
    ## ... Similar to previous best
    ## Run 491 stress 9.805955e-05 
    ## ... Procrustes: rmse 0.0001841507  max resid 0.0003682389 
    ## ... Similar to previous best
    ## Run 492 stress 9.671964e-05 
    ## ... Procrustes: rmse 0.0001584293  max resid 0.0002753415 
    ## ... Similar to previous best
    ## Run 493 stress 9.291128e-05 
    ## ... Procrustes: rmse 0.0001423624  max resid 0.0002680019 
    ## ... Similar to previous best
    ## Run 494 stress 9.657442e-05 
    ## ... Procrustes: rmse 0.0001888026  max resid 0.0003532323 
    ## ... Similar to previous best
    ## Run 495 stress 9.098805e-05 
    ## ... Procrustes: rmse 0.0001667946  max resid 0.0003609218 
    ## ... Similar to previous best
    ## Run 496 stress 9.365533e-05 
    ## ... Procrustes: rmse 0.0001840282  max resid 0.0002816126 
    ## ... Similar to previous best
    ## Run 497 stress 7.702302e-05 
    ## ... Procrustes: rmse 0.0001301232  max resid 0.0002328202 
    ## ... Similar to previous best
    ## Run 498 stress 9.51923e-05 
    ## ... Procrustes: rmse 0.0001787679  max resid 0.0003309244 
    ## ... Similar to previous best
    ## Run 499 stress 9.395884e-05 
    ## ... Procrustes: rmse 0.0001651712  max resid 0.0003411318 
    ## ... Similar to previous best
    ## Run 500 stress 9.789302e-05 
    ## ... Procrustes: rmse 0.0001775938  max resid 0.0003540753 
    ## ... Similar to previous best
    ## *** Best solution repeated 485 times

    ## Warning in metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Mixed lakes
SD_beta_env_M_NMDS <- metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.001279543 
    ## Run 1 stress 0.2842805 
    ## Run 2 stress 0.1990774 
    ## Run 3 stress 0.1776015 
    ## Run 4 stress 0.001409997 
    ## ... Procrustes: rmse 0.002164193  max resid 0.003001915 
    ## ... Similar to previous best
    ## Run 5 stress 0.001223744 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0009758979  max resid 0.001354174 
    ## ... Similar to previous best
    ## Run 6 stress 0.0008072729 
    ## ... New best solution
    ## ... Procrustes: rmse 0.007977041  max resid 0.01104417 
    ## Run 7 stress 0.001281659 
    ## ... Procrustes: rmse 0.008986766  max resid 0.01246101 
    ## Run 8 stress 0.001336399 
    ## Run 9 stress 9.727294e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03430608  max resid 0.04710037 
    ## Run 10 stress 0.0001467305 
    ## ... Procrustes: rmse 0.008775868  max resid 0.01195483 
    ## Run 11 stress 0.2373687 
    ## Run 12 stress 9.913329e-05 
    ## ... Procrustes: rmse 0.0003270031  max resid 0.0006522846 
    ## ... Similar to previous best
    ## Run 13 stress 0.2440272 
    ## Run 14 stress 9.126881e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001332339  max resid 0.0002200589 
    ## ... Similar to previous best
    ## Run 15 stress 0.0004488074 
    ## ... Procrustes: rmse 0.01552889  max resid 0.02114723 
    ## Run 16 stress 9.716506e-05 
    ## ... Procrustes: rmse 2.189263e-05  max resid 3.311202e-05 
    ## ... Similar to previous best
    ## Run 17 stress 0.1491957 
    ## Run 18 stress 0.001244237 
    ## Run 19 stress 0.2570118 
    ## Run 20 stress 0.1776016 
    ## Run 21 stress 0.0002627008 
    ## ... Procrustes: rmse 0.01186997  max resid 0.01609783 
    ## Run 22 stress 0.2528294 
    ## Run 23 stress 8.535957e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 6.797152e-05  max resid 0.0001311498 
    ## ... Similar to previous best
    ## Run 24 stress 6.35848e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003973411  max resid 0.0007577442 
    ## ... Similar to previous best
    ## Run 25 stress 0.001377978 
    ## Run 26 stress 0.177601 
    ## Run 27 stress 0.000291301 
    ## ... Procrustes: rmse 0.02054366  max resid 0.02907201 
    ## Run 28 stress 9.292755e-05 
    ## ... Procrustes: rmse 0.0003193904  max resid 0.0005127447 
    ## ... Similar to previous best
    ## Run 29 stress 0.001256674 
    ## Run 30 stress 8.926218e-05 
    ## ... Procrustes: rmse 0.0004155936  max resid 0.0007728275 
    ## ... Similar to previous best
    ## Run 31 stress 0.1491957 
    ## Run 32 stress 0.0007363164 
    ## Run 33 stress 0.001293895 
    ## Run 34 stress 0.2848524 
    ## Run 35 stress 0.001318415 
    ## Run 36 stress 0.001342449 
    ## Run 37 stress 0.001153073 
    ## Run 38 stress 0.2797325 
    ## Run 39 stress 9.473623e-05 
    ## ... Procrustes: rmse 0.000366354  max resid 0.0007647748 
    ## ... Similar to previous best
    ## Run 40 stress 9.159509e-05 
    ## ... Procrustes: rmse 0.0003536579  max resid 0.0007292846 
    ## ... Similar to previous best
    ## Run 41 stress 0.1491957 
    ## Run 42 stress 0.0002710317 
    ## ... Procrustes: rmse 0.01169194  max resid 0.01610863 
    ## Run 43 stress 8.963699e-05 
    ## ... Procrustes: rmse 0.0003489526  max resid 0.0007116498 
    ## ... Similar to previous best
    ## Run 44 stress 8.267266e-05 
    ## ... Procrustes: rmse 0.0004128155  max resid 0.0007567826 
    ## ... Similar to previous best
    ## Run 45 stress 0.0005052248 
    ## ... Procrustes: rmse 0.01610341  max resid 0.02219551 
    ## Run 46 stress 7.012661e-05 
    ## ... Procrustes: rmse 0.0003946102  max resid 0.0007003885 
    ## ... Similar to previous best
    ## Run 47 stress 0.2846414 
    ## Run 48 stress 9.310177e-05 
    ## ... Procrustes: rmse 0.0003379199  max resid 0.0005971217 
    ## ... Similar to previous best
    ## Run 49 stress 0.1990774 
    ## Run 50 stress 0.0005013673 
    ## ... Procrustes: rmse 0.01601375  max resid 0.02207202 
    ## Run 51 stress 0.001236049 
    ## Run 52 stress 0.001317905 
    ## Run 53 stress 0.001291433 
    ## Run 54 stress 0.001384294 
    ## Run 55 stress 0.1491957 
    ## Run 56 stress 8.40947e-05 
    ## ... Procrustes: rmse 0.000122657  max resid 0.0001728445 
    ## ... Similar to previous best
    ## Run 57 stress 9.219582e-05 
    ## ... Procrustes: rmse 0.0003768513  max resid 0.0005688966 
    ## ... Similar to previous best
    ## Run 58 stress 0.001221696 
    ## Run 59 stress 0.001225508 
    ## Run 60 stress 8.2806e-05 
    ## ... Procrustes: rmse 0.0003365182  max resid 0.0005629733 
    ## ... Similar to previous best
    ## Run 61 stress 0.001282568 
    ## Run 62 stress 0.001291201 
    ## Run 63 stress 0.1990774 
    ## Run 64 stress 0.00143095 
    ## Run 65 stress 9.624835e-05 
    ## ... Procrustes: rmse 0.004796238  max resid 0.006617524 
    ## ... Similar to previous best
    ## Run 66 stress 9.232946e-05 
    ## ... Procrustes: rmse 0.0003565175  max resid 0.0004711877 
    ## ... Similar to previous best
    ## Run 67 stress 8.911158e-05 
    ## ... Procrustes: rmse 0.0003843086  max resid 0.0005509902 
    ## ... Similar to previous best
    ## Run 68 stress 9.541624e-05 
    ## ... Procrustes: rmse 0.0004009893  max resid 0.0007602606 
    ## ... Similar to previous best
    ## Run 69 stress 8.546782e-05 
    ## ... Procrustes: rmse 0.0001315918  max resid 0.0001815375 
    ## ... Similar to previous best
    ## Run 70 stress 8.891072e-05 
    ## ... Procrustes: rmse 0.0004032796  max resid 0.0005987993 
    ## ... Similar to previous best
    ## Run 71 stress 0.001131735 
    ## Run 72 stress 8.8327e-05 
    ## ... Procrustes: rmse 0.0004180124  max resid 0.0006112094 
    ## ... Similar to previous best
    ## Run 73 stress 9.6488e-05 
    ## ... Procrustes: rmse 0.0003492915  max resid 0.0006401847 
    ## ... Similar to previous best
    ## Run 74 stress 9.192478e-05 
    ## ... Procrustes: rmse 0.000361071  max resid 0.000732505 
    ## ... Similar to previous best
    ## Run 75 stress 0.0004733297 
    ## ... Procrustes: rmse 0.01557544  max resid 0.02146687 
    ## Run 76 stress 0.00137287 
    ## Run 77 stress 8.624956e-05 
    ## ... Procrustes: rmse 0.000375401  max resid 0.0005824715 
    ## ... Similar to previous best
    ## Run 78 stress 0.001227137 
    ## Run 79 stress 0.001272732 
    ## Run 80 stress 9.438478e-05 
    ## ... Procrustes: rmse 0.0004072622  max resid 0.0006780625 
    ## ... Similar to previous best
    ## Run 81 stress 0.001419993 
    ## Run 82 stress 7.482961e-05 
    ## ... Procrustes: rmse 0.0003287369  max resid 0.0005915774 
    ## ... Similar to previous best
    ## Run 83 stress 0.0004946903 
    ## ... Procrustes: rmse 0.01593014  max resid 0.0219567 
    ## Run 84 stress 8.915789e-05 
    ## ... Procrustes: rmse 0.0003531883  max resid 0.0005661395 
    ## ... Similar to previous best
    ## Run 85 stress 0.1776012 
    ## Run 86 stress 9.311098e-05 
    ## ... Procrustes: rmse 0.0003854262  max resid 0.0006109867 
    ## ... Similar to previous best
    ## Run 87 stress 8.401862e-05 
    ## ... Procrustes: rmse 0.0003174872  max resid 0.0006024593 
    ## ... Similar to previous best
    ## Run 88 stress 6.627058e-05 
    ## ... Procrustes: rmse 0.001053272  max resid 0.002073468 
    ## ... Similar to previous best
    ## Run 89 stress 9.056588e-05 
    ## ... Procrustes: rmse 0.0004182693  max resid 0.0007771194 
    ## ... Similar to previous best
    ## Run 90 stress 7.094965e-05 
    ## ... Procrustes: rmse 0.0003901272  max resid 0.0006807883 
    ## ... Similar to previous best
    ## Run 91 stress 0.001267328 
    ## Run 92 stress 0.001249268 
    ## Run 93 stress 0.1990774 
    ## Run 94 stress 0.1990774 
    ## Run 95 stress 0.0004869816 
    ## ... Procrustes: rmse 0.0158042  max resid 0.02178257 
    ## Run 96 stress 0.001241854 
    ## Run 97 stress 0.0005000826 
    ## ... Procrustes: rmse 0.01601853  max resid 0.02207895 
    ## Run 98 stress 9.681695e-05 
    ## ... Procrustes: rmse 0.0004208734  max resid 0.0007908608 
    ## ... Similar to previous best
    ## Run 99 stress 0.1990774 
    ## Run 100 stress 0.001377772 
    ## Run 101 stress 9.790241e-05 
    ## ... Procrustes: rmse 0.000420467  max resid 0.0007791527 
    ## ... Similar to previous best
    ## Run 102 stress 0.0004823745 
    ## ... Procrustes: rmse 0.01572645  max resid 0.02167505 
    ## Run 103 stress 0.000460078 
    ## ... Procrustes: rmse 0.01535027  max resid 0.02115645 
    ## Run 104 stress 9.125964e-05 
    ## ... Procrustes: rmse 0.0003636494  max resid 0.0007537927 
    ## ... Similar to previous best
    ## Run 105 stress 0.0008729882 
    ## Run 106 stress 9.964833e-05 
    ## ... Procrustes: rmse 0.0004182876  max resid 0.0007706255 
    ## ... Similar to previous best
    ## Run 107 stress 2.73566e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000367523  max resid 0.0005178491 
    ## ... Similar to previous best
    ## Run 108 stress 9.275637e-05 
    ## ... Procrustes: rmse 0.0001595538  max resid 0.0002285451 
    ## ... Similar to previous best
    ## Run 109 stress 0.1491957 
    ## Run 110 stress 0.1990774 
    ## Run 111 stress 0.001146426 
    ## Run 112 stress 0.0004250854 
    ## ... Procrustes: rmse 0.01509941  max resid 0.02081857 
    ## Run 113 stress 9.583707e-05 
    ## ... Procrustes: rmse 0.0001688349  max resid 0.0002554234 
    ## ... Similar to previous best
    ## Run 114 stress 9.077286e-05 
    ## ... Procrustes: rmse 0.002658642  max resid 0.003683795 
    ## ... Similar to previous best
    ## Run 115 stress 0.00130299 
    ## Run 116 stress 0.001301979 
    ## Run 117 stress 9.930088e-05 
    ## ... Procrustes: rmse 0.0001805159  max resid 0.0002778385 
    ## ... Similar to previous best
    ## Run 118 stress 0.2520602 
    ## Run 119 stress 0.001016616 
    ## Run 120 stress 9.522259e-05 
    ## ... Procrustes: rmse 0.0001715423  max resid 0.0002490565 
    ## ... Similar to previous best
    ## Run 121 stress 0.0004785691 
    ## ... Procrustes: rmse 0.01601521  max resid 0.02208425 
    ## Run 122 stress 0.0002746403 
    ## ... Procrustes: rmse 0.02004023  max resid 0.02761623 
    ## Run 123 stress 9.296065e-05 
    ## ... Procrustes: rmse 0.0001245213  max resid 0.0002055039 
    ## ... Similar to previous best
    ## Run 124 stress 0.001186317 
    ## Run 125 stress 0.001224092 
    ## Run 126 stress 0.1491957 
    ## Run 127 stress 0.0004855976 
    ## ... Procrustes: rmse 0.01613826  max resid 0.02225277 
    ## Run 128 stress 0.1990774 
    ## Run 129 stress 0.177602 
    ## Run 130 stress 0.0004790093 
    ## ... Procrustes: rmse 0.01597309  max resid 0.02202699 
    ## Run 131 stress 0.2440272 
    ## Run 132 stress 0.0004000949 
    ## ... Procrustes: rmse 0.01462926  max resid 0.02017263 
    ## Run 133 stress 0.001412729 
    ## Run 134 stress 0.2852167 
    ## Run 135 stress 9.184696e-05 
    ## ... Procrustes: rmse 0.0001828052  max resid 0.000265438 
    ## ... Similar to previous best
    ## Run 136 stress 0.0004713403 
    ## ... Procrustes: rmse 0.01590122  max resid 0.02192393 
    ## Run 137 stress 0.001271566 
    ## Run 138 stress 9.143021e-05 
    ## ... Procrustes: rmse 0.0001761074  max resid 0.0002263264 
    ## ... Similar to previous best
    ## Run 139 stress 9.14706e-05 
    ## ... Procrustes: rmse 0.0001503752  max resid 0.0002325708 
    ## ... Similar to previous best
    ## Run 140 stress 0.0003696101 
    ## ... Procrustes: rmse 0.02323024  max resid 0.03203073 
    ## Run 141 stress 0.0004435857 
    ## ... Procrustes: rmse 0.01542482  max resid 0.02126783 
    ## Run 142 stress 0.0001863611 
    ## ... Procrustes: rmse 0.009987616  max resid 0.01376505 
    ## Run 143 stress 9.460753e-05 
    ## ... Procrustes: rmse 0.0001716239  max resid 0.0002482196 
    ## ... Similar to previous best
    ## Run 144 stress 0.2806642 
    ## Run 145 stress 9.672985e-05 
    ## ... Procrustes: rmse 0.0001820565  max resid 0.0002629836 
    ## ... Similar to previous best
    ## Run 146 stress 0.2696744 
    ## Run 147 stress 0.1776015 
    ## Run 148 stress 0.1491957 
    ## Run 149 stress 8.667608e-05 
    ## ... Procrustes: rmse 0.0001679678  max resid 0.0002355054 
    ## ... Similar to previous best
    ## Run 150 stress 9.158728e-05 
    ## ... Procrustes: rmse 0.0001744584  max resid 0.0002382894 
    ## ... Similar to previous best
    ## Run 151 stress 0.1491957 
    ## Run 152 stress 0.0006905809 
    ## Run 153 stress 8.095071e-05 
    ## ... Procrustes: rmse 0.001730609  max resid 0.002396886 
    ## ... Similar to previous best
    ## Run 154 stress 8.979141e-05 
    ## ... Procrustes: rmse 0.0001786981  max resid 0.0002592166 
    ## ... Similar to previous best
    ## Run 155 stress 0.001182153 
    ## Run 156 stress 0.001284588 
    ## Run 157 stress 0.1776021 
    ## Run 158 stress 0.177601 
    ## Run 159 stress 0.001169798 
    ## Run 160 stress 0.1776017 
    ## Run 161 stress 0.0004758042 
    ## ... Procrustes: rmse 0.01596959  max resid 0.02201901 
    ## Run 162 stress 0.0004869711 
    ## ... Procrustes: rmse 0.01613521  max resid 0.02225023 
    ## Run 163 stress 0.0008222494 
    ## Run 164 stress 9.357061e-05 
    ## ... Procrustes: rmse 0.000184639  max resid 0.0002740987 
    ## ... Similar to previous best
    ## Run 165 stress 9.011638e-05 
    ## ... Procrustes: rmse 0.0001629228  max resid 0.0002305053 
    ## ... Similar to previous best
    ## Run 166 stress 0.1776011 
    ## Run 167 stress 8.122209e-05 
    ## ... Procrustes: rmse 0.000313554  max resid 0.000447935 
    ## ... Similar to previous best
    ## Run 168 stress 0.001341386 
    ## Run 169 stress 0.0004824213 
    ## ... Procrustes: rmse 0.01604965  max resid 0.02212898 
    ## Run 170 stress 9.102489e-05 
    ## ... Procrustes: rmse 0.0001781916  max resid 0.0002608748 
    ## ... Similar to previous best
    ## Run 171 stress 0.0001382801 
    ## ... Procrustes: rmse 0.008598132  max resid 0.01184731 
    ## Run 172 stress 0.001323451 
    ## Run 173 stress 8.83506e-05 
    ## ... Procrustes: rmse 0.0001088962  max resid 0.0001400125 
    ## ... Similar to previous best
    ## Run 174 stress 9.664966e-05 
    ## ... Procrustes: rmse 0.0001066263  max resid 0.0001674736 
    ## ... Similar to previous best
    ## Run 175 stress 9.844743e-05 
    ## ... Procrustes: rmse 0.0001625694  max resid 0.0002490435 
    ## ... Similar to previous best
    ## Run 176 stress 0.1491957 
    ## Run 177 stress 0.0007007811 
    ## Run 178 stress 7.72596e-05 
    ## ... Procrustes: rmse 0.0008243749  max resid 0.001136551 
    ## ... Similar to previous best
    ## Run 179 stress 9.74303e-05 
    ## ... Procrustes: rmse 0.0001222747  max resid 0.0002004235 
    ## ... Similar to previous best
    ## Run 180 stress 0.0004055908 
    ## ... Procrustes: rmse 0.01474859  max resid 0.02033436 
    ## Run 181 stress 0.0001086765 
    ## ... Procrustes: rmse 0.0126055  max resid 0.0173391 
    ## Run 182 stress 0.1990774 
    ## Run 183 stress 9.900955e-05 
    ## ... Procrustes: rmse 0.0001849045  max resid 0.0002721703 
    ## ... Similar to previous best
    ## Run 184 stress 0.001357652 
    ## Run 185 stress 0.0005716985 
    ## Run 186 stress 8.830679e-05 
    ## ... Procrustes: rmse 0.0001457302  max resid 0.0002227654 
    ## ... Similar to previous best
    ## Run 187 stress 0.001280014 
    ## Run 188 stress 0.0002531788 
    ## ... Procrustes: rmse 0.01163691  max resid 0.01604059 
    ## Run 189 stress 9.441646e-05 
    ## ... Procrustes: rmse 0.0001729976  max resid 0.0002513548 
    ## ... Similar to previous best
    ## Run 190 stress 9.243171e-05 
    ## ... Procrustes: rmse 0.0001586238  max resid 0.0002411232 
    ## ... Similar to previous best
    ## Run 191 stress 8.679013e-05 
    ## ... Procrustes: rmse 0.003962848  max resid 0.00547194 
    ## ... Similar to previous best
    ## Run 192 stress 0.001186556 
    ## Run 193 stress 0.001172894 
    ## Run 194 stress 0.001061791 
    ## Run 195 stress 0.0005124515 
    ## ... Procrustes: rmse 0.01657878  max resid 0.0228602 
    ## Run 196 stress 0.2696744 
    ## Run 197 stress 8.718951e-05 
    ## ... Procrustes: rmse 0.0002398223  max resid 0.0003440897 
    ## ... Similar to previous best
    ## Run 198 stress 0.0004785667 
    ## ... Procrustes: rmse 0.01602292  max resid 0.02209277 
    ## Run 199 stress 9.221392e-05 
    ## ... Procrustes: rmse 0.0002723298  max resid 0.0003619797 
    ## ... Similar to previous best
    ## Run 200 stress 0.2842805 
    ## Run 201 stress 7.288763e-05 
    ## ... Procrustes: rmse 0.0001131581  max resid 0.0001963942 
    ## ... Similar to previous best
    ## Run 202 stress 0.1491957 
    ## Run 203 stress 0.2528294 
    ## Run 204 stress 0.001423362 
    ## Run 205 stress 8.208649e-05 
    ## ... Procrustes: rmse 0.0002127954  max resid 0.0003155192 
    ## ... Similar to previous best
    ## Run 206 stress 0.1491957 
    ## Run 207 stress 9.685374e-05 
    ## ... Procrustes: rmse 0.000167458  max resid 0.0002431225 
    ## ... Similar to previous best
    ## Run 208 stress 0.001095389 
    ## Run 209 stress 9.082448e-05 
    ## ... Procrustes: rmse 0.0001482451  max resid 0.0002106236 
    ## ... Similar to previous best
    ## Run 210 stress 0.0007086961 
    ## Run 211 stress 0.0007252699 
    ## Run 212 stress 0.001275321 
    ## Run 213 stress 8.58923e-05 
    ## ... Procrustes: rmse 0.0001405151  max resid 0.0002139501 
    ## ... Similar to previous best
    ## Run 214 stress 9.320422e-05 
    ## ... Procrustes: rmse 0.0001836963  max resid 0.0002671133 
    ## ... Similar to previous best
    ## Run 215 stress 0.2842805 
    ## Run 216 stress 0.0004880058 
    ## ... Procrustes: rmse 0.016178  max resid 0.02230859 
    ## Run 217 stress 9.143483e-05 
    ## ... Procrustes: rmse 0.0001530151  max resid 0.0002325621 
    ## ... Similar to previous best
    ## Run 218 stress 0.3083098 
    ## Run 219 stress 0.001301267 
    ## Run 220 stress 0.1491957 
    ## Run 221 stress 0.0004667002 
    ## ... Procrustes: rmse 0.01575629  max resid 0.02172549 
    ## Run 222 stress 9.221493e-05 
    ## ... Procrustes: rmse 0.000146875  max resid 0.0002266223 
    ## ... Similar to previous best
    ## Run 223 stress 9.638175e-05 
    ## ... Procrustes: rmse 0.0001675245  max resid 0.0002532298 
    ## ... Similar to previous best
    ## Run 224 stress 9.707954e-05 
    ## ... Procrustes: rmse 0.000109477  max resid 0.0001846372 
    ## ... Similar to previous best
    ## Run 225 stress 0.2494706 
    ## Run 226 stress 7.957985e-05 
    ## ... Procrustes: rmse 0.0002864551  max resid 0.0004027272 
    ## ... Similar to previous best
    ## Run 227 stress 0.001223301 
    ## Run 228 stress 0.001341207 
    ## Run 229 stress 9.499888e-05 
    ## ... Procrustes: rmse 0.0001843106  max resid 0.0002470373 
    ## ... Similar to previous best
    ## Run 230 stress 0.001342305 
    ## Run 231 stress 9.638937e-05 
    ## ... Procrustes: rmse 0.0001884419  max resid 0.0002618133 
    ## ... Similar to previous best
    ## Run 232 stress 9.405861e-05 
    ## ... Procrustes: rmse 0.0001546504  max resid 0.000251447 
    ## ... Similar to previous best
    ## Run 233 stress 0.0004749772 
    ## ... Procrustes: rmse 0.01595285  max resid 0.02199664 
    ## Run 234 stress 0.0007433255 
    ## Run 235 stress 0.001192856 
    ## Run 236 stress 8.231062e-05 
    ## ... Procrustes: rmse 0.0001550153  max resid 0.0002079952 
    ## ... Similar to previous best
    ## Run 237 stress 0.001515933 
    ## Run 238 stress 0.1776017 
    ## Run 239 stress 9.693964e-05 
    ## ... Procrustes: rmse 0.0001822971  max resid 0.0002458443 
    ## ... Similar to previous best
    ## Run 240 stress 0.3083099 
    ## Run 241 stress 0.001478522 
    ## Run 242 stress 0.001428722 
    ## Run 243 stress 9.976753e-05 
    ## ... Procrustes: rmse 0.0001895721  max resid 0.0002546085 
    ## ... Similar to previous best
    ## Run 244 stress 9.58015e-05 
    ## ... Procrustes: rmse 0.0001801228  max resid 0.0002717532 
    ## ... Similar to previous best
    ## Run 245 stress 0.0004960835 
    ## ... Procrustes: rmse 0.02693301  max resid 0.03715882 
    ## Run 246 stress 0.0004677976 
    ## ... Procrustes: rmse 0.01582769  max resid 0.02182424 
    ## Run 247 stress 9.840958e-05 
    ## ... Procrustes: rmse 0.0001926824  max resid 0.0002533963 
    ## ... Similar to previous best
    ## Run 248 stress 9.736587e-05 
    ## ... Procrustes: rmse 0.0001843868  max resid 0.0002732192 
    ## ... Similar to previous best
    ## Run 249 stress 0.001315367 
    ## Run 250 stress 0.0004722165 
    ## ... Procrustes: rmse 0.01591429  max resid 0.02194325 
    ## Run 251 stress 9.379352e-05 
    ## ... Procrustes: rmse 0.001015192  max resid 0.001360057 
    ## ... Similar to previous best
    ## Run 252 stress 0.2373687 
    ## Run 253 stress 0.2842805 
    ## Run 254 stress 0.0008451178 
    ## Run 255 stress 0.0004761494 
    ## ... Procrustes: rmse 0.01598203  max resid 0.02203631 
    ## Run 256 stress 0.0008536464 
    ## Run 257 stress 0.1491957 
    ## Run 258 stress 9.356406e-05 
    ## ... Procrustes: rmse 0.0001800646  max resid 0.0002406919 
    ## ... Similar to previous best
    ## Run 259 stress 0.001224171 
    ## Run 260 stress 0.0006592655 
    ## Run 261 stress 9.972549e-05 
    ## ... Procrustes: rmse 0.0001838878  max resid 0.00026998 
    ## ... Similar to previous best
    ## Run 262 stress 0.001422189 
    ## Run 263 stress 7.944996e-05 
    ## ... Procrustes: rmse 0.0001008094  max resid 0.0001809249 
    ## ... Similar to previous best
    ## Run 264 stress 0.001291334 
    ## Run 265 stress 0.001261122 
    ## Run 266 stress 0.2520602 
    ## Run 267 stress 9.490316e-05 
    ## ... Procrustes: rmse 0.001430726  max resid 0.001973101 
    ## ... Similar to previous best
    ## Run 268 stress 0.1990774 
    ## Run 269 stress 0.001213787 
    ## Run 270 stress 9.963025e-05 
    ## ... Procrustes: rmse 0.0001832876  max resid 0.0002685201 
    ## ... Similar to previous best
    ## Run 271 stress 0.2842805 
    ## Run 272 stress 0.001343197 
    ## Run 273 stress 0.2848517 
    ## Run 274 stress 0.177601 
    ## Run 275 stress 9.592097e-05 
    ## ... Procrustes: rmse 0.0001764589  max resid 0.0002551367 
    ## ... Similar to previous best
    ## Run 276 stress 0.1491957 
    ## Run 277 stress 8.57519e-05 
    ## ... Procrustes: rmse 0.0001083949  max resid 0.0001524613 
    ## ... Similar to previous best
    ## Run 278 stress 9.125924e-05 
    ## ... Procrustes: rmse 0.0001748716  max resid 0.0002542546 
    ## ... Similar to previous best
    ## Run 279 stress 0.0004959265 
    ## ... Procrustes: rmse 0.01630109  max resid 0.02247636 
    ## Run 280 stress 9.065515e-05 
    ## ... Procrustes: rmse 0.0001097858  max resid 0.000144856 
    ## ... Similar to previous best
    ## Run 281 stress 9.552746e-05 
    ## ... Procrustes: rmse 0.0001858193  max resid 0.0002495856 
    ## ... Similar to previous best
    ## Run 282 stress 0.001365749 
    ## Run 283 stress 0.0001404189 
    ## ... Procrustes: rmse 0.01433009  max resid 0.0197223 
    ## Run 284 stress 9.876082e-05 
    ## ... Procrustes: rmse 0.000193345  max resid 0.0002693919 
    ## ... Similar to previous best
    ## Run 285 stress 0.177601 
    ## Run 286 stress 0.1491957 
    ## Run 287 stress 9.821116e-05 
    ## ... Procrustes: rmse 0.0001879211  max resid 0.0002763187 
    ## ... Similar to previous best
    ## Run 288 stress 0.001154451 
    ## Run 289 stress 7.867395e-05 
    ## ... Procrustes: rmse 0.0001372143  max resid 0.0002159007 
    ## ... Similar to previous best
    ## Run 290 stress 9.32215e-05 
    ## ... Procrustes: rmse 0.0001432631  max resid 0.0002526609 
    ## ... Similar to previous best
    ## Run 291 stress 9.953439e-05 
    ## ... Procrustes: rmse 0.0001469152  max resid 0.0002219789 
    ## ... Similar to previous best
    ## Run 292 stress 0.0004218018 
    ## ... Procrustes: rmse 0.01501511  max resid 0.02071213 
    ## Run 293 stress 0.001023512 
    ## Run 294 stress 0.0009602313 
    ## Run 295 stress 9.783978e-05 
    ## ... Procrustes: rmse 0.0001718113  max resid 0.0002606394 
    ## ... Similar to previous best
    ## Run 296 stress 0.0001287727 
    ## ... Procrustes: rmse 0.008297325  max resid 0.0114321 
    ## Run 297 stress 9.856314e-05 
    ## ... Procrustes: rmse 0.0001969297  max resid 0.0002794727 
    ## ... Similar to previous best
    ## Run 298 stress 0.001050904 
    ## Run 299 stress 0.001309556 
    ## Run 300 stress 8.583871e-05 
    ## ... Procrustes: rmse 0.0001067498  max resid 0.0001817685 
    ## ... Similar to previous best
    ## Run 301 stress 0.0006974875 
    ## Run 302 stress 0.001412254 
    ## Run 303 stress 0.1776012 
    ## Run 304 stress 0.001379462 
    ## Run 305 stress 0.0001508883 
    ## ... Procrustes: rmse 0.01485444  max resid 0.02044684 
    ## Run 306 stress 0.1491957 
    ## Run 307 stress 0.1990774 
    ## Run 308 stress 0.2528294 
    ## Run 309 stress 7.662872e-05 
    ## ... Procrustes: rmse 0.001566714  max resid 0.002173106 
    ## ... Similar to previous best
    ## Run 310 stress 9.446672e-05 
    ## ... Procrustes: rmse 0.0001726387  max resid 0.0002513015 
    ## ... Similar to previous best
    ## Run 311 stress 0.001148508 
    ## Run 312 stress 0.1491957 
    ## Run 313 stress 0.0004953347 
    ## ... Procrustes: rmse 0.01612413  max resid 0.02225816 
    ## Run 314 stress 0.001137389 
    ## Run 315 stress 0.001280926 
    ## Run 316 stress 9.747652e-05 
    ## ... Procrustes: rmse 0.0001587441  max resid 0.000247429 
    ## ... Similar to previous best
    ## Run 317 stress 9.579191e-05 
    ## ... Procrustes: rmse 0.0001866887  max resid 0.0002582075 
    ## ... Similar to previous best
    ## Run 318 stress 0.001227528 
    ## Run 319 stress 0.001329034 
    ## Run 320 stress 0.0004587541 
    ## ... Procrustes: rmse 0.01568616  max resid 0.02162814 
    ## Run 321 stress 0.001384448 
    ## Run 322 stress 0.1990774 
    ## Run 323 stress 0.001419015 
    ## Run 324 stress 9.83667e-05 
    ## ... Procrustes: rmse 0.0001796598  max resid 0.0002578636 
    ## ... Similar to previous best
    ## Run 325 stress 9.461751e-05 
    ## ... Procrustes: rmse 0.0001615243  max resid 0.0002586896 
    ## ... Similar to previous best
    ## Run 326 stress 0.001186719 
    ## Run 327 stress 0.3078457 
    ## Run 328 stress 8.401871e-05 
    ## ... Procrustes: rmse 0.0001152762  max resid 0.0001908309 
    ## ... Similar to previous best
    ## Run 329 stress 0.1491957 
    ## Run 330 stress 0.001348726 
    ## Run 331 stress 9.763303e-05 
    ## ... Procrustes: rmse 0.0001310332  max resid 0.0001778264 
    ## ... Similar to previous best
    ## Run 332 stress 9.467583e-05 
    ## ... Procrustes: rmse 0.0001598755  max resid 0.0002451379 
    ## ... Similar to previous best
    ## Run 333 stress 0.1990774 
    ## Run 334 stress 0.001210291 
    ## Run 335 stress 9.917749e-05 
    ## ... Procrustes: rmse 0.0001632246  max resid 0.0002512964 
    ## ... Similar to previous best
    ## Run 336 stress 0.1990774 
    ## Run 337 stress 9.720972e-05 
    ## ... Procrustes: rmse 0.0001939016  max resid 0.0002754243 
    ## ... Similar to previous best
    ## Run 338 stress 0.1491957 
    ## Run 339 stress 0.001451406 
    ## Run 340 stress 8.607938e-05 
    ## ... Procrustes: rmse 0.0001506222  max resid 0.0002464854 
    ## ... Similar to previous best
    ## Run 341 stress 0.001315496 
    ## Run 342 stress 0.2842598 
    ## Run 343 stress 9.692631e-05 
    ## ... Procrustes: rmse 0.002028075  max resid 0.002788623 
    ## ... Similar to previous best
    ## Run 344 stress 0.0003201685 
    ## ... Procrustes: rmse 0.02163331  max resid 0.02982038 
    ## Run 345 stress 0.1491957 
    ## Run 346 stress 9.97252e-05 
    ## ... Procrustes: rmse 0.0001408296  max resid 0.0002220951 
    ## ... Similar to previous best
    ## Run 347 stress 0.0004788439 
    ## ... Procrustes: rmse 0.01602712  max resid 0.02209874 
    ## Run 348 stress 0.001069476 
    ## Run 349 stress 9.526312e-05 
    ## ... Procrustes: rmse 0.0001351252  max resid 0.0001907533 
    ## ... Similar to previous best
    ## Run 350 stress 0.1491957 
    ## Run 351 stress 0.001329087 
    ## Run 352 stress 0.000488688 
    ## ... Procrustes: rmse 0.026731  max resid 0.03687893 
    ## Run 353 stress 0.001396423 
    ## Run 354 stress 8.741753e-05 
    ## ... Procrustes: rmse 0.0001547116  max resid 0.0002158763 
    ## ... Similar to previous best
    ## Run 355 stress 0.001273503 
    ## Run 356 stress 8.991518e-05 
    ## ... Procrustes: rmse 0.0001524093  max resid 0.0002314459 
    ## ... Similar to previous best
    ## Run 357 stress 8.994125e-05 
    ## ... Procrustes: rmse 0.007624829  max resid 0.01049303 
    ## Run 358 stress 0.001162533 
    ## Run 359 stress 9.748058e-05 
    ## ... Procrustes: rmse 0.0001663124  max resid 0.0002474436 
    ## ... Similar to previous best
    ## Run 360 stress 9.956073e-05 
    ## ... Procrustes: rmse 0.006231297  max resid 0.008598396 
    ## Run 361 stress 0.1491957 
    ## Run 362 stress 0.0008377162 
    ## Run 363 stress 0.001249442 
    ## Run 364 stress 0.001226867 
    ## Run 365 stress 0.0013808 
    ## Run 366 stress 0.2570108 
    ## Run 367 stress 0.001212485 
    ## Run 368 stress 0.0003742912 
    ## ... Procrustes: rmse 0.01373179  max resid 0.0189324 
    ## Run 369 stress 9.789094e-05 
    ## ... Procrustes: rmse 0.0001961805  max resid 0.0002798237 
    ## ... Similar to previous best
    ## Run 370 stress 9.846111e-05 
    ## ... Procrustes: rmse 0.0001646157  max resid 0.0002707054 
    ## ... Similar to previous best
    ## Run 371 stress 8.446965e-05 
    ## ... Procrustes: rmse 0.0001369292  max resid 0.0002075714 
    ## ... Similar to previous best
    ## Run 372 stress 0.001372761 
    ## Run 373 stress 9.234325e-05 
    ## ... Procrustes: rmse 0.0001537785  max resid 0.0002357096 
    ## ... Similar to previous best
    ## Run 374 stress 9.561294e-05 
    ## ... Procrustes: rmse 0.0001794119  max resid 0.0002428567 
    ## ... Similar to previous best
    ## Run 375 stress 9.149965e-05 
    ## ... Procrustes: rmse 0.0001266379  max resid 0.0001764232 
    ## ... Similar to previous best
    ## Run 376 stress 0.0008614516 
    ## Run 377 stress 0.00116986 
    ## Run 378 stress 9.676155e-05 
    ## ... Procrustes: rmse 0.0001855447  max resid 0.0002648787 
    ## ... Similar to previous best
    ## Run 379 stress 0.0007080939 
    ## Run 380 stress 9.850351e-05 
    ## ... Procrustes: rmse 0.0001333217  max resid 0.0002200018 
    ## ... Similar to previous best
    ## Run 381 stress 9.096452e-05 
    ## ... Procrustes: rmse 0.0001732134  max resid 0.0002294734 
    ## ... Similar to previous best
    ## Run 382 stress 9.954982e-05 
    ## ... Procrustes: rmse 0.0003645718  max resid 0.0005152986 
    ## ... Similar to previous best
    ## Run 383 stress 0.0003386613 
    ## ... Procrustes: rmse 0.02225183  max resid 0.03067614 
    ## Run 384 stress 0.000505806 
    ## ... Procrustes: rmse 0.01647068  max resid 0.0227105 
    ## Run 385 stress 9.084241e-05 
    ## ... Procrustes: rmse 0.0001654756  max resid 0.000239741 
    ## ... Similar to previous best
    ## Run 386 stress 9.799882e-05 
    ## ... Procrustes: rmse 0.0001902501  max resid 0.0002518169 
    ## ... Similar to previous best
    ## Run 387 stress 0.0002872317 
    ## ... Procrustes: rmse 0.01240712  max resid 0.01710394 
    ## Run 388 stress 0.2528294 
    ## Run 389 stress 0.001319558 
    ## Run 390 stress 9.352185e-05 
    ## ... Procrustes: rmse 0.0001217971  max resid 0.0001974064 
    ## ... Similar to previous best
    ## Run 391 stress 0.001335177 
    ## Run 392 stress 8.30299e-05 
    ## ... Procrustes: rmse 0.0001028844  max resid 0.0001482573 
    ## ... Similar to previous best
    ## Run 393 stress 0.1776017 
    ## Run 394 stress 9.432402e-05 
    ## ... Procrustes: rmse 0.0001689717  max resid 0.0002433197 
    ## ... Similar to previous best
    ## Run 395 stress 8.659746e-05 
    ## ... Procrustes: rmse 0.0001712308  max resid 0.0002503361 
    ## ... Similar to previous best
    ## Run 396 stress 0.001193078 
    ## Run 397 stress 8.289359e-05 
    ## ... Procrustes: rmse 0.0001140937  max resid 0.0001720684 
    ## ... Similar to previous best
    ## Run 398 stress 0.1491957 
    ## Run 399 stress 0.0004588797 
    ## ... Procrustes: rmse 0.01568593  max resid 0.02162786 
    ## Run 400 stress 9.633355e-05 
    ## ... Procrustes: rmse 0.0001840631  max resid 0.0002437031 
    ## ... Similar to previous best
    ## Run 401 stress 0.001133838 
    ## Run 402 stress 0.1491957 
    ## Run 403 stress 8.783859e-05 
    ## ... Procrustes: rmse 0.000145908  max resid 0.0001971131 
    ## ... Similar to previous best
    ## Run 404 stress 0.001222635 
    ## Run 405 stress 0.001199983 
    ## Run 406 stress 0.0004673417 
    ## ... Procrustes: rmse 0.01583103  max resid 0.02183231 
    ## Run 407 stress 0.0001146686 
    ## ... Procrustes: rmse 0.007827996  max resid 0.01078443 
    ## Run 408 stress 0.257011 
    ## Run 409 stress 0.001329542 
    ## Run 410 stress 0.1990774 
    ## Run 411 stress 0.001157197 
    ## Run 412 stress 8.703943e-05 
    ## ... Procrustes: rmse 0.00015915  max resid 0.0002246091 
    ## ... Similar to previous best
    ## Run 413 stress 0.1491957 
    ## Run 414 stress 0.001466126 
    ## Run 415 stress 0.1491957 
    ## Run 416 stress 0.1990774 
    ## Run 417 stress 0.0009853795 
    ## Run 418 stress 0.0004547233 
    ## ... Procrustes: rmse 0.01561764  max resid 0.02153474 
    ## Run 419 stress 8.52097e-05 
    ## ... Procrustes: rmse 0.0001558469  max resid 0.0002491388 
    ## ... Similar to previous best
    ## Run 420 stress 0.1491957 
    ## Run 421 stress 0.001301379 
    ## Run 422 stress 0.000445333 
    ## ... Procrustes: rmse 0.01545495  max resid 0.02130894 
    ## Run 423 stress 0.001244153 
    ## Run 424 stress 9.806088e-05 
    ## ... Procrustes: rmse 0.0001888373  max resid 0.0002504999 
    ## ... Similar to previous best
    ## Run 425 stress 8.868508e-05 
    ## ... Procrustes: rmse 0.0001575432  max resid 0.0002241007 
    ## ... Similar to previous best
    ## Run 426 stress 0.001302195 
    ## Run 427 stress 7.654464e-05 
    ## ... Procrustes: rmse 0.0001319577  max resid 0.0001835546 
    ## ... Similar to previous best
    ## Run 428 stress 0.001347578 
    ## Run 429 stress 0.0009856261 
    ## Run 430 stress 0.00144269 
    ## Run 431 stress 0.001128011 
    ## Run 432 stress 0.001243104 
    ## Run 433 stress 9.132315e-05 
    ## ... Procrustes: rmse 0.007814233  max resid 0.01068541 
    ## Run 434 stress 9.380979e-05 
    ## ... Procrustes: rmse 0.0001618273  max resid 0.0002611849 
    ## ... Similar to previous best
    ## Run 435 stress 0.0004717093 
    ## ... Procrustes: rmse 0.01590577  max resid 0.0219312 
    ## Run 436 stress 9.676293e-05 
    ## ... Procrustes: rmse 0.0001361189  max resid 0.0001772978 
    ## ... Similar to previous best
    ## Run 437 stress 0.1990774 
    ## Run 438 stress 0.001431827 
    ## Run 439 stress 0.001171392 
    ## Run 440 stress 0.2440272 
    ## Run 441 stress 0.0004645623 
    ## ... Procrustes: rmse 0.01572699  max resid 0.02168949 
    ## Run 442 stress 0.001053629 
    ## Run 443 stress 0.001370219 
    ## Run 444 stress 0.0004505958 
    ## ... Procrustes: rmse 0.01554225  max resid 0.02143016 
    ## Run 445 stress 0.2494706 
    ## Run 446 stress 9.944657e-05 
    ## ... Procrustes: rmse 0.0001956608  max resid 0.0002811879 
    ## ... Similar to previous best
    ## Run 447 stress 8.799994e-05 
    ## ... Procrustes: rmse 0.001938834  max resid 0.002681386 
    ## ... Similar to previous best
    ## Run 448 stress 8.598856e-05 
    ## ... Procrustes: rmse 0.000161036  max resid 0.0002113588 
    ## ... Similar to previous best
    ## Run 449 stress 0.1491957 
    ## Run 450 stress 0.2848533 
    ## Run 451 stress 7.841433e-05 
    ## ... Procrustes: rmse 0.0009097815  max resid 0.001256018 
    ## ... Similar to previous best
    ## Run 452 stress 9.412664e-05 
    ## ... Procrustes: rmse 0.0001283051  max resid 0.0002103865 
    ## ... Similar to previous best
    ## Run 453 stress 0.001285485 
    ## Run 454 stress 8.415953e-05 
    ## ... Procrustes: rmse 9.577651e-05  max resid 0.0001756973 
    ## ... Similar to previous best
    ## Run 455 stress 9.354784e-05 
    ## ... Procrustes: rmse 0.003074929  max resid 0.004292791 
    ## ... Similar to previous best
    ## Run 456 stress 0.001185495 
    ## Run 457 stress 9.516336e-05 
    ## ... Procrustes: rmse 0.0001681486  max resid 0.0002562034 
    ## ... Similar to previous best
    ## Run 458 stress 8.763721e-05 
    ## ... Procrustes: rmse 0.0001515904  max resid 0.0002147925 
    ## ... Similar to previous best
    ## Run 459 stress 0.177601 
    ## Run 460 stress 0.001216609 
    ## Run 461 stress 0.2612458 
    ## Run 462 stress 0.0007422376 
    ## Run 463 stress 7.279707e-05 
    ## ... Procrustes: rmse 0.0001247492  max resid 0.0001702573 
    ## ... Similar to previous best
    ## Run 464 stress 0.1491957 
    ## Run 465 stress 0.2842805 
    ## Run 466 stress 9.072009e-05 
    ## ... Procrustes: rmse 0.0001651512  max resid 0.0002206956 
    ## ... Similar to previous best
    ## Run 467 stress 0.001107864 
    ## Run 468 stress 8.338216e-05 
    ## ... Procrustes: rmse 0.002960644  max resid 0.004060362 
    ## ... Similar to previous best
    ## Run 469 stress 0.001371228 
    ## Run 470 stress 9.036667e-05 
    ## ... Procrustes: rmse 0.0001157091  max resid 0.0001875975 
    ## ... Similar to previous best
    ## Run 471 stress 0.001189058 
    ## Run 472 stress 0.2570113 
    ## Run 473 stress 0.001249021 
    ## Run 474 stress 0.1990774 
    ## Run 475 stress 0.001172674 
    ## Run 476 stress 9.284612e-05 
    ## ... Procrustes: rmse 0.0001590668  max resid 0.0002381009 
    ## ... Similar to previous best
    ## Run 477 stress 9.035571e-05 
    ## ... Procrustes: rmse 0.0005732302  max resid 0.0008015653 
    ## ... Similar to previous best
    ## Run 478 stress 9.220296e-05 
    ## ... Procrustes: rmse 0.0009473269  max resid 0.001304991 
    ## ... Similar to previous best
    ## Run 479 stress 0.001324034 
    ## Run 480 stress 0.001193568 
    ## Run 481 stress 0.1990774 
    ## Run 482 stress 9.757447e-05 
    ## ... Procrustes: rmse 0.000194264  max resid 0.0002819581 
    ## ... Similar to previous best
    ## Run 483 stress 0.1776012 
    ## Run 484 stress 0.177601 
    ## Run 485 stress 5.829907e-05 
    ## ... Procrustes: rmse 0.0003950919  max resid 0.0005404687 
    ## ... Similar to previous best
    ## Run 486 stress 0.1776012 
    ## Run 487 stress 9.671879e-05 
    ## ... Procrustes: rmse 0.000159938  max resid 0.0002515901 
    ## ... Similar to previous best
    ## Run 488 stress 9.614396e-05 
    ## ... Procrustes: rmse 0.0001829525  max resid 0.0002511895 
    ## ... Similar to previous best
    ## Run 489 stress 0.2842805 
    ## Run 490 stress 0.001291947 
    ## Run 491 stress 0.1491957 
    ## Run 492 stress 0.1990774 
    ## Run 493 stress 0.0004644148 
    ## ... Procrustes: rmse 0.01575435  max resid 0.02172493 
    ## Run 494 stress 9.846575e-05 
    ## ... Procrustes: rmse 0.0001958898  max resid 0.0002834104 
    ## ... Similar to previous best
    ## Run 495 stress 9.916406e-05 
    ## ... Procrustes: rmse 0.008105611  max resid 0.01116007 
    ## Run 496 stress 9.946675e-05 
    ## ... Procrustes: rmse 0.0001880225  max resid 0.0002725407 
    ## ... Similar to previous best
    ## Run 497 stress 0.001337302 
    ## Run 498 stress 9.423271e-05 
    ## ... Procrustes: rmse 0.0001720151  max resid 0.0002488259 
    ## ... Similar to previous best
    ## Run 499 stress 0.001412459 
    ## Run 500 stress 9.914512e-05 
    ## ... Procrustes: rmse 0.0001688964  max resid 0.0002595224 
    ## ... Similar to previous best
    ## *** Best solution repeated 139 times

    ## Warning in metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
SD_beta_geo_NMDS <- metaMDS(SD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09021173 
    ## Run 1 stress 0.1065375 
    ## Run 2 stress 0.08926088 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04369267  max resid 0.1170552 
    ## Run 3 stress 0.1101446 
    ## Run 4 stress 0.1071321 
    ## Run 5 stress 0.08938974 
    ## ... Procrustes: rmse 0.03596652  max resid 0.1183055 
    ## Run 6 stress 0.09503441 
    ## Run 7 stress 0.08938538 
    ## ... Procrustes: rmse 0.01176607  max resid 0.0397294 
    ## Run 8 stress 0.1092093 
    ## Run 9 stress 0.08946334 
    ## ... Procrustes: rmse 0.03759537  max resid 0.1178957 
    ## Run 10 stress 0.09087341 
    ## Run 11 stress 0.08926112 
    ## ... Procrustes: rmse 0.0001243549  max resid 0.0003372689 
    ## ... Similar to previous best
    ## Run 12 stress 0.08946652 
    ## ... Procrustes: rmse 0.03326447  max resid 0.1164009 
    ## Run 13 stress 0.0892609 
    ## ... Procrustes: rmse 2.540788e-05  max resid 6.168535e-05 
    ## ... Similar to previous best
    ## Run 14 stress 0.08938537 
    ## ... Procrustes: rmse 0.0117884  max resid 0.03974529 
    ## Run 15 stress 0.1061304 
    ## Run 16 stress 0.08926106 
    ## ... Procrustes: rmse 9.477263e-05  max resid 0.0003183089 
    ## ... Similar to previous best
    ## Run 17 stress 0.09099578 
    ## Run 18 stress 0.09088601 
    ## Run 19 stress 0.0950342 
    ## Run 20 stress 0.110145 
    ## Run 21 stress 0.0950345 
    ## Run 22 stress 0.1074434 
    ## Run 23 stress 0.08946338 
    ## ... Procrustes: rmse 0.03753063  max resid 0.1177963 
    ## Run 24 stress 0.1079015 
    ## Run 25 stress 0.1075421 
    ## Run 26 stress 0.1091829 
    ## Run 27 stress 0.109213 
    ## Run 28 stress 0.1092559 
    ## Run 29 stress 0.1087583 
    ## Run 30 stress 0.1087548 
    ## Run 31 stress 0.1056903 
    ## Run 32 stress 0.1091239 
    ## Run 33 stress 0.08938541 
    ## ... Procrustes: rmse 0.0118134  max resid 0.03971795 
    ## Run 34 stress 0.08926097 
    ## ... Procrustes: rmse 5.266639e-05  max resid 0.0001745832 
    ## ... Similar to previous best
    ## Run 35 stress 0.09228493 
    ## Run 36 stress 0.08946339 
    ## ... Procrustes: rmse 0.03752639  max resid 0.1177874 
    ## Run 37 stress 0.08938541 
    ## ... Procrustes: rmse 0.01175702  max resid 0.03975768 
    ## Run 38 stress 0.1075421 
    ## Run 39 stress 0.08938963 
    ## ... Procrustes: rmse 0.03599057  max resid 0.1183433 
    ## Run 40 stress 0.09109116 
    ## Run 41 stress 0.1067503 
    ## Run 42 stress 0.08926099 
    ## ... Procrustes: rmse 5.594497e-05  max resid 0.0001814376 
    ## ... Similar to previous best
    ## Run 43 stress 0.09262381 
    ## Run 44 stress 0.1101454 
    ## Run 45 stress 0.09503442 
    ## Run 46 stress 0.0893855 
    ## ... Procrustes: rmse 0.01184042  max resid 0.03969907 
    ## Run 47 stress 0.08926071 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001423592  max resid 0.000469945 
    ## ... Similar to previous best
    ## Run 48 stress 0.107108 
    ## Run 49 stress 0.09088608 
    ## Run 50 stress 0.08938538 
    ## ... Procrustes: rmse 0.01179681  max resid 0.03972055 
    ## Run 51 stress 0.08938537 
    ## ... Procrustes: rmse 0.01181289  max resid 0.03973031 
    ## Run 52 stress 0.08926076 
    ## ... Procrustes: rmse 6.143416e-05  max resid 0.0001766378 
    ## ... Similar to previous best
    ## Run 53 stress 0.08946652 
    ## ... Procrustes: rmse 0.03326774  max resid 0.1164287 
    ## Run 54 stress 0.1092127 
    ## Run 55 stress 0.08946651 
    ## ... Procrustes: rmse 0.03326071  max resid 0.1164097 
    ## Run 56 stress 0.0892607 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002891899  max resid 0.0007167418 
    ## ... Similar to previous best
    ## Run 57 stress 0.09039132 
    ## Run 58 stress 0.1080638 
    ## Run 59 stress 0.09612805 
    ## Run 60 stress 0.1074434 
    ## Run 61 stress 0.1092559 
    ## Run 62 stress 0.0950345 
    ## Run 63 stress 0.1066195 
    ## Run 64 stress 0.08951716 
    ## ... Procrustes: rmse 0.0350062  max resid 0.1156143 
    ## Run 65 stress 0.08938964 
    ## ... Procrustes: rmse 0.03590231  max resid 0.1182356 
    ## Run 66 stress 0.09044609 
    ## Run 67 stress 0.08938971 
    ## ... Procrustes: rmse 0.03595116  max resid 0.1183156 
    ## Run 68 stress 0.09592149 
    ## Run 69 stress 0.08938538 
    ## ... Procrustes: rmse 0.01179749  max resid 0.0394932 
    ## Run 70 stress 0.09228498 
    ## Run 71 stress 0.0959433 
    ## Run 72 stress 0.1071322 
    ## Run 73 stress 0.1067589 
    ## Run 74 stress 0.1074433 
    ## Run 75 stress 0.08926072 
    ## ... Procrustes: rmse 6.075315e-05  max resid 0.0001863616 
    ## ... Similar to previous best
    ## Run 76 stress 0.09503443 
    ## Run 77 stress 0.08946656 
    ## ... Procrustes: rmse 0.0331562  max resid 0.1162753 
    ## Run 78 stress 0.09503434 
    ## Run 79 stress 0.09503436 
    ## Run 80 stress 0.1113255 
    ## Run 81 stress 0.0892607 
    ## ... Procrustes: rmse 0.0002743143  max resid 0.0006592852 
    ## ... Similar to previous best
    ## Run 82 stress 0.09099546 
    ## Run 83 stress 0.1076302 
    ## Run 84 stress 0.09018708 
    ## Run 85 stress 0.108839 
    ## Run 86 stress 0.092285 
    ## Run 87 stress 0.08926093 
    ## ... Procrustes: rmse 0.0003391906  max resid 0.00107403 
    ## ... Similar to previous best
    ## Run 88 stress 0.08938538 
    ## ... Procrustes: rmse 0.0117789  max resid 0.03946771 
    ## Run 89 stress 0.1074432 
    ## Run 90 stress 0.1067389 
    ## Run 91 stress 0.1086561 
    ## Run 92 stress 0.0894633 
    ## ... Procrustes: rmse 0.03745372  max resid 0.1177175 
    ## Run 93 stress 0.09099569 
    ## Run 94 stress 0.1092131 
    ## Run 95 stress 0.08946654 
    ## ... Procrustes: rmse 0.03316185  max resid 0.116281 
    ## Run 96 stress 0.1076303 
    ## Run 97 stress 0.09130085 
    ## Run 98 stress 0.08951723 
    ## ... Procrustes: rmse 0.03498012  max resid 0.11558 
    ## Run 99 stress 0.09503455 
    ## Run 100 stress 0.08926108 
    ## ... Procrustes: rmse 0.0005317937  max resid 0.001508867 
    ## ... Similar to previous best
    ## Run 101 stress 0.0950343 
    ## Run 102 stress 0.08938542 
    ## ... Procrustes: rmse 0.01179802  max resid 0.03944465 
    ## Run 103 stress 0.1066199 
    ## Run 104 stress 0.08938966 
    ## ... Procrustes: rmse 0.03589343  max resid 0.11822 
    ## Run 105 stress 0.08938566 
    ## ... Procrustes: rmse 0.01184265  max resid 0.03938704 
    ## Run 106 stress 0.1075421 
    ## Run 107 stress 0.1056908 
    ## Run 108 stress 0.08926079 
    ## ... Procrustes: rmse 0.000301661  max resid 0.0009297185 
    ## ... Similar to previous best
    ## Run 109 stress 0.08946334 
    ## ... Procrustes: rmse 0.0374408  max resid 0.1176955 
    ## Run 110 stress 0.1101454 
    ## Run 111 stress 0.08946656 
    ## ... Procrustes: rmse 0.03315549  max resid 0.1162741 
    ## Run 112 stress 0.1052647 
    ## Run 113 stress 0.1063194 
    ## Run 114 stress 0.0894665 
    ## ... Procrustes: rmse 0.03318432  max resid 0.1163228 
    ## Run 115 stress 0.1074432 
    ## Run 116 stress 0.1101451 
    ## Run 117 stress 0.1067373 
    ## Run 118 stress 0.08946653 
    ## ... Procrustes: rmse 0.03319699  max resid 0.1163619 
    ## Run 119 stress 0.1092947 
    ## Run 120 stress 0.08946328 
    ## ... Procrustes: rmse 0.03746591  max resid 0.1177365 
    ## Run 121 stress 0.09088616 
    ## Run 122 stress 0.1103705 
    ## Run 123 stress 0.08926086 
    ## ... Procrustes: rmse 0.0002488807  max resid 0.0008216613 
    ## ... Similar to previous best
    ## Run 124 stress 0.1105339 
    ## Run 125 stress 0.09018711 
    ## Run 126 stress 0.0892607 
    ## ... New best solution
    ## ... Procrustes: rmse 1.149269e-05  max resid 2.854617e-05 
    ## ... Similar to previous best
    ## Run 127 stress 0.09178281 
    ## Run 128 stress 0.08938962 
    ## ... Procrustes: rmse 0.035909  max resid 0.1182449 
    ## Run 129 stress 0.08946655 
    ## ... Procrustes: rmse 0.03315732  max resid 0.1162779 
    ## Run 130 stress 0.09099539 
    ## Run 131 stress 0.09262339 
    ## Run 132 stress 0.1056892 
    ## Run 133 stress 0.08938542 
    ## ... Procrustes: rmse 0.011792  max resid 0.03951656 
    ## Run 134 stress 0.09178295 
    ## Run 135 stress 0.0892611 
    ## ... Procrustes: rmse 0.0005327705  max resid 0.001578593 
    ## ... Similar to previous best
    ## Run 136 stress 0.1087868 
    ## Run 137 stress 0.1063189 
    ## Run 138 stress 0.08938541 
    ## ... Procrustes: rmse 0.01179451  max resid 0.03944915 
    ## Run 139 stress 0.08938559 
    ## ... Procrustes: rmse 0.01169945  max resid 0.03952251 
    ## Run 140 stress 0.0894665 
    ## ... Procrustes: rmse 0.03317732  max resid 0.1163118 
    ## Run 141 stress 0.09021167 
    ## Run 142 stress 0.1087578 
    ## Run 143 stress 0.08926092 
    ## ... Procrustes: rmse 0.0004304308  max resid 0.001246242 
    ## ... Similar to previous best
    ## Run 144 stress 0.09130089 
    ## Run 145 stress 0.08926077 
    ## ... Procrustes: rmse 0.0003342156  max resid 0.0008825609 
    ## ... Similar to previous best
    ## Run 146 stress 0.1091234 
    ## Run 147 stress 0.08926067 
    ## ... New best solution
    ## ... Procrustes: rmse 9.204279e-05  max resid 0.0002678846 
    ## ... Similar to previous best
    ## Run 148 stress 0.1074431 
    ## Run 149 stress 0.1085127 
    ## Run 150 stress 0.09592133 
    ## Run 151 stress 0.0908734 
    ## Run 152 stress 0.08926076 
    ## ... Procrustes: rmse 0.0002435216  max resid 0.00059482 
    ## ... Similar to previous best
    ## Run 153 stress 0.09592164 
    ## Run 154 stress 0.08938549 
    ## ... Procrustes: rmse 0.01170177  max resid 0.03945627 
    ## Run 155 stress 0.1085126 
    ## Run 156 stress 0.08946668 
    ## ... Procrustes: rmse 0.03315896  max resid 0.1162697 
    ## Run 157 stress 0.09503454 
    ## Run 158 stress 0.09099537 
    ## Run 159 stress 0.1104472 
    ## Run 160 stress 0.1056903 
    ## Run 161 stress 0.09099568 
    ## Run 162 stress 0.1075421 
    ## Run 163 stress 0.08926078 
    ## ... Procrustes: rmse 0.0002587886  max resid 0.0005757462 
    ## ... Similar to previous best
    ## Run 164 stress 0.1071323 
    ## Run 165 stress 0.08926084 
    ## ... Procrustes: rmse 0.0003022067  max resid 0.0008010168 
    ## ... Similar to previous best
    ## Run 166 stress 0.08926069 
    ## ... Procrustes: rmse 0.0001189428  max resid 0.0003464405 
    ## ... Similar to previous best
    ## Run 167 stress 0.1052646 
    ## Run 168 stress 0.1074433 
    ## Run 169 stress 0.08938539 
    ## ... Procrustes: rmse 0.01178251  max resid 0.03955424 
    ## Run 170 stress 0.08951719 
    ## ... Procrustes: rmse 0.03505203  max resid 0.1156589 
    ## Run 171 stress 0.0892609 
    ## ... Procrustes: rmse 0.0003408022  max resid 0.0009051312 
    ## ... Similar to previous best
    ## Run 172 stress 0.08938539 
    ## ... Procrustes: rmse 0.01179315  max resid 0.03951044 
    ## Run 173 stress 0.08946651 
    ## ... Procrustes: rmse 0.03319574  max resid 0.1163316 
    ## Run 174 stress 0.09503441 
    ## Run 175 stress 0.1067594 
    ## Run 176 stress 0.09262318 
    ## Run 177 stress 0.09178283 
    ## Run 178 stress 0.1076304 
    ## Run 179 stress 0.08938568 
    ## ... Procrustes: rmse 0.01166547  max resid 0.03940878 
    ## Run 180 stress 0.1092094 
    ## Run 181 stress 0.08938538 
    ## ... Procrustes: rmse 0.01174853  max resid 0.03950485 
    ## Run 182 stress 0.09109108 
    ## Run 183 stress 0.1071276 
    ## Run 184 stress 0.1056906 
    ## Run 185 stress 0.09262372 
    ## Run 186 stress 0.09021165 
    ## Run 187 stress 0.08946345 
    ## ... Procrustes: rmse 0.03744794  max resid 0.1177007 
    ## Run 188 stress 0.08926103 
    ## ... Procrustes: rmse 0.0003839538  max resid 0.001103222 
    ## ... Similar to previous best
    ## Run 189 stress 0.1074433 
    ## Run 190 stress 0.09099529 
    ## Run 191 stress 0.09087358 
    ## Run 192 stress 0.1092128 
    ## Run 193 stress 0.08938966 
    ## ... Procrustes: rmse 0.03595724  max resid 0.1183241 
    ## Run 194 stress 0.1056895 
    ## Run 195 stress 0.1092093 
    ## Run 196 stress 0.09087353 
    ## Run 197 stress 0.1075422 
    ## Run 198 stress 0.08926083 
    ## ... Procrustes: rmse 0.0002822608  max resid 0.0008089493 
    ## ... Similar to previous best
    ## Run 199 stress 0.08926104 
    ## ... Procrustes: rmse 0.0004224985  max resid 0.001110308 
    ## ... Similar to previous best
    ## Run 200 stress 0.1071322 
    ## Run 201 stress 0.09099552 
    ## Run 202 stress 0.1092129 
    ## Run 203 stress 0.1101443 
    ## Run 204 stress 0.09088606 
    ## Run 205 stress 0.1076303 
    ## Run 206 stress 0.1108231 
    ## Run 207 stress 0.1074432 
    ## Run 208 stress 0.08926077 
    ## ... Procrustes: rmse 0.0002455338  max resid 0.0006530326 
    ## ... Similar to previous best
    ## Run 209 stress 0.08938542 
    ## ... Procrustes: rmse 0.01179441  max resid 0.03957507 
    ## Run 210 stress 0.08946652 
    ## ... Procrustes: rmse 0.03318792  max resid 0.1163182 
    ## Run 211 stress 0.08946663 
    ## ... Procrustes: rmse 0.03316676  max resid 0.1162861 
    ## Run 212 stress 0.105265 
    ## Run 213 stress 0.105265 
    ## Run 214 stress 0.1091449 
    ## Run 215 stress 0.09612787 
    ## Run 216 stress 0.1074433 
    ## Run 217 stress 0.1064429 
    ## Run 218 stress 0.1071323 
    ## Run 219 stress 0.09087353 
    ## Run 220 stress 0.09592151 
    ## Run 221 stress 0.09503459 
    ## Run 222 stress 0.09018708 
    ## Run 223 stress 0.1096133 
    ## Run 224 stress 0.1079014 
    ## Run 225 stress 0.1071036 
    ## Run 226 stress 0.08938542 
    ## ... Procrustes: rmse 0.01172907  max resid 0.03953338 
    ## Run 227 stress 0.09503413 
    ## Run 228 stress 0.09039134 
    ## Run 229 stress 0.09503451 
    ## Run 230 stress 0.09039129 
    ## Run 231 stress 0.09109108 
    ## Run 232 stress 0.08938569 
    ## ... Procrustes: rmse 0.01184475  max resid 0.0394626 
    ## Run 233 stress 0.09592136 
    ## Run 234 stress 0.0892608 
    ## ... Procrustes: rmse 0.000297911  max resid 0.001002904 
    ## ... Similar to previous best
    ## Run 235 stress 0.09087366 
    ## Run 236 stress 0.0910911 
    ## Run 237 stress 0.1075422 
    ## Run 238 stress 0.08938537 
    ## ... Procrustes: rmse 0.01177195  max resid 0.03953434 
    ## Run 239 stress 0.1092093 
    ## Run 240 stress 0.1067379 
    ## Run 241 stress 0.0903914 
    ## Run 242 stress 0.08946651 
    ## ... Procrustes: rmse 0.03319383  max resid 0.1163272 
    ## Run 243 stress 0.08946333 
    ## ... Procrustes: rmse 0.0374654  max resid 0.1177269 
    ## Run 244 stress 0.08926091 
    ## ... Procrustes: rmse 0.0003340678  max resid 0.0009427162 
    ## ... Similar to previous best
    ## Run 245 stress 0.1105338 
    ## Run 246 stress 0.08938965 
    ## ... Procrustes: rmse 0.03595823  max resid 0.1183238 
    ## Run 247 stress 0.08926076 
    ## ... Procrustes: rmse 0.0002399818  max resid 0.0005734425 
    ## ... Similar to previous best
    ## Run 248 stress 0.08926068 
    ## ... Procrustes: rmse 0.0001666015  max resid 0.0004252386 
    ## ... Similar to previous best
    ## Run 249 stress 0.1068392 
    ## Run 250 stress 0.09503465 
    ## Run 251 stress 0.08938966 
    ## ... Procrustes: rmse 0.03595712  max resid 0.1183203 
    ## Run 252 stress 0.1075421 
    ## Run 253 stress 0.0893898 
    ## ... Procrustes: rmse 0.03589217  max resid 0.1182125 
    ## Run 254 stress 0.08946659 
    ## ... Procrustes: rmse 0.03317634  max resid 0.1163027 
    ## Run 255 stress 0.09039136 
    ## Run 256 stress 0.1074432 
    ## Run 257 stress 0.1091825 
    ## Run 258 stress 0.08946328 
    ## ... Procrustes: rmse 0.03749753  max resid 0.1177807 
    ## Run 259 stress 0.08946652 
    ## ... Procrustes: rmse 0.03321462  max resid 0.1163712 
    ## Run 260 stress 0.1052649 
    ## Run 261 stress 0.08938961 
    ## ... Procrustes: rmse 0.03593919  max resid 0.1182858 
    ## Run 262 stress 0.1075421 
    ## Run 263 stress 0.1076305 
    ## Run 264 stress 0.08938964 
    ## ... Procrustes: rmse 0.03591857  max resid 0.1182576 
    ## Run 265 stress 0.09592142 
    ## Run 266 stress 0.1092127 
    ## Run 267 stress 0.09775589 
    ## Run 268 stress 0.08926096 
    ## ... Procrustes: rmse 0.0003399627  max resid 0.000993494 
    ## ... Similar to previous best
    ## Run 269 stress 0.09021165 
    ## Run 270 stress 0.1092097 
    ## Run 271 stress 0.1075422 
    ## Run 272 stress 0.09018707 
    ## Run 273 stress 0.09130088 
    ## Run 274 stress 0.08926097 
    ## ... Procrustes: rmse 0.0003815153  max resid 0.001037379 
    ## ... Similar to previous best
    ## Run 275 stress 0.1071029 
    ## Run 276 stress 0.0892609 
    ## ... Procrustes: rmse 0.0002365032  max resid 0.0005446289 
    ## ... Similar to previous best
    ## Run 277 stress 0.1079017 
    ## Run 278 stress 0.08926065 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001305301  max resid 0.000355764 
    ## ... Similar to previous best
    ## Run 279 stress 0.10613 
    ## Run 280 stress 0.1071321 
    ## Run 281 stress 0.1108321 
    ## Run 282 stress 0.08926066 
    ## ... Procrustes: rmse 0.0001727361  max resid 0.0005646264 
    ## ... Similar to previous best
    ## Run 283 stress 0.09503431 
    ## Run 284 stress 0.1080638 
    ## Run 285 stress 0.1052649 
    ## Run 286 stress 0.1066195 
    ## Run 287 stress 0.08938558 
    ## ... Procrustes: rmse 0.01177454  max resid 0.03968064 
    ## Run 288 stress 0.08946651 
    ## ... Procrustes: rmse 0.03324886  max resid 0.1163992 
    ## Run 289 stress 0.08938544 
    ## ... Procrustes: rmse 0.01185371  max resid 0.03976408 
    ## Run 290 stress 0.08946335 
    ## ... Procrustes: rmse 0.03750173  max resid 0.1177508 
    ## Run 291 stress 0.1056901 
    ## Run 292 stress 0.09775573 
    ## Run 293 stress 0.08926079 
    ## ... Procrustes: rmse 0.0001601573  max resid 0.0005223506 
    ## ... Similar to previous best
    ## Run 294 stress 0.1092127 
    ## Run 295 stress 0.1091452 
    ## Run 296 stress 0.1085129 
    ## Run 297 stress 0.1060433 
    ## Run 298 stress 0.09018709 
    ## Run 299 stress 0.1063192 
    ## Run 300 stress 0.08938548 
    ## ... Procrustes: rmse 0.01184099  max resid 0.0396396 
    ## Run 301 stress 0.1052647 
    ## Run 302 stress 0.08938981 
    ## ... Procrustes: rmse 0.03592104  max resid 0.1182358 
    ## Run 303 stress 0.08926106 
    ## ... Procrustes: rmse 0.0002929072  max resid 0.0009596558 
    ## ... Similar to previous best
    ## Run 304 stress 0.09109109 
    ## Run 305 stress 0.09503413 
    ## Run 306 stress 0.105265 
    ## Run 307 stress 0.09178282 
    ## Run 308 stress 0.08938545 
    ## ... Procrustes: rmse 0.01185185  max resid 0.03976886 
    ## Run 309 stress 0.08946335 
    ## ... Procrustes: rmse 0.03750263  max resid 0.117754 
    ## Run 310 stress 0.08938961 
    ## ... Procrustes: rmse 0.0359666  max resid 0.1183139 
    ## Run 311 stress 0.09088601 
    ## Run 312 stress 0.1063198 
    ## Run 313 stress 0.08946334 
    ## ... Procrustes: rmse 0.03750386  max resid 0.1177526 
    ## Run 314 stress 0.1074432 
    ## Run 315 stress 0.1071322 
    ## Run 316 stress 0.08938544 
    ## ... Procrustes: rmse 0.01178249  max resid 0.03970106 
    ## Run 317 stress 0.1062565 
    ## Run 318 stress 0.1092127 
    ## Run 319 stress 0.1065383 
    ## Run 320 stress 0.1080637 
    ## Run 321 stress 0.09087349 
    ## Run 322 stress 0.08926118 
    ## ... Procrustes: rmse 0.0003819259  max resid 0.001270943 
    ## ... Similar to previous best
    ## Run 323 stress 0.10569 
    ## Run 324 stress 0.1062567 
    ## Run 325 stress 0.08938962 
    ## ... Procrustes: rmse 0.03597976  max resid 0.1183292 
    ## Run 326 stress 0.1056895 
    ## Run 327 stress 0.0903913 
    ## Run 328 stress 0.1109748 
    ## Run 329 stress 0.09503417 
    ## Run 330 stress 0.0918089 
    ## Run 331 stress 0.1092361 
    ## Run 332 stress 0.1071322 
    ## Run 333 stress 0.09503429 
    ## Run 334 stress 0.1056893 
    ## Run 335 stress 0.08926084 
    ## ... Procrustes: rmse 0.0001951754  max resid 0.0006165859 
    ## ... Similar to previous best
    ## Run 336 stress 0.1096139 
    ## Run 337 stress 0.0893854 
    ## ... Procrustes: rmse 0.01183289  max resid 0.03974482 
    ## Run 338 stress 0.08946654 
    ## ... Procrustes: rmse 0.03321875  max resid 0.1163396 
    ## Run 339 stress 0.08946654 
    ## ... Procrustes: rmse 0.03325344  max resid 0.1164152 
    ## Run 340 stress 0.0913009 
    ## Run 341 stress 0.105265 
    ## Run 342 stress 0.09592134 
    ## Run 343 stress 0.109561 
    ## Run 344 stress 0.09503445 
    ## Run 345 stress 0.09503417 
    ## Run 346 stress 0.08946342 
    ## ... Procrustes: rmse 0.03748922  max resid 0.1177326 
    ## Run 347 stress 0.08926099 
    ## ... Procrustes: rmse 0.0002756386  max resid 0.0008178076 
    ## ... Similar to previous best
    ## Run 348 stress 0.08938548 
    ## ... Procrustes: rmse 0.01175666  max resid 0.03963526 
    ## Run 349 stress 0.1056895 
    ## Run 350 stress 0.08938962 
    ## ... Procrustes: rmse 0.03595886  max resid 0.1182976 
    ## Run 351 stress 0.1080637 
    ## Run 352 stress 0.08938543 
    ## ... Procrustes: rmse 0.01182035  max resid 0.03975308 
    ## Run 353 stress 0.08938972 
    ## ... Procrustes: rmse 0.03593224  max resid 0.1182562 
    ## Run 354 stress 0.08926105 
    ## ... Procrustes: rmse 0.0003162579  max resid 0.0009866743 
    ## ... Similar to previous best
    ## Run 355 stress 0.08938548 
    ## ... Procrustes: rmse 0.01175031  max resid 0.03964099 
    ## Run 356 stress 0.1101455 
    ## Run 357 stress 0.1075423 
    ## Run 358 stress 0.1066196 
    ## Run 359 stress 0.09109113 
    ## Run 360 stress 0.09503452 
    ## Run 361 stress 0.09503433 
    ## Run 362 stress 0.08946651 
    ## ... Procrustes: rmse 0.03324554  max resid 0.1163868 
    ## Run 363 stress 0.08938961 
    ## ... Procrustes: rmse 0.0359748  max resid 0.118323 
    ## Run 364 stress 0.1056898 
    ## Run 365 stress 0.08926078 
    ## ... Procrustes: rmse 0.0001526577  max resid 0.0005092658 
    ## ... Similar to previous best
    ## Run 366 stress 0.09503416 
    ## Run 367 stress 0.09130098 
    ## Run 368 stress 0.09087348 
    ## Run 369 stress 0.08926068 
    ## ... Procrustes: rmse 0.0001282327  max resid 0.0003620019 
    ## ... Similar to previous best
    ## Run 370 stress 0.09503428 
    ## Run 371 stress 0.09130087 
    ## Run 372 stress 0.08946329 
    ## ... Procrustes: rmse 0.03752219  max resid 0.1177802 
    ## Run 373 stress 0.08938546 
    ## ... Procrustes: rmse 0.01187289  max resid 0.03968169 
    ## Run 374 stress 0.1076305 
    ## Run 375 stress 0.1074432 
    ## Run 376 stress 0.08946655 
    ## ... Procrustes: rmse 0.03325468  max resid 0.116419 
    ## Run 377 stress 0.09039132 
    ## Run 378 stress 0.08926071 
    ## ... Procrustes: rmse 0.0002618557  max resid 0.0006628543 
    ## ... Similar to previous best
    ## Run 379 stress 0.09021166 
    ## Run 380 stress 0.1105345 
    ## Run 381 stress 0.08946652 
    ## ... Procrustes: rmse 0.03322279  max resid 0.1163466 
    ## Run 382 stress 0.1079017 
    ## Run 383 stress 0.1092127 
    ## Run 384 stress 0.1056902 
    ## Run 385 stress 0.1060431 
    ## Run 386 stress 0.1065366 
    ## Run 387 stress 0.08938964 
    ## ... Procrustes: rmse 0.03594913  max resid 0.1182823 
    ## Run 388 stress 0.1056907 
    ## Run 389 stress 0.0892608 
    ## ... Procrustes: rmse 0.0001634694  max resid 0.0005117664 
    ## ... Similar to previous best
    ## Run 390 stress 0.1071322 
    ## Run 391 stress 0.08946335 
    ## ... Procrustes: rmse 0.03756453  max resid 0.1178485 
    ## Run 392 stress 0.08938975 
    ## ... Procrustes: rmse 0.0359327  max resid 0.1182582 
    ## Run 393 stress 0.0893854 
    ## ... Procrustes: rmse 0.01184293  max resid 0.0397458 
    ## Run 394 stress 0.09039131 
    ## Run 395 stress 0.08938985 
    ## ... Procrustes: rmse 0.03592191  max resid 0.1182363 
    ## Run 396 stress 0.1111026 
    ## Run 397 stress 0.09039131 
    ## Run 398 stress 0.1075422 
    ## Run 399 stress 0.1101448 
    ## Run 400 stress 0.1056903 
    ## Run 401 stress 0.08938548 
    ## ... Procrustes: rmse 0.01176924  max resid 0.03973408 
    ## Run 402 stress 0.1091827 
    ## Run 403 stress 0.1060434 
    ## Run 404 stress 0.1056902 
    ## Run 405 stress 0.1079018 
    ## Run 406 stress 0.1074434 
    ## Run 407 stress 0.0892607 
    ## ... Procrustes: rmse 7.092829e-05  max resid 0.000223247 
    ## ... Similar to previous best
    ## Run 408 stress 0.08938549 
    ## ... Procrustes: rmse 0.01177025  max resid 0.0397385 
    ## Run 409 stress 0.08946338 
    ## ... Procrustes: rmse 0.03749606  max resid 0.1177371 
    ## Run 410 stress 0.1062568 
    ## Run 411 stress 0.1086565 
    ## Run 412 stress 0.09503414 
    ## Run 413 stress 0.1076305 
    ## Run 414 stress 0.09180884 
    ## Run 415 stress 0.1075421 
    ## Run 416 stress 0.1061304 
    ## Run 417 stress 0.1087573 
    ## Run 418 stress 0.1092096 
    ## Run 419 stress 0.1056907 
    ## Run 420 stress 0.09021165 
    ## Run 421 stress 0.08926101 
    ## ... Procrustes: rmse 0.0002934628  max resid 0.0009795639 
    ## ... Similar to previous best
    ## Run 422 stress 0.107112 
    ## Run 423 stress 0.1061303 
    ## Run 424 stress 0.0895172 
    ## ... Procrustes: rmse 0.03505351  max resid 0.1156294 
    ## Run 425 stress 0.08938975 
    ## ... Procrustes: rmse 0.03592897  max resid 0.1182514 
    ## Run 426 stress 0.1056903 
    ## Run 427 stress 0.1065358 
    ## Run 428 stress 0.08938963 
    ## ... Procrustes: rmse 0.03598265  max resid 0.1183375 
    ## Run 429 stress 0.08938539 
    ## ... Procrustes: rmse 0.01181823  max resid 0.03971542 
    ## Run 430 stress 0.09018707 
    ## Run 431 stress 0.1091827 
    ## Run 432 stress 0.08946332 
    ## ... Procrustes: rmse 0.03755022  max resid 0.1178243 
    ## Run 433 stress 0.09503451 
    ## Run 434 stress 0.08926074 
    ## ... Procrustes: rmse 0.0002022587  max resid 0.0005320922 
    ## ... Similar to previous best
    ## Run 435 stress 0.08946328 
    ## ... Procrustes: rmse 0.03752994  max resid 0.1177971 
    ## Run 436 stress 0.0908734 
    ## Run 437 stress 0.09503437 
    ## Run 438 stress 0.08938538 
    ## ... Procrustes: rmse 0.01182668  max resid 0.03972473 
    ## Run 439 stress 0.09503461 
    ## Run 440 stress 0.1124789 
    ## Run 441 stress 0.09039132 
    ## Run 442 stress 0.106443 
    ## Run 443 stress 0.1052646 
    ## Run 444 stress 0.09087358 
    ## Run 445 stress 0.1087862 
    ## Run 446 stress 0.08926083 
    ## ... Procrustes: rmse 0.0001764419  max resid 0.0005775209 
    ## ... Similar to previous best
    ## Run 447 stress 0.08951721 
    ## ... Procrustes: rmse 0.03510164  max resid 0.1156917 
    ## Run 448 stress 0.08926113 
    ## ... Procrustes: rmse 0.0003295303  max resid 0.001092904 
    ## ... Similar to previous best
    ## Run 449 stress 0.08946665 
    ## ... Procrustes: rmse 0.03319826  max resid 0.1163082 
    ## Run 450 stress 0.09099542 
    ## Run 451 stress 0.1056901 
    ## Run 452 stress 0.1056892 
    ## Run 453 stress 0.1076303 
    ## Run 454 stress 0.08926116 
    ## ... Procrustes: rmse 0.0003657158  max resid 0.001149314 
    ## ... Similar to previous best
    ## Run 455 stress 0.111895 
    ## Run 456 stress 0.106009 
    ## Run 457 stress 0.1075421 
    ## Run 458 stress 0.1080639 
    ## Run 459 stress 0.09039136 
    ## Run 460 stress 0.0950347 
    ## Run 461 stress 0.08951718 
    ## ... Procrustes: rmse 0.03506349  max resid 0.1156485 
    ## Run 462 stress 0.08938539 
    ## ... Procrustes: rmse 0.0118222  max resid 0.03973671 
    ## Run 463 stress 0.08926086 
    ## ... Procrustes: rmse 0.0001506426  max resid 0.0004700172 
    ## ... Similar to previous best
    ## Run 464 stress 0.1067511 
    ## Run 465 stress 0.09087362 
    ## Run 466 stress 0.0903913 
    ## Run 467 stress 0.09178283 
    ## Run 468 stress 0.09087339 
    ## Run 469 stress 0.1060428 
    ## Run 470 stress 0.08938963 
    ## ... Procrustes: rmse 0.03595232  max resid 0.1182898 
    ## Run 471 stress 0.08938539 
    ## ... Procrustes: rmse 0.01180596  max resid 0.03972636 
    ## Run 472 stress 0.08926076 
    ## ... Procrustes: rmse 0.0001280225  max resid 0.0004026126 
    ## ... Similar to previous best
    ## Run 473 stress 0.09021168 
    ## Run 474 stress 0.08938962 
    ## ... Procrustes: rmse 0.03598134  max resid 0.118333 
    ## Run 475 stress 0.1104177 
    ## Run 476 stress 0.1065382 
    ## Run 477 stress 0.0893854 
    ## ... Procrustes: rmse 0.01178669  max resid 0.03968158 
    ## Run 478 stress 0.08938964 
    ## ... Procrustes: rmse 0.03594971  max resid 0.1182823 
    ## Run 479 stress 0.1074434 
    ## Run 480 stress 0.08938538 
    ## ... Procrustes: rmse 0.01182626  max resid 0.03972385 
    ## Run 481 stress 0.08938964 
    ## ... Procrustes: rmse 0.03595259  max resid 0.1182844 
    ## Run 482 stress 0.09018707 
    ## Run 483 stress 0.1103711 
    ## Run 484 stress 0.08946332 
    ## ... Procrustes: rmse 0.03750916  max resid 0.1177603 
    ## Run 485 stress 0.08938965 
    ## ... Procrustes: rmse 0.03594729  max resid 0.1182802 
    ## Run 486 stress 0.1074433 
    ## Run 487 stress 0.0961283 
    ## Run 488 stress 0.08926077 
    ## ... Procrustes: rmse 0.000138886  max resid 0.0004415683 
    ## ... Similar to previous best
    ## Run 489 stress 0.08926074 
    ## ... Procrustes: rmse 0.0001178802  max resid 0.0003825041 
    ## ... Similar to previous best
    ## Run 490 stress 0.1060429 
    ## Run 491 stress 0.08926073 
    ## ... Procrustes: rmse 0.0001022436  max resid 0.0003219412 
    ## ... Similar to previous best
    ## Run 492 stress 0.08946336 
    ## ... Procrustes: rmse 0.0375008  max resid 0.1177468 
    ## Run 493 stress 0.1076303 
    ## Run 494 stress 0.1066201 
    ## Run 495 stress 0.09039146 
    ## Run 496 stress 0.1066192 
    ## Run 497 stress 0.08938965 
    ## ... Procrustes: rmse 0.03598976  max resid 0.1183518 
    ## Run 498 stress 0.1075421 
    ## Run 499 stress 0.1101453 
    ## Run 500 stress 0.1086567 
    ## *** Best solution repeated 23 times

``` r
# Mixed and stratified lakes
SD_beta_geo_MS_NMDS <- metaMDS(SD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09407452 
    ## Run 2 stress 0.09030398 
    ## Run 3 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01442582  max resid 0.04279712 
    ## Run 4 stress 0.0940797 
    ## Run 5 stress 0.09145321 
    ## Run 6 stress 0.09159089 
    ## Run 7 stress 0.09407975 
    ## Run 8 stress 0.09159086 
    ## Run 9 stress 0.0946591 
    ## Run 10 stress 0.09381855 
    ## Run 11 stress 0.09969992 
    ## Run 12 stress 0.09308987 
    ## Run 13 stress 0.09159086 
    ## Run 14 stress 0.09969991 
    ## Run 15 stress 0.09030395 
    ## Run 16 stress 0.09407972 
    ## Run 17 stress 0.0926837 
    ## Run 18 stress 0.1054535 
    ## Run 19 stress 0.09337231 
    ## Run 20 stress 0.08503551 
    ## Run 21 stress 0.09159099 
    ## Run 22 stress 0.09308949 
    ## Run 23 stress 0.09145328 
    ## Run 24 stress 0.09404394 
    ## Run 25 stress 0.09408005 
    ## Run 26 stress 0.09159086 
    ## Run 27 stress 0.09159093 
    ## Run 28 stress 0.09969965 
    ## Run 29 stress 0.09145337 
    ## Run 30 stress 0.09145323 
    ## Run 31 stress 0.09308958 
    ## Run 32 stress 0.09286084 
    ## Run 33 stress 0.09407986 
    ## Run 34 stress 0.09969967 
    ## Run 35 stress 0.09590559 
    ## Run 36 stress 0.09030418 
    ## Run 37 stress 0.09969995 
    ## Run 38 stress 0.08440265 
    ## ... Procrustes: rmse 0.000267192  max resid 0.0004766824 
    ## ... Similar to previous best
    ## Run 39 stress 0.08973874 
    ## Run 40 stress 0.09286093 
    ## Run 41 stress 0.09407979 
    ## Run 42 stress 0.09030393 
    ## Run 43 stress 0.09445477 
    ## Run 44 stress 0.09030393 
    ## Run 45 stress 0.09308995 
    ## Run 46 stress 0.09969983 
    ## Run 47 stress 0.0877347 
    ## Run 48 stress 0.09400534 
    ## Run 49 stress 0.09030397 
    ## Run 50 stress 0.08503487 
    ## Run 51 stress 0.09400516 
    ## Run 52 stress 0.09380575 
    ## Run 53 stress 0.09030412 
    ## Run 54 stress 0.09416133 
    ## Run 55 stress 0.08503472 
    ## Run 56 stress 0.08440262 
    ## ... Procrustes: rmse 0.0002501382  max resid 0.0004903714 
    ## ... Similar to previous best
    ## Run 57 stress 0.08503478 
    ## Run 58 stress 0.0930897 
    ## Run 59 stress 0.09669868 
    ## Run 60 stress 0.09308946 
    ## Run 61 stress 0.09374381 
    ## Run 62 stress 0.09030398 
    ## Run 63 stress 0.08503592 
    ## Run 64 stress 0.09407416 
    ## Run 65 stress 0.09159083 
    ## Run 66 stress 0.08973892 
    ## Run 67 stress 0.09030401 
    ## Run 68 stress 0.08503469 
    ## Run 69 stress 0.09374281 
    ## Run 70 stress 0.08503511 
    ## Run 71 stress 0.09030393 
    ## Run 72 stress 0.0940798 
    ## Run 73 stress 0.09268349 
    ## Run 74 stress 0.09145338 
    ## Run 75 stress 0.09417784 
    ## Run 76 stress 0.08973865 
    ## Run 77 stress 0.08773489 
    ## Run 78 stress 0.09030397 
    ## Run 79 stress 0.09030394 
    ## Run 80 stress 0.09408 
    ## Run 81 stress 0.09337234 
    ## Run 82 stress 0.09416139 
    ## Run 83 stress 0.09407992 
    ## Run 84 stress 0.09308975 
    ## Run 85 stress 0.09969982 
    ## Run 86 stress 0.08440254 
    ## ... Procrustes: rmse 0.0001600029  max resid 0.0002988514 
    ## ... Similar to previous best
    ## Run 87 stress 0.09380585 
    ## Run 88 stress 0.09407389 
    ## Run 89 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001836494  max resid 0.0003603279 
    ## ... Similar to previous best
    ## Run 90 stress 0.09969988 
    ## Run 91 stress 0.09308979 
    ## Run 92 stress 0.09760824 
    ## Run 93 stress 0.1052848 
    ## Run 94 stress 0.08440261 
    ## ... Procrustes: rmse 0.0001427644  max resid 0.0002832826 
    ## ... Similar to previous best
    ## Run 95 stress 0.09407966 
    ## Run 96 stress 0.1038043 
    ## Run 97 stress 0.09416131 
    ## Run 98 stress 0.09445463 
    ## Run 99 stress 0.08503458 
    ## Run 100 stress 0.09969964 
    ## Run 101 stress 0.08503682 
    ## Run 102 stress 0.09168944 
    ## Run 103 stress 0.09539192 
    ## Run 104 stress 0.09030407 
    ## Run 105 stress 0.09030414 
    ## Run 106 stress 0.08773469 
    ## Run 107 stress 0.09465907 
    ## Run 108 stress 0.09535535 
    ## Run 109 stress 0.08973875 
    ## Run 110 stress 0.09464484 
    ## Run 111 stress 0.1113784 
    ## Run 112 stress 0.09464404 
    ## Run 113 stress 0.09159083 
    ## Run 114 stress 0.09381854 
    ## Run 115 stress 0.100461 
    ## Run 116 stress 0.1052849 
    ## Run 117 stress 0.09168927 
    ## Run 118 stress 0.09445491 
    ## Run 119 stress 0.09168942 
    ## Run 120 stress 0.08973896 
    ## Run 121 stress 0.09407143 
    ## Run 122 stress 0.09969982 
    ## Run 123 stress 0.09416131 
    ## Run 124 stress 0.09465912 
    ## Run 125 stress 0.09337277 
    ## Run 126 stress 0.09416139 
    ## Run 127 stress 0.09268319 
    ## Run 128 stress 0.3125725 
    ## Run 129 stress 0.0937423 
    ## Run 130 stress 0.09337258 
    ## Run 131 stress 0.09286085 
    ## Run 132 stress 0.09337298 
    ## Run 133 stress 0.09268392 
    ## Run 134 stress 0.08503527 
    ## Run 135 stress 0.09381847 
    ## Run 136 stress 0.09407968 
    ## Run 137 stress 0.0940799 
    ## Run 138 stress 0.09308947 
    ## Run 139 stress 0.09168936 
    ## Run 140 stress 0.08773478 
    ## Run 141 stress 0.09416131 
    ## Run 142 stress 0.09159084 
    ## Run 143 stress 0.09268347 
    ## Run 144 stress 0.0850349 
    ## Run 145 stress 0.09286084 
    ## Run 146 stress 0.1038044 
    ## Run 147 stress 0.09268412 
    ## Run 148 stress 0.09030403 
    ## Run 149 stress 0.09268343 
    ## Run 150 stress 0.09465907 
    ## Run 151 stress 0.09760813 
    ## Run 152 stress 0.09145325 
    ## Run 153 stress 0.08503529 
    ## Run 154 stress 0.09268364 
    ## Run 155 stress 0.08440257 
    ## ... Procrustes: rmse 0.0001894823  max resid 0.0003526914 
    ## ... Similar to previous best
    ## Run 156 stress 0.09400526 
    ## Run 157 stress 0.09030404 
    ## Run 158 stress 0.09145322 
    ## Run 159 stress 0.09030414 
    ## Run 160 stress 0.08503472 
    ## Run 161 stress 0.09030395 
    ## Run 162 stress 0.09030411 
    ## Run 163 stress 0.09321833 
    ## Run 164 stress 0.09721201 
    ## Run 165 stress 0.09535524 
    ## Run 166 stress 0.08440253 
    ## ... Procrustes: rmse 3.495558e-05  max resid 7.550402e-05 
    ## ... Similar to previous best
    ## Run 167 stress 0.09337273 
    ## Run 168 stress 0.09159096 
    ## Run 169 stress 0.09030403 
    ## Run 170 stress 0.09286085 
    ## Run 171 stress 0.0937425 
    ## Run 172 stress 0.08773469 
    ## Run 173 stress 0.3047121 
    ## Run 174 stress 0.09407986 
    ## Run 175 stress 0.09381846 
    ## Run 176 stress 0.09168931 
    ## Run 177 stress 0.08773474 
    ## Run 178 stress 0.09969996 
    ## Run 179 stress 0.09590892 
    ## Run 180 stress 0.09407978 
    ## Run 181 stress 0.09539201 
    ## Run 182 stress 0.09669829 
    ## Run 183 stress 0.08973888 
    ## Run 184 stress 0.08440272 
    ## ... Procrustes: rmse 0.0003173712  max resid 0.0005796054 
    ## ... Similar to previous best
    ## Run 185 stress 0.09308955 
    ## Run 186 stress 0.09286083 
    ## Run 187 stress 0.1052851 
    ## Run 188 stress 0.09416138 
    ## Run 189 stress 0.09286086 
    ## Run 190 stress 0.09268326 
    ## Run 191 stress 0.09337241 
    ## Run 192 stress 0.09308966 
    ## Run 193 stress 0.09159087 
    ## Run 194 stress 0.09380578 
    ## Run 195 stress 0.08503473 
    ## Run 196 stress 0.09407968 
    ## Run 197 stress 0.09337286 
    ## Run 198 stress 0.09407984 
    ## Run 199 stress 0.08773479 
    ## Run 200 stress 0.0946591 
    ## Run 201 stress 0.09286083 
    ## Run 202 stress 0.09168939 
    ## Run 203 stress 0.09407985 
    ## Run 204 stress 0.0926838 
    ## Run 205 stress 0.09465904 
    ## Run 206 stress 0.08440262 
    ## ... Procrustes: rmse 0.0002529243  max resid 0.000464642 
    ## ... Similar to previous best
    ## Run 207 stress 0.08773469 
    ## Run 208 stress 0.095392 
    ## Run 209 stress 0.08503494 
    ## Run 210 stress 0.09030398 
    ## Run 211 stress 0.09407996 
    ## Run 212 stress 0.09286088 
    ## Run 213 stress 0.08503662 
    ## Run 214 stress 0.08503512 
    ## Run 215 stress 0.09407986 
    ## Run 216 stress 0.09268392 
    ## Run 217 stress 0.0903041 
    ## Run 218 stress 0.0877347 
    ## Run 219 stress 0.08973862 
    ## Run 220 stress 0.09030399 
    ## Run 221 stress 0.09407966 
    ## Run 222 stress 0.09465905 
    ## Run 223 stress 0.09268373 
    ## Run 224 stress 0.1053432 
    ## Run 225 stress 0.09374303 
    ## Run 226 stress 0.08503507 
    ## Run 227 stress 0.1038049 
    ## Run 228 stress 0.1038035 
    ## Run 229 stress 0.09407966 
    ## Run 230 stress 0.09030402 
    ## Run 231 stress 0.08973886 
    ## Run 232 stress 0.09416137 
    ## Run 233 stress 0.093373 
    ## Run 234 stress 0.09407969 
    ## Run 235 stress 0.0844028 
    ## ... Procrustes: rmse 0.0003726764  max resid 0.0006710347 
    ## ... Similar to previous best
    ## Run 236 stress 0.08773468 
    ## Run 237 stress 0.095355 
    ## Run 238 stress 0.100461 
    ## Run 239 stress 0.09337275 
    ## Run 240 stress 0.09268387 
    ## Run 241 stress 0.09969968 
    ## Run 242 stress 0.1038041 
    ## Run 243 stress 0.08973877 
    ## Run 244 stress 0.09168942 
    ## Run 245 stress 0.09337275 
    ## Run 246 stress 0.09030403 
    ## Run 247 stress 0.08973871 
    ## Run 248 stress 0.09669869 
    ## Run 249 stress 0.09535426 
    ## Run 250 stress 0.08503677 
    ## Run 251 stress 0.08503461 
    ## Run 252 stress 0.09969972 
    ## Run 253 stress 0.09308946 
    ## Run 254 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001626402  max resid 0.0003420487 
    ## ... Similar to previous best
    ## Run 255 stress 0.09030397 
    ## Run 256 stress 0.09416142 
    ## Run 257 stress 0.0844026 
    ## ... Procrustes: rmse 0.0002060907  max resid 0.0004246158 
    ## ... Similar to previous best
    ## Run 258 stress 0.341173 
    ## Run 259 stress 0.08973877 
    ## Run 260 stress 0.09374294 
    ## Run 261 stress 0.0850348 
    ## Run 262 stress 0.0953921 
    ## Run 263 stress 0.09030402 
    ## Run 264 stress 0.09407976 
    ## Run 265 stress 0.08440259 
    ## ... Procrustes: rmse 0.0002081767  max resid 0.0003811008 
    ## ... Similar to previous best
    ## Run 266 stress 0.09268349 
    ## Run 267 stress 0.09168954 
    ## Run 268 stress 0.09145322 
    ## Run 269 stress 0.08773469 
    ## Run 270 stress 0.09407974 
    ## Run 271 stress 0.09168931 
    ## Run 272 stress 0.08503573 
    ## Run 273 stress 0.09168934 
    ## Run 274 stress 0.09286098 
    ## Run 275 stress 0.09286087 
    ## Run 276 stress 0.09168926 
    ## Run 277 stress 0.09337272 
    ## Run 278 stress 0.08503485 
    ## Run 279 stress 0.08503541 
    ## Run 280 stress 0.08440262 
    ## ... Procrustes: rmse 0.0002446126  max resid 0.0004962473 
    ## ... Similar to previous best
    ## Run 281 stress 0.08503547 
    ## Run 282 stress 0.09308946 
    ## Run 283 stress 0.09145322 
    ## Run 284 stress 0.09400525 
    ## Run 285 stress 0.2461313 
    ## Run 286 stress 0.08440252 
    ## ... Procrustes: rmse 7.962986e-05  max resid 0.0001486167 
    ## ... Similar to previous best
    ## Run 287 stress 0.09969967 
    ## Run 288 stress 0.09969972 
    ## Run 289 stress 0.09407988 
    ## Run 290 stress 0.08773467 
    ## Run 291 stress 0.0916893 
    ## Run 292 stress 0.09407112 
    ## Run 293 stress 0.09374324 
    ## Run 294 stress 0.09407999 
    ## Run 295 stress 0.09407974 
    ## Run 296 stress 0.08973874 
    ## Run 297 stress 0.09539207 
    ## Run 298 stress 0.09407989 
    ## Run 299 stress 0.09760835 
    ## Run 300 stress 0.09168933 
    ## Run 301 stress 0.0937424 
    ## Run 302 stress 0.09268384 
    ## Run 303 stress 0.09407987 
    ## Run 304 stress 0.1038035 
    ## Run 305 stress 0.09308973 
    ## Run 306 stress 0.09403396 
    ## Run 307 stress 0.08973865 
    ## Run 308 stress 0.09465905 
    ## Run 309 stress 0.0941614 
    ## Run 310 stress 0.09403399 
    ## Run 311 stress 0.1052851 
    ## Run 312 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001726998  max resid 0.0003279366 
    ## ... Similar to previous best
    ## Run 313 stress 0.0940739 
    ## Run 314 stress 0.09408012 
    ## Run 315 stress 0.08503467 
    ## Run 316 stress 0.0850364 
    ## Run 317 stress 0.09408008 
    ## Run 318 stress 0.08973875 
    ## Run 319 stress 0.0897388 
    ## Run 320 stress 0.09969988 
    ## Run 321 stress 0.08973863 
    ## Run 322 stress 0.09416147 
    ## Run 323 stress 0.08503518 
    ## Run 324 stress 0.09721207 
    ## Run 325 stress 0.09321708 
    ## Run 326 stress 0.0933725 
    ## Run 327 stress 0.09286086 
    ## Run 328 stress 0.09308981 
    ## Run 329 stress 0.08773469 
    ## Run 330 stress 0.09337302 
    ## Run 331 stress 0.09590877 
    ## Run 332 stress 0.08503456 
    ## Run 333 stress 0.09407138 
    ## Run 334 stress 0.0844026 
    ## ... Procrustes: rmse 9.321933e-05  max resid 0.000175767 
    ## ... Similar to previous best
    ## Run 335 stress 0.09407964 
    ## Run 336 stress 0.09374324 
    ## Run 337 stress 0.08773465 
    ## Run 338 stress 0.09145332 
    ## Run 339 stress 0.09286084 
    ## Run 340 stress 0.09168959 
    ## Run 341 stress 0.09380588 
    ## Run 342 stress 0.09374377 
    ## Run 343 stress 0.08973867 
    ## Run 344 stress 0.0877349 
    ## Run 345 stress 0.09268366 
    ## Run 346 stress 0.09760821 
    ## Run 347 stress 0.09286089 
    ## Run 348 stress 0.09417799 
    ## Run 349 stress 0.08503465 
    ## Run 350 stress 0.0897387 
    ## Run 351 stress 0.08503466 
    ## Run 352 stress 0.09337276 
    ## Run 353 stress 0.09590579 
    ## Run 354 stress 0.09408002 
    ## Run 355 stress 0.08773478 
    ## Run 356 stress 0.08503591 
    ## Run 357 stress 0.08773485 
    ## Run 358 stress 0.08973884 
    ## Run 359 stress 0.09145321 
    ## Run 360 stress 0.09590539 
    ## Run 361 stress 0.09374324 
    ## Run 362 stress 0.0915909 
    ## Run 363 stress 0.09286087 
    ## Run 364 stress 0.09337325 
    ## Run 365 stress 0.09445482 
    ## Run 366 stress 0.09168949 
    ## Run 367 stress 0.08440254 
    ## ... Procrustes: rmse 0.0001422611  max resid 0.0002697313 
    ## ... Similar to previous best
    ## Run 368 stress 0.09030401 
    ## Run 369 stress 0.09407969 
    ## Run 370 stress 0.09268403 
    ## Run 371 stress 0.09030393 
    ## Run 372 stress 0.08773477 
    ## Run 373 stress 0.09416149 
    ## Run 374 stress 0.09416135 
    ## Run 375 stress 0.09030394 
    ## Run 376 stress 0.09407134 
    ## Run 377 stress 0.09408013 
    ## Run 378 stress 0.09030403 
    ## Run 379 stress 0.09416139 
    ## Run 380 stress 0.09381854 
    ## Run 381 stress 0.08503615 
    ## Run 382 stress 0.0877347 
    ## Run 383 stress 0.08440272 
    ## ... Procrustes: rmse 0.0003252325  max resid 0.000588525 
    ## ... Similar to previous best
    ## Run 384 stress 0.09308993 
    ## Run 385 stress 0.08503557 
    ## Run 386 stress 0.1038043 
    ## Run 387 stress 0.09407975 
    ## Run 388 stress 0.08503515 
    ## Run 389 stress 0.09465907 
    ## Run 390 stress 0.09721196 
    ## Run 391 stress 0.09408013 
    ## Run 392 stress 0.08503535 
    ## Run 393 stress 0.09286084 
    ## Run 394 stress 0.08503459 
    ## Run 395 stress 0.09407965 
    ## Run 396 stress 0.08973885 
    ## Run 397 stress 0.08503457 
    ## Run 398 stress 0.09268366 
    ## Run 399 stress 0.08973881 
    ## Run 400 stress 0.09268372 
    ## Run 401 stress 0.09286084 
    ## Run 402 stress 0.103805 
    ## Run 403 stress 0.090304 
    ## Run 404 stress 0.09159085 
    ## Run 405 stress 0.09145321 
    ## Run 406 stress 0.08503517 
    ## Run 407 stress 0.09030395 
    ## Run 408 stress 0.09145321 
    ## Run 409 stress 0.09030396 
    ## Run 410 stress 0.09145322 
    ## Run 411 stress 0.09969988 
    ## Run 412 stress 0.09445473 
    ## Run 413 stress 0.0953548 
    ## Run 414 stress 0.09380586 
    ## Run 415 stress 0.09465904 
    ## Run 416 stress 0.09308953 
    ## Run 417 stress 0.09337267 
    ## Run 418 stress 0.09030395 
    ## Run 419 stress 0.09030403 
    ## Run 420 stress 0.09159086 
    ## Run 421 stress 0.09159084 
    ## Run 422 stress 0.09447317 
    ## Run 423 stress 0.09337254 
    ## Run 424 stress 0.09286092 
    ## Run 425 stress 0.09145323 
    ## Run 426 stress 0.09030402 
    ## Run 427 stress 0.09407967 
    ## Run 428 stress 0.08440253 
    ## ... Procrustes: rmse 4.322859e-05  max resid 8.815631e-05 
    ## ... Similar to previous best
    ## Run 429 stress 0.09168946 
    ## Run 430 stress 0.08773468 
    ## Run 431 stress 0.09286087 
    ## Run 432 stress 0.0933725 
    ## Run 433 stress 0.08503482 
    ## Run 434 stress 0.09381848 
    ## Run 435 stress 0.08773468 
    ## Run 436 stress 0.0930899 
    ## Run 437 stress 0.09030398 
    ## Run 438 stress 0.09539195 
    ## Run 439 stress 0.09268424 
    ## Run 440 stress 0.1019602 
    ## Run 441 stress 0.09381855 
    ## Run 442 stress 0.09308969 
    ## Run 443 stress 0.09168943 
    ## Run 444 stress 0.08503597 
    ## Run 445 stress 0.09268337 
    ## Run 446 stress 0.09159086 
    ## Run 447 stress 0.08773467 
    ## Run 448 stress 0.09268403 
    ## Run 449 stress 0.09403404 
    ## Run 450 stress 0.09590632 
    ## Run 451 stress 0.09760821 
    ## Run 452 stress 0.09159083 
    ## Run 453 stress 0.09337244 
    ## Run 454 stress 0.08503462 
    ## Run 455 stress 0.09308946 
    ## Run 456 stress 0.09337286 
    ## Run 457 stress 0.09030393 
    ## Run 458 stress 0.0915909 
    ## Run 459 stress 0.1038044 
    ## Run 460 stress 0.08773469 
    ## Run 461 stress 0.08503474 
    ## Run 462 stress 0.08503507 
    ## Run 463 stress 0.09337244 
    ## Run 464 stress 0.1038048 
    ## Run 465 stress 0.08973874 
    ## Run 466 stress 0.09760814 
    ## Run 467 stress 0.09159083 
    ## Run 468 stress 0.08503456 
    ## Run 469 stress 0.09145321 
    ## Run 470 stress 0.091591 
    ## Run 471 stress 0.09380576 
    ## Run 472 stress 0.09168952 
    ## Run 473 stress 0.09286084 
    ## Run 474 stress 0.08773487 
    ## Run 475 stress 0.09286083 
    ## Run 476 stress 0.08973863 
    ## Run 477 stress 0.0926834 
    ## Run 478 stress 0.09030404 
    ## Run 479 stress 0.08773466 
    ## Run 480 stress 0.08503584 
    ## Run 481 stress 0.09286099 
    ## Run 482 stress 0.09416135 
    ## Run 483 stress 0.08773465 
    ## Run 484 stress 0.09145324 
    ## Run 485 stress 0.3126038 
    ## Run 486 stress 0.09030398 
    ## Run 487 stress 0.09760832 
    ## Run 488 stress 0.1072865 
    ## Run 489 stress 0.09145327 
    ## Run 490 stress 0.09969972 
    ## Run 491 stress 0.09268327 
    ## Run 492 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 6.298335e-05  max resid 0.0001055546 
    ## ... Similar to previous best
    ## Run 493 stress 0.09159093 
    ## Run 494 stress 0.08773484 
    ## Run 495 stress 0.09465905 
    ## Run 496 stress 0.09030393 
    ## Run 497 stress 0.09407995 
    ## Run 498 stress 0.09168928 
    ## Run 499 stress 0.0933723 
    ## Run 500 stress 0.09030394 
    ## *** Best solution repeated 1 times

``` r
# Ocean sites and mixed lakes
SD_beta_geo_OM_NMDS <- metaMDS(SD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07365783 
    ## Run 1 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744355  max resid 0.05105504 
    ## Run 2 stress 0.07365783 
    ## ... Procrustes: rmse 2.656544e-05  max resid 6.14824e-05 
    ## ... Similar to previous best
    ## Run 3 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744353  max resid 0.0510383 
    ## Run 4 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743902  max resid 0.05102231 
    ## Run 5 stress 0.08000473 
    ## Run 6 stress 0.08233751 
    ## Run 7 stress 0.07629238 
    ## Run 8 stress 0.07365783 
    ## ... Procrustes: rmse 3.34402e-05  max resid 7.724538e-05 
    ## ... Similar to previous best
    ## Run 9 stress 0.07365784 
    ## ... Procrustes: rmse 6.148183e-05  max resid 0.0001438381 
    ## ... Similar to previous best
    ## Run 10 stress 0.07732935 
    ## Run 11 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001592566  max resid 0.0003755693 
    ## ... Similar to previous best
    ## Run 12 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744274  max resid 0.05102854 
    ## Run 13 stress 0.07365784 
    ## ... Procrustes: rmse 5.043155e-05  max resid 0.0001189144 
    ## ... Similar to previous best
    ## Run 14 stress 0.08000468 
    ## Run 15 stress 0.0762924 
    ## Run 16 stress 0.07629242 
    ## Run 17 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174326  max resid 0.0509812 
    ## Run 18 stress 0.07365783 
    ## ... Procrustes: rmse 1.855696e-05  max resid 4.173907e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743963  max resid 0.05101747 
    ## Run 20 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744412  max resid 0.05105408 
    ## Run 21 stress 0.08233761 
    ## Run 22 stress 0.07629236 
    ## Run 23 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174434  max resid 0.05104925 
    ## Run 24 stress 0.07365785 
    ## ... Procrustes: rmse 4.670104e-05  max resid 9.920924e-05 
    ## ... Similar to previous best
    ## Run 25 stress 0.07378231 
    ## ... Procrustes: rmse 0.01743867  max resid 0.05104694 
    ## Run 26 stress 0.07629235 
    ## Run 27 stress 0.07732936 
    ## Run 28 stress 0.08000467 
    ## Run 29 stress 0.08233759 
    ## Run 30 stress 0.07629246 
    ## Run 31 stress 0.07629239 
    ## Run 32 stress 0.07629238 
    ## Run 33 stress 0.0800048 
    ## Run 34 stress 0.08233755 
    ## Run 35 stress 0.07629245 
    ## Run 36 stress 0.07629241 
    ## Run 37 stress 0.07365785 
    ## ... Procrustes: rmse 6.841234e-05  max resid 0.0001615611 
    ## ... Similar to previous best
    ## Run 38 stress 0.07629236 
    ## Run 39 stress 0.08000467 
    ## Run 40 stress 0.07629238 
    ## Run 41 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001087112  max resid 0.0002534332 
    ## ... Similar to previous best
    ## Run 42 stress 0.07732933 
    ## Run 43 stress 0.08288034 
    ## Run 44 stress 0.07629234 
    ## Run 45 stress 0.08288044 
    ## Run 46 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744921  max resid 0.05105779 
    ## Run 47 stress 0.08000468 
    ## Run 48 stress 0.07378232 
    ## ... Procrustes: rmse 0.01745738  max resid 0.05106721 
    ## Run 49 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744259  max resid 0.05104918 
    ## Run 50 stress 0.07365784 
    ## ... Procrustes: rmse 5.778325e-05  max resid 0.0001362481 
    ## ... Similar to previous best
    ## Run 51 stress 0.07629233 
    ## Run 52 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744873  max resid 0.05105305 
    ## Run 53 stress 0.07365787 
    ## ... Procrustes: rmse 9.591438e-05  max resid 0.0002242868 
    ## ... Similar to previous best
    ## Run 54 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744282  max resid 0.05103357 
    ## Run 55 stress 0.07378229 
    ## ... Procrustes: rmse 0.01743882  max resid 0.05104501 
    ## Run 56 stress 0.07629245 
    ## Run 57 stress 0.08000468 
    ## Run 58 stress 0.08000468 
    ## Run 59 stress 0.0762924 
    ## Run 60 stress 0.07365784 
    ## ... Procrustes: rmse 1.570792e-05  max resid 3.214413e-05 
    ## ... Similar to previous best
    ## Run 61 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744349  max resid 0.05103471 
    ## Run 62 stress 0.08000468 
    ## Run 63 stress 0.08000479 
    ## Run 64 stress 0.07365789 
    ## ... Procrustes: rmse 3.872671e-05  max resid 7.473095e-05 
    ## ... Similar to previous best
    ## Run 65 stress 0.08000471 
    ## Run 66 stress 0.07378232 
    ## ... Procrustes: rmse 0.01743453  max resid 0.05103708 
    ## Run 67 stress 0.07629241 
    ## Run 68 stress 0.07629247 
    ## Run 69 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001244367  max resid 0.0002950518 
    ## ... Similar to previous best
    ## Run 70 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174428  max resid 0.05105371 
    ## Run 71 stress 0.08233768 
    ## Run 72 stress 0.08000479 
    ## Run 73 stress 0.07365784 
    ## ... Procrustes: rmse 3.181729e-05  max resid 7.434474e-05 
    ## ... Similar to previous best
    ## Run 74 stress 0.07629238 
    ## Run 75 stress 0.07629234 
    ## Run 76 stress 0.07629242 
    ## Run 77 stress 0.07629239 
    ## Run 78 stress 0.07365783 
    ## ... Procrustes: rmse 1.241059e-05  max resid 2.492758e-05 
    ## ... Similar to previous best
    ## Run 79 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744331  max resid 0.05103303 
    ## Run 80 stress 0.07365794 
    ## ... Procrustes: rmse 0.0001561526  max resid 0.0003731935 
    ## ... Similar to previous best
    ## Run 81 stress 0.07629233 
    ## Run 82 stress 0.0773293 
    ## Run 83 stress 0.07365784 
    ## ... Procrustes: rmse 1.873393e-05  max resid 3.525786e-05 
    ## ... Similar to previous best
    ## Run 84 stress 0.07629241 
    ## Run 85 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744275  max resid 0.05103734 
    ## Run 86 stress 0.08233775 
    ## Run 87 stress 0.07365783 
    ## ... Procrustes: rmse 7.275414e-06  max resid 1.62581e-05 
    ## ... Similar to previous best
    ## Run 88 stress 0.07629234 
    ## Run 89 stress 0.0800048 
    ## Run 90 stress 0.07629242 
    ## Run 91 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744251  max resid 0.05106221 
    ## Run 92 stress 0.07365784 
    ## ... Procrustes: rmse 2.110189e-05  max resid 3.329333e-05 
    ## ... Similar to previous best
    ## Run 93 stress 0.07629238 
    ## Run 94 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745434  max resid 0.05109556 
    ## Run 95 stress 0.0800048 
    ## Run 96 stress 0.0800047 
    ## Run 97 stress 0.08000469 
    ## Run 98 stress 0.07732939 
    ## Run 99 stress 0.08000469 
    ## Run 100 stress 0.0737823 
    ## ... Procrustes: rmse 0.01743005  max resid 0.05095211 
    ## Run 101 stress 0.07365783 
    ## ... Procrustes: rmse 1.273237e-05  max resid 2.868715e-05 
    ## ... Similar to previous best
    ## Run 102 stress 0.07629234 
    ## Run 103 stress 0.07365784 
    ## ... Procrustes: rmse 4.064943e-05  max resid 9.409481e-05 
    ## ... Similar to previous best
    ## Run 104 stress 0.07629239 
    ## Run 105 stress 0.07629235 
    ## Run 106 stress 0.08233758 
    ## Run 107 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001721814  max resid 0.0004040143 
    ## ... Similar to previous best
    ## Run 108 stress 0.07365785 
    ## ... Procrustes: rmse 5.230781e-05  max resid 0.0001249957 
    ## ... Similar to previous best
    ## Run 109 stress 0.07365791 
    ## ... Procrustes: rmse 4.460609e-05  max resid 6.777737e-05 
    ## ... Similar to previous best
    ## Run 110 stress 0.08000475 
    ## Run 111 stress 0.2612248 
    ## Run 112 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001151704  max resid 0.0002684963 
    ## ... Similar to previous best
    ## Run 113 stress 0.07365786 
    ## ... Procrustes: rmse 3.760051e-05  max resid 8.738887e-05 
    ## ... Similar to previous best
    ## Run 114 stress 0.07629235 
    ## Run 115 stress 0.07629239 
    ## Run 116 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001192608  max resid 0.0002757149 
    ## ... Similar to previous best
    ## Run 117 stress 0.08000473 
    ## Run 118 stress 0.07365794 
    ## ... Procrustes: rmse 4.650791e-05  max resid 7.27549e-05 
    ## ... Similar to previous best
    ## Run 119 stress 0.08288041 
    ## Run 120 stress 0.0773294 
    ## Run 121 stress 0.07365784 
    ## ... Procrustes: rmse 3.190733e-05  max resid 7.314903e-05 
    ## ... Similar to previous best
    ## Run 122 stress 0.07629237 
    ## Run 123 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174429  max resid 0.05102501 
    ## Run 124 stress 0.08233753 
    ## Run 125 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744433  max resid 0.05104262 
    ## Run 126 stress 0.08233754 
    ## Run 127 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744435  max resid 0.05105333 
    ## Run 128 stress 0.07629234 
    ## Run 129 stress 0.0828805 
    ## Run 130 stress 0.08000471 
    ## Run 131 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001325452  max resid 0.0003149138 
    ## ... Similar to previous best
    ## Run 132 stress 0.07629248 
    ## Run 133 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744253  max resid 0.05101687 
    ## Run 134 stress 0.07365784 
    ## ... Procrustes: rmse 5.014111e-05  max resid 0.0001198634 
    ## ... Similar to previous best
    ## Run 135 stress 0.0800047 
    ## Run 136 stress 0.07629238 
    ## Run 137 stress 0.07378229 
    ## ... Procrustes: rmse 0.01741974  max resid 0.05095258 
    ## Run 138 stress 0.07629235 
    ## Run 139 stress 0.07629245 
    ## Run 140 stress 0.07365783 
    ## ... Procrustes: rmse 1.520035e-05  max resid 3.512463e-05 
    ## ... Similar to previous best
    ## Run 141 stress 0.08000469 
    ## Run 142 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744781  max resid 0.05106422 
    ## Run 143 stress 0.07629242 
    ## Run 144 stress 0.07365786 
    ## ... Procrustes: rmse 6.796543e-05  max resid 0.000162454 
    ## ... Similar to previous best
    ## Run 145 stress 0.07629237 
    ## Run 146 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744787  max resid 0.05104249 
    ## Run 147 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744498  max resid 0.05104697 
    ## Run 148 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001484871  max resid 0.0003516231 
    ## ... Similar to previous best
    ## Run 149 stress 0.07629243 
    ## Run 150 stress 0.08233752 
    ## Run 151 stress 0.07365784 
    ## ... Procrustes: rmse 3.759841e-05  max resid 8.479601e-05 
    ## ... Similar to previous best
    ## Run 152 stress 0.08000469 
    ## Run 153 stress 0.07629235 
    ## Run 154 stress 0.07732927 
    ## Run 155 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744613  max resid 0.05105655 
    ## Run 156 stress 0.0762924 
    ## Run 157 stress 0.08000467 
    ## Run 158 stress 0.07365785 
    ## ... Procrustes: rmse 6.907927e-05  max resid 0.0001651099 
    ## ... Similar to previous best
    ## Run 159 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744742  max resid 0.05106214 
    ## Run 160 stress 0.07629234 
    ## Run 161 stress 0.07365785 
    ## ... Procrustes: rmse 7.063352e-05  max resid 0.0001684058 
    ## ... Similar to previous best
    ## Run 162 stress 0.08000472 
    ## Run 163 stress 0.0736579 
    ## ... Procrustes: rmse 0.000145111  max resid 0.0003441414 
    ## ... Similar to previous best
    ## Run 164 stress 0.07732935 
    ## Run 165 stress 0.08000468 
    ## Run 166 stress 0.07365787 
    ## ... Procrustes: rmse 9.852193e-05  max resid 0.0002333599 
    ## ... Similar to previous best
    ## Run 167 stress 0.0762924 
    ## Run 168 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744669  max resid 0.05106971 
    ## Run 169 stress 0.07629241 
    ## Run 170 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001180344  max resid 0.0002785066 
    ## ... Similar to previous best
    ## Run 171 stress 0.07365783 
    ## ... Procrustes: rmse 1.507482e-05  max resid 3.275631e-05 
    ## ... Similar to previous best
    ## Run 172 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001485381  max resid 0.0003499968 
    ## ... Similar to previous best
    ## Run 173 stress 0.08000476 
    ## Run 174 stress 0.07365786 
    ## ... Procrustes: rmse 7.935489e-05  max resid 0.0001896532 
    ## ... Similar to previous best
    ## Run 175 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001145642  max resid 0.0002717561 
    ## ... Similar to previous best
    ## Run 176 stress 0.07365785 
    ## ... Procrustes: rmse 7.203416e-05  max resid 0.0001706884 
    ## ... Similar to previous best
    ## Run 177 stress 0.3407909 
    ## Run 178 stress 0.08233752 
    ## Run 179 stress 0.07365785 
    ## ... Procrustes: rmse 7.617627e-05  max resid 0.0001802531 
    ## ... Similar to previous best
    ## Run 180 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174404  max resid 0.05102482 
    ## Run 181 stress 0.07365787 
    ## ... Procrustes: rmse 9.390593e-05  max resid 0.000222234 
    ## ... Similar to previous best
    ## Run 182 stress 0.08288041 
    ## Run 183 stress 0.07365784 
    ## ... Procrustes: rmse 4.834859e-05  max resid 0.0001135281 
    ## ... Similar to previous best
    ## Run 184 stress 0.07629235 
    ## Run 185 stress 0.07629237 
    ## Run 186 stress 0.2829247 
    ## Run 187 stress 0.07732929 
    ## Run 188 stress 0.08000467 
    ## Run 189 stress 0.07629236 
    ## Run 190 stress 0.07629233 
    ## Run 191 stress 0.08233771 
    ## Run 192 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744254  max resid 0.0510402 
    ## Run 193 stress 0.08000471 
    ## Run 194 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744512  max resid 0.05104806 
    ## Run 195 stress 0.07365783 
    ## ... Procrustes: rmse 6.550094e-06  max resid 1.289413e-05 
    ## ... Similar to previous best
    ## Run 196 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744743  max resid 0.051053 
    ## Run 197 stress 0.07365784 
    ## ... Procrustes: rmse 3.731913e-05  max resid 8.925129e-05 
    ## ... Similar to previous best
    ## Run 198 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001237724  max resid 0.0002948863 
    ## ... Similar to previous best
    ## Run 199 stress 0.07629237 
    ## Run 200 stress 0.07629237 
    ## Run 201 stress 0.07629236 
    ## Run 202 stress 0.08233752 
    ## Run 203 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744065  max resid 0.05102596 
    ## Run 204 stress 0.07629239 
    ## Run 205 stress 0.07629242 
    ## Run 206 stress 0.0762924 
    ## Run 207 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745976  max resid 0.05112262 
    ## Run 208 stress 0.07365786 
    ## ... Procrustes: rmse 8.608075e-05  max resid 0.0001984912 
    ## ... Similar to previous best
    ## Run 209 stress 0.08000467 
    ## Run 210 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744184  max resid 0.05103129 
    ## Run 211 stress 0.07629239 
    ## Run 212 stress 0.08233753 
    ## Run 213 stress 0.07365783 
    ## ... Procrustes: rmse 1.285393e-05  max resid 2.966714e-05 
    ## ... Similar to previous best
    ## Run 214 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744306  max resid 0.0510337 
    ## Run 215 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745484  max resid 0.05110408 
    ## Run 216 stress 0.08000469 
    ## Run 217 stress 0.08000474 
    ## Run 218 stress 0.07629237 
    ## Run 219 stress 0.07629239 
    ## Run 220 stress 0.08000474 
    ## Run 221 stress 0.08000469 
    ## Run 222 stress 0.07365787 
    ## ... Procrustes: rmse 9.016765e-05  max resid 0.0002156814 
    ## ... Similar to previous best
    ## Run 223 stress 0.349572 
    ## Run 224 stress 0.07378232 
    ## ... Procrustes: rmse 0.01743582  max resid 0.05104114 
    ## Run 225 stress 0.07365784 
    ## ... Procrustes: rmse 2.1014e-05  max resid 4.364949e-05 
    ## ... Similar to previous best
    ## Run 226 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744332  max resid 0.05103193 
    ## Run 227 stress 0.08000468 
    ## Run 228 stress 0.07378229 
    ## ... Procrustes: rmse 0.0174407  max resid 0.05105035 
    ## Run 229 stress 0.07629236 
    ## Run 230 stress 0.0800047 
    ## Run 231 stress 0.3582285 
    ## Run 232 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744618  max resid 0.05106052 
    ## Run 233 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 1.44101e-05  max resid 3.210968e-05 
    ## ... Similar to previous best
    ## Run 234 stress 0.08233752 
    ## Run 235 stress 0.07365786 
    ## ... Procrustes: rmse 9.870077e-05  max resid 0.0002332712 
    ## ... Similar to previous best
    ## Run 236 stress 0.07365784 
    ## ... Procrustes: rmse 6.714262e-05  max resid 0.000156394 
    ## ... Similar to previous best
    ## Run 237 stress 0.07365795 
    ## ... Procrustes: rmse 0.000186237  max resid 0.0004357293 
    ## ... Similar to previous best
    ## Run 238 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744986  max resid 0.05106238 
    ## Run 239 stress 0.07629244 
    ## Run 240 stress 0.07365783 
    ## ... Procrustes: rmse 3.690795e-05  max resid 8.643456e-05 
    ## ... Similar to previous best
    ## Run 241 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001027814  max resid 0.0002353417 
    ## ... Similar to previous best
    ## Run 242 stress 0.07365784 
    ## ... Procrustes: rmse 4.5232e-05  max resid 0.0001069533 
    ## ... Similar to previous best
    ## Run 243 stress 0.07365783 
    ## ... Procrustes: rmse 2.808255e-05  max resid 6.514201e-05 
    ## ... Similar to previous best
    ## Run 244 stress 0.08000474 
    ## Run 245 stress 0.2754056 
    ## Run 246 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744703  max resid 0.05104384 
    ## Run 247 stress 0.08000471 
    ## Run 248 stress 0.07629236 
    ## Run 249 stress 0.07629234 
    ## Run 250 stress 0.07629242 
    ## Run 251 stress 0.07365783 
    ## ... Procrustes: rmse 5.568335e-06  max resid 1.018337e-05 
    ## ... Similar to previous best
    ## Run 252 stress 0.07365783 
    ## ... Procrustes: rmse 7.670252e-06  max resid 1.547031e-05 
    ## ... Similar to previous best
    ## Run 253 stress 0.08288041 
    ## Run 254 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745063  max resid 0.05107311 
    ## Run 255 stress 0.08288043 
    ## Run 256 stress 0.07629233 
    ## Run 257 stress 0.07629244 
    ## Run 258 stress 0.08288051 
    ## Run 259 stress 0.07365784 
    ## ... Procrustes: rmse 5.481712e-05  max resid 0.0001297983 
    ## ... Similar to previous best
    ## Run 260 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744939  max resid 0.05105102 
    ## Run 261 stress 0.08000467 
    ## Run 262 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745216  max resid 0.05106123 
    ## Run 263 stress 0.07732931 
    ## Run 264 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744289  max resid 0.05106018 
    ## Run 265 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 4.642708e-06  max resid 8.479016e-06 
    ## ... Similar to previous best
    ## Run 266 stress 0.07378232 
    ## ... Procrustes: rmse 0.0174585  max resid 0.05112166 
    ## Run 267 stress 0.08288057 
    ## Run 268 stress 0.07629234 
    ## Run 269 stress 0.07629233 
    ## Run 270 stress 0.07365784 
    ## ... Procrustes: rmse 2.903442e-05  max resid 6.895597e-05 
    ## ... Similar to previous best
    ## Run 271 stress 0.07629236 
    ## Run 272 stress 0.08233762 
    ## Run 273 stress 0.07629235 
    ## Run 274 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744746  max resid 0.05105615 
    ## Run 275 stress 0.07732929 
    ## Run 276 stress 0.07629236 
    ## Run 277 stress 0.07629237 
    ## Run 278 stress 0.07365784 
    ## ... Procrustes: rmse 3.795691e-05  max resid 8.793167e-05 
    ## ... Similar to previous best
    ## Run 279 stress 0.2763959 
    ## Run 280 stress 0.08000472 
    ## Run 281 stress 0.07365784 
    ## ... Procrustes: rmse 5.329101e-05  max resid 0.0001250191 
    ## ... Similar to previous best
    ## Run 282 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744566  max resid 0.05104449 
    ## Run 283 stress 0.08000468 
    ## Run 284 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744667  max resid 0.0510473 
    ## Run 285 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744852  max resid 0.05104214 
    ## Run 286 stress 0.08233759 
    ## Run 287 stress 0.08288039 
    ## Run 288 stress 0.07629234 
    ## Run 289 stress 0.07629238 
    ## Run 290 stress 0.08288041 
    ## Run 291 stress 0.07629236 
    ## Run 292 stress 0.07629234 
    ## Run 293 stress 0.0800047 
    ## Run 294 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001505056  max resid 0.0003557934 
    ## ... Similar to previous best
    ## Run 295 stress 0.07629237 
    ## Run 296 stress 0.07629234 
    ## Run 297 stress 0.08000469 
    ## Run 298 stress 0.08288054 
    ## Run 299 stress 0.3578478 
    ## Run 300 stress 0.07378229 
    ## ... Procrustes: rmse 0.01743739  max resid 0.05098694 
    ## Run 301 stress 0.08000468 
    ## Run 302 stress 0.07629241 
    ## Run 303 stress 0.07365783 
    ## ... Procrustes: rmse 3.761079e-05  max resid 8.749342e-05 
    ## ... Similar to previous best
    ## Run 304 stress 0.07378231 
    ## ... Procrustes: rmse 0.01744731  max resid 0.05104729 
    ## Run 305 stress 0.07365784 
    ## ... Procrustes: rmse 5.420798e-05  max resid 0.0001276701 
    ## ... Similar to previous best
    ## Run 306 stress 0.07732932 
    ## Run 307 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744671  max resid 0.05104131 
    ## Run 308 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744795  max resid 0.05105723 
    ## Run 309 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744184  max resid 0.05105208 
    ## Run 310 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744887  max resid 0.05106245 
    ## Run 311 stress 0.07629237 
    ## Run 312 stress 0.3578489 
    ## Run 313 stress 0.07365785 
    ## ... Procrustes: rmse 6.541036e-05  max resid 0.0001509008 
    ## ... Similar to previous best
    ## Run 314 stress 0.07629235 
    ## Run 315 stress 0.07629239 
    ## Run 316 stress 0.07378234 
    ## ... Procrustes: rmse 0.01745794  max resid 0.05113178 
    ## Run 317 stress 0.07365784 
    ## ... Procrustes: rmse 6.427951e-05  max resid 0.0001526469 
    ## ... Similar to previous best
    ## Run 318 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744567  max resid 0.05104393 
    ## Run 319 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001404712  max resid 0.000332026 
    ## ... Similar to previous best
    ## Run 320 stress 0.07629243 
    ## Run 321 stress 0.07629239 
    ## Run 322 stress 0.07732936 
    ## Run 323 stress 0.08233752 
    ## Run 324 stress 0.07629237 
    ## Run 325 stress 0.08000474 
    ## Run 326 stress 0.07629234 
    ## Run 327 stress 0.07629236 
    ## Run 328 stress 0.07629234 
    ## Run 329 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744394  max resid 0.05103539 
    ## Run 330 stress 0.08288037 
    ## Run 331 stress 0.08000469 
    ## Run 332 stress 0.07629234 
    ## Run 333 stress 0.07365784 
    ## ... Procrustes: rmse 6.281208e-05  max resid 0.0001483441 
    ## ... Similar to previous best
    ## Run 334 stress 0.07629233 
    ## Run 335 stress 0.07365785 
    ## ... Procrustes: rmse 8.292026e-05  max resid 0.0001933315 
    ## ... Similar to previous best
    ## Run 336 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744588  max resid 0.05104182 
    ## Run 337 stress 0.07629235 
    ## Run 338 stress 0.07378232 
    ## ... Procrustes: rmse 0.01743898  max resid 0.05105139 
    ## Run 339 stress 0.07732934 
    ## Run 340 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744615  max resid 0.05103571 
    ## Run 341 stress 0.08000467 
    ## Run 342 stress 0.07732935 
    ## Run 343 stress 0.08000474 
    ## Run 344 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745841  max resid 0.0511124 
    ## Run 345 stress 0.08000474 
    ## Run 346 stress 0.07629233 
    ## Run 347 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001194054  max resid 0.0002834737 
    ## ... Similar to previous best
    ## Run 348 stress 0.07365794 
    ## ... Procrustes: rmse 0.0001702534  max resid 0.0003952368 
    ## ... Similar to previous best
    ## Run 349 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 6.70418e-06  max resid 1.263117e-05 
    ## ... Similar to previous best
    ## Run 350 stress 0.07365785 
    ## ... Procrustes: rmse 5.577951e-05  max resid 0.0001278889 
    ## ... Similar to previous best
    ## Run 351 stress 0.07365783 
    ## ... Procrustes: rmse 6.340622e-06  max resid 1.45654e-05 
    ## ... Similar to previous best
    ## Run 352 stress 0.07365783 
    ## ... Procrustes: rmse 3.460309e-05  max resid 7.912173e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.07365786 
    ## ... Procrustes: rmse 8.727226e-05  max resid 0.0002061646 
    ## ... Similar to previous best
    ## Run 354 stress 0.3582284 
    ## Run 355 stress 0.07629242 
    ## Run 356 stress 0.07629247 
    ## Run 357 stress 0.08000475 
    ## Run 358 stress 0.07629234 
    ## Run 359 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744331  max resid 0.05103458 
    ## Run 360 stress 0.07365786 
    ## ... Procrustes: rmse 7.917426e-05  max resid 0.0001873229 
    ## ... Similar to previous best
    ## Run 361 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001418844  max resid 0.0003379328 
    ## ... Similar to previous best
    ## Run 362 stress 0.07365786 
    ## ... Procrustes: rmse 8.612436e-05  max resid 0.0002051414 
    ## ... Similar to previous best
    ## Run 363 stress 0.07629235 
    ## Run 364 stress 0.07629243 
    ## Run 365 stress 0.07629244 
    ## Run 366 stress 0.07732936 
    ## Run 367 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001300785  max resid 0.0003132721 
    ## ... Similar to previous best
    ## Run 368 stress 0.07629247 
    ## Run 369 stress 0.07378232 
    ## ... Procrustes: rmse 0.01744981  max resid 0.05102743 
    ## Run 370 stress 0.07365783 
    ## ... Procrustes: rmse 1.135305e-05  max resid 2.34135e-05 
    ## ... Similar to previous best
    ## Run 371 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744289  max resid 0.05104829 
    ## Run 372 stress 0.08000469 
    ## Run 373 stress 0.07629245 
    ## Run 374 stress 0.08233772 
    ## Run 375 stress 0.07629233 
    ## Run 376 stress 0.08000467 
    ## Run 377 stress 0.07365784 
    ## ... Procrustes: rmse 3.215788e-05  max resid 7.204257e-05 
    ## ... Similar to previous best
    ## Run 378 stress 0.08233774 
    ## Run 379 stress 0.07365783 
    ## ... Procrustes: rmse 2.453539e-05  max resid 6.049106e-05 
    ## ... Similar to previous best
    ## Run 380 stress 0.0762924 
    ## Run 381 stress 0.07365784 
    ## ... Procrustes: rmse 3.792003e-05  max resid 8.693712e-05 
    ## ... Similar to previous best
    ## Run 382 stress 0.07365783 
    ## ... Procrustes: rmse 2.203014e-05  max resid 5.009242e-05 
    ## ... Similar to previous best
    ## Run 383 stress 0.07629238 
    ## Run 384 stress 0.07365786 
    ## ... Procrustes: rmse 9.003273e-05  max resid 0.000215519 
    ## ... Similar to previous best
    ## Run 385 stress 0.07378231 
    ## ... Procrustes: rmse 0.01743913  max resid 0.05105103 
    ## Run 386 stress 0.07732927 
    ## Run 387 stress 0.07365787 
    ## ... Procrustes: rmse 9.969249e-05  max resid 0.0002376989 
    ## ... Similar to previous best
    ## Run 388 stress 0.07365785 
    ## ... Procrustes: rmse 6.332539e-05  max resid 0.0001457117 
    ## ... Similar to previous best
    ## Run 389 stress 0.07378235 
    ## ... Procrustes: rmse 0.0174432  max resid 0.0510296 
    ## Run 390 stress 0.07732929 
    ## Run 391 stress 0.07732928 
    ## Run 392 stress 0.07629235 
    ## Run 393 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743699  max resid 0.05100947 
    ## Run 394 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744447  max resid 0.05104646 
    ## Run 395 stress 0.08000476 
    ## Run 396 stress 0.07365783 
    ## ... Procrustes: rmse 8.549629e-06  max resid 1.677713e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.07732935 
    ## Run 398 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744271  max resid 0.05102285 
    ## Run 399 stress 0.07365786 
    ## ... Procrustes: rmse 8.439818e-05  max resid 0.0002022216 
    ## ... Similar to previous best
    ## Run 400 stress 0.07365784 
    ## ... Procrustes: rmse 6.233004e-05  max resid 0.0001458607 
    ## ... Similar to previous best
    ## Run 401 stress 0.07629252 
    ## Run 402 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745774  max resid 0.05108736 
    ## Run 403 stress 0.07629242 
    ## Run 404 stress 0.08000472 
    ## Run 405 stress 0.08000468 
    ## Run 406 stress 0.07732926 
    ## Run 407 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001417025  max resid 0.0003401614 
    ## ... Similar to previous best
    ## Run 408 stress 0.08233781 
    ## Run 409 stress 0.08000469 
    ## Run 410 stress 0.07365784 
    ## ... Procrustes: rmse 4.369629e-05  max resid 0.0001060968 
    ## ... Similar to previous best
    ## Run 411 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744815  max resid 0.05106958 
    ## Run 412 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001117094  max resid 0.0002689531 
    ## ... Similar to previous best
    ## Run 413 stress 0.08000467 
    ## Run 414 stress 0.08000468 
    ## Run 415 stress 0.08000474 
    ## Run 416 stress 0.08233766 
    ## Run 417 stress 0.08000469 
    ## Run 418 stress 0.07629233 
    ## Run 419 stress 0.0823377 
    ## Run 420 stress 0.07629234 
    ## Run 421 stress 0.07365783 
    ## ... Procrustes: rmse 4.760889e-06  max resid 1.195214e-05 
    ## ... Similar to previous best
    ## Run 422 stress 0.0800047 
    ## Run 423 stress 0.07365787 
    ## ... Procrustes: rmse 9.95104e-05  max resid 0.0002389075 
    ## ... Similar to previous best
    ## Run 424 stress 0.08000475 
    ## Run 425 stress 0.07629243 
    ## Run 426 stress 0.3495546 
    ## Run 427 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744083  max resid 0.05104695 
    ## Run 428 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745731  max resid 0.05111559 
    ## Run 429 stress 0.07732946 
    ## Run 430 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001693646  max resid 0.0004003113 
    ## ... Similar to previous best
    ## Run 431 stress 0.07629243 
    ## Run 432 stress 0.08000468 
    ## Run 433 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001535137  max resid 0.0003679579 
    ## ... Similar to previous best
    ## Run 434 stress 0.07629239 
    ## Run 435 stress 0.0800047 
    ## Run 436 stress 0.08233759 
    ## Run 437 stress 0.07629234 
    ## Run 438 stress 0.07629238 
    ## Run 439 stress 0.07629246 
    ## Run 440 stress 0.07365784 
    ## ... Procrustes: rmse 3.06893e-05  max resid 7.573051e-05 
    ## ... Similar to previous best
    ## Run 441 stress 0.08000468 
    ## Run 442 stress 0.0773294 
    ## Run 443 stress 0.07629238 
    ## Run 444 stress 0.08000467 
    ## Run 445 stress 0.07629242 
    ## Run 446 stress 0.07365783 
    ## ... Procrustes: rmse 1.695012e-05  max resid 4.027012e-05 
    ## ... Similar to previous best
    ## Run 447 stress 0.07629245 
    ## Run 448 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745512  max resid 0.05109978 
    ## Run 449 stress 0.08000468 
    ## Run 450 stress 0.08000469 
    ## Run 451 stress 0.07378232 
    ## ... Procrustes: rmse 0.01743342  max resid 0.05103518 
    ## Run 452 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744503  max resid 0.05104757 
    ## Run 453 stress 0.07365785 
    ## ... Procrustes: rmse 7.115099e-05  max resid 0.0001663709 
    ## ... Similar to previous best
    ## Run 454 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745132  max resid 0.05107682 
    ## Run 455 stress 0.07629236 
    ## Run 456 stress 0.08000469 
    ## Run 457 stress 0.08000472 
    ## Run 458 stress 0.07365786 
    ## ... Procrustes: rmse 9.067159e-05  max resid 0.0002164657 
    ## ... Similar to previous best
    ## Run 459 stress 0.0800047 
    ## Run 460 stress 0.08000476 
    ## Run 461 stress 0.07629242 
    ## Run 462 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744544  max resid 0.05104459 
    ## Run 463 stress 0.0737823 
    ## ... Procrustes: rmse 0.01745894  max resid 0.05112297 
    ## Run 464 stress 0.08288055 
    ## Run 465 stress 0.07629237 
    ## Run 466 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743988  max resid 0.05102362 
    ## Run 467 stress 0.07365785 
    ## ... Procrustes: rmse 7.041209e-05  max resid 0.0001603291 
    ## ... Similar to previous best
    ## Run 468 stress 0.08000483 
    ## Run 469 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744291  max resid 0.05103217 
    ## Run 470 stress 0.07629249 
    ## Run 471 stress 0.08000468 
    ## Run 472 stress 0.07629247 
    ## Run 473 stress 0.07365783 
    ## ... Procrustes: rmse 7.177866e-06  max resid 1.579151e-05 
    ## ... Similar to previous best
    ## Run 474 stress 0.3578479 
    ## Run 475 stress 0.07629238 
    ## Run 476 stress 0.07365783 
    ## ... Procrustes: rmse 5.22689e-06  max resid 1.030471e-05 
    ## ... Similar to previous best
    ## Run 477 stress 0.07365783 
    ## ... Procrustes: rmse 1.492249e-05  max resid 3.10197e-05 
    ## ... Similar to previous best
    ## Run 478 stress 0.07629248 
    ## Run 479 stress 0.07732945 
    ## Run 480 stress 0.08000467 
    ## Run 481 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001024586  max resid 0.000242706 
    ## ... Similar to previous best
    ## Run 482 stress 0.07629244 
    ## Run 483 stress 0.07365783 
    ## ... Procrustes: rmse 2.203048e-05  max resid 5.43653e-05 
    ## ... Similar to previous best
    ## Run 484 stress 0.07629241 
    ## Run 485 stress 0.0762924 
    ## Run 486 stress 0.08000474 
    ## Run 487 stress 0.07365785 
    ## ... Procrustes: rmse 7.747758e-05  max resid 0.0001860031 
    ## ... Similar to previous best
    ## Run 488 stress 0.08000467 
    ## Run 489 stress 0.07365783 
    ## ... Procrustes: rmse 2.825454e-05  max resid 6.500697e-05 
    ## ... Similar to previous best
    ## Run 490 stress 0.08000472 
    ## Run 491 stress 0.07629235 
    ## Run 492 stress 0.07629235 
    ## Run 493 stress 0.07365785 
    ## ... Procrustes: rmse 6.264216e-05  max resid 0.0001523141 
    ## ... Similar to previous best
    ## Run 494 stress 0.07365784 
    ## ... Procrustes: rmse 3.145497e-05  max resid 7.738602e-05 
    ## ... Similar to previous best
    ## Run 495 stress 0.07629237 
    ## Run 496 stress 0.0800047 
    ## Run 497 stress 0.07365785 
    ## ... Procrustes: rmse 6.184987e-05  max resid 0.0001444328 
    ## ... Similar to previous best
    ## Run 498 stress 0.07378231 
    ## ... Procrustes: rmse 0.01746016  max resid 0.0511365 
    ## Run 499 stress 0.07365783 
    ## ... Procrustes: rmse 2.660446e-05  max resid 6.569419e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.3495701 
    ## *** Best solution repeated 43 times

``` r
# Stratified lakes and ocean sites
SD_beta_geo_SO_NMDS <- metaMDS(SD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07970526 
    ## Run 1 stress 0.08340296 
    ## Run 2 stress 0.06978197 
    ## ... New best solution
    ## ... Procrustes: rmse 0.09973819  max resid 0.2618589 
    ## Run 3 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01327206  max resid 0.03360458 
    ## Run 4 stress 0.08448454 
    ## Run 5 stress 0.07428316 
    ## Run 6 stress 0.07970526 
    ## Run 7 stress 0.07250812 
    ## Run 8 stress 0.08448452 
    ## Run 9 stress 0.06942776 
    ## ... Procrustes: rmse 2.382562e-05  max resid 6.285738e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.07970525 
    ## Run 11 stress 0.07250812 
    ## Run 12 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322559  max resid 0.03327007 
    ## Run 13 stress 0.08448443 
    ## Run 14 stress 0.07970525 
    ## Run 15 stress 0.07250812 
    ## Run 16 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326008  max resid 0.03334163 
    ## Run 17 stress 0.0742832 
    ## Run 18 stress 0.07428314 
    ## Run 19 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 9.901914e-06  max resid 2.732135e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.07428315 
    ## Run 21 stress 0.06942776 
    ## ... Procrustes: rmse 6.181147e-06  max resid 1.83011e-05 
    ## ... Similar to previous best
    ## Run 22 stress 0.07970525 
    ## Run 23 stress 0.06942778 
    ## ... Procrustes: rmse 7.573131e-05  max resid 0.000194448 
    ## ... Similar to previous best
    ## Run 24 stress 0.07428315 
    ## Run 25 stress 0.07970525 
    ## Run 26 stress 0.07970526 
    ## Run 27 stress 0.06942777 
    ## ... Procrustes: rmse 4.447555e-05  max resid 0.0001132618 
    ## ... Similar to previous best
    ## Run 28 stress 0.07970525 
    ## Run 29 stress 0.06978202 
    ## ... Procrustes: rmse 0.01328534  max resid 0.03339393 
    ## Run 30 stress 0.06942776 
    ## ... Procrustes: rmse 3.110085e-05  max resid 8.03547e-05 
    ## ... Similar to previous best
    ## Run 31 stress 0.07970526 
    ## Run 32 stress 0.06942777 
    ## ... Procrustes: rmse 4.324978e-05  max resid 0.0001143909 
    ## ... Similar to previous best
    ## Run 33 stress 0.08448438 
    ## Run 34 stress 0.07844926 
    ## Run 35 stress 0.08340295 
    ## Run 36 stress 0.06978199 
    ## ... Procrustes: rmse 0.01328034  max resid 0.03337741 
    ## Run 37 stress 0.07970525 
    ## Run 38 stress 0.0784496 
    ## Run 39 stress 0.08448437 
    ## Run 40 stress 0.08340287 
    ## Run 41 stress 0.07428313 
    ## Run 42 stress 0.08448441 
    ## Run 43 stress 0.07428316 
    ## Run 44 stress 0.07428313 
    ## Run 45 stress 0.08448447 
    ## Run 46 stress 0.08340287 
    ## Run 47 stress 0.07428313 
    ## Run 48 stress 0.07428319 
    ## Run 49 stress 0.07970525 
    ## Run 50 stress 0.07250812 
    ## Run 51 stress 0.07844952 
    ## Run 52 stress 0.06942776 
    ## ... Procrustes: rmse 1.83125e-05  max resid 4.713165e-05 
    ## ... Similar to previous best
    ## Run 53 stress 0.07428312 
    ## Run 54 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 6.039994e-06  max resid 1.555395e-05 
    ## ... Similar to previous best
    ## Run 55 stress 0.08340294 
    ## Run 56 stress 0.08340287 
    ## Run 57 stress 0.06942778 
    ## ... Procrustes: rmse 5.72981e-05  max resid 0.0001479011 
    ## ... Similar to previous best
    ## Run 58 stress 0.07428313 
    ## Run 59 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326433  max resid 0.03334649 
    ## Run 60 stress 0.07970525 
    ## Run 61 stress 0.07970525 
    ## Run 62 stress 0.06942777 
    ## ... Procrustes: rmse 3.613207e-05  max resid 9.328035e-05 
    ## ... Similar to previous best
    ## Run 63 stress 0.08340293 
    ## Run 64 stress 0.06942776 
    ## ... Procrustes: rmse 1.734603e-05  max resid 4.494677e-05 
    ## ... Similar to previous best
    ## Run 65 stress 0.07970526 
    ## Run 66 stress 0.07428317 
    ## Run 67 stress 0.06942776 
    ## ... Procrustes: rmse 7.164762e-06  max resid 1.814341e-05 
    ## ... Similar to previous best
    ## Run 68 stress 0.07970525 
    ## Run 69 stress 0.07250814 
    ## Run 70 stress 0.08448436 
    ## Run 71 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326921  max resid 0.033355 
    ## Run 72 stress 0.06942776 
    ## ... Procrustes: rmse 9.910228e-06  max resid 2.488447e-05 
    ## ... Similar to previous best
    ## Run 73 stress 0.07428313 
    ## Run 74 stress 0.07428314 
    ## Run 75 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325211  max resid 0.03332065 
    ## Run 76 stress 0.08340287 
    ## Run 77 stress 0.07844925 
    ## Run 78 stress 0.08340288 
    ## Run 79 stress 0.07428314 
    ## Run 80 stress 0.08448449 
    ## Run 81 stress 0.08340288 
    ## Run 82 stress 0.06942776 
    ## ... Procrustes: rmse 1.512451e-05  max resid 3.875315e-05 
    ## ... Similar to previous best
    ## Run 83 stress 0.07250812 
    ## Run 84 stress 0.07844917 
    ## Run 85 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323964  max resid 0.03329623 
    ## Run 86 stress 0.07970525 
    ## Run 87 stress 0.08448441 
    ## Run 88 stress 0.07428313 
    ## Run 89 stress 0.07970525 
    ## Run 90 stress 0.07428313 
    ## Run 91 stress 0.07428315 
    ## Run 92 stress 0.07250815 
    ## Run 93 stress 0.07970526 
    ## Run 94 stress 0.07250813 
    ## Run 95 stress 0.06942776 
    ## ... Procrustes: rmse 2.38561e-05  max resid 6.1393e-05 
    ## ... Similar to previous best
    ## Run 96 stress 0.08448447 
    ## Run 97 stress 0.07428316 
    ## Run 98 stress 0.07428319 
    ## Run 99 stress 0.07250812 
    ## Run 100 stress 0.07428313 
    ## Run 101 stress 0.07970526 
    ## Run 102 stress 0.08340289 
    ## Run 103 stress 0.07250816 
    ## Run 104 stress 0.07428313 
    ## Run 105 stress 0.06942778 
    ## ... Procrustes: rmse 6.505836e-05  max resid 0.0001672878 
    ## ... Similar to previous best
    ## Run 106 stress 0.07250814 
    ## Run 107 stress 0.07428314 
    ## Run 108 stress 0.07970525 
    ## Run 109 stress 0.06942776 
    ## ... Procrustes: rmse 1.676937e-05  max resid 4.285743e-05 
    ## ... Similar to previous best
    ## Run 110 stress 0.07428313 
    ## Run 111 stress 0.07428314 
    ## Run 112 stress 0.07250812 
    ## Run 113 stress 0.06942776 
    ## ... Procrustes: rmse 3.102506e-05  max resid 7.973731e-05 
    ## ... Similar to previous best
    ## Run 114 stress 0.07250813 
    ## Run 115 stress 0.07250815 
    ## Run 116 stress 0.08340286 
    ## Run 117 stress 0.07428313 
    ## Run 118 stress 0.08340289 
    ## Run 119 stress 0.07250812 
    ## Run 120 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 6.318946e-06  max resid 1.594576e-05 
    ## ... Similar to previous best
    ## Run 121 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324891  max resid 0.03331928 
    ## Run 122 stress 0.06942776 
    ## ... Procrustes: rmse 6.622405e-06  max resid 1.595537e-05 
    ## ... Similar to previous best
    ## Run 123 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325328  max resid 0.03332713 
    ## Run 124 stress 0.07970526 
    ## Run 125 stress 0.069782 
    ## ... Procrustes: rmse 0.01328116  max resid 0.03338205 
    ## Run 126 stress 0.07428314 
    ## Run 127 stress 0.07428317 
    ## Run 128 stress 0.07428317 
    ## Run 129 stress 0.07250812 
    ## Run 130 stress 0.07428316 
    ## Run 131 stress 0.06978194 
    ## ... Procrustes: rmse 0.01324363  max resid 0.03331028 
    ## Run 132 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324523  max resid 0.03331064 
    ## Run 133 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323751  max resid 0.03329467 
    ## Run 134 stress 0.07428315 
    ## Run 135 stress 0.08340286 
    ## Run 136 stress 0.06942776 
    ## ... Procrustes: rmse 8.703485e-06  max resid 2.368016e-05 
    ## ... Similar to previous best
    ## Run 137 stress 0.06978199 
    ## ... Procrustes: rmse 0.01327898  max resid 0.03337681 
    ## Run 138 stress 0.07428316 
    ## Run 139 stress 0.07250815 
    ## Run 140 stress 0.07250812 
    ## Run 141 stress 0.07428312 
    ## Run 142 stress 0.07428313 
    ## Run 143 stress 0.07428315 
    ## Run 144 stress 0.07428313 
    ## Run 145 stress 0.07970525 
    ## Run 146 stress 0.06942777 
    ## ... Procrustes: rmse 3.417725e-05  max resid 8.817812e-05 
    ## ... Similar to previous best
    ## Run 147 stress 0.06942778 
    ## ... Procrustes: rmse 5.115797e-05  max resid 0.0001316019 
    ## ... Similar to previous best
    ## Run 148 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323127  max resid 0.03328032 
    ## Run 149 stress 0.07970525 
    ## Run 150 stress 0.08340286 
    ## Run 151 stress 0.07428313 
    ## Run 152 stress 0.06978199 
    ## ... Procrustes: rmse 0.01319525  max resid 0.03320935 
    ## Run 153 stress 0.07428315 
    ## Run 154 stress 0.08448447 
    ## Run 155 stress 0.06978192 
    ## ... Procrustes: rmse 0.01323153  max resid 0.03328355 
    ## Run 156 stress 0.08448446 
    ## Run 157 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319291  max resid 0.03319966 
    ## Run 158 stress 0.06942776 
    ## ... Procrustes: rmse 9.539237e-07  max resid 2.518352e-06 
    ## ... Similar to previous best
    ## Run 159 stress 0.06978194 
    ## ... Procrustes: rmse 0.01319143  max resid 0.03319192 
    ## Run 160 stress 0.07428315 
    ## Run 161 stress 0.07970525 
    ## Run 162 stress 0.06978201 
    ## ... Procrustes: rmse 0.01328197  max resid 0.03338114 
    ## Run 163 stress 0.06978192 
    ## ... Procrustes: rmse 0.01320106  max resid 0.03321869 
    ## Run 164 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323972  max resid 0.03329828 
    ## Run 165 stress 0.06942778 
    ## ... Procrustes: rmse 2.424833e-05  max resid 6.783554e-05 
    ## ... Similar to previous best
    ## Run 166 stress 0.07970525 
    ## Run 167 stress 0.07250812 
    ## Run 168 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321182  max resid 0.03323878 
    ## Run 169 stress 0.08448436 
    ## Run 170 stress 0.08340291 
    ## Run 171 stress 0.07428314 
    ## Run 172 stress 0.07428314 
    ## Run 173 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324956  max resid 0.03331605 
    ## Run 174 stress 0.06942776 
    ## ... Procrustes: rmse 2.473809e-05  max resid 6.433903e-05 
    ## ... Similar to previous best
    ## Run 175 stress 0.08340293 
    ## Run 176 stress 0.08340287 
    ## Run 177 stress 0.07970525 
    ## Run 178 stress 0.07970525 
    ## Run 179 stress 0.07250813 
    ## Run 180 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322711  max resid 0.03327143 
    ## Run 181 stress 0.07970526 
    ## Run 182 stress 0.07428315 
    ## Run 183 stress 0.07844922 
    ## Run 184 stress 0.08340287 
    ## Run 185 stress 0.07428315 
    ## Run 186 stress 0.08340287 
    ## Run 187 stress 0.07428315 
    ## Run 188 stress 0.08340288 
    ## Run 189 stress 0.0742832 
    ## Run 190 stress 0.08448442 
    ## Run 191 stress 0.0834029 
    ## Run 192 stress 0.08340298 
    ## Run 193 stress 0.07970526 
    ## Run 194 stress 0.07250812 
    ## Run 195 stress 0.08340293 
    ## Run 196 stress 0.07428313 
    ## Run 197 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323949  max resid 0.03329852 
    ## Run 198 stress 0.07428316 
    ## Run 199 stress 0.08340288 
    ## Run 200 stress 0.08340287 
    ## Run 201 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324497  max resid 0.03330947 
    ## Run 202 stress 0.07428322 
    ## Run 203 stress 0.06942776 
    ## ... Procrustes: rmse 1.03625e-05  max resid 2.754822e-05 
    ## ... Similar to previous best
    ## Run 204 stress 0.07970525 
    ## Run 205 stress 0.06942776 
    ## ... Procrustes: rmse 1.120408e-06  max resid 2.499331e-06 
    ## ... Similar to previous best
    ## Run 206 stress 0.08340287 
    ## Run 207 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327472  max resid 0.03336852 
    ## Run 208 stress 0.08448437 
    ## Run 209 stress 0.08340288 
    ## Run 210 stress 0.08340287 
    ## Run 211 stress 0.08340288 
    ## Run 212 stress 0.07250814 
    ## Run 213 stress 0.08340291 
    ## Run 214 stress 0.07428313 
    ## Run 215 stress 0.07844962 
    ## Run 216 stress 0.06978197 
    ## ... Procrustes: rmse 0.01317038  max resid 0.03314872 
    ## Run 217 stress 0.07970526 
    ## Run 218 stress 0.07428313 
    ## Run 219 stress 0.08340289 
    ## Run 220 stress 0.08340287 
    ## Run 221 stress 0.0834029 
    ## Run 222 stress 0.07970525 
    ## Run 223 stress 0.07250812 
    ## Run 224 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326127  max resid 0.03334263 
    ## Run 225 stress 0.06978191 
    ## ... Procrustes: rmse 0.01320143  max resid 0.03321741 
    ## Run 226 stress 0.07428313 
    ## Run 227 stress 0.08448452 
    ## Run 228 stress 0.07250812 
    ## Run 229 stress 0.08448435 
    ## Run 230 stress 0.07428313 
    ## Run 231 stress 0.08448438 
    ## Run 232 stress 0.08340289 
    ## Run 233 stress 0.07428314 
    ## Run 234 stress 0.0834029 
    ## Run 235 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324457  max resid 0.0333084 
    ## Run 236 stress 0.07428314 
    ## Run 237 stress 0.07250812 
    ## Run 238 stress 0.06942776 
    ## ... Procrustes: rmse 7.201613e-06  max resid 1.705305e-05 
    ## ... Similar to previous best
    ## Run 239 stress 0.07250816 
    ## Run 240 stress 0.07428313 
    ## Run 241 stress 0.07428312 
    ## Run 242 stress 0.07428313 
    ## Run 243 stress 0.07970525 
    ## Run 244 stress 0.06942778 
    ## ... Procrustes: rmse 5.867361e-05  max resid 0.0001509305 
    ## ... Similar to previous best
    ## Run 245 stress 0.07428315 
    ## Run 246 stress 0.08340287 
    ## Run 247 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324993  max resid 0.03331902 
    ## Run 248 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324794  max resid 0.03331692 
    ## Run 249 stress 0.06978202 
    ## ... Procrustes: rmse 0.01315754  max resid 0.03312157 
    ## Run 250 stress 0.07428313 
    ## Run 251 stress 0.0834029 
    ## Run 252 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326195  max resid 0.03334488 
    ## Run 253 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325439  max resid 0.03332367 
    ## Run 254 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327151  max resid 0.03336204 
    ## Run 255 stress 0.06942776 
    ## ... Procrustes: rmse 3.596879e-06  max resid 9.059661e-06 
    ## ... Similar to previous best
    ## Run 256 stress 0.07970525 
    ## Run 257 stress 0.07970525 
    ## Run 258 stress 0.08340289 
    ## Run 259 stress 0.08340287 
    ## Run 260 stress 0.07428314 
    ## Run 261 stress 0.07428313 
    ## Run 262 stress 0.08340294 
    ## Run 263 stress 0.08340287 
    ## Run 264 stress 0.07428313 
    ## Run 265 stress 0.07428315 
    ## Run 266 stress 0.07970525 
    ## Run 267 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327431  max resid 0.033368 
    ## Run 268 stress 0.06942777 
    ## ... Procrustes: rmse 3.765217e-05  max resid 9.681495e-05 
    ## ... Similar to previous best
    ## Run 269 stress 0.07428313 
    ## Run 270 stress 0.06942777 
    ## ... Procrustes: rmse 3.550942e-05  max resid 9.427211e-05 
    ## ... Similar to previous best
    ## Run 271 stress 0.08340287 
    ## Run 272 stress 0.08340288 
    ## Run 273 stress 0.06942776 
    ## ... Procrustes: rmse 1.309592e-05  max resid 3.674549e-05 
    ## ... Similar to previous best
    ## Run 274 stress 0.07250812 
    ## Run 275 stress 0.06942776 
    ## ... Procrustes: rmse 1.022464e-05  max resid 2.36521e-05 
    ## ... Similar to previous best
    ## Run 276 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321843  max resid 0.03325293 
    ## Run 277 stress 0.08448454 
    ## Run 278 stress 0.06942778 
    ## ... Procrustes: rmse 5.892604e-05  max resid 0.0001517952 
    ## ... Similar to previous best
    ## Run 279 stress 0.07970525 
    ## Run 280 stress 0.07970526 
    ## Run 281 stress 0.07250814 
    ## Run 282 stress 0.07250813 
    ## Run 283 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323072  max resid 0.03327924 
    ## Run 284 stress 0.07970526 
    ## Run 285 stress 0.08448436 
    ## Run 286 stress 0.07428313 
    ## Run 287 stress 0.06942776 
    ## ... Procrustes: rmse 1.391118e-05  max resid 3.589061e-05 
    ## ... Similar to previous best
    ## Run 288 stress 0.07428314 
    ## Run 289 stress 0.07428319 
    ## Run 290 stress 0.07428314 
    ## Run 291 stress 0.07250812 
    ## Run 292 stress 0.06942776 
    ## ... Procrustes: rmse 3.374046e-05  max resid 8.727082e-05 
    ## ... Similar to previous best
    ## Run 293 stress 0.07970525 
    ## Run 294 stress 0.08340295 
    ## Run 295 stress 0.07428314 
    ## Run 296 stress 0.06942777 
    ## ... Procrustes: rmse 3.962836e-05  max resid 0.0001026724 
    ## ... Similar to previous best
    ## Run 297 stress 0.08340287 
    ## Run 298 stress 0.06942776 
    ## ... Procrustes: rmse 2.215109e-05  max resid 5.751679e-05 
    ## ... Similar to previous best
    ## Run 299 stress 0.07250812 
    ## Run 300 stress 0.08340289 
    ## Run 301 stress 0.07428319 
    ## Run 302 stress 0.08340296 
    ## Run 303 stress 0.07428317 
    ## Run 304 stress 0.07250813 
    ## Run 305 stress 0.07970525 
    ## Run 306 stress 0.06942776 
    ## ... Procrustes: rmse 4.385734e-06  max resid 1.128346e-05 
    ## ... Similar to previous best
    ## Run 307 stress 0.07250812 
    ## Run 308 stress 0.06978189 
    ## ... Procrustes: rmse 0.013224  max resid 0.03326448 
    ## Run 309 stress 0.07970525 
    ## Run 310 stress 0.07428312 
    ## Run 311 stress 0.07428315 
    ## Run 312 stress 0.08340287 
    ## Run 313 stress 0.07428316 
    ## Run 314 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326086  max resid 0.03334199 
    ## Run 315 stress 0.06978196 
    ## ... Procrustes: rmse 0.01327037  max resid 0.0333625 
    ## Run 316 stress 0.07844933 
    ## Run 317 stress 0.07428313 
    ## Run 318 stress 0.07428313 
    ## Run 319 stress 0.06942776 
    ## ... Procrustes: rmse 1.222754e-05  max resid 2.833861e-05 
    ## ... Similar to previous best
    ## Run 320 stress 0.08340287 
    ## Run 321 stress 0.06942776 
    ## ... Procrustes: rmse 1.374573e-05  max resid 3.586922e-05 
    ## ... Similar to previous best
    ## Run 322 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323828  max resid 0.03329675 
    ## Run 323 stress 0.07428316 
    ## Run 324 stress 0.07428317 
    ## Run 325 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327044  max resid 0.03336036 
    ## Run 326 stress 0.06978189 
    ## ... Procrustes: rmse 0.0132226  max resid 0.03326174 
    ## Run 327 stress 0.07428314 
    ## Run 328 stress 0.069782 
    ## ... Procrustes: rmse 0.01316301  max resid 0.03313182 
    ## Run 329 stress 0.07428313 
    ## Run 330 stress 0.07970525 
    ## Run 331 stress 0.07428315 
    ## Run 332 stress 0.08340288 
    ## Run 333 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326962  max resid 0.03336043 
    ## Run 334 stress 0.08448453 
    ## Run 335 stress 0.06942776 
    ## ... Procrustes: rmse 3.498231e-06  max resid 8.739199e-06 
    ## ... Similar to previous best
    ## Run 336 stress 0.07970525 
    ## Run 337 stress 0.07428314 
    ## Run 338 stress 0.08340287 
    ## Run 339 stress 0.069782 
    ## ... Procrustes: rmse 0.0132802  max resid 0.03338483 
    ## Run 340 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.822746e-06  max resid 6.789046e-06 
    ## ... Similar to previous best
    ## Run 341 stress 0.06942776 
    ## ... Procrustes: rmse 4.797833e-06  max resid 1.223265e-05 
    ## ... Similar to previous best
    ## Run 342 stress 0.08340293 
    ## Run 343 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326867  max resid 0.03335808 
    ## Run 344 stress 0.06942777 
    ## ... Procrustes: rmse 2.514124e-05  max resid 6.82778e-05 
    ## ... Similar to previous best
    ## Run 345 stress 0.07428315 
    ## Run 346 stress 0.07970526 
    ## Run 347 stress 0.07428313 
    ## Run 348 stress 0.07250813 
    ## Run 349 stress 0.07428313 
    ## Run 350 stress 0.07250812 
    ## Run 351 stress 0.07428315 
    ## Run 352 stress 0.07970525 
    ## Run 353 stress 0.06978196 
    ## ... Procrustes: rmse 0.01318569  max resid 0.0331845 
    ## Run 354 stress 0.08340286 
    ## Run 355 stress 0.069782 
    ## ... Procrustes: rmse 0.01328249  max resid 0.03338503 
    ## Run 356 stress 0.07250816 
    ## Run 357 stress 0.06942776 
    ## ... Procrustes: rmse 6.979662e-06  max resid 1.667426e-05 
    ## ... Similar to previous best
    ## Run 358 stress 0.06942777 
    ## ... Procrustes: rmse 1.594998e-05  max resid 4.559016e-05 
    ## ... Similar to previous best
    ## Run 359 stress 0.07250814 
    ## Run 360 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326159  max resid 0.03334128 
    ## Run 361 stress 0.08448456 
    ## Run 362 stress 0.08448451 
    ## Run 363 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325377  max resid 0.03332545 
    ## Run 364 stress 0.08340294 
    ## Run 365 stress 0.07428314 
    ## Run 366 stress 0.08340292 
    ## Run 367 stress 0.07970525 
    ## Run 368 stress 0.07970525 
    ## Run 369 stress 0.07970525 
    ## Run 370 stress 0.07970525 
    ## Run 371 stress 0.08340289 
    ## Run 372 stress 0.07428312 
    ## Run 373 stress 0.0834029 
    ## Run 374 stress 0.07250813 
    ## Run 375 stress 0.06942777 
    ## ... Procrustes: rmse 1.27941e-05  max resid 4.104888e-05 
    ## ... Similar to previous best
    ## Run 376 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322693  max resid 0.03326988 
    ## Run 377 stress 0.08340291 
    ## Run 378 stress 0.06978196 
    ## ... Procrustes: rmse 0.01317467  max resid 0.03315856 
    ## Run 379 stress 0.07428321 
    ## Run 380 stress 0.07428313 
    ## Run 381 stress 0.06978192 
    ## ... Procrustes: rmse 0.01322713  max resid 0.03326479 
    ## Run 382 stress 0.07970526 
    ## Run 383 stress 0.07250812 
    ## Run 384 stress 0.07428314 
    ## Run 385 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324315  max resid 0.03330429 
    ## Run 386 stress 0.06942776 
    ## ... Procrustes: rmse 1.088145e-05  max resid 3.200779e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.07970526 
    ## Run 388 stress 0.07428313 
    ## Run 389 stress 0.07970525 
    ## Run 390 stress 0.07428316 
    ## Run 391 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323469  max resid 0.03328352 
    ## Run 392 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327496  max resid 0.03336837 
    ## Run 393 stress 0.07970525 
    ## Run 394 stress 0.08340288 
    ## Run 395 stress 0.07970525 
    ## Run 396 stress 0.07428339 
    ## Run 397 stress 0.07428313 
    ## Run 398 stress 0.07428313 
    ## Run 399 stress 0.07428315 
    ## Run 400 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321746  max resid 0.03325078 
    ## Run 401 stress 0.07250812 
    ## Run 402 stress 0.08340291 
    ## Run 403 stress 0.07428314 
    ## Run 404 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321479  max resid 0.03324353 
    ## Run 405 stress 0.07250813 
    ## Run 406 stress 0.06978196 
    ## ... Procrustes: rmse 0.01324141  max resid 0.03329421 
    ## Run 407 stress 0.07428315 
    ## Run 408 stress 0.08340293 
    ## Run 409 stress 0.06942776 
    ## ... Procrustes: rmse 1.624667e-05  max resid 4.202314e-05 
    ## ... Similar to previous best
    ## Run 410 stress 0.08448451 
    ## Run 411 stress 0.08340287 
    ## Run 412 stress 0.08340291 
    ## Run 413 stress 0.08448458 
    ## Run 414 stress 0.08340294 
    ## Run 415 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323035  max resid 0.03327545 
    ## Run 416 stress 0.07428313 
    ## Run 417 stress 0.07970525 
    ## Run 418 stress 0.07428313 
    ## Run 419 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327169  max resid 0.03336499 
    ## Run 420 stress 0.06942776 
    ## ... Procrustes: rmse 1.875934e-05  max resid 4.856335e-05 
    ## ... Similar to previous best
    ## Run 421 stress 0.08340286 
    ## Run 422 stress 0.07250812 
    ## Run 423 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326333  max resid 0.03334317 
    ## Run 424 stress 0.08340287 
    ## Run 425 stress 0.06978192 
    ## ... Procrustes: rmse 0.0131938  max resid 0.0331993 
    ## Run 426 stress 0.08340288 
    ## Run 427 stress 0.07844951 
    ## Run 428 stress 0.08340287 
    ## Run 429 stress 0.08448442 
    ## Run 430 stress 0.07970526 
    ## Run 431 stress 0.08340287 
    ## Run 432 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322809  max resid 0.03327113 
    ## Run 433 stress 0.08340286 
    ## Run 434 stress 0.07970526 
    ## Run 435 stress 0.07428313 
    ## Run 436 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323405  max resid 0.03328789 
    ## Run 437 stress 0.07428318 
    ## Run 438 stress 0.2768842 
    ## Run 439 stress 0.06978197 
    ## ... Procrustes: rmse 0.01317087  max resid 0.03314767 
    ## Run 440 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324115  max resid 0.03329855 
    ## Run 441 stress 0.07250812 
    ## Run 442 stress 0.06942776 
    ## ... Procrustes: rmse 2.184312e-05  max resid 5.600371e-05 
    ## ... Similar to previous best
    ## Run 443 stress 0.08340288 
    ## Run 444 stress 0.07250815 
    ## Run 445 stress 0.08340294 
    ## Run 446 stress 0.07970527 
    ## Run 447 stress 0.07428313 
    ## Run 448 stress 0.06942776 
    ## ... Procrustes: rmse 7.833245e-06  max resid 1.868645e-05 
    ## ... Similar to previous best
    ## Run 449 stress 0.06942776 
    ## ... Procrustes: rmse 2.065562e-05  max resid 5.357121e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320408  max resid 0.03321972 
    ## Run 451 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324223  max resid 0.03330026 
    ## Run 452 stress 0.07428315 
    ## Run 453 stress 0.0834029 
    ## Run 454 stress 0.06942776 
    ## ... Procrustes: rmse 5.444957e-06  max resid 1.400084e-05 
    ## ... Similar to previous best
    ## Run 455 stress 0.08340296 
    ## Run 456 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323204  max resid 0.03328222 
    ## Run 457 stress 0.07428313 
    ## Run 458 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324061  max resid 0.03330038 
    ## Run 459 stress 0.07428312 
    ## Run 460 stress 0.0844844 
    ## Run 461 stress 0.07428312 
    ## Run 462 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324513  max resid 0.03330743 
    ## Run 463 stress 0.08340298 
    ## Run 464 stress 0.07970525 
    ## Run 465 stress 0.06942776 
    ## ... Procrustes: rmse 1.994037e-05  max resid 5.129603e-05 
    ## ... Similar to previous best
    ## Run 466 stress 0.07428313 
    ## Run 467 stress 0.07970526 
    ## Run 468 stress 0.07970525 
    ## Run 469 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325592  max resid 0.03333149 
    ## Run 470 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326648  max resid 0.03335453 
    ## Run 471 stress 0.069782 
    ## ... Procrustes: rmse 0.01316734  max resid 0.03314322 
    ## Run 472 stress 0.06942776 
    ## ... Procrustes: rmse 1.062081e-05  max resid 2.740935e-05 
    ## ... Similar to previous best
    ## Run 473 stress 0.07428315 
    ## Run 474 stress 0.06942776 
    ## ... Procrustes: rmse 1.833167e-05  max resid 4.716065e-05 
    ## ... Similar to previous best
    ## Run 475 stress 0.0844844 
    ## Run 476 stress 0.07970525 
    ## Run 477 stress 0.06942776 
    ## ... Procrustes: rmse 3.197518e-05  max resid 8.276376e-05 
    ## ... Similar to previous best
    ## Run 478 stress 0.08448446 
    ## Run 479 stress 0.07970526 
    ## Run 480 stress 0.07428315 
    ## Run 481 stress 0.06978195 
    ## ... Procrustes: rmse 0.01325835  max resid 0.03333216 
    ## Run 482 stress 0.07250813 
    ## Run 483 stress 0.07970525 
    ## Run 484 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326675  max resid 0.03335135 
    ## Run 485 stress 0.08448435 
    ## Run 486 stress 0.07970527 
    ## Run 487 stress 0.07428313 
    ## Run 488 stress 0.06942776 
    ## ... Procrustes: rmse 6.045086e-06  max resid 1.55693e-05 
    ## ... Similar to previous best
    ## Run 489 stress 0.0844845 
    ## Run 490 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.258139e-06  max resid 6.017379e-06 
    ## ... Similar to previous best
    ## Run 491 stress 0.0784492 
    ## Run 492 stress 0.07970526 
    ## Run 493 stress 0.08340287 
    ## Run 494 stress 0.08448435 
    ## Run 495 stress 0.08340292 
    ## Run 496 stress 0.08340288 
    ## Run 497 stress 0.06942778 
    ## ... Procrustes: rmse 1.664605e-05  max resid 5.374701e-05 
    ## ... Similar to previous best
    ## Run 498 stress 0.07250812 
    ## Run 499 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321394  max resid 0.03324481 
    ## Run 500 stress 0.07428313 
    ## *** Best solution repeated 2 times

``` r
# Mixed lakes
SD_beta_geo_M_NMDS <- metaMDS(SD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.001279543 
    ## Run 1 stress 0.00143331 
    ## ... Procrustes: rmse 0.002545359  max resid 0.003532005 
    ## ... Similar to previous best
    ## Run 2 stress 0.0004151886 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01869439  max resid 0.02584874 
    ## Run 3 stress 0.0004467355 
    ## ... Procrustes: rmse 0.02506942  max resid 0.05533164 
    ## Run 4 stress 0.1990774 
    ## Run 5 stress 9.017742e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.02460335  max resid 0.03371954 
    ## Run 6 stress 9.192602e-05 
    ## ... Procrustes: rmse 9.079513e-05  max resid 0.0002119056 
    ## ... Similar to previous best
    ## Run 7 stress 0.001334274 
    ## Run 8 stress 9.74241e-05 
    ## ... Procrustes: rmse 0.0006840046  max resid 0.0009773625 
    ## ... Similar to previous best
    ## Run 9 stress 0.001275227 
    ## Run 10 stress 9.991878e-05 
    ## ... Procrustes: rmse 0.007224472  max resid 0.009834956 
    ## Run 11 stress 0.0004716211 
    ## ... Procrustes: rmse 0.02621442  max resid 0.03599011 
    ## Run 12 stress 8.883685e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002367064  max resid 0.0004772825 
    ## ... Similar to previous best
    ## Run 13 stress 8.288232e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 7.107263e-05  max resid 0.0001393292 
    ## ... Similar to previous best
    ## Run 14 stress 0.0004757201 
    ## ... Procrustes: rmse 0.01597285  max resid 0.02182225 
    ## Run 15 stress 0.001196335 
    ## Run 16 stress 0.000464923 
    ## ... Procrustes: rmse 0.0157919  max resid 0.02157226 
    ## Run 17 stress 0.001241285 
    ## Run 18 stress 0.2842805 
    ## Run 19 stress 0.0006918359 
    ## Run 20 stress 0.0004735054 
    ## ... Procrustes: rmse 0.0159372  max resid 0.02177269 
    ## Run 21 stress 9.103486e-05 
    ## ... Procrustes: rmse 0.0001947639  max resid 0.0004128252 
    ## ... Similar to previous best
    ## Run 22 stress 9.494891e-05 
    ## ... Procrustes: rmse 0.0001429357  max resid 0.0002566049 
    ## ... Similar to previous best
    ## Run 23 stress 0.0005101747 
    ## ... Procrustes: rmse 0.01653606  max resid 0.02259891 
    ## Run 24 stress 0.0004809768 
    ## ... Procrustes: rmse 0.01604741  max resid 0.0219348 
    ## Run 25 stress 0.0001326132 
    ## ... Procrustes: rmse 0.008420785  max resid 0.01140037 
    ## Run 26 stress 0.001244038 
    ## Run 27 stress 0.00125852 
    ## Run 28 stress 0.2842598 
    ## Run 29 stress 7.427978e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0006991487  max resid 0.001015572 
    ## ... Similar to previous best
    ## Run 30 stress 9.748204e-05 
    ## ... Procrustes: rmse 0.0008031622  max resid 0.001571288 
    ## ... Similar to previous best
    ## Run 31 stress 0.2842805 
    ## Run 32 stress 0.1491957 
    ## Run 33 stress 0.1491957 
    ## Run 34 stress 0.001183699 
    ## Run 35 stress 9.754784e-05 
    ## ... Procrustes: rmse 0.0008074972  max resid 0.001152068 
    ## ... Similar to previous best
    ## Run 36 stress 0.00011999 
    ## ... Procrustes: rmse 0.007890218  max resid 0.01226729 
    ## Run 37 stress 0.0004582786 
    ## ... Procrustes: rmse 0.01553509  max resid 0.02285089 
    ## Run 38 stress 0.001059736 
    ## Run 39 stress 8.37101e-05 
    ## ... Procrustes: rmse 0.001945209  max resid 0.00263469 
    ## ... Similar to previous best
    ## Run 40 stress 0.001357978 
    ## Run 41 stress 0.1491957 
    ## Run 42 stress 0.0009132642 
    ## Run 43 stress 9.879091e-05 
    ## ... Procrustes: rmse 0.000648086  max resid 0.0009342831 
    ## ... Similar to previous best
    ## Run 44 stress 9.327121e-05 
    ## ... Procrustes: rmse 0.0007166804  max resid 0.001103979 
    ## ... Similar to previous best
    ## Run 45 stress 8.162552e-05 
    ## ... Procrustes: rmse 0.00075981  max resid 0.001106457 
    ## ... Similar to previous best
    ## Run 46 stress 0.0004336955 
    ## ... Procrustes: rmse 0.02439032  max resid 0.03362787 
    ## Run 47 stress 0.001224023 
    ## Run 48 stress 0.2696744 
    ## Run 49 stress 0.001340604 
    ## Run 50 stress 9.344033e-05 
    ## ... Procrustes: rmse 0.0007578422  max resid 0.001137718 
    ## ... Similar to previous best
    ## Run 51 stress 7.839292e-05 
    ## ... Procrustes: rmse 0.00087711  max resid 0.001929225 
    ## ... Similar to previous best
    ## Run 52 stress 0.001392915 
    ## Run 53 stress 0.001279995 
    ## Run 54 stress 9.257639e-05 
    ## ... Procrustes: rmse 0.0007329866  max resid 0.001171608 
    ## ... Similar to previous best
    ## Run 55 stress 0.1990774 
    ## Run 56 stress 9.032931e-05 
    ## ... Procrustes: rmse 0.002373701  max resid 0.004537678 
    ## ... Similar to previous best
    ## Run 57 stress 5.163138e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0007273262  max resid 0.001021869 
    ## ... Similar to previous best
    ## Run 58 stress 0.0004769795 
    ## ... Procrustes: rmse 0.01598632  max resid 0.02192063 
    ## Run 59 stress 0.001248091 
    ## Run 60 stress 0.001329074 
    ## Run 61 stress 0.0005061149 
    ## ... Procrustes: rmse 0.01650218  max resid 0.02263345 
    ## Run 62 stress 0.001357126 
    ## Run 63 stress 0.1491957 
    ## Run 64 stress 0.2520602 
    ## Run 65 stress 9.684231e-05 
    ## ... Procrustes: rmse 0.0002103861  max resid 0.0004082363 
    ## ... Similar to previous best
    ## Run 66 stress 0.1990774 
    ## Run 67 stress 0.001255163 
    ## Run 68 stress 0.1491957 
    ## Run 69 stress 9.112266e-05 
    ## ... Procrustes: rmse 0.0001014281  max resid 0.0001401837 
    ## ... Similar to previous best
    ## Run 70 stress 8.984321e-05 
    ## ... Procrustes: rmse 0.0001921728  max resid 0.0003489132 
    ## ... Similar to previous best
    ## Run 71 stress 0.2520602 
    ## Run 72 stress 9.967668e-05 
    ## ... Procrustes: rmse 0.0001088125  max resid 0.0001712048 
    ## ... Similar to previous best
    ## Run 73 stress 0.1990774 
    ## Run 74 stress 0.0001912901 
    ## ... Procrustes: rmse 0.0158737  max resid 0.02186311 
    ## Run 75 stress 9.999674e-05 
    ## ... Procrustes: rmse 0.0001791218  max resid 0.0003090205 
    ## ... Similar to previous best
    ## Run 76 stress 0.1491957 
    ## Run 77 stress 0.1491957 
    ## Run 78 stress 9.987657e-05 
    ## ... Procrustes: rmse 0.0001198786  max resid 0.0001862291 
    ## ... Similar to previous best
    ## Run 79 stress 9.282643e-05 
    ## ... Procrustes: rmse 8.911508e-05  max resid 0.0001264995 
    ## ... Similar to previous best
    ## Run 80 stress 0.2440272 
    ## Run 81 stress 0.001267273 
    ## Run 82 stress 0.001190309 
    ## Run 83 stress 8.249246e-05 
    ## ... Procrustes: rmse 8.431106e-05  max resid 0.0001326329 
    ## ... Similar to previous best
    ## Run 84 stress 0.0005814743 
    ## Run 85 stress 9.130431e-05 
    ## ... Procrustes: rmse 0.0001969094  max resid 0.0003433233 
    ## ... Similar to previous best
    ## Run 86 stress 0.0002424192 
    ## ... Procrustes: rmse 0.0187589  max resid 0.02585484 
    ## Run 87 stress 9.698771e-05 
    ## ... Procrustes: rmse 0.0001236549  max resid 0.0001660125 
    ## ... Similar to previous best
    ## Run 88 stress 9.779913e-05 
    ## ... Procrustes: rmse 8.648926e-05  max resid 0.0001307898 
    ## ... Similar to previous best
    ## Run 89 stress 0.001055295 
    ## Run 90 stress 9.317788e-05 
    ## ... Procrustes: rmse 0.0001974682  max resid 0.0003506246 
    ## ... Similar to previous best
    ## Run 91 stress 0.1990774 
    ## Run 92 stress 0.1491957 
    ## Run 93 stress 0.3083098 
    ## Run 94 stress 0.001206729 
    ## Run 95 stress 0.0005052209 
    ## ... Procrustes: rmse 0.01648262  max resid 0.02260767 
    ## Run 96 stress 0.0003339895 
    ## ... Procrustes: rmse 0.01340164  max resid 0.01835653 
    ## Run 97 stress 0.001371459 
    ## Run 98 stress 0.1776017 
    ## Run 99 stress 0.001284754 
    ## Run 100 stress 0.000576128 
    ## Run 101 stress 0.1776012 
    ## Run 102 stress 0.001315731 
    ## Run 103 stress 0.1776019 
    ## Run 104 stress 0.2848515 
    ## Run 105 stress 0.1990774 
    ## Run 106 stress 9.151121e-05 
    ## ... Procrustes: rmse 0.0001242491  max resid 0.0002312954 
    ## ... Similar to previous best
    ## Run 107 stress 9.807255e-05 
    ## ... Procrustes: rmse 0.000211945  max resid 0.0003707386 
    ## ... Similar to previous best
    ## Run 108 stress 9.926716e-05 
    ## ... Procrustes: rmse 0.0005655656  max resid 0.0008150626 
    ## ... Similar to previous best
    ## Run 109 stress 0.0004718299 
    ## ... Procrustes: rmse 0.01592078  max resid 0.02182973 
    ## Run 110 stress 0.001243029 
    ## Run 111 stress 0.2578642 
    ## Run 112 stress 0.2795549 
    ## Run 113 stress 8.849682e-05 
    ## ... Procrustes: rmse 0.0001172438  max resid 0.0002144889 
    ## ... Similar to previous best
    ## Run 114 stress 0.001110317 
    ## Run 115 stress 0.001453578 
    ## Run 116 stress 0.2528294 
    ## Run 117 stress 0.001203199 
    ## Run 118 stress 9.07987e-05 
    ## ... Procrustes: rmse 9.614649e-05  max resid 0.0001488274 
    ## ... Similar to previous best
    ## Run 119 stress 9.679915e-05 
    ## ... Procrustes: rmse 0.0001735267  max resid 0.0003136558 
    ## ... Similar to previous best
    ## Run 120 stress 9.29638e-05 
    ## ... Procrustes: rmse 0.0001070378  max resid 0.0001680941 
    ## ... Similar to previous best
    ## Run 121 stress 0.0004425972 
    ## ... Procrustes: rmse 0.01543316  max resid 0.02115686 
    ## Run 122 stress 9.418308e-05 
    ## ... Procrustes: rmse 0.0001879162  max resid 0.0003232363 
    ## ... Similar to previous best
    ## Run 123 stress 0.1491957 
    ## Run 124 stress 0.1491957 
    ## Run 125 stress 0.0004539769 
    ## ... Procrustes: rmse 0.01561026  max resid 0.02139821 
    ## Run 126 stress 0.001272327 
    ## Run 127 stress 0.3083098 
    ## Run 128 stress 9.451584e-05 
    ## ... Procrustes: rmse 0.0001088598  max resid 0.0001712019 
    ## ... Similar to previous best
    ## Run 129 stress 9.981519e-05 
    ## ... Procrustes: rmse 0.00016704  max resid 0.0003043985 
    ## ... Similar to previous best
    ## Run 130 stress 0.000173458 
    ## ... Procrustes: rmse 0.009657885  max resid 0.01318775 
    ## Run 131 stress 0.001117001 
    ## Run 132 stress 0.2520602 
    ## Run 133 stress 7.841734e-05 
    ## ... Procrustes: rmse 0.001407501  max resid 0.001860696 
    ## ... Similar to previous best
    ## Run 134 stress 0.0004202588 
    ## ... Procrustes: rmse 0.0150209  max resid 0.02058806 
    ## Run 135 stress 9.951158e-05 
    ## ... Procrustes: rmse 0.0001286894  max resid 0.0001735131 
    ## ... Similar to previous best
    ## Run 136 stress 8.457592e-05 
    ## ... Procrustes: rmse 0.0002967479  max resid 0.000501023 
    ## ... Similar to previous best
    ## Run 137 stress 0.1491957 
    ## Run 138 stress 0.001419836 
    ## Run 139 stress 9.838033e-05 
    ## ... Procrustes: rmse 0.0004156539  max resid 0.0006122663 
    ## ... Similar to previous best
    ## Run 140 stress 0.00132958 
    ## Run 141 stress 0.001247743 
    ## Run 142 stress 9.453657e-05 
    ## ... Procrustes: rmse 0.0001738855  max resid 0.0003004344 
    ## ... Similar to previous best
    ## Run 143 stress 9.172597e-05 
    ## ... Procrustes: rmse 0.0001509554  max resid 0.0003081693 
    ## ... Similar to previous best
    ## Run 144 stress 0.0004948679 
    ## ... Procrustes: rmse 0.01631738  max resid 0.02237724 
    ## Run 145 stress 0.001277623 
    ## Run 146 stress 0.001379681 
    ## Run 147 stress 0.001454221 
    ## Run 148 stress 0.3083099 
    ## Run 149 stress 9.892812e-05 
    ## ... Procrustes: rmse 0.0001163777  max resid 0.0001494353 
    ## ... Similar to previous best
    ## Run 150 stress 9.199075e-05 
    ## ... Procrustes: rmse 0.0001986061  max resid 0.0003496705 
    ## ... Similar to previous best
    ## Run 151 stress 0.000381296 
    ## ... Procrustes: rmse 0.01432412  max resid 0.01962676 
    ## Run 152 stress 0.0003646566 
    ## ... Procrustes: rmse 0.01400751  max resid 0.01918989 
    ## Run 153 stress 9.699217e-05 
    ## ... Procrustes: rmse 0.004329409  max resid 0.005914444 
    ## ... Similar to previous best
    ## Run 154 stress 0.001290109 
    ## Run 155 stress 0.2528294 
    ## Run 156 stress 0.0004070333 
    ## ... Procrustes: rmse 0.01454421  max resid 0.01993129 
    ## Run 157 stress 0.001021978 
    ## Run 158 stress 0.1990774 
    ## Run 159 stress 0.0008805882 
    ## Run 160 stress 0.1491957 
    ## Run 161 stress 0.0001027462 
    ## ... Procrustes: rmse 0.00743356  max resid 0.0101183 
    ## Run 162 stress 0.001372964 
    ## Run 163 stress 7.166126e-05 
    ## ... Procrustes: rmse 0.0001374359  max resid 0.0002608351 
    ## ... Similar to previous best
    ## Run 164 stress 8.786894e-05 
    ## ... Procrustes: rmse 0.0001398023  max resid 0.000286013 
    ## ... Similar to previous best
    ## Run 165 stress 0.1491957 
    ## Run 166 stress 8.059975e-05 
    ## ... Procrustes: rmse 0.0001442582  max resid 0.0002209068 
    ## ... Similar to previous best
    ## Run 167 stress 0.001276808 
    ## Run 168 stress 9.32939e-05 
    ## ... Procrustes: rmse 0.0001383338  max resid 0.0002686482 
    ## ... Similar to previous best
    ## Run 169 stress 8.530542e-05 
    ## ... Procrustes: rmse 0.0001569469  max resid 0.0003156753 
    ## ... Similar to previous best
    ## Run 170 stress 0.000446545 
    ## ... Procrustes: rmse 0.01547839  max resid 0.02122017 
    ## Run 171 stress 9.962544e-05 
    ## ... Procrustes: rmse 0.0001312126  max resid 0.0001745649 
    ## ... Similar to previous best
    ## Run 172 stress 0.1990774 
    ## Run 173 stress 0.001434121 
    ## Run 174 stress 9.854081e-05 
    ## ... Procrustes: rmse 0.0001133227  max resid 0.0001635652 
    ## ... Similar to previous best
    ## Run 175 stress 9.019459e-05 
    ## ... Procrustes: rmse 0.0002015125  max resid 0.0003500347 
    ## ... Similar to previous best
    ## Run 176 stress 0.0004690447 
    ## ... Procrustes: rmse 0.01588622  max resid 0.0217821 
    ## Run 177 stress 0.00140339 
    ## Run 178 stress 0.001101004 
    ## Run 179 stress 9.767514e-05 
    ## ... Procrustes: rmse 0.0001479944  max resid 0.0002886174 
    ## ... Similar to previous best
    ## Run 180 stress 0.0005905388 
    ## Run 181 stress 9.475238e-05 
    ## ... Procrustes: rmse 0.0001156326  max resid 0.0001844591 
    ## ... Similar to previous best
    ## Run 182 stress 0.1491957 
    ## Run 183 stress 7.123072e-05 
    ## ... Procrustes: rmse 0.001697231  max resid 0.002325826 
    ## ... Similar to previous best
    ## Run 184 stress 0.001441923 
    ## Run 185 stress 9.084272e-05 
    ## ... Procrustes: rmse 9.966424e-05  max resid 0.0001670361 
    ## ... Similar to previous best
    ## Run 186 stress 9.334119e-05 
    ## ... Procrustes: rmse 0.0001723933  max resid 0.0003402055 
    ## ... Similar to previous best
    ## Run 187 stress 9.211693e-05 
    ## ... Procrustes: rmse 0.0001011984  max resid 0.0001591369 
    ## ... Similar to previous best
    ## Run 188 stress 9.109952e-05 
    ## ... Procrustes: rmse 0.0001261156  max resid 0.0002135764 
    ## ... Similar to previous best
    ## Run 189 stress 7.600325e-05 
    ## ... Procrustes: rmse 6.901196e-05  max resid 0.000109171 
    ## ... Similar to previous best
    ## Run 190 stress 0.00134914 
    ## Run 191 stress 0.001332345 
    ## Run 192 stress 0.00134671 
    ## Run 193 stress 8.410415e-05 
    ## ... Procrustes: rmse 0.0001558004  max resid 0.0002924836 
    ## ... Similar to previous best
    ## Run 194 stress 8.992868e-05 
    ## ... Procrustes: rmse 0.000125679  max resid 0.0002179096 
    ## ... Similar to previous best
    ## Run 195 stress 0.0002156587 
    ## ... Procrustes: rmse 0.01768773  max resid 0.02437367 
    ## Run 196 stress 0.0004595252 
    ## ... Procrustes: rmse 0.01572367  max resid 0.02155927 
    ## Run 197 stress 7.401646e-05 
    ## ... Procrustes: rmse 0.0001568031  max resid 0.0002808045 
    ## ... Similar to previous best
    ## Run 198 stress 6.154078e-05 
    ## ... Procrustes: rmse 0.002359437  max resid 0.003262137 
    ## ... Similar to previous best
    ## Run 199 stress 0.001338118 
    ## Run 200 stress 0.001483278 
    ## Run 201 stress 0.001431781 
    ## Run 202 stress 9.540211e-05 
    ## ... Procrustes: rmse 0.0001098922  max resid 0.0001554603 
    ## ... Similar to previous best
    ## Run 203 stress 9.454958e-05 
    ## ... Procrustes: rmse 0.000132318  max resid 0.0002023251 
    ## ... Similar to previous best
    ## Run 204 stress 0.001129277 
    ## Run 205 stress 0.177602 
    ## Run 206 stress 9.865746e-05 
    ## ... Procrustes: rmse 0.0001776067  max resid 0.0002625443 
    ## ... Similar to previous best
    ## Run 207 stress 0.001191765 
    ## Run 208 stress 0.001193217 
    ## Run 209 stress 0.2440272 
    ## Run 210 stress 7.279002e-05 
    ## ... Procrustes: rmse 6.169568e-05  max resid 9.690555e-05 
    ## ... Similar to previous best
    ## Run 211 stress 0.1491957 
    ## Run 212 stress 9.270226e-05 
    ## ... Procrustes: rmse 9.741673e-05  max resid 0.0001484307 
    ## ... Similar to previous best
    ## Run 213 stress 9.224009e-05 
    ## ... Procrustes: rmse 0.000199313  max resid 0.0003396638 
    ## ... Similar to previous best
    ## Run 214 stress 0.0004547118 
    ## ... Procrustes: rmse 0.01564195  max resid 0.02144484 
    ## Run 215 stress 0.0004417644 
    ## ... Procrustes: rmse 0.02534556  max resid 0.03497096 
    ## Run 216 stress 9.2701e-05 
    ## ... Procrustes: rmse 0.0001664701  max resid 0.000329982 
    ## ... Similar to previous best
    ## Run 217 stress 9.821318e-05 
    ## ... Procrustes: rmse 0.0001182007  max resid 0.0001835052 
    ## ... Similar to previous best
    ## Run 218 stress 9.955088e-05 
    ## ... Procrustes: rmse 0.0001042118  max resid 0.0001640854 
    ## ... Similar to previous best
    ## Run 219 stress 5.041912e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000189604  max resid 0.0002962005 
    ## ... Similar to previous best
    ## Run 220 stress 0.2797319 
    ## Run 221 stress 9.615865e-05 
    ## ... Procrustes: rmse 0.0002298991  max resid 0.0004727348 
    ## ... Similar to previous best
    ## Run 222 stress 7.605326e-05 
    ## ... Procrustes: rmse 0.0002414776  max resid 0.0004549068 
    ## ... Similar to previous best
    ## Run 223 stress 0.3083098 
    ## Run 224 stress 0.001203268 
    ## Run 225 stress 0.00137825 
    ## Run 226 stress 0.001256125 
    ## Run 227 stress 9.016557e-05 
    ## ... Procrustes: rmse 0.0002170577  max resid 0.0004927731 
    ## ... Similar to previous best
    ## Run 228 stress 0.1491957 
    ## Run 229 stress 0.001276227 
    ## Run 230 stress 0.1776017 
    ## Run 231 stress 0.0012856 
    ## Run 232 stress 9.417755e-05 
    ## ... Procrustes: rmse 0.0002513605  max resid 0.0004949309 
    ## ... Similar to previous best
    ## Run 233 stress 9.381913e-05 
    ## ... Procrustes: rmse 0.0001496492  max resid 0.0003029233 
    ## ... Similar to previous best
    ## Run 234 stress 0.001236779 
    ## Run 235 stress 0.0004630759 
    ## ... Procrustes: rmse 0.01574724  max resid 0.02219075 
    ## Run 236 stress 0.1491957 
    ## Run 237 stress 8.818446e-05 
    ## ... Procrustes: rmse 0.0002870724  max resid 0.0004851135 
    ## ... Similar to previous best
    ## Run 238 stress 9.803264e-05 
    ## ... Procrustes: rmse 0.0002334657  max resid 0.0004791361 
    ## ... Similar to previous best
    ## Run 239 stress 0.0005048962 
    ## ... Procrustes: rmse 0.01644609  max resid 0.02315502 
    ## Run 240 stress 0.001240294 
    ## Run 241 stress 8.255019e-05 
    ## ... Procrustes: rmse 0.0002538671  max resid 0.0004325999 
    ## ... Similar to previous best
    ## Run 242 stress 0.001309481 
    ## Run 243 stress 0.0005380835 
    ## ... Procrustes: rmse 0.02777289  max resid 0.03833816 
    ## Run 244 stress 0.001319467 
    ## Run 245 stress 9.232138e-05 
    ## ... Procrustes: rmse 0.0002252335  max resid 0.0004912498 
    ## ... Similar to previous best
    ## Run 246 stress 0.00133373 
    ## Run 247 stress 0.00131892 
    ## Run 248 stress 0.001158544 
    ## Run 249 stress 0.1491957 
    ## Run 250 stress 9.270255e-05 
    ## ... Procrustes: rmse 0.0002397758  max resid 0.0004334478 
    ## ... Similar to previous best
    ## Run 251 stress 9.622957e-05 
    ## ... Procrustes: rmse 0.0002161946  max resid 0.0004798153 
    ## ... Similar to previous best
    ## Run 252 stress 0.00132306 
    ## Run 253 stress 0.001345746 
    ## Run 254 stress 0.001436594 
    ## Run 255 stress 0.001137153 
    ## Run 256 stress 0.0005070488 
    ## ... Procrustes: rmse 0.01641356  max resid 0.02311085 
    ## Run 257 stress 9.672795e-05 
    ## ... Procrustes: rmse 0.0002799821  max resid 0.0005372981 
    ## ... Similar to previous best
    ## Run 258 stress 0.0004912414 
    ## ... Procrustes: rmse 0.01622165  max resid 0.02284602 
    ## Run 259 stress 0.001343233 
    ## Run 260 stress 9.556704e-05 
    ## ... Procrustes: rmse 0.0001900163  max resid 0.0003878827 
    ## ... Similar to previous best
    ## Run 261 stress 0.001260526 
    ## Run 262 stress 0.0004511515 
    ## ... Procrustes: rmse 0.01554285  max resid 0.02191417 
    ## Run 263 stress 0.00134308 
    ## Run 264 stress 9.144873e-05 
    ## ... Procrustes: rmse 0.0001951195  max resid 0.0003418384 
    ## ... Similar to previous best
    ## Run 265 stress 9.998318e-05 
    ## ... Procrustes: rmse 0.000340926  max resid 0.0005817063 
    ## ... Similar to previous best
    ## Run 266 stress 0.1491957 
    ## Run 267 stress 0.001399319 
    ## Run 268 stress 0.001354975 
    ## Run 269 stress 0.0013086 
    ## Run 270 stress 0.000479151 
    ## ... Procrustes: rmse 0.02621671  max resid 0.03618106 
    ## Run 271 stress 0.1491957 
    ## Run 272 stress 0.177602 
    ## Run 273 stress 0.2570116 
    ## Run 274 stress 9.305465e-05 
    ## ... Procrustes: rmse 0.0001907428  max resid 0.0003567631 
    ## ... Similar to previous best
    ## Run 275 stress 0.0004599027 
    ## ... Procrustes: rmse 0.02568027  max resid 0.03543782 
    ## Run 276 stress 9.882665e-05 
    ## ... Procrustes: rmse 0.0001502369  max resid 0.00030584 
    ## ... Similar to previous best
    ## Run 277 stress 8.674921e-05 
    ## ... Procrustes: rmse 0.0002209214  max resid 0.0004370076 
    ## ... Similar to previous best
    ## Run 278 stress 0.0002480416 
    ## ... Procrustes: rmse 0.0187804  max resid 0.02588804 
    ## Run 279 stress 0.0004821446 
    ## ... Procrustes: rmse 0.01598275  max resid 0.02251527 
    ## Run 280 stress 0.2494706 
    ## Run 281 stress 9.718291e-05 
    ## ... Procrustes: rmse 0.0002304908  max resid 0.0004762433 
    ## ... Similar to previous best
    ## Run 282 stress 0.1776011 
    ## Run 283 stress 0.1491957 
    ## Run 284 stress 0.0004590575 
    ## ... Procrustes: rmse 0.02565444  max resid 0.03540213 
    ## Run 285 stress 0.1491957 
    ## Run 286 stress 9.545668e-05 
    ## ... Procrustes: rmse 0.0001843998  max resid 0.000331394 
    ## ... Similar to previous best
    ## Run 287 stress 0.001261558 
    ## Run 288 stress 9.230862e-05 
    ## ... Procrustes: rmse 0.0002298771  max resid 0.0004661734 
    ## ... Similar to previous best
    ## Run 289 stress 9.542827e-05 
    ## ... Procrustes: rmse 0.0002588467  max resid 0.0004939621 
    ## ... Similar to previous best
    ## Run 290 stress 0.2528294 
    ## Run 291 stress 0.2520602 
    ## Run 292 stress 9.366583e-05 
    ## ... Procrustes: rmse 0.007438133  max resid 0.01025125 
    ## Run 293 stress 0.177602 
    ## Run 294 stress 0.1491957 
    ## Run 295 stress 0.1491957 
    ## Run 296 stress 9.696143e-05 
    ## ... Procrustes: rmse 0.0002315146  max resid 0.0004992836 
    ## ... Similar to previous best
    ## Run 297 stress 9.192674e-05 
    ## ... Procrustes: rmse 0.0002203761  max resid 0.0004854497 
    ## ... Similar to previous best
    ## Run 298 stress 0.0001127045 
    ## ... Procrustes: rmse 0.007741839  max resid 0.0111476 
    ## Run 299 stress 7.433508e-05 
    ## ... Procrustes: rmse 0.0002215332  max resid 0.0004699447 
    ## ... Similar to previous best
    ## Run 300 stress 0.1491957 
    ## Run 301 stress 0.0004859408 
    ## ... Procrustes: rmse 0.01613237  max resid 0.02272207 
    ## Run 302 stress 0.2933522 
    ## Run 303 stress 8.309965e-05 
    ## ... Procrustes: rmse 0.0002172149  max resid 0.0004799517 
    ## ... Similar to previous best
    ## Run 304 stress 0.001161856 
    ## Run 305 stress 0.0004848169 
    ## ... Procrustes: rmse 0.01607162  max resid 0.02263887 
    ## Run 306 stress 9.45593e-05 
    ## ... Procrustes: rmse 0.000220013  max resid 0.0005057304 
    ## ... Similar to previous best
    ## Run 307 stress 0.001304193 
    ## Run 308 stress 6.897358e-05 
    ## ... Procrustes: rmse 0.0004057039  max resid 0.0005754897 
    ## ... Similar to previous best
    ## Run 309 stress 0.2842805 
    ## Run 310 stress 0.2795549 
    ## Run 311 stress 0.1491957 
    ## Run 312 stress 9.26429e-05 
    ## ... Procrustes: rmse 0.0002203249  max resid 0.0004886729 
    ## ... Similar to previous best
    ## Run 313 stress 0.0002382558 
    ## ... Procrustes: rmse 0.01841378  max resid 0.0253809 
    ## Run 314 stress 0.001199716 
    ## Run 315 stress 0.1491957 
    ## Run 316 stress 0.001275978 
    ## Run 317 stress 0.1491957 
    ## Run 318 stress 9.535837e-05 
    ## ... Procrustes: rmse 0.0002718349  max resid 0.0004388457 
    ## ... Similar to previous best
    ## Run 319 stress 8.810282e-05 
    ## ... Procrustes: rmse 0.0001507869  max resid 0.0003030898 
    ## ... Similar to previous best
    ## Run 320 stress 0.001306919 
    ## Run 321 stress 0.2797319 
    ## Run 322 stress 8.920959e-05 
    ## ... Procrustes: rmse 0.0002919978  max resid 0.0005180909 
    ## ... Similar to previous best
    ## Run 323 stress 9.140985e-05 
    ## ... Procrustes: rmse 0.0002478912  max resid 0.0004886732 
    ## ... Similar to previous best
    ## Run 324 stress 9.488497e-05 
    ## ... Procrustes: rmse 0.0001573882  max resid 0.0002981452 
    ## ... Similar to previous best
    ## Run 325 stress 9.904017e-05 
    ## ... Procrustes: rmse 0.0001972704  max resid 0.0003601112 
    ## ... Similar to previous best
    ## Run 326 stress 0.1491957 
    ## Run 327 stress 0.2520602 
    ## Run 328 stress 0.2842805 
    ## Run 329 stress 0.001381842 
    ## Run 330 stress 6.107709e-05 
    ## ... Procrustes: rmse 0.0002178514  max resid 0.0004009549 
    ## ... Similar to previous best
    ## Run 331 stress 0.0004627572 
    ## ... Procrustes: rmse 0.01573214  max resid 0.02217847 
    ## Run 332 stress 9.045028e-05 
    ## ... Procrustes: rmse 0.0001535203  max resid 0.0003005938 
    ## ... Similar to previous best
    ## Run 333 stress 9.024668e-05 
    ## ... Procrustes: rmse 0.0002475772  max resid 0.0004837087 
    ## ... Similar to previous best
    ## Run 334 stress 9.374091e-05 
    ## ... Procrustes: rmse 0.0002642793  max resid 0.0004534515 
    ## ... Similar to previous best
    ## Run 335 stress 0.0004979684 
    ## ... Procrustes: rmse 0.01632182  max resid 0.02299048 
    ## Run 336 stress 0.001442298 
    ## Run 337 stress 0.2842805 
    ## Run 338 stress 0.0008663837 
    ## Run 339 stress 0.2646972 
    ## Run 340 stress 0.0004840951 
    ## ... Procrustes: rmse 0.016099  max resid 0.0226764 
    ## Run 341 stress 0.001325625 
    ## Run 342 stress 0.001310695 
    ## Run 343 stress 0.0001184157 
    ## ... Procrustes: rmse 0.007944327  max resid 0.01142113 
    ## Run 344 stress 0.001238972 
    ## Run 345 stress 9.99912e-05 
    ## ... Procrustes: rmse 0.0002701247  max resid 0.0004365215 
    ## ... Similar to previous best
    ## Run 346 stress 9.320211e-05 
    ## ... Procrustes: rmse 0.0001496495  max resid 0.0003030354 
    ## ... Similar to previous best
    ## Run 347 stress 9.342925e-05 
    ## ... Procrustes: rmse 0.0001839664  max resid 0.0003750565 
    ## ... Similar to previous best
    ## Run 348 stress 0.001351141 
    ## Run 349 stress 5.362599e-05 
    ## ... Procrustes: rmse 0.0001073664  max resid 0.0001624957 
    ## ... Similar to previous best
    ## Run 350 stress 0.3083098 
    ## Run 351 stress 0.1491957 
    ## Run 352 stress 0.1776012 
    ## Run 353 stress 8.044276e-05 
    ## ... Procrustes: rmse 0.0003450025  max resid 0.0005185404 
    ## ... Similar to previous best
    ## Run 354 stress 0.001210642 
    ## Run 355 stress 0.001071592 
    ## Run 356 stress 0.0004266654 
    ## ... Procrustes: rmse 0.01511428  max resid 0.02131728 
    ## Run 357 stress 0.0001118226 
    ## ... Procrustes: rmse 0.007719146  max resid 0.01111012 
    ## Run 358 stress 9.17694e-05 
    ## ... Procrustes: rmse 0.0002465515  max resid 0.00047859 
    ## ... Similar to previous best
    ## Run 359 stress 9.714131e-05 
    ## ... Procrustes: rmse 0.000194988  max resid 0.0003986982 
    ## ... Similar to previous best
    ## Run 360 stress 9.794167e-05 
    ## ... Procrustes: rmse 0.007114841  max resid 0.01029674 
    ## Run 361 stress 0.001197389 
    ## Run 362 stress 0.001270963 
    ## Run 363 stress 9.370408e-05 
    ## ... Procrustes: rmse 0.0002182089  max resid 0.0004748591 
    ## ... Similar to previous best
    ## Run 364 stress 0.00048744 
    ## ... Procrustes: rmse 0.01615793  max resid 0.02275815 
    ## Run 365 stress 0.0006278368 
    ## Run 366 stress 8.4702e-05 
    ## ... Procrustes: rmse 0.001692197  max resid 0.002309748 
    ## ... Similar to previous best
    ## Run 367 stress 0.0004573328 
    ## ... Procrustes: rmse 0.01547672  max resid 0.02181867 
    ## Run 368 stress 8.686405e-05 
    ## ... Procrustes: rmse 0.0008406287  max resid 0.001573548 
    ## ... Similar to previous best
    ## Run 369 stress 0.0004398438 
    ## ... Procrustes: rmse 0.0153475  max resid 0.02163943 
    ## Run 370 stress 0.001268195 
    ## Run 371 stress 0.0003379358 
    ## ... Procrustes: rmse 0.01344768  max resid 0.01901781 
    ## Run 372 stress 0.001235286 
    ## Run 373 stress 0.3083098 
    ## Run 374 stress 0.1491957 
    ## Run 375 stress 9.625649e-05 
    ## ... Procrustes: rmse 0.0002183126  max resid 0.0004853704 
    ## ... Similar to previous best
    ## Run 376 stress 8.942137e-05 
    ## ... Procrustes: rmse 0.0008469309  max resid 0.001634664 
    ## ... Similar to previous best
    ## Run 377 stress 0.001498572 
    ## Run 378 stress 0.001296121 
    ## Run 379 stress 9.741777e-05 
    ## ... Procrustes: rmse 0.0001926098  max resid 0.0002453076 
    ## ... Similar to previous best
    ## Run 380 stress 0.001279447 
    ## Run 381 stress 0.001089638 
    ## Run 382 stress 9.557095e-05 
    ## ... Procrustes: rmse 0.000252971  max resid 0.0004921381 
    ## ... Similar to previous best
    ## Run 383 stress 9.56586e-05 
    ## ... Procrustes: rmse 0.0002121319  max resid 0.0004635914 
    ## ... Similar to previous best
    ## Run 384 stress 0.001438408 
    ## Run 385 stress 0.00113315 
    ## Run 386 stress 0.0007982967 
    ## Run 387 stress 8.58359e-05 
    ## ... Procrustes: rmse 0.0001579339  max resid 0.0003189727 
    ## ... Similar to previous best
    ## Run 388 stress 9.144537e-05 
    ## ... Procrustes: rmse 0.0002224816  max resid 0.0004868654 
    ## ... Similar to previous best
    ## Run 389 stress 9.377827e-05 
    ## ... Procrustes: rmse 0.0002299855  max resid 0.0005042912 
    ## ... Similar to previous best
    ## Run 390 stress 8.328682e-05 
    ## ... Procrustes: rmse 0.0002596347  max resid 0.0004373191 
    ## ... Similar to previous best
    ## Run 391 stress 0.1491957 
    ## Run 392 stress 0.0002330567 
    ## ... Procrustes: rmse 0.01116145  max resid 0.01586239 
    ## Run 393 stress 9.087972e-05 
    ## ... Procrustes: rmse 0.0002276813  max resid 0.0004935284 
    ## ... Similar to previous best
    ## Run 394 stress 0.0003382173 
    ## ... Procrustes: rmse 0.01334638  max resid 0.01887649 
    ## Run 395 stress 9.738206e-05 
    ## ... Procrustes: rmse 0.000235674  max resid 0.000511834 
    ## ... Similar to previous best
    ## Run 396 stress 0.001273042 
    ## Run 397 stress 0.0004788614 
    ## ... Procrustes: rmse 0.01601387  max resid 0.02255875 
    ## Run 398 stress 0.1491957 
    ## Run 399 stress 0.1491957 
    ## Run 400 stress 0.3083098 
    ## Run 401 stress 0.1990774 
    ## Run 402 stress 9.954942e-05 
    ## ... Procrustes: rmse 0.0001970266  max resid 0.0003696588 
    ## ... Similar to previous best
    ## Run 403 stress 9.885852e-05 
    ## ... Procrustes: rmse 0.0002863285  max resid 0.0004401249 
    ## ... Similar to previous best
    ## Run 404 stress 7.729512e-05 
    ## ... Procrustes: rmse 0.002294502  max resid 0.003667913 
    ## ... Similar to previous best
    ## Run 405 stress 8.622008e-05 
    ## ... Procrustes: rmse 0.0001522284  max resid 0.000300912 
    ## ... Similar to previous best
    ## Run 406 stress 9.399551e-05 
    ## ... Procrustes: rmse 0.0002682297  max resid 0.0004769708 
    ## ... Similar to previous best
    ## Run 407 stress 0.1491957 
    ## Run 408 stress 9.682378e-05 
    ## ... Procrustes: rmse 0.0001493008  max resid 0.0003031542 
    ## ... Similar to previous best
    ## Run 409 stress 9.56795e-05 
    ## ... Procrustes: rmse 0.0001492293  max resid 0.0003036849 
    ## ... Similar to previous best
    ## Run 410 stress 0.001430273 
    ## Run 411 stress 0.1491957 
    ## Run 412 stress 0.2852153 
    ## Run 413 stress 0.0004598444 
    ## ... Procrustes: rmse 0.01569348  max resid 0.02211659 
    ## Run 414 stress 9.452122e-05 
    ## ... Procrustes: rmse 0.0002760943  max resid 0.0005224456 
    ## ... Similar to previous best
    ## Run 415 stress 9.095921e-05 
    ## ... Procrustes: rmse 0.0002301917  max resid 0.0004954447 
    ## ... Similar to previous best
    ## Run 416 stress 0.0005947222 
    ## Run 417 stress 9.38987e-05 
    ## ... Procrustes: rmse 0.000192836  max resid 0.0003288141 
    ## ... Similar to previous best
    ## Run 418 stress 0.001201478 
    ## Run 419 stress 0.001342371 
    ## Run 420 stress 0.0004784044 
    ## ... Procrustes: rmse 0.01600151  max resid 0.0225426 
    ## Run 421 stress 8.414929e-05 
    ## ... Procrustes: rmse 0.001344859  max resid 0.002297191 
    ## ... Similar to previous best
    ## Run 422 stress 0.001284104 
    ## Run 423 stress 0.2854395 
    ## Run 424 stress 9.69938e-05 
    ## ... Procrustes: rmse 0.0002169285  max resid 0.0004819245 
    ## ... Similar to previous best
    ## Run 425 stress 8.620635e-05 
    ## ... Procrustes: rmse 0.0002904614  max resid 0.0004303025 
    ## ... Similar to previous best
    ## Run 426 stress 0.001321007 
    ## Run 427 stress 0.1491957 
    ## Run 428 stress 0.0004003309 
    ## ... Procrustes: rmse 0.01462097  max resid 0.02064431 
    ## Run 429 stress 0.1491957 
    ## Run 430 stress 0.1491957 
    ## Run 431 stress 9.422691e-05 
    ## ... Procrustes: rmse 0.0001919397  max resid 0.0003220852 
    ## ... Similar to previous best
    ## Run 432 stress 9.82596e-05 
    ## ... Procrustes: rmse 0.0001500118  max resid 0.0003002558 
    ## ... Similar to previous best
    ## Run 433 stress 0.0004854957 
    ## ... Procrustes: rmse 0.01612581  max resid 0.02271349 
    ## Run 434 stress 0.0004595965 
    ## ... Procrustes: rmse 0.01568762  max resid 0.0221088 
    ## Run 435 stress 0.177601 
    ## Run 436 stress 9.188043e-05 
    ## ... Procrustes: rmse 0.0002227936  max resid 0.0003428494 
    ## ... Similar to previous best
    ## Run 437 stress 0.0004870758 
    ## ... Procrustes: rmse 0.01610788  max resid 0.02268921 
    ## Run 438 stress 0.001160908 
    ## Run 439 stress 8.465301e-05 
    ## ... Procrustes: rmse 0.0001538465  max resid 0.0002911633 
    ## ... Similar to previous best
    ## Run 440 stress 6.790397e-05 
    ## ... Procrustes: rmse 0.0008546822  max resid 0.001594745 
    ## ... Similar to previous best
    ## Run 441 stress 0.0004426504 
    ## ... Procrustes: rmse 0.01539663  max resid 0.02170768 
    ## Run 442 stress 9.579313e-05 
    ## ... Procrustes: rmse 0.0077021  max resid 0.01055098 
    ## Run 443 stress 0.001001011 
    ## Run 444 stress 0.0004546699 
    ## ... Procrustes: rmse 0.01560422  max resid 0.02199351 
    ## Run 445 stress 0.0004487775 
    ## ... Procrustes: rmse 0.01549652  max resid 0.02184443 
    ## Run 446 stress 0.001136673 
    ## Run 447 stress 9.283832e-05 
    ## ... Procrustes: rmse 0.0002238593  max resid 0.0004600915 
    ## ... Similar to previous best
    ## Run 448 stress 0.0003059276 
    ## ... Procrustes: rmse 0.02089899  max resid 0.02881847 
    ## Run 449 stress 0.001355242 
    ## Run 450 stress 0.001359837 
    ## Run 451 stress 7.761507e-05 
    ## ... Procrustes: rmse 0.0001882458  max resid 0.000294189 
    ## ... Similar to previous best
    ## Run 452 stress 0.001253959 
    ## Run 453 stress 9.614233e-05 
    ## ... Procrustes: rmse 0.005137119  max resid 0.007560495 
    ## Run 454 stress 8.994664e-05 
    ## ... Procrustes: rmse 0.0001501513  max resid 0.0003003267 
    ## ... Similar to previous best
    ## Run 455 stress 0.2842805 
    ## Run 456 stress 9.019983e-05 
    ## ... Procrustes: rmse 0.0001914882  max resid 0.000395007 
    ## ... Similar to previous best
    ## Run 457 stress 0.001239418 
    ## Run 458 stress 8.620133e-05 
    ## ... Procrustes: rmse 0.0002251706  max resid 0.000489601 
    ## ... Similar to previous best
    ## Run 459 stress 0.1776011 
    ## Run 460 stress 0.177601 
    ## Run 461 stress 9.848145e-05 
    ## ... Procrustes: rmse 0.0002552959  max resid 0.0003877223 
    ## ... Similar to previous best
    ## Run 462 stress 0.001283521 
    ## Run 463 stress 0.1491957 
    ## Run 464 stress 0.001144293 
    ## Run 465 stress 0.2509018 
    ## Run 466 stress 0.0010726 
    ## Run 467 stress 0.1776011 
    ## Run 468 stress 8.787064e-05 
    ## ... Procrustes: rmse 0.0001511305  max resid 0.0003025702 
    ## ... Similar to previous best
    ## Run 469 stress 0.001386708 
    ## Run 470 stress 9.764084e-05 
    ## ... Procrustes: rmse 0.0002566658  max resid 0.0005109681 
    ## ... Similar to previous best
    ## Run 471 stress 0.001354871 
    ## Run 472 stress 0.1491957 
    ## Run 473 stress 0.2494706 
    ## Run 474 stress 0.2373687 
    ## Run 475 stress 9.351753e-05 
    ## ... Procrustes: rmse 0.0002765981  max resid 0.0004526614 
    ## ... Similar to previous best
    ## Run 476 stress 0.0005122485 
    ## ... Procrustes: rmse 0.01656557  max resid 0.0233201 
    ## Run 477 stress 0.00110187 
    ## Run 478 stress 0.001408881 
    ## Run 479 stress 9.825589e-05 
    ## ... Procrustes: rmse 0.000216757  max resid 0.0005001045 
    ## ... Similar to previous best
    ## Run 480 stress 9.997437e-05 
    ## ... Procrustes: rmse 0.007297066  max resid 0.01052768 
    ## Run 481 stress 8.880811e-05 
    ## ... Procrustes: rmse 0.000154997  max resid 0.0002941248 
    ## ... Similar to previous best
    ## Run 482 stress 0.001372135 
    ## Run 483 stress 8.568316e-05 
    ## ... Procrustes: rmse 0.002498428  max resid 0.003909312 
    ## ... Similar to previous best
    ## Run 484 stress 0.0004415108 
    ## ... Procrustes: rmse 0.01536269  max resid 0.0216607 
    ## Run 485 stress 8.678046e-05 
    ## ... Procrustes: rmse 0.0002270361  max resid 0.0004570055 
    ## ... Similar to previous best
    ## Run 486 stress 0.001354251 
    ## Run 487 stress 8.983715e-05 
    ## ... Procrustes: rmse 0.0005929931  max resid 0.001230797 
    ## ... Similar to previous best
    ## Run 488 stress 0.0005920803 
    ## Run 489 stress 9.49707e-05 
    ## ... Procrustes: rmse 0.0002780426  max resid 0.0005305866 
    ## ... Similar to previous best
    ## Run 490 stress 0.2570113 
    ## Run 491 stress 0.001159566 
    ## Run 492 stress 0.001336881 
    ## Run 493 stress 9.77521e-05 
    ## ... Procrustes: rmse 0.0002529268  max resid 0.0005005927 
    ## ... Similar to previous best
    ## Run 494 stress 0.0004958572 
    ## ... Procrustes: rmse 0.01629599  max resid 0.02294772 
    ## Run 495 stress 0.001251824 
    ## Run 496 stress 8.171629e-05 
    ## ... Procrustes: rmse 0.0001945964  max resid 0.0003664614 
    ## ... Similar to previous best
    ## Run 497 stress 0.1990774 
    ## Run 498 stress 0.1491957 
    ## Run 499 stress 7.383914e-05 
    ## ... Procrustes: rmse 0.0002096645  max resid 0.0003997575 
    ## ... Similar to previous best
    ## Run 500 stress 0.001202987 
    ## *** Best solution repeated 97 times

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

### SD beta with env and geo correlated variables using envfit

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
    ## env[surveyed_sites_env, c(34)]  0.99783 -0.06581 0.7388  0.029 *
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
    ## env[surveyed_sites_env, c(34)]  0.99783 -0.06581 0.7388  0.029 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by Site_type
(SD_beta_env_A_ef <- envfit(SD_beta_env_NMDS, env[surveyed_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.48817 -0.87275 0.0585  0.796  
    ## salinity_median     0.99783 -0.06581 0.7388  0.040 *
    ## oxygen_median       0.79324  0.60891 0.5361  0.810  
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
    ## temperature_median -0.48817 -0.87275 0.0585   1.00
    ## salinity_median     0.99783 -0.06581 0.7388   0.12
    ## oxygen_median       0.79324  0.60891 0.5361   1.00
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
    ## temperature_median -0.56548 -0.82476 0.0778  0.812  
    ## salinity_median     0.99173  0.12831 0.7458  0.033 *
    ## oxygen_median       0.94285  0.33321 0.3552  0.913  
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
    ## temperature_median -0.56548 -0.82476 0.0778  1.000  
    ## salinity_median     0.99173  0.12831 0.7458  0.099 .
    ## oxygen_median       0.94285  0.33321 0.3552  1.000  
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
    ## temperature_median -0.10614  0.99435 0.0705  0.699  
    ## salinity_median    -0.78732 -0.61654 0.7842  0.011 *
    ## oxygen_median      -0.96712  0.25431 0.7188  0.044 *
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
    ## salinity_median    -0.78732 -0.61654 0.7842  0.033 *
    ## oxygen_median      -0.96712  0.25431 0.7188  0.132  
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
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## temperature_median  0.0002052 -1.0000000 0.0842  0.554
    ## salinity_median    -0.0008765  1.0000000 0.5308  0.403
    ## oxygen_median      -0.0101115  0.9999500 0.7117  0.924
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
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## temperature_median  0.0002052 -1.0000000 0.0842      1
    ## salinity_median    -0.0008765  1.0000000 0.5308      1
    ## oxygen_median      -0.0101115  0.9999500 0.7117      1
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
    ## temperature_median -0.00002655 -1.00000000 0.3848  0.303
    ## salinity_median     0.00061114 -1.00000000 0.3915  0.311
    ## oxygen_median       0.00152040 -1.00000000 0.5842  0.123
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
    ## temperature_median -0.00002655 -1.00000000 0.3848  0.909
    ## salinity_median     0.00061114 -1.00000000 0.3915  0.933
    ## oxygen_median       0.00152040 -1.00000000 0.5842  0.369
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
    ## temperature_median       -0.35881 -0.93341 0.1964  0.606
    ## salinity_median          -0.70036  0.71379 0.5991  0.112
    ## oxygen_median             0.98543 -0.17006 0.5185  0.157
    ## distance_to_ocean_mean_m  0.01326 -0.99991 0.1726  0.611
    ## max_depth                -0.49998 -0.86604 0.0383  0.901
    ## logArea                  -0.26513 -0.96421 0.2144  0.585
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
    ## salinity_median          -0.70036  0.71379 0.5991  0.672
    ## oxygen_median             0.98543 -0.17006 0.5185  0.942
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
    ## distance_to_ocean_min_m -0.99979  0.02070 0.6139  0.118
    ## max_depth               -0.17731 -0.98416 0.1184  0.784
    ## logArea                  0.26569 -0.96406 0.2080  0.113
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
    ## distance_to_ocean_min_m -0.99979  0.02070 0.6139  0.354
    ## max_depth               -0.17731 -0.98416 0.1184  1.000
    ## logArea                  0.26569 -0.96406 0.2080  0.339
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by Site_type
(SD_beta_geo_A_ef <- envfit(SD_beta_geo_NMDS, env[surveyed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,19]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.99979  0.02070 0.6139  0.152
    ## max_depth               -0.17731 -0.98416 0.1184  0.778
    ## logArea                  0.26569 -0.96406 0.2080  0.119
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
    ## distance_to_ocean_min_m -0.99979  0.02070 0.6139  0.456
    ## max_depth               -0.17731 -0.98416 0.1184  1.000
    ## logArea                  0.26569 -0.96406 0.2080  0.357
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
    ## distance_to_ocean_min_m -0.94033  0.34026 0.5171  0.167
    ## max_depth               -0.21481 -0.97665 0.1851  0.662
    ## logArea                 -0.06539  0.99786 0.0118  0.955
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
    ## distance_to_ocean_min_m -0.94033  0.34026 0.5171  0.501
    ## max_depth               -0.21481 -0.97665 0.1851  1.000
    ## logArea                 -0.06539  0.99786 0.0118  1.000
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
    ## distance_to_ocean_min_m  0.86175 -0.50733 0.3198  0.070 .
    ## max_depth               -0.68923 -0.72454 0.5689  0.023 *
    ## logArea                 -0.84908  0.52827 0.2993  0.206  
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
    ## distance_to_ocean_min_m  0.86175 -0.50733 0.3198  0.210  
    ## max_depth               -0.68923 -0.72454 0.5689  0.069 .
    ## logArea                 -0.84908  0.52827 0.2993  0.618  
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
    ## distance_to_ocean_min_m  0.86963  0.49371 0.6793  0.205
    ## max_depth                0.68533  0.72824 0.0510  0.973
    ## logArea                 -0.52042  0.85391 0.2187  0.523
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
    ## distance_to_ocean_min_m  0.86963  0.49371 0.6793  0.615
    ## max_depth                0.68533  0.72824 0.0510  1.000
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
    ##                               NMDS1       NMDS2     r2 Pr(>r)   
    ## distance_to_ocean_min_m -0.00068819  1.00000000 0.6622  0.079 . 
    ## max_depth                0.00177299 -1.00000000 0.8227  0.015 * 
    ## logArea                  0.00144800  1.00000000 0.8168  0.007 **
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
    ##                               NMDS1       NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.00068819  1.00000000 0.6622  0.237  
    ## max_depth                0.00177299 -1.00000000 0.8227  0.045 *
    ## logArea                  0.00144800  1.00000000 0.8168  0.021 *
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
    ##       Significance: 0.239 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.289 0.317 0.347 0.367 
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
    ## 0.576 0.601 0.627 0.661 
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
    ##       Significance: 0.466 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.593 0.628 0.664 0.688 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_mant_pv <- rbind(SD_beta_env_mant_t$signif, SD_beta_env_mant_s$signif, SD_beta_env_mant_o$signif)
SD_beta_env_mant_pv <- SD_beta_env_mant_pv[,1]
(SD_beta_env_mant_pv <- p.adjust(SD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.717 0.018 1.000

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
    ##       Significance: 0.331 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.331 0.365 0.404 0.432 
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
    ##       Significance: 0.012 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.567 0.601 0.640 0.671 
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
    ##       Significance: 0.605 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.491 0.525 0.569 0.620 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_MS_mant_pv <- rbind(SD_beta_env_MS_mant_t$signif, SD_beta_env_MS_mant_s$signif, SD_beta_env_MS_mant_o$signif)
SD_beta_env_MS_mant_pv <- SD_beta_env_MS_mant_pv[,1]
(SD_beta_env_MS_mant_pv <- p.adjust(SD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 0.993 0.036 1.000

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
    ##       Significance: 0.76 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.131 0.171 0.225 0.400 
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
    ##       Significance: 0.038 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.314 0.396 0.466 0.541 
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
    ##       Significance: 0.023 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.319 0.443 0.505 0.626 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_OM_mant_pv <- rbind(SD_beta_env_OM_mant_t$signif, SD_beta_env_OM_mant_s$signif, SD_beta_env_OM_mant_o$signif)
SD_beta_env_OM_mant_pv <- SD_beta_env_OM_mant_pv[,1]
(SD_beta_env_OM_mant_pv <- p.adjust(SD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.114 0.069

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
    ##       Significance: 0.214 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0254 0.0516 0.0655 0.0803 
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
    ##       Significance: 0.01 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.467 0.503 0.534 0.552 
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
    ##       Significance: 0.38 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.707 0.735 0.753 0.774 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_SO_mant_pv <- rbind(SD_beta_env_SO_mant_t$signif, SD_beta_env_SO_mant_s$signif, SD_beta_env_SO_mant_o$signif)
SD_beta_env_SO_mant_pv <- SD_beta_env_SO_mant_pv[,1]
(SD_beta_env_SO_mant_pv <- p.adjust(SD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 0.642 0.030 1.000

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
    ##       Significance: 0.771 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.219 0.320 0.418 0.678 
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
    ##       Significance: 0.186 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.284 0.343 0.397 0.437 
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
    ##       Significance: 0.102 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.281 0.436 0.517 0.623 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_M_mant_pv <- rbind(SD_beta_env_M_mant_t$signif, SD_beta_env_M_mant_s$signif, SD_beta_env_M_mant_o$signif)
SD_beta_env_M_mant_pv <- SD_beta_env_M_mant_pv[,1]
(SD_beta_env_M_mant_pv <- p.adjust(SD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.558 0.306

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
    ##       Significance: 0.235 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.225 0.290 0.331 0.374 
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
    ##       Significance: 0.485 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.548 0.577 0.600 0.622 
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
    ##       Significance: 0.64 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.197 0.227 0.258 0.283 
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
    ##       Significance: 0.394 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0793 0.0932 0.1044 0.1114 
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
    ##       Significance: 0.409 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.433 0.470 0.511 0.565 
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
    ##       Significance: 0.767 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.216 0.257 0.282 0.316 
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
    ##       Significance: 0.765 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0705 0.0914 0.1113 0.1370 
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
    ##       Significance: 0.618 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.245 0.281 0.308 0.355 
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
    ##       Significance: 0.06 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.174 0.244 0.317 0.369 
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
    ##       Significance: 0.089 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.195 0.241 0.277 0.311 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_OM_mant_pv <- rbind( SD_beta_geo_OM_mant_dmean$signif, SD_beta_geo_OM_mant_md$signif, SD_beta_geo_OM_mant_la$signif)
SD_beta_geo_OM_mant_pv <- SD_beta_geo_OM_mant_pv[,1]
(SD_beta_geo_OM_mant_pv <- p.adjust(SD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.180 0.267

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
    ##       Significance: 0.475 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.684 0.709 0.732 0.744 
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
    ##       Significance: 0.6 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.141 0.182 0.212 0.244 
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
    ##       Significance: 0.232 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.273 0.290 0.302 0.314 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_SO_mant_pv <- rbind( SD_beta_geo_SO_mant_dmean$signif, SD_beta_geo_SO_mant_md$signif, SD_beta_geo_SO_mant_la$signif)
SD_beta_geo_SO_mant_pv <- SD_beta_geo_SO_mant_pv[,1]
(SD_beta_geo_SO_mant_pv <- p.adjust(SD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 1.000 0.696

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
    ##       Significance: 0.489 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.233 0.363 0.476 0.762 
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
    ## 0.260 0.411 0.510 0.573 
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
    ##       Significance: 0.041 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.262 0.383 0.469 0.550 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_M_mant_pv <- rbind( SD_beta_geo_M_mant_dmean$signif, SD_beta_geo_M_mant_md$signif, SD_beta_geo_M_mant_la$signif)
SD_beta_geo_M_mant_pv <- SD_beta_geo_M_mant_pv[,1]
(SD_beta_geo_M_mant_pv <- p.adjust(SD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.006 0.123

### SD beta NMDS ordination plots

``` r
# SD beta total NMDS scores
SD_beta_NMDS_data.scores <- as.data.frame(scores(SD_beta_NMDS))
SD_beta_NMDS_data.scores$Site_type <- env[surveyed_sites,19]
SD_beta_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
SD_beta_NMDS_data.scores$Site_type <- factor(SD_beta_NMDS_data.scores$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

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
SD_beta_ef_plot <- ggplot(data = SD_beta_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = SD_beta_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
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
SD_beta_ef_plot_CI90 <- ggplot(data = SD_beta_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.90) +
  geom_point(data = SD_beta_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
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
SD_beta_rep_NMDS_data.scores$Site_type <- env[surveyed_sites,19]
SD_beta_rep_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
SD_beta_rep_NMDS_data.scores$Site_type <- factor(SD_beta_rep_NMDS_data.scores$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of SD beta replacement with CI = 0.95
SD_beta_rep_plot <- ggplot(data = SD_beta_rep_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = SD_beta_rep_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
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
SD_beta_ric_NMDS_data.scores$Site_type <- env[surveyed_sites,19]
SD_beta_ric_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
SD_beta_ric_NMDS_data.scores$Site_type <- factor(SD_beta_ric_NMDS_data.scores$Site_type, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of SD beta richness with CI = 0.95
SD_beta_ric_plot <- ggplot(data = SD_beta_ric_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = SD_beta_ric_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
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
SD_beta_ref_NMDS_data.scores$Site_type <- env[,19]
SD_beta_ref_NMDS_data.scores$Lakes <- env[,1]
SD_beta_ref_NMDS_data.scores$Site_type <- factor(SD_beta_ref_NMDS_data.scores$Site_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of SD beta ref with CI = 0.95
SD_beta_ref_plot <- ggplot(data = SD_beta_ref_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Site_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Site_type), type = "t", level = 0.95) +
  geom_point(data = SD_beta_ref_NMDS_data.scores, aes(color = Site_type, fill = Site_type), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = SD_beta_ref_NMDS_data.scores, label = SD_beta_ref_NMDS_data.scores$Lakes, 
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
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/SD_beta_ref_NMDS.jpg", SD_beta_ref_plot, width = 6.26, height = 6, units = "in")
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

### SD beta NMDS and distance outliers

``` r
# Determine NMDS1 outliers
ordered_SD_beta_NMDS1_data.scores <- SD_beta_NMDS_data.scores[order(SD_beta_NMDS_data.scores$NMDS1), ]

outlier_SD_beta_NMDS1 <- ordered_SD_beta_NMDS1_data.scores %>%
  group_by(Site_type) %>%
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
(outlier_SD_beta_NMDS1_plot <- ggplot(outlier_SD_beta_NMDS1, aes(x = Site_type, y = NMDS1, fill = Site_type)) +
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

![](SD_analyses_files/figure-gfm/SD%20beta%20outliers-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/SD/outlier_SD_beta_NMDS1.jpg", outlier_SD_beta_NMDS1_plot, width = 6.26, height = 6, units = "in")


# Determine NMDS2 outliers
ordered_SD_beta_NMDS2_data.scores <- SD_beta_NMDS_data.scores[order(SD_beta_NMDS_data.scores$NMDS2), ]

outlier_SD_beta_NMDS2 <- ordered_SD_beta_NMDS2_data.scores %>%
  group_by(Site_type) %>%
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
(outlier_SD_beta_NMDS2_plot <- ggplot(outlier_SD_beta_NMDS2, aes(x = Site_type, y = NMDS2, fill = Site_type)) +
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
SD_beta_mean_dist <- aggregate(SD_beta_dist_matrix, by = list(Site_type_group), FUN = mean, na.rm = TRUE)

# Rename rows, delete column, and transpose matrix and make it a data frame
row.names(SD_beta_mean_dist) <- SD_beta_mean_dist$Group.1
SD_beta_mean_dist <- SD_beta_mean_dist[,-1] 
SD_beta_mean_dist <- as.data.frame(t(SD_beta_mean_dist))

# Add Site_type column
SD_beta_mean_dist$Site_type <- env[surveyed_sites,19]
SD_beta_mean_dist$Site_type <- factor(SD_beta_mean_dist$Site_type, levels = c("Ocean", "Mixed", "Stratified"))


# Determine ocean sites outliers based on distance from other ocean sites
outlier_SD_beta_mean_dist_ocean <- SD_beta_mean_dist[ocean_sites,] %>%
  group_by(Site_type) %>%
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
  group_by(Site_type) %>%
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
  group_by(Site_type) %>%
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

(outlier_SD_beta_mean_dist_plot <- ggplot(outlier_SD_beta_mean_dist, aes(x = Site_type, y = Distances, fill = Site_type)) +
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
