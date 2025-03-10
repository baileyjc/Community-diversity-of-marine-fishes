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

### SD alpha outliers for each site type

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

### SD alpha & site type

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

### SD alpha with env and geo linear models & ANOVAs

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

### SD alpha with env and geo correlated variables

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
    ## 3 Stratified vs Reference  1 0.6261283 1.909357 0.2143091   0.121      0.726
    ## 4          Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.043      0.258
    ## 5      Mixed vs Reference  1 0.5979203 1.966332 0.2193017   0.104      0.624
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
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.041      0.123

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
    ## 3      Mixed vs Ocean  1 0.6396174 2.135322 0.1625633   0.004      0.012   .

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
    ## 1 Stratified vs Mixed  1 1.4219030 4.731564 0.2827927   0.002      0.006   *
    ## 2 Stratified vs Ocean  1 1.3260066 4.147587 0.2931657   0.006      0.018   .
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.068      0.204

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
    ## Run 1 stress 0.1051402 
    ## ... Procrustes: rmse 0.02306848  max resid 0.0784661 
    ## Run 2 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 0.006704424  max resid 0.02560846 
    ## Run 3 stress 0.1049903 
    ## ... Procrustes: rmse 0.006753004  max resid 0.02588609 
    ## Run 4 stress 0.1135598 
    ## Run 5 stress 0.113594 
    ## Run 6 stress 0.1084964 
    ## Run 7 stress 0.1084964 
    ## Run 8 stress 0.1079382 
    ## Run 9 stress 0.1050144 
    ## ... Procrustes: rmse 0.009670564  max resid 0.0261064 
    ## Run 10 stress 0.1051924 
    ## ... Procrustes: rmse 0.02283927  max resid 0.07721692 
    ## Run 11 stress 0.1084964 
    ## Run 12 stress 0.1079382 
    ## Run 13 stress 0.1063129 
    ## Run 14 stress 0.1138394 
    ## Run 15 stress 0.106313 
    ## Run 16 stress 0.1088547 
    ## Run 17 stress 0.1309383 
    ## Run 18 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179807  max resid 0.07684065 
    ## Run 19 stress 0.1082722 
    ## Run 20 stress 0.1050144 
    ## ... Procrustes: rmse 0.009670469  max resid 0.02610666 
    ## Run 21 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179705  max resid 0.07683455 
    ## Run 22 stress 0.1165301 
    ## Run 23 stress 0.1051923 
    ## ... Procrustes: rmse 0.022837  max resid 0.07728249 
    ## Run 24 stress 0.10835 
    ## Run 25 stress 0.1138391 
    ## Run 26 stress 0.1082722 
    ## Run 27 stress 0.1135597 
    ## Run 28 stress 0.1160572 
    ## Run 29 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179883  max resid 0.07683994 
    ## Run 30 stress 0.1084965 
    ## Run 31 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283639  max resid 0.0772341 
    ## Run 32 stress 0.1084964 
    ## Run 33 stress 0.1083503 
    ## Run 34 stress 0.1063128 
    ## Run 35 stress 0.1083992 
    ## Run 36 stress 0.1084965 
    ## Run 37 stress 0.1088547 
    ## Run 38 stress 0.1084965 
    ## Run 39 stress 0.1079382 
    ## Run 40 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179764  max resid 0.07683638 
    ## Run 41 stress 0.1058706 
    ## Run 42 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283634  max resid 0.07726923 
    ## Run 43 stress 0.10835 
    ## Run 44 stress 0.1135599 
    ## Run 45 stress 0.1079384 
    ## Run 46 stress 0.1058423 
    ## Run 47 stress 0.1063128 
    ## Run 48 stress 0.1049903 
    ## ... Procrustes: rmse 0.0067341  max resid 0.02580083 
    ## Run 49 stress 0.1174541 
    ## Run 50 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180227  max resid 0.0768528 
    ## Run 51 stress 0.1088547 
    ## Run 52 stress 0.1136647 
    ## Run 53 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179228  max resid 0.07681658 
    ## Run 54 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.834503e-05  max resid 5.379292e-05 
    ## ... Similar to previous best
    ## Run 55 stress 0.1135598 
    ## Run 56 stress 0.1084964 
    ## Run 57 stress 0.1083499 
    ## Run 58 stress 0.1136646 
    ## Run 59 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283953  max resid 0.07729072 
    ## Run 60 stress 0.1050144 
    ## ... Procrustes: rmse 0.009670585  max resid 0.02611666 
    ## Run 61 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179961  max resid 0.07684819 
    ## Run 62 stress 0.1063128 
    ## Run 63 stress 0.1160548 
    ## Run 64 stress 0.1135598 
    ## Run 65 stress 0.1049903 
    ## ... Procrustes: rmse 0.006806239  max resid 0.02616487 
    ## Run 66 stress 0.1165311 
    ## Run 67 stress 0.1050144 
    ## ... Procrustes: rmse 0.009670649  max resid 0.02613296 
    ## Run 68 stress 0.1083499 
    ## Run 69 stress 0.1135596 
    ## Run 70 stress 0.1135597 
    ## Run 71 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283378  max resid 0.07728365 
    ## Run 72 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180056  max resid 0.07685662 
    ## Run 73 stress 0.1138395 
    ## Run 74 stress 0.1063128 
    ## Run 75 stress 0.1165296 
    ## Run 76 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228364  max resid 0.07727486 
    ## Run 77 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283802  max resid 0.07729571 
    ## Run 78 stress 0.1309387 
    ## Run 79 stress 0.1058423 
    ## Run 80 stress 0.1135598 
    ## Run 81 stress 0.10835 
    ## Run 82 stress 0.1084964 
    ## Run 83 stress 0.1088547 
    ## Run 84 stress 0.10835 
    ## Run 85 stress 0.1160538 
    ## Run 86 stress 0.1063128 
    ## Run 87 stress 0.1138397 
    ## Run 88 stress 0.1138394 
    ## Run 89 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179773  max resid 0.07684207 
    ## Run 90 stress 0.1058423 
    ## Run 91 stress 0.1063128 
    ## Run 92 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180533  max resid 0.07687297 
    ## Run 93 stress 0.1165301 
    ## Run 94 stress 0.1049903 
    ## ... Procrustes: rmse 0.006730759  max resid 0.02581232 
    ## Run 95 stress 0.1160543 
    ## Run 96 stress 0.1051924 
    ## ... Procrustes: rmse 0.02283861  max resid 0.0772838 
    ## Run 97 stress 0.1049903 
    ## ... Procrustes: rmse 0.006733725  max resid 0.02581908 
    ## Run 98 stress 0.1135597 
    ## Run 99 stress 0.1050338 
    ## ... Procrustes: rmse 0.006253721  max resid 0.02184196 
    ## Run 100 stress 0.1084964 
    ## Run 101 stress 0.1136646 
    ## Run 102 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228399  max resid 0.07726715 
    ## Run 103 stress 0.1136646 
    ## Run 104 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284021  max resid 0.07731139 
    ## Run 105 stress 0.1058423 
    ## Run 106 stress 0.1088547 
    ## Run 107 stress 0.107938 
    ## Run 108 stress 0.1063128 
    ## Run 109 stress 0.1084965 
    ## Run 110 stress 0.1136646 
    ## Run 111 stress 0.1138392 
    ## Run 112 stress 0.1063128 
    ## Run 113 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283957  max resid 0.07727297 
    ## Run 114 stress 0.1084965 
    ## Run 115 stress 0.113594 
    ## Run 116 stress 0.1083996 
    ## Run 117 stress 0.1050144 
    ## ... Procrustes: rmse 0.009669495  max resid 0.02613537 
    ## Run 118 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218033  max resid 0.07687519 
    ## Run 119 stress 0.109799 
    ## Run 120 stress 0.1084966 
    ## Run 121 stress 0.1063128 
    ## Run 122 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283874  max resid 0.07731349 
    ## Run 123 stress 0.1050144 
    ## ... Procrustes: rmse 0.009653657  max resid 0.02602546 
    ## Run 124 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283784  max resid 0.0772857 
    ## Run 125 stress 0.1049903 
    ## ... Procrustes: rmse 0.006734461  max resid 0.02582733 
    ## Run 126 stress 0.1135596 
    ## Run 127 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283745  max resid 0.07730357 
    ## Run 128 stress 0.113594 
    ## Run 129 stress 0.108399 
    ## Run 130 stress 0.1138391 
    ## Run 131 stress 0.1063128 
    ## Run 132 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283781  max resid 0.07730313 
    ## Run 133 stress 0.1136646 
    ## Run 134 stress 0.1088547 
    ## Run 135 stress 0.1088547 
    ## Run 136 stress 0.1102721 
    ## Run 137 stress 0.1082722 
    ## Run 138 stress 0.116055 
    ## Run 139 stress 0.1049903 
    ## ... Procrustes: rmse 0.006718159  max resid 0.02575312 
    ## Run 140 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283832  max resid 0.07729884 
    ## Run 141 stress 0.1079381 
    ## Run 142 stress 0.1102723 
    ## Run 143 stress 0.1084965 
    ## Run 144 stress 0.1084965 
    ## Run 145 stress 0.1049903 
    ## ... Procrustes: rmse 0.006701041  max resid 0.0256749 
    ## Run 146 stress 0.1050144 
    ## ... Procrustes: rmse 0.009677336  max resid 0.02612533 
    ## Run 147 stress 0.1097989 
    ## Run 148 stress 0.1049903 
    ## ... Procrustes: rmse 0.006720292  max resid 0.025761 
    ## Run 149 stress 0.11356 
    ## Run 150 stress 0.1097987 
    ## Run 151 stress 0.1049903 
    ## ... Procrustes: rmse 0.006767049  max resid 0.02597724 
    ## Run 152 stress 0.1079382 
    ## Run 153 stress 0.1050144 
    ## ... Procrustes: rmse 0.009664472  max resid 0.02610087 
    ## Run 154 stress 0.1082722 
    ## Run 155 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228372  max resid 0.07729805 
    ## Run 156 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228385  max resid 0.07729609 
    ## Run 157 stress 0.1165301 
    ## Run 158 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180624  max resid 0.07687123 
    ## Run 159 stress 0.1084964 
    ## Run 160 stress 0.1050144 
    ## ... Procrustes: rmse 0.009659717  max resid 0.02606611 
    ## Run 161 stress 0.1049791 
    ## ... Procrustes: rmse 8.192712e-06  max resid 1.600257e-05 
    ## ... Similar to previous best
    ## Run 162 stress 0.1088548 
    ## Run 163 stress 0.1084965 
    ## Run 164 stress 0.1135597 
    ## Run 165 stress 0.1102722 
    ## Run 166 stress 0.1079382 
    ## Run 167 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179904  max resid 0.07684609 
    ## Run 168 stress 0.113594 
    ## Run 169 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180172  max resid 0.07686116 
    ## Run 170 stress 0.1136647 
    ## Run 171 stress 0.1136646 
    ## Run 172 stress 0.1050338 
    ## ... Procrustes: rmse 0.006255986  max resid 0.02185937 
    ## Run 173 stress 0.1063128 
    ## Run 174 stress 0.1135941 
    ## Run 175 stress 0.1058423 
    ## Run 176 stress 0.1079382 
    ## Run 177 stress 0.1050144 
    ## ... Procrustes: rmse 0.009662862  max resid 0.02608238 
    ## Run 178 stress 0.1083991 
    ## Run 179 stress 0.1063128 
    ## Run 180 stress 0.1063128 
    ## Run 181 stress 0.1135597 
    ## Run 182 stress 0.1135597 
    ## Run 183 stress 0.1084965 
    ## Run 184 stress 0.1102722 
    ## Run 185 stress 0.1082722 
    ## Run 186 stress 0.1084965 
    ## Run 187 stress 0.1050144 
    ## ... Procrustes: rmse 0.00969756  max resid 0.02627074 
    ## Run 188 stress 0.1136646 
    ## Run 189 stress 0.1136646 
    ## Run 190 stress 0.1088548 
    ## Run 191 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179915  max resid 0.07683902 
    ## Run 192 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180106  max resid 0.07685441 
    ## Run 193 stress 0.1050338 
    ## ... Procrustes: rmse 0.006252978  max resid 0.02184729 
    ## Run 194 stress 0.1063128 
    ## Run 195 stress 0.1138392 
    ## Run 196 stress 0.1084965 
    ## Run 197 stress 0.113594 
    ## Run 198 stress 0.1082722 
    ## Run 199 stress 0.1063128 
    ## Run 200 stress 0.1160546 
    ## Run 201 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180003  max resid 0.07682206 
    ## Run 202 stress 0.1063128 
    ## Run 203 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283879  max resid 0.07730088 
    ## Run 204 stress 0.1063128 
    ## Run 205 stress 0.1135941 
    ## Run 206 stress 0.1138392 
    ## Run 207 stress 0.1084965 
    ## Run 208 stress 0.1049903 
    ## ... Procrustes: rmse 0.006772648  max resid 0.02600642 
    ## Run 209 stress 0.1049904 
    ## ... Procrustes: rmse 0.006791106  max resid 0.02609318 
    ## Run 210 stress 0.1088547 
    ## Run 211 stress 0.1063128 
    ## Run 212 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283821  max resid 0.07729933 
    ## Run 213 stress 0.1063128 
    ## Run 214 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180145  max resid 0.07685244 
    ## Run 215 stress 0.1102721 
    ## Run 216 stress 0.1088548 
    ## Run 217 stress 0.1084965 
    ## Run 218 stress 0.113594 
    ## Run 219 stress 0.1160542 
    ## Run 220 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283531  max resid 0.07726355 
    ## Run 221 stress 0.1136646 
    ## Run 222 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180333  max resid 0.07687337 
    ## Run 223 stress 0.1051402 
    ## ... Procrustes: rmse 0.0217997  max resid 0.0768438 
    ## Run 224 stress 0.1050144 
    ## ... Procrustes: rmse 0.009656429  max resid 0.0260564 
    ## Run 225 stress 0.1135598 
    ## Run 226 stress 0.1049903 
    ## ... Procrustes: rmse 0.00672508  max resid 0.02578434 
    ## Run 227 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283922  max resid 0.07729014 
    ## Run 228 stress 0.1063128 
    ## Run 229 stress 0.1136647 
    ## Run 230 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228406  max resid 0.07728918 
    ## Run 231 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284138  max resid 0.07723829 
    ## Run 232 stress 0.1063128 
    ## Run 233 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180169  max resid 0.07685632 
    ## Run 234 stress 0.1049903 
    ## ... Procrustes: rmse 0.006726558  max resid 0.02579199 
    ## Run 235 stress 0.1084966 
    ## Run 236 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283603  max resid 0.07726774 
    ## Run 237 stress 0.1135598 
    ## Run 238 stress 0.1136647 
    ## Run 239 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284407  max resid 0.07727735 
    ## Run 240 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180201  max resid 0.07685785 
    ## Run 241 stress 0.1084965 
    ## Run 242 stress 0.1082722 
    ## Run 243 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180269  max resid 0.07686139 
    ## Run 244 stress 0.1051403 
    ## ... Procrustes: rmse 0.02179978  max resid 0.07681505 
    ## Run 245 stress 0.1083501 
    ## Run 246 stress 0.1160548 
    ## Run 247 stress 0.1102725 
    ## Run 248 stress 0.1136646 
    ## Run 249 stress 0.1051923 
    ## ... Procrustes: rmse 0.022834  max resid 0.07728896 
    ## Run 250 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180362  max resid 0.07686225 
    ## Run 251 stress 0.1135598 
    ## Run 252 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181093  max resid 0.07688588 
    ## Run 253 stress 0.1136648 
    ## Run 254 stress 0.1051924 
    ## ... Procrustes: rmse 0.02283982  max resid 0.07725475 
    ## Run 255 stress 0.1063129 
    ## Run 256 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180618  max resid 0.07687358 
    ## Run 257 stress 0.1160594 
    ## Run 258 stress 0.1135599 
    ## Run 259 stress 0.1102721 
    ## Run 260 stress 0.1063128 
    ## Run 261 stress 0.1049903 
    ## ... Procrustes: rmse 0.006778037  max resid 0.02603004 
    ## Run 262 stress 0.1084965 
    ## Run 263 stress 0.1063128 
    ## Run 264 stress 0.1082722 
    ## Run 265 stress 0.1050144 
    ## ... Procrustes: rmse 0.009662703  max resid 0.02609515 
    ## Run 266 stress 0.1088547 
    ## Run 267 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218035  max resid 0.07686848 
    ## Run 268 stress 0.1097987 
    ## Run 269 stress 0.1084964 
    ## Run 270 stress 0.1083501 
    ## Run 271 stress 0.10835 
    ## Run 272 stress 0.1079381 
    ## Run 273 stress 0.1136647 
    ## Run 274 stress 0.10835 
    ## Run 275 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283764  max resid 0.07728639 
    ## Run 276 stress 0.1084965 
    ## Run 277 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180294  max resid 0.07687732 
    ## Run 278 stress 0.1083999 
    ## Run 279 stress 0.1136647 
    ## Run 280 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180635  max resid 0.076872 
    ## Run 281 stress 0.1082722 
    ## Run 282 stress 0.1174541 
    ## Run 283 stress 0.1136646 
    ## Run 284 stress 0.1079382 
    ## Run 285 stress 0.1160542 
    ## Run 286 stress 0.1136647 
    ## Run 287 stress 0.1160555 
    ## Run 288 stress 0.1097988 
    ## Run 289 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179964  max resid 0.07684834 
    ## Run 290 stress 0.1079384 
    ## Run 291 stress 0.1050338 
    ## ... Procrustes: rmse 0.006254991  max resid 0.0218496 
    ## Run 292 stress 0.1135597 
    ## Run 293 stress 0.1063128 
    ## Run 294 stress 0.1088547 
    ## Run 295 stress 0.113594 
    ## Run 296 stress 0.1049903 
    ## ... Procrustes: rmse 0.006708908  max resid 0.02571151 
    ## Run 297 stress 0.1058423 
    ## Run 298 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284291  max resid 0.07734176 
    ## Run 299 stress 0.1135598 
    ## Run 300 stress 0.1136647 
    ## Run 301 stress 0.1136648 
    ## Run 302 stress 0.1138394 
    ## Run 303 stress 0.1135598 
    ## Run 304 stress 0.1309381 
    ## Run 305 stress 0.1136647 
    ## Run 306 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180128  max resid 0.07685057 
    ## Run 307 stress 0.1088547 
    ## Run 308 stress 0.1050144 
    ## ... Procrustes: rmse 0.009681378  max resid 0.02617481 
    ## Run 309 stress 0.1084965 
    ## Run 310 stress 0.1050144 
    ## ... Procrustes: rmse 0.00965867  max resid 0.02605794 
    ## Run 311 stress 0.1160556 
    ## Run 312 stress 0.1049903 
    ## ... Procrustes: rmse 0.006758432  max resid 0.02593823 
    ## Run 313 stress 0.1138393 
    ## Run 314 stress 0.1102722 
    ## Run 315 stress 0.1084964 
    ## Run 316 stress 0.1058423 
    ## Run 317 stress 0.1049903 
    ## ... Procrustes: rmse 0.006719947  max resid 0.02576673 
    ## Run 318 stress 0.1079382 
    ## Run 319 stress 0.1088548 
    ## Run 320 stress 0.1138394 
    ## Run 321 stress 0.1138393 
    ## Run 322 stress 0.1051402 
    ## ... Procrustes: rmse 0.021802  max resid 0.07685645 
    ## Run 323 stress 0.1063128 
    ## Run 324 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284076  max resid 0.077287 
    ## Run 325 stress 0.1135598 
    ## Run 326 stress 0.1063128 
    ## Run 327 stress 0.10835 
    ## Run 328 stress 0.1174547 
    ## Run 329 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283921  max resid 0.07728747 
    ## Run 330 stress 0.1049791 
    ## ... Procrustes: rmse 2.898565e-05  max resid 8.461822e-05 
    ## ... Similar to previous best
    ## Run 331 stress 0.1079382 
    ## Run 332 stress 0.1136647 
    ## Run 333 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180081  max resid 0.07685367 
    ## Run 334 stress 0.1050338 
    ## ... Procrustes: rmse 0.006257979  max resid 0.02186938 
    ## Run 335 stress 0.1050144 
    ## ... Procrustes: rmse 0.00966562  max resid 0.02609802 
    ## Run 336 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283822  max resid 0.07729165 
    ## Run 337 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283556  max resid 0.07730108 
    ## Run 338 stress 0.1082723 
    ## Run 339 stress 0.1174531 
    ## Run 340 stress 0.1082723 
    ## Run 341 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180145  max resid 0.07685473 
    ## Run 342 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180242  max resid 0.07685761 
    ## Run 343 stress 0.1082723 
    ## Run 344 stress 0.1050144 
    ## ... Procrustes: rmse 0.009665878  max resid 0.02611084 
    ## Run 345 stress 0.1138397 
    ## Run 346 stress 0.1083999 
    ## Run 347 stress 0.1084966 
    ## Run 348 stress 0.1050144 
    ## ... Procrustes: rmse 0.00966087  max resid 0.02607637 
    ## Run 349 stress 0.113594 
    ## Run 350 stress 0.1136647 
    ## Run 351 stress 0.1138986 
    ## Run 352 stress 0.1058423 
    ## Run 353 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180097  max resid 0.07684197 
    ## Run 354 stress 0.1058423 
    ## Run 355 stress 0.1084965 
    ## Run 356 stress 0.1160544 
    ## Run 357 stress 0.1135597 
    ## Run 358 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283596  max resid 0.07729363 
    ## Run 359 stress 0.1063128 
    ## Run 360 stress 0.1063128 
    ## Run 361 stress 0.113594 
    ## Run 362 stress 0.1136646 
    ## Run 363 stress 0.1084965 
    ## Run 364 stress 0.1051403 
    ## ... Procrustes: rmse 0.02179716  max resid 0.07680828 
    ## Run 365 stress 0.1083503 
    ## Run 366 stress 0.1050144 
    ## ... Procrustes: rmse 0.009659084  max resid 0.0260606 
    ## Run 367 stress 0.1082722 
    ## Run 368 stress 0.1083996 
    ## Run 369 stress 0.1136647 
    ## Run 370 stress 0.1138393 
    ## Run 371 stress 0.1165304 
    ## Run 372 stress 0.1135597 
    ## Run 373 stress 0.10835 
    ## Run 374 stress 0.1082722 
    ## Run 375 stress 0.1050338 
    ## ... Procrustes: rmse 0.006252142  max resid 0.02183516 
    ## Run 376 stress 0.1097989 
    ## Run 377 stress 0.1082722 
    ## Run 378 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218007  max resid 0.07685308 
    ## Run 379 stress 0.1049791 
    ## ... Procrustes: rmse 3.474883e-05  max resid 0.0001104696 
    ## ... Similar to previous best
    ## Run 380 stress 0.1049792 
    ## ... Procrustes: rmse 2.563025e-05  max resid 8.986985e-05 
    ## ... Similar to previous best
    ## Run 381 stress 0.1084965 
    ## Run 382 stress 0.1050338 
    ## ... Procrustes: rmse 0.006246801  max resid 0.02182938 
    ## Run 383 stress 0.1050144 
    ## ... Procrustes: rmse 0.009685844  max resid 0.0262297 
    ## Run 384 stress 0.1063128 
    ## Run 385 stress 0.1309381 
    ## Run 386 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179983  max resid 0.07685083 
    ## Run 387 stress 0.1136646 
    ## Run 388 stress 0.1138397 
    ## Run 389 stress 0.113594 
    ## Run 390 stress 0.1063128 
    ## Run 391 stress 0.1082722 
    ## Run 392 stress 0.1049791 
    ## ... Procrustes: rmse 2.035327e-05  max resid 6.299709e-05 
    ## ... Similar to previous best
    ## Run 393 stress 0.1136646 
    ## Run 394 stress 0.1050144 
    ## ... Procrustes: rmse 0.009667684  max resid 0.02613524 
    ## Run 395 stress 0.1136647 
    ## Run 396 stress 0.1084964 
    ## Run 397 stress 0.1063128 
    ## Run 398 stress 0.1063128 
    ## Run 399 stress 0.1083501 
    ## Run 400 stress 0.1051403 
    ## ... Procrustes: rmse 0.02179234  max resid 0.07682599 
    ## Run 401 stress 0.1136646 
    ## Run 402 stress 0.1079381 
    ## Run 403 stress 0.1097987 
    ## Run 404 stress 0.1084965 
    ## Run 405 stress 0.1135597 
    ## Run 406 stress 0.1063129 
    ## Run 407 stress 0.1138394 
    ## Run 408 stress 0.1088548 
    ## Run 409 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283778  max resid 0.07729993 
    ## Run 410 stress 0.1097989 
    ## Run 411 stress 0.1049903 
    ## ... Procrustes: rmse 0.006688667  max resid 0.02561864 
    ## Run 412 stress 0.1050144 
    ## ... Procrustes: rmse 0.009668509  max resid 0.02612679 
    ## Run 413 stress 0.1097987 
    ## Run 414 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180212  max resid 0.07686432 
    ## Run 415 stress 0.1049903 
    ## ... Procrustes: rmse 0.006767158  max resid 0.02597885 
    ## Run 416 stress 0.1063128 
    ## Run 417 stress 0.1063128 
    ## Run 418 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284141  max resid 0.07733332 
    ## Run 419 stress 0.1084964 
    ## Run 420 stress 0.10835 
    ## Run 421 stress 0.1135597 
    ## Run 422 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283816  max resid 0.07729411 
    ## Run 423 stress 0.1084966 
    ## Run 424 stress 0.1051924 
    ## ... Procrustes: rmse 0.02283736  max resid 0.0773008 
    ## Run 425 stress 0.1051924 
    ## ... Procrustes: rmse 0.02283868  max resid 0.07734744 
    ## Run 426 stress 0.10835 
    ## Run 427 stress 0.1084964 
    ## Run 428 stress 0.1136646 
    ## Run 429 stress 0.1083502 
    ## Run 430 stress 0.1084964 
    ## Run 431 stress 0.1138393 
    ## Run 432 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 1.419122e-06  max resid 3.384575e-06 
    ## ... Similar to previous best
    ## Run 433 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180108  max resid 0.07685007 
    ## Run 434 stress 0.1102724 
    ## Run 435 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283809  max resid 0.07728965 
    ## Run 436 stress 0.1079382 
    ## Run 437 stress 0.1050144 
    ## ... Procrustes: rmse 0.009650618  max resid 0.02602701 
    ## Run 438 stress 0.1083501 
    ## Run 439 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180099  max resid 0.07685485 
    ## Run 440 stress 0.1063128 
    ## Run 441 stress 0.1082722 
    ## Run 442 stress 0.1050144 
    ## ... Procrustes: rmse 0.009674232  max resid 0.02612797 
    ## Run 443 stress 0.1084966 
    ## Run 444 stress 0.1050144 
    ## ... Procrustes: rmse 0.009666509  max resid 0.02610646 
    ## Run 445 stress 0.1082722 
    ## Run 446 stress 0.1088547 
    ## Run 447 stress 0.1049903 
    ## ... Procrustes: rmse 0.006730065  max resid 0.02580724 
    ## Run 448 stress 0.1063128 
    ## Run 449 stress 0.1063128 
    ## Run 450 stress 0.1050144 
    ## ... Procrustes: rmse 0.00966302  max resid 0.02608714 
    ## Run 451 stress 0.1138392 
    ## Run 452 stress 0.1050144 
    ## ... Procrustes: rmse 0.009659451  max resid 0.02608373 
    ## Run 453 stress 0.1063128 
    ## Run 454 stress 0.1079382 
    ## Run 455 stress 0.1135598 
    ## Run 456 stress 0.1063128 
    ## Run 457 stress 0.1097989 
    ## Run 458 stress 0.1082722 
    ## Run 459 stress 0.1050144 
    ## ... Procrustes: rmse 0.009664439  max resid 0.02610145 
    ## Run 460 stress 0.1063128 
    ## Run 461 stress 0.1050144 
    ## ... Procrustes: rmse 0.009677667  max resid 0.02612631 
    ## Run 462 stress 0.1135599 
    ## Run 463 stress 0.1136646 
    ## Run 464 stress 0.1084964 
    ## Run 465 stress 0.1102722 
    ## Run 466 stress 0.1136646 
    ## Run 467 stress 0.1049903 
    ## ... Procrustes: rmse 0.006694513  max resid 0.02564551 
    ## Run 468 stress 0.1136646 
    ## Run 469 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283888  max resid 0.0772945 
    ## Run 470 stress 0.1136646 
    ## Run 471 stress 0.1063128 
    ## Run 472 stress 0.1138396 
    ## Run 473 stress 0.1088547 
    ## Run 474 stress 0.1174546 
    ## Run 475 stress 0.1084965 
    ## Run 476 stress 0.1063128 
    ## Run 477 stress 0.10835 
    ## Run 478 stress 0.1051402 
    ## ... Procrustes: rmse 0.021802  max resid 0.07686289 
    ## Run 479 stress 0.1138394 
    ## Run 480 stress 0.1084965 
    ## Run 481 stress 0.1049903 
    ## ... Procrustes: rmse 0.006722709  max resid 0.02577573 
    ## Run 482 stress 0.1160546 
    ## Run 483 stress 0.1051924 
    ## ... Procrustes: rmse 0.02283767  max resid 0.07734632 
    ## Run 484 stress 0.1063128 
    ## Run 485 stress 0.1136646 
    ## Run 486 stress 0.1082722 
    ## Run 487 stress 0.1138394 
    ## Run 488 stress 0.1174546 
    ## Run 489 stress 0.1135597 
    ## Run 490 stress 0.1063128 
    ## Run 491 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180176  max resid 0.07686616 
    ## Run 492 stress 0.1084965 
    ## Run 493 stress 0.1136647 
    ## Run 494 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 4.013099e-06  max resid 9.746237e-06 
    ## ... Similar to previous best
    ## Run 495 stress 0.1138392 
    ## Run 496 stress 0.1050338 
    ## ... Procrustes: rmse 0.006255377  max resid 0.021859 
    ## Run 497 stress 0.1050338 
    ## ... Procrustes: rmse 0.006267415  max resid 0.02191224 
    ## Run 498 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284054  max resid 0.07732273 
    ## Run 499 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180844  max resid 0.07687917 
    ## Run 500 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283954  max resid 0.07730875 
    ## *** Best solution repeated 1 times

``` r
# Surveyed sites 
SD_beta_NMDS <- metaMDS(SD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09021173 
    ## Run 1 stress 0.1071321 
    ## Run 2 stress 0.08946332 
    ## ... New best solution
    ## ... Procrustes: rmse 0.02372455  max resid 0.07469596 
    ## Run 3 stress 0.1056896 
    ## Run 4 stress 0.08946651 
    ## ... Procrustes: rmse 0.01224025  max resid 0.03878368 
    ## Run 5 stress 0.08926109 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03756121  max resid 0.1170845 
    ## Run 6 stress 0.08926066 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003147625  max resid 0.001044397 
    ## ... Similar to previous best
    ## Run 7 stress 0.09130085 
    ## Run 8 stress 0.08926068 
    ## ... Procrustes: rmse 1.424551e-05  max resid 5.125359e-05 
    ## ... Similar to previous best
    ## Run 9 stress 0.1104177 
    ## Run 10 stress 0.1091456 
    ## Run 11 stress 0.1071321 
    ## Run 12 stress 0.09021172 
    ## Run 13 stress 0.09088614 
    ## Run 14 stress 0.109614 
    ## Run 15 stress 0.0893854 
    ## ... Procrustes: rmse 0.01184195  max resid 0.03974659 
    ## Run 16 stress 0.09178286 
    ## Run 17 stress 0.08946651 
    ## ... Procrustes: rmse 0.0332414  max resid 0.1163796 
    ## Run 18 stress 0.1052649 
    ## Run 19 stress 0.0895172 
    ## ... Procrustes: rmse 0.03506311  max resid 0.1156585 
    ## Run 20 stress 0.08946328 
    ## ... Procrustes: rmse 0.03752812  max resid 0.1177916 
    ## Run 21 stress 0.1052647 
    ## Run 22 stress 0.09087346 
    ## Run 23 stress 0.1092126 
    ## Run 24 stress 0.1056895 
    ## Run 25 stress 0.08938549 
    ## ... Procrustes: rmse 0.0118224  max resid 0.03962769 
    ## Run 26 stress 0.08938541 
    ## ... Procrustes: rmse 0.01180295  max resid 0.03973215 
    ## Run 27 stress 0.1061301 
    ## Run 28 stress 0.09088614 
    ## Run 29 stress 0.1052647 
    ## Run 30 stress 0.1074435 
    ## Run 31 stress 0.08938545 
    ## ... Procrustes: rmse 0.01185347  max resid 0.03976987 
    ## Run 32 stress 0.1074431 
    ## Run 33 stress 0.1105336 
    ## Run 34 stress 0.1062566 
    ## Run 35 stress 0.08926081 
    ## ... Procrustes: rmse 0.0001391409  max resid 0.0004275671 
    ## ... Similar to previous best
    ## Run 36 stress 0.1066193 
    ## Run 37 stress 0.1052647 
    ## Run 38 stress 0.09592127 
    ## Run 39 stress 0.08926065 
    ## ... New best solution
    ## ... Procrustes: rmse 3.369398e-05  max resid 0.0001081154 
    ## ... Similar to previous best
    ## Run 40 stress 0.08946658 
    ## ... Procrustes: rmse 0.03322135  max resid 0.1163408 
    ## Run 41 stress 0.1060095 
    ## Run 42 stress 0.1062561 
    ## Run 43 stress 0.1056903 
    ## Run 44 stress 0.08926064 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001572446  max resid 0.0004198568 
    ## ... Similar to previous best
    ## Run 45 stress 0.08946329 
    ## ... Procrustes: rmse 0.0374895  max resid 0.1177482 
    ## Run 46 stress 0.1065336 
    ## Run 47 stress 0.08926072 
    ## ... Procrustes: rmse 0.0001035353  max resid 0.0002532578 
    ## ... Similar to previous best
    ## Run 48 stress 0.08926082 
    ## ... Procrustes: rmse 0.000339306  max resid 0.001010941 
    ## ... Similar to previous best
    ## Run 49 stress 0.08938966 
    ## ... Procrustes: rmse 0.03590923  max resid 0.1182281 
    ## Run 50 stress 0.1056901 
    ## Run 51 stress 0.1080637 
    ## Run 52 stress 0.1061303 
    ## Run 53 stress 0.09087351 
    ## Run 54 stress 0.09088616 
    ## Run 55 stress 0.08946328 
    ## ... Procrustes: rmse 0.03749525  max resid 0.1177586 
    ## Run 56 stress 0.09612854 
    ## Run 57 stress 0.1065376 
    ## Run 58 stress 0.08938963 
    ## ... Procrustes: rmse 0.03592022  max resid 0.1182455 
    ## Run 59 stress 0.1071322 
    ## Run 60 stress 0.1076302 
    ## Run 61 stress 0.08938542 
    ## ... Procrustes: rmse 0.01186258  max resid 0.03960645 
    ## Run 62 stress 0.08938541 
    ## ... Procrustes: rmse 0.0117741  max resid 0.03954977 
    ## Run 63 stress 0.1071321 
    ## Run 64 stress 0.0893854 
    ## ... Procrustes: rmse 0.01177749  max resid 0.03957441 
    ## Run 65 stress 0.09503435 
    ## Run 66 stress 0.1061304 
    ## Run 67 stress 0.110866 
    ## Run 68 stress 0.106419 
    ## Run 69 stress 0.1085124 
    ## Run 70 stress 0.08926075 
    ## ... Procrustes: rmse 0.0001115497  max resid 0.0003130775 
    ## ... Similar to previous best
    ## Run 71 stress 0.08938539 
    ## ... Procrustes: rmse 0.01181317  max resid 0.03962427 
    ## Run 72 stress 0.1064427 
    ## Run 73 stress 0.09088603 
    ## Run 74 stress 0.08938538 
    ## ... Procrustes: rmse 0.01180139  max resid 0.03959194 
    ## Run 75 stress 0.1061299 
    ## Run 76 stress 0.1110588 
    ## Run 77 stress 0.09099561 
    ## Run 78 stress 0.1089202 
    ## Run 79 stress 0.1091827 
    ## Run 80 stress 0.08938549 
    ## ... Procrustes: rmse 0.01188873  max resid 0.0396254 
    ## Run 81 stress 0.08926095 
    ## ... Procrustes: rmse 0.0004303957  max resid 0.001339725 
    ## ... Similar to previous best
    ## Run 82 stress 0.09178297 
    ## Run 83 stress 0.09044616 
    ## Run 84 stress 0.1074432 
    ## Run 85 stress 0.08946332 
    ## ... Procrustes: rmse 0.03751712  max resid 0.1177811 
    ## Run 86 stress 0.1052652 
    ## Run 87 stress 0.09099578 
    ## Run 88 stress 0.08938963 
    ## ... Procrustes: rmse 0.03594914  max resid 0.1182939 
    ## Run 89 stress 0.08938546 
    ## ... Procrustes: rmse 0.01175808  max resid 0.039609 
    ## Run 90 stress 0.1071321 
    ## Run 91 stress 0.1095625 
    ## Run 92 stress 0.08938971 
    ## ... Procrustes: rmse 0.03590269  max resid 0.1182167 
    ## Run 93 stress 0.08938966 
    ## ... Procrustes: rmse 0.03591075  max resid 0.1182322 
    ## Run 94 stress 0.1071322 
    ## Run 95 stress 0.08946654 
    ## ... Procrustes: rmse 0.03318261  max resid 0.1163013 
    ## Run 96 stress 0.08938963 
    ## ... Procrustes: rmse 0.03592682  max resid 0.1182541 
    ## Run 97 stress 0.1074434 
    ## Run 98 stress 0.1075422 
    ## Run 99 stress 0.1056907 
    ## Run 100 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 4.53799e-05  max resid 0.0001142424 
    ## ... Similar to previous best
    ## Run 101 stress 0.0892607 
    ## ... Procrustes: rmse 0.0001405004  max resid 0.0004674123 
    ## ... Similar to previous best
    ## Run 102 stress 0.1056898 
    ## Run 103 stress 0.1067481 
    ## Run 104 stress 0.1087579 
    ## Run 105 stress 0.1074431 
    ## Run 106 stress 0.09503444 
    ## Run 107 stress 0.08926064 
    ## ... Procrustes: rmse 9.654144e-05  max resid 0.0002236477 
    ## ... Similar to previous best
    ## Run 108 stress 0.1052649 
    ## Run 109 stress 0.08946665 
    ## ... Procrustes: rmse 0.03317379  max resid 0.1162765 
    ## Run 110 stress 0.08938559 
    ## ... Procrustes: rmse 0.01177499  max resid 0.03970233 
    ## Run 111 stress 0.1065378 
    ## Run 112 stress 0.09178294 
    ## Run 113 stress 0.09039132 
    ## Run 114 stress 0.1088389 
    ## Run 115 stress 0.1079014 
    ## Run 116 stress 0.08938557 
    ## ... Procrustes: rmse 0.01190488  max resid 0.0396001 
    ## Run 117 stress 0.09018712 
    ## Run 118 stress 0.08926092 
    ## ... Procrustes: rmse 0.0003699585  max resid 0.001196281 
    ## ... Similar to previous best
    ## Run 119 stress 0.0926237 
    ## Run 120 stress 0.09612817 
    ## Run 121 stress 0.09018711 
    ## Run 122 stress 0.09109108 
    ## Run 123 stress 0.090886 
    ## Run 124 stress 0.1096132 
    ## Run 125 stress 0.1108316 
    ## Run 126 stress 0.09503467 
    ## Run 127 stress 0.1064194 
    ## Run 128 stress 0.08946657 
    ## ... Procrustes: rmse 0.03318656  max resid 0.1162962 
    ## Run 129 stress 0.106043 
    ## Run 130 stress 0.0910911 
    ## Run 131 stress 0.1071322 
    ## Run 132 stress 0.1067902 
    ## Run 133 stress 0.08926077 
    ## ... Procrustes: rmse 0.0002758554  max resid 0.0008459118 
    ## ... Similar to previous best
    ## Run 134 stress 0.08946331 
    ## ... Procrustes: rmse 0.03750004  max resid 0.1177556 
    ## Run 135 stress 0.1056905 
    ## Run 136 stress 0.1067871 
    ## Run 137 stress 0.1065363 
    ## Run 138 stress 0.1071322 
    ## Run 139 stress 0.09109108 
    ## Run 140 stress 0.08946651 
    ## ... Procrustes: rmse 0.03320613  max resid 0.1163307 
    ## Run 141 stress 0.1088391 
    ## Run 142 stress 0.1103715 
    ## Run 143 stress 0.1065361 
    ## Run 144 stress 0.09087351 
    ## Run 145 stress 0.1056899 
    ## Run 146 stress 0.08926078 
    ## ... Procrustes: rmse 0.0002795074  max resid 0.0008563457 
    ## ... Similar to previous best
    ## Run 147 stress 0.09503413 
    ## Run 148 stress 0.08938551 
    ## ... Procrustes: rmse 0.0118299  max resid 0.03967629 
    ## Run 149 stress 0.09018707 
    ## Run 150 stress 0.1076304 
    ## Run 151 stress 0.08938541 
    ## ... Procrustes: rmse 0.01181157  max resid 0.03967735 
    ## Run 152 stress 0.1075422 
    ## Run 153 stress 0.1092209 
    ## Run 154 stress 0.08938965 
    ## ... Procrustes: rmse 0.03596622  max resid 0.118321 
    ## Run 155 stress 0.1060432 
    ## Run 156 stress 0.08938544 
    ## ... Procrustes: rmse 0.01185961  max resid 0.0396983 
    ## Run 157 stress 0.08951718 
    ## ... Procrustes: rmse 0.03503455  max resid 0.1156154 
    ## Run 158 stress 0.1099387 
    ## Run 159 stress 0.1071322 
    ## Run 160 stress 0.08938969 
    ## ... Procrustes: rmse 0.03591698  max resid 0.1182274 
    ## Run 161 stress 0.08926072 
    ## ... Procrustes: rmse 0.0001411896  max resid 0.0003793719 
    ## ... Similar to previous best
    ## Run 162 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183551  max resid 0.03966676 
    ## Run 163 stress 0.08926083 
    ## ... Procrustes: rmse 0.000316003  max resid 0.0009948822 
    ## ... Similar to previous best
    ## Run 164 stress 0.09503432 
    ## Run 165 stress 0.08946333 
    ## ... Procrustes: rmse 0.03748227  max resid 0.1177236 
    ## Run 166 stress 0.1056902 
    ## Run 167 stress 0.09021169 
    ## Run 168 stress 0.09044612 
    ## Run 169 stress 0.09503433 
    ## Run 170 stress 0.09503457 
    ## Run 171 stress 0.1092127 
    ## Run 172 stress 0.09109108 
    ## Run 173 stress 0.08946662 
    ## ... Procrustes: rmse 0.0331751  max resid 0.1162631 
    ## Run 174 stress 0.1074431 
    ## Run 175 stress 0.1085144 
    ## Run 176 stress 0.1087572 
    ## Run 177 stress 0.08946654 
    ## ... Procrustes: rmse 0.03323059  max resid 0.1163896 
    ## Run 178 stress 0.08938975 
    ## ... Procrustes: rmse 0.03590782  max resid 0.1182234 
    ## Run 179 stress 0.08938968 
    ## ... Procrustes: rmse 0.0359701  max resid 0.1183307 
    ## Run 180 stress 0.1074416 
    ## Run 181 stress 0.1056902 
    ## Run 182 stress 0.08938542 
    ## ... Procrustes: rmse 0.0117987  max resid 0.03967157 
    ## Run 183 stress 0.09503418 
    ## Run 184 stress 0.08926104 
    ## ... Procrustes: rmse 0.0004252201  max resid 0.001372046 
    ## ... Similar to previous best
    ## Run 185 stress 0.08938557 
    ## ... Procrustes: rmse 0.01189619  max resid 0.03958718 
    ## Run 186 stress 0.1071321 
    ## Run 187 stress 0.09503419 
    ## Run 188 stress 0.08926066 
    ## ... Procrustes: rmse 0.0001218844  max resid 0.0003657289 
    ## ... Similar to previous best
    ## Run 189 stress 0.1071322 
    ## Run 190 stress 0.09180904 
    ## Run 191 stress 0.09503446 
    ## Run 192 stress 0.1075421 
    ## Run 193 stress 0.1074434 
    ## Run 194 stress 0.08938975 
    ## ... Procrustes: rmse 0.03590668  max resid 0.118217 
    ## Run 195 stress 0.1091508 
    ## Run 196 stress 0.08926109 
    ## ... Procrustes: rmse 0.0004657189  max resid 0.001438146 
    ## ... Similar to previous best
    ## Run 197 stress 0.09039131 
    ## Run 198 stress 0.08926075 
    ## ... Procrustes: rmse 0.0002543459  max resid 0.0007820651 
    ## ... Similar to previous best
    ## Run 199 stress 0.1076306 
    ## Run 200 stress 0.1062557 
    ## Run 201 stress 0.1071322 
    ## Run 202 stress 0.09503439 
    ## Run 203 stress 0.09087344 
    ## Run 204 stress 0.1091877 
    ## Run 205 stress 0.0893854 
    ## ... Procrustes: rmse 0.01179694  max resid 0.03962322 
    ## Run 206 stress 0.1062594 
    ## Run 207 stress 0.1056897 
    ## Run 208 stress 0.08946329 
    ## ... Procrustes: rmse 0.03749887  max resid 0.1177608 
    ## Run 209 stress 0.09109109 
    ## Run 210 stress 0.09021171 
    ## Run 211 stress 0.09503415 
    ## Run 212 stress 0.1067503 
    ## Run 213 stress 0.1075421 
    ## Run 214 stress 0.09109114 
    ## Run 215 stress 0.1064198 
    ## Run 216 stress 0.09087358 
    ## Run 217 stress 0.1092092 
    ## Run 218 stress 0.0893897 
    ## ... Procrustes: rmse 0.03591502  max resid 0.118232 
    ## Run 219 stress 0.08951728 
    ## ... Procrustes: rmse 0.0350128  max resid 0.1155916 
    ## Run 220 stress 0.09109108 
    ## Run 221 stress 0.1052651 
    ## Run 222 stress 0.09088601 
    ## Run 223 stress 0.1063193 
    ## Run 224 stress 0.1103727 
    ## Run 225 stress 0.08926095 
    ## ... Procrustes: rmse 0.0003808486  max resid 0.001174132 
    ## ... Similar to previous best
    ## Run 226 stress 0.1103727 
    ## Run 227 stress 0.08926089 
    ## ... Procrustes: rmse 0.0003568229  max resid 0.001084342 
    ## ... Similar to previous best
    ## Run 228 stress 0.08926083 
    ## ... Procrustes: rmse 0.0003228534  max resid 0.0009991023 
    ## ... Similar to previous best
    ## Run 229 stress 0.09180878 
    ## Run 230 stress 0.08946329 
    ## ... Procrustes: rmse 0.03749725  max resid 0.1177432 
    ## Run 231 stress 0.09592157 
    ## Run 232 stress 0.08946328 
    ## ... Procrustes: rmse 0.03751223  max resid 0.1177679 
    ## Run 233 stress 0.105265 
    ## Run 234 stress 0.106419 
    ## Run 235 stress 0.09018707 
    ## Run 236 stress 0.1060429 
    ## Run 237 stress 0.09178289 
    ## Run 238 stress 0.08926098 
    ## ... Procrustes: rmse 0.000413267  max resid 0.001314518 
    ## ... Similar to previous best
    ## Run 239 stress 0.105691 
    ## Run 240 stress 0.08946654 
    ## ... Procrustes: rmse 0.03319293  max resid 0.116308 
    ## Run 241 stress 0.1067899 
    ## Run 242 stress 0.09503442 
    ## Run 243 stress 0.0893898 
    ## ... Procrustes: rmse 0.03590056  max resid 0.1182063 
    ## Run 244 stress 0.08926069 
    ## ... Procrustes: rmse 0.0001351478  max resid 0.0004229896 
    ## ... Similar to previous best
    ## Run 245 stress 0.08946331 
    ## ... Procrustes: rmse 0.03753016  max resid 0.1178054 
    ## Run 246 stress 0.08938964 
    ## ... Procrustes: rmse 0.03593376  max resid 0.1182623 
    ## Run 247 stress 0.1056896 
    ## Run 248 stress 0.1060096 
    ## Run 249 stress 0.08926063 
    ## ... Procrustes: rmse 1.197196e-05  max resid 2.405276e-05 
    ## ... Similar to previous best
    ## Run 250 stress 0.1063197 
    ## Run 251 stress 0.1067905 
    ## Run 252 stress 0.09044616 
    ## Run 253 stress 0.09044625 
    ## Run 254 stress 0.1108642 
    ## Run 255 stress 0.09021165 
    ## Run 256 stress 0.09130091 
    ## Run 257 stress 0.1091879 
    ## Run 258 stress 0.08938979 
    ## ... Procrustes: rmse 0.03590198  max resid 0.1182141 
    ## Run 259 stress 0.08938961 
    ## ... Procrustes: rmse 0.0359501  max resid 0.118287 
    ## Run 260 stress 0.0908734 
    ## Run 261 stress 0.09087347 
    ## Run 262 stress 0.09503434 
    ## Run 263 stress 0.09503431 
    ## Run 264 stress 0.08946657 
    ## ... Procrustes: rmse 0.0331891  max resid 0.1162984 
    ## Run 265 stress 0.1101443 
    ## Run 266 stress 0.1064195 
    ## Run 267 stress 0.1056893 
    ## Run 268 stress 0.08926073 
    ## ... Procrustes: rmse 0.0001372507  max resid 0.0003984349 
    ## ... Similar to previous best
    ## Run 269 stress 0.09018708 
    ## Run 270 stress 0.08938564 
    ## ... Procrustes: rmse 0.01174317  max resid 0.03952877 
    ## Run 271 stress 0.09039134 
    ## Run 272 stress 0.08946651 
    ## ... Procrustes: rmse 0.03320755  max resid 0.1163307 
    ## Run 273 stress 0.09039137 
    ## Run 274 stress 0.1076302 
    ## Run 275 stress 0.1074434 
    ## Run 276 stress 0.09088601 
    ## Run 277 stress 0.08946651 
    ## ... Procrustes: rmse 0.03322154  max resid 0.1163614 
    ## Run 278 stress 0.1071109 
    ## Run 279 stress 0.09130092 
    ## Run 280 stress 0.09180904 
    ## Run 281 stress 0.09503445 
    ## Run 282 stress 0.1074432 
    ## Run 283 stress 0.1087863 
    ## Run 284 stress 0.0894666 
    ## ... Procrustes: rmse 0.03318204  max resid 0.116291 
    ## Run 285 stress 0.08946344 
    ## ... Procrustes: rmse 0.037463  max resid 0.1176879 
    ## Run 286 stress 0.1071323 
    ## Run 287 stress 0.1075422 
    ## Run 288 stress 0.08938554 
    ## ... Procrustes: rmse 0.01188471  max resid 0.03958335 
    ## Run 289 stress 0.0908734 
    ## Run 290 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183584  max resid 0.03966437 
    ## Run 291 stress 0.08926088 
    ## ... Procrustes: rmse 0.0003525421  max resid 0.001089533 
    ## ... Similar to previous best
    ## Run 292 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595217  max resid 0.1182911 
    ## Run 293 stress 0.08946334 
    ## ... Procrustes: rmse 0.03748148  max resid 0.117723 
    ## Run 294 stress 0.1079016 
    ## Run 295 stress 0.1071323 
    ## Run 296 stress 0.08938543 
    ## ... Procrustes: rmse 0.01182855  max resid 0.03958246 
    ## Run 297 stress 0.09503454 
    ## Run 298 stress 0.08938549 
    ## ... Procrustes: rmse 0.01182115  max resid 0.0397129 
    ## Run 299 stress 0.1087565 
    ## Run 300 stress 0.105265 
    ## Run 301 stress 0.08938542 
    ## ... Procrustes: rmse 0.01180751  max resid 0.03959303 
    ## Run 302 stress 0.09109108 
    ## Run 303 stress 0.108839 
    ## Run 304 stress 0.08946656 
    ## ... Procrustes: rmse 0.03319594  max resid 0.1163165 
    ## Run 305 stress 0.1092127 
    ## Run 306 stress 0.1075422 
    ## Run 307 stress 0.08926063 
    ## ... Procrustes: rmse 2.614489e-05  max resid 7.072608e-05 
    ## ... Similar to previous best
    ## Run 308 stress 0.0893854 
    ## ... Procrustes: rmse 0.01185878  max resid 0.03963509 
    ## Run 309 stress 0.1060094 
    ## Run 310 stress 0.08926082 
    ## ... Procrustes: rmse 0.0003119779  max resid 0.0009754106 
    ## ... Similar to previous best
    ## Run 311 stress 0.09592174 
    ## Run 312 stress 0.08926087 
    ## ... Procrustes: rmse 0.0003478658  max resid 0.00110681 
    ## ... Similar to previous best
    ## Run 313 stress 0.09039153 
    ## Run 314 stress 0.09503437 
    ## Run 315 stress 0.08946663 
    ## ... Procrustes: rmse 0.03324353  max resid 0.1163925 
    ## Run 316 stress 0.08946328 
    ## ... Procrustes: rmse 0.03750309  max resid 0.1177579 
    ## Run 317 stress 0.1065366 
    ## Run 318 stress 0.1071321 
    ## Run 319 stress 0.08926104 
    ## ... Procrustes: rmse 0.0003901347  max resid 0.00116508 
    ## ... Similar to previous best
    ## Run 320 stress 0.1074433 
    ## Run 321 stress 0.1108643 
    ## Run 322 stress 0.09178298 
    ## Run 323 stress 0.08938538 
    ## ... Procrustes: rmse 0.01184278  max resid 0.03966897 
    ## Run 324 stress 0.1061301 
    ## Run 325 stress 0.08926088 
    ## ... Procrustes: rmse 0.0003054578  max resid 0.0009524554 
    ## ... Similar to previous best
    ## Run 326 stress 0.09087338 
    ## Run 327 stress 0.08946328 
    ## ... Procrustes: rmse 0.03750847  max resid 0.1177657 
    ## Run 328 stress 0.08946651 
    ## ... Procrustes: rmse 0.03322284  max resid 0.1163661 
    ## Run 329 stress 0.1067588 
    ## Run 330 stress 0.1052649 
    ## Run 331 stress 0.1085126 
    ## Run 332 stress 0.09503418 
    ## Run 333 stress 0.09088622 
    ## Run 334 stress 0.1056901 
    ## Run 335 stress 0.1065356 
    ## Run 336 stress 0.09018707 
    ## Run 337 stress 0.09178285 
    ## Run 338 stress 0.08926068 
    ## ... Procrustes: rmse 0.0001754314  max resid 0.0005287691 
    ## ... Similar to previous best
    ## Run 339 stress 0.08946333 
    ## ... Procrustes: rmse 0.03748432  max resid 0.1177282 
    ## Run 340 stress 0.1052649 
    ## Run 341 stress 0.08926075 
    ## ... Procrustes: rmse 0.0002504595  max resid 0.0007818547 
    ## ... Similar to previous best
    ## Run 342 stress 0.0893855 
    ## ... Procrustes: rmse 0.01188462  max resid 0.03960272 
    ## Run 343 stress 0.1060095 
    ## Run 344 stress 0.1076305 
    ## Run 345 stress 0.1108321 
    ## Run 346 stress 0.1074432 
    ## Run 347 stress 0.1087862 
    ## Run 348 stress 0.090886 
    ## Run 349 stress 0.1105344 
    ## Run 350 stress 0.1061301 
    ## Run 351 stress 0.1076303 
    ## Run 352 stress 0.08926096 
    ## ... Procrustes: rmse 0.0004021363  max resid 0.00128584 
    ## ... Similar to previous best
    ## Run 353 stress 0.08938541 
    ## ... Procrustes: rmse 0.01185225  max resid 0.03963203 
    ## Run 354 stress 0.1089907 
    ## Run 355 stress 0.1071321 
    ## Run 356 stress 0.09592177 
    ## Run 357 stress 0.1071067 
    ## Run 358 stress 0.08938963 
    ## ... Procrustes: rmse 0.03593666  max resid 0.1182636 
    ## Run 359 stress 0.110371 
    ## Run 360 stress 0.08926093 
    ## ... Procrustes: rmse 0.0003893577  max resid 0.001204516 
    ## ... Similar to previous best
    ## Run 361 stress 0.08926081 
    ## ... Procrustes: rmse 0.0002833934  max resid 0.0008613953 
    ## ... Similar to previous best
    ## Run 362 stress 0.09503438 
    ## Run 363 stress 0.08946331 
    ## ... Procrustes: rmse 0.03751223  max resid 0.1177907 
    ## Run 364 stress 0.08938964 
    ## ... Procrustes: rmse 0.03592943  max resid 0.1182495 
    ## Run 365 stress 0.1108464 
    ## Run 366 stress 0.08938542 
    ## ... Procrustes: rmse 0.0118049  max resid 0.03967733 
    ## Run 367 stress 0.1065364 
    ## Run 368 stress 0.08938538 
    ## ... Procrustes: rmse 0.0118307  max resid 0.0396521 
    ## Run 369 stress 0.1071322 
    ## Run 370 stress 0.08938966 
    ## ... Procrustes: rmse 0.0359244  max resid 0.1182455 
    ## Run 371 stress 0.1099393 
    ## Run 372 stress 0.1087573 
    ## Run 373 stress 0.1110405 
    ## Run 374 stress 0.09044618 
    ## Run 375 stress 0.1061302 
    ## Run 376 stress 0.1087864 
    ## Run 377 stress 0.09109108 
    ## Run 378 stress 0.08946655 
    ## ... Procrustes: rmse 0.03319079  max resid 0.1163047 
    ## Run 379 stress 0.1111379 
    ## Run 380 stress 0.1063195 
    ## Run 381 stress 0.08938568 
    ## ... Procrustes: rmse 0.01175463  max resid 0.03971583 
    ## Run 382 stress 0.1086568 
    ## Run 383 stress 0.08946651 
    ## ... Procrustes: rmse 0.03320162  max resid 0.1163212 
    ## Run 384 stress 0.1071321 
    ## Run 385 stress 0.08938551 
    ## ... Procrustes: rmse 0.01176676  max resid 0.03967609 
    ## Run 386 stress 0.08946334 
    ## ... Procrustes: rmse 0.03753054  max resid 0.1178026 
    ## Run 387 stress 0.08926067 
    ## ... Procrustes: rmse 3.019667e-05  max resid 8.32972e-05 
    ## ... Similar to previous best
    ## Run 388 stress 0.1087865 
    ## Run 389 stress 0.1062566 
    ## Run 390 stress 0.1067895 
    ## Run 391 stress 0.0892608 
    ## ... Procrustes: rmse 0.0002634078  max resid 0.0008283245 
    ## ... Similar to previous best
    ## Run 392 stress 0.1092363 
    ## Run 393 stress 0.1087567 
    ## Run 394 stress 0.08938546 
    ## ... Procrustes: rmse 0.01184808  max resid 0.03958106 
    ## Run 395 stress 0.08938538 
    ## ... Procrustes: rmse 0.01184595  max resid 0.03965882 
    ## Run 396 stress 0.1074434 
    ## Run 397 stress 0.1092093 
    ## Run 398 stress 0.08946651 
    ## ... Procrustes: rmse 0.03320388  max resid 0.1163262 
    ## Run 399 stress 0.09180905 
    ## Run 400 stress 0.08938961 
    ## ... Procrustes: rmse 0.03593929  max resid 0.1182677 
    ## Run 401 stress 0.08946339 
    ## ... Procrustes: rmse 0.03747055  max resid 0.1177081 
    ## Run 402 stress 0.1071321 
    ## Run 403 stress 0.1089907 
    ## Run 404 stress 0.09503418 
    ## Run 405 stress 0.08938964 
    ## ... Procrustes: rmse 0.03592781  max resid 0.1182519 
    ## Run 406 stress 0.1063192 
    ## Run 407 stress 0.08938966 
    ## ... Procrustes: rmse 0.03596681  max resid 0.1183225 
    ## Run 408 stress 0.08938964 
    ## ... Procrustes: rmse 0.03592841  max resid 0.1182533 
    ## Run 409 stress 0.08946652 
    ## ... Procrustes: rmse 0.03320052  max resid 0.1163199 
    ## Run 410 stress 0.08938964 
    ## ... Procrustes: rmse 0.03596412  max resid 0.1183155 
    ## Run 411 stress 0.1061307 
    ## Run 412 stress 0.1056902 
    ## Run 413 stress 0.0895172 
    ## ... Procrustes: rmse 0.03502965  max resid 0.1156085 
    ## Run 414 stress 0.1065381 
    ## Run 415 stress 0.09021183 
    ## Run 416 stress 0.08946337 
    ## ... Procrustes: rmse 0.03747396  max resid 0.1177115 
    ## Run 417 stress 0.1071323 
    ## Run 418 stress 0.1076305 
    ## Run 419 stress 0.1080638 
    ## Run 420 stress 0.1093479 
    ## Run 421 stress 0.08946329 
    ## ... Procrustes: rmse 0.03750008  max resid 0.1177321 
    ## Run 422 stress 0.0913009 
    ## Run 423 stress 0.09044613 
    ## Run 424 stress 0.08946331 
    ## ... Procrustes: rmse 0.03748884  max resid 0.1177349 
    ## Run 425 stress 0.1056897 
    ## Run 426 stress 0.106443 
    ## Run 427 stress 0.08938975 
    ## ... Procrustes: rmse 0.03591098  max resid 0.118224 
    ## Run 428 stress 0.1074434 
    ## Run 429 stress 0.09503415 
    ## Run 430 stress 0.09088605 
    ## Run 431 stress 0.08926088 
    ## ... Procrustes: rmse 0.0003416354  max resid 0.00106935 
    ## ... Similar to previous best
    ## Run 432 stress 0.09503418 
    ## Run 433 stress 0.08938962 
    ## ... Procrustes: rmse 0.03595921  max resid 0.1183026 
    ## Run 434 stress 0.08926077 
    ## ... Procrustes: rmse 0.0002692555  max resid 0.0008400816 
    ## ... Similar to previous best
    ## Run 435 stress 0.1086562 
    ## Run 436 stress 0.1092126 
    ## Run 437 stress 0.08938965 
    ## ... Procrustes: rmse 0.03596588  max resid 0.1183189 
    ## Run 438 stress 0.1076302 
    ## Run 439 stress 0.08938962 
    ## ... Procrustes: rmse 0.03593552  max resid 0.1182645 
    ## Run 440 stress 0.08926109 
    ## ... Procrustes: rmse 0.0004676632  max resid 0.001523334 
    ## ... Similar to previous best
    ## Run 441 stress 0.1065373 
    ## Run 442 stress 0.1067599 
    ## Run 443 stress 0.1067512 
    ## Run 444 stress 0.08951717 
    ## ... Procrustes: rmse 0.0350551  max resid 0.1156387 
    ## Run 445 stress 0.1093316 
    ## Run 446 stress 0.08938963 
    ## ... Procrustes: rmse 0.03595768  max resid 0.1183021 
    ## Run 447 stress 0.08926066 
    ## ... Procrustes: rmse 9.187132e-05  max resid 0.0002883283 
    ## ... Similar to previous best
    ## Run 448 stress 0.1071322 
    ## Run 449 stress 0.1067586 
    ## Run 450 stress 0.09262364 
    ## Run 451 stress 0.1056903 
    ## Run 452 stress 0.1092126 
    ## Run 453 stress 0.1087865 
    ## Run 454 stress 0.1103719 
    ## Run 455 stress 0.08946651 
    ## ... Procrustes: rmse 0.03322027  max resid 0.1163588 
    ## Run 456 stress 0.09592171 
    ## Run 457 stress 0.1064197 
    ## Run 458 stress 0.1076302 
    ## Run 459 stress 0.1064193 
    ## Run 460 stress 0.09228491 
    ## Run 461 stress 0.09503443 
    ## Run 462 stress 0.0894633 
    ## ... Procrustes: rmse 0.03749321  max resid 0.1177424 
    ## Run 463 stress 0.08938543 
    ## ... Procrustes: rmse 0.0117963  max resid 0.03967816 
    ## Run 464 stress 0.08938542 
    ## ... Procrustes: rmse 0.01184989  max resid 0.03960323 
    ## Run 465 stress 0.09503415 
    ## Run 466 stress 0.08926063 
    ## ... Procrustes: rmse 6.811197e-05  max resid 0.0001786333 
    ## ... Similar to previous best
    ## Run 467 stress 0.1063191 
    ## Run 468 stress 0.08946328 
    ## ... Procrustes: rmse 0.03751185  max resid 0.1177736 
    ## Run 469 stress 0.08946329 
    ## ... Procrustes: rmse 0.0375002  max resid 0.1177503 
    ## Run 470 stress 0.1071322 
    ## Run 471 stress 0.08946651 
    ## ... Procrustes: rmse 0.03321802  max resid 0.1163488 
    ## Run 472 stress 0.08926075 
    ## ... Procrustes: rmse 0.0002478676  max resid 0.0007565734 
    ## ... Similar to previous best
    ## Run 473 stress 0.1056906 
    ## Run 474 stress 0.1060432 
    ## Run 475 stress 0.1060094 
    ## Run 476 stress 0.08926066 
    ## ... Procrustes: rmse 0.000102744  max resid 0.0003457622 
    ## ... Similar to previous best
    ## Run 477 stress 0.09021166 
    ## Run 478 stress 0.09099572 
    ## Run 479 stress 0.1056902 
    ## Run 480 stress 0.1065369 
    ## Run 481 stress 0.09503424 
    ## Run 482 stress 0.09262344 
    ## Run 483 stress 0.08938966 
    ## ... Procrustes: rmse 0.03596755  max resid 0.1183236 
    ## Run 484 stress 0.09109108 
    ## Run 485 stress 0.09503425 
    ## Run 486 stress 0.08938542 
    ## ... Procrustes: rmse 0.01185311  max resid 0.03969467 
    ## Run 487 stress 0.09088613 
    ## Run 488 stress 0.08926063 
    ## ... Procrustes: rmse 6.327128e-05  max resid 0.0001861823 
    ## ... Similar to previous best
    ## Run 489 stress 0.1052651 
    ## Run 490 stress 0.08926112 
    ## ... Procrustes: rmse 0.0004838758  max resid 0.001506809 
    ## ... Similar to previous best
    ## Run 491 stress 0.08926067 
    ## ... Procrustes: rmse 0.0001544813  max resid 0.0004660102 
    ## ... Similar to previous best
    ## Run 492 stress 0.1056902 
    ## Run 493 stress 0.1064427 
    ## Run 494 stress 0.09612868 
    ## Run 495 stress 0.08938966 
    ## ... Procrustes: rmse 0.03592346  max resid 0.1182464 
    ## Run 496 stress 0.1074432 
    ## Run 497 stress 0.08926071 
    ## ... Procrustes: rmse 0.0001608664  max resid 0.0005038377 
    ## ... Similar to previous best
    ## Run 498 stress 0.08938976 
    ## ... Procrustes: rmse 0.0359161  max resid 0.1182293 
    ## Run 499 stress 0.09018707 
    ## Run 500 stress 0.09018707 
    ## *** Best solution repeated 43 times

``` r
round(SD_beta_NMDS$stress, digits = 2)
```

    ## [1] 0.09

``` r
SD_beta_rep_NMDS <- metaMDS(SD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.3323292 
    ## Run 1 stress 0.3477205 
    ## Run 2 stress 0.3393891 
    ## Run 3 stress 0.3310174 
    ## ... New best solution
    ## ... Procrustes: rmse 0.2051922  max resid 0.3125521 
    ## Run 4 stress 0.3397407 
    ## Run 5 stress 0.3328103 
    ## Run 6 stress 0.3348521 
    ## Run 7 stress 0.3329455 
    ## Run 8 stress 0.3274821 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1887699  max resid 0.2903252 
    ## Run 9 stress 0.3279689 
    ## ... Procrustes: rmse 0.125435  max resid 0.2992949 
    ## Run 10 stress 0.344441 
    ## Run 11 stress 0.3311506 
    ## Run 12 stress 0.3377399 
    ## Run 13 stress 0.3380822 
    ## Run 14 stress 0.3257292 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1631847  max resid 0.3013942 
    ## Run 15 stress 0.3302954 
    ## Run 16 stress 0.345307 
    ## Run 17 stress 0.3342958 
    ## Run 18 stress 0.345834 
    ## Run 19 stress 0.3362337 
    ## Run 20 stress 0.3296661 
    ## Run 21 stress 0.3448237 
    ## Run 22 stress 0.3338715 
    ## Run 23 stress 0.3367941 
    ## Run 24 stress 0.3311099 
    ## Run 25 stress 0.3347155 
    ## Run 26 stress 0.3281724 
    ## Run 27 stress 0.3464853 
    ## Run 28 stress 0.3304485 
    ## Run 29 stress 0.3304801 
    ## Run 30 stress 0.3448738 
    ## Run 31 stress 0.3362703 
    ## Run 32 stress 0.3326862 
    ## Run 33 stress 0.3422811 
    ## Run 34 stress 0.3302607 
    ## Run 35 stress 0.3378078 
    ## Run 36 stress 0.3504092 
    ## Run 37 stress 0.3329337 
    ## Run 38 stress 0.3345621 
    ## Run 39 stress 0.3321225 
    ## Run 40 stress 0.3361582 
    ## Run 41 stress 0.3295194 
    ## Run 42 stress 0.3400534 
    ## Run 43 stress 0.3377178 
    ## Run 44 stress 0.341957 
    ## Run 45 stress 0.3400975 
    ## Run 46 stress 0.3399593 
    ## Run 47 stress 0.343949 
    ## Run 48 stress 0.3395899 
    ## Run 49 stress 0.3310406 
    ## Run 50 stress 0.3272015 
    ## Run 51 stress 0.3295699 
    ## Run 52 stress 0.342644 
    ## Run 53 stress 0.3331379 
    ## Run 54 stress 0.343699 
    ## Run 55 stress 0.3342502 
    ## Run 56 stress 0.3277876 
    ## Run 57 stress 0.3327845 
    ## Run 58 stress 0.3361677 
    ## Run 59 stress 0.3302233 
    ## Run 60 stress 0.3378371 
    ## Run 61 stress 0.3338811 
    ## Run 62 stress 0.3427105 
    ## Run 63 stress 0.3473636 
    ## Run 64 stress 0.3362589 
    ## Run 65 stress 0.3479628 
    ## Run 66 stress 0.339936 
    ## Run 67 stress 0.3367427 
    ## Run 68 stress 0.3336442 
    ## Run 69 stress 0.3285008 
    ## Run 70 stress 0.3368508 
    ## Run 71 stress 0.3382128 
    ## Run 72 stress 0.3403189 
    ## Run 73 stress 0.3419922 
    ## Run 74 stress 0.3311181 
    ## Run 75 stress 0.3322371 
    ## Run 76 stress 0.3455809 
    ## Run 77 stress 0.3343225 
    ## Run 78 stress 0.3355978 
    ## Run 79 stress 0.325401 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1043995  max resid 0.2185487 
    ## Run 80 stress 0.3459697 
    ## Run 81 stress 0.3362896 
    ## Run 82 stress 0.3317101 
    ## Run 83 stress 0.3310327 
    ## Run 84 stress 0.3444612 
    ## Run 85 stress 0.3340817 
    ## Run 86 stress 0.3445708 
    ## Run 87 stress 0.338198 
    ## Run 88 stress 0.3467837 
    ## Run 89 stress 0.3563306 
    ## Run 90 stress 0.342636 
    ## Run 91 stress 0.3395253 
    ## Run 92 stress 0.329879 
    ## Run 93 stress 0.3439435 
    ## Run 94 stress 0.3441523 
    ## Run 95 stress 0.3461067 
    ## Run 96 stress 0.3361784 
    ## Run 97 stress 0.3330091 
    ## Run 98 stress 0.3408097 
    ## Run 99 stress 0.3417494 
    ## Run 100 stress 0.3351007 
    ## Run 101 stress 0.3377691 
    ## Run 102 stress 0.3448168 
    ## Run 103 stress 0.3368338 
    ## Run 104 stress 0.3557296 
    ## Run 105 stress 0.3536736 
    ## Run 106 stress 0.3274071 
    ## Run 107 stress 0.3279411 
    ## Run 108 stress 0.3477768 
    ## Run 109 stress 0.3265958 
    ## Run 110 stress 0.3342363 
    ## Run 111 stress 0.3402139 
    ## Run 112 stress 0.3484128 
    ## Run 113 stress 0.3483411 
    ## Run 114 stress 0.3376883 
    ## Run 115 stress 0.3378223 
    ## Run 116 stress 0.3459105 
    ## Run 117 stress 0.3319752 
    ## Run 118 stress 0.3370624 
    ## Run 119 stress 0.3296484 
    ## Run 120 stress 0.3320944 
    ## Run 121 stress 0.3296876 
    ## Run 122 stress 0.3299268 
    ## Run 123 stress 0.3349951 
    ## Run 124 stress 0.3318048 
    ## Run 125 stress 0.3352392 
    ## Run 126 stress 0.3337996 
    ## Run 127 stress 0.3341507 
    ## Run 128 stress 0.3370334 
    ## Run 129 stress 0.3345234 
    ## Run 130 stress 0.3376151 
    ## Run 131 stress 0.3469585 
    ## Run 132 stress 0.343131 
    ## Run 133 stress 0.3365962 
    ## Run 134 stress 0.3346546 
    ## Run 135 stress 0.3325195 
    ## Run 136 stress 0.3341644 
    ## Run 137 stress 0.3379499 
    ## Run 138 stress 0.3305007 
    ## Run 139 stress 0.3454731 
    ## Run 140 stress 0.3378797 
    ## Run 141 stress 0.3413805 
    ## Run 142 stress 0.3330843 
    ## Run 143 stress 0.3359614 
    ## Run 144 stress 0.3451361 
    ## Run 145 stress 0.3374943 
    ## Run 146 stress 0.3272611 
    ## Run 147 stress 0.3419364 
    ## Run 148 stress 0.3352123 
    ## Run 149 stress 0.3440531 
    ## Run 150 stress 0.3331528 
    ## Run 151 stress 0.3309173 
    ## Run 152 stress 0.3443137 
    ## Run 153 stress 0.3456009 
    ## Run 154 stress 0.3461754 
    ## Run 155 stress 0.3335619 
    ## Run 156 stress 0.3374247 
    ## Run 157 stress 0.3295081 
    ## Run 158 stress 0.3348225 
    ## Run 159 stress 0.3295527 
    ## Run 160 stress 0.3436688 
    ## Run 161 stress 0.345884 
    ## Run 162 stress 0.3390804 
    ## Run 163 stress 0.3405281 
    ## Run 164 stress 0.3354104 
    ## Run 165 stress 0.336927 
    ## Run 166 stress 0.3338617 
    ## Run 167 stress 0.3281862 
    ## Run 168 stress 0.3361767 
    ## Run 169 stress 0.3451394 
    ## Run 170 stress 0.3363158 
    ## Run 171 stress 0.3486268 
    ## Run 172 stress 0.3392873 
    ## Run 173 stress 0.3352619 
    ## Run 174 stress 0.3367235 
    ## Run 175 stress 0.3321536 
    ## Run 176 stress 0.3351152 
    ## Run 177 stress 0.3311728 
    ## Run 178 stress 0.3395447 
    ## Run 179 stress 0.3383175 
    ## Run 180 stress 0.3303096 
    ## Run 181 stress 0.3455987 
    ## Run 182 stress 0.3284972 
    ## Run 183 stress 0.3369741 
    ## Run 184 stress 0.3344346 
    ## Run 185 stress 0.3394776 
    ## Run 186 stress 0.3390534 
    ## Run 187 stress 0.344019 
    ## Run 188 stress 0.3347983 
    ## Run 189 stress 0.3456333 
    ## Run 190 stress 0.3395976 
    ## Run 191 stress 0.3357594 
    ## Run 192 stress 0.3391454 
    ## Run 193 stress 0.3343868 
    ## Run 194 stress 0.343363 
    ## Run 195 stress 0.3494829 
    ## Run 196 stress 0.3406012 
    ## Run 197 stress 0.3457875 
    ## Run 198 stress 0.3361102 
    ## Run 199 stress 0.3417356 
    ## Run 200 stress 0.3334557 
    ## Run 201 stress 0.3342421 
    ## Run 202 stress 0.3499989 
    ## Run 203 stress 0.34232 
    ## Run 204 stress 0.343715 
    ## Run 205 stress 0.3292742 
    ## Run 206 stress 0.3333191 
    ## Run 207 stress 0.329451 
    ## Run 208 stress 0.333867 
    ## Run 209 stress 0.3363882 
    ## Run 210 stress 0.3315401 
    ## Run 211 stress 0.3497717 
    ## Run 212 stress 0.329641 
    ## Run 213 stress 0.3328555 
    ## Run 214 stress 0.3380201 
    ## Run 215 stress 0.3301022 
    ## Run 216 stress 0.3341346 
    ## Run 217 stress 0.3322433 
    ## Run 218 stress 0.3304443 
    ## Run 219 stress 0.3489742 
    ## Run 220 stress 0.3329645 
    ## Run 221 stress 0.3366984 
    ## Run 222 stress 0.3334305 
    ## Run 223 stress 0.3441065 
    ## Run 224 stress 0.3379053 
    ## Run 225 stress 0.3450252 
    ## Run 226 stress 0.3342912 
    ## Run 227 stress 0.3389479 
    ## Run 228 stress 0.330955 
    ## Run 229 stress 0.3283885 
    ## Run 230 stress 0.3373691 
    ## Run 231 stress 0.3448828 
    ## Run 232 stress 0.3346012 
    ## Run 233 stress 0.3430089 
    ## Run 234 stress 0.3380529 
    ## Run 235 stress 0.3318835 
    ## Run 236 stress 0.338848 
    ## Run 237 stress 0.3337156 
    ## Run 238 stress 0.3505367 
    ## Run 239 stress 0.3375125 
    ## Run 240 stress 0.3438984 
    ## Run 241 stress 0.3394731 
    ## Run 242 stress 0.3342346 
    ## Run 243 stress 0.327141 
    ## Run 244 stress 0.3363486 
    ## Run 245 stress 0.3357564 
    ## Run 246 stress 0.3353031 
    ## Run 247 stress 0.345168 
    ## Run 248 stress 0.3365049 
    ## Run 249 stress 0.3266709 
    ## Run 250 stress 0.3411907 
    ## Run 251 stress 0.3520086 
    ## Run 252 stress 0.3446231 
    ## Run 253 stress 0.3427463 
    ## Run 254 stress 0.3280338 
    ## Run 255 stress 0.3348124 
    ## Run 256 stress 0.3421409 
    ## Run 257 stress 0.3334386 
    ## Run 258 stress 0.3340035 
    ## Run 259 stress 0.3266756 
    ## Run 260 stress 0.3438989 
    ## Run 261 stress 0.3296173 
    ## Run 262 stress 0.3485204 
    ## Run 263 stress 0.3337494 
    ## Run 264 stress 0.3465632 
    ## Run 265 stress 0.3332962 
    ## Run 266 stress 0.335043 
    ## Run 267 stress 0.3318217 
    ## Run 268 stress 0.3359686 
    ## Run 269 stress 0.3388657 
    ## Run 270 stress 0.340153 
    ## Run 271 stress 0.3287211 
    ## Run 272 stress 0.3356784 
    ## Run 273 stress 0.3447692 
    ## Run 274 stress 0.3371359 
    ## Run 275 stress 0.3273408 
    ## Run 276 stress 0.3333296 
    ## Run 277 stress 0.3324332 
    ## Run 278 stress 0.3341343 
    ## Run 279 stress 0.3428412 
    ## Run 280 stress 0.3338432 
    ## Run 281 stress 0.342408 
    ## Run 282 stress 0.3349024 
    ## Run 283 stress 0.3456266 
    ## Run 284 stress 0.3255777 
    ## ... Procrustes: rmse 0.1875827  max resid 0.3473997 
    ## Run 285 stress 0.3269593 
    ## Run 286 stress 0.3385033 
    ## Run 287 stress 0.3311583 
    ## Run 288 stress 0.3304721 
    ## Run 289 stress 0.3339118 
    ## Run 290 stress 0.3306412 
    ## Run 291 stress 0.334619 
    ## Run 292 stress 0.333688 
    ## Run 293 stress 0.3438459 
    ## Run 294 stress 0.333679 
    ## Run 295 stress 0.3445737 
    ## Run 296 stress 0.3441017 
    ## Run 297 stress 0.3438261 
    ## Run 298 stress 0.330924 
    ## Run 299 stress 0.3414375 
    ## Run 300 stress 0.3326393 
    ## Run 301 stress 0.3303866 
    ## Run 302 stress 0.3335337 
    ## Run 303 stress 0.3346094 
    ## Run 304 stress 0.3405141 
    ## Run 305 stress 0.3261708 
    ## Run 306 stress 0.3500453 
    ## Run 307 stress 0.3373598 
    ## Run 308 stress 0.3359387 
    ## Run 309 stress 0.3378156 
    ## Run 310 stress 0.3365442 
    ## Run 311 stress 0.3348251 
    ## Run 312 stress 0.337027 
    ## Run 313 stress 0.3289294 
    ## Run 314 stress 0.3464468 
    ## Run 315 stress 0.335308 
    ## Run 316 stress 0.3343991 
    ## Run 317 stress 0.337967 
    ## Run 318 stress 0.3401782 
    ## Run 319 stress 0.3306574 
    ## Run 320 stress 0.3266842 
    ## Run 321 stress 0.3424873 
    ## Run 322 stress 0.3424747 
    ## Run 323 stress 0.3312634 
    ## Run 324 stress 0.3374495 
    ## Run 325 stress 0.3378243 
    ## Run 326 stress 0.344489 
    ## Run 327 stress 0.3370302 
    ## Run 328 stress 0.3286507 
    ## Run 329 stress 0.3414318 
    ## Run 330 stress 0.3350659 
    ## Run 331 stress 0.3349956 
    ## Run 332 stress 0.3337889 
    ## Run 333 stress 0.3307617 
    ## Run 334 stress 0.3277174 
    ## Run 335 stress 0.3488369 
    ## Run 336 stress 0.3358406 
    ## Run 337 stress 0.3310145 
    ## Run 338 stress 0.3419236 
    ## Run 339 stress 0.3300401 
    ## Run 340 stress 0.326369 
    ## Run 341 stress 0.3270633 
    ## Run 342 stress 0.3335387 
    ## Run 343 stress 0.3478992 
    ## Run 344 stress 0.3426288 
    ## Run 345 stress 0.3358551 
    ## Run 346 stress 0.3377315 
    ## Run 347 stress 0.3408469 
    ## Run 348 stress 0.3422312 
    ## Run 349 stress 0.3382376 
    ## Run 350 stress 0.3354356 
    ## Run 351 stress 0.3307856 
    ## Run 352 stress 0.3337567 
    ## Run 353 stress 0.3322548 
    ## Run 354 stress 0.3301324 
    ## Run 355 stress 0.3401181 
    ## Run 356 stress 0.3460957 
    ## Run 357 stress 0.3407849 
    ## Run 358 stress 0.3375581 
    ## Run 359 stress 0.3265269 
    ## Run 360 stress 0.3429794 
    ## Run 361 stress 0.3321661 
    ## Run 362 stress 0.3434485 
    ## Run 363 stress 0.3471126 
    ## Run 364 stress 0.3316828 
    ## Run 365 stress 0.3373528 
    ## Run 366 stress 0.3274787 
    ## Run 367 stress 0.3343282 
    ## Run 368 stress 0.3349542 
    ## Run 369 stress 0.3394684 
    ## Run 370 stress 0.3293595 
    ## Run 371 stress 0.3507815 
    ## Run 372 stress 0.3437978 
    ## Run 373 stress 0.334377 
    ## Run 374 stress 0.3299606 
    ## Run 375 stress 0.3329241 
    ## Run 376 stress 0.3417735 
    ## Run 377 stress 0.3363746 
    ## Run 378 stress 0.3338037 
    ## Run 379 stress 0.3300163 
    ## Run 380 stress 0.3392632 
    ## Run 381 stress 0.3307801 
    ## Run 382 stress 0.3327167 
    ## Run 383 stress 0.3388389 
    ## Run 384 stress 0.3419264 
    ## Run 385 stress 0.3412711 
    ## Run 386 stress 0.3511142 
    ## Run 387 stress 0.345985 
    ## Run 388 stress 0.3353473 
    ## Run 389 stress 0.337323 
    ## Run 390 stress 0.3450473 
    ## Run 391 stress 0.3275857 
    ## Run 392 stress 0.3405165 
    ## Run 393 stress 0.3426935 
    ## Run 394 stress 0.3310973 
    ## Run 395 stress 0.3393022 
    ## Run 396 stress 0.330863 
    ## Run 397 stress 0.3464534 
    ## Run 398 stress 0.3318489 
    ## Run 399 stress 0.3415152 
    ## Run 400 stress 0.347048 
    ## Run 401 stress 0.3357362 
    ## Run 402 stress 0.3344668 
    ## Run 403 stress 0.3376654 
    ## Run 404 stress 0.3361544 
    ## Run 405 stress 0.3351647 
    ## Run 406 stress 0.3410473 
    ## Run 407 stress 0.3277172 
    ## Run 408 stress 0.3343576 
    ## Run 409 stress 0.3466937 
    ## Run 410 stress 0.3329528 
    ## Run 411 stress 0.3286907 
    ## Run 412 stress 0.3348991 
    ## Run 413 stress 0.3318497 
    ## Run 414 stress 0.3437579 
    ## Run 415 stress 0.3435923 
    ## Run 416 stress 0.3443866 
    ## Run 417 stress 0.347116 
    ## Run 418 stress 0.3360803 
    ## Run 419 stress 0.3360713 
    ## Run 420 stress 0.3439601 
    ## Run 421 stress 0.3326257 
    ## Run 422 stress 0.3380095 
    ## Run 423 stress 0.3320352 
    ## Run 424 stress 0.3297284 
    ## Run 425 stress 0.3374106 
    ## Run 426 stress 0.3350132 
    ## Run 427 stress 0.3436127 
    ## Run 428 stress 0.3408108 
    ## Run 429 stress 0.3366362 
    ## Run 430 stress 0.3302753 
    ## Run 431 stress 0.3462494 
    ## Run 432 stress 0.3314348 
    ## Run 433 stress 0.336229 
    ## Run 434 stress 0.3431712 
    ## Run 435 stress 0.3446121 
    ## Run 436 stress 0.3399871 
    ## Run 437 stress 0.3317857 
    ## Run 438 stress 0.3369619 
    ## Run 439 stress 0.3412335 
    ## Run 440 stress 0.3334051 
    ## Run 441 stress 0.3590658 
    ## Run 442 stress 0.3258176 
    ## ... Procrustes: rmse 0.184165  max resid 0.3703082 
    ## Run 443 stress 0.3300407 
    ## Run 444 stress 0.3359614 
    ## Run 445 stress 0.3272196 
    ## Run 446 stress 0.3350114 
    ## Run 447 stress 0.3339274 
    ## Run 448 stress 0.3317298 
    ## Run 449 stress 0.332162 
    ## Run 450 stress 0.3532004 
    ## Run 451 stress 0.3457028 
    ## Run 452 stress 0.3354388 
    ## Run 453 stress 0.3342147 
    ## Run 454 stress 0.3291457 
    ## Run 455 stress 0.3367427 
    ## Run 456 stress 0.332602 
    ## Run 457 stress 0.3289894 
    ## Run 458 stress 0.3307239 
    ## Run 459 stress 0.3346727 
    ## Run 460 stress 0.3488789 
    ## Run 461 stress 0.3276049 
    ## Run 462 stress 0.3359852 
    ## Run 463 stress 0.3394742 
    ## Run 464 stress 0.3398216 
    ## Run 465 stress 0.3290607 
    ## Run 466 stress 0.3329517 
    ## Run 467 stress 0.343687 
    ## Run 468 stress 0.3474706 
    ## Run 469 stress 0.3331451 
    ## Run 470 stress 0.3378798 
    ## Run 471 stress 0.3604261 
    ## Run 472 stress 0.333802 
    ## Run 473 stress 0.330971 
    ## Run 474 stress 0.3344493 
    ## Run 475 stress 0.3238982 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1785097  max resid 0.3167929 
    ## Run 476 stress 0.3289958 
    ## Run 477 stress 0.3330246 
    ## Run 478 stress 0.3457383 
    ## Run 479 stress 0.3378241 
    ## Run 480 stress 0.3367103 
    ## Run 481 stress 0.3372988 
    ## Run 482 stress 0.3350357 
    ## Run 483 stress 0.3377461 
    ## Run 484 stress 0.3347461 
    ## Run 485 stress 0.3419975 
    ## Run 486 stress 0.3452214 
    ## Run 487 stress 0.3394911 
    ## Run 488 stress 0.3339448 
    ## Run 489 stress 0.3365251 
    ## Run 490 stress 0.3407958 
    ## Run 491 stress 0.344803 
    ## Run 492 stress 0.3385303 
    ## Run 493 stress 0.3347435 
    ## Run 494 stress 0.3418132 
    ## Run 495 stress 0.3450853 
    ## Run 496 stress 0.3440126 
    ## Run 497 stress 0.330025 
    ## Run 498 stress 0.3304545 
    ## Run 499 stress 0.3332359 
    ## Run 500 stress 0.3441669 
    ## *** Best solution was not repeated -- monoMDS stopping criteria:
    ##    500: stress ratio > sratmax

``` r
SD_beta_ric_NMDS <- metaMDS(SD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.01879552 
    ## Run 1 stress 0.02492209 
    ## Run 2 stress 0.0187957 
    ## ... Procrustes: rmse 8.810404e-05  max resid 0.000180803 
    ## ... Similar to previous best
    ## Run 3 stress 0.02509457 
    ## Run 4 stress 0.01879574 
    ## ... Procrustes: rmse 0.001235931  max resid 0.002545466 
    ## ... Similar to previous best
    ## Run 5 stress 0.01879552 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001132813  max resid 0.002333155 
    ## ... Similar to previous best
    ## Run 6 stress 0.01879581 
    ## ... Procrustes: rmse 0.0001351062  max resid 0.0002783315 
    ## ... Similar to previous best
    ## Run 7 stress 0.0187956 
    ## ... Procrustes: rmse 4.170714e-05  max resid 8.588439e-05 
    ## ... Similar to previous best
    ## Run 8 stress 0.01886178 
    ## ... Procrustes: rmse 0.004807406  max resid 0.009559606 
    ## ... Similar to previous best
    ## Run 9 stress 0.01905542 
    ## ... Procrustes: rmse 0.01189604  max resid 0.02427352 
    ## Run 10 stress 0.01887066 
    ## ... Procrustes: rmse 0.005344815  max resid 0.0106888 
    ## Run 11 stress 0.02492203 
    ## Run 12 stress 0.01879658 
    ## ... Procrustes: rmse 0.0004122275  max resid 0.0008493134 
    ## ... Similar to previous best
    ## Run 13 stress 0.01879538 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001050669  max resid 0.002166579 
    ## ... Similar to previous best
    ## Run 14 stress 0.02486132 
    ## Run 15 stress 0.02509421 
    ## Run 16 stress 0.02492165 
    ## Run 17 stress 0.01885858 
    ## ... Procrustes: rmse 0.005490172  max resid 0.01100597 
    ## Run 18 stress 0.02492171 
    ## Run 19 stress 0.02520899 
    ## Run 20 stress 0.0249233 
    ## Run 21 stress 0.02509453 
    ## Run 22 stress 0.01879573 
    ## ... Procrustes: rmse 0.001139182  max resid 0.002346821 
    ## ... Similar to previous best
    ## Run 23 stress 0.02509451 
    ## Run 24 stress 0.01879591 
    ## ... Procrustes: rmse 0.001217231  max resid 0.002507399 
    ## ... Similar to previous best
    ## Run 25 stress 0.02509455 
    ## Run 26 stress 0.01879969 
    ## ... Procrustes: rmse 0.002163914  max resid 0.004456473 
    ## ... Similar to previous best
    ## Run 27 stress 0.02520905 
    ## Run 28 stress 0.01945635 
    ## Run 29 stress 0.02509463 
    ## Run 30 stress 0.01903551 
    ## ... Procrustes: rmse 0.01271166  max resid 0.02604341 
    ## Run 31 stress 0.02509445 
    ## Run 32 stress 0.02492166 
    ## Run 33 stress 0.0250943 
    ## Run 34 stress 0.01879567 
    ## ... Procrustes: rmse 0.001122414  max resid 0.002312181 
    ## ... Similar to previous best
    ## Run 35 stress 0.0188348 
    ## ... Procrustes: rmse 0.005406611  max resid 0.01112074 
    ## Run 36 stress 0.01889755 
    ## ... Procrustes: rmse 0.007734734  max resid 0.01567371 
    ## Run 37 stress 0.3710573 
    ## Run 38 stress 0.02509454 
    ## Run 39 stress 0.0187959 
    ## ... Procrustes: rmse 0.001218087  max resid 0.002509174 
    ## ... Similar to previous best
    ## Run 40 stress 0.02520876 
    ## Run 41 stress 0.01885518 
    ## ... Procrustes: rmse 0.006501779  max resid 0.01336102 
    ## Run 42 stress 0.01879596 
    ## ... Procrustes: rmse 0.001236709  max resid 0.002547518 
    ## ... Similar to previous best
    ## Run 43 stress 0.01892383 
    ## ... Procrustes: rmse 0.008879615  max resid 0.01805021 
    ## Run 44 stress 0.02520888 
    ## Run 45 stress 0.01879606 
    ## ... Procrustes: rmse 0.001275333  max resid 0.002626617 
    ## ... Similar to previous best
    ## Run 46 stress 0.02520888 
    ## Run 47 stress 0.01883157 
    ## ... Procrustes: rmse 0.003119606  max resid 0.006003497 
    ## ... Similar to previous best
    ## Run 48 stress 0.0188221 
    ## ... Procrustes: rmse 0.002506864  max resid 0.004925652 
    ## ... Similar to previous best
    ## Run 49 stress 0.0249219 
    ## Run 50 stress 0.02509469 
    ## Run 51 stress 0.01879613 
    ## ... Procrustes: rmse 0.001312036  max resid 0.002702339 
    ## ... Similar to previous best
    ## Run 52 stress 0.01918577 
    ## ... Procrustes: rmse 0.01531891  max resid 0.03134752 
    ## Run 53 stress 0.01879581 
    ## ... Procrustes: rmse 0.001188991  max resid 0.002449135 
    ## ... Similar to previous best
    ## Run 54 stress 0.02509451 
    ## Run 55 stress 0.01908948 
    ## ... Procrustes: rmse 0.01360652  max resid 0.02776617 
    ## Run 56 stress 0.02520902 
    ## Run 57 stress 0.01883618 
    ## ... Procrustes: rmse 0.004999351  max resid 0.01027152 
    ## Run 58 stress 0.02520905 
    ## Run 59 stress 0.01879574 
    ## ... Procrustes: rmse 0.001145205  max resid 0.002359185 
    ## ... Similar to previous best
    ## Run 60 stress 0.01879575 
    ## ... Procrustes: rmse 0.001157053  max resid 0.002383519 
    ## ... Similar to previous best
    ## Run 61 stress 0.02486126 
    ## Run 62 stress 0.02509438 
    ## Run 63 stress 0.0188989 
    ## ... Procrustes: rmse 0.007785729  max resid 0.01578403 
    ## Run 64 stress 0.01880228 
    ## ... Procrustes: rmse 0.002339901  max resid 0.004789874 
    ## ... Similar to previous best
    ## Run 65 stress 0.02492164 
    ## Run 66 stress 0.01879811 
    ## ... Procrustes: rmse 0.001854267  max resid 0.003818514 
    ## ... Similar to previous best
    ## Run 67 stress 0.02493108 
    ## Run 68 stress 0.0250946 
    ## Run 69 stress 0.02493132 
    ## Run 70 stress 0.0195818 
    ## Run 71 stress 0.02509447 
    ## Run 72 stress 0.01879584 
    ## ... Procrustes: rmse 0.001201629  max resid 0.002475157 
    ## ... Similar to previous best
    ## Run 73 stress 0.01884947 
    ## ... Procrustes: rmse 0.004953412  max resid 0.009881429 
    ## ... Similar to previous best
    ## Run 74 stress 0.02492188 
    ## Run 75 stress 0.01884415 
    ## ... Procrustes: rmse 0.004490021  max resid 0.008906998 
    ## ... Similar to previous best
    ## Run 76 stress 0.01879595 
    ## ... Procrustes: rmse 0.001234294  max resid 0.00254254 
    ## ... Similar to previous best
    ## Run 77 stress 0.01879598 
    ## ... Procrustes: rmse 0.001247593  max resid 0.00256992 
    ## ... Similar to previous best
    ## Run 78 stress 0.01879585 
    ## ... Procrustes: rmse 0.0002365273  max resid 0.0004857876 
    ## ... Similar to previous best
    ## Run 79 stress 0.0250947 
    ## Run 80 stress 0.01879574 
    ## ... Procrustes: rmse 0.001154245  max resid 0.002377719 
    ## ... Similar to previous best
    ## Run 81 stress 0.01879526 
    ## ... New best solution
    ## ... Procrustes: rmse 7.493256e-05  max resid 0.0001545125 
    ## ... Similar to previous best
    ## Run 82 stress 0.01884474 
    ## ... Procrustes: rmse 0.004472696  max resid 0.008867211 
    ## ... Similar to previous best
    ## Run 83 stress 0.01879572 
    ## ... Procrustes: rmse 0.000242637  max resid 0.0004989792 
    ## ... Similar to previous best
    ## Run 84 stress 0.02520913 
    ## Run 85 stress 0.02486107 
    ## Run 86 stress 0.01879586 
    ## ... Procrustes: rmse 0.0003133931  max resid 0.0006444154 
    ## ... Similar to previous best
    ## Run 87 stress 0.01879714 
    ## ... Procrustes: rmse 0.001476077  max resid 0.003040155 
    ## ... Similar to previous best
    ## Run 88 stress 0.02509426 
    ## Run 89 stress 0.01884749 
    ## ... Procrustes: rmse 0.006048286  max resid 0.01243211 
    ## Run 90 stress 0.02486117 
    ## Run 91 stress 0.02492188 
    ## Run 92 stress 0.01879589 
    ## ... Procrustes: rmse 0.001144814  max resid 0.002357934 
    ## ... Similar to previous best
    ## Run 93 stress 0.01879741 
    ## ... Procrustes: rmse 0.001616141  max resid 0.003328153 
    ## ... Similar to previous best
    ## Run 94 stress 0.0249234 
    ## Run 95 stress 0.01885353 
    ## ... Procrustes: rmse 0.006220364  max resid 0.01278676 
    ## Run 96 stress 0.01879588 
    ## ... Procrustes: rmse 0.0003236369  max resid 0.0006654459 
    ## ... Similar to previous best
    ## Run 97 stress 0.02486101 
    ## Run 98 stress 0.02486109 
    ## Run 99 stress 0.02520903 
    ## Run 100 stress 0.0249219 
    ## Run 101 stress 0.02492188 
    ## Run 102 stress 0.02495517 
    ## Run 103 stress 0.02520881 
    ## Run 104 stress 0.01879577 
    ## ... Procrustes: rmse 0.0002767888  max resid 0.0005692618 
    ## ... Similar to previous best
    ## Run 105 stress 0.02520875 
    ## Run 106 stress 0.01906834 
    ## ... Procrustes: rmse 0.01236941  max resid 0.02521604 
    ## Run 107 stress 0.02492167 
    ## Run 108 stress 0.01880238 
    ## ... Procrustes: rmse 0.00226936  max resid 0.0046338 
    ## ... Similar to previous best
    ## Run 109 stress 0.01879575 
    ## ... Procrustes: rmse 0.001087037  max resid 0.002239101 
    ## ... Similar to previous best
    ## Run 110 stress 0.0252089 
    ## Run 111 stress 0.01880743 
    ## ... Procrustes: rmse 0.002283914  max resid 0.004703103 
    ## ... Similar to previous best
    ## Run 112 stress 0.02486113 
    ## Run 113 stress 0.01879977 
    ## ... Procrustes: rmse 0.002074599  max resid 0.004273138 
    ## ... Similar to previous best
    ## Run 114 stress 0.02509443 
    ## Run 115 stress 0.02509456 
    ## Run 116 stress 0.02509431 
    ## Run 117 stress 0.01882157 
    ## ... Procrustes: rmse 0.0044071  max resid 0.009068751 
    ## ... Similar to previous best
    ## Run 118 stress 0.01879596 
    ## ... Procrustes: rmse 0.001173541  max resid 0.002416533 
    ## ... Similar to previous best
    ## Run 119 stress 0.01880906 
    ## ... Procrustes: rmse 0.003324937  max resid 0.006846107 
    ## ... Similar to previous best
    ## Run 120 stress 0.01879585 
    ## ... Procrustes: rmse 0.001121295  max resid 0.00230941 
    ## ... Similar to previous best
    ## Run 121 stress 0.01886534 
    ## ... Procrustes: rmse 0.005967007  max resid 0.01200054 
    ## Run 122 stress 0.01879572 
    ## ... Procrustes: rmse 0.001071552  max resid 0.002207203 
    ## ... Similar to previous best
    ## Run 123 stress 0.02486121 
    ## Run 124 stress 0.01879539 
    ## ... Procrustes: rmse 0.0009014229  max resid 0.001856901 
    ## ... Similar to previous best
    ## Run 125 stress 0.02492192 
    ## Run 126 stress 0.2237017 
    ## Run 127 stress 0.01879673 
    ## ... Procrustes: rmse 0.001276582  max resid 0.002628641 
    ## ... Similar to previous best
    ## Run 128 stress 0.01879585 
    ## ... Procrustes: rmse 0.001118294  max resid 0.002303221 
    ## ... Similar to previous best
    ## Run 129 stress 0.02509473 
    ## Run 130 stress 0.02487666 
    ## Run 131 stress 0.01889774 
    ## ... Procrustes: rmse 0.007692674  max resid 0.0155888 
    ## Run 132 stress 0.01883636 
    ## ... Procrustes: rmse 0.003665042  max resid 0.007158692 
    ## ... Similar to previous best
    ## Run 133 stress 0.01899638 
    ## ... Procrustes: rmse 0.01131787  max resid 0.02307501 
    ## Run 134 stress 0.02520893 
    ## Run 135 stress 0.02509461 
    ## Run 136 stress 0.01879685 
    ## ... Procrustes: rmse 0.001282466  max resid 0.002641408 
    ## ... Similar to previous best
    ## Run 137 stress 0.0187958 
    ## ... Procrustes: rmse 0.001107977  max resid 0.002282125 
    ## ... Similar to previous best
    ## Run 138 stress 0.02486131 
    ## Run 139 stress 0.01879593 
    ## ... Procrustes: rmse 0.0003416045  max resid 0.0007023374 
    ## ... Similar to previous best
    ## Run 140 stress 0.01879547 
    ## ... Procrustes: rmse 0.0001292767  max resid 0.0002664257 
    ## ... Similar to previous best
    ## Run 141 stress 0.02492212 
    ## Run 142 stress 0.02509444 
    ## Run 143 stress 0.01879787 
    ## ... Procrustes: rmse 0.001704315  max resid 0.003509155 
    ## ... Similar to previous best
    ## Run 144 stress 0.01879576 
    ## ... Procrustes: rmse 0.001090397  max resid 0.002245987 
    ## ... Similar to previous best
    ## Run 145 stress 0.01879604 
    ## ... Procrustes: rmse 0.001187703  max resid 0.002445967 
    ## ... Similar to previous best
    ## Run 146 stress 0.0250945 
    ## Run 147 stress 0.01891505 
    ## ... Procrustes: rmse 0.008794437  max resid 0.01805938 
    ## Run 148 stress 0.02509441 
    ## Run 149 stress 0.01879593 
    ## ... Procrustes: rmse 0.001161977  max resid 0.00239336 
    ## ... Similar to previous best
    ## Run 150 stress 0.01881731 
    ## ... Procrustes: rmse 0.002367408  max resid 0.004753639 
    ## ... Similar to previous best
    ## Run 151 stress 0.02486123 
    ## Run 152 stress 0.01879875 
    ## ... Procrustes: rmse 0.001896878  max resid 0.003906 
    ## ... Similar to previous best
    ## Run 153 stress 0.02486115 
    ## Run 154 stress 0.0187959 
    ## ... Procrustes: rmse 0.001142713  max resid 0.002353786 
    ## ... Similar to previous best
    ## Run 155 stress 0.01879874 
    ## ... Procrustes: rmse 0.001818683  max resid 0.003745856 
    ## ... Similar to previous best
    ## Run 156 stress 0.02509436 
    ## Run 157 stress 0.01883746 
    ## ... Procrustes: rmse 0.003787593  max resid 0.007417805 
    ## ... Similar to previous best
    ## Run 158 stress 0.01879574 
    ## ... Procrustes: rmse 0.001084392  max resid 0.002233623 
    ## ... Similar to previous best
    ## Run 159 stress 0.02492189 
    ## Run 160 stress 0.0249218 
    ## Run 161 stress 0.01886705 
    ## ... Procrustes: rmse 0.00704181  max resid 0.01446921 
    ## Run 162 stress 0.02509459 
    ## Run 163 stress 0.02486122 
    ## Run 164 stress 0.01879564 
    ## ... Procrustes: rmse 0.00103556  max resid 0.002133082 
    ## ... Similar to previous best
    ## Run 165 stress 0.02509471 
    ## Run 166 stress 0.01925728 
    ## ... Procrustes: rmse 0.01750205  max resid 0.0359668 
    ## Run 167 stress 0.01879579 
    ## ... Procrustes: rmse 0.000285453  max resid 0.0005870492 
    ## ... Similar to previous best
    ## Run 168 stress 0.02509467 
    ## Run 169 stress 0.02486135 
    ## Run 170 stress 0.02000181 
    ## Run 171 stress 0.01879605 
    ## ... Procrustes: rmse 0.0003763299  max resid 0.000773547 
    ## ... Similar to previous best
    ## Run 172 stress 0.01879599 
    ## ... Procrustes: rmse 0.001169801  max resid 0.002409559 
    ## ... Similar to previous best
    ## Run 173 stress 0.01879557 
    ## ... Procrustes: rmse 0.001002443  max resid 0.002064908 
    ## ... Similar to previous best
    ## Run 174 stress 0.02486121 
    ## Run 175 stress 0.02509476 
    ## Run 176 stress 0.0187954 
    ## ... Procrustes: rmse 0.0009071969  max resid 0.001868918 
    ## ... Similar to previous best
    ## Run 177 stress 0.01879582 
    ## ... Procrustes: rmse 0.001114687  max resid 0.002296015 
    ## ... Similar to previous best
    ## Run 178 stress 0.02520876 
    ## Run 179 stress 0.01879709 
    ## ... Procrustes: rmse 0.001532162  max resid 0.003155473 
    ## ... Similar to previous best
    ## Run 180 stress 0.02509462 
    ## Run 181 stress 0.01879546 
    ## ... Procrustes: rmse 0.0001248778  max resid 0.0002572438 
    ## ... Similar to previous best
    ## Run 182 stress 0.01879572 
    ## ... Procrustes: rmse 0.001071771  max resid 0.002207639 
    ## ... Similar to previous best
    ## Run 183 stress 0.01883151 
    ## ... Procrustes: rmse 0.003018448  max resid 0.00578543 
    ## ... Similar to previous best
    ## Run 184 stress 0.01879559 
    ## ... Procrustes: rmse 0.00101195  max resid 0.002084527 
    ## ... Similar to previous best
    ## Run 185 stress 0.02509473 
    ## Run 186 stress 0.3743875 
    ## Run 187 stress 0.02486121 
    ## Run 188 stress 0.01879799 
    ## ... Procrustes: rmse 0.00175127  max resid 0.00360638 
    ## ... Similar to previous best
    ## Run 189 stress 0.01879769 
    ## ... Procrustes: rmse 0.001682449  max resid 0.003464705 
    ## ... Similar to previous best
    ## Run 190 stress 0.01879572 
    ## ... Procrustes: rmse 0.0002489124  max resid 0.000511998 
    ## ... Similar to previous best
    ## Run 191 stress 0.02492182 
    ## Run 192 stress 0.01879597 
    ## ... Procrustes: rmse 0.001161187  max resid 0.002391411 
    ## ... Similar to previous best
    ## Run 193 stress 0.02492178 
    ## Run 194 stress 0.01879583 
    ## ... Procrustes: rmse 0.001112304  max resid 0.002290887 
    ## ... Similar to previous best
    ## Run 195 stress 0.02486124 
    ## Run 196 stress 0.01879586 
    ## ... Procrustes: rmse 0.001132079  max resid 0.002331821 
    ## ... Similar to previous best
    ## Run 197 stress 0.02494974 
    ## Run 198 stress 0.01879592 
    ## ... Procrustes: rmse 0.001148058  max resid 0.002364807 
    ## ... Similar to previous best
    ## Run 199 stress 0.02509472 
    ## Run 200 stress 0.0187956 
    ## ... Procrustes: rmse 0.001016366  max resid 0.002093737 
    ## ... Similar to previous best
    ## Run 201 stress 0.01879567 
    ## ... Procrustes: rmse 0.0002326255  max resid 0.0004785847 
    ## ... Similar to previous best
    ## Run 202 stress 0.01879575 
    ## ... Procrustes: rmse 0.0002580647  max resid 0.000530669 
    ## ... Similar to previous best
    ## Run 203 stress 0.02492186 
    ## Run 204 stress 0.01880077 
    ## ... Procrustes: rmse 0.002273829  max resid 0.004682926 
    ## ... Similar to previous best
    ## Run 205 stress 0.01879546 
    ## ... Procrustes: rmse 0.0001200228  max resid 0.0002473292 
    ## ... Similar to previous best
    ## Run 206 stress 0.01879567 
    ## ... Procrustes: rmse 0.001049229  max resid 0.002161236 
    ## ... Similar to previous best
    ## Run 207 stress 0.02509463 
    ## Run 208 stress 0.01884492 
    ## ... Procrustes: rmse 0.004504951  max resid 0.008935824 
    ## ... Similar to previous best
    ## Run 209 stress 0.01896211 
    ## ... Procrustes: rmse 0.0104498  max resid 0.02142736 
    ## Run 210 stress 0.02509481 
    ## Run 211 stress 0.0187958 
    ## ... Procrustes: rmse 0.001108119  max resid 0.002282473 
    ## ... Similar to previous best
    ## Run 212 stress 0.01879601 
    ## ... Procrustes: rmse 0.001172456  max resid 0.002415044 
    ## ... Similar to previous best
    ## Run 213 stress 0.01879573 
    ## ... Procrustes: rmse 0.001061927  max resid 0.002187533 
    ## ... Similar to previous best
    ## Run 214 stress 0.01879566 
    ## ... Procrustes: rmse 0.001036741  max resid 0.002135342 
    ## ... Similar to previous best
    ## Run 215 stress 0.01879593 
    ## ... Procrustes: rmse 0.001162938  max resid 0.002395273 
    ## ... Similar to previous best
    ## Run 216 stress 0.01879588 
    ## ... Procrustes: rmse 0.001131585  max resid 0.002330554 
    ## ... Similar to previous best
    ## Run 217 stress 0.02509467 
    ## Run 218 stress 0.01879567 
    ## ... Procrustes: rmse 0.0002300798  max resid 0.0004733676 
    ## ... Similar to previous best
    ## Run 219 stress 0.02509483 
    ## Run 220 stress 0.01904821 
    ## ... Procrustes: rmse 0.01279352  max resid 0.02611358 
    ## Run 221 stress 0.01879581 
    ## ... Procrustes: rmse 0.001096039  max resid 0.002257737 
    ## ... Similar to previous best
    ## Run 222 stress 0.02492195 
    ## Run 223 stress 0.01880952 
    ## ... Procrustes: rmse 0.003372416  max resid 0.00694361 
    ## ... Similar to previous best
    ## Run 224 stress 0.018796 
    ## ... Procrustes: rmse 0.001182393  max resid 0.002435149 
    ## ... Similar to previous best
    ## Run 225 stress 0.01879575 
    ## ... Procrustes: rmse 0.001083092  max resid 0.002231038 
    ## ... Similar to previous best
    ## Run 226 stress 0.02509452 
    ## Run 227 stress 0.01879574 
    ## ... Procrustes: rmse 0.001081963  max resid 0.002228702 
    ## ... Similar to previous best
    ## Run 228 stress 0.01879595 
    ## ... Procrustes: rmse 0.001152873  max resid 0.00237474 
    ## ... Similar to previous best
    ## Run 229 stress 0.0192251 
    ## ... Procrustes: rmse 0.01607498  max resid 0.03304461 
    ## Run 230 stress 0.01881028 
    ## ... Procrustes: rmse 0.003449768  max resid 0.007102548 
    ## ... Similar to previous best
    ## Run 231 stress 0.02492192 
    ## Run 232 stress 0.01879563 
    ## ... Procrustes: rmse 0.001032843  max resid 0.002127545 
    ## ... Similar to previous best
    ## Run 233 stress 0.01894914 
    ## ... Procrustes: rmse 0.009761996  max resid 0.01987118 
    ## Run 234 stress 0.02509473 
    ## Run 235 stress 0.01891782 
    ## ... Procrustes: rmse 0.008415826  max resid 0.01708329 
    ## Run 236 stress 0.02520898 
    ## Run 237 stress 0.02492159 
    ## Run 238 stress 0.01895785 
    ## ... Procrustes: rmse 0.01044259  max resid 0.02142023 
    ## Run 239 stress 0.02486145 
    ## Run 240 stress 0.01885653 
    ## ... Procrustes: rmse 0.005401038  max resid 0.01081877 
    ## Run 241 stress 0.01892552 
    ## ... Procrustes: rmse 0.009363657  max resid 0.01921976 
    ## Run 242 stress 0.02486131 
    ## Run 243 stress 0.01883452 
    ## ... Procrustes: rmse 0.003395547  max resid 0.006584374 
    ## ... Similar to previous best
    ## Run 244 stress 0.0250947 
    ## Run 245 stress 0.01879575 
    ## ... Procrustes: rmse 0.001086555  max resid 0.002238134 
    ## ... Similar to previous best
    ## Run 246 stress 0.01879574 
    ## ... Procrustes: rmse 0.001079588  max resid 0.002223664 
    ## ... Similar to previous best
    ## Run 247 stress 0.0249219 
    ## Run 248 stress 0.02509472 
    ## Run 249 stress 0.01879583 
    ## ... Procrustes: rmse 0.00112255  max resid 0.002312187 
    ## ... Similar to previous best
    ## Run 250 stress 0.02509491 
    ## Run 251 stress 0.02509456 
    ## Run 252 stress 0.01879721 
    ## ... Procrustes: rmse 0.001358536  max resid 0.002798306 
    ## ... Similar to previous best
    ## Run 253 stress 0.02509444 
    ## Run 254 stress 0.0250945 
    ## Run 255 stress 0.01893422 
    ## ... Procrustes: rmse 0.009676581  max resid 0.01985735 
    ## Run 256 stress 0.01879552 
    ## ... Procrustes: rmse 0.0009751967  max resid 0.002008786 
    ## ... Similar to previous best
    ## Run 257 stress 0.01889803 
    ## ... Procrustes: rmse 0.008358704  max resid 0.01716556 
    ## Run 258 stress 0.01967164 
    ## Run 259 stress 0.01883271 
    ## ... Procrustes: rmse 0.003213268  max resid 0.006197865 
    ## ... Similar to previous best
    ## Run 260 stress 0.01879551 
    ## ... Procrustes: rmse 0.0009514458  max resid 0.00195972 
    ## ... Similar to previous best
    ## Run 261 stress 0.01999811 
    ## Run 262 stress 0.01879559 
    ## ... Procrustes: rmse 0.001000438  max resid 0.00206093 
    ## ... Similar to previous best
    ## Run 263 stress 0.01879856 
    ## ... Procrustes: rmse 0.001871539  max resid 0.003854372 
    ## ... Similar to previous best
    ## Run 264 stress 0.01879544 
    ## ... Procrustes: rmse 0.0009266859  max resid 0.001909873 
    ## ... Similar to previous best
    ## Run 265 stress 0.02509463 
    ## Run 266 stress 0.02509448 
    ## Run 267 stress 0.01879585 
    ## ... Procrustes: rmse 0.001124903  max resid 0.002316875 
    ## ... Similar to previous best
    ## Run 268 stress 0.02509476 
    ## Run 269 stress 0.01883611 
    ## ... Procrustes: rmse 0.00362298  max resid 0.007068282 
    ## ... Similar to previous best
    ## Run 270 stress 0.01879553 
    ## ... Procrustes: rmse 0.0009806001  max resid 0.002020034 
    ## ... Similar to previous best
    ## Run 271 stress 0.02520917 
    ## Run 272 stress 0.02492197 
    ## Run 273 stress 0.02509438 
    ## Run 274 stress 0.02509472 
    ## Run 275 stress 0.0188151 
    ## ... Procrustes: rmse 0.003898568  max resid 0.008025206 
    ## ... Similar to previous best
    ## Run 276 stress 0.01879597 
    ## ... Procrustes: rmse 0.001170505  max resid 0.002410976 
    ## ... Similar to previous best
    ## Run 277 stress 0.01879583 
    ## ... Procrustes: rmse 0.001101178  max resid 0.002267903 
    ## ... Similar to previous best
    ## Run 278 stress 0.3166934 
    ## Run 279 stress 0.01879891 
    ## ... Procrustes: rmse 0.001940884  max resid 0.003996955 
    ## ... Similar to previous best
    ## Run 280 stress 0.02509466 
    ## Run 281 stress 0.01883498 
    ## ... Procrustes: rmse 0.003503136  max resid 0.006812948 
    ## ... Similar to previous best
    ## Run 282 stress 0.01879583 
    ## ... Procrustes: rmse 0.001122303  max resid 0.002311679 
    ## ... Similar to previous best
    ## Run 283 stress 0.02494681 
    ## Run 284 stress 0.01879642 
    ## ... Procrustes: rmse 0.001337423  max resid 0.002754393 
    ## ... Similar to previous best
    ## Run 285 stress 0.01883086 
    ## ... Procrustes: rmse 0.002945099  max resid 0.005631266 
    ## ... Similar to previous best
    ## Run 286 stress 0.0188006 
    ## ... Procrustes: rmse 0.002233915  max resid 0.004600405 
    ## ... Similar to previous best
    ## Run 287 stress 0.02509476 
    ## Run 288 stress 0.223643 
    ## Run 289 stress 0.01879565 
    ## ... Procrustes: rmse 0.001040101  max resid 0.002142486 
    ## ... Similar to previous best
    ## Run 290 stress 0.02102324 
    ## Run 291 stress 0.0188117 
    ## ... Procrustes: rmse 0.003168277  max resid 0.006525049 
    ## ... Similar to previous best
    ## Run 292 stress 0.01879598 
    ## ... Procrustes: rmse 0.001169487  max resid 0.002408871 
    ## ... Similar to previous best
    ## Run 293 stress 0.02492447 
    ## Run 294 stress 0.02509474 
    ## Run 295 stress 0.01910931 
    ## ... Procrustes: rmse 0.01431527  max resid 0.02927604 
    ## Run 296 stress 0.02520876 
    ## Run 297 stress 0.02509442 
    ## Run 298 stress 0.01879593 
    ## ... Procrustes: rmse 0.0003423252  max resid 0.0007037935 
    ## ... Similar to previous best
    ## Run 299 stress 0.02509455 
    ## Run 300 stress 0.02520894 
    ## Run 301 stress 0.01887138 
    ## ... Procrustes: rmse 0.00622697  max resid 0.01254052 
    ## Run 302 stress 0.01879715 
    ## ... Procrustes: rmse 0.001548595  max resid 0.003189048 
    ## ... Similar to previous best
    ## Run 303 stress 0.02509473 
    ## Run 304 stress 0.01879611 
    ## ... Procrustes: rmse 0.001217427  max resid 0.002507594 
    ## ... Similar to previous best
    ## Run 305 stress 0.02509449 
    ## Run 306 stress 0.02492169 
    ## Run 307 stress 0.01879584 
    ## ... Procrustes: rmse 0.001094666  max resid 0.002255669 
    ## ... Similar to previous best
    ## Run 308 stress 0.02520861 
    ## Run 309 stress 0.0188573 
    ## ... Procrustes: rmse 0.006572481  max resid 0.01350999 
    ## Run 310 stress 0.02486123 
    ## Run 311 stress 0.02492205 
    ## Run 312 stress 0.01880086 
    ## ... Procrustes: rmse 0.001901379  max resid 0.003916143 
    ## ... Similar to previous best
    ## Run 313 stress 0.01879594 
    ## ... Procrustes: rmse 0.001150499  max resid 0.002369849 
    ## ... Similar to previous best
    ## Run 314 stress 0.01879567 
    ## ... Procrustes: rmse 0.001047808  max resid 0.002158408 
    ## ... Similar to previous best
    ## Run 315 stress 0.01879561 
    ## ... Procrustes: rmse 0.001018387  max resid 0.002098064 
    ## ... Similar to previous best
    ## Run 316 stress 0.0187959 
    ## ... Procrustes: rmse 0.0003313197  max resid 0.0006812553 
    ## ... Similar to previous best
    ## Run 317 stress 0.01879575 
    ## ... Procrustes: rmse 0.001087494  max resid 0.002239994 
    ## ... Similar to previous best
    ## Run 318 stress 0.01883142 
    ## ... Procrustes: rmse 0.00302451  max resid 0.005797008 
    ## ... Similar to previous best
    ## Run 319 stress 0.02520892 
    ## Run 320 stress 0.02486137 
    ## Run 321 stress 0.02509472 
    ## Run 322 stress 0.0187957 
    ## ... Procrustes: rmse 0.0002460806  max resid 0.0005062149 
    ## ... Similar to previous best
    ## Run 323 stress 0.02492194 
    ## Run 324 stress 0.2262188 
    ## Run 325 stress 0.02492191 
    ## Run 326 stress 0.01883388 
    ## ... Procrustes: rmse 0.005272717  max resid 0.01084296 
    ## Run 327 stress 0.02492193 
    ## Run 328 stress 0.02509461 
    ## Run 329 stress 0.02509476 
    ## Run 330 stress 0.0188746 
    ## ... Procrustes: rmse 0.006369285  max resid 0.01283349 
    ## Run 331 stress 0.01884563 
    ## ... Procrustes: rmse 0.004559885  max resid 0.009053416 
    ## ... Similar to previous best
    ## Run 332 stress 0.01882516 
    ## ... Procrustes: rmse 0.004686985  max resid 0.009641759 
    ## ... Similar to previous best
    ## Run 333 stress 0.02509448 
    ## Run 334 stress 0.0187963 
    ## ... Procrustes: rmse 0.001275263  max resid 0.002626675 
    ## ... Similar to previous best
    ## Run 335 stress 0.02492192 
    ## Run 336 stress 0.02492193 
    ## Run 337 stress 0.02492213 
    ## Run 338 stress 0.02509464 
    ## Run 339 stress 0.01883015 
    ## ... Procrustes: rmse 0.002804549  max resid 0.005343492 
    ## ... Similar to previous best
    ## Run 340 stress 0.0250948 
    ## Run 341 stress 0.02486128 
    ## Run 342 stress 0.01879586 
    ## ... Procrustes: rmse 0.00113148  max resid 0.002330582 
    ## ... Similar to previous best
    ## Run 343 stress 0.01901419 
    ## ... Procrustes: rmse 0.0116515  max resid 0.02375415 
    ## Run 344 stress 0.0188127 
    ## ... Procrustes: rmse 0.003682415  max resid 0.007580364 
    ## ... Similar to previous best
    ## Run 345 stress 0.02509454 
    ## Run 346 stress 0.01879537 
    ## ... Procrustes: rmse 0.0008885717  max resid 0.001830509 
    ## ... Similar to previous best
    ## Run 347 stress 0.01879563 
    ## ... Procrustes: rmse 0.001029729  max resid 0.002121132 
    ## ... Similar to previous best
    ## Run 348 stress 0.01879538 
    ## ... Procrustes: rmse 0.0008926807  max resid 0.001838945 
    ## ... Similar to previous best
    ## Run 349 stress 0.0188317 
    ## ... Procrustes: rmse 0.003070782  max resid 0.005895036 
    ## ... Similar to previous best
    ## Run 350 stress 0.0187955 
    ## ... Procrustes: rmse 0.0009529272  max resid 0.001962792 
    ## ... Similar to previous best
    ## Run 351 stress 0.02486124 
    ## Run 352 stress 0.02509481 
    ## Run 353 stress 0.02486128 
    ## Run 354 stress 0.01886294 
    ## ... Procrustes: rmse 0.006850696  max resid 0.01407957 
    ## Run 355 stress 0.02509442 
    ## Run 356 stress 0.02492166 
    ## Run 357 stress 0.01879577 
    ## ... Procrustes: rmse 0.0002738874  max resid 0.0005632293 
    ## ... Similar to previous best
    ## Run 358 stress 0.0187956 
    ## ... Procrustes: rmse 0.001014937  max resid 0.00209073 
    ## ... Similar to previous best
    ## Run 359 stress 0.02509489 
    ## Run 360 stress 0.01880428 
    ## ... Procrustes: rmse 0.002775769  max resid 0.005716573 
    ## ... Similar to previous best
    ## Run 361 stress 0.02509444 
    ## Run 362 stress 0.02509487 
    ## Run 363 stress 0.02520904 
    ## Run 364 stress 0.01881087 
    ## ... Procrustes: rmse 0.003506081  max resid 0.007218036 
    ## ... Similar to previous best
    ## Run 365 stress 0.02509474 
    ## Run 366 stress 0.01879537 
    ## ... Procrustes: rmse 0.0008891614  max resid 0.00183182 
    ## ... Similar to previous best
    ## Run 367 stress 0.01879566 
    ## ... Procrustes: rmse 0.001044239  max resid 0.002150995 
    ## ... Similar to previous best
    ## Run 368 stress 0.0250514 
    ## Run 369 stress 0.01879849 
    ## ... Procrustes: rmse 0.001853232  max resid 0.003816694 
    ## ... Similar to previous best
    ## Run 370 stress 0.0187956 
    ## ... Procrustes: rmse 0.001015813  max resid 0.002092434 
    ## ... Similar to previous best
    ## Run 371 stress 0.02509435 
    ## Run 372 stress 0.01879564 
    ## ... Procrustes: rmse 0.001034532  max resid 0.002131009 
    ## ... Similar to previous best
    ## Run 373 stress 0.02520886 
    ## Run 374 stress 0.02509477 
    ## Run 375 stress 0.02509443 
    ## Run 376 stress 0.01879595 
    ## ... Procrustes: rmse 0.00116103  max resid 0.002391496 
    ## ... Similar to previous best
    ## Run 377 stress 0.01882463 
    ## ... Procrustes: rmse 0.004636909  max resid 0.009539971 
    ## ... Similar to previous best
    ## Run 378 stress 0.0188316 
    ## ... Procrustes: rmse 0.003036169  max resid 0.005823067 
    ## ... Similar to previous best
    ## Run 379 stress 0.02509467 
    ## Run 380 stress 0.02492204 
    ## Run 381 stress 0.01879593 
    ## ... Procrustes: rmse 0.00114882  max resid 0.002366009 
    ## ... Similar to previous best
    ## Run 382 stress 0.01885996 
    ## ... Procrustes: rmse 0.006709203  max resid 0.01378597 
    ## Run 383 stress 0.0187958 
    ## ... Procrustes: rmse 0.000289419  max resid 0.0005951709 
    ## ... Similar to previous best
    ## Run 384 stress 0.01883407 
    ## ... Procrustes: rmse 0.003390619  max resid 0.006573974 
    ## ... Similar to previous best
    ## Run 385 stress 0.0187959 
    ## ... Procrustes: rmse 0.0003175415  max resid 0.0006527951 
    ## ... Similar to previous best
    ## Run 386 stress 0.02509443 
    ## Run 387 stress 0.01879568 
    ## ... Procrustes: rmse 0.0002338223  max resid 0.0004818298 
    ## ... Similar to previous best
    ## Run 388 stress 0.01879546 
    ## ... Procrustes: rmse 0.0009436359  max resid 0.001943877 
    ## ... Similar to previous best
    ## Run 389 stress 0.02509448 
    ## Run 390 stress 0.01879547 
    ## ... Procrustes: rmse 0.0009496124  max resid 0.001956211 
    ## ... Similar to previous best
    ## Run 391 stress 0.02509466 
    ## Run 392 stress 0.01879582 
    ## ... Procrustes: rmse 0.0002997182  max resid 0.000616335 
    ## ... Similar to previous best
    ## Run 393 stress 0.02509472 
    ## Run 394 stress 0.02520883 
    ## Run 395 stress 0.02509454 
    ## Run 396 stress 0.0250947 
    ## Run 397 stress 0.0187966 
    ## ... Procrustes: rmse 0.001392781  max resid 0.002868345 
    ## ... Similar to previous best
    ## Run 398 stress 0.02492179 
    ## Run 399 stress 0.01884857 
    ## ... Procrustes: rmse 0.004814304  max resid 0.009588342 
    ## ... Similar to previous best
    ## Run 400 stress 0.01894971 
    ## ... Procrustes: rmse 0.009769035  max resid 0.01988604 
    ## Run 401 stress 0.02509453 
    ## Run 402 stress 0.02486104 
    ## Run 403 stress 0.0248613 
    ## Run 404 stress 0.01885485 
    ## ... Procrustes: rmse 0.005271542  max resid 0.01054867 
    ## Run 405 stress 0.01879676 
    ## ... Procrustes: rmse 0.00144175  max resid 0.002969116 
    ## ... Similar to previous best
    ## Run 406 stress 0.01883897 
    ## ... Procrustes: rmse 0.003946538  max resid 0.007754369 
    ## ... Similar to previous best
    ## Run 407 stress 0.0249474 
    ## Run 408 stress 0.01879579 
    ## ... Procrustes: rmse 0.001106098  max resid 0.002278337 
    ## ... Similar to previous best
    ## Run 409 stress 0.01879581 
    ## ... Procrustes: rmse 0.001106101  max resid 0.002278186 
    ## ... Similar to previous best
    ## Run 410 stress 0.02509469 
    ## Run 411 stress 0.02509435 
    ## Run 412 stress 0.02492181 
    ## Run 413 stress 0.02492191 
    ## Run 414 stress 0.01879609 
    ## ... Procrustes: rmse 0.001212372  max resid 0.002497178 
    ## ... Similar to previous best
    ## Run 415 stress 0.01879575 
    ## ... Procrustes: rmse 0.0002641763  max resid 0.0005430518 
    ## ... Similar to previous best
    ## Run 416 stress 0.02492199 
    ## Run 417 stress 0.02486144 
    ## Run 418 stress 0.01884356 
    ## ... Procrustes: rmse 0.005851281  max resid 0.01203032 
    ## Run 419 stress 0.02509482 
    ## Run 420 stress 0.02492184 
    ## Run 421 stress 0.01879587 
    ## ... Procrustes: rmse 0.001138388  max resid 0.00234474 
    ## ... Similar to previous best
    ## Run 422 stress 0.02509468 
    ## Run 423 stress 0.01880517 
    ## ... Procrustes: rmse 0.002837373  max resid 0.005841878 
    ## ... Similar to previous best
    ## Run 424 stress 0.01883545 
    ## ... Procrustes: rmse 0.003559825  max resid 0.00693312 
    ## ... Similar to previous best
    ## Run 425 stress 0.02486126 
    ## Run 426 stress 0.02486138 
    ## Run 427 stress 0.02492205 
    ## Run 428 stress 0.02509432 
    ## Run 429 stress 0.02492197 
    ## Run 430 stress 0.02492573 
    ## Run 431 stress 0.0249218 
    ## Run 432 stress 0.01879571 
    ## ... Procrustes: rmse 0.001070946  max resid 0.002205964 
    ## ... Similar to previous best
    ## Run 433 stress 0.02509451 
    ## Run 434 stress 0.02486109 
    ## Run 435 stress 0.01883751 
    ## ... Procrustes: rmse 0.003789817  max resid 0.007423113 
    ## ... Similar to previous best
    ## Run 436 stress 0.0187957 
    ## ... Procrustes: rmse 0.0002454042  max resid 0.00050483 
    ## ... Similar to previous best
    ## Run 437 stress 0.02509459 
    ## Run 438 stress 0.01883648 
    ## ... Procrustes: rmse 0.003679165  max resid 0.007186953 
    ## ... Similar to previous best
    ## Run 439 stress 0.02509442 
    ## Run 440 stress 0.01885784 
    ## ... Procrustes: rmse 0.005495891  max resid 0.0110162 
    ## Run 441 stress 0.01908905 
    ## ... Procrustes: rmse 0.0136702  max resid 0.02798187 
    ## Run 442 stress 0.02492174 
    ## Run 443 stress 0.02509435 
    ## Run 444 stress 0.02509444 
    ## Run 445 stress 0.02509465 
    ## Run 446 stress 0.01884275 
    ## ... Procrustes: rmse 0.004314566  max resid 0.008534042 
    ## ... Similar to previous best
    ## Run 447 stress 0.01901759 
    ## ... Procrustes: rmse 0.0120254  max resid 0.02465492 
    ## Run 448 stress 0.02492183 
    ## Run 449 stress 0.01880445 
    ## ... Procrustes: rmse 0.002798305  max resid 0.00576287 
    ## ... Similar to previous best
    ## Run 450 stress 0.01879592 
    ## ... Procrustes: rmse 0.001138051  max resid 0.002343793 
    ## ... Similar to previous best
    ## Run 451 stress 0.01888063 
    ## ... Procrustes: rmse 0.006683243  max resid 0.01348952 
    ## Run 452 stress 0.02492194 
    ## Run 453 stress 0.01892925 
    ## ... Procrustes: rmse 0.009022313  max resid 0.01834424 
    ## Run 454 stress 0.01884167 
    ## ... Procrustes: rmse 0.004022314  max resid 0.007913555 
    ## ... Similar to previous best
    ## Run 455 stress 0.02486123 
    ## Run 456 stress 0.01879598 
    ## ... Procrustes: rmse 0.001173053  max resid 0.002416182 
    ## ... Similar to previous best
    ## Run 457 stress 0.02509475 
    ## Run 458 stress 0.01879534 
    ## ... Procrustes: rmse 0.000873417  max resid 0.001799301 
    ## ... Similar to previous best
    ## Run 459 stress 0.02509459 
    ## Run 460 stress 0.01879587 
    ## ... Procrustes: rmse 0.001131044  max resid 0.002329476 
    ## ... Similar to previous best
    ## Run 461 stress 0.02486139 
    ## Run 462 stress 0.01885661 
    ## ... Procrustes: rmse 0.004959676  max resid 0.0101671 
    ## Run 463 stress 0.02492213 
    ## Run 464 stress 0.01879544 
    ## ... Procrustes: rmse 0.0009335955  max resid 0.001923173 
    ## ... Similar to previous best
    ## Run 465 stress 0.02486128 
    ## Run 466 stress 0.0188315 
    ## ... Procrustes: rmse 0.003040145  max resid 0.005830413 
    ## ... Similar to previous best
    ## Run 467 stress 0.02486123 
    ## Run 468 stress 0.02492182 
    ## Run 469 stress 0.01879589 
    ## ... Procrustes: rmse 0.001130141  max resid 0.00232755 
    ## ... Similar to previous best
    ## Run 470 stress 0.01879583 
    ## ... Procrustes: rmse 0.001121673  max resid 0.002310304 
    ## ... Similar to previous best
    ## Run 471 stress 0.01881732 
    ## ... Procrustes: rmse 0.004067582  max resid 0.008369645 
    ## ... Similar to previous best
    ## Run 472 stress 0.02509468 
    ## Run 473 stress 0.01880136 
    ## ... Procrustes: rmse 0.002367882  max resid 0.004876762 
    ## ... Similar to previous best
    ## Run 474 stress 0.01879585 
    ## ... Procrustes: rmse 0.001131673  max resid 0.002330945 
    ## ... Similar to previous best
    ## Run 475 stress 0.02509424 
    ## Run 476 stress 0.024861 
    ## Run 477 stress 0.02509439 
    ## Run 478 stress 0.02520892 
    ## Run 479 stress 0.01900852 
    ## ... Procrustes: rmse 0.01191783  max resid 0.02442603 
    ## Run 480 stress 0.01889462 
    ## ... Procrustes: rmse 0.006893373  max resid 0.01413445 
    ## Run 481 stress 0.01880605 
    ## ... Procrustes: rmse 0.002987959  max resid 0.006152646 
    ## ... Similar to previous best
    ## Run 482 stress 0.01879584 
    ## ... Procrustes: rmse 0.001126429  max resid 0.002320132 
    ## ... Similar to previous best
    ## Run 483 stress 0.01882744 
    ## ... Procrustes: rmse 0.004841875  max resid 0.009960145 
    ## ... Similar to previous best
    ## Run 484 stress 0.02509484 
    ## Run 485 stress 0.01879601 
    ## ... Procrustes: rmse 0.001184781  max resid 0.002440375 
    ## ... Similar to previous best
    ## Run 486 stress 0.02486124 
    ## Run 487 stress 0.01879578 
    ## ... Procrustes: rmse 0.001100331  max resid 0.002266459 
    ## ... Similar to previous best
    ## Run 488 stress 0.0249218 
    ## Run 489 stress 0.02492203 
    ## Run 490 stress 0.0250943 
    ## Run 491 stress 0.01879564 
    ## ... Procrustes: rmse 0.001032574  max resid 0.002126861 
    ## ... Similar to previous best
    ## Run 492 stress 0.02509449 
    ## Run 493 stress 0.01913885 
    ## ... Procrustes: rmse 0.01497499  max resid 0.03074519 
    ## Run 494 stress 0.01879572 
    ## ... Procrustes: rmse 0.0002537762  max resid 0.0005220225 
    ## ... Similar to previous best
    ## Run 495 stress 0.02509443 
    ## Run 496 stress 0.02492191 
    ## Run 497 stress 0.02486128 
    ## Run 498 stress 0.01885401 
    ## ... Procrustes: rmse 0.005204498  max resid 0.0104062 
    ## Run 499 stress 0.01879997 
    ## ... Procrustes: rmse 0.00207654  max resid 0.004277297 
    ## ... Similar to previous best
    ## Run 500 stress 0.0248612 
    ## *** Best solution repeated 183 times

``` r
# Mixed and stratified lakes
SD_beta_MS_NMDS <- metaMDS(SD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09416137 
    ## Run 2 stress 0.09969988 
    ## Run 3 stress 0.08773465 
    ## Run 4 stress 0.09268347 
    ## Run 5 stress 0.09337284 
    ## Run 6 stress 0.09030398 
    ## Run 7 stress 0.09712942 
    ## Run 8 stress 0.09407998 
    ## Run 9 stress 0.09337246 
    ## Run 10 stress 0.08503618 
    ## ... Procrustes: rmse 0.0001449861  max resid 0.0003546182 
    ## ... Similar to previous best
    ## Run 11 stress 0.0940742 
    ## Run 12 stress 0.09030412 
    ## Run 13 stress 0.09321616 
    ## Run 14 stress 0.08503538 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003625486  max resid 0.001028975 
    ## ... Similar to previous best
    ## Run 15 stress 0.08503484 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001273448  max resid 0.003063582 
    ## ... Similar to previous best
    ## Run 16 stress 0.09760831 
    ## Run 17 stress 0.09407968 
    ## Run 18 stress 0.09159099 
    ## Run 19 stress 0.08503731 
    ## ... Procrustes: rmse 0.002061867  max resid 0.00513876 
    ## ... Similar to previous best
    ## Run 20 stress 0.09030394 
    ## Run 21 stress 0.08773467 
    ## Run 22 stress 0.08503481 
    ## ... New best solution
    ## ... Procrustes: rmse 6.053401e-05  max resid 0.0001166795 
    ## ... Similar to previous best
    ## Run 23 stress 0.09407994 
    ## Run 24 stress 0.09400519 
    ## Run 25 stress 0.09268386 
    ## Run 26 stress 0.09407982 
    ## Run 27 stress 0.09030403 
    ## Run 28 stress 0.3156497 
    ## Run 29 stress 0.09159088 
    ## Run 30 stress 0.09145322 
    ## Run 31 stress 0.09417766 
    ## Run 32 stress 0.09447297 
    ## Run 33 stress 0.0914534 
    ## Run 34 stress 0.0897387 
    ## Run 35 stress 0.09337301 
    ## Run 36 stress 0.08773489 
    ## Run 37 stress 0.09159083 
    ## Run 38 stress 0.09030399 
    ## Run 39 stress 0.08503467 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001383789  max resid 0.0002802602 
    ## ... Similar to previous best
    ## Run 40 stress 0.09407998 
    ## Run 41 stress 0.08503639 
    ## ... Procrustes: rmse 0.001706559  max resid 0.004304975 
    ## ... Similar to previous best
    ## Run 42 stress 0.09308946 
    ## Run 43 stress 0.09407968 
    ## Run 44 stress 0.08503463 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003233067  max resid 0.0007878932 
    ## ... Similar to previous best
    ## Run 45 stress 0.2567613 
    ## Run 46 stress 0.09465904 
    ## Run 47 stress 0.09337267 
    ## Run 48 stress 0.09465905 
    ## Run 49 stress 0.08503481 
    ## ... Procrustes: rmse 0.0004445492  max resid 0.001010372 
    ## ... Similar to previous best
    ## Run 50 stress 0.09380599 
    ## Run 51 stress 0.08440261 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01417467  max resid 0.04272789 
    ## Run 52 stress 0.09407993 
    ## Run 53 stress 0.08503508 
    ## Run 54 stress 0.09337261 
    ## Run 55 stress 0.09030409 
    ## Run 56 stress 0.0996997 
    ## Run 57 stress 0.1052851 
    ## Run 58 stress 0.09760813 
    ## Run 59 stress 0.08773466 
    ## Run 60 stress 0.1072866 
    ## Run 61 stress 0.08973865 
    ## Run 62 stress 0.08973863 
    ## Run 63 stress 0.09286083 
    ## Run 64 stress 0.09407978 
    ## Run 65 stress 0.08503521 
    ## Run 66 stress 0.08503544 
    ## Run 67 stress 0.09030412 
    ## Run 68 stress 0.09286088 
    ## Run 69 stress 0.09030406 
    ## Run 70 stress 0.08973903 
    ## Run 71 stress 0.08773468 
    ## Run 72 stress 0.09407973 
    ## Run 73 stress 0.09381853 
    ## Run 74 stress 0.09408004 
    ## Run 75 stress 0.08773485 
    ## Run 76 stress 0.08440272 
    ## ... Procrustes: rmse 8.079507e-05  max resid 0.0001448407 
    ## ... Similar to previous best
    ## Run 77 stress 0.08973883 
    ## Run 78 stress 0.09374409 
    ## Run 79 stress 0.09286083 
    ## Run 80 stress 0.09268328 
    ## Run 81 stress 0.09400538 
    ## Run 82 stress 0.08503471 
    ## Run 83 stress 0.09407983 
    ## Run 84 stress 0.09465905 
    ## Run 85 stress 0.09268385 
    ## Run 86 stress 0.09030411 
    ## Run 87 stress 0.08503486 
    ## Run 88 stress 0.09030405 
    ## Run 89 stress 0.09416134 
    ## Run 90 stress 0.09465905 
    ## Run 91 stress 0.09268345 
    ## Run 92 stress 0.09969985 
    ## Run 93 stress 0.09030405 
    ## Run 94 stress 0.09268336 
    ## Run 95 stress 0.09969991 
    ## Run 96 stress 0.0941614 
    ## Run 97 stress 0.09145325 
    ## Run 98 stress 0.09268419 
    ## Run 99 stress 0.09465911 
    ## Run 100 stress 0.08973877 
    ## Run 101 stress 0.09286086 
    ## Run 102 stress 0.09268421 
    ## Run 103 stress 0.09308972 
    ## Run 104 stress 0.1038038 
    ## Run 105 stress 0.09535501 
    ## Run 106 stress 0.3186041 
    ## Run 107 stress 0.09030396 
    ## Run 108 stress 0.09145321 
    ## Run 109 stress 0.08973863 
    ## Run 110 stress 0.08440252 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001197079  max resid 0.0002141481 
    ## ... Similar to previous best
    ## Run 111 stress 0.09416151 
    ## Run 112 stress 0.09374186 
    ## Run 113 stress 0.08973867 
    ## Run 114 stress 0.09145327 
    ## Run 115 stress 0.08503473 
    ## Run 116 stress 0.09465905 
    ## Run 117 stress 0.09590947 
    ## Run 118 stress 0.08773478 
    ## Run 119 stress 0.08973878 
    ## Run 120 stress 0.09030401 
    ## Run 121 stress 0.09030394 
    ## Run 122 stress 0.09416149 
    ## Run 123 stress 0.08503498 
    ## Run 124 stress 0.0996998 
    ## Run 125 stress 0.08773489 
    ## Run 126 stress 0.08503471 
    ## Run 127 stress 0.0877349 
    ## Run 128 stress 0.0940344 
    ## Run 129 stress 0.08773475 
    ## Run 130 stress 0.09407991 
    ## Run 131 stress 0.09030396 
    ## Run 132 stress 0.08773466 
    ## Run 133 stress 0.09145321 
    ## Run 134 stress 0.08973883 
    ## Run 135 stress 0.09416134 
    ## Run 136 stress 0.08440255 
    ## ... Procrustes: rmse 3.034026e-05  max resid 5.725993e-05 
    ## ... Similar to previous best
    ## Run 137 stress 0.1072864 
    ## Run 138 stress 0.09268381 
    ## Run 139 stress 0.08503467 
    ## Run 140 stress 0.08973884 
    ## Run 141 stress 0.09407133 
    ## Run 142 stress 0.09168935 
    ## Run 143 stress 0.1038046 
    ## Run 144 stress 0.09268422 
    ## Run 145 stress 0.09407973 
    ## Run 146 stress 0.09030401 
    ## Run 147 stress 0.09408 
    ## Run 148 stress 0.08973891 
    ## Run 149 stress 0.08440257 
    ## ... Procrustes: rmse 7.230956e-05  max resid 0.0001322316 
    ## ... Similar to previous best
    ## Run 150 stress 0.09286083 
    ## Run 151 stress 0.09308947 
    ## Run 152 stress 0.08773477 
    ## Run 153 stress 0.09447297 
    ## Run 154 stress 0.0940799 
    ## Run 155 stress 0.0940343 
    ## Run 156 stress 0.1038033 
    ## Run 157 stress 0.08503511 
    ## Run 158 stress 0.09374243 
    ## Run 159 stress 0.08503462 
    ## Run 160 stress 0.09308964 
    ## Run 161 stress 0.09417792 
    ## Run 162 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 9.708313e-05  max resid 0.0001956451 
    ## ... Similar to previous best
    ## Run 163 stress 0.09407984 
    ## Run 164 stress 0.09030396 
    ## Run 165 stress 0.09159087 
    ## Run 166 stress 0.08503511 
    ## Run 167 stress 0.08973877 
    ## Run 168 stress 0.09308978 
    ## Run 169 stress 0.1038036 
    ## Run 170 stress 0.08773485 
    ## Run 171 stress 0.1019603 
    ## Run 172 stress 0.09268329 
    ## Run 173 stress 0.09400541 
    ## Run 174 stress 0.0940799 
    ## Run 175 stress 0.09381858 
    ## Run 176 stress 0.08440261 
    ## ... Procrustes: rmse 0.0002135226  max resid 0.000413882 
    ## ... Similar to previous best
    ## Run 177 stress 0.0850348 
    ## Run 178 stress 0.08503567 
    ## Run 179 stress 0.09308969 
    ## Run 180 stress 0.0897388 
    ## Run 181 stress 0.08440264 
    ## ... Procrustes: rmse 0.0001940217  max resid 0.000368827 
    ## ... Similar to previous best
    ## Run 182 stress 0.1052852 
    ## Run 183 stress 0.09337308 
    ## Run 184 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 4.552517e-05  max resid 9.621819e-05 
    ## ... Similar to previous best
    ## Run 185 stress 0.09030405 
    ## Run 186 stress 0.2573973 
    ## Run 187 stress 0.09407977 
    ## Run 188 stress 0.09030397 
    ## Run 189 stress 0.09416136 
    ## Run 190 stress 0.09465915 
    ## Run 191 stress 0.0897388 
    ## Run 192 stress 0.09403432 
    ## Run 193 stress 0.09407988 
    ## Run 194 stress 0.1038047 
    ## Run 195 stress 0.09159085 
    ## Run 196 stress 0.09374323 
    ## Run 197 stress 0.08503465 
    ## Run 198 stress 0.08440251 
    ## ... Procrustes: rmse 3.435069e-05  max resid 7.566648e-05 
    ## ... Similar to previous best
    ## Run 199 stress 0.09159099 
    ## Run 200 stress 0.08503512 
    ## Run 201 stress 0.090304 
    ## Run 202 stress 0.2574213 
    ## Run 203 stress 0.08973878 
    ## Run 204 stress 0.09168947 
    ## Run 205 stress 0.09407976 
    ## Run 206 stress 0.09539201 
    ## Run 207 stress 0.0940053 
    ## Run 208 stress 0.09168953 
    ## Run 209 stress 0.09268348 
    ## Run 210 stress 0.0940053 
    ## Run 211 stress 0.09159098 
    ## Run 212 stress 0.08773475 
    ## Run 213 stress 0.0915909 
    ## Run 214 stress 0.09168932 
    ## Run 215 stress 0.09407998 
    ## Run 216 stress 0.09374275 
    ## Run 217 stress 0.09535473 
    ## Run 218 stress 0.09168958 
    ## Run 219 stress 0.09407969 
    ## Run 220 stress 0.09969967 
    ## Run 221 stress 0.08773468 
    ## Run 222 stress 0.09268332 
    ## Run 223 stress 0.1026215 
    ## Run 224 stress 0.09337275 
    ## Run 225 stress 0.09337229 
    ## Run 226 stress 0.09374268 
    ## Run 227 stress 0.0959082 
    ## Run 228 stress 0.09464481 
    ## Run 229 stress 0.0930897 
    ## Run 230 stress 0.09407986 
    ## Run 231 stress 0.09610743 
    ## Run 232 stress 0.09407975 
    ## Run 233 stress 0.0930897 
    ## Run 234 stress 0.09407967 
    ## Run 235 stress 0.09407997 
    ## Run 236 stress 0.08503503 
    ## Run 237 stress 0.0938059 
    ## Run 238 stress 0.09407984 
    ## Run 239 stress 0.09168947 
    ## Run 240 stress 0.0941614 
    ## Run 241 stress 0.08773466 
    ## Run 242 stress 0.09337284 
    ## Run 243 stress 0.08773481 
    ## Run 244 stress 0.08773475 
    ## Run 245 stress 0.09159084 
    ## Run 246 stress 0.09445476 
    ## Run 247 stress 0.093806 
    ## Run 248 stress 0.09969965 
    ## Run 249 stress 0.09416153 
    ## Run 250 stress 0.09465905 
    ## Run 251 stress 0.08973864 
    ## Run 252 stress 0.09286088 
    ## Run 253 stress 0.09407987 
    ## Run 254 stress 0.09030393 
    ## Run 255 stress 0.09159095 
    ## Run 256 stress 0.09030402 
    ## Run 257 stress 0.09407992 
    ## Run 258 stress 0.09268341 
    ## Run 259 stress 0.09145326 
    ## Run 260 stress 0.09407984 
    ## Run 261 stress 0.09464402 
    ## Run 262 stress 0.08973866 
    ## Run 263 stress 0.09416132 
    ## Run 264 stress 0.08973867 
    ## Run 265 stress 0.09308963 
    ## Run 266 stress 0.08973863 
    ## Run 267 stress 0.0940798 
    ## Run 268 stress 0.09407105 
    ## Run 269 stress 0.08773468 
    ## Run 270 stress 0.0850348 
    ## Run 271 stress 0.08973875 
    ## Run 272 stress 0.09407985 
    ## Run 273 stress 0.09286084 
    ## Run 274 stress 0.09407977 
    ## Run 275 stress 0.09445476 
    ## Run 276 stress 0.08440258 
    ## ... Procrustes: rmse 0.0001378724  max resid 0.0002544098 
    ## ... Similar to previous best
    ## Run 277 stress 0.0877349 
    ## Run 278 stress 0.0941613 
    ## Run 279 stress 0.09030393 
    ## Run 280 stress 0.09168952 
    ## Run 281 stress 0.09286089 
    ## Run 282 stress 0.09030396 
    ## Run 283 stress 0.1038041 
    ## Run 284 stress 0.09337283 
    ## Run 285 stress 0.09030401 
    ## Run 286 stress 0.09337253 
    ## Run 287 stress 0.09380584 
    ## Run 288 stress 0.08440251 
    ## ... Procrustes: rmse 4.365571e-05  max resid 9.300955e-05 
    ## ... Similar to previous best
    ## Run 289 stress 0.09145321 
    ## Run 290 stress 0.09308972 
    ## Run 291 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001509288  max resid 0.0002799796 
    ## ... Similar to previous best
    ## Run 292 stress 0.09030409 
    ## Run 293 stress 0.09308957 
    ## Run 294 stress 0.09159085 
    ## Run 295 stress 0.09760824 
    ## Run 296 stress 0.09381856 
    ## Run 297 stress 0.09721205 
    ## Run 298 stress 0.08503676 
    ## Run 299 stress 0.09168942 
    ## Run 300 stress 0.09168946 
    ## Run 301 stress 0.09381845 
    ## Run 302 stress 0.09321704 
    ## Run 303 stress 0.09321831 
    ## Run 304 stress 0.08440266 
    ## ... Procrustes: rmse 0.0001923824  max resid 0.0003835864 
    ## ... Similar to previous best
    ## Run 305 stress 0.0940712 
    ## Run 306 stress 0.09308973 
    ## Run 307 stress 0.09286093 
    ## Run 308 stress 0.08973881 
    ## Run 309 stress 0.09168961 
    ## Run 310 stress 0.092684 
    ## Run 311 stress 0.0850348 
    ## Run 312 stress 0.08773466 
    ## Run 313 stress 0.08503582 
    ## Run 314 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001191454  max resid 0.0002448238 
    ## ... Similar to previous best
    ## Run 315 stress 0.08973867 
    ## Run 316 stress 0.1052851 
    ## Run 317 stress 0.0946591 
    ## Run 318 stress 0.09030408 
    ## Run 319 stress 0.09337245 
    ## Run 320 stress 0.0946591 
    ## Run 321 stress 0.08440253 
    ## ... Procrustes: rmse 7.734382e-05  max resid 0.0001690926 
    ## ... Similar to previous best
    ## Run 322 stress 0.08973879 
    ## Run 323 stress 0.08973883 
    ## Run 324 stress 0.09337285 
    ## Run 325 stress 0.08973887 
    ## Run 326 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001944901  max resid 0.0003847078 
    ## ... Similar to previous best
    ## Run 327 stress 0.09465904 
    ## Run 328 stress 0.09030393 
    ## Run 329 stress 0.09407998 
    ## Run 330 stress 0.09308957 
    ## Run 331 stress 0.09407986 
    ## Run 332 stress 0.08440254 
    ## ... Procrustes: rmse 8.226952e-05  max resid 0.0001524174 
    ## ... Similar to previous best
    ## Run 333 stress 0.09407969 
    ## Run 334 stress 0.09416142 
    ## Run 335 stress 0.09159086 
    ## Run 336 stress 0.09268352 
    ## Run 337 stress 0.08503456 
    ## Run 338 stress 0.09337241 
    ## Run 339 stress 0.08973868 
    ## Run 340 stress 0.08440262 
    ## ... Procrustes: rmse 0.0001788196  max resid 0.0003269969 
    ## ... Similar to previous best
    ## Run 341 stress 0.08503466 
    ## Run 342 stress 0.08503504 
    ## Run 343 stress 0.09168954 
    ## Run 344 stress 0.09145334 
    ## Run 345 stress 0.09159085 
    ## Run 346 stress 0.0938058 
    ## Run 347 stress 0.09308974 
    ## Run 348 stress 0.09337254 
    ## Run 349 stress 0.09308946 
    ## Run 350 stress 0.09969979 
    ## Run 351 stress 0.08503527 
    ## Run 352 stress 0.09145323 
    ## Run 353 stress 0.09159087 
    ## Run 354 stress 0.1038034 
    ## Run 355 stress 0.08773496 
    ## Run 356 stress 0.09590516 
    ## Run 357 stress 0.09145321 
    ## Run 358 stress 0.08503489 
    ## Run 359 stress 0.09465905 
    ## Run 360 stress 0.09417777 
    ## Run 361 stress 0.08503507 
    ## Run 362 stress 0.09286088 
    ## Run 363 stress 0.09407969 
    ## Run 364 stress 0.09159088 
    ## Run 365 stress 0.09268398 
    ## Run 366 stress 0.08503669 
    ## Run 367 stress 0.09407975 
    ## Run 368 stress 0.08440251 
    ## ... Procrustes: rmse 1.925116e-05  max resid 3.648878e-05 
    ## ... Similar to previous best
    ## Run 369 stress 0.09445488 
    ## Run 370 stress 0.08503563 
    ## Run 371 stress 0.08503513 
    ## Run 372 stress 0.09337275 
    ## Run 373 stress 0.09159085 
    ## Run 374 stress 0.09374264 
    ## Run 375 stress 0.09465905 
    ## Run 376 stress 0.2574686 
    ## Run 377 stress 0.08503502 
    ## Run 378 stress 0.0940801 
    ## Run 379 stress 0.09407457 
    ## Run 380 stress 0.09969988 
    ## Run 381 stress 0.09407977 
    ## Run 382 stress 0.09286083 
    ## Run 383 stress 0.09145326 
    ## Run 384 stress 0.09030393 
    ## Run 385 stress 0.09268336 
    ## Run 386 stress 0.09030403 
    ## Run 387 stress 0.09760821 
    ## Run 388 stress 0.08503609 
    ## Run 389 stress 0.1038039 
    ## Run 390 stress 0.09286084 
    ## Run 391 stress 0.09145335 
    ## Run 392 stress 0.1038037 
    ## Run 393 stress 0.08503462 
    ## Run 394 stress 0.1026215 
    ## Run 395 stress 0.09286083 
    ## Run 396 stress 0.09145334 
    ## Run 397 stress 0.09969989 
    ## Run 398 stress 0.2456114 
    ## Run 399 stress 0.08440262 
    ## ... Procrustes: rmse 0.0001657915  max resid 0.000318133 
    ## ... Similar to previous best
    ## Run 400 stress 0.09407993 
    ## Run 401 stress 0.09030403 
    ## Run 402 stress 0.09030407 
    ## Run 403 stress 0.08440265 
    ## ... Procrustes: rmse 0.000186071  max resid 0.0003520616 
    ## ... Similar to previous best
    ## Run 404 stress 0.09969963 
    ## Run 405 stress 0.08503508 
    ## Run 406 stress 0.2511907 
    ## Run 407 stress 0.09969984 
    ## Run 408 stress 0.09610744 
    ## Run 409 stress 0.08503475 
    ## Run 410 stress 0.08503516 
    ## Run 411 stress 0.09286097 
    ## Run 412 stress 0.09908186 
    ## Run 413 stress 0.09308947 
    ## Run 414 stress 0.09337272 
    ## Run 415 stress 0.08773483 
    ## Run 416 stress 0.09268363 
    ## Run 417 stress 0.09669869 
    ## Run 418 stress 0.09407105 
    ## Run 419 stress 0.09145332 
    ## Run 420 stress 0.09308975 
    ## Run 421 stress 0.09268404 
    ## Run 422 stress 0.09030401 
    ## Run 423 stress 0.08773473 
    ## Run 424 stress 0.09321916 
    ## Run 425 stress 0.09969994 
    ## Run 426 stress 0.08503572 
    ## Run 427 stress 0.3178306 
    ## Run 428 stress 0.08440258 
    ## ... Procrustes: rmse 0.0001405435  max resid 0.0002619542 
    ## ... Similar to previous best
    ## Run 429 stress 0.09465906 
    ## Run 430 stress 0.09417798 
    ## Run 431 stress 0.09168947 
    ## Run 432 stress 0.09400561 
    ## Run 433 stress 0.0930896 
    ## Run 434 stress 0.09416142 
    ## Run 435 stress 0.09168933 
    ## Run 436 stress 0.08773476 
    ## Run 437 stress 0.09403406 
    ## Run 438 stress 0.09159084 
    ## Run 439 stress 0.09969971 
    ## Run 440 stress 0.09321612 
    ## Run 441 stress 0.08503601 
    ## Run 442 stress 0.08503567 
    ## Run 443 stress 0.0850361 
    ## Run 444 stress 0.08440271 
    ## ... Procrustes: rmse 0.0002451742  max resid 0.0004435302 
    ## ... Similar to previous best
    ## Run 445 stress 0.08440276 
    ## ... Procrustes: rmse 0.0003207772  max resid 0.0006302176 
    ## ... Similar to previous best
    ## Run 446 stress 0.08503668 
    ## Run 447 stress 0.08973878 
    ## Run 448 stress 0.0940797 
    ## Run 449 stress 0.0946591 
    ## Run 450 stress 0.09465909 
    ## Run 451 stress 0.09380576 
    ## Run 452 stress 0.0940798 
    ## Run 453 stress 0.090304 
    ## Run 454 stress 0.08773481 
    ## Run 455 stress 0.09308964 
    ## Run 456 stress 0.09535514 
    ## Run 457 stress 0.09159083 
    ## Run 458 stress 0.08440253 
    ## ... Procrustes: rmse 7.700429e-05  max resid 0.0001547187 
    ## ... Similar to previous best
    ## Run 459 stress 0.08973898 
    ## Run 460 stress 0.09159084 
    ## Run 461 stress 0.08973864 
    ## Run 462 stress 0.08973893 
    ## Run 463 stress 0.09416137 
    ## Run 464 stress 0.08503478 
    ## Run 465 stress 0.09268378 
    ## Run 466 stress 0.09407994 
    ## Run 467 stress 0.09374256 
    ## Run 468 stress 0.08503501 
    ## Run 469 stress 0.3198888 
    ## Run 470 stress 0.0940739 
    ## Run 471 stress 0.09590912 
    ## Run 472 stress 0.0937419 
    ## Run 473 stress 0.09416142 
    ## Run 474 stress 0.09168959 
    ## Run 475 stress 0.09168947 
    ## Run 476 stress 0.09374222 
    ## Run 477 stress 0.09337298 
    ## Run 478 stress 0.09590683 
    ## Run 479 stress 0.309074 
    ## Run 480 stress 0.08440271 
    ## ... Procrustes: rmse 0.0002501736  max resid 0.0004511562 
    ## ... Similar to previous best
    ## Run 481 stress 0.0903041 
    ## Run 482 stress 0.09969969 
    ## Run 483 stress 0.0940798 
    ## Run 484 stress 0.09535438 
    ## Run 485 stress 0.08503656 
    ## Run 486 stress 0.09030411 
    ## Run 487 stress 0.09977534 
    ## Run 488 stress 0.08773466 
    ## Run 489 stress 0.09145325 
    ## Run 490 stress 0.08503539 
    ## Run 491 stress 0.09969995 
    ## Run 492 stress 0.08440273 
    ## ... Procrustes: rmse 0.000252804  max resid 0.0004616358 
    ## ... Similar to previous best
    ## Run 493 stress 0.08440272 
    ## ... Procrustes: rmse 0.0002415455  max resid 0.000413082 
    ## ... Similar to previous best
    ## Run 494 stress 0.09416131 
    ## Run 495 stress 0.08503511 
    ## Run 496 stress 0.08773472 
    ## Run 497 stress 0.09159083 
    ## Run 498 stress 0.08773465 
    ## Run 499 stress 0.08503479 
    ## Run 500 stress 0.08503462 
    ## *** Best solution repeated 21 times

``` r
# Ocean sites and mixed lakes
SD_beta_OM_NMDS <- metaMDS(SD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07365783 
    ## Run 1 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743652  max resid 0.05101221 
    ## Run 2 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001145876  max resid 0.0002699915 
    ## ... Similar to previous best
    ## Run 3 stress 0.07365786 
    ## ... Procrustes: rmse 6.756139e-05  max resid 0.0001616859 
    ## ... Similar to previous best
    ## Run 4 stress 0.07732931 
    ## Run 5 stress 0.07732925 
    ## Run 6 stress 0.3495557 
    ## Run 7 stress 0.07378232 
    ## ... Procrustes: rmse 0.01746284  max resid 0.05114882 
    ## Run 8 stress 0.07365783 
    ## ... Procrustes: rmse 2.419051e-05  max resid 5.734775e-05 
    ## ... Similar to previous best
    ## Run 9 stress 0.07365783 
    ## ... Procrustes: rmse 1.464699e-05  max resid 3.235132e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.08000467 
    ## Run 11 stress 0.07629235 
    ## Run 12 stress 0.07629239 
    ## Run 13 stress 0.07732926 
    ## Run 14 stress 0.07629242 
    ## Run 15 stress 0.08000467 
    ## Run 16 stress 0.07365785 
    ## ... Procrustes: rmse 8.150716e-05  max resid 0.000194126 
    ## ... Similar to previous best
    ## Run 17 stress 0.07378228 
    ## ... Procrustes: rmse 0.01743948  max resid 0.05100027 
    ## Run 18 stress 0.08000473 
    ## Run 19 stress 0.08000474 
    ## Run 20 stress 0.07365787 
    ## ... Procrustes: rmse 9.345249e-05  max resid 0.0002165002 
    ## ... Similar to previous best
    ## Run 21 stress 0.07629233 
    ## Run 22 stress 0.07365786 
    ## ... Procrustes: rmse 7.671448e-05  max resid 0.0001789183 
    ## ... Similar to previous best
    ## Run 23 stress 0.0762924 
    ## Run 24 stress 0.07732927 
    ## Run 25 stress 0.08000471 
    ## Run 26 stress 0.07629234 
    ## Run 27 stress 0.07629239 
    ## Run 28 stress 0.07365784 
    ## ... Procrustes: rmse 5.329867e-05  max resid 0.000123712 
    ## ... Similar to previous best
    ## Run 29 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745536  max resid 0.05108356 
    ## Run 30 stress 0.08000468 
    ## Run 31 stress 0.07365784 
    ## ... Procrustes: rmse 4.870508e-05  max resid 0.0001149839 
    ## ... Similar to previous best
    ## Run 32 stress 0.0737823 
    ## ... Procrustes: rmse 0.01747049  max resid 0.05117443 
    ## Run 33 stress 0.07629244 
    ## Run 34 stress 0.07629233 
    ## Run 35 stress 0.07629252 
    ## Run 36 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174441  max resid 0.05104209 
    ## Run 37 stress 0.07365784 
    ## ... Procrustes: rmse 4.743888e-05  max resid 0.0001098647 
    ## ... Similar to previous best
    ## Run 38 stress 0.08000472 
    ## Run 39 stress 0.07629234 
    ## Run 40 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001656469  max resid 0.0003897693 
    ## ... Similar to previous best
    ## Run 41 stress 0.0762924 
    ## Run 42 stress 0.07629238 
    ## Run 43 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745197  max resid 0.05106697 
    ## Run 44 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744194  max resid 0.0510285 
    ## Run 45 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001209524  max resid 0.0002889866 
    ## ... Similar to previous best
    ## Run 46 stress 0.07732926 
    ## Run 47 stress 0.07629236 
    ## Run 48 stress 0.07629235 
    ## Run 49 stress 0.07732932 
    ## Run 50 stress 0.07365785 
    ## ... Procrustes: rmse 6.898509e-05  max resid 0.0001639082 
    ## ... Similar to previous best
    ## Run 51 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744278  max resid 0.05102261 
    ## Run 52 stress 0.07365784 
    ## ... Procrustes: rmse 4.022724e-05  max resid 9.580554e-05 
    ## ... Similar to previous best
    ## Run 53 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745532  max resid 0.05109929 
    ## Run 54 stress 0.07365783 
    ## ... Procrustes: rmse 2.059119e-05  max resid 4.851631e-05 
    ## ... Similar to previous best
    ## Run 55 stress 0.07365784 
    ## ... Procrustes: rmse 6.072646e-05  max resid 0.0001430002 
    ## ... Similar to previous best
    ## Run 56 stress 0.0737823 
    ## ... Procrustes: rmse 0.01743904  max resid 0.05104797 
    ## Run 57 stress 0.07378231 
    ## ... Procrustes: rmse 0.01744456  max resid 0.0510117 
    ## Run 58 stress 0.08000467 
    ## Run 59 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745005  max resid 0.05105769 
    ## Run 60 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 1.248161e-05  max resid 2.875832e-05 
    ## ... Similar to previous best
    ## Run 61 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174416  max resid 0.05102279 
    ## Run 62 stress 0.3097962 
    ## Run 63 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744162  max resid 0.0510568 
    ## Run 64 stress 0.07365792 
    ## ... Procrustes: rmse 0.000162715  max resid 0.0003854743 
    ## ... Similar to previous best
    ## Run 65 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001630934  max resid 0.0003845393 
    ## ... Similar to previous best
    ## Run 66 stress 0.07365791 
    ## ... Procrustes: rmse 0.000140388  max resid 0.0003355135 
    ## ... Similar to previous best
    ## Run 67 stress 0.08000472 
    ## Run 68 stress 0.07629244 
    ## Run 69 stress 0.07629234 
    ## Run 70 stress 0.07732925 
    ## Run 71 stress 0.07732931 
    ## Run 72 stress 0.07629239 
    ## Run 73 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744248  max resid 0.05102096 
    ## Run 74 stress 0.08233759 
    ## Run 75 stress 0.07365783 
    ## ... Procrustes: rmse 1.89592e-05  max resid 4.663835e-05 
    ## ... Similar to previous best
    ## Run 76 stress 0.07365783 
    ## ... Procrustes: rmse 1.61782e-05  max resid 3.793023e-05 
    ## ... Similar to previous best
    ## Run 77 stress 0.08233758 
    ## Run 78 stress 0.07629243 
    ## Run 79 stress 0.08288044 
    ## Run 80 stress 0.0773293 
    ## Run 81 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174475  max resid 0.05104337 
    ## Run 82 stress 0.07629237 
    ## Run 83 stress 0.07732944 
    ## Run 84 stress 0.08288052 
    ## Run 85 stress 0.08233758 
    ## Run 86 stress 0.08000469 
    ## Run 87 stress 0.07629236 
    ## Run 88 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744485  max resid 0.05103868 
    ## Run 89 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744343  max resid 0.05103932 
    ## Run 90 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744859  max resid 0.05106585 
    ## Run 91 stress 0.0823376 
    ## Run 92 stress 0.08000469 
    ## Run 93 stress 0.07629236 
    ## Run 94 stress 0.07629246 
    ## Run 95 stress 0.07629234 
    ## Run 96 stress 0.07629235 
    ## Run 97 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001360771  max resid 0.0003250775 
    ## ... Similar to previous best
    ## Run 98 stress 0.07629243 
    ## Run 99 stress 0.07629241 
    ## Run 100 stress 0.08000468 
    ## Run 101 stress 0.07365784 
    ## ... Procrustes: rmse 4.496595e-05  max resid 0.0001052816 
    ## ... Similar to previous best
    ## Run 102 stress 0.07732928 
    ## Run 103 stress 0.08233753 
    ## Run 104 stress 0.08000467 
    ## Run 105 stress 0.07365786 
    ## ... Procrustes: rmse 9.460288e-05  max resid 0.0002211075 
    ## ... Similar to previous best
    ## Run 106 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744566  max resid 0.05106214 
    ## Run 107 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744281  max resid 0.05101435 
    ## Run 108 stress 0.07365783 
    ## ... Procrustes: rmse 3.496777e-05  max resid 8.176406e-05 
    ## ... Similar to previous best
    ## Run 109 stress 0.07629237 
    ## Run 110 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744908  max resid 0.05107279 
    ## Run 111 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744801  max resid 0.05106931 
    ## Run 112 stress 0.08000485 
    ## Run 113 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 6.162905e-06  max resid 1.319751e-05 
    ## ... Similar to previous best
    ## Run 114 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744836  max resid 0.05103788 
    ## Run 115 stress 0.07365783 
    ## ... Procrustes: rmse 9.192189e-06  max resid 1.963196e-05 
    ## ... Similar to previous best
    ## Run 116 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744564  max resid 0.05105816 
    ## Run 117 stress 0.08000475 
    ## Run 118 stress 0.0823377 
    ## Run 119 stress 0.07365784 
    ## ... Procrustes: rmse 4.497347e-05  max resid 0.0001062784 
    ## ... Similar to previous best
    ## Run 120 stress 0.07365783 
    ## ... Procrustes: rmse 3.365649e-05  max resid 7.893323e-05 
    ## ... Similar to previous best
    ## Run 121 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744605  max resid 0.05105322 
    ## Run 122 stress 0.07732948 
    ## Run 123 stress 0.07365788 
    ## ... Procrustes: rmse 5.918329e-05  max resid 0.0001187954 
    ## ... Similar to previous best
    ## Run 124 stress 0.08288055 
    ## Run 125 stress 0.07629233 
    ## Run 126 stress 0.07629243 
    ## Run 127 stress 0.08000477 
    ## Run 128 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001299809  max resid 0.0003041402 
    ## ... Similar to previous best
    ## Run 129 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001135459  max resid 0.0002683627 
    ## ... Similar to previous best
    ## Run 130 stress 0.07378231 
    ## ... Procrustes: rmse 0.0174434  max resid 0.05099868 
    ## Run 131 stress 0.07378232 
    ## ... Procrustes: rmse 0.01743332  max resid 0.05095861 
    ## Run 132 stress 0.07629234 
    ## Run 133 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744937  max resid 0.05104269 
    ## Run 134 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001112903  max resid 0.0002539344 
    ## ... Similar to previous best
    ## Run 135 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744795  max resid 0.0510518 
    ## Run 136 stress 0.07629239 
    ## Run 137 stress 0.07378227 
    ## ... Procrustes: rmse 0.01742742  max resid 0.05096964 
    ## Run 138 stress 0.07629238 
    ## Run 139 stress 0.07629234 
    ## Run 140 stress 0.07629237 
    ## Run 141 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001555131  max resid 0.0003654147 
    ## ... Similar to previous best
    ## Run 142 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744343  max resid 0.05102208 
    ## Run 143 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001210497  max resid 0.0002852945 
    ## ... Similar to previous best
    ## Run 144 stress 0.07365784 
    ## ... Procrustes: rmse 4.146557e-05  max resid 9.878609e-05 
    ## ... Similar to previous best
    ## Run 145 stress 0.07365785 
    ## ... Procrustes: rmse 6.698548e-05  max resid 0.0001564063 
    ## ... Similar to previous best
    ## Run 146 stress 0.07732942 
    ## Run 147 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743676  max resid 0.05100814 
    ## Run 148 stress 0.07365784 
    ## ... Procrustes: rmse 3.575052e-05  max resid 7.735919e-05 
    ## ... Similar to previous best
    ## Run 149 stress 0.07365785 
    ## ... Procrustes: rmse 6.736129e-05  max resid 0.0001574323 
    ## ... Similar to previous best
    ## Run 150 stress 0.08000472 
    ## Run 151 stress 0.07378227 
    ## ... Procrustes: rmse 0.01742135  max resid 0.05093606 
    ## Run 152 stress 0.08000472 
    ## Run 153 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001367238  max resid 0.0003187801 
    ## ... Similar to previous best
    ## Run 154 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001477535  max resid 0.0003471194 
    ## ... Similar to previous best
    ## Run 155 stress 0.07629235 
    ## Run 156 stress 0.07629234 
    ## Run 157 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001181495  max resid 0.0002796368 
    ## ... Similar to previous best
    ## Run 158 stress 0.07629245 
    ## Run 159 stress 0.07378233 
    ## ... Procrustes: rmse 0.01744934  max resid 0.05103031 
    ## Run 160 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744459  max resid 0.05102641 
    ## Run 161 stress 0.07732939 
    ## Run 162 stress 0.07365784 
    ## ... Procrustes: rmse 6.300052e-05  max resid 0.0001491374 
    ## ... Similar to previous best
    ## Run 163 stress 0.07629234 
    ## Run 164 stress 0.3036923 
    ## Run 165 stress 0.08000469 
    ## Run 166 stress 0.07629239 
    ## Run 167 stress 0.08288034 
    ## Run 168 stress 0.08000472 
    ## Run 169 stress 0.07732934 
    ## Run 170 stress 0.07365784 
    ## ... Procrustes: rmse 3.745855e-05  max resid 8.541974e-05 
    ## ... Similar to previous best
    ## Run 171 stress 0.07365783 
    ## ... Procrustes: rmse 8.436653e-06  max resid 1.51422e-05 
    ## ... Similar to previous best
    ## Run 172 stress 0.07365784 
    ## ... Procrustes: rmse 3.64983e-05  max resid 8.789025e-05 
    ## ... Similar to previous best
    ## Run 173 stress 0.07732933 
    ## Run 174 stress 0.07365784 
    ## ... Procrustes: rmse 3.624997e-05  max resid 8.587269e-05 
    ## ... Similar to previous best
    ## Run 175 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001346711  max resid 0.0003186657 
    ## ... Similar to previous best
    ## Run 176 stress 0.07629239 
    ## Run 177 stress 0.07365783 
    ## ... Procrustes: rmse 1.444679e-05  max resid 3.343075e-05 
    ## ... Similar to previous best
    ## Run 178 stress 0.07365784 
    ## ... Procrustes: rmse 6.343275e-05  max resid 0.0001500743 
    ## ... Similar to previous best
    ## Run 179 stress 0.2574761 
    ## Run 180 stress 0.07629233 
    ## Run 181 stress 0.07629234 
    ## Run 182 stress 0.3407914 
    ## Run 183 stress 0.07629246 
    ## Run 184 stress 0.07365783 
    ## ... Procrustes: rmse 1.811823e-05  max resid 4.298472e-05 
    ## ... Similar to previous best
    ## Run 185 stress 0.08233765 
    ## Run 186 stress 0.07365783 
    ## ... Procrustes: rmse 1.20334e-05  max resid 2.158426e-05 
    ## ... Similar to previous best
    ## Run 187 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744384  max resid 0.05101767 
    ## Run 188 stress 0.07732936 
    ## Run 189 stress 0.08000476 
    ## Run 190 stress 0.07732943 
    ## Run 191 stress 0.0773294 
    ## Run 192 stress 0.07365783 
    ## ... Procrustes: rmse 4.099015e-06  max resid 8.641282e-06 
    ## ... Similar to previous best
    ## Run 193 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744791  max resid 0.05104749 
    ## Run 194 stress 0.07365784 
    ## ... Procrustes: rmse 4.318434e-05  max resid 0.0001027151 
    ## ... Similar to previous best
    ## Run 195 stress 0.07365784 
    ## ... Procrustes: rmse 4.654598e-05  max resid 0.0001074245 
    ## ... Similar to previous best
    ## Run 196 stress 0.0737823 
    ## ... Procrustes: rmse 0.01746248  max resid 0.05113753 
    ## Run 197 stress 0.08000468 
    ## Run 198 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744665  max resid 0.0510513 
    ## Run 199 stress 0.07378228 
    ## ... Procrustes: rmse 0.01743063  max resid 0.05095907 
    ## Run 200 stress 0.07629236 
    ## Run 201 stress 0.08288048 
    ## Run 202 stress 0.07365784 
    ## ... Procrustes: rmse 3.850049e-05  max resid 8.999335e-05 
    ## ... Similar to previous best
    ## Run 203 stress 0.07365783 
    ## ... Procrustes: rmse 1.499931e-05  max resid 3.066611e-05 
    ## ... Similar to previous best
    ## Run 204 stress 0.07629242 
    ## Run 205 stress 0.08000467 
    ## Run 206 stress 0.07365784 
    ## ... Procrustes: rmse 3.170628e-05  max resid 7.191258e-05 
    ## ... Similar to previous best
    ## Run 207 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744329  max resid 0.05102593 
    ## Run 208 stress 0.07629235 
    ## Run 209 stress 0.08288046 
    ## Run 210 stress 0.07629243 
    ## Run 211 stress 0.08000468 
    ## Run 212 stress 0.08288048 
    ## Run 213 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744679  max resid 0.05105355 
    ## Run 214 stress 0.07629241 
    ## Run 215 stress 0.07365783 
    ## ... Procrustes: rmse 2.738987e-05  max resid 6.525943e-05 
    ## ... Similar to previous best
    ## Run 216 stress 0.0762924 
    ## Run 217 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745027  max resid 0.05104591 
    ## Run 218 stress 0.07365784 
    ## ... Procrustes: rmse 5.643965e-05  max resid 0.0001358145 
    ## ... Similar to previous best
    ## Run 219 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001315098  max resid 0.0003087572 
    ## ... Similar to previous best
    ## Run 220 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001040592  max resid 0.0002432956 
    ## ... Similar to previous best
    ## Run 221 stress 0.07365783 
    ## ... Procrustes: rmse 2.292895e-06  max resid 4.079156e-06 
    ## ... Similar to previous best
    ## Run 222 stress 0.07629234 
    ## Run 223 stress 0.07629245 
    ## Run 224 stress 0.07378231 
    ## ... Procrustes: rmse 0.01745953  max resid 0.05112432 
    ## Run 225 stress 0.07365783 
    ## ... Procrustes: rmse 4.15036e-06  max resid 7.18169e-06 
    ## ... Similar to previous best
    ## Run 226 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744791  max resid 0.05104943 
    ## Run 227 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744841  max resid 0.05105336 
    ## Run 228 stress 0.0823376 
    ## Run 229 stress 0.07629239 
    ## Run 230 stress 0.07629236 
    ## Run 231 stress 0.07365786 
    ## ... Procrustes: rmse 8.722331e-05  max resid 0.0002068863 
    ## ... Similar to previous best
    ## Run 232 stress 0.07629243 
    ## Run 233 stress 0.07629246 
    ## Run 234 stress 0.07629239 
    ## Run 235 stress 0.07365786 
    ## ... Procrustes: rmse 9.212586e-05  max resid 0.0002174219 
    ## ... Similar to previous best
    ## Run 236 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744981  max resid 0.05107144 
    ## Run 237 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745146  max resid 0.05105366 
    ## Run 238 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001061974  max resid 0.0002503329 
    ## ... Similar to previous best
    ## Run 239 stress 0.07732926 
    ## Run 240 stress 0.07365784 
    ## ... Procrustes: rmse 5.924461e-05  max resid 0.0001367758 
    ## ... Similar to previous best
    ## Run 241 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744843  max resid 0.05106123 
    ## Run 242 stress 0.07378234 
    ## ... Procrustes: rmse 0.01743809  max resid 0.05096492 
    ## Run 243 stress 0.07629233 
    ## Run 244 stress 0.07365783 
    ## ... Procrustes: rmse 2.523912e-05  max resid 6.143135e-05 
    ## ... Similar to previous best
    ## Run 245 stress 0.07629241 
    ## Run 246 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001116302  max resid 0.0002635658 
    ## ... Similar to previous best
    ## Run 247 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744493  max resid 0.05106291 
    ## Run 248 stress 0.07629236 
    ## Run 249 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744568  max resid 0.05103959 
    ## Run 250 stress 0.07365784 
    ## ... Procrustes: rmse 4.610709e-05  max resid 0.000109044 
    ## ... Similar to previous best
    ## Run 251 stress 0.07732926 
    ## Run 252 stress 0.08233762 
    ## Run 253 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 2.420304e-06  max resid 4.368735e-06 
    ## ... Similar to previous best
    ## Run 254 stress 0.08000477 
    ## Run 255 stress 0.07629239 
    ## Run 256 stress 0.07629249 
    ## Run 257 stress 0.08000469 
    ## Run 258 stress 0.2901992 
    ## Run 259 stress 0.07629244 
    ## Run 260 stress 0.07629237 
    ## Run 261 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744803  max resid 0.05106153 
    ## Run 262 stress 0.07732933 
    ## Run 263 stress 0.07732932 
    ## Run 264 stress 0.07629235 
    ## Run 265 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744636  max resid 0.05105237 
    ## Run 266 stress 0.07378232 
    ## ... Procrustes: rmse 0.01743894  max resid 0.05105122 
    ## Run 267 stress 0.07365785 
    ## ... Procrustes: rmse 7.123107e-05  max resid 0.0001681343 
    ## ... Similar to previous best
    ## Run 268 stress 0.07365785 
    ## ... Procrustes: rmse 6.441122e-05  max resid 0.0001516688 
    ## ... Similar to previous best
    ## Run 269 stress 0.0828804 
    ## Run 270 stress 0.08000468 
    ## Run 271 stress 0.08000469 
    ## Run 272 stress 0.07365784 
    ## ... Procrustes: rmse 2.118492e-05  max resid 4.772039e-05 
    ## ... Similar to previous best
    ## Run 273 stress 0.3416215 
    ## Run 274 stress 0.08000468 
    ## Run 275 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744495  max resid 0.05103943 
    ## Run 276 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745602  max resid 0.05110557 
    ## Run 277 stress 0.07629234 
    ## Run 278 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744485  max resid 0.0510244 
    ## Run 279 stress 0.07629244 
    ## Run 280 stress 0.07629242 
    ## Run 281 stress 0.07629233 
    ## Run 282 stress 0.08233763 
    ## Run 283 stress 0.07629234 
    ## Run 284 stress 0.07732935 
    ## Run 285 stress 0.07365784 
    ## ... Procrustes: rmse 3.404597e-05  max resid 8.132954e-05 
    ## ... Similar to previous best
    ## Run 286 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744592  max resid 0.05104681 
    ## Run 287 stress 0.08000467 
    ## Run 288 stress 0.08000467 
    ## Run 289 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001404513  max resid 0.00033567 
    ## ... Similar to previous best
    ## Run 290 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745645  max resid 0.05108798 
    ## Run 291 stress 0.07365783 
    ## ... Procrustes: rmse 2.243066e-05  max resid 5.265191e-05 
    ## ... Similar to previous best
    ## Run 292 stress 0.07365784 
    ## ... Procrustes: rmse 3.391016e-05  max resid 7.829938e-05 
    ## ... Similar to previous best
    ## Run 293 stress 0.07365789 
    ## ... Procrustes: rmse 0.000121309  max resid 0.000282444 
    ## ... Similar to previous best
    ## Run 294 stress 0.07365784 
    ## ... Procrustes: rmse 6.232419e-05  max resid 0.0001460051 
    ## ... Similar to previous best
    ## Run 295 stress 0.07629244 
    ## Run 296 stress 0.08233752 
    ## Run 297 stress 0.2584051 
    ## Run 298 stress 0.07365783 
    ## ... Procrustes: rmse 2.376738e-05  max resid 5.360976e-05 
    ## ... Similar to previous best
    ## Run 299 stress 0.07365783 
    ## ... Procrustes: rmse 3.025242e-05  max resid 6.944577e-05 
    ## ... Similar to previous best
    ## Run 300 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744681  max resid 0.05106143 
    ## Run 301 stress 0.08288057 
    ## Run 302 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745159  max resid 0.05106488 
    ## Run 303 stress 0.07378227 
    ## ... Procrustes: rmse 0.01746232  max resid 0.05111576 
    ## Run 304 stress 0.07629235 
    ## Run 305 stress 0.07629245 
    ## Run 306 stress 0.07365783 
    ## ... Procrustes: rmse 1.740012e-06  max resid 3.878521e-06 
    ## ... Similar to previous best
    ## Run 307 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744708  max resid 0.05105545 
    ## Run 308 stress 0.07365783 
    ## ... Procrustes: rmse 4.957015e-06  max resid 9.914937e-06 
    ## ... Similar to previous best
    ## Run 309 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744596  max resid 0.05104483 
    ## Run 310 stress 0.08000467 
    ## Run 311 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001091052  max resid 0.0002578006 
    ## ... Similar to previous best
    ## Run 312 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743606  max resid 0.05100704 
    ## Run 313 stress 0.07732939 
    ## Run 314 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744422  max resid 0.05105839 
    ## Run 315 stress 0.08233755 
    ## Run 316 stress 0.07378234 
    ## ... Procrustes: rmse 0.01747006  max resid 0.05118247 
    ## Run 317 stress 0.3495478 
    ## Run 318 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744562  max resid 0.05105801 
    ## Run 319 stress 0.07629236 
    ## Run 320 stress 0.0828804 
    ## Run 321 stress 0.0800047 
    ## Run 322 stress 0.07629234 
    ## Run 323 stress 0.08000475 
    ## Run 324 stress 0.08233753 
    ## Run 325 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744625  max resid 0.05106508 
    ## Run 326 stress 0.08000467 
    ## Run 327 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744488  max resid 0.05104602 
    ## Run 328 stress 0.07365786 
    ## ... Procrustes: rmse 3.462797e-05  max resid 6.29198e-05 
    ## ... Similar to previous best
    ## Run 329 stress 0.08233754 
    ## Run 330 stress 0.07365785 
    ## ... Procrustes: rmse 6.234563e-05  max resid 0.0001490423 
    ## ... Similar to previous best
    ## Run 331 stress 0.08000471 
    ## Run 332 stress 0.08000473 
    ## Run 333 stress 0.07629239 
    ## Run 334 stress 0.07629236 
    ## Run 335 stress 0.08000468 
    ## Run 336 stress 0.07365783 
    ## ... Procrustes: rmse 2.828472e-05  max resid 6.693658e-05 
    ## ... Similar to previous best
    ## Run 337 stress 0.08000471 
    ## Run 338 stress 0.07365783 
    ## ... Procrustes: rmse 1.703915e-05  max resid 3.43229e-05 
    ## ... Similar to previous best
    ## Run 339 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744272  max resid 0.05101369 
    ## Run 340 stress 0.07378234 
    ## ... Procrustes: rmse 0.01743469  max resid 0.05104206 
    ## Run 341 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745657  max resid 0.05108957 
    ## Run 342 stress 0.07365785 
    ## ... Procrustes: rmse 6.857751e-05  max resid 0.0001593145 
    ## ... Similar to previous best
    ## Run 343 stress 0.07365785 
    ## ... Procrustes: rmse 7.376529e-05  max resid 0.0001741287 
    ## ... Similar to previous best
    ## Run 344 stress 0.07732931 
    ## Run 345 stress 0.07629242 
    ## Run 346 stress 0.0762925 
    ## Run 347 stress 0.07365783 
    ## ... Procrustes: rmse 6.018298e-06  max resid 1.267659e-05 
    ## ... Similar to previous best
    ## Run 348 stress 0.08288056 
    ## Run 349 stress 0.08000469 
    ## Run 350 stress 0.07365785 
    ## ... Procrustes: rmse 6.972633e-05  max resid 0.0001654482 
    ## ... Similar to previous best
    ## Run 351 stress 0.07732936 
    ## Run 352 stress 0.0800047 
    ## Run 353 stress 0.07629243 
    ## Run 354 stress 0.07365784 
    ## ... Procrustes: rmse 1.864548e-05  max resid 2.769305e-05 
    ## ... Similar to previous best
    ## Run 355 stress 0.0823376 
    ## Run 356 stress 0.07629236 
    ## Run 357 stress 0.07629247 
    ## Run 358 stress 0.08233759 
    ## Run 359 stress 0.07365783 
    ## ... Procrustes: rmse 2.1061e-05  max resid 4.939558e-05 
    ## ... Similar to previous best
    ## Run 360 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001022675  max resid 0.0002392158 
    ## ... Similar to previous best
    ## Run 361 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744538  max resid 0.05103905 
    ## Run 362 stress 0.08000468 
    ## Run 363 stress 0.08233752 
    ## Run 364 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744434  max resid 0.05103307 
    ## Run 365 stress 0.07629238 
    ## Run 366 stress 0.07629236 
    ## Run 367 stress 0.07365784 
    ## ... Procrustes: rmse 3.053094e-05  max resid 6.80556e-05 
    ## ... Similar to previous best
    ## Run 368 stress 0.07365783 
    ## ... Procrustes: rmse 7.78512e-06  max resid 1.825027e-05 
    ## ... Similar to previous best
    ## Run 369 stress 0.07365783 
    ## ... Procrustes: rmse 1.244136e-05  max resid 2.823853e-05 
    ## ... Similar to previous best
    ## Run 370 stress 0.0773294 
    ## Run 371 stress 0.08233772 
    ## Run 372 stress 0.07732926 
    ## Run 373 stress 0.07378227 
    ## ... Procrustes: rmse 0.01742534  max resid 0.05096361 
    ## Run 374 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744149  max resid 0.05105155 
    ## Run 375 stress 0.08233757 
    ## Run 376 stress 0.08000469 
    ## Run 377 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744989  max resid 0.05105323 
    ## Run 378 stress 0.08233759 
    ## Run 379 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744537  max resid 0.05103789 
    ## Run 380 stress 0.08000477 
    ## Run 381 stress 0.07732925 
    ## Run 382 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001636422  max resid 0.000381257 
    ## ... Similar to previous best
    ## Run 383 stress 0.07365784 
    ## ... Procrustes: rmse 5.165175e-05  max resid 0.0001225693 
    ## ... Similar to previous best
    ## Run 384 stress 0.07365785 
    ## ... Procrustes: rmse 7.185669e-05  max resid 0.0001686949 
    ## ... Similar to previous best
    ## Run 385 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001500991  max resid 0.0003512107 
    ## ... Similar to previous best
    ## Run 386 stress 0.07365784 
    ## ... Procrustes: rmse 3.672454e-05  max resid 8.618607e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.08000467 
    ## Run 388 stress 0.08000468 
    ## Run 389 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744524  max resid 0.05103465 
    ## Run 390 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001465998  max resid 0.0003451991 
    ## ... Similar to previous best
    ## Run 391 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745614  max resid 0.05110103 
    ## Run 392 stress 0.07629238 
    ## Run 393 stress 0.07365785 
    ## ... Procrustes: rmse 8.083485e-05  max resid 0.0001920075 
    ## ... Similar to previous best
    ## Run 394 stress 0.07629241 
    ## Run 395 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174438  max resid 0.05103244 
    ## Run 396 stress 0.07365784 
    ## ... Procrustes: rmse 5.378527e-05  max resid 0.0001266416 
    ## ... Similar to previous best
    ## Run 397 stress 0.08288065 
    ## Run 398 stress 0.08233755 
    ## Run 399 stress 0.08000467 
    ## Run 400 stress 0.07629233 
    ## Run 401 stress 0.07365783 
    ## ... Procrustes: rmse 7.041027e-06  max resid 1.654254e-05 
    ## ... Similar to previous best
    ## Run 402 stress 0.07365783 
    ## ... Procrustes: rmse 1.653494e-05  max resid 3.947652e-05 
    ## ... Similar to previous best
    ## Run 403 stress 0.07732937 
    ## Run 404 stress 0.0823376 
    ## Run 405 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744962  max resid 0.05105963 
    ## Run 406 stress 0.07732932 
    ## Run 407 stress 0.07365786 
    ## ... Procrustes: rmse 8.707754e-05  max resid 0.0002071476 
    ## ... Similar to previous best
    ## Run 408 stress 0.08233754 
    ## Run 409 stress 0.07629233 
    ## Run 410 stress 0.08233775 
    ## Run 411 stress 0.07378231 
    ## ... Procrustes: rmse 0.01746245  max resid 0.05114023 
    ## Run 412 stress 0.08000467 
    ## Run 413 stress 0.07629238 
    ## Run 414 stress 0.07365783 
    ## ... Procrustes: rmse 2.03458e-06  max resid 3.819842e-06 
    ## ... Similar to previous best
    ## Run 415 stress 0.07365785 
    ## ... Procrustes: rmse 7.176659e-05  max resid 0.0001685789 
    ## ... Similar to previous best
    ## Run 416 stress 0.07629234 
    ## Run 417 stress 0.08233757 
    ## Run 418 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745397  max resid 0.05105652 
    ## Run 419 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744801  max resid 0.05106116 
    ## Run 420 stress 0.07629244 
    ## Run 421 stress 0.08000478 
    ## Run 422 stress 0.07629237 
    ## Run 423 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744619  max resid 0.05105514 
    ## Run 424 stress 0.07378228 
    ## ... Procrustes: rmse 0.01746169  max resid 0.05108693 
    ## Run 425 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744624  max resid 0.05104806 
    ## Run 426 stress 0.07629244 
    ## Run 427 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001239916  max resid 0.0002963486 
    ## ... Similar to previous best
    ## Run 428 stress 0.08000467 
    ## Run 429 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001534498  max resid 0.0003635598 
    ## ... Similar to previous best
    ## Run 430 stress 0.07629234 
    ## Run 431 stress 0.0823377 
    ## Run 432 stress 0.07365783 
    ## ... Procrustes: rmse 2.882564e-05  max resid 6.608023e-05 
    ## ... Similar to previous best
    ## Run 433 stress 0.07365784 
    ## ... Procrustes: rmse 1.377898e-05  max resid 2.71211e-05 
    ## ... Similar to previous best
    ## Run 434 stress 0.08000468 
    ## Run 435 stress 0.07629233 
    ## Run 436 stress 0.07365784 
    ## ... Procrustes: rmse 3.58999e-05  max resid 8.440805e-05 
    ## ... Similar to previous best
    ## Run 437 stress 0.07629243 
    ## Run 438 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744497  max resid 0.05105862 
    ## Run 439 stress 0.07629235 
    ## Run 440 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744503  max resid 0.05104167 
    ## Run 441 stress 0.08000472 
    ## Run 442 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174452  max resid 0.05103167 
    ## Run 443 stress 0.07629237 
    ## Run 444 stress 0.07365785 
    ## ... Procrustes: rmse 7.110571e-05  max resid 0.0001674313 
    ## ... Similar to previous best
    ## Run 445 stress 0.07732931 
    ## Run 446 stress 0.08288047 
    ## Run 447 stress 0.08000471 
    ## Run 448 stress 0.07732926 
    ## Run 449 stress 0.07629235 
    ## Run 450 stress 0.07629243 
    ## Run 451 stress 0.07732928 
    ## Run 452 stress 0.08233755 
    ## Run 453 stress 0.07365783 
    ## ... Procrustes: rmse 3.09684e-05  max resid 7.325817e-05 
    ## ... Similar to previous best
    ## Run 454 stress 0.07629234 
    ## Run 455 stress 0.0762924 
    ## Run 456 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744383  max resid 0.05102244 
    ## Run 457 stress 0.08000471 
    ## Run 458 stress 0.07629238 
    ## Run 459 stress 0.0800047 
    ## Run 460 stress 0.07629233 
    ## Run 461 stress 0.07629235 
    ## Run 462 stress 0.08233756 
    ## Run 463 stress 0.07365784 
    ## ... Procrustes: rmse 2.953125e-05  max resid 6.923587e-05 
    ## ... Similar to previous best
    ## Run 464 stress 0.08233754 
    ## Run 465 stress 0.0737823 
    ## ... Procrustes: rmse 0.01746318  max resid 0.05114493 
    ## Run 466 stress 0.07629234 
    ## Run 467 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744731  max resid 0.05105607 
    ## Run 468 stress 0.07378231 
    ## ... Procrustes: rmse 0.01743682  max resid 0.05104223 
    ## Run 469 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001058017  max resid 0.0002514839 
    ## ... Similar to previous best
    ## Run 470 stress 0.07378228 
    ## ... Procrustes: rmse 0.01741758  max resid 0.05094829 
    ## Run 471 stress 0.08000467 
    ## Run 472 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744132  max resid 0.05105167 
    ## Run 473 stress 0.07629238 
    ## Run 474 stress 0.07365784 
    ## ... Procrustes: rmse 3.595133e-05  max resid 8.400811e-05 
    ## ... Similar to previous best
    ## Run 475 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745455  max resid 0.05109875 
    ## Run 476 stress 0.07629237 
    ## Run 477 stress 0.07365783 
    ## ... Procrustes: rmse 9.974601e-06  max resid 2.350766e-05 
    ## ... Similar to previous best
    ## Run 478 stress 0.07378229 
    ## ... Procrustes: rmse 0.0174476  max resid 0.05102852 
    ## Run 479 stress 0.07629245 
    ## Run 480 stress 0.08000468 
    ## Run 481 stress 0.07629235 
    ## Run 482 stress 0.07365786 
    ## ... Procrustes: rmse 9.654734e-05  max resid 0.0002275064 
    ## ... Similar to previous best
    ## Run 483 stress 0.07629235 
    ## Run 484 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744282  max resid 0.05102228 
    ## Run 485 stress 0.07629239 
    ## Run 486 stress 0.07365783 
    ## ... Procrustes: rmse 3.20369e-05  max resid 7.551237e-05 
    ## ... Similar to previous best
    ## Run 487 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174416  max resid 0.05100835 
    ## Run 488 stress 0.07629243 
    ## Run 489 stress 0.07365786 
    ## ... Procrustes: rmse 3.801549e-05  max resid 6.935055e-05 
    ## ... Similar to previous best
    ## Run 490 stress 0.07629237 
    ## Run 491 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744363  max resid 0.05103136 
    ## Run 492 stress 0.07365783 
    ## ... Procrustes: rmse 1.907446e-05  max resid 4.201158e-05 
    ## ... Similar to previous best
    ## Run 493 stress 0.08233769 
    ## Run 494 stress 0.0828805 
    ## Run 495 stress 0.0823376 
    ## Run 496 stress 0.08288054 
    ## Run 497 stress 0.08000477 
    ## Run 498 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001670829  max resid 0.0003931981 
    ## ... Similar to previous best
    ## Run 499 stress 0.07378235 
    ## ... Procrustes: rmse 0.01743664  max resid 0.05104825 
    ## Run 500 stress 0.07365783 
    ## ... Procrustes: rmse 7.21369e-06  max resid 1.683139e-05 
    ## ... Similar to previous best
    ## *** Best solution repeated 59 times

``` r
# Stratified lakes and ocean sites
SD_beta_SO_NMDS <- metaMDS(SD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07970526 
    ## Run 1 stress 0.07250812 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0899475  max resid 0.2576263 
    ## Run 2 stress 0.07428313 
    ## Run 3 stress 0.07250814 
    ## ... Procrustes: rmse 4.199043e-05  max resid 0.0001182409 
    ## ... Similar to previous best
    ## Run 4 stress 0.07250815 
    ## ... Procrustes: rmse 4.347845e-05  max resid 0.0001231815 
    ## ... Similar to previous best
    ## Run 5 stress 0.06978191 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04137704  max resid 0.1124036 
    ## Run 6 stress 0.08340287 
    ## Run 7 stress 0.06978196 
    ## ... Procrustes: rmse 4.697293e-05  max resid 0.0001285422 
    ## ... Similar to previous best
    ## Run 8 stress 0.07250812 
    ## Run 9 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01324925  max resid 0.03355772 
    ## Run 10 stress 0.083403 
    ## Run 11 stress 0.07428313 
    ## Run 12 stress 0.08340287 
    ## Run 13 stress 0.07428313 
    ## Run 14 stress 0.07970525 
    ## Run 15 stress 0.07970525 
    ## Run 16 stress 0.07428318 
    ## Run 17 stress 0.07250812 
    ## Run 18 stress 0.07844938 
    ## Run 19 stress 0.06942776 
    ## ... Procrustes: rmse 7.963591e-06  max resid 2.05168e-05 
    ## ... Similar to previous best
    ## Run 20 stress 0.08448443 
    ## Run 21 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324159  max resid 0.03329861 
    ## Run 22 stress 0.07428313 
    ## Run 23 stress 0.07428316 
    ## Run 24 stress 0.06942776 
    ## ... Procrustes: rmse 1.324091e-05  max resid 3.495705e-05 
    ## ... Similar to previous best
    ## Run 25 stress 0.0742832 
    ## Run 26 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324839  max resid 0.03331588 
    ## Run 27 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326639  max resid 0.03335204 
    ## Run 28 stress 0.08340286 
    ## Run 29 stress 0.07970526 
    ## Run 30 stress 0.07428314 
    ## Run 31 stress 0.08340286 
    ## Run 32 stress 0.07844915 
    ## Run 33 stress 0.06942777 
    ## ... Procrustes: rmse 5.115083e-05  max resid 0.0001332449 
    ## ... Similar to previous best
    ## Run 34 stress 0.07970525 
    ## Run 35 stress 0.08340288 
    ## Run 36 stress 0.07250812 
    ## Run 37 stress 0.07428313 
    ## Run 38 stress 0.06942777 
    ## ... Procrustes: rmse 4.817537e-05  max resid 0.0001242163 
    ## ... Similar to previous best
    ## Run 39 stress 0.07970525 
    ## Run 40 stress 0.08340289 
    ## Run 41 stress 0.06942776 
    ## ... Procrustes: rmse 9.924423e-06  max resid 2.556365e-05 
    ## ... Similar to previous best
    ## Run 42 stress 0.07250812 
    ## Run 43 stress 0.06942776 
    ## ... Procrustes: rmse 1.13887e-05  max resid 2.942742e-05 
    ## ... Similar to previous best
    ## Run 44 stress 0.08340287 
    ## Run 45 stress 0.07428312 
    ## Run 46 stress 0.06942776 
    ## ... Procrustes: rmse 1.015404e-05  max resid 2.630904e-05 
    ## ... Similar to previous best
    ## Run 47 stress 0.07250813 
    ## Run 48 stress 0.07970526 
    ## Run 49 stress 0.06978197 
    ## ... Procrustes: rmse 0.01317551  max resid 0.03316274 
    ## Run 50 stress 0.07428313 
    ## Run 51 stress 0.07428314 
    ## Run 52 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324473  max resid 0.03330758 
    ## Run 53 stress 0.07428319 
    ## Run 54 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325703  max resid 0.03333388 
    ## Run 55 stress 0.07970526 
    ## Run 56 stress 0.07970525 
    ## Run 57 stress 0.07428313 
    ## Run 58 stress 0.08448436 
    ## Run 59 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132308  max resid 0.03328036 
    ## Run 60 stress 0.07428313 
    ## Run 61 stress 0.06978196 
    ## ... Procrustes: rmse 0.0132647  max resid 0.03335245 
    ## Run 62 stress 0.06942777 
    ## ... Procrustes: rmse 3.251345e-05  max resid 8.303736e-05 
    ## ... Similar to previous best
    ## Run 63 stress 0.07970525 
    ## Run 64 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325111  max resid 0.03331983 
    ## Run 65 stress 0.0742832 
    ## Run 66 stress 0.07428316 
    ## Run 67 stress 0.07970526 
    ## Run 68 stress 0.08340286 
    ## Run 69 stress 0.08340287 
    ## Run 70 stress 0.07250812 
    ## Run 71 stress 0.06942776 
    ## ... Procrustes: rmse 6.066136e-06  max resid 1.670241e-05 
    ## ... Similar to previous best
    ## Run 72 stress 0.07428314 
    ## Run 73 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325829  max resid 0.03333578 
    ## Run 74 stress 0.07428316 
    ## Run 75 stress 0.06942776 
    ## ... Procrustes: rmse 1.023847e-05  max resid 3.107356e-05 
    ## ... Similar to previous best
    ## Run 76 stress 0.07970526 
    ## Run 77 stress 0.07250813 
    ## Run 78 stress 0.07970526 
    ## Run 79 stress 0.07250812 
    ## Run 80 stress 0.07250812 
    ## Run 81 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319437  max resid 0.03320194 
    ## Run 82 stress 0.08340289 
    ## Run 83 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322226  max resid 0.03325471 
    ## Run 84 stress 0.08340292 
    ## Run 85 stress 0.06942776 
    ## ... Procrustes: rmse 9.359232e-06  max resid 2.407648e-05 
    ## ... Similar to previous best
    ## Run 86 stress 0.06942776 
    ## ... Procrustes: rmse 3.286806e-05  max resid 8.504663e-05 
    ## ... Similar to previous best
    ## Run 87 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324906  max resid 0.03331623 
    ## Run 88 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325892  max resid 0.03333745 
    ## Run 89 stress 0.07428313 
    ## Run 90 stress 0.06978199 
    ## ... Procrustes: rmse 0.01328166  max resid 0.03338172 
    ## Run 91 stress 0.07844928 
    ## Run 92 stress 0.07970525 
    ## Run 93 stress 0.07428313 
    ## Run 94 stress 0.07428314 
    ## Run 95 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324846  max resid 0.03331464 
    ## Run 96 stress 0.06942779 
    ## ... Procrustes: rmse 2.53364e-05  max resid 7.585218e-05 
    ## ... Similar to previous best
    ## Run 97 stress 0.08448439 
    ## Run 98 stress 0.08340288 
    ## Run 99 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321003  max resid 0.03323276 
    ## Run 100 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323599  max resid 0.03329008 
    ## Run 101 stress 0.07428315 
    ## Run 102 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327532  max resid 0.03336953 
    ## Run 103 stress 0.07970525 
    ## Run 104 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319073  max resid 0.03319243 
    ## Run 105 stress 0.07250814 
    ## Run 106 stress 0.07428312 
    ## Run 107 stress 0.07250813 
    ## Run 108 stress 0.07428316 
    ## Run 109 stress 0.06978193 
    ## ... Procrustes: rmse 0.01318959  max resid 0.03319079 
    ## Run 110 stress 0.08340286 
    ## Run 111 stress 0.07970525 
    ## Run 112 stress 0.07428312 
    ## Run 113 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322164  max resid 0.03325758 
    ## Run 114 stress 0.07250813 
    ## Run 115 stress 0.0834029 
    ## Run 116 stress 0.08448443 
    ## Run 117 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320646  max resid 0.03322672 
    ## Run 118 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 3.330061e-06  max resid 6.743018e-06 
    ## ... Similar to previous best
    ## Run 119 stress 0.08340288 
    ## Run 120 stress 0.08340289 
    ## Run 121 stress 0.08340296 
    ## Run 122 stress 0.08340289 
    ## Run 123 stress 0.07970525 
    ## Run 124 stress 0.08340287 
    ## Run 125 stress 0.08340289 
    ## Run 126 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324048  max resid 0.03329915 
    ## Run 127 stress 0.07250814 
    ## Run 128 stress 0.07970525 
    ## Run 129 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325113  max resid 0.033319 
    ## Run 130 stress 0.07428313 
    ## Run 131 stress 0.07970525 
    ## Run 132 stress 0.06942777 
    ## ... Procrustes: rmse 3.878972e-05  max resid 0.0001011624 
    ## ... Similar to previous best
    ## Run 133 stress 0.07844956 
    ## Run 134 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320694  max resid 0.03322713 
    ## Run 135 stress 0.08340286 
    ## Run 136 stress 0.06942776 
    ## ... Procrustes: rmse 1.273111e-05  max resid 3.350004e-05 
    ## ... Similar to previous best
    ## Run 137 stress 0.07250812 
    ## Run 138 stress 0.07428314 
    ## Run 139 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322594  max resid 0.03326565 
    ## Run 140 stress 0.07428313 
    ## Run 141 stress 0.07250813 
    ## Run 142 stress 0.07428313 
    ## Run 143 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324198  max resid 0.03329941 
    ## Run 144 stress 0.06942776 
    ## ... Procrustes: rmse 8.707786e-06  max resid 2.258349e-05 
    ## ... Similar to previous best
    ## Run 145 stress 0.07250812 
    ## Run 146 stress 0.07250813 
    ## Run 147 stress 0.07250815 
    ## Run 148 stress 0.07970526 
    ## Run 149 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321758  max resid 0.03324678 
    ## Run 150 stress 0.08448451 
    ## Run 151 stress 0.07250813 
    ## Run 152 stress 0.08340293 
    ## Run 153 stress 0.07428317 
    ## Run 154 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325882  max resid 0.03333536 
    ## Run 155 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323819  max resid 0.03329235 
    ## Run 156 stress 0.07250812 
    ## Run 157 stress 0.07250812 
    ## Run 158 stress 0.07970525 
    ## Run 159 stress 0.06942776 
    ## ... Procrustes: rmse 9.896426e-06  max resid 2.636967e-05 
    ## ... Similar to previous best
    ## Run 160 stress 0.07428317 
    ## Run 161 stress 0.08340287 
    ## Run 162 stress 0.0834029 
    ## Run 163 stress 0.08448439 
    ## Run 164 stress 0.07428313 
    ## Run 165 stress 0.08340294 
    ## Run 166 stress 0.08340289 
    ## Run 167 stress 0.07428314 
    ## Run 168 stress 0.07428313 
    ## Run 169 stress 0.07428315 
    ## Run 170 stress 0.07428314 
    ## Run 171 stress 0.08340293 
    ## Run 172 stress 0.07970526 
    ## Run 173 stress 0.08448448 
    ## Run 174 stress 0.07970527 
    ## Run 175 stress 0.07250812 
    ## Run 176 stress 0.06942776 
    ## ... Procrustes: rmse 3.813708e-06  max resid 9.90938e-06 
    ## ... Similar to previous best
    ## Run 177 stress 0.07428315 
    ## Run 178 stress 0.07428313 
    ## Run 179 stress 0.06942776 
    ## ... Procrustes: rmse 2.14404e-05  max resid 5.557281e-05 
    ## ... Similar to previous best
    ## Run 180 stress 0.08340286 
    ## Run 181 stress 0.07428313 
    ## Run 182 stress 0.07970525 
    ## Run 183 stress 0.07970525 
    ## Run 184 stress 0.07428316 
    ## Run 185 stress 0.08340289 
    ## Run 186 stress 0.07428314 
    ## Run 187 stress 0.08340287 
    ## Run 188 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.947386e-06  max resid 6.196982e-06 
    ## ... Similar to previous best
    ## Run 189 stress 0.06942778 
    ## ... Procrustes: rmse 2.224989e-05  max resid 6.593271e-05 
    ## ... Similar to previous best
    ## Run 190 stress 0.07428314 
    ## Run 191 stress 0.07428313 
    ## Run 192 stress 0.07250812 
    ## Run 193 stress 0.06942776 
    ## ... Procrustes: rmse 1.485416e-05  max resid 3.881194e-05 
    ## ... Similar to previous best
    ## Run 194 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321381  max resid 0.03324399 
    ## Run 195 stress 0.07970525 
    ## Run 196 stress 0.07250813 
    ## Run 197 stress 0.07844943 
    ## Run 198 stress 0.07970525 
    ## Run 199 stress 0.07250812 
    ## Run 200 stress 0.07250812 
    ## Run 201 stress 0.07844964 
    ## Run 202 stress 0.07428314 
    ## Run 203 stress 0.07250812 
    ## Run 204 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326644  max resid 0.03335282 
    ## Run 205 stress 0.08448447 
    ## Run 206 stress 0.07844929 
    ## Run 207 stress 0.07250814 
    ## Run 208 stress 0.083403 
    ## Run 209 stress 0.07250814 
    ## Run 210 stress 0.069782 
    ## ... Procrustes: rmse 0.01328035  max resid 0.03338213 
    ## Run 211 stress 0.08340287 
    ## Run 212 stress 0.08340292 
    ## Run 213 stress 0.08340295 
    ## Run 214 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326498  max resid 0.03335238 
    ## Run 215 stress 0.07428313 
    ## Run 216 stress 0.06942776 
    ## ... Procrustes: rmse 3.243671e-05  max resid 8.37713e-05 
    ## ... Similar to previous best
    ## Run 217 stress 0.08340288 
    ## Run 218 stress 0.07428313 
    ## Run 219 stress 0.06978197 
    ## ... Procrustes: rmse 0.01321225  max resid 0.0332334 
    ## Run 220 stress 0.07250812 
    ## Run 221 stress 0.07428313 
    ## Run 222 stress 0.06978195 
    ## ... Procrustes: rmse 0.01317811  max resid 0.03316471 
    ## Run 223 stress 0.06942776 
    ## ... Procrustes: rmse 1.228096e-05  max resid 3.2947e-05 
    ## ... Similar to previous best
    ## Run 224 stress 0.07428313 
    ## Run 225 stress 0.0834029 
    ## Run 226 stress 0.07428315 
    ## Run 227 stress 0.07250812 
    ## Run 228 stress 0.08340293 
    ## Run 229 stress 0.07428312 
    ## Run 230 stress 0.07970525 
    ## Run 231 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322512  max resid 0.03326661 
    ## Run 232 stress 0.08340289 
    ## Run 233 stress 0.07250813 
    ## Run 234 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322724  max resid 0.03327171 
    ## Run 235 stress 0.08448437 
    ## Run 236 stress 0.08448451 
    ## Run 237 stress 0.06942776 
    ## ... Procrustes: rmse 7.853729e-06  max resid 2.043931e-05 
    ## ... Similar to previous best
    ## Run 238 stress 0.07250813 
    ## Run 239 stress 0.07428315 
    ## Run 240 stress 0.06942776 
    ## ... Procrustes: rmse 2.655938e-06  max resid 8.036295e-06 
    ## ... Similar to previous best
    ## Run 241 stress 0.07250812 
    ## Run 242 stress 0.07428314 
    ## Run 243 stress 0.07970526 
    ## Run 244 stress 0.07428319 
    ## Run 245 stress 0.08448439 
    ## Run 246 stress 0.08448444 
    ## Run 247 stress 0.06942777 
    ## ... Procrustes: rmse 4.007037e-05  max resid 0.0001033002 
    ## ... Similar to previous best
    ## Run 248 stress 0.07250812 
    ## Run 249 stress 0.07250812 
    ## Run 250 stress 0.07428312 
    ## Run 251 stress 0.06978191 
    ## ... Procrustes: rmse 0.01320198  max resid 0.03321606 
    ## Run 252 stress 0.08340289 
    ## Run 253 stress 0.06978196 
    ## ... Procrustes: rmse 0.01318113  max resid 0.03316925 
    ## Run 254 stress 0.07428319 
    ## Run 255 stress 0.07250819 
    ## Run 256 stress 0.08340295 
    ## Run 257 stress 0.0834029 
    ## Run 258 stress 0.07970525 
    ## Run 259 stress 0.06942776 
    ## ... Procrustes: rmse 1.484648e-05  max resid 3.84656e-05 
    ## ... Similar to previous best
    ## Run 260 stress 0.06978191 
    ## ... Procrustes: rmse 0.01320664  max resid 0.03322913 
    ## Run 261 stress 0.08340287 
    ## Run 262 stress 0.07250814 
    ## Run 263 stress 0.07970526 
    ## Run 264 stress 0.07428313 
    ## Run 265 stress 0.08448436 
    ## Run 266 stress 0.08340289 
    ## Run 267 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326411  max resid 0.03334657 
    ## Run 268 stress 0.08340287 
    ## Run 269 stress 0.07250812 
    ## Run 270 stress 0.06942777 
    ## ... Procrustes: rmse 5.344851e-05  max resid 0.0001376598 
    ## ... Similar to previous best
    ## Run 271 stress 0.08340287 
    ## Run 272 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320903  max resid 0.03323286 
    ## Run 273 stress 0.06978193 
    ## ... Procrustes: rmse 0.01324734  max resid 0.03331673 
    ## Run 274 stress 0.08340287 
    ## Run 275 stress 0.07970526 
    ## Run 276 stress 0.08340291 
    ## Run 277 stress 0.07428313 
    ## Run 278 stress 0.07970525 
    ## Run 279 stress 0.07428313 
    ## Run 280 stress 0.07428313 
    ## Run 281 stress 0.08340289 
    ## Run 282 stress 0.07970525 
    ## Run 283 stress 0.07428313 
    ## Run 284 stress 0.07250813 
    ## Run 285 stress 0.07250813 
    ## Run 286 stress 0.07250812 
    ## Run 287 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323897  max resid 0.03329516 
    ## Run 288 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325843  max resid 0.03333673 
    ## Run 289 stress 0.06978203 
    ## ... Procrustes: rmse 0.0132847  max resid 0.03339316 
    ## Run 290 stress 0.07250814 
    ## Run 291 stress 0.0784497 
    ## Run 292 stress 0.06978191 
    ## ... Procrustes: rmse 0.01321187  max resid 0.033242 
    ## Run 293 stress 0.07428314 
    ## Run 294 stress 0.07428315 
    ## Run 295 stress 0.07970525 
    ## Run 296 stress 0.07428314 
    ## Run 297 stress 0.06942776 
    ## ... Procrustes: rmse 6.044251e-06  max resid 1.537847e-05 
    ## ... Similar to previous best
    ## Run 298 stress 0.07428313 
    ## Run 299 stress 0.07844963 
    ## Run 300 stress 0.08448439 
    ## Run 301 stress 0.08448454 
    ## Run 302 stress 0.06978201 
    ## ... Procrustes: rmse 0.01328452  max resid 0.03338803 
    ## Run 303 stress 0.07970525 
    ## Run 304 stress 0.07250813 
    ## Run 305 stress 0.07428313 
    ## Run 306 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326479  max resid 0.03334555 
    ## Run 307 stress 0.07250812 
    ## Run 308 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323671  max resid 0.03329269 
    ## Run 309 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323877  max resid 0.03329728 
    ## Run 310 stress 0.06978194 
    ## ... Procrustes: rmse 0.01319506  max resid 0.03320625 
    ## Run 311 stress 0.08340291 
    ## Run 312 stress 0.07970525 
    ## Run 313 stress 0.07250812 
    ## Run 314 stress 0.06942776 
    ## ... Procrustes: rmse 1.86697e-06  max resid 3.430295e-06 
    ## ... Similar to previous best
    ## Run 315 stress 0.07428314 
    ## Run 316 stress 0.08340291 
    ## Run 317 stress 0.07250812 
    ## Run 318 stress 0.07428313 
    ## Run 319 stress 0.069782 
    ## ... Procrustes: rmse 0.01328416  max resid 0.03338885 
    ## Run 320 stress 0.06942776 
    ## ... Procrustes: rmse 2.987465e-05  max resid 7.71965e-05 
    ## ... Similar to previous best
    ## Run 321 stress 0.08340287 
    ## Run 322 stress 0.07428317 
    ## Run 323 stress 0.07428312 
    ## Run 324 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325844  max resid 0.03333803 
    ## Run 325 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325185  max resid 0.03332267 
    ## Run 326 stress 0.07250812 
    ## Run 327 stress 0.07970526 
    ## Run 328 stress 0.06942776 
    ## ... Procrustes: rmse 9.401718e-06  max resid 2.545446e-05 
    ## ... Similar to previous best
    ## Run 329 stress 0.07428313 
    ## Run 330 stress 0.08340289 
    ## Run 331 stress 0.07250812 
    ## Run 332 stress 0.06942776 
    ## ... Procrustes: rmse 2.82674e-05  max resid 7.34207e-05 
    ## ... Similar to previous best
    ## Run 333 stress 0.07428313 
    ## Run 334 stress 0.08340293 
    ## Run 335 stress 0.07428313 
    ## Run 336 stress 0.07970525 
    ## Run 337 stress 0.07250816 
    ## Run 338 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132116  max resid 0.03323857 
    ## Run 339 stress 0.0834029 
    ## Run 340 stress 0.07428317 
    ## Run 341 stress 0.07250813 
    ## Run 342 stress 0.07970525 
    ## Run 343 stress 0.08340297 
    ## Run 344 stress 0.07428316 
    ## Run 345 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322373  max resid 0.03326278 
    ## Run 346 stress 0.07428313 
    ## Run 347 stress 0.06942776 
    ## ... Procrustes: rmse 5.601619e-06  max resid 1.457747e-05 
    ## ... Similar to previous best
    ## Run 348 stress 0.06978199 
    ## ... Procrustes: rmse 0.01326726  max resid 0.03335919 
    ## Run 349 stress 0.07428319 
    ## Run 350 stress 0.07844943 
    ## Run 351 stress 0.08340286 
    ## Run 352 stress 0.08340289 
    ## Run 353 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327588  max resid 0.03337356 
    ## Run 354 stress 0.0834029 
    ## Run 355 stress 0.08340287 
    ## Run 356 stress 0.07428313 
    ## Run 357 stress 0.07250814 
    ## Run 358 stress 0.08448451 
    ## Run 359 stress 0.06942777 
    ## ... Procrustes: rmse 4.762339e-05  max resid 0.0001228909 
    ## ... Similar to previous best
    ## Run 360 stress 0.07250813 
    ## Run 361 stress 0.06942776 
    ## ... Procrustes: rmse 1.010745e-05  max resid 2.654519e-05 
    ## ... Similar to previous best
    ## Run 362 stress 0.06942776 
    ## ... Procrustes: rmse 4.657937e-06  max resid 1.196351e-05 
    ## ... Similar to previous best
    ## Run 363 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321521  max resid 0.0332464 
    ## Run 364 stress 0.06942777 
    ## ... Procrustes: rmse 3.288453e-05  max resid 8.633357e-05 
    ## ... Similar to previous best
    ## Run 365 stress 0.07970525 
    ## Run 366 stress 0.07250813 
    ## Run 367 stress 0.07970525 
    ## Run 368 stress 0.08448437 
    ## Run 369 stress 0.07970525 
    ## Run 370 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325303  max resid 0.03332664 
    ## Run 371 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322749  max resid 0.03326749 
    ## Run 372 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326193  max resid 0.03334391 
    ## Run 373 stress 0.08340286 
    ## Run 374 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325027  max resid 0.03331997 
    ## Run 375 stress 0.07250812 
    ## Run 376 stress 0.07428316 
    ## Run 377 stress 0.07428314 
    ## Run 378 stress 0.07250812 
    ## Run 379 stress 0.06942776 
    ## ... Procrustes: rmse 1.923011e-05  max resid 5.172693e-05 
    ## ... Similar to previous best
    ## Run 380 stress 0.07970525 
    ## Run 381 stress 0.07428314 
    ## Run 382 stress 0.07428314 
    ## Run 383 stress 0.07428314 
    ## Run 384 stress 0.08340288 
    ## Run 385 stress 0.0697819 
    ## ... Procrustes: rmse 0.01324131  max resid 0.03330186 
    ## Run 386 stress 0.06978196 
    ## ... Procrustes: rmse 0.01317892  max resid 0.03317008 
    ## Run 387 stress 0.07428313 
    ## Run 388 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323491  max resid 0.03328682 
    ## Run 389 stress 0.07250812 
    ## Run 390 stress 0.07428313 
    ## Run 391 stress 0.07428318 
    ## Run 392 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326233  max resid 0.0333438 
    ## Run 393 stress 0.07970526 
    ## Run 394 stress 0.07970525 
    ## Run 395 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325614  max resid 0.03333202 
    ## Run 396 stress 0.07428323 
    ## Run 397 stress 0.07250814 
    ## Run 398 stress 0.06942777 
    ## ... Procrustes: rmse 4.403181e-05  max resid 0.0001134433 
    ## ... Similar to previous best
    ## Run 399 stress 0.06978195 
    ## ... Procrustes: rmse 0.0131778  max resid 0.03316569 
    ## Run 400 stress 0.07428314 
    ## Run 401 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323245  max resid 0.03327971 
    ## Run 402 stress 0.07250814 
    ## Run 403 stress 0.07250813 
    ## Run 404 stress 0.08448441 
    ## Run 405 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322532  max resid 0.0332692 
    ## Run 406 stress 0.06978205 
    ## ... Procrustes: rmse 0.01328869  max resid 0.0333941 
    ## Run 407 stress 0.08340288 
    ## Run 408 stress 0.07428313 
    ## Run 409 stress 0.07250815 
    ## Run 410 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323581  max resid 0.03328973 
    ## Run 411 stress 0.06978192 
    ## ... Procrustes: rmse 0.013252  max resid 0.03332279 
    ## Run 412 stress 0.07428316 
    ## Run 413 stress 0.07844915 
    ## Run 414 stress 0.08340287 
    ## Run 415 stress 0.06942777 
    ## ... Procrustes: rmse 3.427908e-05  max resid 8.89867e-05 
    ## ... Similar to previous best
    ## Run 416 stress 0.08340289 
    ## Run 417 stress 0.07428315 
    ## Run 418 stress 0.08340289 
    ## Run 419 stress 0.07250813 
    ## Run 420 stress 0.08340289 
    ## Run 421 stress 0.07250825 
    ## Run 422 stress 0.07970526 
    ## Run 423 stress 0.07970525 
    ## Run 424 stress 0.07970528 
    ## Run 425 stress 0.07250813 
    ## Run 426 stress 0.07970525 
    ## Run 427 stress 0.07428314 
    ## Run 428 stress 0.07970525 
    ## Run 429 stress 0.07970525 
    ## Run 430 stress 0.08340287 
    ## Run 431 stress 0.08340294 
    ## Run 432 stress 0.07428313 
    ## Run 433 stress 0.0834029 
    ## Run 434 stress 0.06978203 
    ## ... Procrustes: rmse 0.01315537  max resid 0.03311407 
    ## Run 435 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.465427e-06  max resid 6.378422e-06 
    ## ... Similar to previous best
    ## Run 436 stress 0.07428313 
    ## Run 437 stress 0.08340287 
    ## Run 438 stress 0.07970525 
    ## Run 439 stress 0.07970525 
    ## Run 440 stress 0.08340293 
    ## Run 441 stress 0.06942777 
    ## ... Procrustes: rmse 3.846893e-05  max resid 9.955981e-05 
    ## ... Similar to previous best
    ## Run 442 stress 0.07250816 
    ## Run 443 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322465  max resid 0.03326583 
    ## Run 444 stress 0.06942776 
    ## ... Procrustes: rmse 2.340873e-05  max resid 6.090124e-05 
    ## ... Similar to previous best
    ## Run 445 stress 0.07970525 
    ## Run 446 stress 0.07970526 
    ## Run 447 stress 0.07250812 
    ## Run 448 stress 0.07250812 
    ## Run 449 stress 0.07844918 
    ## Run 450 stress 0.07970525 
    ## Run 451 stress 0.08448439 
    ## Run 452 stress 0.06942777 
    ## ... Procrustes: rmse 3.969914e-05  max resid 0.0001021573 
    ## ... Similar to previous best
    ## Run 453 stress 0.07250812 
    ## Run 454 stress 0.07428313 
    ## Run 455 stress 0.07970525 
    ## Run 456 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327216  max resid 0.03336577 
    ## Run 457 stress 0.07970525 
    ## Run 458 stress 0.07428314 
    ## Run 459 stress 0.07970525 
    ## Run 460 stress 0.07428318 
    ## Run 461 stress 0.07428313 
    ## Run 462 stress 0.06942777 
    ## ... Procrustes: rmse 4.391223e-05  max resid 0.0001131922 
    ## ... Similar to previous best
    ## Run 463 stress 0.06942776 
    ## ... Procrustes: rmse 5.32646e-06  max resid 1.404398e-05 
    ## ... Similar to previous best
    ## Run 464 stress 0.07970525 
    ## Run 465 stress 0.07428314 
    ## Run 466 stress 0.08340287 
    ## Run 467 stress 0.07428317 
    ## Run 468 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327823  max resid 0.0333764 
    ## Run 469 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322718  max resid 0.03327119 
    ## Run 470 stress 0.08448442 
    ## Run 471 stress 0.08340287 
    ## Run 472 stress 0.06978196 
    ## ... Procrustes: rmse 0.01317614  max resid 0.0331617 
    ## Run 473 stress 0.08448453 
    ## Run 474 stress 0.07970526 
    ## Run 475 stress 0.07428313 
    ## Run 476 stress 0.07428315 
    ## Run 477 stress 0.07970525 
    ## Run 478 stress 0.08340293 
    ## Run 479 stress 0.07428315 
    ## Run 480 stress 0.07250812 
    ## Run 481 stress 0.07428314 
    ## Run 482 stress 0.07970525 
    ## Run 483 stress 0.06942776 
    ## ... Procrustes: rmse 1.174934e-05  max resid 3.040674e-05 
    ## ... Similar to previous best
    ## Run 484 stress 0.06942777 
    ## ... Procrustes: rmse 4.176892e-05  max resid 0.0001075386 
    ## ... Similar to previous best
    ## Run 485 stress 0.08340287 
    ## Run 486 stress 0.07428321 
    ## Run 487 stress 0.08340293 
    ## Run 488 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326977  max resid 0.03336292 
    ## Run 489 stress 0.07970525 
    ## Run 490 stress 0.07428319 
    ## Run 491 stress 0.08448451 
    ## Run 492 stress 0.07428313 
    ## Run 493 stress 0.07970525 
    ## Run 494 stress 0.06942776 
    ## ... Procrustes: rmse 3.721135e-06  max resid 1.072049e-05 
    ## ... Similar to previous best
    ## Run 495 stress 0.07844959 
    ## Run 496 stress 0.07428316 
    ## Run 497 stress 0.07250814 
    ## Run 498 stress 0.0844844 
    ## Run 499 stress 0.06942776 
    ## ... Procrustes: rmse 2.415923e-05  max resid 6.226481e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.07250815 
    ## *** Best solution repeated 10 times

``` r
# Stratified lakes
SD_beta_S_NMDS <- metaMDS(SD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1410449 
    ## Run 1 stress 0.1551863 
    ## Run 2 stress 0.1628606 
    ## Run 3 stress 0.1572668 
    ## Run 4 stress 0.1572668 
    ## Run 5 stress 0.1415299 
    ## ... Procrustes: rmse 0.1479839  max resid 0.2321187 
    ## Run 6 stress 0.1410448 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003239049  max resid 0.0007641685 
    ## ... Similar to previous best
    ## Run 7 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 0.2673437  max resid 0.5426798 
    ## Run 8 stress 0.1320258 
    ## ... Procrustes: rmse 3.59266e-06  max resid 7.699449e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.1383681 
    ## Run 10 stress 0.1572668 
    ## Run 11 stress 0.1970303 
    ## Run 12 stress 0.1693717 
    ## Run 13 stress 0.249638 
    ## Run 14 stress 0.1445811 
    ## Run 15 stress 0.1693717 
    ## Run 16 stress 0.1320258 
    ## ... Procrustes: rmse 1.516948e-06  max resid 2.786313e-06 
    ## ... Similar to previous best
    ## Run 17 stress 0.2224298 
    ## Run 18 stress 0.1407298 
    ## Run 19 stress 0.1320258 
    ## ... Procrustes: rmse 2.412527e-06  max resid 4.565481e-06 
    ## ... Similar to previous best
    ## Run 20 stress 0.1693717 
    ## Run 21 stress 0.1407206 
    ## Run 22 stress 0.1383681 
    ## Run 23 stress 0.1415299 
    ## Run 24 stress 0.1320258 
    ## ... Procrustes: rmse 2.546922e-06  max resid 5.522397e-06 
    ## ... Similar to previous best
    ## Run 25 stress 0.1693717 
    ## Run 26 stress 0.1410448 
    ## Run 27 stress 0.1693717 
    ## Run 28 stress 0.1415299 
    ## Run 29 stress 0.141045 
    ## Run 30 stress 0.1383682 
    ## Run 31 stress 0.2422742 
    ## Run 32 stress 0.1410448 
    ## Run 33 stress 0.1830332 
    ## Run 34 stress 0.1410448 
    ## Run 35 stress 0.2085297 
    ## Run 36 stress 0.2224298 
    ## Run 37 stress 0.2842805 
    ## Run 38 stress 0.1407298 
    ## Run 39 stress 0.1598154 
    ## Run 40 stress 0.1776793 
    ## Run 41 stress 0.1410448 
    ## Run 42 stress 0.1434197 
    ## Run 43 stress 0.1407298 
    ## Run 44 stress 0.1434197 
    ## Run 45 stress 0.3072832 
    ## Run 46 stress 0.1407208 
    ## Run 47 stress 0.1383682 
    ## Run 48 stress 0.1415299 
    ## Run 49 stress 0.2950733 
    ## Run 50 stress 0.1410449 
    ## Run 51 stress 0.1410448 
    ## Run 52 stress 0.1320258 
    ## ... Procrustes: rmse 6.568987e-07  max resid 1.11291e-06 
    ## ... Similar to previous best
    ## Run 53 stress 0.1771966 
    ## Run 54 stress 0.1771966 
    ## Run 55 stress 0.1407298 
    ## Run 56 stress 0.1434197 
    ## Run 57 stress 0.2562346 
    ## Run 58 stress 0.1407207 
    ## Run 59 stress 0.1693717 
    ## Run 60 stress 0.1572668 
    ## Run 61 stress 0.1407298 
    ## Run 62 stress 0.1383681 
    ## Run 63 stress 0.1410451 
    ## Run 64 stress 0.1410449 
    ## Run 65 stress 0.1415299 
    ## Run 66 stress 0.1320258 
    ## ... Procrustes: rmse 7.259252e-07  max resid 1.607847e-06 
    ## ... Similar to previous best
    ## Run 67 stress 0.1407206 
    ## Run 68 stress 0.1410448 
    ## Run 69 stress 0.1407298 
    ## Run 70 stress 0.1583016 
    ## Run 71 stress 0.1415299 
    ## Run 72 stress 0.1410448 
    ## Run 73 stress 0.2385029 
    ## Run 74 stress 0.1445811 
    ## Run 75 stress 0.1410449 
    ## Run 76 stress 0.1572668 
    ## Run 77 stress 0.1383681 
    ## Run 78 stress 0.1415299 
    ## Run 79 stress 0.1407208 
    ## Run 80 stress 0.1771969 
    ## Run 81 stress 0.1410449 
    ## Run 82 stress 0.1415299 
    ## Run 83 stress 0.1551863 
    ## Run 84 stress 0.1383681 
    ## Run 85 stress 0.1410455 
    ## Run 86 stress 0.1320258 
    ## ... Procrustes: rmse 1.001653e-06  max resid 1.736392e-06 
    ## ... Similar to previous best
    ## Run 87 stress 0.141045 
    ## Run 88 stress 0.1445811 
    ## Run 89 stress 0.1415299 
    ## Run 90 stress 0.1587623 
    ## Run 91 stress 0.1383682 
    ## Run 92 stress 0.1415299 
    ## Run 93 stress 0.1410448 
    ## Run 94 stress 0.1415299 
    ## Run 95 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 5.422313e-07  max resid 9.163433e-07 
    ## ... Similar to previous best
    ## Run 96 stress 0.1796865 
    ## Run 97 stress 0.1415299 
    ## Run 98 stress 0.1407298 
    ## Run 99 stress 0.1415299 
    ## Run 100 stress 0.1663544 
    ## Run 101 stress 0.1320258 
    ## ... Procrustes: rmse 2.671851e-06  max resid 5.500391e-06 
    ## ... Similar to previous best
    ## Run 102 stress 0.2422668 
    ## Run 103 stress 0.1776793 
    ## Run 104 stress 0.1383681 
    ## Run 105 stress 0.1628606 
    ## Run 106 stress 0.1998955 
    ## Run 107 stress 0.1410448 
    ## Run 108 stress 0.1415299 
    ## Run 109 stress 0.1407298 
    ## Run 110 stress 0.1320258 
    ## ... Procrustes: rmse 1.833585e-06  max resid 3.862342e-06 
    ## ... Similar to previous best
    ## Run 111 stress 0.1410448 
    ## Run 112 stress 0.1415299 
    ## Run 113 stress 0.1383681 
    ## Run 114 stress 0.1551863 
    ## Run 115 stress 0.1410454 
    ## Run 116 stress 0.2008158 
    ## Run 117 stress 0.2629507 
    ## Run 118 stress 0.1407207 
    ## Run 119 stress 0.1320258 
    ## ... Procrustes: rmse 1.838452e-06  max resid 3.702733e-06 
    ## ... Similar to previous best
    ## Run 120 stress 0.1410452 
    ## Run 121 stress 0.1407298 
    ## Run 122 stress 0.1410449 
    ## Run 123 stress 0.1598154 
    ## Run 124 stress 0.1572668 
    ## Run 125 stress 0.2008162 
    ## Run 126 stress 0.1572668 
    ## Run 127 stress 0.1383681 
    ## Run 128 stress 0.1320258 
    ## ... Procrustes: rmse 1.071797e-06  max resid 1.856955e-06 
    ## ... Similar to previous best
    ## Run 129 stress 0.1320258 
    ## ... Procrustes: rmse 1.438811e-06  max resid 2.93209e-06 
    ## ... Similar to previous best
    ## Run 130 stress 0.1771969 
    ## Run 131 stress 0.2615852 
    ## Run 132 stress 0.1320258 
    ## ... Procrustes: rmse 9.230304e-07  max resid 1.923927e-06 
    ## ... Similar to previous best
    ## Run 133 stress 0.1771966 
    ## Run 134 stress 0.1383682 
    ## Run 135 stress 0.1383682 
    ## Run 136 stress 0.2373121 
    ## Run 137 stress 0.1407298 
    ## Run 138 stress 0.1771969 
    ## Run 139 stress 0.1415299 
    ## Run 140 stress 0.1415299 
    ## Run 141 stress 0.1415299 
    ## Run 142 stress 0.1320258 
    ## ... Procrustes: rmse 2.153778e-06  max resid 4.387058e-06 
    ## ... Similar to previous best
    ## Run 143 stress 0.1320258 
    ## ... Procrustes: rmse 9.388027e-07  max resid 1.854531e-06 
    ## ... Similar to previous best
    ## Run 144 stress 0.1407209 
    ## Run 145 stress 0.1830332 
    ## Run 146 stress 0.1415299 
    ## Run 147 stress 0.2422667 
    ## Run 148 stress 0.1383682 
    ## Run 149 stress 0.2224298 
    ## Run 150 stress 0.1434197 
    ## Run 151 stress 0.1383682 
    ## Run 152 stress 0.1320258 
    ## ... Procrustes: rmse 1.802363e-06  max resid 3.599801e-06 
    ## ... Similar to previous best
    ## Run 153 stress 0.1407206 
    ## Run 154 stress 0.1410449 
    ## Run 155 stress 0.1383682 
    ## Run 156 stress 0.1551863 
    ## Run 157 stress 0.1320258 
    ## ... Procrustes: rmse 6.695167e-07  max resid 1.114489e-06 
    ## ... Similar to previous best
    ## Run 158 stress 0.1383681 
    ## Run 159 stress 0.1383681 
    ## Run 160 stress 0.1383682 
    ## Run 161 stress 0.1383682 
    ## Run 162 stress 0.1434197 
    ## Run 163 stress 0.1383681 
    ## Run 164 stress 0.2501063 
    ## Run 165 stress 0.1320258 
    ## ... Procrustes: rmse 1.082182e-06  max resid 2.182385e-06 
    ## ... Similar to previous best
    ## Run 166 stress 0.2711712 
    ## Run 167 stress 0.2775528 
    ## Run 168 stress 0.1415299 
    ## Run 169 stress 0.2578639 
    ## Run 170 stress 0.2422667 
    ## Run 171 stress 0.1693717 
    ## Run 172 stress 0.141045 
    ## Run 173 stress 0.1572668 
    ## Run 174 stress 0.1776793 
    ## Run 175 stress 0.1572668 
    ## Run 176 stress 0.1410449 
    ## Run 177 stress 0.2416997 
    ## Run 178 stress 0.1410454 
    ## Run 179 stress 0.1929519 
    ## Run 180 stress 0.1572668 
    ## Run 181 stress 0.1415299 
    ## Run 182 stress 0.1407298 
    ## Run 183 stress 0.1320258 
    ## ... Procrustes: rmse 1.650854e-06  max resid 3.361897e-06 
    ## ... Similar to previous best
    ## Run 184 stress 0.1320258 
    ## ... Procrustes: rmse 2.226078e-06  max resid 4.589255e-06 
    ## ... Similar to previous best
    ## Run 185 stress 0.1410449 
    ## Run 186 stress 0.1796864 
    ## Run 187 stress 0.1572668 
    ## Run 188 stress 0.1410448 
    ## Run 189 stress 0.1771969 
    ## Run 190 stress 0.2385029 
    ## Run 191 stress 0.1598154 
    ## Run 192 stress 0.1407298 
    ## Run 193 stress 0.1320258 
    ## ... Procrustes: rmse 2.164903e-06  max resid 4.423368e-06 
    ## ... Similar to previous best
    ## Run 194 stress 0.2496383 
    ## Run 195 stress 0.1693717 
    ## Run 196 stress 0.1628606 
    ## Run 197 stress 0.1551863 
    ## Run 198 stress 0.1415299 
    ## Run 199 stress 0.1693717 
    ## Run 200 stress 0.1415299 
    ## Run 201 stress 0.1320258 
    ## ... Procrustes: rmse 1.710018e-06  max resid 3.459675e-06 
    ## ... Similar to previous best
    ## Run 202 stress 0.2085297 
    ## Run 203 stress 0.1320258 
    ## ... Procrustes: rmse 2.111973e-06  max resid 4.202859e-06 
    ## ... Similar to previous best
    ## Run 204 stress 0.1415299 
    ## Run 205 stress 0.1383681 
    ## Run 206 stress 0.1802587 
    ## Run 207 stress 0.1410449 
    ## Run 208 stress 0.1383681 
    ## Run 209 stress 0.1830332 
    ## Run 210 stress 0.1320258 
    ## ... Procrustes: rmse 4.504367e-07  max resid 6.560955e-07 
    ## ... Similar to previous best
    ## Run 211 stress 0.2501064 
    ## Run 212 stress 0.1802587 
    ## Run 213 stress 0.2385029 
    ## Run 214 stress 0.2562346 
    ## Run 215 stress 0.1802587 
    ## Run 216 stress 0.1320258 
    ## ... Procrustes: rmse 7.290918e-07  max resid 1.365046e-06 
    ## ... Similar to previous best
    ## Run 217 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 3.483495e-07  max resid 5.298438e-07 
    ## ... Similar to previous best
    ## Run 218 stress 0.1628606 
    ## Run 219 stress 0.1970303 
    ## Run 220 stress 0.2373121 
    ## Run 221 stress 0.1962476 
    ## Run 222 stress 0.1410451 
    ## Run 223 stress 0.1434197 
    ## Run 224 stress 0.1410453 
    ## Run 225 stress 0.1572668 
    ## Run 226 stress 0.1383682 
    ## Run 227 stress 0.1320258 
    ## ... Procrustes: rmse 9.624235e-07  max resid 1.847014e-06 
    ## ... Similar to previous best
    ## Run 228 stress 0.2385029 
    ## Run 229 stress 0.2005461 
    ## Run 230 stress 0.1320258 
    ## ... Procrustes: rmse 2.433644e-06  max resid 4.772946e-06 
    ## ... Similar to previous best
    ## Run 231 stress 0.1693717 
    ## Run 232 stress 0.1383681 
    ## Run 233 stress 0.1693718 
    ## Run 234 stress 0.1415299 
    ## Run 235 stress 0.2373121 
    ## Run 236 stress 0.1434197 
    ## Run 237 stress 0.2384387 
    ## Run 238 stress 0.1551863 
    ## Run 239 stress 0.1383681 
    ## Run 240 stress 0.1693717 
    ## Run 241 stress 0.1407298 
    ## Run 242 stress 0.1551863 
    ## Run 243 stress 0.1572668 
    ## Run 244 stress 0.1320258 
    ## ... Procrustes: rmse 2.43175e-06  max resid 5.028085e-06 
    ## ... Similar to previous best
    ## Run 245 stress 0.1410448 
    ## Run 246 stress 0.2848515 
    ## Run 247 stress 0.1598154 
    ## Run 248 stress 0.1587623 
    ## Run 249 stress 0.1598154 
    ## Run 250 stress 0.1415299 
    ## Run 251 stress 0.1581642 
    ## Run 252 stress 0.1771969 
    ## Run 253 stress 0.1415299 
    ## Run 254 stress 0.1551863 
    ## Run 255 stress 0.2422742 
    ## Run 256 stress 0.2085297 
    ## Run 257 stress 0.1320258 
    ## ... Procrustes: rmse 1.294108e-06  max resid 2.578885e-06 
    ## ... Similar to previous best
    ## Run 258 stress 0.2501063 
    ## Run 259 stress 0.2373121 
    ## Run 260 stress 0.1383681 
    ## Run 261 stress 0.141045 
    ## Run 262 stress 0.1407298 
    ## Run 263 stress 0.1320258 
    ## ... Procrustes: rmse 1.49471e-06  max resid 2.94546e-06 
    ## ... Similar to previous best
    ## Run 264 stress 0.1771966 
    ## Run 265 stress 0.1583016 
    ## Run 266 stress 0.1551863 
    ## Run 267 stress 0.3083098 
    ## Run 268 stress 0.1320258 
    ## ... Procrustes: rmse 1.435209e-06  max resid 3.056381e-06 
    ## ... Similar to previous best
    ## Run 269 stress 0.1830332 
    ## Run 270 stress 0.1410449 
    ## Run 271 stress 0.1407207 
    ## Run 272 stress 0.1415299 
    ## Run 273 stress 0.1693717 
    ## Run 274 stress 0.2496389 
    ## Run 275 stress 0.1771969 
    ## Run 276 stress 0.1693718 
    ## Run 277 stress 0.2422742 
    ## Run 278 stress 0.1407298 
    ## Run 279 stress 0.1551863 
    ## Run 280 stress 0.1415299 
    ## Run 281 stress 0.1410448 
    ## Run 282 stress 0.1551863 
    ## Run 283 stress 0.1320258 
    ## ... Procrustes: rmse 1.590267e-06  max resid 2.937537e-06 
    ## ... Similar to previous best
    ## Run 284 stress 0.1551863 
    ## Run 285 stress 0.1551863 
    ## Run 286 stress 0.141045 
    ## Run 287 stress 0.1776793 
    ## Run 288 stress 0.1410448 
    ## Run 289 stress 0.2578644 
    ## Run 290 stress 0.1551863 
    ## Run 291 stress 0.2501065 
    ## Run 292 stress 0.1320258 
    ## ... Procrustes: rmse 1.530572e-05  max resid 3.160138e-05 
    ## ... Similar to previous best
    ## Run 293 stress 0.2015531 
    ## Run 294 stress 0.1383682 
    ## Run 295 stress 0.1320258 
    ## ... Procrustes: rmse 1.343782e-06  max resid 2.663951e-06 
    ## ... Similar to previous best
    ## Run 296 stress 0.3083098 
    ## Run 297 stress 0.1407298 
    ## Run 298 stress 0.1551863 
    ## Run 299 stress 0.2496385 
    ## Run 300 stress 0.1572668 
    ## Run 301 stress 0.1583016 
    ## Run 302 stress 0.1583016 
    ## Run 303 stress 0.1415299 
    ## Run 304 stress 0.1320258 
    ## ... Procrustes: rmse 9.345357e-07  max resid 1.84274e-06 
    ## ... Similar to previous best
    ## Run 305 stress 0.1383681 
    ## Run 306 stress 0.1320258 
    ## ... Procrustes: rmse 1.790288e-06  max resid 3.70563e-06 
    ## ... Similar to previous best
    ## Run 307 stress 0.1693717 
    ## Run 308 stress 0.2385029 
    ## Run 309 stress 0.1572668 
    ## Run 310 stress 0.2118435 
    ## Run 311 stress 0.2385029 
    ## Run 312 stress 0.1551863 
    ## Run 313 stress 0.1572668 
    ## Run 314 stress 0.1551863 
    ## Run 315 stress 0.1410448 
    ## Run 316 stress 0.2074935 
    ## Run 317 stress 0.1693717 
    ## Run 318 stress 0.1693717 
    ## Run 319 stress 0.1572668 
    ## Run 320 stress 0.1572668 
    ## Run 321 stress 0.1434197 
    ## Run 322 stress 0.1572668 
    ## Run 323 stress 0.1410452 
    ## Run 324 stress 0.2842805 
    ## Run 325 stress 0.2015531 
    ## Run 326 stress 0.1383682 
    ## Run 327 stress 0.1572668 
    ## Run 328 stress 0.1383681 
    ## Run 329 stress 0.1407298 
    ## Run 330 stress 0.1693717 
    ## Run 331 stress 0.1572668 
    ## Run 332 stress 0.2848522 
    ## Run 333 stress 0.1410449 
    ## Run 334 stress 0.1410449 
    ## Run 335 stress 0.2550045 
    ## Run 336 stress 0.1410449 
    ## Run 337 stress 0.2224298 
    ## Run 338 stress 0.1693717 
    ## Run 339 stress 0.1320258 
    ## ... Procrustes: rmse 4.230667e-06  max resid 8.504216e-06 
    ## ... Similar to previous best
    ## Run 340 stress 0.1581642 
    ## Run 341 stress 0.1407207 
    ## Run 342 stress 0.1410449 
    ## Run 343 stress 0.1407298 
    ## Run 344 stress 0.1320258 
    ## ... Procrustes: rmse 4.385412e-07  max resid 7.115145e-07 
    ## ... Similar to previous best
    ## Run 345 stress 0.2005461 
    ## Run 346 stress 0.2008159 
    ## Run 347 stress 0.1407207 
    ## Run 348 stress 0.1410449 
    ## Run 349 stress 0.1970303 
    ## Run 350 stress 0.1581642 
    ## Run 351 stress 0.1771966 
    ## Run 352 stress 0.1572668 
    ## Run 353 stress 0.1598154 
    ## Run 354 stress 0.1410452 
    ## Run 355 stress 0.1415299 
    ## Run 356 stress 0.1415299 
    ## Run 357 stress 0.1551863 
    ## Run 358 stress 0.2385029 
    ## Run 359 stress 0.1407298 
    ## Run 360 stress 0.1407298 
    ## Run 361 stress 0.1771969 
    ## Run 362 stress 0.1587623 
    ## Run 363 stress 0.1407298 
    ## Run 364 stress 0.2422667 
    ## Run 365 stress 0.1410449 
    ## Run 366 stress 0.1410449 
    ## Run 367 stress 0.1415299 
    ## Run 368 stress 0.257865 
    ## Run 369 stress 0.2008163 
    ## Run 370 stress 0.1415299 
    ## Run 371 stress 0.1693717 
    ## Run 372 stress 0.2951814 
    ## Run 373 stress 0.1320258 
    ## ... Procrustes: rmse 6.109703e-07  max resid 1.312944e-06 
    ## ... Similar to previous best
    ## Run 374 stress 0.3083098 
    ## Run 375 stress 0.1572668 
    ## Run 376 stress 0.1771969 
    ## Run 377 stress 0.1551863 
    ## Run 378 stress 0.1320258 
    ## ... Procrustes: rmse 1.705109e-06  max resid 3.586657e-06 
    ## ... Similar to previous best
    ## Run 379 stress 0.1320258 
    ## ... Procrustes: rmse 1.722519e-07  max resid 2.442231e-07 
    ## ... Similar to previous best
    ## Run 380 stress 0.1320258 
    ## ... Procrustes: rmse 1.635755e-06  max resid 2.774688e-06 
    ## ... Similar to previous best
    ## Run 381 stress 0.1693717 
    ## Run 382 stress 0.2373121 
    ## Run 383 stress 0.1320258 
    ## ... Procrustes: rmse 2.904771e-07  max resid 6.346359e-07 
    ## ... Similar to previous best
    ## Run 384 stress 0.1410449 
    ## Run 385 stress 0.1415299 
    ## Run 386 stress 0.1415299 
    ## Run 387 stress 0.1383682 
    ## Run 388 stress 0.1693717 
    ## Run 389 stress 0.1962476 
    ## Run 390 stress 0.1572668 
    ## Run 391 stress 0.2456775 
    ## Run 392 stress 0.2227526 
    ## Run 393 stress 0.1587623 
    ## Run 394 stress 0.1598154 
    ## Run 395 stress 0.1551863 
    ## Run 396 stress 0.2385029 
    ## Run 397 stress 0.1663544 
    ## Run 398 stress 0.1410449 
    ## Run 399 stress 0.1407298 
    ## Run 400 stress 0.1410449 
    ## Run 401 stress 0.2538018 
    ## Run 402 stress 0.1572668 
    ## Run 403 stress 0.1415299 
    ## Run 404 stress 0.1320258 
    ## ... Procrustes: rmse 2.248625e-06  max resid 4.642025e-06 
    ## ... Similar to previous best
    ## Run 405 stress 0.1598154 
    ## Run 406 stress 0.1320258 
    ## ... Procrustes: rmse 1.523657e-06  max resid 3.16616e-06 
    ## ... Similar to previous best
    ## Run 407 stress 0.1771969 
    ## Run 408 stress 0.1383681 
    ## Run 409 stress 0.1410448 
    ## Run 410 stress 0.1320258 
    ## ... Procrustes: rmse 7.53688e-07  max resid 1.56911e-06 
    ## ... Similar to previous best
    ## Run 411 stress 0.1415299 
    ## Run 412 stress 0.1320258 
    ## ... Procrustes: rmse 2.050751e-06  max resid 4.282026e-06 
    ## ... Similar to previous best
    ## Run 413 stress 0.1410449 
    ## Run 414 stress 0.1572668 
    ## Run 415 stress 0.1410449 
    ## Run 416 stress 0.1415299 
    ## Run 417 stress 0.1383681 
    ## Run 418 stress 0.2384387 
    ## Run 419 stress 0.1410449 
    ## Run 420 stress 0.1410449 
    ## Run 421 stress 0.1572668 
    ## Run 422 stress 0.1410449 
    ## Run 423 stress 0.1410451 
    ## Run 424 stress 0.1407298 
    ## Run 425 stress 0.2008165 
    ## Run 426 stress 0.1410452 
    ## Run 427 stress 0.2842805 
    ## Run 428 stress 0.2422742 
    ## Run 429 stress 0.1320258 
    ## ... Procrustes: rmse 6.552488e-07  max resid 1.304476e-06 
    ## ... Similar to previous best
    ## Run 430 stress 0.2224298 
    ## Run 431 stress 0.1320258 
    ## ... Procrustes: rmse 1.379263e-06  max resid 2.787972e-06 
    ## ... Similar to previous best
    ## Run 432 stress 0.2600627 
    ## Run 433 stress 0.1802587 
    ## Run 434 stress 0.1407206 
    ## Run 435 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 4.150649e-07  max resid 7.572278e-07 
    ## ... Similar to previous best
    ## Run 436 stress 0.1383681 
    ## Run 437 stress 0.1320258 
    ## ... Procrustes: rmse 4.596773e-07  max resid 6.280919e-07 
    ## ... Similar to previous best
    ## Run 438 stress 0.2600627 
    ## Run 439 stress 0.1383681 
    ## Run 440 stress 0.1583016 
    ## Run 441 stress 0.1551863 
    ## Run 442 stress 0.1407298 
    ## Run 443 stress 0.1320258 
    ## ... Procrustes: rmse 2.554794e-06  max resid 5.237957e-06 
    ## ... Similar to previous best
    ## Run 444 stress 0.1320258 
    ## ... Procrustes: rmse 1.060118e-06  max resid 1.8903e-06 
    ## ... Similar to previous best
    ## Run 445 stress 0.1802752 
    ## Run 446 stress 0.1962476 
    ## Run 447 stress 0.1445811 
    ## Run 448 stress 0.1445811 
    ## Run 449 stress 0.1410448 
    ## Run 450 stress 0.1572668 
    ## Run 451 stress 0.3083098 
    ## Run 452 stress 0.2008163 
    ## Run 453 stress 0.1410449 
    ## Run 454 stress 0.1551863 
    ## Run 455 stress 0.1415299 
    ## Run 456 stress 0.1320258 
    ## ... Procrustes: rmse 8.858831e-07  max resid 1.78293e-06 
    ## ... Similar to previous best
    ## Run 457 stress 0.1551863 
    ## Run 458 stress 0.1320258 
    ## ... Procrustes: rmse 1.605846e-06  max resid 3.269041e-06 
    ## ... Similar to previous best
    ## Run 459 stress 0.1407207 
    ## Run 460 stress 0.1551863 
    ## Run 461 stress 0.1415299 
    ## Run 462 stress 0.1383681 
    ## Run 463 stress 0.1407298 
    ## Run 464 stress 0.1693717 
    ## Run 465 stress 0.1415299 
    ## Run 466 stress 0.2514078 
    ## Run 467 stress 0.1320258 
    ## ... Procrustes: rmse 1.836307e-06  max resid 3.773163e-06 
    ## ... Similar to previous best
    ## Run 468 stress 0.1320258 
    ## ... Procrustes: rmse 6.03313e-07  max resid 9.455794e-07 
    ## ... Similar to previous best
    ## Run 469 stress 0.1410452 
    ## Run 470 stress 0.1551863 
    ## Run 471 stress 0.1410449 
    ## Run 472 stress 0.1320258 
    ## ... Procrustes: rmse 2.684945e-06  max resid 3.358554e-06 
    ## ... Similar to previous best
    ## Run 473 stress 0.1320258 
    ## ... Procrustes: rmse 6.072795e-07  max resid 1.253921e-06 
    ## ... Similar to previous best
    ## Run 474 stress 0.1410449 
    ## Run 475 stress 0.2951814 
    ## Run 476 stress 0.1415299 
    ## Run 477 stress 0.2224298 
    ## Run 478 stress 0.2008161 
    ## Run 479 stress 0.1551863 
    ## Run 480 stress 0.1407208 
    ## Run 481 stress 0.1572668 
    ## Run 482 stress 0.1320258 
    ## ... Procrustes: rmse 9.266168e-07  max resid 1.77528e-06 
    ## ... Similar to previous best
    ## Run 483 stress 0.1410448 
    ## Run 484 stress 0.1415299 
    ## Run 485 stress 0.3083098 
    ## Run 486 stress 0.2711712 
    ## Run 487 stress 0.1551863 
    ## Run 488 stress 0.1415299 
    ## Run 489 stress 0.1320258 
    ## ... Procrustes: rmse 2.963889e-07  max resid 4.664831e-07 
    ## ... Similar to previous best
    ## Run 490 stress 0.1771969 
    ## Run 491 stress 0.1407298 
    ## Run 492 stress 0.2612821 
    ## Run 493 stress 0.2612821 
    ## Run 494 stress 0.1572668 
    ## Run 495 stress 0.1415299 
    ## Run 496 stress 0.1771969 
    ## Run 497 stress 0.1407298 
    ## Run 498 stress 0.1551863 
    ## Run 499 stress 0.1410449 
    ## Run 500 stress 0.141045 
    ## *** Best solution repeated 12 times

``` r
### Environmental
# Surveyed sites 
SD_beta_env_NMDS <- metaMDS(SD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07669178 
    ## Run 1 stress 0.08877307 
    ## Run 2 stress 0.08337114 
    ## Run 3 stress 0.07951495 
    ## Run 4 stress 0.08175707 
    ## Run 5 stress 0.08835185 
    ## Run 6 stress 0.0766915 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003129945  max resid 0.0005444043 
    ## ... Similar to previous best
    ## Run 7 stress 0.07669156 
    ## ... Procrustes: rmse 9.147349e-05  max resid 0.0001763011 
    ## ... Similar to previous best
    ## Run 8 stress 0.08459585 
    ## Run 9 stress 0.08182125 
    ## Run 10 stress 0.08063935 
    ## Run 11 stress 0.09070292 
    ## Run 12 stress 0.08265438 
    ## Run 13 stress 0.07669166 
    ## ... Procrustes: rmse 0.0002576301  max resid 0.0005012593 
    ## ... Similar to previous best
    ## Run 14 stress 0.08704848 
    ## Run 15 stress 0.08252153 
    ## Run 16 stress 0.08003182 
    ## Run 17 stress 0.0837573 
    ## Run 18 stress 0.08842187 
    ## Run 19 stress 0.08252117 
    ## Run 20 stress 0.08305281 
    ## Run 21 stress 0.07669172 
    ## ... Procrustes: rmse 0.0002890363  max resid 0.0005246463 
    ## ... Similar to previous best
    ## Run 22 stress 0.08842149 
    ## Run 23 stress 0.08433684 
    ## Run 24 stress 0.08189386 
    ## Run 25 stress 0.08243365 
    ## Run 26 stress 0.08175694 
    ## Run 27 stress 0.08274929 
    ## Run 28 stress 0.07948951 
    ## Run 29 stress 0.08433705 
    ## Run 30 stress 0.07949007 
    ## Run 31 stress 0.07669179 
    ## ... Procrustes: rmse 0.0003311396  max resid 0.0005667345 
    ## ... Similar to previous best
    ## Run 32 stress 0.08443336 
    ## Run 33 stress 0.08305372 
    ## Run 34 stress 0.08217355 
    ## Run 35 stress 0.08217343 
    ## Run 36 stress 0.08360263 
    ## Run 37 stress 0.07960667 
    ## Run 38 stress 0.08379934 
    ## Run 39 stress 0.08493055 
    ## Run 40 stress 0.07871747 
    ## Run 41 stress 0.0766915 
    ## ... Procrustes: rmse 6.80182e-05  max resid 0.0001203415 
    ## ... Similar to previous best
    ## Run 42 stress 0.08485162 
    ## Run 43 stress 0.07951471 
    ## Run 44 stress 0.08175695 
    ## Run 45 stress 0.08265452 
    ## Run 46 stress 0.08968607 
    ## Run 47 stress 0.08485155 
    ## Run 48 stress 0.0825011 
    ## Run 49 stress 0.08957596 
    ## Run 50 stress 0.08252136 
    ## Run 51 stress 0.08013823 
    ## Run 52 stress 0.08097776 
    ## Run 53 stress 0.08175694 
    ## Run 54 stress 0.07669155 
    ## ... Procrustes: rmse 8.341059e-05  max resid 0.0001504676 
    ## ... Similar to previous best
    ## Run 55 stress 0.08249608 
    ## Run 56 stress 0.09001684 
    ## Run 57 stress 0.08957542 
    ## Run 58 stress 0.08265457 
    ## Run 59 stress 0.07960684 
    ## Run 60 stress 0.07669157 
    ## ... Procrustes: rmse 8.59863e-05  max resid 0.0001950563 
    ## ... Similar to previous best
    ## Run 61 stress 0.07868331 
    ## Run 62 stress 0.08175689 
    ## Run 63 stress 0.08305109 
    ## Run 64 stress 0.08421031 
    ## Run 65 stress 0.08012959 
    ## Run 66 stress 0.08317313 
    ## Run 67 stress 0.08175694 
    ## Run 68 stress 0.07669154 
    ## ... Procrustes: rmse 0.0001270917  max resid 0.0002034454 
    ## ... Similar to previous best
    ## Run 69 stress 0.08265407 
    ## Run 70 stress 0.08003174 
    ## Run 71 stress 0.08355808 
    ## Run 72 stress 0.08249609 
    ## Run 73 stress 0.0850702 
    ## Run 74 stress 0.08175716 
    ## Run 75 stress 0.08175689 
    ## Run 76 stress 0.08104924 
    ## Run 77 stress 0.08433686 
    ## Run 78 stress 0.08250125 
    ## Run 79 stress 0.08305321 
    ## Run 80 stress 0.0787173 
    ## Run 81 stress 0.08314388 
    ## Run 82 stress 0.07868328 
    ## Run 83 stress 0.08217357 
    ## Run 84 stress 0.0824179 
    ## Run 85 stress 0.09091502 
    ## Run 86 stress 0.08175692 
    ## Run 87 stress 0.08265392 
    ## Run 88 stress 0.08249608 
    ## Run 89 stress 0.08835182 
    ## Run 90 stress 0.08527033 
    ## Run 91 stress 0.0849934 
    ## Run 92 stress 0.08182121 
    ## Run 93 stress 0.08375734 
    ## Run 94 stress 0.0825012 
    ## Run 95 stress 0.08097755 
    ## Run 96 stress 0.08241494 
    ## Run 97 stress 0.08477454 
    ## Run 98 stress 0.08104901 
    ## Run 99 stress 0.08104927 
    ## Run 100 stress 0.08453212 
    ## Run 101 stress 0.07669159 
    ## ... Procrustes: rmse 0.0002029704  max resid 0.000422746 
    ## ... Similar to previous best
    ## Run 102 stress 0.07949037 
    ## Run 103 stress 0.08957555 
    ## Run 104 stress 0.08329763 
    ## Run 105 stress 0.0826537 
    ## Run 106 stress 0.08097768 
    ## Run 107 stress 0.08314415 
    ## Run 108 stress 0.07669166 
    ## ... Procrustes: rmse 0.0002218623  max resid 0.0004217531 
    ## ... Similar to previous best
    ## Run 109 stress 0.08175707 
    ## Run 110 stress 0.0833716 
    ## Run 111 stress 0.07960669 
    ## Run 112 stress 0.08477421 
    ## Run 113 stress 0.08217344 
    ## Run 114 stress 0.08485118 
    ## Run 115 stress 0.08097788 
    ## Run 116 stress 0.08527029 
    ## Run 117 stress 0.0828359 
    ## Run 118 stress 0.08097773 
    ## Run 119 stress 0.08682648 
    ## Run 120 stress 0.08453197 
    ## Run 121 stress 0.08241789 
    ## Run 122 stress 0.08927487 
    ## Run 123 stress 0.08012965 
    ## Run 124 stress 0.08405262 
    ## Run 125 stress 0.08265496 
    ## Run 126 stress 0.0800322 
    ## Run 127 stress 0.08003181 
    ## Run 128 stress 0.08097775 
    ## Run 129 stress 0.08104904 
    ## Run 130 stress 0.08305496 
    ## Run 131 stress 0.07669164 
    ## ... Procrustes: rmse 0.0001634005  max resid 0.0002942162 
    ## ... Similar to previous best
    ## Run 132 stress 0.07669151 
    ## ... Procrustes: rmse 8.213394e-05  max resid 0.0001410561 
    ## ... Similar to previous best
    ## Run 133 stress 0.08355832 
    ## Run 134 stress 0.08182109 
    ## Run 135 stress 0.0831445 
    ## Run 136 stress 0.08175706 
    ## Run 137 stress 0.07871725 
    ## Run 138 stress 0.08835131 
    ## Run 139 stress 0.08013523 
    ## Run 140 stress 0.09028554 
    ## Run 141 stress 0.08684132 
    ## Run 142 stress 0.07669152 
    ## ... Procrustes: rmse 0.0001127933  max resid 0.000223295 
    ## ... Similar to previous best
    ## Run 143 stress 0.08182131 
    ## Run 144 stress 0.08012947 
    ## Run 145 stress 0.0817569 
    ## Run 146 stress 0.08898817 
    ## Run 147 stress 0.08927447 
    ## Run 148 stress 0.08421018 
    ## Run 149 stress 0.08144832 
    ## Run 150 stress 0.08063931 
    ## Run 151 stress 0.08175696 
    ## Run 152 stress 0.08252787 
    ## Run 153 stress 0.07948971 
    ## Run 154 stress 0.07868345 
    ## Run 155 stress 0.08097757 
    ## Run 156 stress 0.08485126 
    ## Run 157 stress 0.08003173 
    ## Run 158 stress 0.08486998 
    ## Run 159 stress 0.08360271 
    ## Run 160 stress 0.08314383 
    ## Run 161 stress 0.07669154 
    ## ... Procrustes: rmse 0.0001538049  max resid 0.000341435 
    ## ... Similar to previous best
    ## Run 162 stress 0.07871748 
    ## Run 163 stress 0.08421018 
    ## Run 164 stress 0.08876554 
    ## Run 165 stress 0.08104914 
    ## Run 166 stress 0.09243745 
    ## Run 167 stress 0.08530157 
    ## Run 168 stress 0.08929133 
    ## Run 169 stress 0.08265404 
    ## Run 170 stress 0.0787175 
    ## Run 171 stress 0.07669157 
    ## ... Procrustes: rmse 0.0001697018  max resid 0.0003187214 
    ## ... Similar to previous best
    ## Run 172 stress 0.08507064 
    ## Run 173 stress 0.08514254 
    ## Run 174 stress 0.0825212 
    ## Run 175 stress 0.0831443 
    ## Run 176 stress 0.08717155 
    ## Run 177 stress 0.07871763 
    ## Run 178 stress 0.08175689 
    ## Run 179 stress 0.08493134 
    ## Run 180 stress 0.0817569 
    ## Run 181 stress 0.09063309 
    ## Run 182 stress 0.07960675 
    ## Run 183 stress 0.08957505 
    ## Run 184 stress 0.08615932 
    ## Run 185 stress 0.08493065 
    ## Run 186 stress 0.07871754 
    ## Run 187 stress 0.08305383 
    ## Run 188 stress 0.08477442 
    ## Run 189 stress 0.08684209 
    ## Run 190 stress 0.08314425 
    ## Run 191 stress 0.08477446 
    ## Run 192 stress 0.08553053 
    ## Run 193 stress 0.08144959 
    ## Run 194 stress 0.08012959 
    ## Run 195 stress 0.07669174 
    ## ... Procrustes: rmse 0.000292574  max resid 0.0005095758 
    ## ... Similar to previous best
    ## Run 196 stress 0.07960679 
    ## Run 197 stress 0.08175703 
    ## Run 198 stress 0.08957546 
    ## Run 199 stress 0.08252815 
    ## Run 200 stress 0.08217339 
    ## Run 201 stress 0.08097758 
    ## Run 202 stress 0.08929051 
    ## Run 203 stress 0.08362754 
    ## Run 204 stress 0.08556626 
    ## Run 205 stress 0.07871749 
    ## Run 206 stress 0.08175699 
    ## Run 207 stress 0.07960659 
    ## Run 208 stress 0.07868385 
    ## Run 209 stress 0.08012974 
    ## Run 210 stress 0.0766915 
    ## ... New best solution
    ## ... Procrustes: rmse 3.978064e-05  max resid 6.548462e-05 
    ## ... Similar to previous best
    ## Run 211 stress 0.08265426 
    ## Run 212 stress 0.08012939 
    ## Run 213 stress 0.08175721 
    ## Run 214 stress 0.07669173 
    ## ... Procrustes: rmse 0.0003066252  max resid 0.0005258476 
    ## ... Similar to previous best
    ## Run 215 stress 0.09103929 
    ## Run 216 stress 0.08012988 
    ## Run 217 stress 0.08104924 
    ## Run 218 stress 0.0825214 
    ## Run 219 stress 0.08097777 
    ## Run 220 stress 0.08265428 
    ## Run 221 stress 0.08175716 
    ## Run 222 stress 0.08217347 
    ## Run 223 stress 0.0809777 
    ## Run 224 stress 0.07669166 
    ## ... Procrustes: rmse 0.0002282595  max resid 0.0003977994 
    ## ... Similar to previous best
    ## Run 225 stress 0.08252118 
    ## Run 226 stress 0.08547293 
    ## Run 227 stress 0.08012984 
    ## Run 228 stress 0.08175712 
    ## Run 229 stress 0.08682751 
    ## Run 230 stress 0.07669156 
    ## ... Procrustes: rmse 0.0001727769  max resid 0.0002959146 
    ## ... Similar to previous best
    ## Run 231 stress 0.07731519 
    ## Run 232 stress 0.08375706 
    ## Run 233 stress 0.08175699 
    ## Run 234 stress 0.08249603 
    ## Run 235 stress 0.08097768 
    ## Run 236 stress 0.07669179 
    ## ... Procrustes: rmse 0.000311266  max resid 0.0005520507 
    ## ... Similar to previous best
    ## Run 237 stress 0.08252119 
    ## Run 238 stress 0.08144905 
    ## Run 239 stress 0.0831734 
    ## Run 240 stress 0.0794901 
    ## Run 241 stress 0.08252142 
    ## Run 242 stress 0.08012969 
    ## Run 243 stress 0.0766917 
    ## ... Procrustes: rmse 0.000282134  max resid 0.0005051555 
    ## ... Similar to previous best
    ## Run 244 stress 0.07669163 
    ## ... Procrustes: rmse 0.0002316918  max resid 0.0003983574 
    ## ... Similar to previous best
    ## Run 245 stress 0.08672773 
    ## Run 246 stress 0.08104886 
    ## Run 247 stress 0.0787174 
    ## Run 248 stress 0.07669169 
    ## ... Procrustes: rmse 0.0002488789  max resid 0.0005366758 
    ## ... Similar to previous best
    ## Run 249 stress 0.08506931 
    ## Run 250 stress 0.08433709 
    ## Run 251 stress 0.08464781 
    ## Run 252 stress 0.07669173 
    ## ... Procrustes: rmse 0.0003036923  max resid 0.0005238629 
    ## ... Similar to previous best
    ## Run 253 stress 0.07960671 
    ## Run 254 stress 0.08360263 
    ## Run 255 stress 0.08265446 
    ## Run 256 stress 0.0787173 
    ## Run 257 stress 0.08175708 
    ## Run 258 stress 0.08243374 
    ## Run 259 stress 0.08305437 
    ## Run 260 stress 0.08252134 
    ## Run 261 stress 0.08175708 
    ## Run 262 stress 0.08283575 
    ## Run 263 stress 0.08243369 
    ## Run 264 stress 0.07871727 
    ## Run 265 stress 0.08418563 
    ## Run 266 stress 0.07960695 
    ## Run 267 stress 0.08063933 
    ## Run 268 stress 0.08189374 
    ## Run 269 stress 0.07960663 
    ## Run 270 stress 0.08063945 
    ## Run 271 stress 0.08433685 
    ## Run 272 stress 0.07951497 
    ## Run 273 stress 0.08175716 
    ## Run 274 stress 0.08547312 
    ## Run 275 stress 0.07871749 
    ## Run 276 stress 0.08433686 
    ## Run 277 stress 0.08421008 
    ## Run 278 stress 0.08265402 
    ## Run 279 stress 0.08684091 
    ## Run 280 stress 0.07949038 
    ## Run 281 stress 0.07731522 
    ## Run 282 stress 0.08443224 
    ## Run 283 stress 0.07943343 
    ## Run 284 stress 0.08265371 
    ## Run 285 stress 0.08097769 
    ## Run 286 stress 0.08104914 
    ## Run 287 stress 0.08104915 
    ## Run 288 stress 0.08252127 
    ## Run 289 stress 0.08243364 
    ## Run 290 stress 0.08175696 
    ## Run 291 stress 0.08175706 
    ## Run 292 stress 0.08063952 
    ## Run 293 stress 0.08433684 
    ## Run 294 stress 0.08003186 
    ## Run 295 stress 0.08175708 
    ## Run 296 stress 0.07949013 
    ## Run 297 stress 0.08104942 
    ## Run 298 stress 0.07948952 
    ## Run 299 stress 0.08210327 
    ## Run 300 stress 0.07868338 
    ## Run 301 stress 0.08898788 
    ## Run 302 stress 0.08104903 
    ## Run 303 stress 0.08433688 
    ## Run 304 stress 0.07871732 
    ## Run 305 stress 0.07669151 
    ## ... Procrustes: rmse 5.796656e-05  max resid 9.129115e-05 
    ## ... Similar to previous best
    ## Run 306 stress 0.08265355 
    ## Run 307 stress 0.07960669 
    ## Run 308 stress 0.08250149 
    ## Run 309 stress 0.08499314 
    ## Run 310 stress 0.08097757 
    ## Run 311 stress 0.08510822 
    ## Run 312 stress 0.08547277 
    ## Run 313 stress 0.07871753 
    ## Run 314 stress 0.08189397 
    ## Run 315 stress 0.0852702 
    ## Run 316 stress 0.08876472 
    ## Run 317 stress 0.08249601 
    ## Run 318 stress 0.0830542 
    ## Run 319 stress 0.0850703 
    ## Run 320 stress 0.08175691 
    ## Run 321 stress 0.08175719 
    ## Run 322 stress 0.07731526 
    ## Run 323 stress 0.0824339 
    ## Run 324 stress 0.08684314 
    ## Run 325 stress 0.08355837 
    ## Run 326 stress 0.07871748 
    ## Run 327 stress 0.08375712 
    ## Run 328 stress 0.0830551 
    ## Run 329 stress 0.08104896 
    ## Run 330 stress 0.08370327 
    ## Run 331 stress 0.08144932 
    ## Run 332 stress 0.08243379 
    ## Run 333 stress 0.08175689 
    ## Run 334 stress 0.07669156 
    ## ... Procrustes: rmse 0.0001041564  max resid 0.0002121108 
    ## ... Similar to previous best
    ## Run 335 stress 0.08360262 
    ## Run 336 stress 0.08252147 
    ## Run 337 stress 0.08012975 
    ## Run 338 stress 0.07731518 
    ## Run 339 stress 0.08003226 
    ## Run 340 stress 0.07868362 
    ## Run 341 stress 0.08500321 
    ## Run 342 stress 0.07669167 
    ## ... Procrustes: rmse 0.0002295484  max resid 0.0003979641 
    ## ... Similar to previous best
    ## Run 343 stress 0.08418557 
    ## Run 344 stress 0.08104912 
    ## Run 345 stress 0.08684153 
    ## Run 346 stress 0.08012946 
    ## Run 347 stress 0.08553075 
    ## Run 348 stress 0.08175699 
    ## Run 349 stress 0.07669161 
    ## ... Procrustes: rmse 0.0002178621  max resid 0.0003721931 
    ## ... Similar to previous best
    ## Run 350 stress 0.08527044 
    ## Run 351 stress 0.08375743 
    ## Run 352 stress 0.08305498 
    ## Run 353 stress 0.07960677 
    ## Run 354 stress 0.07871745 
    ## Run 355 stress 0.08252126 
    ## Run 356 stress 0.08418586 
    ## Run 357 stress 0.08317353 
    ## Run 358 stress 0.08175703 
    ## Run 359 stress 0.07669168 
    ## ... Procrustes: rmse 0.0002406296  max resid 0.0004261252 
    ## ... Similar to previous best
    ## Run 360 stress 0.08305358 
    ## Run 361 stress 0.09012166 
    ## Run 362 stress 0.08421017 
    ## Run 363 stress 0.08355809 
    ## Run 364 stress 0.08003201 
    ## Run 365 stress 0.08527016 
    ## Run 366 stress 0.07669158 
    ## ... Procrustes: rmse 0.0001884522  max resid 0.0003417876 
    ## ... Similar to previous best
    ## Run 367 stress 0.08063965 
    ## Run 368 stress 0.08362758 
    ## Run 369 stress 0.07949019 
    ## Run 370 stress 0.08842213 
    ## Run 371 stress 0.08175706 
    ## Run 372 stress 0.08097741 
    ## Run 373 stress 0.07960665 
    ## Run 374 stress 0.07951491 
    ## Run 375 stress 0.08443353 
    ## Run 376 stress 0.07943389 
    ## Run 377 stress 0.08252148 
    ## Run 378 stress 0.08314422 
    ## Run 379 stress 0.08421015 
    ## Run 380 stress 0.08189373 
    ## Run 381 stress 0.0766916 
    ## ... Procrustes: rmse 0.0001767694  max resid 0.0002998353 
    ## ... Similar to previous best
    ## Run 382 stress 0.08013735 
    ## Run 383 stress 0.08957635 
    ## Run 384 stress 0.08012956 
    ## Run 385 stress 0.08104901 
    ## Run 386 stress 0.07731525 
    ## Run 387 stress 0.07868353 
    ## Run 388 stress 0.08217345 
    ## Run 389 stress 0.08249616 
    ## Run 390 stress 0.07871745 
    ## Run 391 stress 0.08265315 
    ## Run 392 stress 0.07871743 
    ## Run 393 stress 0.07948998 
    ## Run 394 stress 0.08486851 
    ## Run 395 stress 0.0833711 
    ## Run 396 stress 0.07871757 
    ## Run 397 stress 0.08104893 
    ## Run 398 stress 0.08360257 
    ## Run 399 stress 0.08305334 
    ## Run 400 stress 0.08360275 
    ## Run 401 stress 0.08314422 
    ## Run 402 stress 0.0818211 
    ## Run 403 stress 0.08097762 
    ## Run 404 stress 0.08305318 
    ## Run 405 stress 0.07868376 
    ## Run 406 stress 0.08888032 
    ## Run 407 stress 0.09196136 
    ## Run 408 stress 0.08175704 
    ## Run 409 stress 0.07669158 
    ## ... Procrustes: rmse 0.0001797906  max resid 0.0003091177 
    ## ... Similar to previous best
    ## Run 410 stress 0.08305478 
    ## Run 411 stress 0.07731533 
    ## Run 412 stress 0.08175703 
    ## Run 413 stress 0.08252134 
    ## Run 414 stress 0.08252116 
    ## Run 415 stress 0.08189384 
    ## Run 416 stress 0.08583944 
    ## Run 417 stress 0.07949005 
    ## Run 418 stress 0.0821852 
    ## Run 419 stress 0.08175694 
    ## Run 420 stress 0.07669168 
    ## ... Procrustes: rmse 0.0002685961  max resid 0.0004610657 
    ## ... Similar to previous best
    ## Run 421 stress 0.07669161 
    ## ... Procrustes: rmse 0.0002199179  max resid 0.000377384 
    ## ... Similar to previous best
    ## Run 422 stress 0.08957575 
    ## Run 423 stress 0.08360262 
    ## Run 424 stress 0.08175695 
    ## Run 425 stress 0.08104926 
    ## Run 426 stress 0.07731536 
    ## Run 427 stress 0.08421011 
    ## Run 428 stress 0.07868353 
    ## Run 429 stress 0.08684206 
    ## Run 430 stress 0.08063939 
    ## Run 431 stress 0.07669169 
    ## ... Procrustes: rmse 0.0002535972  max resid 0.0004451577 
    ## ... Similar to previous best
    ## Run 432 stress 0.07871754 
    ## Run 433 stress 0.08250107 
    ## Run 434 stress 0.08249597 
    ## Run 435 stress 0.08252158 
    ## Run 436 stress 0.07871739 
    ## Run 437 stress 0.07943394 
    ## Run 438 stress 0.08314388 
    ## Run 439 stress 0.08250191 
    ## Run 440 stress 0.0824149 
    ## Run 441 stress 0.07868342 
    ## Run 442 stress 0.08305348 
    ## Run 443 stress 0.08283558 
    ## Run 444 stress 0.08175693 
    ## Run 445 stress 0.08003192 
    ## Run 446 stress 0.07868341 
    ## Run 447 stress 0.08339579 
    ## Run 448 stress 0.08189402 
    ## Run 449 stress 0.08144929 
    ## Run 450 stress 0.08283567 
    ## Run 451 stress 0.08317319 
    ## Run 452 stress 0.08355802 
    ## Run 453 stress 0.08012953 
    ## Run 454 stress 0.07669152 
    ## ... Procrustes: rmse 0.0001069811  max resid 0.0002158412 
    ## ... Similar to previous best
    ## Run 455 stress 0.08175711 
    ## Run 456 stress 0.08274924 
    ## Run 457 stress 0.08217348 
    ## Run 458 stress 0.0826538 
    ## Run 459 stress 0.08012942 
    ## Run 460 stress 0.08175694 
    ## Run 461 stress 0.08527037 
    ## Run 462 stress 0.07731535 
    ## Run 463 stress 0.08182117 
    ## Run 464 stress 0.08421018 
    ## Run 465 stress 0.08012953 
    ## Run 466 stress 0.08249626 
    ## Run 467 stress 0.08355867 
    ## Run 468 stress 0.0809774 
    ## Run 469 stress 0.08012959 
    ## Run 470 stress 0.08337166 
    ## Run 471 stress 0.08314413 
    ## Run 472 stress 0.08003214 
    ## Run 473 stress 0.08274923 
    ## Run 474 stress 0.08556639 
    ## Run 475 stress 0.08530175 
    ## Run 476 stress 0.08673788 
    ## Run 477 stress 0.0786836 
    ## Run 478 stress 0.08673752 
    ## Run 479 stress 0.08305338 
    ## Run 480 stress 0.07871727 
    ## Run 481 stress 0.0883514 
    ## Run 482 stress 0.08583951 
    ## Run 483 stress 0.08317351 
    ## Run 484 stress 0.08104902 
    ## Run 485 stress 0.08842083 
    ## Run 486 stress 0.08433698 
    ## Run 487 stress 0.08405257 
    ## Run 488 stress 0.08842128 
    ## Run 489 stress 0.08314417 
    ## Run 490 stress 0.08144938 
    ## Run 491 stress 0.08375743 
    ## Run 492 stress 0.08249618 
    ## Run 493 stress 0.08360272 
    ## Run 494 stress 0.08012995 
    ## Run 495 stress 0.08337167 
    ## Run 496 stress 0.08217351 
    ## Run 497 stress 0.08175705 
    ## Run 498 stress 0.08835134 
    ## Run 499 stress 0.08314414 
    ## Run 500 stress 0.08433714 
    ## *** Best solution repeated 21 times

``` r
# Mixed and stratified lakes
SD_beta_env_MS_NMDS <- metaMDS(SD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.08973862 
    ## Run 2 stress 0.09030414 
    ## Run 3 stress 0.09030419 
    ## Run 4 stress 0.09159095 
    ## Run 5 stress 0.0938185 
    ## Run 6 stress 0.0940798 
    ## Run 7 stress 0.08440266 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01443363  max resid 0.04307049 
    ## Run 8 stress 0.09145329 
    ## Run 9 stress 0.09308946 
    ## Run 10 stress 0.09268327 
    ## Run 11 stress 0.09168936 
    ## Run 12 stress 0.09416155 
    ## Run 13 stress 0.09030405 
    ## Run 14 stress 0.09159084 
    ## Run 15 stress 0.1026215 
    ## Run 16 stress 0.09407972 
    ## Run 17 stress 0.09286089 
    ## Run 18 stress 0.09308959 
    ## Run 19 stress 0.1038049 
    ## Run 20 stress 0.09535447 
    ## Run 21 stress 0.0844027 
    ## ... Procrustes: rmse 3.437305e-05  max resid 6.905761e-05 
    ## ... Similar to previous best
    ## Run 22 stress 0.09464466 
    ## Run 23 stress 0.09030404 
    ## Run 24 stress 0.2643456 
    ## Run 25 stress 0.09030407 
    ## Run 26 stress 0.08440261 
    ## ... New best solution
    ## ... Procrustes: rmse 8.861405e-05  max resid 0.0002436909 
    ## ... Similar to previous best
    ## Run 27 stress 0.09407411 
    ## Run 28 stress 0.09030404 
    ## Run 29 stress 0.09969972 
    ## Run 30 stress 0.0937419 
    ## Run 31 stress 0.09407979 
    ## Run 32 stress 0.08973862 
    ## Run 33 stress 0.08503468 
    ## Run 34 stress 0.09286083 
    ## Run 35 stress 0.08773481 
    ## Run 36 stress 0.09268324 
    ## Run 37 stress 0.09464464 
    ## Run 38 stress 0.08773469 
    ## Run 39 stress 0.09337261 
    ## Run 40 stress 0.08973862 
    ## Run 41 stress 0.09030395 
    ## Run 42 stress 0.09407456 
    ## Run 43 stress 0.08503616 
    ## Run 44 stress 0.09268381 
    ## Run 45 stress 0.09159094 
    ## Run 46 stress 0.08503533 
    ## Run 47 stress 0.097212 
    ## Run 48 stress 0.08503518 
    ## Run 49 stress 0.09337305 
    ## Run 50 stress 0.09168966 
    ## Run 51 stress 0.1026217 
    ## Run 52 stress 0.09447312 
    ## Run 53 stress 0.09374254 
    ## Run 54 stress 0.09030414 
    ## Run 55 stress 0.085035 
    ## Run 56 stress 0.08773465 
    ## Run 57 stress 0.09669859 
    ## Run 58 stress 0.09407987 
    ## Run 59 stress 0.08440263 
    ## ... Procrustes: rmse 7.564596e-05  max resid 0.0002140094 
    ## ... Similar to previous best
    ## Run 60 stress 0.08440258 
    ## ... New best solution
    ## ... Procrustes: rmse 6.290199e-05  max resid 0.0001322527 
    ## ... Similar to previous best
    ## Run 61 stress 0.09168932 
    ## Run 62 stress 0.09464422 
    ## Run 63 stress 0.08503499 
    ## Run 64 stress 0.09030411 
    ## Run 65 stress 0.09400524 
    ## Run 66 stress 0.09539209 
    ## Run 67 stress 0.09168946 
    ## Run 68 stress 0.09145322 
    ## Run 69 stress 0.09145332 
    ## Run 70 stress 0.09374405 
    ## Run 71 stress 0.09760842 
    ## Run 72 stress 0.09268322 
    ## Run 73 stress 0.0897387 
    ## Run 74 stress 0.0928609 
    ## Run 75 stress 0.08503476 
    ## Run 76 stress 0.09374342 
    ## Run 77 stress 0.09030396 
    ## Run 78 stress 0.08503683 
    ## Run 79 stress 0.09030408 
    ## Run 80 stress 0.08973878 
    ## Run 81 stress 0.09337289 
    ## Run 82 stress 0.09268339 
    ## Run 83 stress 0.09669867 
    ## Run 84 stress 0.09407996 
    ## Run 85 stress 0.09159095 
    ## Run 86 stress 0.09465904 
    ## Run 87 stress 0.09380599 
    ## Run 88 stress 0.1026218 
    ## Run 89 stress 0.09308958 
    ## Run 90 stress 0.09286092 
    ## Run 91 stress 0.08503467 
    ## Run 92 stress 0.09465905 
    ## Run 93 stress 0.0914533 
    ## Run 94 stress 0.0940714 
    ## Run 95 stress 0.09159092 
    ## Run 96 stress 0.09969992 
    ## Run 97 stress 0.09407996 
    ## Run 98 stress 0.09760819 
    ## Run 99 stress 0.08973867 
    ## Run 100 stress 0.09416132 
    ## Run 101 stress 0.09464417 
    ## Run 102 stress 0.09445504 
    ## Run 103 stress 0.08773468 
    ## Run 104 stress 0.09030395 
    ## Run 105 stress 0.09465907 
    ## Run 106 stress 0.08973873 
    ## Run 107 stress 0.09337283 
    ## Run 108 stress 0.09416143 
    ## Run 109 stress 0.09407984 
    ## Run 110 stress 0.09407977 
    ## Run 111 stress 0.09286085 
    ## Run 112 stress 0.09407997 
    ## Run 113 stress 0.08503474 
    ## Run 114 stress 0.08773481 
    ## Run 115 stress 0.3315595 
    ## Run 116 stress 0.1038052 
    ## Run 117 stress 0.09407974 
    ## Run 118 stress 0.08503481 
    ## Run 119 stress 0.09308954 
    ## Run 120 stress 0.094455 
    ## Run 121 stress 0.09416134 
    ## Run 122 stress 0.09408017 
    ## Run 123 stress 0.09407979 
    ## Run 124 stress 0.1052848 
    ## Run 125 stress 0.0844026 
    ## ... Procrustes: rmse 7.964058e-05  max resid 0.0001566431 
    ## ... Similar to previous best
    ## Run 126 stress 0.08503481 
    ## Run 127 stress 0.09159092 
    ## Run 128 stress 0.09447304 
    ## Run 129 stress 0.09268341 
    ## Run 130 stress 0.09145321 
    ## Run 131 stress 0.09268327 
    ## Run 132 stress 0.09969977 
    ## Run 133 stress 0.1026215 
    ## Run 134 stress 0.09308948 
    ## Run 135 stress 0.0850362 
    ## Run 136 stress 0.09145325 
    ## Run 137 stress 0.09337294 
    ## Run 138 stress 0.08773483 
    ## Run 139 stress 0.09321661 
    ## Run 140 stress 0.09416131 
    ## Run 141 stress 0.0850346 
    ## Run 142 stress 0.1097011 
    ## Run 143 stress 0.09969983 
    ## Run 144 stress 0.09168954 
    ## Run 145 stress 0.08503501 
    ## Run 146 stress 0.09969976 
    ## Run 147 stress 0.09286083 
    ## Run 148 stress 0.09268357 
    ## Run 149 stress 0.10196 
    ## Run 150 stress 0.09308947 
    ## Run 151 stress 0.09590485 
    ## Run 152 stress 0.09416143 
    ## Run 153 stress 0.08503503 
    ## Run 154 stress 0.09539206 
    ## Run 155 stress 0.09669889 
    ## Run 156 stress 0.09969984 
    ## Run 157 stress 0.08440266 
    ## ... Procrustes: rmse 7.13537e-05  max resid 0.0001288165 
    ## ... Similar to previous best
    ## Run 158 stress 0.08503492 
    ## Run 159 stress 0.09380582 
    ## Run 160 stress 0.08440253 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002397181  max resid 0.000442016 
    ## ... Similar to previous best
    ## Run 161 stress 0.09337272 
    ## Run 162 stress 0.09168928 
    ## Run 163 stress 0.09308973 
    ## Run 164 stress 0.09268332 
    ## Run 165 stress 0.09465906 
    ## Run 166 stress 0.0877347 
    ## Run 167 stress 0.09337292 
    ## Run 168 stress 0.09159089 
    ## Run 169 stress 0.09030414 
    ## Run 170 stress 0.08773476 
    ## Run 171 stress 0.08503627 
    ## Run 172 stress 0.09416143 
    ## Run 173 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001055188  max resid 0.0002588295 
    ## ... Similar to previous best
    ## Run 174 stress 0.09380575 
    ## Run 175 stress 0.0938058 
    ## Run 176 stress 0.09286085 
    ## Run 177 stress 0.09381854 
    ## Run 178 stress 0.09465904 
    ## Run 179 stress 0.09400537 
    ## Run 180 stress 0.08440273 
    ## ... Procrustes: rmse 0.0003620957  max resid 0.0006520358 
    ## ... Similar to previous best
    ## Run 181 stress 0.09268345 
    ## Run 182 stress 0.1038034 
    ## Run 183 stress 0.09407984 
    ## Run 184 stress 0.08503486 
    ## Run 185 stress 0.093373 
    ## Run 186 stress 0.08440256 
    ## ... Procrustes: rmse 0.0002187814  max resid 0.000397911 
    ## ... Similar to previous best
    ## Run 187 stress 0.09535497 
    ## Run 188 stress 0.09159087 
    ## Run 189 stress 0.09407998 
    ## Run 190 stress 0.09969983 
    ## Run 191 stress 0.09464483 
    ## Run 192 stress 0.09308954 
    ## Run 193 stress 0.08503484 
    ## Run 194 stress 0.09268359 
    ## Run 195 stress 0.0916895 
    ## Run 196 stress 0.0926837 
    ## Run 197 stress 0.09159083 
    ## Run 198 stress 0.09268363 
    ## Run 199 stress 0.09535448 
    ## Run 200 stress 0.09465904 
    ## Run 201 stress 0.08503621 
    ## Run 202 stress 0.09407103 
    ## Run 203 stress 0.09030393 
    ## Run 204 stress 0.09445499 
    ## Run 205 stress 0.09374271 
    ## Run 206 stress 0.09145322 
    ## Run 207 stress 0.09030402 
    ## Run 208 stress 0.08503464 
    ## Run 209 stress 0.09159086 
    ## Run 210 stress 0.09590906 
    ## Run 211 stress 0.09159102 
    ## Run 212 stress 0.09159088 
    ## Run 213 stress 0.08440254 
    ## ... Procrustes: rmse 3.820531e-05  max resid 8.016045e-05 
    ## ... Similar to previous best
    ## Run 214 stress 0.09145324 
    ## Run 215 stress 0.08503468 
    ## Run 216 stress 0.09464429 
    ## Run 217 stress 0.09337275 
    ## Run 218 stress 0.09407992 
    ## Run 219 stress 0.09465905 
    ## Run 220 stress 0.09416136 
    ## Run 221 stress 0.1038049 
    ## Run 222 stress 0.09407988 
    ## Run 223 stress 0.09145326 
    ## Run 224 stress 0.08773468 
    ## Run 225 stress 0.08973883 
    ## Run 226 stress 0.09168939 
    ## Run 227 stress 0.09286084 
    ## Run 228 stress 0.09268328 
    ## Run 229 stress 0.09145326 
    ## Run 230 stress 0.09721196 
    ## Run 231 stress 0.09337248 
    ## Run 232 stress 0.09407965 
    ## Run 233 stress 0.08773482 
    ## Run 234 stress 0.09308973 
    ## Run 235 stress 0.08973881 
    ## Run 236 stress 0.09337271 
    ## Run 237 stress 0.08773483 
    ## Run 238 stress 0.0916895 
    ## Run 239 stress 0.09030403 
    ## Run 240 stress 0.09464449 
    ## Run 241 stress 0.09407417 
    ## Run 242 stress 0.0850364 
    ## Run 243 stress 0.09590911 
    ## Run 244 stress 0.09403413 
    ## Run 245 stress 0.08773466 
    ## Run 246 stress 0.08440263 
    ## ... Procrustes: rmse 0.0002738464  max resid 0.0005573499 
    ## ... Similar to previous best
    ## Run 247 stress 0.09168928 
    ## Run 248 stress 0.09286083 
    ## Run 249 stress 0.09145322 
    ## Run 250 stress 0.09539201 
    ## Run 251 stress 0.09030397 
    ## Run 252 stress 0.08503583 
    ## Run 253 stress 0.08773494 
    ## Run 254 stress 0.09145326 
    ## Run 255 stress 0.09268364 
    ## Run 256 stress 0.09969976 
    ## Run 257 stress 0.09407385 
    ## Run 258 stress 0.08503514 
    ## Run 259 stress 0.09337257 
    ## Run 260 stress 0.08440256 
    ## ... Procrustes: rmse 0.0002007861  max resid 0.0003689659 
    ## ... Similar to previous best
    ## Run 261 stress 0.08440263 
    ## ... Procrustes: rmse 0.000277428  max resid 0.0005040673 
    ## ... Similar to previous best
    ## Run 262 stress 0.09407392 
    ## Run 263 stress 0.09337276 
    ## Run 264 stress 0.09590871 
    ## Run 265 stress 0.09030402 
    ## Run 266 stress 0.0940711 
    ## Run 267 stress 0.100461 
    ## Run 268 stress 0.08503472 
    ## Run 269 stress 0.09408008 
    ## Run 270 stress 0.09030396 
    ## Run 271 stress 0.08503706 
    ## Run 272 stress 0.08773469 
    ## Run 273 stress 0.08440255 
    ## ... Procrustes: rmse 5.040325e-05  max resid 0.0001598678 
    ## ... Similar to previous best
    ## Run 274 stress 0.09407996 
    ## Run 275 stress 0.09268398 
    ## Run 276 stress 0.1038045 
    ## Run 277 stress 0.09030399 
    ## Run 278 stress 0.09268355 
    ## Run 279 stress 0.08440258 
    ## ... Procrustes: rmse 0.0002114195  max resid 0.0004341329 
    ## ... Similar to previous best
    ## Run 280 stress 0.08503475 
    ## Run 281 stress 0.09400517 
    ## Run 282 stress 0.08503627 
    ## Run 283 stress 0.09407975 
    ## Run 284 stress 0.09268386 
    ## Run 285 stress 0.08440253 
    ## ... Procrustes: rmse 0.0001537799  max resid 0.0003146199 
    ## ... Similar to previous best
    ## Run 286 stress 0.08973888 
    ## Run 287 stress 0.09447305 
    ## Run 288 stress 0.09030403 
    ## Run 289 stress 0.09308978 
    ## Run 290 stress 0.08773472 
    ## Run 291 stress 0.09969992 
    ## Run 292 stress 0.1038046 
    ## Run 293 stress 0.09145321 
    ## Run 294 stress 0.08773477 
    ## Run 295 stress 0.09286086 
    ## Run 296 stress 0.09374352 
    ## Run 297 stress 0.08503492 
    ## Run 298 stress 0.09268369 
    ## Run 299 stress 0.09407964 
    ## Run 300 stress 0.09760826 
    ## Run 301 stress 0.09969965 
    ## Run 302 stress 0.0850365 
    ## Run 303 stress 0.08773476 
    ## Run 304 stress 0.08773469 
    ## Run 305 stress 0.09374309 
    ## Run 306 stress 0.08973884 
    ## Run 307 stress 0.0941614 
    ## Run 308 stress 0.09535493 
    ## Run 309 stress 0.08773466 
    ## Run 310 stress 0.09030396 
    ## Run 311 stress 0.09030393 
    ## Run 312 stress 0.09407987 
    ## Run 313 stress 0.0928609 
    ## Run 314 stress 0.09286083 
    ## Run 315 stress 0.09374215 
    ## Run 316 stress 0.09535495 
    ## Run 317 stress 0.09969981 
    ## Run 318 stress 0.09445486 
    ## Run 319 stress 0.09337298 
    ## Run 320 stress 0.09268345 
    ## Run 321 stress 0.09380576 
    ## Run 322 stress 0.3245465 
    ## Run 323 stress 0.08503479 
    ## Run 324 stress 0.09159084 
    ## Run 325 stress 0.09539204 
    ## Run 326 stress 0.09030397 
    ## Run 327 stress 0.09400516 
    ## Run 328 stress 0.09407975 
    ## Run 329 stress 0.09416147 
    ## Run 330 stress 0.08773465 
    ## Run 331 stress 0.090304 
    ## Run 332 stress 0.09447322 
    ## Run 333 stress 0.09286083 
    ## Run 334 stress 0.09030404 
    ## Run 335 stress 0.09407983 
    ## Run 336 stress 0.09168933 
    ## Run 337 stress 0.09539195 
    ## Run 338 stress 0.08973865 
    ## Run 339 stress 0.09286084 
    ## Run 340 stress 0.09407978 
    ## Run 341 stress 0.09416131 
    ## Run 342 stress 0.0953547 
    ## Run 343 stress 0.09407989 
    ## Run 344 stress 0.09268367 
    ## Run 345 stress 0.09407974 
    ## Run 346 stress 0.09159098 
    ## Run 347 stress 0.09977532 
    ## Run 348 stress 0.08503472 
    ## Run 349 stress 0.09030395 
    ## Run 350 stress 0.1038047 
    ## Run 351 stress 0.09268318 
    ## Run 352 stress 0.09374297 
    ## Run 353 stress 0.08503503 
    ## Run 354 stress 0.09969965 
    ## Run 355 stress 0.09380577 
    ## Run 356 stress 0.09145321 
    ## Run 357 stress 0.09337241 
    ## Run 358 stress 0.1038084 
    ## Run 359 stress 0.08503491 
    ## Run 360 stress 0.09408009 
    ## Run 361 stress 0.09168943 
    ## Run 362 stress 0.09168926 
    ## Run 363 stress 0.09400542 
    ## Run 364 stress 0.08973862 
    ## Run 365 stress 0.08973866 
    ## Run 366 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 9.516891e-05  max resid 0.000179812 
    ## ... Similar to previous best
    ## Run 367 stress 0.08973875 
    ## Run 368 stress 0.09337309 
    ## Run 369 stress 0.09407968 
    ## Run 370 stress 0.1038037 
    ## Run 371 stress 0.1094097 
    ## Run 372 stress 0.09381855 
    ## Run 373 stress 0.1026216 
    ## Run 374 stress 0.09030399 
    ## Run 375 stress 0.09465905 
    ## Run 376 stress 0.09590922 
    ## Run 377 stress 0.0933729 
    ## Run 378 stress 0.09030406 
    ## Run 379 stress 0.08773478 
    ## Run 380 stress 0.09337289 
    ## Run 381 stress 0.08973875 
    ## Run 382 stress 0.08773481 
    ## Run 383 stress 0.08773482 
    ## Run 384 stress 0.09268347 
    ## Run 385 stress 0.3372334 
    ## Run 386 stress 0.09268322 
    ## Run 387 stress 0.09308948 
    ## Run 388 stress 0.09286084 
    ## Run 389 stress 0.09465907 
    ## Run 390 stress 0.09408001 
    ## Run 391 stress 0.08773485 
    ## Run 392 stress 0.09286092 
    ## Run 393 stress 0.1026218 
    ## Run 394 stress 0.09445475 
    ## Run 395 stress 0.08773472 
    ## Run 396 stress 0.08773468 
    ## Run 397 stress 0.09030407 
    ## Run 398 stress 0.09535503 
    ## Run 399 stress 0.09407983 
    ## Run 400 stress 0.09447308 
    ## Run 401 stress 0.08503518 
    ## Run 402 stress 0.08440251 
    ## ... Procrustes: rmse 5.940697e-06  max resid 1.278219e-05 
    ## ... Similar to previous best
    ## Run 403 stress 0.09465905 
    ## Run 404 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001179397  max resid 0.0002154053 
    ## ... Similar to previous best
    ## Run 405 stress 0.08973868 
    ## Run 406 stress 0.09030399 
    ## Run 407 stress 0.09268318 
    ## Run 408 stress 0.09465906 
    ## Run 409 stress 0.09159093 
    ## Run 410 stress 0.09268341 
    ## Run 411 stress 0.09408 
    ## Run 412 stress 0.09760823 
    ## Run 413 stress 0.09380578 
    ## Run 414 stress 0.0937428 
    ## Run 415 stress 0.09030404 
    ## Run 416 stress 0.09308962 
    ## Run 417 stress 0.1084715 
    ## Run 418 stress 0.09286093 
    ## Run 419 stress 0.09416132 
    ## Run 420 stress 0.08503456 
    ## Run 421 stress 0.09407999 
    ## Run 422 stress 0.08973866 
    ## Run 423 stress 0.09407995 
    ## Run 424 stress 0.09407106 
    ## Run 425 stress 0.09337252 
    ## Run 426 stress 0.08440257 
    ## ... Procrustes: rmse 0.0001230949  max resid 0.0002562109 
    ## ... Similar to previous best
    ## Run 427 stress 0.09969977 
    ## Run 428 stress 0.09286089 
    ## Run 429 stress 0.09445487 
    ## Run 430 stress 0.08773473 
    ## Run 431 stress 0.09407992 
    ## Run 432 stress 0.0850349 
    ## Run 433 stress 0.08440266 
    ## ... Procrustes: rmse 0.0002141818  max resid 0.0004393754 
    ## ... Similar to previous best
    ## Run 434 stress 0.09374218 
    ## Run 435 stress 0.09465905 
    ## Run 436 stress 0.08440261 
    ## ... Procrustes: rmse 0.0001668119  max resid 0.0003059999 
    ## ... Similar to previous best
    ## Run 437 stress 0.08773469 
    ## Run 438 stress 0.09337275 
    ## Run 439 stress 0.09445464 
    ## Run 440 stress 0.09030409 
    ## Run 441 stress 0.09403423 
    ## Run 442 stress 0.08773478 
    ## Run 443 stress 0.09407969 
    ## Run 444 stress 0.09337232 
    ## Run 445 stress 0.08503463 
    ## Run 446 stress 0.0926835 
    ## Run 447 stress 0.08973885 
    ## Run 448 stress 0.09465918 
    ## Run 449 stress 0.09168947 
    ## Run 450 stress 0.0930896 
    ## Run 451 stress 0.09669846 
    ## Run 452 stress 0.08503483 
    ## Run 453 stress 0.09337271 
    ## Run 454 stress 0.09969982 
    ## Run 455 stress 0.09447303 
    ## Run 456 stress 0.09408004 
    ## Run 457 stress 0.09268358 
    ## Run 458 stress 0.09969982 
    ## Run 459 stress 0.08503498 
    ## Run 460 stress 0.08773479 
    ## Run 461 stress 0.09159089 
    ## Run 462 stress 0.09308958 
    ## Run 463 stress 0.09168932 
    ## Run 464 stress 0.08973863 
    ## Run 465 stress 0.09159085 
    ## Run 466 stress 0.09030393 
    ## Run 467 stress 0.08503585 
    ## Run 468 stress 0.09380587 
    ## Run 469 stress 0.09030395 
    ## Run 470 stress 0.09407989 
    ## Run 471 stress 0.09416139 
    ## Run 472 stress 0.0926842 
    ## Run 473 stress 0.1038082 
    ## Run 474 stress 0.08973867 
    ## Run 475 stress 0.09159095 
    ## Run 476 stress 0.0940797 
    ## Run 477 stress 0.09159094 
    ## Run 478 stress 0.08773483 
    ## Run 479 stress 0.0926838 
    ## Run 480 stress 0.09535454 
    ## Run 481 stress 0.09447296 
    ## Run 482 stress 0.08773485 
    ## Run 483 stress 0.09030396 
    ## Run 484 stress 0.08973864 
    ## Run 485 stress 0.09760843 
    ## Run 486 stress 0.08773466 
    ## Run 487 stress 0.09374267 
    ## Run 488 stress 0.09308969 
    ## Run 489 stress 0.0850351 
    ## Run 490 stress 0.09465905 
    ## Run 491 stress 0.08503481 
    ## Run 492 stress 0.09969994 
    ## Run 493 stress 0.09407986 
    ## Run 494 stress 0.08440265 
    ## ... Procrustes: rmse 0.0002076376  max resid 0.0004052121 
    ## ... Similar to previous best
    ## Run 495 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001772595  max resid 0.0003884597 
    ## ... Similar to previous best
    ## Run 496 stress 0.09374282 
    ## Run 497 stress 0.09030412 
    ## Run 498 stress 0.09374258 
    ## Run 499 stress 0.09168933 
    ## Run 500 stress 0.09407976 
    ## *** Best solution repeated 8 times

``` r
# Ocean sites and mixed lakes
SD_beta_env_OM_NMDS <- metaMDS(SD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.06623458 
    ## Run 1 stress 0.06477885 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1518412  max resid 0.2644425 
    ## Run 2 stress 0.08226161 
    ## Run 3 stress 0.07166981 
    ## Run 4 stress 0.06477933 
    ## ... Procrustes: rmse 0.00235608  max resid 0.003964659 
    ## ... Similar to previous best
    ## Run 5 stress 0.06124494 
    ## ... New best solution
    ## ... Procrustes: rmse 0.09181535  max resid 0.2202598 
    ## Run 6 stress 0.06623462 
    ## Run 7 stress 0.07411944 
    ## Run 8 stress 0.06623462 
    ## Run 9 stress 0.07166975 
    ## Run 10 stress 0.06124509 
    ## ... Procrustes: rmse 0.0001347936  max resid 0.0003526973 
    ## ... Similar to previous best
    ## Run 11 stress 0.07166982 
    ## Run 12 stress 0.07411941 
    ## Run 13 stress 0.09469095 
    ## Run 14 stress 0.06623459 
    ## Run 15 stress 0.06477972 
    ## Run 16 stress 0.06124496 
    ## ... Procrustes: rmse 3.713297e-05  max resid 9.592539e-05 
    ## ... Similar to previous best
    ## Run 17 stress 0.08226158 
    ## Run 18 stress 0.06623461 
    ## Run 19 stress 0.07411939 
    ## Run 20 stress 0.06123361 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01500349  max resid 0.04174088 
    ## Run 21 stress 0.06124509 
    ## ... Procrustes: rmse 0.01516359  max resid 0.04198652 
    ## Run 22 stress 0.07166985 
    ## Run 23 stress 0.06623459 
    ## Run 24 stress 0.06623458 
    ## Run 25 stress 0.0822616 
    ## Run 26 stress 0.0612451 
    ## ... Procrustes: rmse 0.01516128  max resid 0.04197897 
    ## Run 27 stress 0.3252911 
    ## Run 28 stress 0.06124488 
    ## ... Procrustes: rmse 0.01487857  max resid 0.04122718 
    ## Run 29 stress 0.06124503 
    ## ... Procrustes: rmse 0.01511067  max resid 0.0418445 
    ## Run 30 stress 0.06623461 
    ## Run 31 stress 0.07411946 
    ## Run 32 stress 0.07166964 
    ## Run 33 stress 0.06477938 
    ## Run 34 stress 0.3121118 
    ## Run 35 stress 0.07411941 
    ## Run 36 stress 0.06124503 
    ## ... Procrustes: rmse 0.01508886  max resid 0.04178489 
    ## Run 37 stress 0.07411945 
    ## Run 38 stress 0.06123362 
    ## ... Procrustes: rmse 6.425744e-05  max resid 0.0001032107 
    ## ... Similar to previous best
    ## Run 39 stress 0.06623457 
    ## Run 40 stress 0.06124495 
    ## ... Procrustes: rmse 0.01452994  max resid 0.04028831 
    ## Run 41 stress 0.06124501 
    ## ... Procrustes: rmse 0.01508544  max resid 0.0417782 
    ## Run 42 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 1.163049e-05  max resid 2.9588e-05 
    ## ... Similar to previous best
    ## Run 43 stress 0.0647789 
    ## Run 44 stress 0.2381016 
    ## Run 45 stress 0.07411938 
    ## Run 46 stress 0.3131458 
    ## Run 47 stress 0.07411962 
    ## Run 48 stress 0.06477928 
    ## Run 49 stress 0.07411939 
    ## Run 50 stress 0.06124503 
    ## ... Procrustes: rmse 0.0151038  max resid 0.04182486 
    ## Run 51 stress 0.06623457 
    ## Run 52 stress 0.08226168 
    ## Run 53 stress 0.08226169 
    ## Run 54 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 6.953719e-06  max resid 8.988161e-06 
    ## ... Similar to previous best
    ## Run 55 stress 0.06124507 
    ## ... Procrustes: rmse 0.01439864  max resid 0.0399345 
    ## Run 56 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 2.727405e-06  max resid 6.059337e-06 
    ## ... Similar to previous best
    ## Run 57 stress 0.0741194 
    ## Run 58 stress 0.07411945 
    ## Run 59 stress 0.0662346 
    ## Run 60 stress 0.06124505 
    ## ... Procrustes: rmse 0.01441904  max resid 0.03998881 
    ## Run 61 stress 0.0741194 
    ## Run 62 stress 0.07411943 
    ## Run 63 stress 0.2345379 
    ## Run 64 stress 0.07166968 
    ## Run 65 stress 0.06124506 
    ## ... Procrustes: rmse 0.01513878  max resid 0.04191922 
    ## Run 66 stress 0.06124497 
    ## ... Procrustes: rmse 0.01504293  max resid 0.0416628 
    ## Run 67 stress 0.07411945 
    ## Run 68 stress 0.07411944 
    ## Run 69 stress 0.0662346 
    ## Run 70 stress 0.07411939 
    ## Run 71 stress 0.07411953 
    ## Run 72 stress 0.06623462 
    ## Run 73 stress 0.0662346 
    ## Run 74 stress 0.06477868 
    ## Run 75 stress 0.0716697 
    ## Run 76 stress 0.07411956 
    ## Run 77 stress 0.3271653 
    ## Run 78 stress 0.0741195 
    ## Run 79 stress 0.06477835 
    ## Run 80 stress 0.06477852 
    ## Run 81 stress 0.06124522 
    ## ... Procrustes: rmse 0.01525575  max resid 0.04223204 
    ## Run 82 stress 0.06124512 
    ## ... Procrustes: rmse 0.01518828  max resid 0.04205142 
    ## Run 83 stress 0.06477904 
    ## Run 84 stress 0.06124518 
    ## ... Procrustes: rmse 0.01523391  max resid 0.04217386 
    ## Run 85 stress 0.0662346 
    ## Run 86 stress 0.3119309 
    ## Run 87 stress 0.06123361 
    ## ... Procrustes: rmse 4.786476e-05  max resid 5.708233e-05 
    ## ... Similar to previous best
    ## Run 88 stress 0.06477881 
    ## Run 89 stress 0.07166983 
    ## Run 90 stress 0.06623458 
    ## Run 91 stress 0.2963972 
    ## Run 92 stress 0.07411941 
    ## Run 93 stress 0.06477943 
    ## Run 94 stress 0.07166978 
    ## Run 95 stress 0.06623457 
    ## Run 96 stress 0.07411945 
    ## Run 97 stress 0.06124524 
    ## ... Procrustes: rmse 0.01526847  max resid 0.04226681 
    ## Run 98 stress 0.08226166 
    ## Run 99 stress 0.07166983 
    ## Run 100 stress 0.0612449 
    ## ... Procrustes: rmse 0.01495145  max resid 0.04141757 
    ## Run 101 stress 0.0716697 
    ## Run 102 stress 0.06623465 
    ## Run 103 stress 0.06623458 
    ## Run 104 stress 0.3258233 
    ## Run 105 stress 0.07411939 
    ## Run 106 stress 0.06124499 
    ## ... Procrustes: rmse 0.01505412  max resid 0.04169063 
    ## Run 107 stress 0.07411948 
    ## Run 108 stress 0.061245 
    ## ... Procrustes: rmse 0.01507063  max resid 0.04173543 
    ## Run 109 stress 0.07166966 
    ## Run 110 stress 0.06623462 
    ## Run 111 stress 0.08226156 
    ## Run 112 stress 0.07411948 
    ## Run 113 stress 0.06623457 
    ## Run 114 stress 0.06623463 
    ## Run 115 stress 0.06124531 
    ## ... Procrustes: rmse 0.01528932  max resid 0.04231897 
    ## Run 116 stress 0.06623458 
    ## Run 117 stress 0.08226158 
    ## Run 118 stress 0.07166967 
    ## Run 119 stress 0.06477846 
    ## Run 120 stress 0.06477952 
    ## Run 121 stress 0.07166965 
    ## Run 122 stress 0.2381016 
    ## Run 123 stress 0.07166986 
    ## Run 124 stress 0.06123364 
    ## ... Procrustes: rmse 5.513539e-05  max resid 7.634229e-05 
    ## ... Similar to previous best
    ## Run 125 stress 0.07166987 
    ## Run 126 stress 0.06123361 
    ## ... Procrustes: rmse 3.573085e-05  max resid 5.302793e-05 
    ## ... Similar to previous best
    ## Run 127 stress 0.3121119 
    ## Run 128 stress 0.08226159 
    ## Run 129 stress 0.3407173 
    ## Run 130 stress 0.08226156 
    ## Run 131 stress 0.06124498 
    ## ... Procrustes: rmse 0.01505091  max resid 0.0416816 
    ## Run 132 stress 0.06623462 
    ## Run 133 stress 0.08226159 
    ## Run 134 stress 0.06477847 
    ## Run 135 stress 0.06623457 
    ## Run 136 stress 0.06124537 
    ## ... Procrustes: rmse 0.0141933  max resid 0.03938121 
    ## Run 137 stress 0.0741195 
    ## Run 138 stress 0.06477845 
    ## Run 139 stress 0.06124507 
    ## ... Procrustes: rmse 0.01514816  max resid 0.0419442 
    ## Run 140 stress 0.06623464 
    ## Run 141 stress 0.06623457 
    ## Run 142 stress 0.2567807 
    ## Run 143 stress 0.07166972 
    ## Run 144 stress 0.06477918 
    ## Run 145 stress 0.06124486 
    ## ... Procrustes: rmse 0.01482504  max resid 0.04107886 
    ## Run 146 stress 0.06124492 
    ## ... Procrustes: rmse 0.01457192  max resid 0.04039997 
    ## Run 147 stress 0.06477861 
    ## Run 148 stress 0.07411939 
    ## Run 149 stress 0.0612452 
    ## ... Procrustes: rmse 0.01524749  max resid 0.04220957 
    ## Run 150 stress 0.07411944 
    ## Run 151 stress 0.07411956 
    ## Run 152 stress 0.07166971 
    ## Run 153 stress 0.0612336 
    ## ... Procrustes: rmse 8.946049e-06  max resid 1.468301e-05 
    ## ... Similar to previous best
    ## Run 154 stress 0.06124496 
    ## ... Procrustes: rmse 0.01503042  max resid 0.04162977 
    ## Run 155 stress 0.07411946 
    ## Run 156 stress 0.06123361 
    ## ... Procrustes: rmse 1.587556e-05  max resid 2.362113e-05 
    ## ... Similar to previous best
    ## Run 157 stress 0.07166969 
    ## Run 158 stress 0.06623457 
    ## Run 159 stress 0.0741194 
    ## Run 160 stress 0.06623463 
    ## Run 161 stress 0.07166964 
    ## Run 162 stress 0.06623463 
    ## Run 163 stress 0.07411942 
    ## Run 164 stress 0.06124511 
    ## ... Procrustes: rmse 0.01518302  max resid 0.04203725 
    ## Run 165 stress 0.06124509 
    ## ... Procrustes: rmse 0.01510039  max resid 0.04181276 
    ## Run 166 stress 0.3142495 
    ## Run 167 stress 0.08226159 
    ## Run 168 stress 0.06477827 
    ## Run 169 stress 0.06623457 
    ## Run 170 stress 0.06124512 
    ## ... Procrustes: rmse 0.01518912  max resid 0.04205295 
    ## Run 171 stress 0.06124507 
    ## ... Procrustes: rmse 0.01512132  max resid 0.04187413 
    ## Run 172 stress 0.07411939 
    ## Run 173 stress 0.06623462 
    ## Run 174 stress 0.09469061 
    ## Run 175 stress 0.3109017 
    ## Run 176 stress 0.06124517 
    ## ... Procrustes: rmse 0.01520809  max resid 0.04210617 
    ## Run 177 stress 0.07411943 
    ## Run 178 stress 0.06623462 
    ## Run 179 stress 0.0741194 
    ## Run 180 stress 0.07411939 
    ## Run 181 stress 0.06623462 
    ## Run 182 stress 0.06623461 
    ## Run 183 stress 0.07411941 
    ## Run 184 stress 0.06623457 
    ## Run 185 stress 0.2592302 
    ## Run 186 stress 0.07411943 
    ## Run 187 stress 0.06124503 
    ## ... Procrustes: rmse 0.01508088  max resid 0.04176176 
    ## Run 188 stress 0.07411947 
    ## Run 189 stress 0.3131327 
    ## Run 190 stress 0.0741194 
    ## Run 191 stress 0.07411942 
    ## Run 192 stress 0.06477843 
    ## Run 193 stress 0.06477988 
    ## Run 194 stress 0.06123361 
    ## ... Procrustes: rmse 3.114222e-05  max resid 4.003063e-05 
    ## ... Similar to previous best
    ## Run 195 stress 0.06477915 
    ## Run 196 stress 0.07166979 
    ## Run 197 stress 0.07411952 
    ## Run 198 stress 0.0716697 
    ## Run 199 stress 0.0647796 
    ## Run 200 stress 0.06124514 
    ## ... Procrustes: rmse 0.01518987  max resid 0.04205699 
    ## Run 201 stress 0.06477863 
    ## Run 202 stress 0.07166969 
    ## Run 203 stress 0.07411951 
    ## Run 204 stress 0.08226167 
    ## Run 205 stress 0.08226162 
    ## Run 206 stress 0.06124495 
    ## ... Procrustes: rmse 0.01502469  max resid 0.04161382 
    ## Run 207 stress 0.07411947 
    ## Run 208 stress 0.06623459 
    ## Run 209 stress 0.3407339 
    ## Run 210 stress 0.07411949 
    ## Run 211 stress 0.06124494 
    ## ... Procrustes: rmse 0.01498622  max resid 0.04151181 
    ## Run 212 stress 0.07411939 
    ## Run 213 stress 0.06477856 
    ## Run 214 stress 0.07411946 
    ## Run 215 stress 0.2345379 
    ## Run 216 stress 0.06623462 
    ## Run 217 stress 0.07411941 
    ## Run 218 stress 0.06124499 
    ## ... Procrustes: rmse 0.01506841  max resid 0.04173151 
    ## Run 219 stress 0.06123361 
    ## ... Procrustes: rmse 2.754135e-05  max resid 6.408148e-05 
    ## ... Similar to previous best
    ## Run 220 stress 0.2324716 
    ## Run 221 stress 0.07411953 
    ## Run 222 stress 0.06477879 
    ## Run 223 stress 0.08226162 
    ## Run 224 stress 0.07411941 
    ## Run 225 stress 0.07411938 
    ## Run 226 stress 0.06477866 
    ## Run 227 stress 0.06477889 
    ## Run 228 stress 0.06623459 
    ## Run 229 stress 0.06623458 
    ## Run 230 stress 0.08226172 
    ## Run 231 stress 0.06623461 
    ## Run 232 stress 0.06477852 
    ## Run 233 stress 0.08226169 
    ## Run 234 stress 0.06623457 
    ## Run 235 stress 0.06623457 
    ## Run 236 stress 0.07411938 
    ## Run 237 stress 0.06623458 
    ## Run 238 stress 0.06124505 
    ## ... Procrustes: rmse 0.01512543  max resid 0.04188257 
    ## Run 239 stress 0.0612336 
    ## ... Procrustes: rmse 7.233169e-06  max resid 1.224294e-05 
    ## ... Similar to previous best
    ## Run 240 stress 0.07411942 
    ## Run 241 stress 0.06477929 
    ## Run 242 stress 0.07166986 
    ## Run 243 stress 0.06477837 
    ## Run 244 stress 0.06477843 
    ## Run 245 stress 0.07411941 
    ## Run 246 stress 0.061245 
    ## ... Procrustes: rmse 0.01508647  max resid 0.04177889 
    ## Run 247 stress 0.07411941 
    ## Run 248 stress 0.0662346 
    ## Run 249 stress 0.06124498 
    ## ... Procrustes: rmse 0.01502971  max resid 0.04162864 
    ## Run 250 stress 0.06124503 
    ## ... Procrustes: rmse 0.01443655  max resid 0.04003583 
    ## Run 251 stress 0.0716697 
    ## Run 252 stress 0.06477858 
    ## Run 253 stress 0.0647791 
    ## Run 254 stress 0.06477835 
    ## Run 255 stress 0.06477957 
    ## Run 256 stress 0.06623459 
    ## Run 257 stress 0.06477884 
    ## Run 258 stress 0.06623458 
    ## Run 259 stress 0.06623464 
    ## Run 260 stress 0.06623459 
    ## Run 261 stress 0.06477844 
    ## Run 262 stress 0.07411942 
    ## Run 263 stress 0.06623468 
    ## Run 264 stress 0.06623458 
    ## Run 265 stress 0.07166965 
    ## Run 266 stress 0.0716697 
    ## Run 267 stress 0.07166964 
    ## Run 268 stress 0.06623467 
    ## Run 269 stress 0.07411942 
    ## Run 270 stress 0.06623457 
    ## Run 271 stress 0.07411941 
    ## Run 272 stress 0.06623457 
    ## Run 273 stress 0.07411943 
    ## Run 274 stress 0.3101175 
    ## Run 275 stress 0.06477837 
    ## Run 276 stress 0.06623458 
    ## Run 277 stress 0.08226156 
    ## Run 278 stress 0.07166969 
    ## Run 279 stress 0.0662346 
    ## Run 280 stress 0.06477914 
    ## Run 281 stress 0.07411943 
    ## Run 282 stress 0.06477855 
    ## Run 283 stress 0.07411938 
    ## Run 284 stress 0.07166968 
    ## Run 285 stress 0.06124498 
    ## ... Procrustes: rmse 0.01504834  max resid 0.04167557 
    ## Run 286 stress 0.06124492 
    ## ... Procrustes: rmse 0.01498523  max resid 0.04150799 
    ## Run 287 stress 0.07411946 
    ## Run 288 stress 0.2548998 
    ## Run 289 stress 0.3284433 
    ## Run 290 stress 0.0647792 
    ## Run 291 stress 0.07411939 
    ## Run 292 stress 0.06124509 
    ## ... Procrustes: rmse 0.01515205  max resid 0.04195268 
    ## Run 293 stress 0.07411947 
    ## Run 294 stress 0.06623458 
    ## Run 295 stress 0.3023568 
    ## Run 296 stress 0.06477854 
    ## Run 297 stress 0.07411941 
    ## Run 298 stress 0.0662346 
    ## Run 299 stress 0.06124521 
    ## ... Procrustes: rmse 0.01524869  max resid 0.0422137 
    ## Run 300 stress 0.07411938 
    ## Run 301 stress 0.06123362 
    ## ... Procrustes: rmse 5.988089e-05  max resid 7.237055e-05 
    ## ... Similar to previous best
    ## Run 302 stress 0.07411944 
    ## Run 303 stress 0.06623461 
    ## Run 304 stress 0.06123362 
    ## ... Procrustes: rmse 5.302403e-05  max resid 0.0001224822 
    ## ... Similar to previous best
    ## Run 305 stress 0.0612452 
    ## ... Procrustes: rmse 0.0152281  max resid 0.04215571 
    ## Run 306 stress 0.06623457 
    ## Run 307 stress 0.06124509 
    ## ... Procrustes: rmse 0.01516513  max resid 0.04198893 
    ## Run 308 stress 0.07411947 
    ## Run 309 stress 0.06477838 
    ## Run 310 stress 0.06623458 
    ## Run 311 stress 0.07166967 
    ## Run 312 stress 0.07166967 
    ## Run 313 stress 0.06124507 
    ## ... Procrustes: rmse 0.015141  max resid 0.04192371 
    ## Run 314 stress 0.07166967 
    ## Run 315 stress 0.06623458 
    ## Run 316 stress 0.06477828 
    ## Run 317 stress 0.06124523 
    ## ... Procrustes: rmse 0.01427887  max resid 0.03961192 
    ## Run 318 stress 0.0612336 
    ## ... Procrustes: rmse 4.824582e-06  max resid 8.748789e-06 
    ## ... Similar to previous best
    ## Run 319 stress 0.0647798 
    ## Run 320 stress 0.06623464 
    ## Run 321 stress 0.0716697 
    ## Run 322 stress 0.06477956 
    ## Run 323 stress 0.07411945 
    ## Run 324 stress 0.08226169 
    ## Run 325 stress 0.07411939 
    ## Run 326 stress 0.06124508 
    ## ... Procrustes: rmse 0.01512499  max resid 0.04187964 
    ## Run 327 stress 0.08226157 
    ## Run 328 stress 0.06124524 
    ## ... Procrustes: rmse 0.01526187  max resid 0.04224961 
    ## Run 329 stress 0.061245 
    ## ... Procrustes: rmse 0.01508312  max resid 0.04176944 
    ## Run 330 stress 0.06623458 
    ## Run 331 stress 0.07166977 
    ## Run 332 stress 0.07411951 
    ## Run 333 stress 0.2324715 
    ## Run 334 stress 0.06477935 
    ## Run 335 stress 0.06623457 
    ## Run 336 stress 0.0612336 
    ## ... Procrustes: rmse 1.892017e-05  max resid 2.757041e-05 
    ## ... Similar to previous best
    ## Run 337 stress 0.0741194 
    ## Run 338 stress 0.06623458 
    ## Run 339 stress 0.251859 
    ## Run 340 stress 0.06623459 
    ## Run 341 stress 0.06477856 
    ## Run 342 stress 0.06477889 
    ## Run 343 stress 0.0647786 
    ## Run 344 stress 0.07411938 
    ## Run 345 stress 0.07166981 
    ## Run 346 stress 0.07166971 
    ## Run 347 stress 0.06123361 
    ## ... Procrustes: rmse 2.16619e-05  max resid 4.681658e-05 
    ## ... Similar to previous best
    ## Run 348 stress 0.2362366 
    ## Run 349 stress 0.06124494 
    ## ... Procrustes: rmse 0.01454269  max resid 0.04032165 
    ## Run 350 stress 0.06623457 
    ## Run 351 stress 0.07411944 
    ## Run 352 stress 0.06623464 
    ## Run 353 stress 0.08226161 
    ## Run 354 stress 0.08226162 
    ## Run 355 stress 0.06477918 
    ## Run 356 stress 0.06477862 
    ## Run 357 stress 0.06124499 
    ## ... Procrustes: rmse 0.01503526  max resid 0.0416395 
    ## Run 358 stress 0.06477844 
    ## Run 359 stress 0.3113442 
    ## Run 360 stress 0.0647783 
    ## Run 361 stress 0.06477986 
    ## Run 362 stress 0.07411943 
    ## Run 363 stress 0.06124504 
    ## ... Procrustes: rmse 0.01511325  max resid 0.04184926 
    ## Run 364 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 2.621854e-06  max resid 4.702836e-06 
    ## ... Similar to previous best
    ## Run 365 stress 0.07166973 
    ## Run 366 stress 0.06123361 
    ## ... Procrustes: rmse 3.692494e-05  max resid 8.464233e-05 
    ## ... Similar to previous best
    ## Run 367 stress 0.06623459 
    ## Run 368 stress 0.06623457 
    ## Run 369 stress 0.07166978 
    ## Run 370 stress 0.07166967 
    ## Run 371 stress 0.2381016 
    ## Run 372 stress 0.07411943 
    ## Run 373 stress 0.07166967 
    ## Run 374 stress 0.06623457 
    ## Run 375 stress 0.07411949 
    ## Run 376 stress 0.07166974 
    ## Run 377 stress 0.3270102 
    ## Run 378 stress 0.2590308 
    ## Run 379 stress 0.3157933 
    ## Run 380 stress 0.0741194 
    ## Run 381 stress 0.06623462 
    ## Run 382 stress 0.07411945 
    ## Run 383 stress 0.07166974 
    ## Run 384 stress 0.09623017 
    ## Run 385 stress 0.06623458 
    ## Run 386 stress 0.0662346 
    ## Run 387 stress 0.07411942 
    ## Run 388 stress 0.06124519 
    ## ... Procrustes: rmse 0.01523193  max resid 0.04216944 
    ## Run 389 stress 0.06623457 
    ## Run 390 stress 0.08226171 
    ## Run 391 stress 0.06123363 
    ## ... Procrustes: rmse 7.777756e-05  max resid 0.0001188124 
    ## ... Similar to previous best
    ## Run 392 stress 0.08226173 
    ## Run 393 stress 0.0741194 
    ## Run 394 stress 0.06124503 
    ## ... Procrustes: rmse 0.01445684  max resid 0.04008937 
    ## Run 395 stress 0.0612336 
    ## ... Procrustes: rmse 3.917122e-06  max resid 8.165702e-06 
    ## ... Similar to previous best
    ## Run 396 stress 0.06477839 
    ## Run 397 stress 0.07411938 
    ## Run 398 stress 0.06124488 
    ## ... Procrustes: rmse 0.01490085  max resid 0.04128226 
    ## Run 399 stress 0.07166966 
    ## Run 400 stress 0.06623457 
    ## Run 401 stress 0.06623457 
    ## Run 402 stress 0.06477909 
    ## Run 403 stress 0.07166975 
    ## Run 404 stress 0.08226158 
    ## Run 405 stress 0.06623465 
    ## Run 406 stress 0.0647797 
    ## Run 407 stress 0.08226166 
    ## Run 408 stress 0.07411939 
    ## Run 409 stress 0.0741194 
    ## Run 410 stress 0.06124489 
    ## ... Procrustes: rmse 0.01489948  max resid 0.04128224 
    ## Run 411 stress 0.06477862 
    ## Run 412 stress 0.0662346 
    ## Run 413 stress 0.2585395 
    ## Run 414 stress 0.0612336 
    ## ... Procrustes: rmse 5.638681e-06  max resid 7.261275e-06 
    ## ... Similar to previous best
    ## Run 415 stress 0.3022789 
    ## Run 416 stress 0.07411938 
    ## Run 417 stress 0.06623457 
    ## Run 418 stress 0.08226159 
    ## Run 419 stress 0.07411939 
    ## Run 420 stress 0.06623459 
    ## Run 421 stress 0.07166965 
    ## Run 422 stress 0.06477935 
    ## Run 423 stress 0.06623458 
    ## Run 424 stress 0.0612336 
    ## ... Procrustes: rmse 1.332151e-05  max resid 1.671631e-05 
    ## ... Similar to previous best
    ## Run 425 stress 0.06623462 
    ## Run 426 stress 0.0716698 
    ## Run 427 stress 0.07166983 
    ## Run 428 stress 0.06623461 
    ## Run 429 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 1.829878e-06  max resid 4.517162e-06 
    ## ... Similar to previous best
    ## Run 430 stress 0.06623468 
    ## Run 431 stress 0.07166972 
    ## Run 432 stress 0.09093542 
    ## Run 433 stress 0.07411941 
    ## Run 434 stress 0.07411941 
    ## Run 435 stress 0.0612336 
    ## ... Procrustes: rmse 6.09775e-06  max resid 1.154597e-05 
    ## ... Similar to previous best
    ## Run 436 stress 0.3266727 
    ## Run 437 stress 0.07411941 
    ## Run 438 stress 0.2362366 
    ## Run 439 stress 0.3270108 
    ## Run 440 stress 0.07166964 
    ## Run 441 stress 0.06623459 
    ## Run 442 stress 0.07166978 
    ## Run 443 stress 0.0822616 
    ## Run 444 stress 0.0822616 
    ## Run 445 stress 0.0716697 
    ## Run 446 stress 0.07411957 
    ## Run 447 stress 0.0822617 
    ## Run 448 stress 0.06124497 
    ## ... Procrustes: rmse 0.01502864  max resid 0.04162247 
    ## Run 449 stress 0.06623462 
    ## Run 450 stress 0.07411947 
    ## Run 451 stress 0.07411939 
    ## Run 452 stress 0.07411941 
    ## Run 453 stress 0.07411951 
    ## Run 454 stress 0.06124495 
    ## ... Procrustes: rmse 0.01500701  max resid 0.04156486 
    ## Run 455 stress 0.06623458 
    ## Run 456 stress 0.06623458 
    ## Run 457 stress 0.06477933 
    ## Run 458 stress 0.0662346 
    ## Run 459 stress 0.08226167 
    ## Run 460 stress 0.06124489 
    ## ... Procrustes: rmse 0.01491437  max resid 0.04131901 
    ## Run 461 stress 0.07411944 
    ## Run 462 stress 0.07166992 
    ## Run 463 stress 0.08226169 
    ## Run 464 stress 0.06124516 
    ## ... Procrustes: rmse 0.01512638  max resid 0.04188321 
    ## Run 465 stress 0.07166973 
    ## Run 466 stress 0.09300798 
    ## Run 467 stress 0.07411939 
    ## Run 468 stress 0.07166969 
    ## Run 469 stress 0.06623466 
    ## Run 470 stress 0.2381016 
    ## Run 471 stress 0.06623463 
    ## Run 472 stress 0.06477843 
    ## Run 473 stress 0.06477892 
    ## Run 474 stress 0.06477912 
    ## Run 475 stress 0.08226166 
    ## Run 476 stress 0.06623457 
    ## Run 477 stress 0.06123363 
    ## ... Procrustes: rmse 7.21643e-05  max resid 8.792706e-05 
    ## ... Similar to previous best
    ## Run 478 stress 0.07411939 
    ## Run 479 stress 0.08226171 
    ## Run 480 stress 0.07411938 
    ## Run 481 stress 0.08226157 
    ## Run 482 stress 0.07166968 
    ## Run 483 stress 0.06477915 
    ## Run 484 stress 0.06623464 
    ## Run 485 stress 0.08226156 
    ## Run 486 stress 0.3157936 
    ## Run 487 stress 0.07411946 
    ## Run 488 stress 0.06623458 
    ## Run 489 stress 0.06623459 
    ## Run 490 stress 0.3022791 
    ## Run 491 stress 0.06124506 
    ## ... Procrustes: rmse 0.01440997  max resid 0.03996453 
    ## Run 492 stress 0.06477897 
    ## Run 493 stress 0.07411959 
    ## Run 494 stress 0.07411945 
    ## Run 495 stress 0.06123361 
    ## ... Procrustes: rmse 1.477055e-05  max resid 2.635249e-05 
    ## ... Similar to previous best
    ## Run 496 stress 0.0647783 
    ## Run 497 stress 0.07411942 
    ## Run 498 stress 0.06124516 
    ## ... Procrustes: rmse 0.01436467  max resid 0.03984033 
    ## Run 499 stress 0.06623458 
    ## Run 500 stress 0.06124522 
    ## ... Procrustes: rmse 0.01525212  max resid 0.04222293 
    ## *** Best solution repeated 4 times

``` r
# Stratified lakes and ocean sites
SD_beta_env_SO_NMDS <- metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.59483e-05 
    ## Run 1 stress 9.74824e-05 
    ## ... Procrustes: rmse 0.0001804319  max resid 0.0003934323 
    ## ... Similar to previous best
    ## Run 2 stress 9.439531e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002135304  max resid 0.0003816813 
    ## ... Similar to previous best
    ## Run 3 stress 9.663334e-05 
    ## ... Procrustes: rmse 0.0002332695  max resid 0.0003434411 
    ## ... Similar to previous best
    ## Run 4 stress 8.5837e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002129844  max resid 0.00030838 
    ## ... Similar to previous best
    ## Run 5 stress 8.655576e-05 
    ## ... Procrustes: rmse 0.0001284006  max resid 0.0002871944 
    ## ... Similar to previous best
    ## Run 6 stress 9.146762e-05 
    ## ... Procrustes: rmse 0.0001808134  max resid 0.0004163403 
    ## ... Similar to previous best
    ## Run 7 stress 9.316035e-05 
    ## ... Procrustes: rmse 0.0002048619  max resid 0.0003678355 
    ## ... Similar to previous best
    ## Run 8 stress 9.677026e-05 
    ## ... Procrustes: rmse 0.0001759588  max resid 0.000341823 
    ## ... Similar to previous best
    ## Run 9 stress 9.599833e-05 
    ## ... Procrustes: rmse 0.0001687351  max resid 0.0002718025 
    ## ... Similar to previous best
    ## Run 10 stress 9.838825e-05 
    ## ... Procrustes: rmse 0.0002075264  max resid 0.0003514261 
    ## ... Similar to previous best
    ## Run 11 stress 9.526343e-05 
    ## ... Procrustes: rmse 0.0001252659  max resid 0.0002626972 
    ## ... Similar to previous best
    ## Run 12 stress 8.49824e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001471383  max resid 0.0002441872 
    ## ... Similar to previous best
    ## Run 13 stress 9.418867e-05 
    ## ... Procrustes: rmse 0.0001563186  max resid 0.0002593675 
    ## ... Similar to previous best
    ## Run 14 stress 9.644937e-05 
    ## ... Procrustes: rmse 0.000100176  max resid 0.0001656249 
    ## ... Similar to previous best
    ## Run 15 stress 9.650065e-05 
    ## ... Procrustes: rmse 0.0001633813  max resid 0.0003633989 
    ## ... Similar to previous best
    ## Run 16 stress 9.552275e-05 
    ## ... Procrustes: rmse 0.0001551469  max resid 0.0002614988 
    ## ... Similar to previous best
    ## Run 17 stress 8.391503e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002189705  max resid 0.0003554991 
    ## ... Similar to previous best
    ## Run 18 stress 9.841953e-05 
    ## ... Procrustes: rmse 0.0002527611  max resid 0.0004105203 
    ## ... Similar to previous best
    ## Run 19 stress 9.730964e-05 
    ## ... Procrustes: rmse 0.0002278784  max resid 0.0003580917 
    ## ... Similar to previous best
    ## Run 20 stress 9.597069e-05 
    ## ... Procrustes: rmse 0.0002591814  max resid 0.0003587168 
    ## ... Similar to previous best
    ## Run 21 stress 7.911804e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001818425  max resid 0.0003037658 
    ## ... Similar to previous best
    ## Run 22 stress 9.800909e-05 
    ## ... Procrustes: rmse 0.0001991257  max resid 0.00033773 
    ## ... Similar to previous best
    ## Run 23 stress 9.780865e-05 
    ## ... Procrustes: rmse 0.0001403024  max resid 0.000317151 
    ## ... Similar to previous best
    ## Run 24 stress 9.903091e-05 
    ## ... Procrustes: rmse 0.0002148164  max resid 0.0004357045 
    ## ... Similar to previous best
    ## Run 25 stress 7.666178e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001488358  max resid 0.0003451044 
    ## ... Similar to previous best
    ## Run 26 stress 9.614364e-05 
    ## ... Procrustes: rmse 0.0001473238  max resid 0.0003230623 
    ## ... Similar to previous best
    ## Run 27 stress 7.875122e-05 
    ## ... Procrustes: rmse 0.0001556626  max resid 0.0003559873 
    ## ... Similar to previous best
    ## Run 28 stress 7.91981e-05 
    ## ... Procrustes: rmse 4.767212e-05  max resid 8.772725e-05 
    ## ... Similar to previous best
    ## Run 29 stress 9.822576e-05 
    ## ... Procrustes: rmse 0.0001862235  max resid 0.0004228676 
    ## ... Similar to previous best
    ## Run 30 stress 9.061029e-05 
    ## ... Procrustes: rmse 0.0001097708  max resid 0.0001827733 
    ## ... Similar to previous best
    ## Run 31 stress 9.780714e-05 
    ## ... Procrustes: rmse 0.0001495202  max resid 0.0003283357 
    ## ... Similar to previous best
    ## Run 32 stress 9.747778e-05 
    ## ... Procrustes: rmse 0.0001787309  max resid 0.0003341559 
    ## ... Similar to previous best
    ## Run 33 stress 9.831901e-05 
    ## ... Procrustes: rmse 0.0001528762  max resid 0.0003203556 
    ## ... Similar to previous best
    ## Run 34 stress 9.562682e-05 
    ## ... Procrustes: rmse 0.0001791315  max resid 0.000413786 
    ## ... Similar to previous best
    ## Run 35 stress 9.445337e-05 
    ## ... Procrustes: rmse 0.0001615007  max resid 0.0003366645 
    ## ... Similar to previous best
    ## Run 36 stress 9.776943e-05 
    ## ... Procrustes: rmse 0.0001901549  max resid 0.0004060719 
    ## ... Similar to previous best
    ## Run 37 stress 8.866084e-05 
    ## ... Procrustes: rmse 0.0001355834  max resid 0.0003066892 
    ## ... Similar to previous best
    ## Run 38 stress 9.726754e-05 
    ## ... Procrustes: rmse 0.0001436267  max resid 0.0003313016 
    ## ... Similar to previous best
    ## Run 39 stress 9.895828e-05 
    ## ... Procrustes: rmse 0.0001862899  max resid 0.0004220628 
    ## ... Similar to previous best
    ## Run 40 stress 9.368197e-05 
    ## ... Procrustes: rmse 9.09292e-05  max resid 0.0002137824 
    ## ... Similar to previous best
    ## Run 41 stress 8.749974e-05 
    ## ... Procrustes: rmse 0.0001603591  max resid 0.0003129428 
    ## ... Similar to previous best
    ## Run 42 stress 9.028617e-05 
    ## ... Procrustes: rmse 0.0001719872  max resid 0.0003996586 
    ## ... Similar to previous best
    ## Run 43 stress 9.659451e-05 
    ## ... Procrustes: rmse 0.0001887398  max resid 0.0004031865 
    ## ... Similar to previous best
    ## Run 44 stress 9.789792e-05 
    ## ... Procrustes: rmse 9.374144e-05  max resid 0.000207358 
    ## ... Similar to previous best
    ## Run 45 stress 9.898833e-05 
    ## ... Procrustes: rmse 0.0001385452  max resid 0.0002836084 
    ## ... Similar to previous best
    ## Run 46 stress 9.557549e-05 
    ## ... Procrustes: rmse 0.0001622288  max resid 0.0003857036 
    ## ... Similar to previous best
    ## Run 47 stress 9.079099e-05 
    ## ... Procrustes: rmse 7.572355e-05  max resid 0.0001339681 
    ## ... Similar to previous best
    ## Run 48 stress 9.758777e-05 
    ## ... Procrustes: rmse 6.980301e-05  max resid 0.0001093259 
    ## ... Similar to previous best
    ## Run 49 stress 9.522635e-05 
    ## ... Procrustes: rmse 0.0001184465  max resid 0.0001855559 
    ## ... Similar to previous best
    ## Run 50 stress 9.21367e-05 
    ## ... Procrustes: rmse 0.000174988  max resid 0.0004030768 
    ## ... Similar to previous best
    ## Run 51 stress 9.774964e-05 
    ## ... Procrustes: rmse 0.0001121692  max resid 0.0001746141 
    ## ... Similar to previous best
    ## Run 52 stress 9.905729e-05 
    ## ... Procrustes: rmse 0.0001526616  max resid 0.0003304671 
    ## ... Similar to previous best
    ## Run 53 stress 9.110898e-05 
    ## ... Procrustes: rmse 0.0001343513  max resid 0.000275665 
    ## ... Similar to previous best
    ## Run 54 stress 9.60828e-05 
    ## ... Procrustes: rmse 8.589539e-05  max resid 0.0001452189 
    ## ... Similar to previous best
    ## Run 55 stress 8.403428e-05 
    ## ... Procrustes: rmse 0.0001403364  max resid 0.0002944154 
    ## ... Similar to previous best
    ## Run 56 stress 9.909249e-05 
    ## ... Procrustes: rmse 6.820663e-05  max resid 0.000114155 
    ## ... Similar to previous best
    ## Run 57 stress 9.849074e-05 
    ## ... Procrustes: rmse 6.109922e-05  max resid 0.0001001582 
    ## ... Similar to previous best
    ## Run 58 stress 9.859927e-05 
    ## ... Procrustes: rmse 8.122323e-05  max resid 0.0001545949 
    ## ... Similar to previous best
    ## Run 59 stress 9.444596e-05 
    ## ... Procrustes: rmse 0.0001775234  max resid 0.0004314701 
    ## ... Similar to previous best
    ## Run 60 stress 9.951749e-05 
    ## ... Procrustes: rmse 0.0001460703  max resid 0.000319991 
    ## ... Similar to previous best
    ## Run 61 stress 9.980675e-05 
    ## ... Procrustes: rmse 9.383129e-05  max resid 0.0001552098 
    ## ... Similar to previous best
    ## Run 62 stress 8.694084e-05 
    ## ... Procrustes: rmse 0.0001520091  max resid 0.0002904619 
    ## ... Similar to previous best
    ## Run 63 stress 9.779628e-05 
    ## ... Procrustes: rmse 0.0001915896  max resid 0.0003290658 
    ## ... Similar to previous best
    ## Run 64 stress 0.2203301 
    ## Run 65 stress 9.981224e-05 
    ## ... Procrustes: rmse 0.0001901178  max resid 0.0003360721 
    ## ... Similar to previous best
    ## Run 66 stress 9.745003e-05 
    ## ... Procrustes: rmse 6.761923e-05  max resid 0.0001063801 
    ## ... Similar to previous best
    ## Run 67 stress 9.379326e-05 
    ## ... Procrustes: rmse 9.824591e-05  max resid 0.0002040353 
    ## ... Similar to previous best
    ## Run 68 stress 8.940828e-05 
    ## ... Procrustes: rmse 0.0001216096  max resid 0.0002808284 
    ## ... Similar to previous best
    ## Run 69 stress 9.870899e-05 
    ## ... Procrustes: rmse 7.122362e-05  max resid 0.0001150868 
    ## ... Similar to previous best
    ## Run 70 stress 8.402064e-05 
    ## ... Procrustes: rmse 0.0001007549  max resid 0.0001874745 
    ## ... Similar to previous best
    ## Run 71 stress 9.049432e-05 
    ## ... Procrustes: rmse 0.0001580667  max resid 0.0003044338 
    ## ... Similar to previous best
    ## Run 72 stress 9.577724e-05 
    ## ... Procrustes: rmse 0.0001639063  max resid 0.0003404474 
    ## ... Similar to previous best
    ## Run 73 stress 9.619225e-05 
    ## ... Procrustes: rmse 0.0001638188  max resid 0.0003400088 
    ## ... Similar to previous best
    ## Run 74 stress 9.77806e-05 
    ## ... Procrustes: rmse 0.0001053491  max resid 0.0001701017 
    ## ... Similar to previous best
    ## Run 75 stress 9.616509e-05 
    ## ... Procrustes: rmse 0.000136036  max resid 0.0003113397 
    ## ... Similar to previous best
    ## Run 76 stress 9.924482e-05 
    ## ... Procrustes: rmse 0.0001244325  max resid 0.0002100687 
    ## ... Similar to previous best
    ## Run 77 stress 9.483449e-05 
    ## ... Procrustes: rmse 0.0001477679  max resid 0.000289595 
    ## ... Similar to previous best
    ## Run 78 stress 9.676646e-05 
    ## ... Procrustes: rmse 0.0001964437  max resid 0.0003345983 
    ## ... Similar to previous best
    ## Run 79 stress 9.129909e-05 
    ## ... Procrustes: rmse 0.0001594662  max resid 0.0003012034 
    ## ... Similar to previous best
    ## Run 80 stress 9.799344e-05 
    ## ... Procrustes: rmse 0.0001545708  max resid 0.0003209194 
    ## ... Similar to previous best
    ## Run 81 stress 9.337034e-05 
    ## ... Procrustes: rmse 9.414793e-05  max resid 0.0001658394 
    ## ... Similar to previous best
    ## Run 82 stress 9.484583e-05 
    ## ... Procrustes: rmse 0.0001256827  max resid 0.0002906524 
    ## ... Similar to previous best
    ## Run 83 stress 8.581257e-05 
    ## ... Procrustes: rmse 0.000151012  max resid 0.0002855072 
    ## ... Similar to previous best
    ## Run 84 stress 9.653284e-05 
    ## ... Procrustes: rmse 0.0001289734  max resid 0.0003073693 
    ## ... Similar to previous best
    ## Run 85 stress 9.954107e-05 
    ## ... Procrustes: rmse 0.0001538143  max resid 0.0003321037 
    ## ... Similar to previous best
    ## Run 86 stress 9.579891e-05 
    ## ... Procrustes: rmse 7.960879e-05  max resid 0.0001164783 
    ## ... Similar to previous best
    ## Run 87 stress 9.857447e-05 
    ## ... Procrustes: rmse 0.0001858916  max resid 0.0003309081 
    ## ... Similar to previous best
    ## Run 88 stress 9.571785e-05 
    ## ... Procrustes: rmse 7.13179e-05  max resid 0.0001094026 
    ## ... Similar to previous best
    ## Run 89 stress 9.155332e-05 
    ## ... Procrustes: rmse 0.0001428538  max resid 0.000289512 
    ## ... Similar to previous best
    ## Run 90 stress 9.727974e-05 
    ## ... Procrustes: rmse 8.62625e-05  max resid 0.0001266545 
    ## ... Similar to previous best
    ## Run 91 stress 9.689318e-05 
    ## ... Procrustes: rmse 0.0001656403  max resid 0.0003112747 
    ## ... Similar to previous best
    ## Run 92 stress 9.606656e-05 
    ## ... Procrustes: rmse 0.0001747289  max resid 0.0004090962 
    ## ... Similar to previous best
    ## Run 93 stress 9.978677e-05 
    ## ... Procrustes: rmse 0.0001012026  max resid 0.0002310807 
    ## ... Similar to previous best
    ## Run 94 stress 8.937696e-05 
    ## ... Procrustes: rmse 0.0001563093  max resid 0.0002993488 
    ## ... Similar to previous best
    ## Run 95 stress 9.461178e-05 
    ## ... Procrustes: rmse 8.181823e-05  max resid 0.0001433617 
    ## ... Similar to previous best
    ## Run 96 stress 9.735016e-05 
    ## ... Procrustes: rmse 0.0001231032  max resid 0.0001904593 
    ## ... Similar to previous best
    ## Run 97 stress 9.711263e-05 
    ## ... Procrustes: rmse 0.0001565056  max resid 0.0003510014 
    ## ... Similar to previous best
    ## Run 98 stress 8.642888e-05 
    ## ... Procrustes: rmse 0.0001480141  max resid 0.0002842569 
    ## ... Similar to previous best
    ## Run 99 stress 9.402915e-05 
    ## ... Procrustes: rmse 6.622526e-05  max resid 0.000100271 
    ## ... Similar to previous best
    ## Run 100 stress 9.737584e-05 
    ## ... Procrustes: rmse 0.0001533581  max resid 0.0003193696 
    ## ... Similar to previous best
    ## Run 101 stress 9.514482e-05 
    ## ... Procrustes: rmse 0.000186574  max resid 0.0003990318 
    ## ... Similar to previous best
    ## Run 102 stress 9.924389e-05 
    ## ... Procrustes: rmse 7.399259e-05  max resid 0.000113891 
    ## ... Similar to previous best
    ## Run 103 stress 8.861528e-05 
    ## ... Procrustes: rmse 0.0001380001  max resid 0.0002629856 
    ## ... Similar to previous best
    ## Run 104 stress 9.804065e-05 
    ## ... Procrustes: rmse 0.000133608  max resid 0.0002432282 
    ## ... Similar to previous best
    ## Run 105 stress 7.942798e-05 
    ## ... Procrustes: rmse 0.0001710465  max resid 0.0002911982 
    ## ... Similar to previous best
    ## Run 106 stress 9.478884e-05 
    ## ... Procrustes: rmse 0.0001731156  max resid 0.0003278797 
    ## ... Similar to previous best
    ## Run 107 stress 9.937065e-05 
    ## ... Procrustes: rmse 0.0001700069  max resid 0.000347491 
    ## ... Similar to previous best
    ## Run 108 stress 9.675649e-05 
    ## ... Procrustes: rmse 0.0001757114  max resid 0.0003984968 
    ## ... Similar to previous best
    ## Run 109 stress 9.945871e-05 
    ## ... Procrustes: rmse 0.0001687268  max resid 0.0003204374 
    ## ... Similar to previous best
    ## Run 110 stress 9.497889e-05 
    ## ... Procrustes: rmse 0.0001442864  max resid 0.0003471551 
    ## ... Similar to previous best
    ## Run 111 stress 0.2386315 
    ## Run 112 stress 9.439408e-05 
    ## ... Procrustes: rmse 0.0001115289  max resid 0.0001811916 
    ## ... Similar to previous best
    ## Run 113 stress 9.643534e-05 
    ## ... Procrustes: rmse 0.0001816631  max resid 0.0003279817 
    ## ... Similar to previous best
    ## Run 114 stress 9.2753e-05 
    ## ... Procrustes: rmse 0.0001733448  max resid 0.0003161019 
    ## ... Similar to previous best
    ## Run 115 stress 9.579063e-05 
    ## ... Procrustes: rmse 6.637636e-05  max resid 0.0001048754 
    ## ... Similar to previous best
    ## Run 116 stress 9.304664e-05 
    ## ... Procrustes: rmse 0.0001763471  max resid 0.0004055909 
    ## ... Similar to previous best
    ## Run 117 stress 9.398837e-05 
    ## ... Procrustes: rmse 0.0001791795  max resid 0.0004114841 
    ## ... Similar to previous best
    ## Run 118 stress 9.556336e-05 
    ## ... Procrustes: rmse 7.074747e-05  max resid 0.0001154704 
    ## ... Similar to previous best
    ## Run 119 stress 9.864627e-05 
    ## ... Procrustes: rmse 0.0001436134  max resid 0.0003187851 
    ## ... Similar to previous best
    ## Run 120 stress 9.344186e-05 
    ## ... Procrustes: rmse 0.0001845151  max resid 0.0003208036 
    ## ... Similar to previous best
    ## Run 121 stress 9.806881e-05 
    ## ... Procrustes: rmse 0.0001558938  max resid 0.0003231656 
    ## ... Similar to previous best
    ## Run 122 stress 9.953079e-05 
    ## ... Procrustes: rmse 0.0001842173  max resid 0.0004150232 
    ## ... Similar to previous best
    ## Run 123 stress 8.895744e-05 
    ## ... Procrustes: rmse 0.0001549461  max resid 0.0002936541 
    ## ... Similar to previous best
    ## Run 124 stress 8.552857e-05 
    ## ... Procrustes: rmse 0.000162811  max resid 0.000384848 
    ## ... Similar to previous best
    ## Run 125 stress 9.711761e-05 
    ## ... Procrustes: rmse 0.0001666496  max resid 0.0003118023 
    ## ... Similar to previous best
    ## Run 126 stress 9.266785e-05 
    ## ... Procrustes: rmse 6.270836e-05  max resid 9.875364e-05 
    ## ... Similar to previous best
    ## Run 127 stress 9.447533e-05 
    ## ... Procrustes: rmse 6.560422e-05  max resid 0.0001033646 
    ## ... Similar to previous best
    ## Run 128 stress 9.645394e-05 
    ## ... Procrustes: rmse 0.0001472494  max resid 0.0003235928 
    ## ... Similar to previous best
    ## Run 129 stress 9.376401e-05 
    ## ... Procrustes: rmse 0.0001156631  max resid 0.0002838369 
    ## ... Similar to previous best
    ## Run 130 stress 9.916659e-05 
    ## ... Procrustes: rmse 0.0001793413  max resid 0.0004066026 
    ## ... Similar to previous best
    ## Run 131 stress 9.633652e-05 
    ## ... Procrustes: rmse 0.0001482523  max resid 0.0003251693 
    ## ... Similar to previous best
    ## Run 132 stress 9.248126e-05 
    ## ... Procrustes: rmse 0.0001422785  max resid 0.0003183772 
    ## ... Similar to previous best
    ## Run 133 stress 9.639767e-05 
    ## ... Procrustes: rmse 0.0001627147  max resid 0.0003523129 
    ## ... Similar to previous best
    ## Run 134 stress 9.635751e-05 
    ## ... Procrustes: rmse 9.005434e-05  max resid 0.0002115561 
    ## ... Similar to previous best
    ## Run 135 stress 9.913346e-05 
    ## ... Procrustes: rmse 0.000130065  max resid 0.000215534 
    ## ... Similar to previous best
    ## Run 136 stress 9.966007e-05 
    ## ... Procrustes: rmse 0.0001830034  max resid 0.0002911348 
    ## ... Similar to previous best
    ## Run 137 stress 9.646866e-05 
    ## ... Procrustes: rmse 0.000183286  max resid 0.0004189581 
    ## ... Similar to previous best
    ## Run 138 stress 9.269151e-05 
    ## ... Procrustes: rmse 8.238981e-05  max resid 0.000138932 
    ## ... Similar to previous best
    ## Run 139 stress 9.352504e-05 
    ## ... Procrustes: rmse 9.989764e-05  max resid 0.0002033471 
    ## ... Similar to previous best
    ## Run 140 stress 9.970433e-05 
    ## ... Procrustes: rmse 0.000185178  max resid 0.0004294725 
    ## ... Similar to previous best
    ## Run 141 stress 8.353403e-05 
    ## ... Procrustes: rmse 0.0001593478  max resid 0.000294193 
    ## ... Similar to previous best
    ## Run 142 stress 9.32179e-05 
    ## ... Procrustes: rmse 0.0001597856  max resid 0.000301388 
    ## ... Similar to previous best
    ## Run 143 stress 9.8668e-05 
    ## ... Procrustes: rmse 7.430983e-05  max resid 0.0001042111 
    ## ... Similar to previous best
    ## Run 144 stress 9.875401e-05 
    ## ... Procrustes: rmse 0.0001563163  max resid 0.0003229152 
    ## ... Similar to previous best
    ## Run 145 stress 0.2974629 
    ## Run 146 stress 9.757191e-05 
    ## ... Procrustes: rmse 7.430164e-05  max resid 0.0001138905 
    ## ... Similar to previous best
    ## Run 147 stress 9.316177e-05 
    ## ... Procrustes: rmse 0.00016169  max resid 0.0003033171 
    ## ... Similar to previous best
    ## Run 148 stress 9.270555e-05 
    ## ... Procrustes: rmse 0.0001616941  max resid 0.0003020037 
    ## ... Similar to previous best
    ## Run 149 stress 9.528847e-05 
    ## ... Procrustes: rmse 0.000187122  max resid 0.0003210144 
    ## ... Similar to previous best
    ## Run 150 stress 9.163091e-05 
    ## ... Procrustes: rmse 0.0001409968  max resid 0.0003154652 
    ## ... Similar to previous best
    ## Run 151 stress 9.767068e-05 
    ## ... Procrustes: rmse 0.0001682789  max resid 0.0003103321 
    ## ... Similar to previous best
    ## Run 152 stress 9.95777e-05 
    ## ... Procrustes: rmse 0.0001684034  max resid 0.0003441681 
    ## ... Similar to previous best
    ## Run 153 stress 9.576208e-05 
    ## ... Procrustes: rmse 0.0001385092  max resid 0.0002983154 
    ## ... Similar to previous best
    ## Run 154 stress 8.856343e-05 
    ## ... Procrustes: rmse 0.000158051  max resid 0.0003417705 
    ## ... Similar to previous best
    ## Run 155 stress 9.788808e-05 
    ## ... Procrustes: rmse 0.0001915835  max resid 0.0002990281 
    ## ... Similar to previous best
    ## Run 156 stress 9.611931e-05 
    ## ... Procrustes: rmse 0.0001411729  max resid 0.0003268855 
    ## ... Similar to previous best
    ## Run 157 stress 9.360955e-05 
    ## ... Procrustes: rmse 8.441186e-05  max resid 0.0001393376 
    ## ... Similar to previous best
    ## Run 158 stress 9.681266e-05 
    ## ... Procrustes: rmse 9.131312e-05  max resid 0.0001509915 
    ## ... Similar to previous best
    ## Run 159 stress 9.38787e-05 
    ## ... Procrustes: rmse 0.0001750858  max resid 0.0004002618 
    ## ... Similar to previous best
    ## Run 160 stress 9.352538e-05 
    ## ... Procrustes: rmse 0.0001161681  max resid 0.0002845526 
    ## ... Similar to previous best
    ## Run 161 stress 9.125644e-05 
    ## ... Procrustes: rmse 5.59091e-05  max resid 8.488013e-05 
    ## ... Similar to previous best
    ## Run 162 stress 7.694886e-05 
    ## ... Procrustes: rmse 0.0001367322  max resid 0.0002499475 
    ## ... Similar to previous best
    ## Run 163 stress 8.868341e-05 
    ## ... Procrustes: rmse 0.0001696743  max resid 0.0003959255 
    ## ... Similar to previous best
    ## Run 164 stress 9.81738e-05 
    ## ... Procrustes: rmse 0.0001527772  max resid 0.0003045348 
    ## ... Similar to previous best
    ## Run 165 stress 9.680254e-05 
    ## ... Procrustes: rmse 8.528836e-05  max resid 0.0001250055 
    ## ... Similar to previous best
    ## Run 166 stress 9.649906e-05 
    ## ... Procrustes: rmse 0.0001507578  max resid 0.0002944058 
    ## ... Similar to previous best
    ## Run 167 stress 9.80466e-05 
    ## ... Procrustes: rmse 0.0001094988  max resid 0.0001713727 
    ## ... Similar to previous best
    ## Run 168 stress 9.592993e-05 
    ## ... Procrustes: rmse 0.0001657715  max resid 0.0003099598 
    ## ... Similar to previous best
    ## Run 169 stress 9.284039e-05 
    ## ... Procrustes: rmse 0.0001693771  max resid 0.0003695523 
    ## ... Similar to previous best
    ## Run 170 stress 8.710353e-05 
    ## ... Procrustes: rmse 0.0001540894  max resid 0.0002911977 
    ## ... Similar to previous best
    ## Run 171 stress 9.229439e-05 
    ## ... Procrustes: rmse 0.0001590766  max resid 0.0002942026 
    ## ... Similar to previous best
    ## Run 172 stress 9.504757e-05 
    ## ... Procrustes: rmse 0.0001350933  max resid 0.0002284745 
    ## ... Similar to previous best
    ## Run 173 stress 9.561774e-05 
    ## ... Procrustes: rmse 0.0001765736  max resid 0.0004096422 
    ## ... Similar to previous best
    ## Run 174 stress 9.910164e-05 
    ## ... Procrustes: rmse 0.0001922162  max resid 0.0004094071 
    ## ... Similar to previous best
    ## Run 175 stress 9.372012e-05 
    ## ... Procrustes: rmse 0.0001563309  max resid 0.0003304887 
    ## ... Similar to previous best
    ## Run 176 stress 8.516679e-05 
    ## ... Procrustes: rmse 0.0001514453  max resid 0.0002899135 
    ## ... Similar to previous best
    ## Run 177 stress 9.59018e-05 
    ## ... Procrustes: rmse 0.0001840194  max resid 0.0004191039 
    ## ... Similar to previous best
    ## Run 178 stress 8.972725e-05 
    ## ... Procrustes: rmse 0.000121882  max resid 0.0002848959 
    ## ... Similar to previous best
    ## Run 179 stress 9.788519e-05 
    ## ... Procrustes: rmse 0.0001695916  max resid 0.0003118437 
    ## ... Similar to previous best
    ## Run 180 stress 9.852976e-05 
    ## ... Procrustes: rmse 0.000166475  max resid 0.0003198635 
    ## ... Similar to previous best
    ## Run 181 stress 9.91601e-05 
    ## ... Procrustes: rmse 7.755555e-05  max resid 0.0001171188 
    ## ... Similar to previous best
    ## Run 182 stress 8.826303e-05 
    ## ... Procrustes: rmse 7.034397e-05  max resid 0.0001096118 
    ## ... Similar to previous best
    ## Run 183 stress 9.603179e-05 
    ## ... Procrustes: rmse 0.0001732138  max resid 0.0003944559 
    ## ... Similar to previous best
    ## Run 184 stress 9.724019e-05 
    ## ... Procrustes: rmse 0.0001972916  max resid 0.0003349334 
    ## ... Similar to previous best
    ## Run 185 stress 9.995321e-05 
    ## ... Procrustes: rmse 0.000202618  max resid 0.0003411748 
    ## ... Similar to previous best
    ## Run 186 stress 9.404085e-05 
    ## ... Procrustes: rmse 0.0001797671  max resid 0.0004111226 
    ## ... Similar to previous best
    ## Run 187 stress 9.964142e-05 
    ## ... Procrustes: rmse 0.0002011344  max resid 0.0003391248 
    ## ... Similar to previous best
    ## Run 188 stress 9.660463e-05 
    ## ... Procrustes: rmse 0.0001644275  max resid 0.0003405204 
    ## ... Similar to previous best
    ## Run 189 stress 9.813512e-05 
    ## ... Procrustes: rmse 0.0001650768  max resid 0.0003103408 
    ## ... Similar to previous best
    ## Run 190 stress 9.92062e-05 
    ## ... Procrustes: rmse 0.0001259788  max resid 0.0001937848 
    ## ... Similar to previous best
    ## Run 191 stress 9.817561e-05 
    ## ... Procrustes: rmse 0.000113417  max resid 0.000174046 
    ## ... Similar to previous best
    ## Run 192 stress 9.653432e-05 
    ## ... Procrustes: rmse 6.344753e-05  max resid 0.0001029728 
    ## ... Similar to previous best
    ## Run 193 stress 9.416141e-05 
    ## ... Procrustes: rmse 0.0001373441  max resid 0.0002705851 
    ## ... Similar to previous best
    ## Run 194 stress 9.621504e-05 
    ## ... Procrustes: rmse 0.0001086263  max resid 0.0001687505 
    ## ... Similar to previous best
    ## Run 195 stress 9.952336e-05 
    ## ... Procrustes: rmse 0.0001866956  max resid 0.0004244783 
    ## ... Similar to previous best
    ## Run 196 stress 9.530513e-05 
    ## ... Procrustes: rmse 0.0001820933  max resid 0.0003266538 
    ## ... Similar to previous best
    ## Run 197 stress 8.892978e-05 
    ## ... Procrustes: rmse 0.0001406224  max resid 0.0003036704 
    ## ... Similar to previous best
    ## Run 198 stress 8.681882e-05 
    ## ... Procrustes: rmse 0.0001617838  max resid 0.0003815478 
    ## ... Similar to previous best
    ## Run 199 stress 9.925682e-05 
    ## ... Procrustes: rmse 0.0001388448  max resid 0.0002792847 
    ## ... Similar to previous best
    ## Run 200 stress 9.791939e-05 
    ## ... Procrustes: rmse 0.000185476  max resid 0.0003308504 
    ## ... Similar to previous best
    ## Run 201 stress 9.935642e-05 
    ## ... Procrustes: rmse 0.0001414041  max resid 0.0002391462 
    ## ... Similar to previous best
    ## Run 202 stress 9.831908e-05 
    ## ... Procrustes: rmse 0.0001396088  max resid 0.0002851728 
    ## ... Similar to previous best
    ## Run 203 stress 8.61748e-05 
    ## ... Procrustes: rmse 0.0001581693  max resid 0.0003705151 
    ## ... Similar to previous best
    ## Run 204 stress 9.046866e-05 
    ## ... Procrustes: rmse 0.0001850367  max resid 0.0003214431 
    ## ... Similar to previous best
    ## Run 205 stress 9.197267e-05 
    ## ... Procrustes: rmse 0.0001509628  max resid 0.0002938688 
    ## ... Similar to previous best
    ## Run 206 stress 9.342399e-05 
    ## ... Procrustes: rmse 0.0001763309  max resid 0.0004048346 
    ## ... Similar to previous best
    ## Run 207 stress 9.245877e-05 
    ## ... Procrustes: rmse 0.000135098  max resid 0.0002765906 
    ## ... Similar to previous best
    ## Run 208 stress 9.644792e-05 
    ## ... Procrustes: rmse 0.0001476127  max resid 0.0003244077 
    ## ... Similar to previous best
    ## Run 209 stress 9.722814e-05 
    ## ... Procrustes: rmse 0.00016788  max resid 0.0003123565 
    ## ... Similar to previous best
    ## Run 210 stress 9.679708e-05 
    ## ... Procrustes: rmse 7.477553e-05  max resid 0.0001173926 
    ## ... Similar to previous best
    ## Run 211 stress 9.64324e-05 
    ## ... Procrustes: rmse 9.022327e-05  max resid 0.0001491569 
    ## ... Similar to previous best
    ## Run 212 stress 9.891189e-05 
    ## ... Procrustes: rmse 0.0001848413  max resid 0.0004204032 
    ## ... Similar to previous best
    ## Run 213 stress 9.138149e-05 
    ## ... Procrustes: rmse 0.0001108404  max resid 0.0001918412 
    ## ... Similar to previous best
    ## Run 214 stress 9.047477e-05 
    ## ... Procrustes: rmse 7.077569e-05  max resid 0.0001141427 
    ## ... Similar to previous best
    ## Run 215 stress 9.799317e-05 
    ## ... Procrustes: rmse 7.398213e-05  max resid 0.0001143834 
    ## ... Similar to previous best
    ## Run 216 stress 9.509092e-05 
    ## ... Procrustes: rmse 0.0001788996  max resid 0.0003836287 
    ## ... Similar to previous best
    ## Run 217 stress 9.84758e-05 
    ## ... Procrustes: rmse 0.0001641777  max resid 0.0003137045 
    ## ... Similar to previous best
    ## Run 218 stress 9.819921e-05 
    ## ... Procrustes: rmse 0.0001865208  max resid 0.0003322824 
    ## ... Similar to previous best
    ## Run 219 stress 9.978654e-05 
    ## ... Procrustes: rmse 9.322458e-05  max resid 0.0001545464 
    ## ... Similar to previous best
    ## Run 220 stress 9.597313e-05 
    ## ... Procrustes: rmse 8.941151e-05  max resid 0.0001486505 
    ## ... Similar to previous best
    ## Run 221 stress 9.910129e-05 
    ## ... Procrustes: rmse 8.875388e-05  max resid 0.0001288347 
    ## ... Similar to previous best
    ## Run 222 stress 0.2582637 
    ## Run 223 stress 9.899859e-05 
    ## ... Procrustes: rmse 6.71936e-05  max resid 0.000124368 
    ## ... Similar to previous best
    ## Run 224 stress 9.081414e-05 
    ## ... Procrustes: rmse 0.0001653712  max resid 0.0003905145 
    ## ... Similar to previous best
    ## Run 225 stress 9.502267e-05 
    ## ... Procrustes: rmse 0.0001353862  max resid 0.0002296198 
    ## ... Similar to previous best
    ## Run 226 stress 9.858709e-05 
    ## ... Procrustes: rmse 0.0001775374  max resid 0.0003890988 
    ## ... Similar to previous best
    ## Run 227 stress 9.614726e-05 
    ## ... Procrustes: rmse 8.957361e-05  max resid 0.0001481193 
    ## ... Similar to previous best
    ## Run 228 stress 4.527529e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001446156  max resid 0.0002236024 
    ## ... Similar to previous best
    ## Run 229 stress 9.300499e-05 
    ## ... Procrustes: rmse 0.0001922366  max resid 0.0002807034 
    ## ... Similar to previous best
    ## Run 230 stress 9.834117e-05 
    ## ... Procrustes: rmse 0.0002191712  max resid 0.0003010396 
    ## ... Similar to previous best
    ## Run 231 stress 8.721007e-05 
    ## ... Procrustes: rmse 0.0001094914  max resid 0.0001845181 
    ## ... Similar to previous best
    ## Run 232 stress 8.194554e-05 
    ## ... Procrustes: rmse 0.0001626929  max resid 0.0002485564 
    ## ... Similar to previous best
    ## Run 233 stress 9.602738e-05 
    ## ... Procrustes: rmse 0.0001939794  max resid 0.0002619648 
    ## ... Similar to previous best
    ## Run 234 stress 8.91316e-05 
    ## ... Procrustes: rmse 0.0001378373  max resid 0.0001874711 
    ## ... Similar to previous best
    ## Run 235 stress 0.2386316 
    ## Run 236 stress 8.739419e-05 
    ## ... Procrustes: rmse 0.0001669794  max resid 0.0002363957 
    ## ... Similar to previous best
    ## Run 237 stress 9.712319e-05 
    ## ... Procrustes: rmse 0.0001985787  max resid 0.0002865217 
    ## ... Similar to previous best
    ## Run 238 stress 9.792845e-05 
    ## ... Procrustes: rmse 0.000133982  max resid 0.000210861 
    ## ... Similar to previous best
    ## Run 239 stress 9.869215e-05 
    ## ... Procrustes: rmse 0.0002028994  max resid 0.000273566 
    ## ... Similar to previous best
    ## Run 240 stress 9.968939e-05 
    ## ... Procrustes: rmse 0.0002148568  max resid 0.0002853396 
    ## ... Similar to previous best
    ## Run 241 stress 9.679532e-05 
    ## ... Procrustes: rmse 0.000209626  max resid 0.0002977195 
    ## ... Similar to previous best
    ## Run 242 stress 9.064749e-05 
    ## ... Procrustes: rmse 0.0001610248  max resid 0.0002555205 
    ## ... Similar to previous best
    ## Run 243 stress 9.462543e-05 
    ## ... Procrustes: rmse 0.0002063579  max resid 0.0002977682 
    ## ... Similar to previous best
    ## Run 244 stress 9.779105e-05 
    ## ... Procrustes: rmse 0.0001302944  max resid 0.0001965203 
    ## ... Similar to previous best
    ## Run 245 stress 9.330487e-05 
    ## ... Procrustes: rmse 0.0001761121  max resid 0.0002136503 
    ## ... Similar to previous best
    ## Run 246 stress 9.485041e-05 
    ## ... Procrustes: rmse 0.0001914652  max resid 0.0002575373 
    ## ... Similar to previous best
    ## Run 247 stress 8.928367e-05 
    ## ... Procrustes: rmse 0.0001577601  max resid 0.0002021973 
    ## ... Similar to previous best
    ## Run 248 stress 9.303298e-05 
    ## ... Procrustes: rmse 0.0001523952  max resid 0.000211409 
    ## ... Similar to previous best
    ## Run 249 stress 9.76661e-05 
    ## ... Procrustes: rmse 9.401662e-05  max resid 0.0001358874 
    ## ... Similar to previous best
    ## Run 250 stress 9.266361e-05 
    ## ... Procrustes: rmse 0.0001571714  max resid 0.00020263 
    ## ... Similar to previous best
    ## Run 251 stress 9.294278e-05 
    ## ... Procrustes: rmse 0.0001610776  max resid 0.0002077154 
    ## ... Similar to previous best
    ## Run 252 stress 9.297652e-05 
    ## ... Procrustes: rmse 0.0001895979  max resid 0.0002682442 
    ## ... Similar to previous best
    ## Run 253 stress 9.451207e-05 
    ## ... Procrustes: rmse 0.0001690505  max resid 0.0002682223 
    ## ... Similar to previous best
    ## Run 254 stress 9.983517e-05 
    ## ... Procrustes: rmse 0.0001903864  max resid 0.0002282318 
    ## ... Similar to previous best
    ## Run 255 stress 9.165034e-05 
    ## ... Procrustes: rmse 0.0001594198  max resid 0.0002107345 
    ## ... Similar to previous best
    ## Run 256 stress 8.994275e-05 
    ## ... Procrustes: rmse 0.0001539972  max resid 0.0001949526 
    ## ... Similar to previous best
    ## Run 257 stress 9.227559e-05 
    ## ... Procrustes: rmse 0.0001965632  max resid 0.0002790587 
    ## ... Similar to previous best
    ## Run 258 stress 9.272514e-05 
    ## ... Procrustes: rmse 0.0001212285  max resid 0.0001880945 
    ## ... Similar to previous best
    ## Run 259 stress 9.74664e-05 
    ## ... Procrustes: rmse 0.0001938004  max resid 0.0002731696 
    ## ... Similar to previous best
    ## Run 260 stress 8.907237e-05 
    ## ... Procrustes: rmse 0.0001370396  max resid 0.0002318804 
    ## ... Similar to previous best
    ## Run 261 stress 9.895223e-05 
    ## ... Procrustes: rmse 0.0001898689  max resid 0.0002649958 
    ## ... Similar to previous best
    ## Run 262 stress 9.840081e-05 
    ## ... Procrustes: rmse 0.000199164  max resid 0.000288868 
    ## ... Similar to previous best
    ## Run 263 stress 9.969799e-05 
    ## ... Procrustes: rmse 0.0002043408  max resid 0.0003007306 
    ## ... Similar to previous best
    ## Run 264 stress 9.352818e-05 
    ## ... Procrustes: rmse 0.0001661713  max resid 0.0002374312 
    ## ... Similar to previous best
    ## Run 265 stress 9.620857e-05 
    ## ... Procrustes: rmse 0.0001655439  max resid 0.0002168506 
    ## ... Similar to previous best
    ## Run 266 stress 8.776962e-05 
    ## ... Procrustes: rmse 0.0001583892  max resid 0.0002510595 
    ## ... Similar to previous best
    ## Run 267 stress 9.925539e-05 
    ## ... Procrustes: rmse 0.0002084361  max resid 0.0003136353 
    ## ... Similar to previous best
    ## Run 268 stress 9.829959e-05 
    ## ... Procrustes: rmse 0.0002157577  max resid 0.0002953498 
    ## ... Similar to previous best
    ## Run 269 stress 9.536688e-05 
    ## ... Procrustes: rmse 0.0002115042  max resid 0.0003066025 
    ## ... Similar to previous best
    ## Run 270 stress 9.441555e-05 
    ## ... Procrustes: rmse 0.0001900742  max resid 0.0002844727 
    ## ... Similar to previous best
    ## Run 271 stress 9.971501e-05 
    ## ... Procrustes: rmse 0.0002102335  max resid 0.0002973912 
    ## ... Similar to previous best
    ## Run 272 stress 9.884305e-05 
    ## ... Procrustes: rmse 0.0002046785  max resid 0.0002966338 
    ## ... Similar to previous best
    ## Run 273 stress 9.160694e-05 
    ## ... Procrustes: rmse 0.0001618102  max resid 0.0002096842 
    ## ... Similar to previous best
    ## Run 274 stress 9.893416e-05 
    ## ... Procrustes: rmse 0.0002222726  max resid 0.0003185853 
    ## ... Similar to previous best
    ## Run 275 stress 8.513557e-05 
    ## ... Procrustes: rmse 0.0001730424  max resid 0.0002543328 
    ## ... Similar to previous best
    ## Run 276 stress 9.923753e-05 
    ## ... Procrustes: rmse 0.0002072803  max resid 0.0003023802 
    ## ... Similar to previous best
    ## Run 277 stress 9.225984e-05 
    ## ... Procrustes: rmse 0.0001639103  max resid 0.0002190977 
    ## ... Similar to previous best
    ## Run 278 stress 9.091083e-05 
    ## ... Procrustes: rmse 0.000182432  max resid 0.0002698132 
    ## ... Similar to previous best
    ## Run 279 stress 9.866144e-05 
    ## ... Procrustes: rmse 0.0001538793  max resid 0.0002399182 
    ## ... Similar to previous best
    ## Run 280 stress 9.669007e-05 
    ## ... Procrustes: rmse 0.0002157084  max resid 0.0003055944 
    ## ... Similar to previous best
    ## Run 281 stress 9.540892e-05 
    ## ... Procrustes: rmse 0.0002130456  max resid 0.0002999259 
    ## ... Similar to previous best
    ## Run 282 stress 9.885286e-05 
    ## ... Procrustes: rmse 0.0001827492  max resid 0.0002634788 
    ## ... Similar to previous best
    ## Run 283 stress 9.92063e-05 
    ## ... Procrustes: rmse 0.0002054979  max resid 0.0002830468 
    ## ... Similar to previous best
    ## Run 284 stress 9.762872e-05 
    ## ... Procrustes: rmse 0.0001814421  max resid 0.0002162912 
    ## ... Similar to previous best
    ## Run 285 stress 9.507658e-05 
    ## ... Procrustes: rmse 0.0001622011  max resid 0.0002355329 
    ## ... Similar to previous best
    ## Run 286 stress 9.498649e-05 
    ## ... Procrustes: rmse 0.000163344  max resid 0.000227737 
    ## ... Similar to previous best
    ## Run 287 stress 8.481274e-05 
    ## ... Procrustes: rmse 0.0001421399  max resid 0.0001821303 
    ## ... Similar to previous best
    ## Run 288 stress 8.357014e-05 
    ## ... Procrustes: rmse 0.0001364682  max resid 0.0001990845 
    ## ... Similar to previous best
    ## Run 289 stress 9.438908e-05 
    ## ... Procrustes: rmse 0.000164079  max resid 0.0002536482 
    ## ... Similar to previous best
    ## Run 290 stress 9.877586e-05 
    ## ... Procrustes: rmse 0.0001748781  max resid 0.0002248069 
    ## ... Similar to previous best
    ## Run 291 stress 9.591911e-05 
    ## ... Procrustes: rmse 0.0002042163  max resid 0.0002731436 
    ## ... Similar to previous best
    ## Run 292 stress 9.430795e-05 
    ## ... Procrustes: rmse 0.0001948095  max resid 0.0002842261 
    ## ... Similar to previous best
    ## Run 293 stress 9.383449e-05 
    ## ... Procrustes: rmse 0.0001257941  max resid 0.0001652476 
    ## ... Similar to previous best
    ## Run 294 stress 8.612421e-05 
    ## ... Procrustes: rmse 0.0001545945  max resid 0.0002417148 
    ## ... Similar to previous best
    ## Run 295 stress 9.412865e-05 
    ## ... Procrustes: rmse 0.0001617579  max resid 0.0002110959 
    ## ... Similar to previous best
    ## Run 296 stress 9.917323e-05 
    ## ... Procrustes: rmse 0.0002101638  max resid 0.0002867361 
    ## ... Similar to previous best
    ## Run 297 stress 8.076948e-05 
    ## ... Procrustes: rmse 9.482469e-05  max resid 0.0001688247 
    ## ... Similar to previous best
    ## Run 298 stress 9.185354e-05 
    ## ... Procrustes: rmse 0.0001825158  max resid 0.0002733124 
    ## ... Similar to previous best
    ## Run 299 stress 9.912885e-05 
    ## ... Procrustes: rmse 0.0001978652  max resid 0.0003024531 
    ## ... Similar to previous best
    ## Run 300 stress 9.991277e-05 
    ## ... Procrustes: rmse 0.000187573  max resid 0.0002691507 
    ## ... Similar to previous best
    ## Run 301 stress 9.477258e-05 
    ## ... Procrustes: rmse 0.0002008606  max resid 0.0002814353 
    ## ... Similar to previous best
    ## Run 302 stress 9.630167e-05 
    ## ... Procrustes: rmse 0.000200589  max resid 0.0002923077 
    ## ... Similar to previous best
    ## Run 303 stress 8.82635e-05 
    ## ... Procrustes: rmse 0.0001540601  max resid 0.0002197842 
    ## ... Similar to previous best
    ## Run 304 stress 9.778541e-05 
    ## ... Procrustes: rmse 0.0001725682  max resid 0.0002210282 
    ## ... Similar to previous best
    ## Run 305 stress 9.78178e-05 
    ## ... Procrustes: rmse 0.0002076395  max resid 0.0002796708 
    ## ... Similar to previous best
    ## Run 306 stress 6.846612e-05 
    ## ... Procrustes: rmse 8.909524e-05  max resid 0.000124052 
    ## ... Similar to previous best
    ## Run 307 stress 9.362516e-05 
    ## ... Procrustes: rmse 0.0001488293  max resid 0.0001954461 
    ## ... Similar to previous best
    ## Run 308 stress 9.825843e-05 
    ## ... Procrustes: rmse 0.0001620255  max resid 0.0002195744 
    ## ... Similar to previous best
    ## Run 309 stress 9.942921e-05 
    ## ... Procrustes: rmse 0.0002232644  max resid 0.0003150786 
    ## ... Similar to previous best
    ## Run 310 stress 9.905366e-05 
    ## ... Procrustes: rmse 0.0002036731  max resid 0.0002982011 
    ## ... Similar to previous best
    ## Run 311 stress 9.451676e-05 
    ## ... Procrustes: rmse 0.0001651728  max resid 0.0002105655 
    ## ... Similar to previous best
    ## Run 312 stress 9.713532e-05 
    ## ... Procrustes: rmse 0.0001656083  max resid 0.0002247788 
    ## ... Similar to previous best
    ## Run 313 stress 9.743924e-05 
    ## ... Procrustes: rmse 0.0002171408  max resid 0.0002967512 
    ## ... Similar to previous best
    ## Run 314 stress 9.270726e-05 
    ## ... Procrustes: rmse 0.0001244967  max resid 0.0001883932 
    ## ... Similar to previous best
    ## Run 315 stress 9.59095e-05 
    ## ... Procrustes: rmse 0.0001677345  max resid 0.0002155987 
    ## ... Similar to previous best
    ## Run 316 stress 9.573093e-05 
    ## ... Procrustes: rmse 0.0001669424  max resid 0.0002128248 
    ## ... Similar to previous best
    ## Run 317 stress 9.892722e-05 
    ## ... Procrustes: rmse 0.00022099  max resid 0.0003034552 
    ## ... Similar to previous best
    ## Run 318 stress 9.615422e-05 
    ## ... Procrustes: rmse 0.0001999224  max resid 0.0002909052 
    ## ... Similar to previous best
    ## Run 319 stress 8.985256e-05 
    ## ... Procrustes: rmse 0.000140907  max resid 0.0002107931 
    ## ... Similar to previous best
    ## Run 320 stress 7.828654e-05 
    ## ... Procrustes: rmse 0.0001545089  max resid 0.0002437452 
    ## ... Similar to previous best
    ## Run 321 stress 9.98954e-05 
    ## ... Procrustes: rmse 0.0001852019  max resid 0.0002969143 
    ## ... Similar to previous best
    ## Run 322 stress 9.34208e-05 
    ## ... Procrustes: rmse 0.0001908242  max resid 0.000266642 
    ## ... Similar to previous best
    ## Run 323 stress 9.083031e-05 
    ## ... Procrustes: rmse 0.0002016715  max resid 0.0002845648 
    ## ... Similar to previous best
    ## Run 324 stress 9.812856e-05 
    ## ... Procrustes: rmse 0.0001857983  max resid 0.0002250835 
    ## ... Similar to previous best
    ## Run 325 stress 9.159478e-05 
    ## ... Procrustes: rmse 0.0001864358  max resid 0.0002744735 
    ## ... Similar to previous best
    ## Run 326 stress 9.40722e-05 
    ## ... Procrustes: rmse 0.0002009128  max resid 0.0002726434 
    ## ... Similar to previous best
    ## Run 327 stress 9.053074e-05 
    ## ... Procrustes: rmse 0.0001813277  max resid 0.0002492131 
    ## ... Similar to previous best
    ## Run 328 stress 9.284527e-05 
    ## ... Procrustes: rmse 0.0001610873  max resid 0.000206041 
    ## ... Similar to previous best
    ## Run 329 stress 8.801981e-05 
    ## ... Procrustes: rmse 0.0001499078  max resid 0.0001906714 
    ## ... Similar to previous best
    ## Run 330 stress 9.1869e-05 
    ## ... Procrustes: rmse 0.0001606909  max resid 0.0002161628 
    ## ... Similar to previous best
    ## Run 331 stress 9.434258e-05 
    ## ... Procrustes: rmse 0.0002064664  max resid 0.0002845534 
    ## ... Similar to previous best
    ## Run 332 stress 9.990704e-05 
    ## ... Procrustes: rmse 0.0001772427  max resid 0.0002266035 
    ## ... Similar to previous best
    ## Run 333 stress 9.902361e-05 
    ## ... Procrustes: rmse 0.0002101536  max resid 0.000282194 
    ## ... Similar to previous best
    ## Run 334 stress 9.620418e-05 
    ## ... Procrustes: rmse 0.0001746117  max resid 0.0002646137 
    ## ... Similar to previous best
    ## Run 335 stress 8.489351e-05 
    ## ... Procrustes: rmse 0.0001173833  max resid 0.0001533821 
    ## ... Similar to previous best
    ## Run 336 stress 9.839956e-05 
    ## ... Procrustes: rmse 0.0001846657  max resid 0.0002666111 
    ## ... Similar to previous best
    ## Run 337 stress 9.846832e-05 
    ## ... Procrustes: rmse 0.0001998271  max resid 0.0002921348 
    ## ... Similar to previous best
    ## Run 338 stress 9.444372e-05 
    ## ... Procrustes: rmse 0.000129729  max resid 0.0001810199 
    ## ... Similar to previous best
    ## Run 339 stress 9.826651e-05 
    ## ... Procrustes: rmse 0.0001874883  max resid 0.0002650981 
    ## ... Similar to previous best
    ## Run 340 stress 8.775411e-05 
    ## ... Procrustes: rmse 0.0001754055  max resid 0.0002381415 
    ## ... Similar to previous best
    ## Run 341 stress 9.143338e-05 
    ## ... Procrustes: rmse 0.0001252835  max resid 0.0001767181 
    ## ... Similar to previous best
    ## Run 342 stress 9.243438e-05 
    ## ... Procrustes: rmse 0.0002047193  max resid 0.0002789557 
    ## ... Similar to previous best
    ## Run 343 stress 9.930679e-05 
    ## ... Procrustes: rmse 0.0002087667  max resid 0.0003041292 
    ## ... Similar to previous best
    ## Run 344 stress 9.635324e-05 
    ## ... Procrustes: rmse 0.0002020852  max resid 0.0003070511 
    ## ... Similar to previous best
    ## Run 345 stress 9.035613e-05 
    ## ... Procrustes: rmse 9.117243e-05  max resid 0.0001583623 
    ## ... Similar to previous best
    ## Run 346 stress 9.541585e-05 
    ## ... Procrustes: rmse 0.0001945493  max resid 0.000259515 
    ## ... Similar to previous best
    ## Run 347 stress 9.682469e-05 
    ## ... Procrustes: rmse 0.0001967855  max resid 0.000294828 
    ## ... Similar to previous best
    ## Run 348 stress 9.777234e-05 
    ## ... Procrustes: rmse 0.0001931203  max resid 0.0002727847 
    ## ... Similar to previous best
    ## Run 349 stress 9.132005e-05 
    ## ... Procrustes: rmse 0.0001852433  max resid 0.000285475 
    ## ... Similar to previous best
    ## Run 350 stress 9.570472e-05 
    ## ... Procrustes: rmse 0.0002074498  max resid 0.0003042197 
    ## ... Similar to previous best
    ## Run 351 stress 9.959942e-05 
    ## ... Procrustes: rmse 0.0001573931  max resid 0.000243231 
    ## ... Similar to previous best
    ## Run 352 stress 9.702057e-05 
    ## ... Procrustes: rmse 0.0001541087  max resid 0.0002117768 
    ## ... Similar to previous best
    ## Run 353 stress 9.884058e-05 
    ## ... Procrustes: rmse 0.0002002422  max resid 0.0002696414 
    ## ... Similar to previous best
    ## Run 354 stress 9.535879e-05 
    ## ... Procrustes: rmse 0.000159132  max resid 0.0002047222 
    ## ... Similar to previous best
    ## Run 355 stress 8.707028e-05 
    ## ... Procrustes: rmse 0.0001430405  max resid 0.000180698 
    ## ... Similar to previous best
    ## Run 356 stress 9.007287e-05 
    ## ... Procrustes: rmse 0.0001811365  max resid 0.0002616127 
    ## ... Similar to previous best
    ## Run 357 stress 9.870264e-05 
    ## ... Procrustes: rmse 0.0001883618  max resid 0.0002238895 
    ## ... Similar to previous best
    ## Run 358 stress 9.66223e-05 
    ## ... Procrustes: rmse 0.0001960941  max resid 0.0002986398 
    ## ... Similar to previous best
    ## Run 359 stress 9.41267e-05 
    ## ... Procrustes: rmse 0.0001608629  max resid 0.0002092103 
    ## ... Similar to previous best
    ## Run 360 stress 7.319723e-05 
    ## ... Procrustes: rmse 0.0001416733  max resid 0.0002274882 
    ## ... Similar to previous best
    ## Run 361 stress 9.574478e-05 
    ## ... Procrustes: rmse 0.0001946762  max resid 0.0002706987 
    ## ... Similar to previous best
    ## Run 362 stress 9.55142e-05 
    ## ... Procrustes: rmse 0.0001285623  max resid 0.0001808458 
    ## ... Similar to previous best
    ## Run 363 stress 9.003463e-05 
    ## ... Procrustes: rmse 0.0001822183  max resid 0.000262147 
    ## ... Similar to previous best
    ## Run 364 stress 9.968852e-05 
    ## ... Procrustes: rmse 0.0002037421  max resid 0.0002983343 
    ## ... Similar to previous best
    ## Run 365 stress 9.578485e-05 
    ## ... Procrustes: rmse 0.0001703965  max resid 0.0002526766 
    ## ... Similar to previous best
    ## Run 366 stress 9.545706e-05 
    ## ... Procrustes: rmse 0.0001923144  max resid 0.000271699 
    ## ... Similar to previous best
    ## Run 367 stress 9.074382e-05 
    ## ... Procrustes: rmse 0.0001799761  max resid 0.0002475841 
    ## ... Similar to previous best
    ## Run 368 stress 9.876056e-05 
    ## ... Procrustes: rmse 0.0002079485  max resid 0.0003160313 
    ## ... Similar to previous best
    ## Run 369 stress 9.121546e-05 
    ## ... Procrustes: rmse 0.0001544587  max resid 0.0002426496 
    ## ... Similar to previous best
    ## Run 370 stress 9.714666e-05 
    ## ... Procrustes: rmse 0.0001705454  max resid 0.0002190039 
    ## ... Similar to previous best
    ## Run 371 stress 9.979851e-05 
    ## ... Procrustes: rmse 0.0002051693  max resid 0.0003159257 
    ## ... Similar to previous best
    ## Run 372 stress 9.78578e-05 
    ## ... Procrustes: rmse 0.000205126  max resid 0.0002998307 
    ## ... Similar to previous best
    ## Run 373 stress 9.955374e-05 
    ## ... Procrustes: rmse 0.0001733356  max resid 0.0002613452 
    ## ... Similar to previous best
    ## Run 374 stress 9.45463e-05 
    ## ... Procrustes: rmse 0.0001935625  max resid 0.000299491 
    ## ... Similar to previous best
    ## Run 375 stress 9.799932e-05 
    ## ... Procrustes: rmse 0.0002028097  max resid 0.0002786448 
    ## ... Similar to previous best
    ## Run 376 stress 9.871139e-05 
    ## ... Procrustes: rmse 0.0001862268  max resid 0.0002959761 
    ## ... Similar to previous best
    ## Run 377 stress 9.845465e-05 
    ## ... Procrustes: rmse 0.0001763688  max resid 0.0002484498 
    ## ... Similar to previous best
    ## Run 378 stress 9.285345e-05 
    ## ... Procrustes: rmse 0.0001614042  max resid 0.0002058438 
    ## ... Similar to previous best
    ## Run 379 stress 9.848429e-05 
    ## ... Procrustes: rmse 0.0001944193  max resid 0.0002981953 
    ## ... Similar to previous best
    ## Run 380 stress 9.751846e-05 
    ## ... Procrustes: rmse 0.0002005314  max resid 0.0002768627 
    ## ... Similar to previous best
    ## Run 381 stress 9.717372e-05 
    ## ... Procrustes: rmse 0.0001831648  max resid 0.0002163477 
    ## ... Similar to previous best
    ## Run 382 stress 9.794632e-05 
    ## ... Procrustes: rmse 0.0001906286  max resid 0.0002692822 
    ## ... Similar to previous best
    ## Run 383 stress 9.478547e-05 
    ## ... Procrustes: rmse 0.0001920128  max resid 0.0002763125 
    ## ... Similar to previous best
    ## Run 384 stress 9.491923e-05 
    ## ... Procrustes: rmse 0.0001320316  max resid 0.0002104374 
    ## ... Similar to previous best
    ## Run 385 stress 9.819523e-05 
    ## ... Procrustes: rmse 0.0001871016  max resid 0.0002623166 
    ## ... Similar to previous best
    ## Run 386 stress 9.42014e-05 
    ## ... Procrustes: rmse 0.0001566449  max resid 0.0002034829 
    ## ... Similar to previous best
    ## Run 387 stress 9.527489e-05 
    ## ... Procrustes: rmse 0.0001915284  max resid 0.0002768001 
    ## ... Similar to previous best
    ## Run 388 stress 9.205102e-05 
    ## ... Procrustes: rmse 0.000153701  max resid 0.0001964131 
    ## ... Similar to previous best
    ## Run 389 stress 8.539255e-05 
    ## ... Procrustes: rmse 0.0001073802  max resid 0.0001594207 
    ## ... Similar to previous best
    ## Run 390 stress 9.980184e-05 
    ## ... Procrustes: rmse 0.0001806621  max resid 0.0002141322 
    ## ... Similar to previous best
    ## Run 391 stress 9.885229e-05 
    ## ... Procrustes: rmse 0.0002164647  max resid 0.0002935451 
    ## ... Similar to previous best
    ## Run 392 stress 9.255174e-05 
    ## ... Procrustes: rmse 0.0001322453  max resid 0.0001876752 
    ## ... Similar to previous best
    ## Run 393 stress 9.336013e-05 
    ## ... Procrustes: rmse 0.0002027734  max resid 0.0002781298 
    ## ... Similar to previous best
    ## Run 394 stress 9.616044e-05 
    ## ... Procrustes: rmse 0.0001970502  max resid 0.0002631571 
    ## ... Similar to previous best
    ## Run 395 stress 9.956875e-05 
    ## ... Procrustes: rmse 0.0002014273  max resid 0.0002926835 
    ## ... Similar to previous best
    ## Run 396 stress 9.696432e-05 
    ## ... Procrustes: rmse 0.0002131345  max resid 0.0002933544 
    ## ... Similar to previous best
    ## Run 397 stress 8.907954e-05 
    ## ... Procrustes: rmse 0.0001524534  max resid 0.0001949273 
    ## ... Similar to previous best
    ## Run 398 stress 9.092804e-05 
    ## ... Procrustes: rmse 0.0001870864  max resid 0.0002736868 
    ## ... Similar to previous best
    ## Run 399 stress 9.165165e-05 
    ## ... Procrustes: rmse 0.0001892837  max resid 0.0002771177 
    ## ... Similar to previous best
    ## Run 400 stress 9.669967e-05 
    ## ... Procrustes: rmse 0.0001962764  max resid 0.0002871735 
    ## ... Similar to previous best
    ## Run 401 stress 9.699622e-05 
    ## ... Procrustes: rmse 0.0001776113  max resid 0.0002691411 
    ## ... Similar to previous best
    ## Run 402 stress 9.664342e-05 
    ## ... Procrustes: rmse 0.0001652424  max resid 0.000227803 
    ## ... Similar to previous best
    ## Run 403 stress 9.550745e-05 
    ## ... Procrustes: rmse 0.0001925591  max resid 0.0002869177 
    ## ... Similar to previous best
    ## Run 404 stress 8.346674e-05 
    ## ... Procrustes: rmse 0.0001771738  max resid 0.0002415221 
    ## ... Similar to previous best
    ## Run 405 stress 9.81208e-05 
    ## ... Procrustes: rmse 0.0002178702  max resid 0.0002995862 
    ## ... Similar to previous best
    ## Run 406 stress 8.946419e-05 
    ## ... Procrustes: rmse 0.0001752386  max resid 0.0002532055 
    ## ... Similar to previous best
    ## Run 407 stress 9.651529e-05 
    ## ... Procrustes: rmse 0.0002113533  max resid 0.0003072847 
    ## ... Similar to previous best
    ## Run 408 stress 9.467429e-05 
    ## ... Procrustes: rmse 0.0001930479  max resid 0.0002697817 
    ## ... Similar to previous best
    ## Run 409 stress 9.436692e-05 
    ## ... Procrustes: rmse 0.0002106301  max resid 0.0003036858 
    ## ... Similar to previous best
    ## Run 410 stress 9.71687e-05 
    ## ... Procrustes: rmse 0.000186433  max resid 0.0002644942 
    ## ... Similar to previous best
    ## Run 411 stress 9.979503e-05 
    ## ... Procrustes: rmse 0.0001261495  max resid 0.0002084886 
    ## ... Similar to previous best
    ## Run 412 stress 9.550905e-05 
    ## ... Procrustes: rmse 0.0001659487  max resid 0.0002175373 
    ## ... Similar to previous best
    ## Run 413 stress 9.373518e-05 
    ## ... Procrustes: rmse 0.0001780054  max resid 0.0002533864 
    ## ... Similar to previous best
    ## Run 414 stress 8.351828e-05 
    ## ... Procrustes: rmse 0.0001389834  max resid 0.0001741848 
    ## ... Similar to previous best
    ## Run 415 stress 9.293433e-05 
    ## ... Procrustes: rmse 0.0001873289  max resid 0.0002637603 
    ## ... Similar to previous best
    ## Run 416 stress 9.584772e-05 
    ## ... Procrustes: rmse 0.0001704581  max resid 0.0002100918 
    ## ... Similar to previous best
    ## Run 417 stress 9.79712e-05 
    ## ... Procrustes: rmse 0.0002178691  max resid 0.00029976 
    ## ... Similar to previous best
    ## Run 418 stress 8.525763e-05 
    ## ... Procrustes: rmse 8.943703e-05  max resid 0.0001763354 
    ## ... Similar to previous best
    ## Run 419 stress 9.687506e-05 
    ## ... Procrustes: rmse 0.0002025717  max resid 0.0002966672 
    ## ... Similar to previous best
    ## Run 420 stress 9.103604e-05 
    ## ... Procrustes: rmse 0.0001940636  max resid 0.0002608609 
    ## ... Similar to previous best
    ## Run 421 stress 8.800832e-05 
    ## ... Procrustes: rmse 0.0001751787  max resid 0.0002377556 
    ## ... Similar to previous best
    ## Run 422 stress 9.187409e-05 
    ## ... Procrustes: rmse 0.0001408426  max resid 0.0001903875 
    ## ... Similar to previous best
    ## Run 423 stress 9.784265e-05 
    ## ... Procrustes: rmse 0.0001955974  max resid 0.0002752975 
    ## ... Similar to previous best
    ## Run 424 stress 9.144881e-05 
    ## ... Procrustes: rmse 0.0001570604  max resid 0.0001968715 
    ## ... Similar to previous best
    ## Run 425 stress 9.841335e-05 
    ## ... Procrustes: rmse 0.0002031145  max resid 0.000310059 
    ## ... Similar to previous best
    ## Run 426 stress 9.189455e-05 
    ## ... Procrustes: rmse 0.0001291741  max resid 0.000175006 
    ## ... Similar to previous best
    ## Run 427 stress 9.945494e-05 
    ## ... Procrustes: rmse 0.0001896571  max resid 0.0002287025 
    ## ... Similar to previous best
    ## Run 428 stress 9.986043e-05 
    ## ... Procrustes: rmse 0.0002102766  max resid 0.0002899591 
    ## ... Similar to previous best
    ## Run 429 stress 7.665617e-05 
    ## ... Procrustes: rmse 0.0001452194  max resid 0.0002097712 
    ## ... Similar to previous best
    ## Run 430 stress 8.695743e-05 
    ## ... Procrustes: rmse 9.657964e-05  max resid 0.000140306 
    ## ... Similar to previous best
    ## Run 431 stress 9.986026e-05 
    ## ... Procrustes: rmse 0.0001781475  max resid 0.0002100772 
    ## ... Similar to previous best
    ## Run 432 stress 9.967241e-05 
    ## ... Procrustes: rmse 0.0002179332  max resid 0.0002913664 
    ## ... Similar to previous best
    ## Run 433 stress 9.584508e-05 
    ## ... Procrustes: rmse 0.0001891045  max resid 0.0002958313 
    ## ... Similar to previous best
    ## Run 434 stress 9.84968e-05 
    ## ... Procrustes: rmse 0.0001580994  max resid 0.0002016526 
    ## ... Similar to previous best
    ## Run 435 stress 9.972184e-05 
    ## ... Procrustes: rmse 0.00020649  max resid 0.0002857796 
    ## ... Similar to previous best
    ## Run 436 stress 9.638042e-05 
    ## ... Procrustes: rmse 0.0001549928  max resid 0.0002007668 
    ## ... Similar to previous best
    ## Run 437 stress 9.836284e-05 
    ## ... Procrustes: rmse 0.0002013356  max resid 0.0002708704 
    ## ... Similar to previous best
    ## Run 438 stress 9.545904e-05 
    ## ... Procrustes: rmse 0.0001986823  max resid 0.0002900948 
    ## ... Similar to previous best
    ## Run 439 stress 9.418382e-05 
    ## ... Procrustes: rmse 0.0001508488  max resid 0.0002189653 
    ## ... Similar to previous best
    ## Run 440 stress 9.109537e-05 
    ## ... Procrustes: rmse 0.0001715961  max resid 0.0002414139 
    ## ... Similar to previous best
    ## Run 441 stress 9.676038e-05 
    ## ... Procrustes: rmse 0.0002156795  max resid 0.0002957415 
    ## ... Similar to previous best
    ## Run 442 stress 9.891461e-05 
    ## ... Procrustes: rmse 0.0001995423  max resid 0.0002642378 
    ## ... Similar to previous best
    ## Run 443 stress 9.793403e-05 
    ## ... Procrustes: rmse 0.00020289  max resid 0.000307588 
    ## ... Similar to previous best
    ## Run 444 stress 9.85287e-05 
    ## ... Procrustes: rmse 0.0001463009  max resid 0.0002195683 
    ## ... Similar to previous best
    ## Run 445 stress 9.845008e-05 
    ## ... Procrustes: rmse 0.0001732584  max resid 0.0002780696 
    ## ... Similar to previous best
    ## Run 446 stress 9.570474e-05 
    ## ... Procrustes: rmse 0.00018275  max resid 0.0002597808 
    ## ... Similar to previous best
    ## Run 447 stress 9.84306e-05 
    ## ... Procrustes: rmse 0.000207369  max resid 0.0002783245 
    ## ... Similar to previous best
    ## Run 448 stress 9.769376e-05 
    ## ... Procrustes: rmse 0.0002006979  max resid 0.0002774497 
    ## ... Similar to previous best
    ## Run 449 stress 9.982677e-05 
    ## ... Procrustes: rmse 0.0001862621  max resid 0.0002261677 
    ## ... Similar to previous best
    ## Run 450 stress 9.758774e-05 
    ## ... Procrustes: rmse 0.0002068509  max resid 0.0002759491 
    ## ... Similar to previous best
    ## Run 451 stress 9.579028e-05 
    ## ... Procrustes: rmse 0.0001579659  max resid 0.0002447331 
    ## ... Similar to previous best
    ## Run 452 stress 9.039315e-05 
    ## ... Procrustes: rmse 0.0001697262  max resid 0.000205704 
    ## ... Similar to previous best
    ## Run 453 stress 9.116327e-05 
    ## ... Procrustes: rmse 0.0001670407  max resid 0.0002553041 
    ## ... Similar to previous best
    ## Run 454 stress 9.871515e-05 
    ## ... Procrustes: rmse 0.0001879378  max resid 0.00022169 
    ## ... Similar to previous best
    ## Run 455 stress 9.082693e-05 
    ## ... Procrustes: rmse 0.000151325  max resid 0.0002454122 
    ## ... Similar to previous best
    ## Run 456 stress 9.978329e-05 
    ## ... Procrustes: rmse 0.0002077672  max resid 0.0003031252 
    ## ... Similar to previous best
    ## Run 457 stress 9.613835e-05 
    ## ... Procrustes: rmse 0.0002133394  max resid 0.0002930763 
    ## ... Similar to previous best
    ## Run 458 stress 9.935153e-05 
    ## ... Procrustes: rmse 0.0002047112  max resid 0.0002739423 
    ## ... Similar to previous best
    ## Run 459 stress 9.865603e-05 
    ## ... Procrustes: rmse 0.0002143032  max resid 0.0002876713 
    ## ... Similar to previous best
    ## Run 460 stress 9.698042e-05 
    ## ... Procrustes: rmse 0.000202672  max resid 0.0003050582 
    ## ... Similar to previous best
    ## Run 461 stress 9.748768e-05 
    ## ... Procrustes: rmse 0.0001312927  max resid 0.0001900826 
    ## ... Similar to previous best
    ## Run 462 stress 8.779163e-05 
    ## ... Procrustes: rmse 0.0001906663  max resid 0.0002581753 
    ## ... Similar to previous best
    ## Run 463 stress 9.553541e-05 
    ## ... Procrustes: rmse 0.0001612925  max resid 0.0002492721 
    ## ... Similar to previous best
    ## Run 464 stress 9.794276e-05 
    ## ... Procrustes: rmse 0.0001736389  max resid 0.0002595075 
    ## ... Similar to previous best
    ## Run 465 stress 9.81641e-05 
    ## ... Procrustes: rmse 0.0001140996  max resid 0.0002203212 
    ## ... Similar to previous best
    ## Run 466 stress 9.317864e-05 
    ## ... Procrustes: rmse 0.0001690214  max resid 0.0002430298 
    ## ... Similar to previous best
    ## Run 467 stress 8.696761e-05 
    ## ... Procrustes: rmse 0.0001423393  max resid 0.0001888702 
    ## ... Similar to previous best
    ## Run 468 stress 9.932236e-05 
    ## ... Procrustes: rmse 0.0001886248  max resid 0.0002281542 
    ## ... Similar to previous best
    ## Run 469 stress 9.932843e-05 
    ## ... Procrustes: rmse 0.0001934698  max resid 0.0002976101 
    ## ... Similar to previous best
    ## Run 470 stress 9.894226e-05 
    ## ... Procrustes: rmse 0.0001739832  max resid 0.0002270456 
    ## ... Similar to previous best
    ## Run 471 stress 9.201806e-05 
    ## ... Procrustes: rmse 0.0001552213  max resid 0.0002434323 
    ## ... Similar to previous best
    ## Run 472 stress 9.889234e-05 
    ## ... Procrustes: rmse 0.0002134242  max resid 0.0002873999 
    ## ... Similar to previous best
    ## Run 473 stress 9.422506e-05 
    ## ... Procrustes: rmse 0.0001571746  max resid 0.0002075629 
    ## ... Similar to previous best
    ## Run 474 stress 9.650269e-05 
    ## ... Procrustes: rmse 0.0001983978  max resid 0.0003041114 
    ## ... Similar to previous best
    ## Run 475 stress 8.41412e-05 
    ## ... Procrustes: rmse 0.0001573198  max resid 0.0002378924 
    ## ... Similar to previous best
    ## Run 476 stress 8.873472e-05 
    ## ... Procrustes: rmse 0.0001679669  max resid 0.0002388335 
    ## ... Similar to previous best
    ## Run 477 stress 9.664687e-05 
    ## ... Procrustes: rmse 0.0001655175  max resid 0.0002134415 
    ## ... Similar to previous best
    ## Run 478 stress 8.807478e-05 
    ## ... Procrustes: rmse 0.0001488597  max resid 0.000197556 
    ## ... Similar to previous best
    ## Run 479 stress 9.50717e-05 
    ## ... Procrustes: rmse 0.0001284939  max resid 0.000188659 
    ## ... Similar to previous best
    ## Run 480 stress 9.396276e-05 
    ## ... Procrustes: rmse 0.0001929996  max resid 0.0002963818 
    ## ... Similar to previous best
    ## Run 481 stress 9.57736e-05 
    ## ... Procrustes: rmse 0.0001598721  max resid 0.0002120317 
    ## ... Similar to previous best
    ## Run 482 stress 9.854277e-05 
    ## ... Procrustes: rmse 0.0002105416  max resid 0.0002921279 
    ## ... Similar to previous best
    ## Run 483 stress 9.413957e-05 
    ## ... Procrustes: rmse 0.0001484312  max resid 0.0001794799 
    ## ... Similar to previous best
    ## Run 484 stress 9.333466e-05 
    ## ... Procrustes: rmse 0.0001283262  max resid 0.0001744776 
    ## ... Similar to previous best
    ## Run 485 stress 9.588942e-05 
    ## ... Procrustes: rmse 0.0002141325  max resid 0.0003012394 
    ## ... Similar to previous best
    ## Run 486 stress 9.257502e-05 
    ## ... Procrustes: rmse 0.0001682994  max resid 0.0002509386 
    ## ... Similar to previous best
    ## Run 487 stress 9.197129e-05 
    ## ... Procrustes: rmse 0.0002010278  max resid 0.0002750164 
    ## ... Similar to previous best
    ## Run 488 stress 9.485568e-05 
    ## ... Procrustes: rmse 0.000197191  max resid 0.0002880937 
    ## ... Similar to previous best
    ## Run 489 stress 9.223261e-05 
    ## ... Procrustes: rmse 0.0001871926  max resid 0.0002682503 
    ## ... Similar to previous best
    ## Run 490 stress 9.269364e-05 
    ## ... Procrustes: rmse 0.0001603381  max resid 0.0002044764 
    ## ... Similar to previous best
    ## Run 491 stress 9.924459e-05 
    ## ... Procrustes: rmse 0.0002058266  max resid 0.0002821152 
    ## ... Similar to previous best
    ## Run 492 stress 8.751201e-05 
    ## ... Procrustes: rmse 0.0001727425  max resid 0.0002645692 
    ## ... Similar to previous best
    ## Run 493 stress 9.658559e-05 
    ## ... Procrustes: rmse 0.0002015591  max resid 0.0002952133 
    ## ... Similar to previous best
    ## Run 494 stress 9.536581e-05 
    ## ... Procrustes: rmse 0.0001452223  max resid 0.0002506188 
    ## ... Similar to previous best
    ## Run 495 stress 9.811144e-05 
    ## ... Procrustes: rmse 0.0001717412  max resid 0.0002271169 
    ## ... Similar to previous best
    ## Run 496 stress 9.699108e-05 
    ## ... Procrustes: rmse 0.0001831356  max resid 0.0002257187 
    ## ... Similar to previous best
    ## Run 497 stress 9.128443e-05 
    ## ... Procrustes: rmse 0.0001726836  max resid 0.0002485372 
    ## ... Similar to previous best
    ## Run 498 stress 9.110785e-05 
    ## ... Procrustes: rmse 0.0001700242  max resid 0.0002060581 
    ## ... Similar to previous best
    ## Run 499 stress 9.926074e-05 
    ## ... Procrustes: rmse 0.0001237194  max resid 0.0001650962 
    ## ... Similar to previous best
    ## Run 500 stress 9.839859e-05 
    ## ... Procrustes: rmse 0.0002143501  max resid 0.0003172108 
    ## ... Similar to previous best
    ## *** Best solution repeated 272 times

    ## Warning in metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Mixed lakes
SD_beta_env_M_NMDS <- metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.001279543 
    ## Run 1 stress 0.001386719 
    ## ... Procrustes: rmse 0.001786825  max resid 0.002478425 
    ## ... Similar to previous best
    ## Run 2 stress 0.00130263 
    ## ... Procrustes: rmse 0.0003907482  max resid 0.0005424721 
    ## ... Similar to previous best
    ## Run 3 stress 8.896918e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04311338  max resid 0.05939338 
    ## Run 4 stress 9.061554e-05 
    ## ... Procrustes: rmse 0.0001900176  max resid 0.0002866281 
    ## ... Similar to previous best
    ## Run 5 stress 8.668926e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001491549  max resid 0.0003542605 
    ## ... Similar to previous best
    ## Run 6 stress 0.001403799 
    ## Run 7 stress 0.1990774 
    ## Run 8 stress 6.390123e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0007354831  max resid 0.001095268 
    ## ... Similar to previous best
    ## Run 9 stress 5.284261e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0007442867  max resid 0.001022954 
    ## ... Similar to previous best
    ## Run 10 stress 0.001369813 
    ## Run 11 stress 0.0001014675 
    ## ... Procrustes: rmse 0.007356313  max resid 0.01008419 
    ## Run 12 stress 8.679138e-05 
    ## ... Procrustes: rmse 8.464228e-05  max resid 0.0001272144 
    ## ... Similar to previous best
    ## Run 13 stress 8.844095e-05 
    ## ... Procrustes: rmse 0.0001180794  max resid 0.0001822175 
    ## ... Similar to previous best
    ## Run 14 stress 0.1491957 
    ## Run 15 stress 0.0003882198 
    ## ... Procrustes: rmse 0.01441669  max resid 0.0198268 
    ## Run 16 stress 9.543293e-05 
    ## ... Procrustes: rmse 0.0001744983  max resid 0.0002447111 
    ## ... Similar to previous best
    ## Run 17 stress 8.810424e-05 
    ## ... Procrustes: rmse 0.0001839579  max resid 0.0002916874 
    ## ... Similar to previous best
    ## Run 18 stress 9.562651e-05 
    ## ... Procrustes: rmse 0.0001442137  max resid 0.0002734732 
    ## ... Similar to previous best
    ## Run 19 stress 9.58221e-05 
    ## ... Procrustes: rmse 0.0001821095  max resid 0.000301296 
    ## ... Similar to previous best
    ## Run 20 stress 9.364018e-05 
    ## ... Procrustes: rmse 0.000174039  max resid 0.000240997 
    ## ... Similar to previous best
    ## Run 21 stress 9.215109e-05 
    ## ... Procrustes: rmse 0.000168851  max resid 0.0002627465 
    ## ... Similar to previous best
    ## Run 22 stress 7.224079e-05 
    ## ... Procrustes: rmse 0.000520504  max resid 0.0007117232 
    ## ... Similar to previous best
    ## Run 23 stress 0.1491957 
    ## Run 24 stress 0.001303238 
    ## Run 25 stress 0.1491957 
    ## Run 26 stress 0.00117672 
    ## Run 27 stress 0.001377832 
    ## Run 28 stress 9.100054e-05 
    ## ... Procrustes: rmse 0.0001737098  max resid 0.0002899501 
    ## ... Similar to previous best
    ## Run 29 stress 0.001206901 
    ## Run 30 stress 0.1990774 
    ## Run 31 stress 0.00107532 
    ## Run 32 stress 0.0004619603 
    ## ... Procrustes: rmse 0.02473885  max resid 0.03408405 
    ## Run 33 stress 0.00102192 
    ## Run 34 stress 0.0014492 
    ## Run 35 stress 8.191206e-05 
    ## ... Procrustes: rmse 0.001976473  max resid 0.002680774 
    ## ... Similar to previous best
    ## Run 36 stress 0.2848512 
    ## Run 37 stress 8.199062e-05 
    ## ... Procrustes: rmse 0.0001043806  max resid 0.0001531507 
    ## ... Similar to previous best
    ## Run 38 stress 0.1491957 
    ## Run 39 stress 9.89629e-05 
    ## ... Procrustes: rmse 0.0001246465  max resid 0.0002027342 
    ## ... Similar to previous best
    ## Run 40 stress 0.0004928961 
    ## ... Procrustes: rmse 0.01625572  max resid 0.02236422 
    ## Run 41 stress 0.1990774 
    ## Run 42 stress 0.001343382 
    ## Run 43 stress 9.687088e-05 
    ## ... Procrustes: rmse 0.0001303346  max resid 0.0001939944 
    ## ... Similar to previous best
    ## Run 44 stress 0.2797322 
    ## Run 45 stress 9.53956e-05 
    ## ... Procrustes: rmse 0.0001246081  max resid 0.0002073314 
    ## ... Similar to previous best
    ## Run 46 stress 0.001204406 
    ## Run 47 stress 8.081589e-05 
    ## ... Procrustes: rmse 0.0004055034  max resid 0.0005183053 
    ## ... Similar to previous best
    ## Run 48 stress 9.46871e-05 
    ## ... Procrustes: rmse 0.0001415629  max resid 0.0002596783 
    ## ... Similar to previous best
    ## Run 49 stress 0.001319374 
    ## Run 50 stress 0.00133837 
    ## Run 51 stress 0.1990774 
    ## Run 52 stress 0.001270678 
    ## Run 53 stress 0.1776012 
    ## Run 54 stress 9.435665e-05 
    ## ... Procrustes: rmse 0.0001252759  max resid 0.0001979716 
    ## ... Similar to previous best
    ## Run 55 stress 0.001663965 
    ## Run 56 stress 0.1491957 
    ## Run 57 stress 0.001334446 
    ## Run 58 stress 9.762397e-05 
    ## ... Procrustes: rmse 0.0001836929  max resid 0.0002971792 
    ## ... Similar to previous best
    ## Run 59 stress 0.001365706 
    ## Run 60 stress 0.001279115 
    ## Run 61 stress 0.1491957 
    ## Run 62 stress 0.1990774 
    ## Run 63 stress 0.001358122 
    ## Run 64 stress 0.001056818 
    ## Run 65 stress 0.000449991 
    ## ... Procrustes: rmse 0.01552985  max resid 0.02136261 
    ## Run 66 stress 0.001347429 
    ## Run 67 stress 0.001265256 
    ## Run 68 stress 0.001271739 
    ## Run 69 stress 9.884329e-05 
    ## ... Procrustes: rmse 0.0001609423  max resid 0.0002233679 
    ## ... Similar to previous best
    ## Run 70 stress 9.033111e-05 
    ## ... Procrustes: rmse 0.0001253044  max resid 0.0001840458 
    ## ... Similar to previous best
    ## Run 71 stress 0.0004894088 
    ## ... Procrustes: rmse 0.01616914  max resid 0.02224531 
    ## Run 72 stress 0.177601 
    ## Run 73 stress 0.1776021 
    ## Run 74 stress 9.377284e-05 
    ## ... Procrustes: rmse 0.00012521  max resid 0.0001927174 
    ## ... Similar to previous best
    ## Run 75 stress 9.543247e-05 
    ## ... Procrustes: rmse 0.0001325214  max resid 0.0002034095 
    ## ... Similar to previous best
    ## Run 76 stress 0.1491957 
    ## Run 77 stress 8.878628e-05 
    ## ... Procrustes: rmse 0.0001256759  max resid 0.0001849767 
    ## ... Similar to previous best
    ## Run 78 stress 0.001287837 
    ## Run 79 stress 0.0004769346 
    ## ... Procrustes: rmse 0.01598754  max resid 0.02199399 
    ## Run 80 stress 0.00141568 
    ## Run 81 stress 0.0004895709 
    ## ... Procrustes: rmse 0.01619568  max resid 0.0222806 
    ## Run 82 stress 0.2848514 
    ## Run 83 stress 0.001199436 
    ## Run 84 stress 9.026434e-05 
    ## ... Procrustes: rmse 0.000127277  max resid 0.0001913886 
    ## ... Similar to previous best
    ## Run 85 stress 9.35111e-05 
    ## ... Procrustes: rmse 0.0001605645  max resid 0.0002107964 
    ## ... Similar to previous best
    ## Run 86 stress 0.0004664731 
    ## ... Procrustes: rmse 0.01566941  max resid 0.021561 
    ## Run 87 stress 0.0004478224 
    ## ... Procrustes: rmse 0.01549325  max resid 0.02131227 
    ## Run 88 stress 0.1776012 
    ## Run 89 stress 9.367177e-05 
    ## ... Procrustes: rmse 0.0001307937  max resid 0.0001930418 
    ## ... Similar to previous best
    ## Run 90 stress 7.739765e-05 
    ## ... Procrustes: rmse 0.0001000991  max resid 0.0001353284 
    ## ... Similar to previous best
    ## Run 91 stress 0.001130517 
    ## Run 92 stress 0.2568282 
    ## Run 93 stress 0.0002375231 
    ## ... Procrustes: rmse 0.01861294  max resid 0.02561745 
    ## Run 94 stress 9.88446e-05 
    ## ... Procrustes: rmse 0.0001260031  max resid 0.0001835632 
    ## ... Similar to previous best
    ## Run 95 stress 0.001297113 
    ## Run 96 stress 0.3083098 
    ## Run 97 stress 4.557763e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0006427156  max resid 0.0008910304 
    ## ... Similar to previous best
    ## Run 98 stress 0.001442272 
    ## Run 99 stress 8.141475e-05 
    ## ... Procrustes: rmse 0.0005070584  max resid 0.0007016265 
    ## ... Similar to previous best
    ## Run 100 stress 9.668125e-05 
    ## ... Procrustes: rmse 0.0005862639  max resid 0.0008300299 
    ## ... Similar to previous best
    ## Run 101 stress 0.0004387551 
    ## ... Procrustes: rmse 0.01521247  max resid 0.02215887 
    ## Run 102 stress 9.520073e-05 
    ## ... Procrustes: rmse 0.0005297499  max resid 0.0007842899 
    ## ... Similar to previous best
    ## Run 103 stress 0.00102999 
    ## Run 104 stress 5.842764e-05 
    ## ... Procrustes: rmse 0.0006523607  max resid 0.0008569265 
    ## ... Similar to previous best
    ## Run 105 stress 0.0002553482 
    ## ... Procrustes: rmse 0.01157473  max resid 0.01713457 
    ## Run 106 stress 9.345019e-05 
    ## ... Procrustes: rmse 0.0005789663  max resid 0.0007875356 
    ## ... Similar to previous best
    ## Run 107 stress 9.730401e-05 
    ## ... Procrustes: rmse 0.006288547  max resid 0.009834102 
    ## Run 108 stress 7.108135e-05 
    ## ... Procrustes: rmse 0.0003941615  max resid 0.0005178931 
    ## ... Similar to previous best
    ## Run 109 stress 0.0004839183 
    ## ... Procrustes: rmse 0.01596048  max resid 0.02319242 
    ## Run 110 stress 9.730037e-05 
    ## ... Procrustes: rmse 0.0006777186  max resid 0.001074892 
    ## ... Similar to previous best
    ## Run 111 stress 0.001388148 
    ## Run 112 stress 9.991264e-05 
    ## ... Procrustes: rmse 0.0005194136  max resid 0.0007934004 
    ## ... Similar to previous best
    ## Run 113 stress 9.659632e-05 
    ## ... Procrustes: rmse 0.0005723353  max resid 0.0008899416 
    ## ... Similar to previous best
    ## Run 114 stress 0.001386671 
    ## Run 115 stress 9.482398e-05 
    ## ... Procrustes: rmse 0.0006667068  max resid 0.0009733497 
    ## ... Similar to previous best
    ## Run 116 stress 0.1491957 
    ## Run 117 stress 9.045647e-05 
    ## ... Procrustes: rmse 0.002834719  max resid 0.005014825 
    ## ... Similar to previous best
    ## Run 118 stress 0.0004669253 
    ## ... Procrustes: rmse 0.0157001  max resid 0.02283224 
    ## Run 119 stress 0.0001115471 
    ## ... Procrustes: rmse 0.007607641  max resid 0.01164742 
    ## Run 120 stress 0.1990774 
    ## Run 121 stress 0.0009759398 
    ## Run 122 stress 0.0008530849 
    ## Run 123 stress 0.1990774 
    ## Run 124 stress 0.0004775509 
    ## ... Procrustes: rmse 0.0158751  max resid 0.02307372 
    ## Run 125 stress 0.0004421719 
    ## ... Procrustes: rmse 0.02439  max resid 0.03366667 
    ## Run 126 stress 0.2646973 
    ## Run 127 stress 7.130512e-05 
    ## ... Procrustes: rmse 5.336094e-05  max resid 9.286648e-05 
    ## ... Similar to previous best
    ## Run 128 stress 0.001325684 
    ## Run 129 stress 0.1990774 
    ## Run 130 stress 0.3082926 
    ## Run 131 stress 0.2806642 
    ## Run 132 stress 0.001356704 
    ## Run 133 stress 0.001136553 
    ## Run 134 stress 0.001377292 
    ## Run 135 stress 0.001327145 
    ## Run 136 stress 0.2373687 
    ## Run 137 stress 0.0004841915 
    ## ... Procrustes: rmse 0.0159249  max resid 0.02314336 
    ## Run 138 stress 8.035339e-05 
    ## ... Procrustes: rmse 0.002873088  max resid 0.003893076 
    ## ... Similar to previous best
    ## Run 139 stress 0.2848523 
    ## Run 140 stress 9.7835e-05 
    ## ... Procrustes: rmse 0.0006057435  max resid 0.001006302 
    ## ... Similar to previous best
    ## Run 141 stress 0.0004609097 
    ## ... Procrustes: rmse 0.01559403  max resid 0.02268591 
    ## Run 142 stress 0.001192863 
    ## Run 143 stress 0.001193467 
    ## Run 144 stress 0.1990774 
    ## Run 145 stress 0.001333772 
    ## Run 146 stress 0.0001361851 
    ## ... Procrustes: rmse 0.008411957  max resid 0.01275978 
    ## Run 147 stress 0.2696744 
    ## Run 148 stress 0.001358236 
    ## Run 149 stress 0.0004565425 
    ## ... Procrustes: rmse 0.01552306  max resid 0.02258781 
    ## Run 150 stress 0.2612473 
    ## Run 151 stress 0.0004889761 
    ## ... Procrustes: rmse 0.01606976  max resid 0.02334344 
    ## Run 152 stress 0.1776015 
    ## Run 153 stress 0.00124175 
    ## Run 154 stress 9.84578e-05 
    ## ... Procrustes: rmse 0.0005725314  max resid 0.0007855554 
    ## ... Similar to previous best
    ## Run 155 stress 9.557882e-05 
    ## ... Procrustes: rmse 0.0005977707  max resid 0.001015658 
    ## ... Similar to previous best
    ## Run 156 stress 0.1491957 
    ## Run 157 stress 0.0004450418 
    ## ... Procrustes: rmse 0.02365007  max resid 0.03263087 
    ## Run 158 stress 9.972449e-05 
    ## ... Procrustes: rmse 0.0006553476  max resid 0.001078471 
    ## ... Similar to previous best
    ## Run 159 stress 9.527679e-05 
    ## ... Procrustes: rmse 0.0006281744  max resid 0.0008572136 
    ## ... Similar to previous best
    ## Run 160 stress 9.071606e-05 
    ## ... Procrustes: rmse 0.0006097233  max resid 0.0009950638 
    ## ... Similar to previous best
    ## Run 161 stress 0.001164779 
    ## Run 162 stress 9.677344e-05 
    ## ... Procrustes: rmse 0.0005236988  max resid 0.0007825516 
    ## ... Similar to previous best
    ## Run 163 stress 0.001325607 
    ## Run 164 stress 0.0002854434 
    ## ... Procrustes: rmse 0.0197692  max resid 0.02726987 
    ## Run 165 stress 9.849753e-05 
    ## ... Procrustes: rmse 0.0006479269  max resid 0.001016417 
    ## ... Similar to previous best
    ## Run 166 stress 8.922607e-05 
    ## ... Procrustes: rmse 0.000631746  max resid 0.0009723158 
    ## ... Similar to previous best
    ## Run 167 stress 0.2373687 
    ## Run 168 stress 0.3082926 
    ## Run 169 stress 9.588017e-05 
    ## ... Procrustes: rmse 0.0006031978  max resid 0.0009932729 
    ## ... Similar to previous best
    ## Run 170 stress 0.1990774 
    ## Run 171 stress 0.2848515 
    ## Run 172 stress 0.1990774 
    ## Run 173 stress 9.418302e-05 
    ## ... Procrustes: rmse 0.0006047873  max resid 0.001027873 
    ## ... Similar to previous best
    ## Run 174 stress 0.0001037599 
    ## ... Procrustes: rmse 0.007332101  max resid 0.011266 
    ## Run 175 stress 8.899905e-05 
    ## ... Procrustes: rmse 0.0005643972  max resid 0.0008874115 
    ## ... Similar to previous best
    ## Run 176 stress 0.1990774 
    ## Run 177 stress 9.58945e-05 
    ## ... Procrustes: rmse 0.0006324168  max resid 0.0009769064 
    ## ... Similar to previous best
    ## Run 178 stress 9.200521e-05 
    ## ... Procrustes: rmse 0.00060932  max resid 0.001033707 
    ## ... Similar to previous best
    ## Run 179 stress 0.001427798 
    ## Run 180 stress 0.1491957 
    ## Run 181 stress 0.00149119 
    ## Run 182 stress 0.00120271 
    ## Run 183 stress 9.449405e-05 
    ## ... Procrustes: rmse 0.0005964085  max resid 0.001012943 
    ## ... Similar to previous best
    ## Run 184 stress 0.001305238 
    ## Run 185 stress 0.001381768 
    ## Run 186 stress 0.1491957 
    ## Run 187 stress 0.001384475 
    ## Run 188 stress 0.1491957 
    ## Run 189 stress 0.0004833752 
    ## ... Procrustes: rmse 0.01597601  max resid 0.02321276 
    ## Run 190 stress 9.645077e-05 
    ## ... Procrustes: rmse 0.000522921  max resid 0.0007882898 
    ## ... Similar to previous best
    ## Run 191 stress 8.628196e-05 
    ## ... Procrustes: rmse 0.0005441733  max resid 0.0008475994 
    ## ... Similar to previous best
    ## Run 192 stress 8.814255e-05 
    ## ... Procrustes: rmse 0.0006217333  max resid 0.0008526297 
    ## ... Similar to previous best
    ## Run 193 stress 0.0003986094 
    ## ... Procrustes: rmse 0.01449545  max resid 0.02116864 
    ## Run 194 stress 0.001394714 
    ## Run 195 stress 0.00140726 
    ## Run 196 stress 0.0001111847 
    ## ... Procrustes: rmse 0.01208922  max resid 0.01665454 
    ## Run 197 stress 0.2696744 
    ## Run 198 stress 8.732192e-05 
    ## ... Procrustes: rmse 0.0006273618  max resid 0.0008527914 
    ## ... Similar to previous best
    ## Run 199 stress 0.1776018 
    ## Run 200 stress 0.0004598042 
    ## ... Procrustes: rmse 0.01545491  max resid 0.02248087 
    ## Run 201 stress 9.745846e-05 
    ## ... Procrustes: rmse 0.006271238  max resid 0.009809051 
    ## Run 202 stress 0.001369611 
    ## Run 203 stress 0.001379033 
    ## Run 204 stress 0.001358533 
    ## Run 205 stress 0.001295469 
    ## Run 206 stress 0.0004837311 
    ## ... Procrustes: rmse 0.01597861  max resid 0.02321664 
    ## Run 207 stress 0.0005021128 
    ## ... Procrustes: rmse 0.01628636  max resid 0.02364157 
    ## Run 208 stress 0.1776021 
    ## Run 209 stress 9.891348e-05 
    ## ... Procrustes: rmse 0.0006096043  max resid 0.001016889 
    ## ... Similar to previous best
    ## Run 210 stress 0.1990774 
    ## Run 211 stress 9.978417e-05 
    ## ... Procrustes: rmse 0.009491299  max resid 0.01314342 
    ## Run 212 stress 0.001169676 
    ## Run 213 stress 0.001360936 
    ## Run 214 stress 9.061728e-05 
    ## ... Procrustes: rmse 0.0006324693  max resid 0.0009751964 
    ## ... Similar to previous best
    ## Run 215 stress 0.0001163342 
    ## ... Procrustes: rmse 0.007771465  max resid 0.01187443 
    ## Run 216 stress 9.193962e-05 
    ## ... Procrustes: rmse 0.0005715997  max resid 0.0008907449 
    ## ... Similar to previous best
    ## Run 217 stress 9.462064e-05 
    ## ... Procrustes: rmse 0.000596959  max resid 0.001012554 
    ## ... Similar to previous best
    ## Run 218 stress 0.001327536 
    ## Run 219 stress 9.433794e-05 
    ## ... Procrustes: rmse 0.0006124865  max resid 0.001003754 
    ## ... Similar to previous best
    ## Run 220 stress 9.732709e-05 
    ## ... Procrustes: rmse 0.0006046955  max resid 0.001036563 
    ## ... Similar to previous best
    ## Run 221 stress 7.72404e-05 
    ## ... Procrustes: rmse 0.0002060046  max resid 0.0003153257 
    ## ... Similar to previous best
    ## Run 222 stress 8.760509e-05 
    ## ... Procrustes: rmse 0.0007322038  max resid 0.001614175 
    ## ... Similar to previous best
    ## Run 223 stress 0.001160924 
    ## Run 224 stress 0.1990774 
    ## Run 225 stress 0.0004754299 
    ## ... Procrustes: rmse 0.01583899  max resid 0.02302309 
    ## Run 226 stress 0.2520602 
    ## Run 227 stress 0.001289183 
    ## Run 228 stress 0.001372672 
    ## Run 229 stress 0.0004829342 
    ## ... Procrustes: rmse 0.01596926  max resid 0.02320271 
    ## Run 230 stress 0.001256771 
    ## Run 231 stress 0.00125186 
    ## Run 232 stress 9.452461e-05 
    ## ... Procrustes: rmse 0.0005937291  max resid 0.001006115 
    ## ... Similar to previous best
    ## Run 233 stress 0.0005028133 
    ## ... Procrustes: rmse 0.01629591  max resid 0.0236558 
    ## Run 234 stress 0.001163117 
    ## Run 235 stress 9.975741e-05 
    ## ... Procrustes: rmse 0.0006331587  max resid 0.000952259 
    ## ... Similar to previous best
    ## Run 236 stress 9.151198e-05 
    ## ... Procrustes: rmse 0.0005627018  max resid 0.0008915845 
    ## ... Similar to previous best
    ## Run 237 stress 0.0004857155 
    ## ... Procrustes: rmse 0.01600193  max resid 0.02325109 
    ## Run 238 stress 0.0004454863 
    ## ... Procrustes: rmse 0.01533103  max resid 0.02232327 
    ## Run 239 stress 9.395683e-05 
    ## ... Procrustes: rmse 0.0005251434  max resid 0.0007916084 
    ## ... Similar to previous best
    ## Run 240 stress 0.001268107 
    ## Run 241 stress 0.0004866228 
    ## ... Procrustes: rmse 0.0155665  max resid 0.02264782 
    ## Run 242 stress 9.340937e-05 
    ## ... Procrustes: rmse 0.0006375073  max resid 0.0009934779 
    ## ... Similar to previous best
    ## Run 243 stress 9.65093e-05 
    ## ... Procrustes: rmse 0.0005900609  max resid 0.0009953088 
    ## ... Similar to previous best
    ## Run 244 stress 9.692748e-05 
    ## ... Procrustes: rmse 0.0006546562  max resid 0.001071912 
    ## ... Similar to previous best
    ## Run 245 stress 0.0001161441 
    ## ... Procrustes: rmse 0.00776567  max resid 0.01186621 
    ## Run 246 stress 0.0008945018 
    ## Run 247 stress 9.298914e-05 
    ## ... Procrustes: rmse 0.0005730134  max resid 0.0008824322 
    ## ... Similar to previous best
    ## Run 248 stress 0.0004897546 
    ## ... Procrustes: rmse 0.01608135  max resid 0.02336273 
    ## Run 249 stress 0.001154277 
    ## Run 250 stress 0.001163963 
    ## Run 251 stress 0.0005005666 
    ## ... Procrustes: rmse 0.01621689  max resid 0.02354888 
    ## Run 252 stress 8.667124e-05 
    ## ... Procrustes: rmse 0.0006202161  max resid 0.001102195 
    ## ... Similar to previous best
    ## Run 253 stress 0.2440272 
    ## Run 254 stress 8.651893e-05 
    ## ... Procrustes: rmse 0.0006516186  max resid 0.001039186 
    ## ... Similar to previous best
    ## Run 255 stress 0.3083098 
    ## Run 256 stress 0.00138135 
    ## Run 257 stress 0.001149217 
    ## Run 258 stress 0.001225986 
    ## Run 259 stress 0.001191279 
    ## Run 260 stress 9.32282e-05 
    ## ... Procrustes: rmse 0.000521656  max resid 0.0007658092 
    ## ... Similar to previous best
    ## Run 261 stress 0.1990774 
    ## Run 262 stress 0.0007389142 
    ## Run 263 stress 0.001320889 
    ## Run 264 stress 0.1776014 
    ## Run 265 stress 0.177601 
    ## Run 266 stress 0.001188027 
    ## Run 267 stress 0.001405947 
    ## Run 268 stress 9.929467e-05 
    ## ... Procrustes: rmse 0.0005848497  max resid 0.000951787 
    ## ... Similar to previous best
    ## Run 269 stress 0.001313251 
    ## Run 270 stress 0.001114746 
    ## Run 271 stress 0.00125951 
    ## Run 272 stress 0.0004823428 
    ## ... Procrustes: rmse 0.01595938  max resid 0.02319011 
    ## Run 273 stress 0.001331219 
    ## Run 274 stress 0.001365031 
    ## Run 275 stress 9.204408e-05 
    ## ... Procrustes: rmse 0.0006489217  max resid 0.001035085 
    ## ... Similar to previous best
    ## Run 276 stress 0.1491957 
    ## Run 277 stress 0.1990774 
    ## Run 278 stress 0.0004682579 
    ## ... Procrustes: rmse 0.01570376  max resid 0.02284307 
    ## Run 279 stress 0.0004742229 
    ## ... Procrustes: rmse 0.01581699  max resid 0.02299351 
    ## Run 280 stress 0.00131104 
    ## Run 281 stress 0.001149804 
    ## Run 282 stress 0.001163151 
    ## Run 283 stress 0.0004856281 
    ## ... Procrustes: rmse 0.0159434  max resid 0.02317245 
    ## Run 284 stress 0.001127786 
    ## Run 285 stress 9.323932e-05 
    ## ... Procrustes: rmse 0.0005263331  max resid 0.0007901883 
    ## ... Similar to previous best
    ## Run 286 stress 9.408299e-05 
    ## ... Procrustes: rmse 0.0006456827  max resid 0.00103653 
    ## ... Similar to previous best
    ## Run 287 stress 0.0002301698 
    ## ... Procrustes: rmse 0.01098239  max resid 0.01631576 
    ## Run 288 stress 9.832023e-05 
    ## ... Procrustes: rmse 0.002017506  max resid 0.003841039 
    ## ... Similar to previous best
    ## Run 289 stress 0.1491957 
    ## Run 290 stress 0.2842805 
    ## Run 291 stress 9.719577e-05 
    ## ... Procrustes: rmse 0.0006312779  max resid 0.0009771636 
    ## ... Similar to previous best
    ## Run 292 stress 9.313137e-05 
    ## ... Procrustes: rmse 0.0002298788  max resid 0.0003695282 
    ## ... Similar to previous best
    ## Run 293 stress 8.379674e-05 
    ## ... Procrustes: rmse 0.0006579144  max resid 0.0009580606 
    ## ... Similar to previous best
    ## Run 294 stress 9.171905e-05 
    ## ... Procrustes: rmse 0.0005632902  max resid 0.0008871587 
    ## ... Similar to previous best
    ## Run 295 stress 9.174157e-05 
    ## ... Procrustes: rmse 0.0005676808  max resid 0.0008861635 
    ## ... Similar to previous best
    ## Run 296 stress 9.250421e-05 
    ## ... Procrustes: rmse 0.0008309432  max resid 0.001149651 
    ## ... Similar to previous best
    ## Run 297 stress 0.0004336684 
    ## ... Procrustes: rmse 0.02452087  max resid 0.03384647 
    ## Run 298 stress 0.1491957 
    ## Run 299 stress 0.00123958 
    ## Run 300 stress 9.658628e-05 
    ## ... Procrustes: rmse 0.0005695864  max resid 0.0008896611 
    ## ... Similar to previous best
    ## Run 301 stress 0.001273588 
    ## Run 302 stress 9.744515e-05 
    ## ... Procrustes: rmse 0.0006073663  max resid 0.001008264 
    ## ... Similar to previous best
    ## Run 303 stress 0.1491957 
    ## Run 304 stress 0.001314297 
    ## Run 305 stress 8.858843e-05 
    ## ... Procrustes: rmse 0.0005984036  max resid 0.0009984737 
    ## ... Similar to previous best
    ## Run 306 stress 0.001208269 
    ## Run 307 stress 0.001336713 
    ## Run 308 stress 0.0005200453 
    ## ... Procrustes: rmse 0.01657039  max resid 0.02404524 
    ## Run 309 stress 8.576812e-05 
    ## ... Procrustes: rmse 0.000535571  max resid 0.0007904822 
    ## ... Similar to previous best
    ## Run 310 stress 9.457613e-05 
    ## ... Procrustes: rmse 0.0006413257  max resid 0.0009924894 
    ## ... Similar to previous best
    ## Run 311 stress 0.001349257 
    ## Run 312 stress 0.00158546 
    ## Run 313 stress 9.869321e-05 
    ## ... Procrustes: rmse 0.0005566246  max resid 0.0008934411 
    ## ... Similar to previous best
    ## Run 314 stress 0.001420889 
    ## Run 315 stress 8.886237e-05 
    ## ... Procrustes: rmse 0.0005837458  max resid 0.0007836415 
    ## ... Similar to previous best
    ## Run 316 stress 0.001262142 
    ## Run 317 stress 9.750197e-05 
    ## ... Procrustes: rmse 0.0006500855  max resid 0.0009329194 
    ## ... Similar to previous best
    ## Run 318 stress 0.1776019 
    ## Run 319 stress 0.0004418579 
    ## ... Procrustes: rmse 0.01526911  max resid 0.02223726 
    ## Run 320 stress 9.004787e-05 
    ## ... Procrustes: rmse 0.0005304877  max resid 0.0007903438 
    ## ... Similar to previous best
    ## Run 321 stress 0.001280916 
    ## Run 322 stress 0.001188644 
    ## Run 323 stress 0.001186967 
    ## Run 324 stress 0.0004783574 
    ## ... Procrustes: rmse 0.01589195  max resid 0.02309719 
    ## Run 325 stress 0.0005055229 
    ## ... Procrustes: rmse 0.0265272  max resid 0.03662579 
    ## Run 326 stress 0.001272573 
    ## Run 327 stress 5.653425e-05 
    ## ... Procrustes: rmse 0.0006509644  max resid 0.0008909108 
    ## ... Similar to previous best
    ## Run 328 stress 8.246879e-05 
    ## ... Procrustes: rmse 0.0006151285  max resid 0.0009879998 
    ## ... Similar to previous best
    ## Run 329 stress 9.28287e-05 
    ## ... Procrustes: rmse 0.0005986291  max resid 0.0009859612 
    ## ... Similar to previous best
    ## Run 330 stress 0.001340164 
    ## Run 331 stress 0.001384436 
    ## Run 332 stress 9.369239e-05 
    ## ... Procrustes: rmse 0.006152385  max resid 0.009643535 
    ## Run 333 stress 9.207215e-05 
    ## ... Procrustes: rmse 0.0006120642  max resid 0.001009321 
    ## ... Similar to previous best
    ## Run 334 stress 0.001247881 
    ## Run 335 stress 9.499389e-05 
    ## ... Procrustes: rmse 0.000596086  max resid 0.001012034 
    ## ... Similar to previous best
    ## Run 336 stress 0.1776022 
    ## Run 337 stress 0.001254704 
    ## Run 338 stress 9.838247e-05 
    ## ... Procrustes: rmse 0.0006473619  max resid 0.001047792 
    ## ... Similar to previous best
    ## Run 339 stress 9.320825e-05 
    ## ... Procrustes: rmse 0.0005965579  max resid 0.001001933 
    ## ... Similar to previous best
    ## Run 340 stress 8.901829e-05 
    ## ... Procrustes: rmse 0.002208855  max resid 0.003052218 
    ## ... Similar to previous best
    ## Run 341 stress 5.257644e-05 
    ## ... Procrustes: rmse 0.0005876088  max resid 0.0007863814 
    ## ... Similar to previous best
    ## Run 342 stress 9.304845e-05 
    ## ... Procrustes: rmse 0.0005349674  max resid 0.0007827934 
    ## ... Similar to previous best
    ## Run 343 stress 9.539588e-05 
    ## ... Procrustes: rmse 0.0006515618  max resid 0.001052869 
    ## ... Similar to previous best
    ## Run 344 stress 0.001298077 
    ## Run 345 stress 0.2494706 
    ## Run 346 stress 0.0003777456 
    ## ... Procrustes: rmse 0.01409847  max resid 0.02062025 
    ## Run 347 stress 0.2842805 
    ## Run 348 stress 0.1491957 
    ## Run 349 stress 9.350834e-05 
    ## ... Procrustes: rmse 0.0005928928  max resid 0.0009856847 
    ## ... Similar to previous best
    ## Run 350 stress 0.0004866398 
    ## ... Procrustes: rmse 0.0160069  max resid 0.02325649 
    ## Run 351 stress 0.0002064806 
    ## ... Procrustes: rmse 0.01039333  max resid 0.01550147 
    ## Run 352 stress 9.070496e-05 
    ## ... Procrustes: rmse 0.0005857713  max resid 0.0008266501 
    ## ... Similar to previous best
    ## Run 353 stress 0.2570116 
    ## Run 354 stress 0.2570117 
    ## Run 355 stress 9.532605e-05 
    ## ... Procrustes: rmse 0.0006074641  max resid 0.001044828 
    ## ... Similar to previous best
    ## Run 356 stress 8.596434e-05 
    ## ... Procrustes: rmse 0.003034337  max resid 0.005334954 
    ## ... Similar to previous best
    ## Run 357 stress 6.882101e-05 
    ## ... Procrustes: rmse 0.000623328  max resid 0.0009513587 
    ## ... Similar to previous best
    ## Run 358 stress 9.492054e-05 
    ## ... Procrustes: rmse 0.0006693604  max resid 0.0009574911 
    ## ... Similar to previous best
    ## Run 359 stress 0.000670305 
    ## Run 360 stress 9.931497e-05 
    ## ... Procrustes: rmse 0.000678428  max resid 0.001081963 
    ## ... Similar to previous best
    ## Run 361 stress 0.1491957 
    ## Run 362 stress 8.73552e-05 
    ## ... Procrustes: rmse 0.0006359989  max resid 0.0009673311 
    ## ... Similar to previous best
    ## Run 363 stress 0.0005004938 
    ## ... Procrustes: rmse 0.01623351  max resid 0.02356894 
    ## Run 364 stress 0.0004873645 
    ## ... Procrustes: rmse 0.01603884  max resid 0.02329981 
    ## Run 365 stress 0.3083098 
    ## Run 366 stress 8.608502e-05 
    ## ... Procrustes: rmse 0.0005792874  max resid 0.0008852587 
    ## ... Similar to previous best
    ## Run 367 stress 8.89628e-05 
    ## ... Procrustes: rmse 0.0005494066  max resid 0.0008025937 
    ## ... Similar to previous best
    ## Run 368 stress 0.001281742 
    ## Run 369 stress 0.000455609 
    ## ... Procrustes: rmse 0.01549363  max resid 0.02255277 
    ## Run 370 stress 8.548161e-05 
    ## ... Procrustes: rmse 0.0008055904  max resid 0.001159342 
    ## ... Similar to previous best
    ## Run 371 stress 0.0004634072 
    ## ... Procrustes: rmse 0.01563671  max resid 0.02274954 
    ## Run 372 stress 9.552891e-05 
    ## ... Procrustes: rmse 0.0005241592  max resid 0.0007883599 
    ## ... Similar to previous best
    ## Run 373 stress 9.539769e-05 
    ## ... Procrustes: rmse 0.0005631123  max resid 0.0008883381 
    ## ... Similar to previous best
    ## Run 374 stress 9.130525e-05 
    ## ... Procrustes: rmse 0.000562208  max resid 0.0008895759 
    ## ... Similar to previous best
    ## Run 375 stress 0.0005189 
    ## ... Procrustes: rmse 0.01655599  max resid 0.02401394 
    ## Run 376 stress 0.001284389 
    ## Run 377 stress 0.1990774 
    ## Run 378 stress 0.1491957 
    ## Run 379 stress 6.208538e-05 
    ## ... Procrustes: rmse 0.0005878124  max resid 0.0007928305 
    ## ... Similar to previous best
    ## Run 380 stress 0.3082926 
    ## Run 381 stress 9.402873e-05 
    ## ... Procrustes: rmse 0.0006511237  max resid 0.001043364 
    ## ... Similar to previous best
    ## Run 382 stress 0.2696744 
    ## Run 383 stress 0.1990774 
    ## Run 384 stress 0.001139338 
    ## Run 385 stress 0.001271494 
    ## Run 386 stress 7.800435e-05 
    ## ... Procrustes: rmse 0.0002523991  max resid 0.0003506221 
    ## ... Similar to previous best
    ## Run 387 stress 0.001193047 
    ## Run 388 stress 0.2440272 
    ## Run 389 stress 0.0004932454 
    ## ... Procrustes: rmse 0.01613977  max resid 0.02344084 
    ## Run 390 stress 9.212236e-05 
    ## ... Procrustes: rmse 0.0005327906  max resid 0.0007974108 
    ## ... Similar to previous best
    ## Run 391 stress 0.00047242 
    ## ... Procrustes: rmse 0.01578387  max resid 0.02294797 
    ## Run 392 stress 0.000519665 
    ## ... Procrustes: rmse 0.01655954  max resid 0.02401898 
    ## Run 393 stress 0.1990774 
    ## Run 394 stress 9.375008e-05 
    ## ... Procrustes: rmse 0.000570948  max resid 0.0008543032 
    ## ... Similar to previous best
    ## Run 395 stress 0.001317369 
    ## Run 396 stress 9.772493e-05 
    ## ... Procrustes: rmse 0.0005747835  max resid 0.0008766222 
    ## ... Similar to previous best
    ## Run 397 stress 8.993531e-05 
    ## ... Procrustes: rmse 0.0005648216  max resid 0.0008872171 
    ## ... Similar to previous best
    ## Run 398 stress 9.853346e-05 
    ## ... Procrustes: rmse 0.000574882  max resid 0.0009015538 
    ## ... Similar to previous best
    ## Run 399 stress 0.0005545485 
    ## Run 400 stress 9.044378e-05 
    ## ... Procrustes: rmse 0.0005894269  max resid 0.0008825989 
    ## ... Similar to previous best
    ## Run 401 stress 9.30066e-05 
    ## ... Procrustes: rmse 0.0005280045  max resid 0.0007953926 
    ## ... Similar to previous best
    ## Run 402 stress 9.458258e-05 
    ## ... Procrustes: rmse 0.0006359683  max resid 0.0009888476 
    ## ... Similar to previous best
    ## Run 403 stress 9.796635e-05 
    ## ... Procrustes: rmse 0.0006096037  max resid 0.001022965 
    ## ... Similar to previous best
    ## Run 404 stress 9.709098e-05 
    ## ... Procrustes: rmse 0.0005874479  max resid 0.0008605598 
    ## ... Similar to previous best
    ## Run 405 stress 0.1990774 
    ## Run 406 stress 8.346169e-05 
    ## ... Procrustes: rmse 0.0006670139  max resid 0.0009700702 
    ## ... Similar to previous best
    ## Run 407 stress 0.001163103 
    ## Run 408 stress 8.536289e-05 
    ## ... Procrustes: rmse 0.001404743  max resid 0.002918045 
    ## ... Similar to previous best
    ## Run 409 stress 0.0004845202 
    ## ... Procrustes: rmse 0.01595674  max resid 0.02318652 
    ## Run 410 stress 9.577255e-05 
    ## ... Procrustes: rmse 0.0005238136  max resid 0.0007935736 
    ## ... Similar to previous best
    ## Run 411 stress 0.1491957 
    ## Run 412 stress 0.001395196 
    ## Run 413 stress 0.001398306 
    ## Run 414 stress 0.1491957 
    ## Run 415 stress 8.930153e-05 
    ## ... Procrustes: rmse 0.0005345759  max resid 0.0007864836 
    ## ... Similar to previous best
    ## Run 416 stress 0.1491957 
    ## Run 417 stress 0.001393743 
    ## Run 418 stress 0.00125099 
    ## Run 419 stress 9.445055e-05 
    ## ... Procrustes: rmse 0.0005791321  max resid 0.0008158135 
    ## ... Similar to previous best
    ## Run 420 stress 9.922233e-05 
    ## ... Procrustes: rmse 0.0005566454  max resid 0.0008781308 
    ## ... Similar to previous best
    ## Run 421 stress 0.1491957 
    ## Run 422 stress 0.001306441 
    ## Run 423 stress 0.1491957 
    ## Run 424 stress 0.1776011 
    ## Run 425 stress 9.607802e-05 
    ## ... Procrustes: rmse 0.000590099  max resid 0.0009962235 
    ## ... Similar to previous best
    ## Run 426 stress 0.001231136 
    ## Run 427 stress 0.0004569993 
    ## ... Procrustes: rmse 0.01553108  max resid 0.02259919 
    ## Run 428 stress 0.1491957 
    ## Run 429 stress 0.00126141 
    ## Run 430 stress 0.0002669531 
    ## ... Procrustes: rmse 0.01909646  max resid 0.02633936 
    ## Run 431 stress 9.9976e-05 
    ## ... Procrustes: rmse 0.007196575  max resid 0.01107806 
    ## Run 432 stress 9.056962e-05 
    ## ... Procrustes: rmse 0.0006307367  max resid 0.0009751972 
    ## ... Similar to previous best
    ## Run 433 stress 9.263764e-05 
    ## ... Procrustes: rmse 0.0006498589  max resid 0.001029452 
    ## ... Similar to previous best
    ## Run 434 stress 7.151392e-05 
    ## ... Procrustes: rmse 0.0006167318  max resid 0.0009737803 
    ## ... Similar to previous best
    ## Run 435 stress 0.1491957 
    ## Run 436 stress 0.0001575604 
    ## ... Procrustes: rmse 0.01451751  max resid 0.02000936 
    ## Run 437 stress 0.001346401 
    ## Run 438 stress 8.628885e-05 
    ## ... Procrustes: rmse 0.0005357266  max resid 0.0007920519 
    ## ... Similar to previous best
    ## Run 439 stress 0.2578649 
    ## Run 440 stress 0.001227478 
    ## Run 441 stress 0.001259591 
    ## Run 442 stress 0.2520602 
    ## Run 443 stress 9.797437e-05 
    ## ... Procrustes: rmse 0.0006359834  max resid 0.0009984593 
    ## ... Similar to previous best
    ## Run 444 stress 8.088147e-05 
    ## ... Procrustes: rmse 0.0005229243  max resid 0.0007180341 
    ## ... Similar to previous best
    ## Run 445 stress 8.430435e-05 
    ## ... Procrustes: rmse 0.0006526143  max resid 0.001003203 
    ## ... Similar to previous best
    ## Run 446 stress 0.1491957 
    ## Run 447 stress 0.1776011 
    ## Run 448 stress 8.584139e-05 
    ## ... Procrustes: rmse 0.0006009068  max resid 0.0009692157 
    ## ... Similar to previous best
    ## Run 449 stress 0.001211073 
    ## Run 450 stress 0.000483219 
    ## ... Procrustes: rmse 0.01597433  max resid 0.02321044 
    ## Run 451 stress 0.001396748 
    ## Run 452 stress 0.1491957 
    ## Run 453 stress 0.00106313 
    ## Run 454 stress 0.1491957 
    ## Run 455 stress 9.059444e-05 
    ## ... Procrustes: rmse 0.0005949957  max resid 0.0009904755 
    ## ... Similar to previous best
    ## Run 456 stress 0.001308796 
    ## Run 457 stress 0.00133497 
    ## Run 458 stress 0.1491957 
    ## Run 459 stress 9.676872e-05 
    ## ... Procrustes: rmse 0.0005562658  max resid 0.0008841287 
    ## ... Similar to previous best
    ## Run 460 stress 0.0001629171 
    ## ... Procrustes: rmse 0.0147733  max resid 0.02036282 
    ## Run 461 stress 0.1491957 
    ## Run 462 stress 0.0004671835 
    ## ... Procrustes: rmse 0.0156996  max resid 0.02283166 
    ## Run 463 stress 0.000433161 
    ## ... Procrustes: rmse 0.01508273  max resid 0.02197937 
    ## Run 464 stress 0.0003993262 
    ## ... Procrustes: rmse 0.01450894  max resid 0.0211874 
    ## Run 465 stress 0.001440189 
    ## Run 466 stress 9.209607e-05 
    ## ... Procrustes: rmse 0.0005770996  max resid 0.00085142 
    ## ... Similar to previous best
    ## Run 467 stress 9.010989e-05 
    ## ... Procrustes: rmse 0.0005345428  max resid 0.0007855372 
    ## ... Similar to previous best
    ## Run 468 stress 0.1491957 
    ## Run 469 stress 9.744524e-05 
    ## ... Procrustes: rmse 0.0005760418  max resid 0.000811907 
    ## ... Similar to previous best
    ## Run 470 stress 0.1491957 
    ## Run 471 stress 0.001362679 
    ## Run 472 stress 0.001171343 
    ## Run 473 stress 8.344909e-05 
    ## ... Procrustes: rmse 0.0006118769  max resid 0.0009961519 
    ## ... Similar to previous best
    ## Run 474 stress 9.474217e-05 
    ## ... Procrustes: rmse 0.0006094232  max resid 0.00100881 
    ## ... Similar to previous best
    ## Run 475 stress 0.1491957 
    ## Run 476 stress 0.2528294 
    ## Run 477 stress 0.001416678 
    ## Run 478 stress 8.633759e-05 
    ## ... Procrustes: rmse 0.0005982529  max resid 0.000988651 
    ## ... Similar to previous best
    ## Run 479 stress 0.001363327 
    ## Run 480 stress 0.001311765 
    ## Run 481 stress 0.00133819 
    ## Run 482 stress 0.0001880048 
    ## ... Procrustes: rmse 0.0159193  max resid 0.02194664 
    ## Run 483 stress 0.00135962 
    ## Run 484 stress 9.022385e-05 
    ## ... Procrustes: rmse 0.0006394403  max resid 0.001011942 
    ## ... Similar to previous best
    ## Run 485 stress 0.2842805 
    ## Run 486 stress 0.0012884 
    ## Run 487 stress 9.645008e-05 
    ## ... Procrustes: rmse 0.0005667542  max resid 0.0008962971 
    ## ... Similar to previous best
    ## Run 488 stress 0.001186421 
    ## Run 489 stress 9.487847e-05 
    ## ... Procrustes: rmse 0.0006117908  max resid 0.001012437 
    ## ... Similar to previous best
    ## Run 490 stress 0.1491957 
    ## Run 491 stress 8.761725e-05 
    ## ... Procrustes: rmse 0.001399305  max resid 0.001992714 
    ## ... Similar to previous best
    ## Run 492 stress 8.62345e-05 
    ## ... Procrustes: rmse 0.0005295765  max resid 0.0007270486 
    ## ... Similar to previous best
    ## Run 493 stress 0.0008213447 
    ## Run 494 stress 9.136179e-05 
    ## ... Procrustes: rmse 0.0005958173  max resid 0.0009998648 
    ## ... Similar to previous best
    ## Run 495 stress 0.0001044246 
    ## ... Procrustes: rmse 0.01169292  max resid 0.01610656 
    ## Run 496 stress 0.0008594191 
    ## Run 497 stress 0.0004534447 
    ## ... Procrustes: rmse 0.02508292  max resid 0.03462478 
    ## Run 498 stress 6.372434e-05 
    ## ... Procrustes: rmse 0.0005659291  max resid 0.0007885791 
    ## ... Similar to previous best
    ## Run 499 stress 0.0004457704 
    ## ... Procrustes: rmse 0.01533683  max resid 0.02233086 
    ## Run 500 stress 0.0005377325 
    ## ... Procrustes: rmse 0.01664951  max resid 0.02414419 
    ## *** Best solution repeated 140 times

    ## Warning in metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
SD_beta_geo_NMDS <- metaMDS(SD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09021173 
    ## Run 1 stress 0.09088611 
    ## Run 2 stress 0.0892608 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04368274  max resid 0.1170431 
    ## Run 3 stress 0.09088603 
    ## Run 4 stress 0.1087569 
    ## Run 5 stress 0.1060095 
    ## Run 6 stress 0.08938554 
    ## ... Procrustes: rmse 0.01174954  max resid 0.03979394 
    ## Run 7 stress 0.09088614 
    ## Run 8 stress 0.1075422 
    ## Run 9 stress 0.106256 
    ## Run 10 stress 0.1060429 
    ## Run 11 stress 0.1092207 
    ## Run 12 stress 0.09018709 
    ## Run 13 stress 0.1075421 
    ## Run 14 stress 0.08946329 
    ## ... Procrustes: rmse 0.03754812  max resid 0.1178263 
    ## Run 15 stress 0.0913009 
    ## Run 16 stress 0.08938537 
    ## ... Procrustes: rmse 0.01179919  max resid 0.03973672 
    ## Run 17 stress 0.08938538 
    ## ... Procrustes: rmse 0.01181823  max resid 0.03976266 
    ## Run 18 stress 0.1079016 
    ## Run 19 stress 0.09088601 
    ## Run 20 stress 0.1056906 
    ## Run 21 stress 0.0959214 
    ## Run 22 stress 0.08926079 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004415205  max resid 0.001152624 
    ## ... Similar to previous best
    ## Run 23 stress 0.08926076 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004082817  max resid 0.001015204 
    ## ... Similar to previous best
    ## Run 24 stress 0.1074435 
    ## Run 25 stress 0.1062564 
    ## Run 26 stress 0.1062566 
    ## Run 27 stress 0.09178302 
    ## Run 28 stress 0.1071321 
    ## Run 29 stress 0.08951719 
    ## ... Procrustes: rmse 0.03511848  max resid 0.1157408 
    ## Run 30 stress 0.08938538 
    ## ... Procrustes: rmse 0.01179211  max resid 0.03973031 
    ## Run 31 stress 0.09775562 
    ## Run 32 stress 0.1074434 
    ## Run 33 stress 0.1074432 
    ## Run 34 stress 0.1074432 
    ## Run 35 stress 0.1110405 
    ## Run 36 stress 0.08926093 
    ## ... Procrustes: rmse 0.0001395361  max resid 0.0004693629 
    ## ... Similar to previous best
    ## Run 37 stress 0.09087357 
    ## Run 38 stress 0.09018707 
    ## Run 39 stress 0.08938538 
    ## ... Procrustes: rmse 0.01180829  max resid 0.03974953 
    ## Run 40 stress 0.08938539 
    ## ... Procrustes: rmse 0.01182873  max resid 0.03977516 
    ## Run 41 stress 0.08938992 
    ## ... Procrustes: rmse 0.03594057  max resid 0.1182713 
    ## Run 42 stress 0.1061299 
    ## Run 43 stress 0.09503457 
    ## Run 44 stress 0.1067382 
    ## Run 45 stress 0.08926072 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003828204  max resid 0.001033976 
    ## ... Similar to previous best
    ## Run 46 stress 0.09503465 
    ## Run 47 stress 0.09262362 
    ## Run 48 stress 0.1104177 
    ## Run 49 stress 0.08938961 
    ## ... Procrustes: rmse 0.03589958  max resid 0.1182269 
    ## Run 50 stress 0.09180901 
    ## Run 51 stress 0.09018707 
    ## Run 52 stress 0.08926087 
    ## ... Procrustes: rmse 0.000470601  max resid 0.001333006 
    ## ... Similar to previous best
    ## Run 53 stress 0.08926077 
    ## ... Procrustes: rmse 0.0004077344  max resid 0.001015006 
    ## ... Similar to previous best
    ## Run 54 stress 0.09039141 
    ## Run 55 stress 0.08938565 
    ## ... Procrustes: rmse 0.01170819  max resid 0.0395279 
    ## Run 56 stress 0.08951734 
    ## ... Procrustes: rmse 0.03495954  max resid 0.1155518 
    ## Run 57 stress 0.09088603 
    ## Run 58 stress 0.08926068 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002055539  max resid 0.0005479854 
    ## ... Similar to previous best
    ## Run 59 stress 0.08926067 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001013846  max resid 0.0002414493 
    ## ... Similar to previous best
    ## Run 60 stress 0.08938539 
    ## ... Procrustes: rmse 0.01188781  max resid 0.03978064 
    ## Run 61 stress 0.09088603 
    ## Run 62 stress 0.08938546 
    ## ... Procrustes: rmse 0.01183519  max resid 0.03980926 
    ## Run 63 stress 0.09594332 
    ## Run 64 stress 0.1052653 
    ## Run 65 stress 0.1074432 
    ## Run 66 stress 0.1091824 
    ## Run 67 stress 0.1076303 
    ## Run 68 stress 0.0894633 
    ## ... Procrustes: rmse 0.03752513  max resid 0.1177767 
    ## Run 69 stress 0.109256 
    ## Run 70 stress 0.08946652 
    ## ... Procrustes: rmse 0.03322894  max resid 0.1163531 
    ## Run 71 stress 0.09087339 
    ## Run 72 stress 0.08951722 
    ## ... Procrustes: rmse 0.03511441  max resid 0.1157031 
    ## Run 73 stress 0.1064188 
    ## Run 74 stress 0.09044609 
    ## Run 75 stress 0.1111371 
    ## Run 76 stress 0.1056902 
    ## Run 77 stress 0.1068388 
    ## Run 78 stress 0.09021168 
    ## Run 79 stress 0.09130086 
    ## Run 80 stress 0.08946336 
    ## ... Procrustes: rmse 0.03750658  max resid 0.1177503 
    ## Run 81 stress 0.09130108 
    ## Run 82 stress 0.08946656 
    ## ... Procrustes: rmse 0.03321873  max resid 0.1163387 
    ## Run 83 stress 0.1075422 
    ## Run 84 stress 0.08938553 
    ## ... Procrustes: rmse 0.01178889  max resid 0.0397222 
    ## Run 85 stress 0.1095612 
    ## Run 86 stress 0.1087561 
    ## Run 87 stress 0.10613 
    ## Run 88 stress 0.09503426 
    ## Run 89 stress 0.1074435 
    ## Run 90 stress 0.08938542 
    ## ... Procrustes: rmse 0.01184161  max resid 0.03972544 
    ## Run 91 stress 0.1080635 
    ## Run 92 stress 0.09262358 
    ## Run 93 stress 0.09039136 
    ## Run 94 stress 0.09021167 
    ## Run 95 stress 0.1060429 
    ## Run 96 stress 0.1074433 
    ## Run 97 stress 0.09021173 
    ## Run 98 stress 0.1056903 
    ## Run 99 stress 0.08926069 
    ## ... Procrustes: rmse 0.0001505875  max resid 0.0003667883 
    ## ... Similar to previous best
    ## Run 100 stress 0.1074433 
    ## Run 101 stress 0.08926065 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001971941  max resid 0.0004896209 
    ## ... Similar to previous best
    ## Run 102 stress 0.09088622 
    ## Run 103 stress 0.08938546 
    ## ... Procrustes: rmse 0.0118811  max resid 0.03958135 
    ## Run 104 stress 0.08938554 
    ## ... Procrustes: rmse 0.01178369  max resid 0.03967632 
    ## Run 105 stress 0.09178295 
    ## Run 106 stress 0.1108636 
    ## Run 107 stress 0.1075423 
    ## Run 108 stress 0.09021165 
    ## Run 109 stress 0.1076304 
    ## Run 110 stress 0.1087865 
    ## Run 111 stress 0.1067909 
    ## Run 112 stress 0.1063192 
    ## Run 113 stress 0.09503421 
    ## Run 114 stress 0.08926065 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002191622  max resid 0.000673147 
    ## ... Similar to previous best
    ## Run 115 stress 0.08946656 
    ## ... Procrustes: rmse 0.03322249  max resid 0.1163483 
    ## Run 116 stress 0.08938542 
    ## ... Procrustes: rmse 0.01177925  max resid 0.03967171 
    ## Run 117 stress 0.09130099 
    ## Run 118 stress 0.1071321 
    ## Run 119 stress 0.1120157 
    ## Run 120 stress 0.08938971 
    ## ... Procrustes: rmse 0.03593409  max resid 0.1182587 
    ## Run 121 stress 0.09109108 
    ## Run 122 stress 0.1060094 
    ## Run 123 stress 0.09130088 
    ## Run 124 stress 0.08926078 
    ## ... Procrustes: rmse 0.0001488863  max resid 0.0004604126 
    ## ... Similar to previous best
    ## Run 125 stress 0.1063193 
    ## Run 126 stress 0.09503436 
    ## Run 127 stress 0.1092126 
    ## Run 128 stress 0.09109109 
    ## Run 129 stress 0.08938962 
    ## ... Procrustes: rmse 0.0359546  max resid 0.1182876 
    ## Run 130 stress 0.1062562 
    ## Run 131 stress 0.0959215 
    ## Run 132 stress 0.1060435 
    ## Run 133 stress 0.090996 
    ## Run 134 stress 0.08938539 
    ## ... Procrustes: rmse 0.01179781  max resid 0.03969175 
    ## Run 135 stress 0.09503461 
    ## Run 136 stress 0.09099561 
    ## Run 137 stress 0.110534 
    ## Run 138 stress 0.09503434 
    ## Run 139 stress 0.1062599 
    ## Run 140 stress 0.09018707 
    ## Run 141 stress 0.08946331 
    ## ... Procrustes: rmse 0.03751292  max resid 0.1177647 
    ## Run 142 stress 0.1091237 
    ## Run 143 stress 0.1071321 
    ## Run 144 stress 0.09109112 
    ## Run 145 stress 0.1087564 
    ## Run 146 stress 0.105265 
    ## Run 147 stress 0.08926084 
    ## ... Procrustes: rmse 0.0003936121  max resid 0.00119054 
    ## ... Similar to previous best
    ## Run 148 stress 0.08926082 
    ## ... Procrustes: rmse 0.0003914526  max resid 0.001219646 
    ## ... Similar to previous best
    ## Run 149 stress 0.08926065 
    ## ... New best solution
    ## ... Procrustes: rmse 7.176882e-05  max resid 0.0001911632 
    ## ... Similar to previous best
    ## Run 150 stress 0.1060088 
    ## Run 151 stress 0.08946652 
    ## ... Procrustes: rmse 0.03320473  max resid 0.1163279 
    ## Run 152 stress 0.09088623 
    ## Run 153 stress 0.09503437 
    ## Run 154 stress 0.09503455 
    ## Run 155 stress 0.09044608 
    ## Run 156 stress 0.1065371 
    ## Run 157 stress 0.08946346 
    ## ... Procrustes: rmse 0.0374618  max resid 0.1176933 
    ## Run 158 stress 0.08926081 
    ## ... Procrustes: rmse 0.0002899436  max resid 0.0008025362 
    ## ... Similar to previous best
    ## Run 159 stress 0.09044616 
    ## Run 160 stress 0.08938976 
    ## ... Procrustes: rmse 0.03591114  max resid 0.1182287 
    ## Run 161 stress 0.08938547 
    ## ... Procrustes: rmse 0.01172277  max resid 0.03953689 
    ## Run 162 stress 0.1091878 
    ## Run 163 stress 0.09088604 
    ## Run 164 stress 0.09503457 
    ## Run 165 stress 0.09021166 
    ## Run 166 stress 0.1074433 
    ## Run 167 stress 0.08946651 
    ## ... Procrustes: rmse 0.03322582  max resid 0.1163683 
    ## Run 168 stress 0.1074432 
    ## Run 169 stress 0.1071322 
    ## Run 170 stress 0.08938972 
    ## ... Procrustes: rmse 0.0359171  max resid 0.1182382 
    ## Run 171 stress 0.1052654 
    ## Run 172 stress 0.091783 
    ## Run 173 stress 0.08938962 
    ## ... Procrustes: rmse 0.03596557  max resid 0.1183194 
    ## Run 174 stress 0.1056908 
    ## Run 175 stress 0.08946652 
    ## ... Procrustes: rmse 0.03320861  max resid 0.1163388 
    ## Run 176 stress 0.1076305 
    ## Run 177 stress 0.09039129 
    ## Run 178 stress 0.1062565 
    ## Run 179 stress 0.09018712 
    ## Run 180 stress 0.08938963 
    ## ... Procrustes: rmse 0.03596752  max resid 0.1183239 
    ## Run 181 stress 0.1091825 
    ## Run 182 stress 0.08951718 
    ## ... Procrustes: rmse 0.03505009  max resid 0.1156675 
    ## Run 183 stress 0.08926064 
    ## ... New best solution
    ## ... Procrustes: rmse 6.976754e-05  max resid 0.0002165311 
    ## ... Similar to previous best
    ## Run 184 stress 0.1088391 
    ## Run 185 stress 0.08946651 
    ## ... Procrustes: rmse 0.03321166  max resid 0.1163542 
    ## Run 186 stress 0.08938965 
    ## ... Procrustes: rmse 0.0359205  max resid 0.1182465 
    ## Run 187 stress 0.08938965 
    ## ... Procrustes: rmse 0.0359239  max resid 0.1182537 
    ## Run 188 stress 0.09018708 
    ## Run 189 stress 0.09018708 
    ## Run 190 stress 0.08938538 
    ## ... Procrustes: rmse 0.01179806  max resid 0.03959828 
    ## Run 191 stress 0.1060433 
    ## Run 192 stress 0.08938972 
    ## ... Procrustes: rmse 0.03590715  max resid 0.1182233 
    ## Run 193 stress 0.09109108 
    ## Run 194 stress 0.09503416 
    ## Run 195 stress 0.09592149 
    ## Run 196 stress 0.08926105 
    ## ... Procrustes: rmse 0.000436532  max resid 0.001288974 
    ## ... Similar to previous best
    ## Run 197 stress 0.1061298 
    ## Run 198 stress 0.09099593 
    ## Run 199 stress 0.1091881 
    ## Run 200 stress 0.1064191 
    ## Run 201 stress 0.09018708 
    ## Run 202 stress 0.08926099 
    ## ... Procrustes: rmse 0.0003166688  max resid 0.0009275163 
    ## ... Similar to previous best
    ## Run 203 stress 0.108657 
    ## Run 204 stress 0.1060437 
    ## Run 205 stress 0.09503458 
    ## Run 206 stress 0.08938963 
    ## ... Procrustes: rmse 0.03595943  max resid 0.1183113 
    ## Run 207 stress 0.09018707 
    ## Run 208 stress 0.08951718 
    ## ... Procrustes: rmse 0.03502474  max resid 0.1156136 
    ## Run 209 stress 0.09503414 
    ## Run 210 stress 0.1105336 
    ## Run 211 stress 0.10613 
    ## Run 212 stress 0.1075421 
    ## Run 213 stress 0.1064187 
    ## Run 214 stress 0.09039137 
    ## Run 215 stress 0.08926076 
    ## ... Procrustes: rmse 0.000204563  max resid 0.0005895915 
    ## ... Similar to previous best
    ## Run 216 stress 0.09039132 
    ## Run 217 stress 0.1104481 
    ## Run 218 stress 0.08938964 
    ## ... Procrustes: rmse 0.03592223  max resid 0.1182497 
    ## Run 219 stress 0.1071321 
    ## Run 220 stress 0.09109113 
    ## Run 221 stress 0.08926098 
    ## ... Procrustes: rmse 0.0004022826  max resid 0.001191782 
    ## ... Similar to previous best
    ## Run 222 stress 0.1067389 
    ## Run 223 stress 0.1064187 
    ## Run 224 stress 0.09503417 
    ## Run 225 stress 0.09018712 
    ## Run 226 stress 0.1074433 
    ## Run 227 stress 0.1092093 
    ## Run 228 stress 0.08938539 
    ## ... Procrustes: rmse 0.01182287  max resid 0.03962976 
    ## Run 229 stress 0.08938538 
    ## ... Procrustes: rmse 0.01181254  max resid 0.03958321 
    ## Run 230 stress 0.09109108 
    ## Run 231 stress 0.1104193 
    ## Run 232 stress 0.09039131 
    ## Run 233 stress 0.1074434 
    ## Run 234 stress 0.09018707 
    ## Run 235 stress 0.09180904 
    ## Run 236 stress 0.08951721 
    ## ... Procrustes: rmse 0.03501767  max resid 0.1156073 
    ## Run 237 stress 0.1067487 
    ## Run 238 stress 0.08938967 
    ## ... Procrustes: rmse 0.03596466  max resid 0.1183268 
    ## Run 239 stress 0.08938539 
    ## ... Procrustes: rmse 0.01181807  max resid 0.03962258 
    ## Run 240 stress 0.09094202 
    ## Run 241 stress 0.1071321 
    ## Run 242 stress 0.09109112 
    ## Run 243 stress 0.1092129 
    ## Run 244 stress 0.1052647 
    ## Run 245 stress 0.0926236 
    ## Run 246 stress 0.1068392 
    ## Run 247 stress 0.1067886 
    ## Run 248 stress 0.1074433 
    ## Run 249 stress 0.1052647 
    ## Run 250 stress 0.1056902 
    ## Run 251 stress 0.09180876 
    ## Run 252 stress 0.08926072 
    ## ... Procrustes: rmse 0.000213553  max resid 0.0005888313 
    ## ... Similar to previous best
    ## Run 253 stress 0.0894666 
    ## ... Procrustes: rmse 0.03317469  max resid 0.1162809 
    ## Run 254 stress 0.1067581 
    ## Run 255 stress 0.08926076 
    ## ... Procrustes: rmse 0.0001624807  max resid 0.0003931574 
    ## ... Similar to previous best
    ## Run 256 stress 0.09039142 
    ## Run 257 stress 0.08926078 
    ## ... Procrustes: rmse 0.0002583172  max resid 0.0007318191 
    ## ... Similar to previous best
    ## Run 258 stress 0.1060435 
    ## Run 259 stress 0.08926064 
    ## ... Procrustes: rmse 9.408127e-05  max resid 0.0002221058 
    ## ... Similar to previous best
    ## Run 260 stress 0.1062567 
    ## Run 261 stress 0.1092093 
    ## Run 262 stress 0.09018713 
    ## Run 263 stress 0.08938968 
    ## ... Procrustes: rmse 0.03591511  max resid 0.118239 
    ## Run 264 stress 0.08938542 
    ## ... Procrustes: rmse 0.01176095  max resid 0.03954643 
    ## Run 265 stress 0.08926065 
    ## ... Procrustes: rmse 6.862937e-05  max resid 0.0002288574 
    ## ... Similar to previous best
    ## Run 266 stress 0.09775593 
    ## Run 267 stress 0.1092093 
    ## Run 268 stress 0.1062597 
    ## Run 269 stress 0.08926071 
    ## ... Procrustes: rmse 0.0001873189  max resid 0.0006351474 
    ## ... Similar to previous best
    ## Run 270 stress 0.08938973 
    ## ... Procrustes: rmse 0.03597901  max resid 0.1183398 
    ## Run 271 stress 0.1065377 
    ## Run 272 stress 0.09021166 
    ## Run 273 stress 0.08938962 
    ## ... Procrustes: rmse 0.03592978  max resid 0.1182608 
    ## Run 274 stress 0.08951724 
    ## ... Procrustes: rmse 0.035012  max resid 0.1156037 
    ## Run 275 stress 0.08946331 
    ## ... Procrustes: rmse 0.03747906  max resid 0.1177287 
    ## Run 276 stress 0.1111368 
    ## Run 277 stress 0.105265 
    ## Run 278 stress 0.08946328 
    ## ... Procrustes: rmse 0.03750355  max resid 0.1177472 
    ## Run 279 stress 0.0894665 
    ## ... Procrustes: rmse 0.03320776  max resid 0.11634 
    ## Run 280 stress 0.08926091 
    ## ... Procrustes: rmse 0.0003602319  max resid 0.001070721 
    ## ... Similar to previous best
    ## Run 281 stress 0.09178293 
    ## Run 282 stress 0.1071321 
    ## Run 283 stress 0.1093311 
    ## Run 284 stress 0.1056901 
    ## Run 285 stress 0.0926236 
    ## Run 286 stress 0.1092126 
    ## Run 287 stress 0.0917829 
    ## Run 288 stress 0.08938561 
    ## ... Procrustes: rmse 0.01186983  max resid 0.03951821 
    ## Run 289 stress 0.1074434 
    ## Run 290 stress 0.09018707 
    ## Run 291 stress 0.09099558 
    ## Run 292 stress 0.08938545 
    ## ... Procrustes: rmse 0.01184431  max resid 0.03956663 
    ## Run 293 stress 0.106043 
    ## Run 294 stress 0.09503414 
    ## Run 295 stress 0.09503428 
    ## Run 296 stress 0.08938541 
    ## ... Procrustes: rmse 0.01177378  max resid 0.03960585 
    ## Run 297 stress 0.09592165 
    ## Run 298 stress 0.09503424 
    ## Run 299 stress 0.08946328 
    ## ... Procrustes: rmse 0.03749763  max resid 0.1177597 
    ## Run 300 stress 0.08926091 
    ## ... Procrustes: rmse 0.0003560971  max resid 0.001078489 
    ## ... Similar to previous best
    ## Run 301 stress 0.08926086 
    ## ... Procrustes: rmse 0.0003296258  max resid 0.0009274905 
    ## ... Similar to previous best
    ## Run 302 stress 0.09109113 
    ## Run 303 stress 0.1071322 
    ## Run 304 stress 0.1091453 
    ## Run 305 stress 0.09228507 
    ## Run 306 stress 0.09044614 
    ## Run 307 stress 0.08938539 
    ## ... Procrustes: rmse 0.0118254  max resid 0.03956925 
    ## Run 308 stress 0.1074434 
    ## Run 309 stress 0.09021167 
    ## Run 310 stress 0.09088609 
    ## Run 311 stress 0.1108451 
    ## Run 312 stress 0.09087344 
    ## Run 313 stress 0.09592158 
    ## Run 314 stress 0.08946341 
    ## ... Procrustes: rmse 0.03745782  max resid 0.1177015 
    ## Run 315 stress 0.08926065 
    ## ... Procrustes: rmse 9.27203e-05  max resid 0.0002262976 
    ## ... Similar to previous best
    ## Run 316 stress 0.1091825 
    ## Run 317 stress 0.1052647 
    ## Run 318 stress 0.0903913 
    ## Run 319 stress 0.08926077 
    ## ... Procrustes: rmse 0.0002653047  max resid 0.0007574933 
    ## ... Similar to previous best
    ## Run 320 stress 0.09018713 
    ## Run 321 stress 0.08938544 
    ## ... Procrustes: rmse 0.01181745  max resid 0.03955729 
    ## Run 322 stress 0.09503452 
    ## Run 323 stress 0.1071321 
    ## Run 324 stress 0.09262357 
    ## Run 325 stress 0.09503422 
    ## Run 326 stress 0.08938966 
    ## ... Procrustes: rmse 0.03596243  max resid 0.1183228 
    ## Run 327 stress 0.09099543 
    ## Run 328 stress 0.09044608 
    ## Run 329 stress 0.09130091 
    ## Run 330 stress 0.08938967 
    ## ... Procrustes: rmse 0.03592622  max resid 0.1182516 
    ## Run 331 stress 0.1086565 
    ## Run 332 stress 0.08938962 
    ## ... Procrustes: rmse 0.03593027  max resid 0.1182596 
    ## Run 333 stress 0.08946651 
    ## ... Procrustes: rmse 0.03319625  max resid 0.1163205 
    ## Run 334 stress 0.1076304 
    ## Run 335 stress 0.09503453 
    ## Run 336 stress 0.1056899 
    ## Run 337 stress 0.0903913 
    ## Run 338 stress 0.09018709 
    ## Run 339 stress 0.09021169 
    ## Run 340 stress 0.1091827 
    ## Run 341 stress 0.1089906 
    ## Run 342 stress 0.08938553 
    ## ... Procrustes: rmse 0.01181766  max resid 0.03949307 
    ## Run 343 stress 0.1076303 
    ## Run 344 stress 0.1108233 
    ## Run 345 stress 0.08938964 
    ## ... Procrustes: rmse 0.03596181  max resid 0.1183155 
    ## Run 346 stress 0.09109108 
    ## Run 347 stress 0.09592146 
    ## Run 348 stress 0.08938961 
    ## ... Procrustes: rmse 0.03593935  max resid 0.1182782 
    ## Run 349 stress 0.1074435 
    ## Run 350 stress 0.08926064 
    ## ... Procrustes: rmse 9.381598e-05  max resid 0.0003002888 
    ## ... Similar to previous best
    ## Run 351 stress 0.08938542 
    ## ... Procrustes: rmse 0.01179865  max resid 0.03958947 
    ## Run 352 stress 0.1067506 
    ## Run 353 stress 0.1061303 
    ## Run 354 stress 0.08946654 
    ## ... Procrustes: rmse 0.03318735  max resid 0.1163033 
    ## Run 355 stress 0.1104482 
    ## Run 356 stress 0.1067906 
    ## Run 357 stress 0.1080639 
    ## Run 358 stress 0.09503419 
    ## Run 359 stress 0.1076304 
    ## Run 360 stress 0.1075421 
    ## Run 361 stress 0.1056905 
    ## Run 362 stress 0.1056897 
    ## Run 363 stress 0.0892607 
    ## ... Procrustes: rmse 0.0001902118  max resid 0.0005021817 
    ## ... Similar to previous best
    ## Run 364 stress 0.1056893 
    ## Run 365 stress 0.08946328 
    ## ... Procrustes: rmse 0.03749928  max resid 0.1177636 
    ## Run 366 stress 0.1079014 
    ## Run 367 stress 0.08926079 
    ## ... Procrustes: rmse 0.0002737098  max resid 0.0007714089 
    ## ... Similar to previous best
    ## Run 368 stress 0.1056897 
    ## Run 369 stress 0.1060096 
    ## Run 370 stress 0.08926073 
    ## ... Procrustes: rmse 0.0002280618  max resid 0.0005205559 
    ## ... Similar to previous best
    ## Run 371 stress 0.08938961 
    ## ... Procrustes: rmse 0.03593432  max resid 0.1182699 
    ## Run 372 stress 0.09592137 
    ## Run 373 stress 0.08938537 
    ## ... Procrustes: rmse 0.01180407  max resid 0.03959461 
    ## Run 374 stress 0.0893854 
    ## ... Procrustes: rmse 0.01183671  max resid 0.0395862 
    ## Run 375 stress 0.1056909 
    ## Run 376 stress 0.08926095 
    ## ... Procrustes: rmse 0.000307268  max resid 0.0008704923 
    ## ... Similar to previous best
    ## Run 377 stress 0.08926077 
    ## ... Procrustes: rmse 0.0002522333  max resid 0.0007219555 
    ## ... Similar to previous best
    ## Run 378 stress 0.1074432 
    ## Run 379 stress 0.09775586 
    ## Run 380 stress 0.09130086 
    ## Run 381 stress 0.08938539 
    ## ... Procrustes: rmse 0.01178363  max resid 0.03959058 
    ## Run 382 stress 0.1071322 
    ## Run 383 stress 0.08946651 
    ## ... Procrustes: rmse 0.03320666  max resid 0.1163369 
    ## Run 384 stress 0.1088368 
    ## Run 385 stress 0.08938542 
    ## ... Procrustes: rmse 0.01182626  max resid 0.03963705 
    ## Run 386 stress 0.09503432 
    ## Run 387 stress 0.08946653 
    ## ... Procrustes: rmse 0.03322028  max resid 0.1163731 
    ## Run 388 stress 0.0893854 
    ## ... Procrustes: rmse 0.01182424  max resid 0.03962696 
    ## Run 389 stress 0.09178284 
    ## Run 390 stress 0.109213 
    ## Run 391 stress 0.08926069 
    ## ... Procrustes: rmse 0.0001629322  max resid 0.0004121342 
    ## ... Similar to previous best
    ## Run 392 stress 0.08938966 
    ## ... Procrustes: rmse 0.03595201  max resid 0.118294 
    ## Run 393 stress 0.08926067 
    ## ... Procrustes: rmse 0.0001544943  max resid 0.0003823147 
    ## ... Similar to previous best
    ## Run 394 stress 0.08946331 
    ## ... Procrustes: rmse 0.03750133  max resid 0.11776 
    ## Run 395 stress 0.09109109 
    ## Run 396 stress 0.09592139 
    ## Run 397 stress 0.09228488 
    ## Run 398 stress 0.0902117 
    ## Run 399 stress 0.1063192 
    ## Run 400 stress 0.08926067 
    ## ... Procrustes: rmse 9.479059e-05  max resid 0.0002555615 
    ## ... Similar to previous best
    ## Run 401 stress 0.1060432 
    ## Run 402 stress 0.09039136 
    ## Run 403 stress 0.09503429 
    ## Run 404 stress 0.1101447 
    ## Run 405 stress 0.1092363 
    ## Run 406 stress 0.08946335 
    ## ... Procrustes: rmse 0.03747086  max resid 0.1177177 
    ## Run 407 stress 0.1104478 
    ## Run 408 stress 0.08926091 
    ## ... Procrustes: rmse 0.0003602086  max resid 0.001055348 
    ## ... Similar to previous best
    ## Run 409 stress 0.1056904 
    ## Run 410 stress 0.1076307 
    ## Run 411 stress 0.1063195 
    ## Run 412 stress 0.09228504 
    ## Run 413 stress 0.1092562 
    ## Run 414 stress 0.1086565 
    ## Run 415 stress 0.09044608 
    ## Run 416 stress 0.09592145 
    ## Run 417 stress 0.09503478 
    ## Run 418 stress 0.1064188 
    ## Run 419 stress 0.1079014 
    ## Run 420 stress 0.1065346 
    ## Run 421 stress 0.1074432 
    ## Run 422 stress 0.08926105 
    ## ... Procrustes: rmse 0.0004096484  max resid 0.001286637 
    ## ... Similar to previous best
    ## Run 423 stress 0.08926064 
    ## ... New best solution
    ## ... Procrustes: rmse 9.392839e-05  max resid 0.0002554057 
    ## ... Similar to previous best
    ## Run 424 stress 0.08946654 
    ## ... Procrustes: rmse 0.0332423  max resid 0.1163949 
    ## Run 425 stress 0.09503424 
    ## Run 426 stress 0.08938542 
    ## ... Procrustes: rmse 0.01189109  max resid 0.03971249 
    ## Run 427 stress 0.1092091 
    ## Run 428 stress 0.1067895 
    ## Run 429 stress 0.08926065 
    ## ... Procrustes: rmse 8.117206e-05  max resid 0.0002281946 
    ## ... Similar to previous best
    ## Run 430 stress 0.08946656 
    ## ... Procrustes: rmse 0.03324487  max resid 0.1164069 
    ## Run 431 stress 0.08938557 
    ## ... Procrustes: rmse 0.01192246  max resid 0.03967838 
    ## Run 432 stress 0.09021166 
    ## Run 433 stress 0.1089907 
    ## Run 434 stress 0.1064187 
    ## Run 435 stress 0.09503461 
    ## Run 436 stress 0.1101451 
    ## Run 437 stress 0.08946664 
    ## ... Procrustes: rmse 0.03318672  max resid 0.1162844 
    ## Run 438 stress 0.1108317 
    ## Run 439 stress 0.1064428 
    ## Run 440 stress 0.1071322 
    ## Run 441 stress 0.1074432 
    ## Run 442 stress 0.09021166 
    ## Run 443 stress 0.1063193 
    ## Run 444 stress 0.08946328 
    ## ... Procrustes: rmse 0.03751448  max resid 0.1177598 
    ## Run 445 stress 0.09109115 
    ## Run 446 stress 0.09178287 
    ## Run 447 stress 0.1087867 
    ## Run 448 stress 0.08946331 
    ## ... Procrustes: rmse 0.03750405  max resid 0.117745 
    ## Run 449 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595774  max resid 0.1182925 
    ## Run 450 stress 0.0894633 
    ## ... Procrustes: rmse 0.03753571  max resid 0.1177921 
    ## Run 451 stress 0.1060433 
    ## Run 452 stress 0.08926113 
    ## ... Procrustes: rmse 0.000418191  max resid 0.001344166 
    ## ... Similar to previous best
    ## Run 453 stress 0.1086563 
    ## Run 454 stress 0.1075421 
    ## Run 455 stress 0.0894634 
    ## ... Procrustes: rmse 0.03748089  max resid 0.1177082 
    ## Run 456 stress 0.09109113 
    ## Run 457 stress 0.08926071 
    ## ... Procrustes: rmse 0.000202335  max resid 0.0005292203 
    ## ... Similar to previous best
    ## Run 458 stress 0.08946652 
    ## ... Procrustes: rmse 0.03323664  max resid 0.1163757 
    ## Run 459 stress 0.09109108 
    ## Run 460 stress 0.08926064 
    ## ... Procrustes: rmse 0.0001171924  max resid 0.0002780896 
    ## ... Similar to previous best
    ## Run 461 stress 0.09130086 
    ## Run 462 stress 0.08938984 
    ## ... Procrustes: rmse 0.03590555  max resid 0.118207 
    ## Run 463 stress 0.09503458 
    ## Run 464 stress 0.1056896 
    ## Run 465 stress 0.09262347 
    ## Run 466 stress 0.1052646 
    ## Run 467 stress 0.08938548 
    ## ... Procrustes: rmse 0.01178841  max resid 0.03971173 
    ## Run 468 stress 0.1095633 
    ## Run 469 stress 0.08926086 
    ## ... Procrustes: rmse 0.0002828505  max resid 0.0008687232 
    ## ... Similar to previous best
    ## Run 470 stress 0.1067899 
    ## Run 471 stress 0.08938538 
    ## ... Procrustes: rmse 0.01185578  max resid 0.03974155 
    ## Run 472 stress 0.08926092 
    ## ... Procrustes: rmse 0.0003465795  max resid 0.0009878832 
    ## ... Similar to previous best
    ## Run 473 stress 0.09087345 
    ## Run 474 stress 0.08946651 
    ## ... Procrustes: rmse 0.03323166  max resid 0.1163615 
    ## Run 475 stress 0.10613 
    ## Run 476 stress 0.1087865 
    ## Run 477 stress 0.1086561 
    ## Run 478 stress 0.1052647 
    ## Run 479 stress 0.1052647 
    ## Run 480 stress 0.09087339 
    ## Run 481 stress 0.0893897 
    ## ... Procrustes: rmse 0.03592367  max resid 0.1182379 
    ## Run 482 stress 0.1052649 
    ## Run 483 stress 0.1074433 
    ## Run 484 stress 0.1067585 
    ## Run 485 stress 0.09592144 
    ## Run 486 stress 0.09087342 
    ## Run 487 stress 0.1087864 
    ## Run 488 stress 0.08938541 
    ## ... Procrustes: rmse 0.0118709  max resid 0.03976677 
    ## Run 489 stress 0.0894665 
    ## ... Procrustes: rmse 0.03322183  max resid 0.1163403 
    ## Run 490 stress 0.09087354 
    ## Run 491 stress 0.1067369 
    ## Run 492 stress 0.0893854 
    ## ... Procrustes: rmse 0.01185826  max resid 0.03974801 
    ## Run 493 stress 0.08938554 
    ## ... Procrustes: rmse 0.01178555  max resid 0.03963833 
    ## Run 494 stress 0.08938548 
    ## ... Procrustes: rmse 0.01184927  max resid 0.03964886 
    ## Run 495 stress 0.09178303 
    ## Run 496 stress 0.1076303 
    ## Run 497 stress 0.1071321 
    ## Run 498 stress 0.09018712 
    ## Run 499 stress 0.09130086 
    ## Run 500 stress 0.08938967 
    ## ... Procrustes: rmse 0.03598158  max resid 0.1183326 
    ## *** Best solution repeated 7 times

``` r
# Mixed and stratified lakes
SD_beta_geo_MS_NMDS <- metaMDS(SD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09407972 
    ## Run 2 stress 0.09145326 
    ## Run 3 stress 0.09286098 
    ## Run 4 stress 0.09030395 
    ## Run 5 stress 0.09159091 
    ## Run 6 stress 0.09030393 
    ## Run 7 stress 0.09374329 
    ## Run 8 stress 0.08973867 
    ## Run 9 stress 0.09030397 
    ## Run 10 stress 0.09381844 
    ## Run 11 stress 0.08440261 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01443498  max resid 0.04304662 
    ## Run 12 stress 0.08440256 
    ## ... New best solution
    ## ... Procrustes: rmse 5.467607e-05  max resid 0.0001130468 
    ## ... Similar to previous best
    ## Run 13 stress 0.09308949 
    ## Run 14 stress 0.08503496 
    ## Run 15 stress 0.08503602 
    ## Run 16 stress 0.09286086 
    ## Run 17 stress 0.09407988 
    ## Run 18 stress 0.2465516 
    ## Run 19 stress 0.09286083 
    ## Run 20 stress 0.1006369 
    ## Run 21 stress 0.09590556 
    ## Run 22 stress 0.08773493 
    ## Run 23 stress 0.09030397 
    ## Run 24 stress 0.08440255 
    ## ... New best solution
    ## ... Procrustes: rmse 2.584927e-05  max resid 4.748479e-05 
    ## ... Similar to previous best
    ## Run 25 stress 0.08773467 
    ## Run 26 stress 0.09159084 
    ## Run 27 stress 0.08773479 
    ## Run 28 stress 0.09286088 
    ## Run 29 stress 0.08773484 
    ## Run 30 stress 0.08773474 
    ## Run 31 stress 0.0914533 
    ## Run 32 stress 0.08503687 
    ## Run 33 stress 0.09407986 
    ## Run 34 stress 0.08773466 
    ## Run 35 stress 0.08440263 
    ## ... Procrustes: rmse 0.0003288718  max resid 0.0006488735 
    ## ... Similar to previous best
    ## Run 36 stress 0.09445468 
    ## Run 37 stress 0.09268389 
    ## Run 38 stress 0.08503466 
    ## Run 39 stress 0.08503456 
    ## Run 40 stress 0.08973872 
    ## Run 41 stress 0.09159095 
    ## Run 42 stress 0.1026218 
    ## Run 43 stress 0.09308966 
    ## Run 44 stress 0.09286083 
    ## Run 45 stress 0.09168938 
    ## Run 46 stress 0.09168936 
    ## Run 47 stress 0.09030407 
    ## Run 48 stress 0.09407976 
    ## Run 49 stress 0.09535493 
    ## Run 50 stress 0.09168929 
    ## Run 51 stress 0.2574702 
    ## Run 52 stress 0.08503491 
    ## Run 53 stress 0.08503457 
    ## Run 54 stress 0.1038047 
    ## Run 55 stress 0.09145337 
    ## Run 56 stress 0.09159087 
    ## Run 57 stress 0.08440255 
    ## ... Procrustes: rmse 1.367014e-05  max resid 3.000905e-05 
    ## ... Similar to previous best
    ## Run 58 stress 0.09030415 
    ## Run 59 stress 0.09400525 
    ## Run 60 stress 0.09308968 
    ## Run 61 stress 0.08773484 
    ## Run 62 stress 0.1005995 
    ## Run 63 stress 0.09337294 
    ## Run 64 stress 0.09286084 
    ## Run 65 stress 0.09407136 
    ## Run 66 stress 0.08440263 
    ## ... Procrustes: rmse 8.555553e-05  max resid 0.000145968 
    ## ... Similar to previous best
    ## Run 67 stress 0.313882 
    ## Run 68 stress 0.08973863 
    ## Run 69 stress 0.09760825 
    ## Run 70 stress 0.08973863 
    ## Run 71 stress 0.09168935 
    ## Run 72 stress 0.08440258 
    ## ... Procrustes: rmse 3.530791e-05  max resid 6.098926e-05 
    ## ... Similar to previous best
    ## Run 73 stress 0.09145321 
    ## Run 74 stress 0.09407976 
    ## Run 75 stress 0.08503704 
    ## Run 76 stress 0.09407964 
    ## Run 77 stress 0.1038037 
    ## Run 78 stress 0.0850362 
    ## Run 79 stress 0.09337239 
    ## Run 80 stress 0.08773466 
    ## Run 81 stress 0.08440257 
    ## ... Procrustes: rmse 3.239912e-05  max resid 5.505026e-05 
    ## ... Similar to previous best
    ## Run 82 stress 0.09407965 
    ## Run 83 stress 0.09145325 
    ## Run 84 stress 0.09407407 
    ## Run 85 stress 0.08773473 
    ## Run 86 stress 0.09445519 
    ## Run 87 stress 0.09268336 
    ## Run 88 stress 0.08503525 
    ## Run 89 stress 0.09416138 
    ## Run 90 stress 0.09268395 
    ## Run 91 stress 0.09168936 
    ## Run 92 stress 0.09268352 
    ## Run 93 stress 0.09159086 
    ## Run 94 stress 0.09539192 
    ## Run 95 stress 0.09407976 
    ## Run 96 stress 0.2923523 
    ## Run 97 stress 0.09465906 
    ## Run 98 stress 0.3312755 
    ## Run 99 stress 0.08503484 
    ## Run 100 stress 0.09286092 
    ## Run 101 stress 0.0903041 
    ## Run 102 stress 0.08440253 
    ## ... New best solution
    ## ... Procrustes: rmse 5.286644e-05  max resid 0.000100027 
    ## ... Similar to previous best
    ## Run 103 stress 0.09465906 
    ## Run 104 stress 0.09407971 
    ## Run 105 stress 0.09465905 
    ## Run 106 stress 0.09447313 
    ## Run 107 stress 0.09145323 
    ## Run 108 stress 0.09159096 
    ## Run 109 stress 0.09030396 
    ## Run 110 stress 0.0916894 
    ## Run 111 stress 0.09445506 
    ## Run 112 stress 0.09268359 
    ## Run 113 stress 0.08773476 
    ## Run 114 stress 0.09168956 
    ## Run 115 stress 0.08440258 
    ## ... Procrustes: rmse 9.25953e-05  max resid 0.0001906279 
    ## ... Similar to previous best
    ## Run 116 stress 0.08973863 
    ## Run 117 stress 0.08440268 
    ## ... Procrustes: rmse 0.000178765  max resid 0.0003087487 
    ## ... Similar to previous best
    ## Run 118 stress 0.09416142 
    ## Run 119 stress 0.0915909 
    ## Run 120 stress 0.09268342 
    ## Run 121 stress 0.09381846 
    ## Run 122 stress 0.09407112 
    ## Run 123 stress 0.09416133 
    ## Run 124 stress 0.08973863 
    ## Run 125 stress 0.09286085 
    ## Run 126 stress 0.09337251 
    ## Run 127 stress 0.09407989 
    ## Run 128 stress 0.09308951 
    ## Run 129 stress 0.0941777 
    ## Run 130 stress 0.09760822 
    ## Run 131 stress 0.09969975 
    ## Run 132 stress 0.09030409 
    ## Run 133 stress 0.09168942 
    ## Run 134 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 4.294909e-05  max resid 8.971112e-05 
    ## ... Similar to previous best
    ## Run 135 stress 0.08503505 
    ## Run 136 stress 0.09145321 
    ## Run 137 stress 0.09464482 
    ## Run 138 stress 0.08973865 
    ## Run 139 stress 0.08503659 
    ## Run 140 stress 0.08773467 
    ## Run 141 stress 0.09416137 
    ## Run 142 stress 0.09145324 
    ## Run 143 stress 0.09337271 
    ## Run 144 stress 0.09030405 
    ## Run 145 stress 0.09030394 
    ## Run 146 stress 0.09969965 
    ## Run 147 stress 0.09337281 
    ## Run 148 stress 0.09374279 
    ## Run 149 stress 0.09417816 
    ## Run 150 stress 0.09407991 
    ## Run 151 stress 0.09381847 
    ## Run 152 stress 0.0996999 
    ## Run 153 stress 0.08973867 
    ## Run 154 stress 0.1052849 
    ## Run 155 stress 0.0903041 
    ## Run 156 stress 0.09145322 
    ## Run 157 stress 0.09286084 
    ## Run 158 stress 0.09030411 
    ## Run 159 stress 0.09030398 
    ## Run 160 stress 0.09268327 
    ## Run 161 stress 0.09969985 
    ## Run 162 stress 0.09145331 
    ## Run 163 stress 0.0897388 
    ## Run 164 stress 0.09535418 
    ## Run 165 stress 0.1038038 
    ## Run 166 stress 0.09030404 
    ## Run 167 stress 0.09168961 
    ## Run 168 stress 0.08503523 
    ## Run 169 stress 0.09404402 
    ## Run 170 stress 0.09030401 
    ## Run 171 stress 0.08773479 
    ## Run 172 stress 0.08440251 
    ## ... Procrustes: rmse 2.098908e-05  max resid 4.090056e-05 
    ## ... Similar to previous best
    ## Run 173 stress 0.09308952 
    ## Run 174 stress 0.09308946 
    ## Run 175 stress 0.09969969 
    ## Run 176 stress 0.09407108 
    ## Run 177 stress 0.3126486 
    ## Run 178 stress 0.0940797 
    ## Run 179 stress 0.0933729 
    ## Run 180 stress 0.09268355 
    ## Run 181 stress 0.09145324 
    ## Run 182 stress 0.09416151 
    ## Run 183 stress 0.09268339 
    ## Run 184 stress 0.0938058 
    ## Run 185 stress 0.09321772 
    ## Run 186 stress 0.09535442 
    ## Run 187 stress 0.08773471 
    ## Run 188 stress 0.08973883 
    ## Run 189 stress 0.09030394 
    ## Run 190 stress 0.08973875 
    ## Run 191 stress 0.08973866 
    ## Run 192 stress 0.09268391 
    ## Run 193 stress 0.1038034 
    ## Run 194 stress 0.09030396 
    ## Run 195 stress 0.09286088 
    ## Run 196 stress 0.0996998 
    ## Run 197 stress 0.09159086 
    ## Run 198 stress 0.09268393 
    ## Run 199 stress 0.09400516 
    ## Run 200 stress 0.091591 
    ## Run 201 stress 0.09465905 
    ## Run 202 stress 0.0850353 
    ## Run 203 stress 0.09374375 
    ## Run 204 stress 0.09268413 
    ## Run 205 stress 0.08503475 
    ## Run 206 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 4.052723e-05  max resid 0.0001011257 
    ## ... Similar to previous best
    ## Run 207 stress 0.08503458 
    ## Run 208 stress 0.09407977 
    ## Run 209 stress 0.08503502 
    ## Run 210 stress 0.0916895 
    ## Run 211 stress 0.09407331 
    ## Run 212 stress 0.09610756 
    ## Run 213 stress 0.09286092 
    ## Run 214 stress 0.09337276 
    ## Run 215 stress 0.09145329 
    ## Run 216 stress 0.09535407 
    ## Run 217 stress 0.08773466 
    ## Run 218 stress 0.2456114 
    ## Run 219 stress 0.09030397 
    ## Run 220 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001384769  max resid 0.0002628066 
    ## ... Similar to previous best
    ## Run 221 stress 0.1038039 
    ## Run 222 stress 0.0916895 
    ## Run 223 stress 0.09159084 
    ## Run 224 stress 0.09337282 
    ## Run 225 stress 0.09030394 
    ## Run 226 stress 0.09381849 
    ## Run 227 stress 0.091591 
    ## Run 228 stress 0.09380582 
    ## Run 229 stress 0.09159084 
    ## Run 230 stress 0.09400543 
    ## Run 231 stress 0.09337248 
    ## Run 232 stress 0.09969964 
    ## Run 233 stress 0.09286087 
    ## Run 234 stress 0.09400515 
    ## Run 235 stress 0.09535467 
    ## Run 236 stress 0.09145324 
    ## Run 237 stress 0.09590814 
    ## Run 238 stress 0.09337243 
    ## Run 239 stress 0.09535405 
    ## Run 240 stress 0.09308959 
    ## Run 241 stress 0.0953553 
    ## Run 242 stress 0.09030406 
    ## Run 243 stress 0.09465904 
    ## Run 244 stress 0.09969968 
    ## Run 245 stress 0.09535496 
    ## Run 246 stress 0.09308982 
    ## Run 247 stress 0.09760813 
    ## Run 248 stress 0.09337229 
    ## Run 249 stress 0.09407985 
    ## Run 250 stress 0.09969972 
    ## Run 251 stress 0.09969966 
    ## Run 252 stress 0.09030405 
    ## Run 253 stress 0.09030413 
    ## Run 254 stress 0.09268388 
    ## Run 255 stress 0.0946445 
    ## Run 256 stress 0.0976083 
    ## Run 257 stress 0.08503461 
    ## Run 258 stress 0.09145323 
    ## Run 259 stress 0.09535482 
    ## Run 260 stress 0.09268341 
    ## Run 261 stress 0.09030411 
    ## Run 262 stress 0.09321801 
    ## Run 263 stress 0.09407996 
    ## Run 264 stress 0.1053421 
    ## Run 265 stress 0.09168946 
    ## Run 266 stress 0.09286083 
    ## Run 267 stress 0.08773471 
    ## Run 268 stress 0.09407996 
    ## Run 269 stress 0.08773469 
    ## Run 270 stress 0.09145331 
    ## Run 271 stress 0.09030397 
    ## Run 272 stress 0.08503499 
    ## Run 273 stress 0.08440252 
    ## ... Procrustes: rmse 6.217341e-05  max resid 0.000117276 
    ## ... Similar to previous best
    ## Run 274 stress 0.09309007 
    ## Run 275 stress 0.09407128 
    ## Run 276 stress 0.09159094 
    ## Run 277 stress 0.09374185 
    ## Run 278 stress 0.09268386 
    ## Run 279 stress 0.09308974 
    ## Run 280 stress 0.1054536 
    ## Run 281 stress 0.09407983 
    ## Run 282 stress 0.08773479 
    ## Run 283 stress 0.3254873 
    ## Run 284 stress 0.09465904 
    ## Run 285 stress 0.09535433 
    ## Run 286 stress 0.09145321 
    ## Run 287 stress 0.09374308 
    ## Run 288 stress 0.09969972 
    ## Run 289 stress 0.09539218 
    ## Run 290 stress 0.09030399 
    ## Run 291 stress 0.09308979 
    ## Run 292 stress 0.09464422 
    ## Run 293 stress 0.1038032 
    ## Run 294 stress 0.08503467 
    ## Run 295 stress 0.09030394 
    ## Run 296 stress 0.09286083 
    ## Run 297 stress 0.09145327 
    ## Run 298 stress 0.09465911 
    ## Run 299 stress 0.09168932 
    ## Run 300 stress 0.09308962 
    ## Run 301 stress 0.08773481 
    ## Run 302 stress 0.09030401 
    ## Run 303 stress 0.09416158 
    ## Run 304 stress 0.08503511 
    ## Run 305 stress 0.09407969 
    ## Run 306 stress 0.09374331 
    ## Run 307 stress 0.09464451 
    ## Run 308 stress 0.09145333 
    ## Run 309 stress 0.09168939 
    ## Run 310 stress 0.09445474 
    ## Run 311 stress 0.09030396 
    ## Run 312 stress 0.08440252 
    ## ... Procrustes: rmse 5.122528e-05  max resid 9.884495e-05 
    ## ... Similar to previous best
    ## Run 313 stress 0.09030399 
    ## Run 314 stress 0.09407988 
    ## Run 315 stress 0.08503502 
    ## Run 316 stress 0.09465909 
    ## Run 317 stress 0.08503491 
    ## Run 318 stress 0.09159088 
    ## Run 319 stress 0.08440264 
    ## ... Procrustes: rmse 0.000163271  max resid 0.0003107752 
    ## ... Similar to previous best
    ## Run 320 stress 0.08440275 
    ## ... Procrustes: rmse 0.0002879083  max resid 0.0005887757 
    ## ... Similar to previous best
    ## Run 321 stress 0.09760815 
    ## Run 322 stress 0.09308949 
    ## Run 323 stress 0.08503486 
    ## Run 324 stress 0.09374267 
    ## Run 325 stress 0.08773469 
    ## Run 326 stress 0.1038088 
    ## Run 327 stress 0.09030409 
    ## Run 328 stress 0.09539193 
    ## Run 329 stress 0.08973887 
    ## Run 330 stress 0.09168929 
    ## Run 331 stress 0.09380596 
    ## Run 332 stress 0.09447296 
    ## Run 333 stress 0.09159091 
    ## Run 334 stress 0.09969995 
    ## Run 335 stress 0.09337257 
    ## Run 336 stress 0.09159099 
    ## Run 337 stress 0.09030398 
    ## Run 338 stress 0.09268384 
    ## Run 339 stress 0.09465911 
    ## Run 340 stress 0.0877349 
    ## Run 341 stress 0.08773467 
    ## Run 342 stress 0.08503486 
    ## Run 343 stress 0.08973882 
    ## Run 344 stress 0.09030396 
    ## Run 345 stress 0.09286103 
    ## Run 346 stress 0.09286087 
    ## Run 347 stress 0.1072862 
    ## Run 348 stress 0.09407355 
    ## Run 349 stress 0.09145343 
    ## Run 350 stress 0.08503586 
    ## Run 351 stress 0.09407967 
    ## Run 352 stress 0.09969976 
    ## Run 353 stress 0.08773485 
    ## Run 354 stress 0.08773469 
    ## Run 355 stress 0.1019597 
    ## Run 356 stress 0.09380585 
    ## Run 357 stress 0.09400521 
    ## Run 358 stress 0.09535496 
    ## Run 359 stress 0.09408013 
    ## Run 360 stress 0.09337268 
    ## Run 361 stress 0.08503458 
    ## Run 362 stress 0.09268354 
    ## Run 363 stress 0.09168945 
    ## Run 364 stress 0.09407975 
    ## Run 365 stress 0.09977536 
    ## Run 366 stress 0.08503494 
    ## Run 367 stress 0.09159084 
    ## Run 368 stress 0.09030403 
    ## Run 369 stress 0.09721199 
    ## Run 370 stress 0.08973871 
    ## Run 371 stress 0.09465913 
    ## Run 372 stress 0.09408008 
    ## Run 373 stress 0.09590923 
    ## Run 374 stress 0.09337269 
    ## Run 375 stress 0.0930898 
    ## Run 376 stress 0.09030393 
    ## Run 377 stress 0.09407975 
    ## Run 378 stress 0.09159088 
    ## Run 379 stress 0.09380593 
    ## Run 380 stress 0.09337232 
    ## Run 381 stress 0.09590578 
    ## Run 382 stress 0.1038042 
    ## Run 383 stress 0.09286086 
    ## Run 384 stress 0.08503467 
    ## Run 385 stress 0.08503465 
    ## Run 386 stress 0.0915909 
    ## Run 387 stress 0.08973883 
    ## Run 388 stress 0.09168936 
    ## Run 389 stress 0.08773466 
    ## Run 390 stress 0.09417786 
    ## Run 391 stress 0.09268327 
    ## Run 392 stress 0.09308955 
    ## Run 393 stress 0.08503512 
    ## Run 394 stress 0.09030394 
    ## Run 395 stress 0.09145327 
    ## Run 396 stress 0.0997 
    ## Run 397 stress 0.09407984 
    ## Run 398 stress 0.0903041 
    ## Run 399 stress 0.09286083 
    ## Run 400 stress 0.09268376 
    ## Run 401 stress 0.09407992 
    ## Run 402 stress 0.09030411 
    ## Run 403 stress 0.0996999 
    ## Run 404 stress 0.0938185 
    ## Run 405 stress 0.08503494 
    ## Run 406 stress 0.1038044 
    ## Run 407 stress 0.09286083 
    ## Run 408 stress 0.08973873 
    ## Run 409 stress 0.1026217 
    ## Run 410 stress 0.09268324 
    ## Run 411 stress 0.09407992 
    ## Run 412 stress 0.09168944 
    ## Run 413 stress 0.0938058 
    ## Run 414 stress 0.09407985 
    ## Run 415 stress 0.09416141 
    ## Run 416 stress 0.09447313 
    ## Run 417 stress 0.08503486 
    ## Run 418 stress 0.09337268 
    ## Run 419 stress 0.09760826 
    ## Run 420 stress 0.08503471 
    ## Run 421 stress 0.09168929 
    ## Run 422 stress 0.09407965 
    ## Run 423 stress 0.09337245 
    ## Run 424 stress 0.09337235 
    ## Run 425 stress 0.0940798 
    ## Run 426 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001129702  max resid 0.0002734122 
    ## ... Similar to previous best
    ## Run 427 stress 0.09374337 
    ## Run 428 stress 0.09030395 
    ## Run 429 stress 0.08973874 
    ## Run 430 stress 0.09465912 
    ## Run 431 stress 0.09969975 
    ## Run 432 stress 0.1026216 
    ## Run 433 stress 0.09308985 
    ## Run 434 stress 0.09159086 
    ## Run 435 stress 0.09407981 
    ## Run 436 stress 0.08773478 
    ## Run 437 stress 0.09408017 
    ## Run 438 stress 0.0915909 
    ## Run 439 stress 0.09408005 
    ## Run 440 stress 0.08973865 
    ## Run 441 stress 0.08973872 
    ## Run 442 stress 0.09286084 
    ## Run 443 stress 0.09337277 
    ## Run 444 stress 0.09539215 
    ## Run 445 stress 0.1038045 
    ## Run 446 stress 0.09337291 
    ## Run 447 stress 0.08973873 
    ## Run 448 stress 0.09030395 
    ## Run 449 stress 0.09168955 
    ## Run 450 stress 0.09286085 
    ## Run 451 stress 0.09145323 
    ## Run 452 stress 0.09268376 
    ## Run 453 stress 0.3399088 
    ## Run 454 stress 0.1038047 
    ## Run 455 stress 0.09159084 
    ## Run 456 stress 0.09400524 
    ## Run 457 stress 0.09268412 
    ## Run 458 stress 0.0940437 
    ## Run 459 stress 0.08773468 
    ## Run 460 stress 0.09535432 
    ## Run 461 stress 0.09030401 
    ## Run 462 stress 0.09030393 
    ## Run 463 stress 0.1026217 
    ## Run 464 stress 0.09416141 
    ## Run 465 stress 0.09407992 
    ## Run 466 stress 0.2455402 
    ## Run 467 stress 0.08440252 
    ## ... Procrustes: rmse 5.859726e-05  max resid 0.0001409224 
    ## ... Similar to previous best
    ## Run 468 stress 0.09465904 
    ## Run 469 stress 0.09145321 
    ## Run 470 stress 0.09030409 
    ## Run 471 stress 0.08773468 
    ## Run 472 stress 0.08773467 
    ## Run 473 stress 0.100461 
    ## Run 474 stress 0.08503481 
    ## Run 475 stress 0.08773467 
    ## Run 476 stress 0.09721197 
    ## Run 477 stress 0.09337273 
    ## Run 478 stress 0.09286092 
    ## Run 479 stress 0.09407997 
    ## Run 480 stress 0.09030393 
    ## Run 481 stress 0.09337293 
    ## Run 482 stress 0.08503504 
    ## Run 483 stress 0.09030406 
    ## Run 484 stress 0.09145326 
    ## Run 485 stress 0.09407999 
    ## Run 486 stress 0.08503627 
    ## Run 487 stress 0.08503477 
    ## Run 488 stress 0.1053431 
    ## Run 489 stress 0.09407979 
    ## Run 490 stress 0.08440266 
    ## ... Procrustes: rmse 0.0002372421  max resid 0.0004361442 
    ## ... Similar to previous best
    ## Run 491 stress 0.09464483 
    ## Run 492 stress 0.09337285 
    ## Run 493 stress 0.091591 
    ## Run 494 stress 0.09407965 
    ## Run 495 stress 0.0996998 
    ## Run 496 stress 0.1052847 
    ## Run 497 stress 0.08503511 
    ## Run 498 stress 0.09374302 
    ## Run 499 stress 0.08973882 
    ## Run 500 stress 0.09416131 
    ## *** Best solution repeated 9 times

``` r
# Ocean sites and mixed lakes
SD_beta_geo_OM_NMDS <- metaMDS(SD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07365783 
    ## Run 1 stress 0.07732944 
    ## Run 2 stress 0.07365783 
    ## ... Procrustes: rmse 1.240521e-05  max resid 2.51165e-05 
    ## ... Similar to previous best
    ## Run 3 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745458  max resid 0.05109044 
    ## Run 4 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744303  max resid 0.05103098 
    ## Run 5 stress 0.08000468 
    ## Run 6 stress 0.07365787 
    ## ... Procrustes: rmse 6.578808e-05  max resid 0.0001500132 
    ## ... Similar to previous best
    ## Run 7 stress 0.07629237 
    ## Run 8 stress 0.07629245 
    ## Run 9 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001329751  max resid 0.0003147872 
    ## ... Similar to previous best
    ## Run 10 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744594  max resid 0.0510469 
    ## Run 11 stress 0.07629235 
    ## Run 12 stress 0.07629238 
    ## Run 13 stress 0.07365787 
    ## ... Procrustes: rmse 9.727162e-05  max resid 0.0002283726 
    ## ... Similar to previous best
    ## Run 14 stress 0.07629248 
    ## Run 15 stress 0.0737823 
    ## ... Procrustes: rmse 0.0174518  max resid 0.0510976 
    ## Run 16 stress 0.08000467 
    ## Run 17 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 7.121399e-06  max resid 1.412816e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174457  max resid 0.05105682 
    ## Run 19 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174464  max resid 0.05103926 
    ## Run 20 stress 0.08000468 
    ## Run 21 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743566  max resid 0.0509937 
    ## Run 22 stress 0.08000469 
    ## Run 23 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001450131  max resid 0.0003378915 
    ## ... Similar to previous best
    ## Run 24 stress 0.08000467 
    ## Run 25 stress 0.07365786 
    ## ... Procrustes: rmse 9.595532e-05  max resid 0.0002250618 
    ## ... Similar to previous best
    ## Run 26 stress 0.07365783 
    ## ... Procrustes: rmse 1.224483e-05  max resid 2.804818e-05 
    ## ... Similar to previous best
    ## Run 27 stress 0.08000473 
    ## Run 28 stress 0.08000467 
    ## Run 29 stress 0.07629241 
    ## Run 30 stress 0.07365785 
    ## ... Procrustes: rmse 6.995468e-05  max resid 0.0001653048 
    ## ... Similar to previous best
    ## Run 31 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744555  max resid 0.05103974 
    ## Run 32 stress 0.08000468 
    ## Run 33 stress 0.07629242 
    ## Run 34 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744666  max resid 0.05105907 
    ## Run 35 stress 0.08000472 
    ## Run 36 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745613  max resid 0.0511059 
    ## Run 37 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001220725  max resid 0.0002890259 
    ## ... Similar to previous best
    ## Run 38 stress 0.07732928 
    ## Run 39 stress 0.07629236 
    ## Run 40 stress 0.07732932 
    ## Run 41 stress 0.07365783 
    ## ... Procrustes: rmse 1.790435e-05  max resid 4.132132e-05 
    ## ... Similar to previous best
    ## Run 42 stress 0.07629243 
    ## Run 43 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744807  max resid 0.05106101 
    ## Run 44 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744624  max resid 0.05104374 
    ## Run 45 stress 0.08000468 
    ## Run 46 stress 0.08000471 
    ## Run 47 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744972  max resid 0.05106431 
    ## Run 48 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001581068  max resid 0.0003716465 
    ## ... Similar to previous best
    ## Run 49 stress 0.07365784 
    ## ... Procrustes: rmse 4.333736e-05  max resid 0.0001021249 
    ## ... Similar to previous best
    ## Run 50 stress 0.07629235 
    ## Run 51 stress 0.07365787 
    ## ... Procrustes: rmse 7.333256e-05  max resid 0.0001603162 
    ## ... Similar to previous best
    ## Run 52 stress 0.07365784 
    ## ... Procrustes: rmse 4.084194e-05  max resid 9.483533e-05 
    ## ... Similar to previous best
    ## Run 53 stress 0.07629236 
    ## Run 54 stress 0.07629234 
    ## Run 55 stress 0.08233765 
    ## Run 56 stress 0.0762925 
    ## Run 57 stress 0.07629246 
    ## Run 58 stress 0.08288042 
    ## Run 59 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744359  max resid 0.05101917 
    ## Run 60 stress 0.07629235 
    ## Run 61 stress 0.08288058 
    ## Run 62 stress 0.08000468 
    ## Run 63 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744474  max resid 0.05103212 
    ## Run 64 stress 0.07365787 
    ## ... Procrustes: rmse 9.909029e-05  max resid 0.0002308364 
    ## ... Similar to previous best
    ## Run 65 stress 0.07629236 
    ## Run 66 stress 0.07629233 
    ## Run 67 stress 0.0800047 
    ## Run 68 stress 0.07629237 
    ## Run 69 stress 0.07629238 
    ## Run 70 stress 0.0800047 
    ## Run 71 stress 0.08000467 
    ## Run 72 stress 0.07629236 
    ## Run 73 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001348969  max resid 0.0003167145 
    ## ... Similar to previous best
    ## Run 74 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744525  max resid 0.05103237 
    ## Run 75 stress 0.08000467 
    ## Run 76 stress 0.08000468 
    ## Run 77 stress 0.2676491 
    ## Run 78 stress 0.07365785 
    ## ... Procrustes: rmse 7.140924e-05  max resid 0.0001677925 
    ## ... Similar to previous best
    ## Run 79 stress 0.07629238 
    ## Run 80 stress 0.07629243 
    ## Run 81 stress 0.08000476 
    ## Run 82 stress 0.07365784 
    ## ... Procrustes: rmse 5.989628e-05  max resid 0.0001399533 
    ## ... Similar to previous best
    ## Run 83 stress 0.08000476 
    ## Run 84 stress 0.07365783 
    ## ... Procrustes: rmse 2.77812e-05  max resid 6.483295e-05 
    ## ... Similar to previous best
    ## Run 85 stress 0.07629234 
    ## Run 86 stress 0.07629241 
    ## Run 87 stress 0.07732927 
    ## Run 88 stress 0.07629234 
    ## Run 89 stress 0.08000476 
    ## Run 90 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001420796  max resid 0.000332789 
    ## ... Similar to previous best
    ## Run 91 stress 0.08000469 
    ## Run 92 stress 0.07629239 
    ## Run 93 stress 0.07732929 
    ## Run 94 stress 0.3578479 
    ## Run 95 stress 0.07732925 
    ## Run 96 stress 0.07365783 
    ## ... Procrustes: rmse 5.285976e-06  max resid 1.124195e-05 
    ## ... Similar to previous best
    ## Run 97 stress 0.08000487 
    ## Run 98 stress 0.0773294 
    ## Run 99 stress 0.07365784 
    ## ... Procrustes: rmse 3.172703e-05  max resid 7.182573e-05 
    ## ... Similar to previous best
    ## Run 100 stress 0.07629235 
    ## Run 101 stress 0.08000472 
    ## Run 102 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744498  max resid 0.05106098 
    ## Run 103 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001459715  max resid 0.0003440823 
    ## ... Similar to previous best
    ## Run 104 stress 0.07629242 
    ## Run 105 stress 0.0762924 
    ## Run 106 stress 0.07732932 
    ## Run 107 stress 0.07365787 
    ## ... Procrustes: rmse 9.42362e-05  max resid 0.0002252607 
    ## ... Similar to previous best
    ## Run 108 stress 0.07365783 
    ## ... Procrustes: rmse 2.224504e-05  max resid 5.269707e-05 
    ## ... Similar to previous best
    ## Run 109 stress 0.07365785 
    ## ... Procrustes: rmse 7.337589e-05  max resid 0.0001724558 
    ## ... Similar to previous best
    ## Run 110 stress 0.0737823 
    ## ... Procrustes: rmse 0.01746409  max resid 0.05109114 
    ## Run 111 stress 0.07629235 
    ## Run 112 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745004  max resid 0.05107574 
    ## Run 113 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744702  max resid 0.05105475 
    ## Run 114 stress 0.07365786 
    ## ... Procrustes: rmse 8.406197e-05  max resid 0.0001992503 
    ## ... Similar to previous best
    ## Run 115 stress 0.08000472 
    ## Run 116 stress 0.08000467 
    ## Run 117 stress 0.07629239 
    ## Run 118 stress 0.07365783 
    ## ... Procrustes: rmse 1.865643e-05  max resid 3.720001e-05 
    ## ... Similar to previous best
    ## Run 119 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744681  max resid 0.05104669 
    ## Run 120 stress 0.08000468 
    ## Run 121 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744915  max resid 0.05104866 
    ## Run 122 stress 0.07629234 
    ## Run 123 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745476  max resid 0.05106025 
    ## Run 124 stress 0.08000467 
    ## Run 125 stress 0.08000468 
    ## Run 126 stress 0.07629239 
    ## Run 127 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743725  max resid 0.05100704 
    ## Run 128 stress 0.08000471 
    ## Run 129 stress 0.07629235 
    ## Run 130 stress 0.07629241 
    ## Run 131 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743372  max resid 0.05100185 
    ## Run 132 stress 0.07378227 
    ## ... Procrustes: rmse 0.01746021  max resid 0.05111234 
    ## Run 133 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001284  max resid 0.0003028373 
    ## ... Similar to previous best
    ## Run 134 stress 0.07365784 
    ## ... Procrustes: rmse 3.932038e-05  max resid 8.708794e-05 
    ## ... Similar to previous best
    ## Run 135 stress 0.07629256 
    ## Run 136 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745763  max resid 0.05108246 
    ## Run 137 stress 0.07365783 
    ## ... Procrustes: rmse 1.039748e-05  max resid 2.486715e-05 
    ## ... Similar to previous best
    ## Run 138 stress 0.08233763 
    ## Run 139 stress 0.07732927 
    ## Run 140 stress 0.08000473 
    ## Run 141 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744578  max resid 0.05104123 
    ## Run 142 stress 0.08233765 
    ## Run 143 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744691  max resid 0.05103155 
    ## Run 144 stress 0.07629241 
    ## Run 145 stress 0.2747107 
    ## Run 146 stress 0.07365785 
    ## ... Procrustes: rmse 7.294741e-05  max resid 0.0001705515 
    ## ... Similar to previous best
    ## Run 147 stress 0.07629239 
    ## Run 148 stress 0.07629233 
    ## Run 149 stress 0.07378231 
    ## ... Procrustes: rmse 0.01744804  max resid 0.05102149 
    ## Run 150 stress 0.07629244 
    ## Run 151 stress 0.08000475 
    ## Run 152 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744735  max resid 0.05105234 
    ## Run 153 stress 0.07365785 
    ## ... Procrustes: rmse 6.657452e-05  max resid 0.0001583416 
    ## ... Similar to previous best
    ## Run 154 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744523  max resid 0.05103474 
    ## Run 155 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174496  max resid 0.05104388 
    ## Run 156 stress 0.07629234 
    ## Run 157 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744658  max resid 0.05101857 
    ## Run 158 stress 0.08233762 
    ## Run 159 stress 0.07629237 
    ## Run 160 stress 0.2560005 
    ## Run 161 stress 0.08288054 
    ## Run 162 stress 0.07365783 
    ## ... Procrustes: rmse 7.253167e-06  max resid 1.611675e-05 
    ## ... Similar to previous best
    ## Run 163 stress 0.07629238 
    ## Run 164 stress 0.08233752 
    ## Run 165 stress 0.07629236 
    ## Run 166 stress 0.07732933 
    ## Run 167 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174437  max resid 0.0510204 
    ## Run 168 stress 0.07629234 
    ## Run 169 stress 0.07365784 
    ## ... Procrustes: rmse 5.007033e-05  max resid 0.0001167126 
    ## ... Similar to previous best
    ## Run 170 stress 0.08288058 
    ## Run 171 stress 0.07629238 
    ## Run 172 stress 0.07365784 
    ## ... Procrustes: rmse 4.350482e-05  max resid 0.0001022001 
    ## ... Similar to previous best
    ## Run 173 stress 0.08000471 
    ## Run 174 stress 0.07629233 
    ## Run 175 stress 0.07732928 
    ## Run 176 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744592  max resid 0.05104368 
    ## Run 177 stress 0.07629244 
    ## Run 178 stress 0.07365783 
    ## ... Procrustes: rmse 1.954998e-05  max resid 4.070751e-05 
    ## ... Similar to previous best
    ## Run 179 stress 0.07378227 
    ## ... Procrustes: rmse 0.01742089  max resid 0.05094064 
    ## Run 180 stress 0.07629234 
    ## Run 181 stress 0.07629238 
    ## Run 182 stress 0.08288051 
    ## Run 183 stress 0.0762924 
    ## Run 184 stress 0.07629233 
    ## Run 185 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744171  max resid 0.05105478 
    ## Run 186 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174445  max resid 0.05103998 
    ## Run 187 stress 0.07732931 
    ## Run 188 stress 0.07378233 
    ## ... Procrustes: rmse 0.01743809  max resid 0.05104929 
    ## Run 189 stress 0.07365783 
    ## ... Procrustes: rmse 9.437243e-06  max resid 2.177822e-05 
    ## ... Similar to previous best
    ## Run 190 stress 0.08233767 
    ## Run 191 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001183219  max resid 0.0002798391 
    ## ... Similar to previous best
    ## Run 192 stress 0.349555 
    ## Run 193 stress 0.08000467 
    ## Run 194 stress 0.08288053 
    ## Run 195 stress 0.07732943 
    ## Run 196 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744468  max resid 0.05103658 
    ## Run 197 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744619  max resid 0.05103963 
    ## Run 198 stress 0.08000468 
    ## Run 199 stress 0.07629236 
    ## Run 200 stress 0.08288046 
    ## Run 201 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744684  max resid 0.05106135 
    ## Run 202 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744472  max resid 0.0510275 
    ## Run 203 stress 0.07365785 
    ## ... Procrustes: rmse 6.917074e-05  max resid 0.0001608254 
    ## ... Similar to previous best
    ## Run 204 stress 0.08233766 
    ## Run 205 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001410217  max resid 0.0003323551 
    ## ... Similar to previous best
    ## Run 206 stress 0.0773293 
    ## Run 207 stress 0.07629242 
    ## Run 208 stress 0.0800047 
    ## Run 209 stress 0.07365786 
    ## ... Procrustes: rmse 8.765354e-05  max resid 0.0002026915 
    ## ... Similar to previous best
    ## Run 210 stress 0.07378229 
    ## ... Procrustes: rmse 0.01746015  max resid 0.05107853 
    ## Run 211 stress 0.07378231 
    ## ... Procrustes: rmse 0.01744391  max resid 0.05106657 
    ## Run 212 stress 0.07629235 
    ## Run 213 stress 0.07732931 
    ## Run 214 stress 0.07365784 
    ## ... Procrustes: rmse 2.681825e-05  max resid 5.584553e-05 
    ## ... Similar to previous best
    ## Run 215 stress 0.0762924 
    ## Run 216 stress 0.08233764 
    ## Run 217 stress 0.07365783 
    ## ... Procrustes: rmse 2.383251e-06  max resid 4.741703e-06 
    ## ... Similar to previous best
    ## Run 218 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001014504  max resid 0.0002387424 
    ## ... Similar to previous best
    ## Run 219 stress 0.07365785 
    ## ... Procrustes: rmse 7.540365e-05  max resid 0.0001781296 
    ## ... Similar to previous best
    ## Run 220 stress 0.08000474 
    ## Run 221 stress 0.07365785 
    ## ... Procrustes: rmse 6.935438e-05  max resid 0.0001637965 
    ## ... Similar to previous best
    ## Run 222 stress 0.07732928 
    ## Run 223 stress 0.07629238 
    ## Run 224 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744437  max resid 0.05103572 
    ## Run 225 stress 0.07629246 
    ## Run 226 stress 0.2918423 
    ## Run 227 stress 0.07629236 
    ## Run 228 stress 0.08233753 
    ## Run 229 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001399277  max resid 0.0003320941 
    ## ... Similar to previous best
    ## Run 230 stress 0.07629233 
    ## Run 231 stress 0.08000467 
    ## Run 232 stress 0.08000468 
    ## Run 233 stress 0.07629235 
    ## Run 234 stress 0.07365784 
    ## ... Procrustes: rmse 3.413181e-05  max resid 8.134993e-05 
    ## ... Similar to previous best
    ## Run 235 stress 0.3495553 
    ## Run 236 stress 0.0762924 
    ## Run 237 stress 0.07365783 
    ## ... Procrustes: rmse 1.660568e-05  max resid 3.610489e-05 
    ## ... Similar to previous best
    ## Run 238 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744673  max resid 0.05105035 
    ## Run 239 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744712  max resid 0.05104027 
    ## Run 240 stress 0.07365784 
    ## ... Procrustes: rmse 3.810032e-05  max resid 8.441811e-05 
    ## ... Similar to previous best
    ## Run 241 stress 0.07365784 
    ## ... Procrustes: rmse 5.282715e-05  max resid 0.0001254305 
    ## ... Similar to previous best
    ## Run 242 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001463754  max resid 0.0003422195 
    ## ... Similar to previous best
    ## Run 243 stress 0.08000468 
    ## Run 244 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001309809  max resid 0.0003132605 
    ## ... Similar to previous best
    ## Run 245 stress 0.07365786 
    ## ... Procrustes: rmse 3.353201e-05  max resid 6.506484e-05 
    ## ... Similar to previous best
    ## Run 246 stress 0.08288039 
    ## Run 247 stress 0.08000467 
    ## Run 248 stress 0.07629238 
    ## Run 249 stress 0.08000467 
    ## Run 250 stress 0.08233758 
    ## Run 251 stress 0.07378231 
    ## ... Procrustes: rmse 0.01744415  max resid 0.0510286 
    ## Run 252 stress 0.07365784 
    ## ... Procrustes: rmse 4.74605e-05  max resid 0.0001134649 
    ## ... Similar to previous best
    ## Run 253 stress 0.07365783 
    ## ... Procrustes: rmse 1.00071e-05  max resid 2.260616e-05 
    ## ... Similar to previous best
    ## Run 254 stress 0.07732941 
    ## Run 255 stress 0.07629235 
    ## Run 256 stress 0.07732928 
    ## Run 257 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001113493  max resid 0.0002658267 
    ## ... Similar to previous best
    ## Run 258 stress 0.3135765 
    ## Run 259 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744958  max resid 0.05107159 
    ## Run 260 stress 0.07365783 
    ## ... Procrustes: rmse 7.447834e-06  max resid 1.62045e-05 
    ## ... Similar to previous best
    ## Run 261 stress 0.07629254 
    ## Run 262 stress 0.07629237 
    ## Run 263 stress 0.07365784 
    ## ... Procrustes: rmse 6.293226e-05  max resid 0.0001484646 
    ## ... Similar to previous best
    ## Run 264 stress 0.07732925 
    ## Run 265 stress 0.0737823 
    ## ... Procrustes: rmse 0.01745659  max resid 0.05110763 
    ## Run 266 stress 0.07365785 
    ## ... Procrustes: rmse 6.793754e-05  max resid 0.0001583319 
    ## ... Similar to previous best
    ## Run 267 stress 0.08233755 
    ## Run 268 stress 0.07629244 
    ## Run 269 stress 0.07629237 
    ## Run 270 stress 0.07365785 
    ## ... Procrustes: rmse 6.923427e-05  max resid 0.0001655281 
    ## ... Similar to previous best
    ## Run 271 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744171  max resid 0.05105337 
    ## Run 272 stress 0.07365785 
    ## ... Procrustes: rmse 2.234402e-05  max resid 5.290415e-05 
    ## ... Similar to previous best
    ## Run 273 stress 0.07365784 
    ## ... Procrustes: rmse 5.356364e-05  max resid 0.0001282009 
    ## ... Similar to previous best
    ## Run 274 stress 0.07365784 
    ## ... Procrustes: rmse 4.777245e-05  max resid 0.0001130681 
    ## ... Similar to previous best
    ## Run 275 stress 0.07732946 
    ## Run 276 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744563  max resid 0.05103462 
    ## Run 277 stress 0.07629242 
    ## Run 278 stress 0.08233756 
    ## Run 279 stress 0.07629237 
    ## Run 280 stress 0.07629235 
    ## Run 281 stress 0.07365784 
    ## ... Procrustes: rmse 4.238677e-05  max resid 9.845139e-05 
    ## ... Similar to previous best
    ## Run 282 stress 0.07732941 
    ## Run 283 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744595  max resid 0.05102786 
    ## Run 284 stress 0.08288034 
    ## Run 285 stress 0.07365788 
    ## ... Procrustes: rmse 2.906192e-05  max resid 6.096063e-05 
    ## ... Similar to previous best
    ## Run 286 stress 0.07629234 
    ## Run 287 stress 0.07629247 
    ## Run 288 stress 0.08233773 
    ## Run 289 stress 0.07365783 
    ## ... Procrustes: rmse 1.582133e-05  max resid 3.630201e-05 
    ## ... Similar to previous best
    ## Run 290 stress 0.07365783 
    ## ... Procrustes: rmse 2.914884e-05  max resid 6.841536e-05 
    ## ... Similar to previous best
    ## Run 291 stress 0.07629237 
    ## Run 292 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745652  max resid 0.05107558 
    ## Run 293 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744853  max resid 0.05105666 
    ## Run 294 stress 0.07629249 
    ## Run 295 stress 0.08000467 
    ## Run 296 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744683  max resid 0.05105702 
    ## Run 297 stress 0.3578479 
    ## Run 298 stress 0.08000477 
    ## Run 299 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745204  max resid 0.05107212 
    ## Run 300 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744687  max resid 0.05102075 
    ## Run 301 stress 0.07365783 
    ## ... Procrustes: rmse 9.268555e-06  max resid 2.041644e-05 
    ## ... Similar to previous best
    ## Run 302 stress 0.08000468 
    ## Run 303 stress 0.07365785 
    ## ... Procrustes: rmse 7.284882e-05  max resid 0.0001686866 
    ## ... Similar to previous best
    ## Run 304 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744694  max resid 0.05104949 
    ## Run 305 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001000275  max resid 0.000237326 
    ## ... Similar to previous best
    ## Run 306 stress 0.07732931 
    ## Run 307 stress 0.07378231 
    ## ... Procrustes: rmse 0.0174394  max resid 0.05105073 
    ## Run 308 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001012374  max resid 0.0002314792 
    ## ... Similar to previous best
    ## Run 309 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743499  max resid 0.05100334 
    ## Run 310 stress 0.08000478 
    ## Run 311 stress 0.08000468 
    ## Run 312 stress 0.07365783 
    ## ... Procrustes: rmse 3.170984e-05  max resid 7.466246e-05 
    ## ... Similar to previous best
    ## Run 313 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744672  max resid 0.05106405 
    ## Run 314 stress 0.07629233 
    ## Run 315 stress 0.07629235 
    ## Run 316 stress 0.07365784 
    ## ... Procrustes: rmse 4.079588e-05  max resid 9.729959e-05 
    ## ... Similar to previous best
    ## Run 317 stress 0.08000476 
    ## Run 318 stress 0.07629249 
    ## Run 319 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001436062  max resid 0.0003406656 
    ## ... Similar to previous best
    ## Run 320 stress 0.08000474 
    ## Run 321 stress 0.08000479 
    ## Run 322 stress 0.07365784 
    ## ... Procrustes: rmse 4.63676e-05  max resid 0.0001042162 
    ## ... Similar to previous best
    ## Run 323 stress 0.0762924 
    ## Run 324 stress 0.08000474 
    ## Run 325 stress 0.07629241 
    ## Run 326 stress 0.0823376 
    ## Run 327 stress 0.08000467 
    ## Run 328 stress 0.08000471 
    ## Run 329 stress 0.08000473 
    ## Run 330 stress 0.07365785 
    ## ... Procrustes: rmse 5.901739e-05  max resid 0.0001412516 
    ## ... Similar to previous best
    ## Run 331 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744715  max resid 0.05105636 
    ## Run 332 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745138  max resid 0.05106944 
    ## Run 333 stress 0.07365785 
    ## ... Procrustes: rmse 5.742525e-05  max resid 0.0001365157 
    ## ... Similar to previous best
    ## Run 334 stress 0.07629246 
    ## Run 335 stress 0.07732928 
    ## Run 336 stress 0.08233769 
    ## Run 337 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744619  max resid 0.05105584 
    ## Run 338 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745071  max resid 0.05107079 
    ## Run 339 stress 0.08000477 
    ## Run 340 stress 0.07629234 
    ## Run 341 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744583  max resid 0.05106223 
    ## Run 342 stress 0.07378229 
    ## ... Procrustes: rmse 0.017443  max resid 0.05105806 
    ## Run 343 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744034  max resid 0.05104951 
    ## Run 344 stress 0.0800047 
    ## Run 345 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744617  max resid 0.05105862 
    ## Run 346 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174519  max resid 0.05108029 
    ## Run 347 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744949  max resid 0.05106549 
    ## Run 348 stress 0.07629246 
    ## Run 349 stress 0.07365784 
    ## ... Procrustes: rmse 1.287915e-05  max resid 2.66696e-05 
    ## ... Similar to previous best
    ## Run 350 stress 0.07732938 
    ## Run 351 stress 0.07365783 
    ## ... Procrustes: rmse 2.293948e-05  max resid 5.189987e-05 
    ## ... Similar to previous best
    ## Run 352 stress 0.07365783 
    ## ... Procrustes: rmse 2.302508e-05  max resid 5.530538e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.07365784 
    ## ... Procrustes: rmse 3.533244e-05  max resid 8.442264e-05 
    ## ... Similar to previous best
    ## Run 354 stress 0.0800047 
    ## Run 355 stress 0.07629239 
    ## Run 356 stress 0.07365787 
    ## ... Procrustes: rmse 8.851741e-05  max resid 0.0002047165 
    ## ... Similar to previous best
    ## Run 357 stress 0.07629234 
    ## Run 358 stress 0.07629243 
    ## Run 359 stress 0.07365785 
    ## ... Procrustes: rmse 4.466862e-05  max resid 0.0001065907 
    ## ... Similar to previous best
    ## Run 360 stress 0.07365783 
    ## ... Procrustes: rmse 6.418846e-06  max resid 1.512943e-05 
    ## ... Similar to previous best
    ## Run 361 stress 0.08000467 
    ## Run 362 stress 0.08233752 
    ## Run 363 stress 0.3578484 
    ## Run 364 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744606  max resid 0.0510389 
    ## Run 365 stress 0.0773294 
    ## Run 366 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001055064  max resid 0.0002392652 
    ## ... Similar to previous best
    ## Run 367 stress 0.07365783 
    ## ... Procrustes: rmse 1.220927e-05  max resid 2.801013e-05 
    ## ... Similar to previous best
    ## Run 368 stress 0.07365787 
    ## ... Procrustes: rmse 9.928963e-05  max resid 0.0002345825 
    ## ... Similar to previous best
    ## Run 369 stress 0.2879036 
    ## Run 370 stress 0.08000469 
    ## Run 371 stress 0.07629234 
    ## Run 372 stress 0.07365783 
    ## ... Procrustes: rmse 2.469487e-05  max resid 5.751712e-05 
    ## ... Similar to previous best
    ## Run 373 stress 0.08233758 
    ## Run 374 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744684  max resid 0.05105625 
    ## Run 375 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744711  max resid 0.05105947 
    ## Run 376 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744985  max resid 0.05106986 
    ## Run 377 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744516  max resid 0.05104106 
    ## Run 378 stress 0.07365786 
    ## ... Procrustes: rmse 8.562468e-05  max resid 0.0002022877 
    ## ... Similar to previous best
    ## Run 379 stress 0.07629244 
    ## Run 380 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745448  max resid 0.05109855 
    ## Run 381 stress 0.07629247 
    ## Run 382 stress 0.08000467 
    ## Run 383 stress 0.07732939 
    ## Run 384 stress 0.07629235 
    ## Run 385 stress 0.07378232 
    ## ... Procrustes: rmse 0.01744129  max resid 0.05098345 
    ## Run 386 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001151655  max resid 0.0002711519 
    ## ... Similar to previous best
    ## Run 387 stress 0.07365784 
    ## ... Procrustes: rmse 1.817409e-05  max resid 3.91469e-05 
    ## ... Similar to previous best
    ## Run 388 stress 0.07365786 
    ## ... Procrustes: rmse 9.832139e-05  max resid 0.0002320902 
    ## ... Similar to previous best
    ## Run 389 stress 0.08233758 
    ## Run 390 stress 0.07629235 
    ## Run 391 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001244197  max resid 0.0002901827 
    ## ... Similar to previous best
    ## Run 392 stress 0.07629233 
    ## Run 393 stress 0.07629242 
    ## Run 394 stress 0.07629236 
    ## Run 395 stress 0.08000471 
    ## Run 396 stress 0.07629233 
    ## Run 397 stress 0.07629241 
    ## Run 398 stress 0.07378229 
    ## ... Procrustes: rmse 0.0174633  max resid 0.05113517 
    ## Run 399 stress 0.08000474 
    ## Run 400 stress 0.0773294 
    ## Run 401 stress 0.07629233 
    ## Run 402 stress 0.08233761 
    ## Run 403 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745943  max resid 0.05112055 
    ## Run 404 stress 0.08000468 
    ## Run 405 stress 0.07365785 
    ## ... Procrustes: rmse 6.726388e-05  max resid 0.0001592586 
    ## ... Similar to previous best
    ## Run 406 stress 0.07732928 
    ## Run 407 stress 0.07378229 
    ## ... Procrustes: rmse 0.0174439  max resid 0.05105647 
    ## Run 408 stress 0.08000467 
    ## Run 409 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745231  max resid 0.05108153 
    ## Run 410 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001117797  max resid 0.0002641493 
    ## ... Similar to previous best
    ## Run 411 stress 0.07365783 
    ## ... Procrustes: rmse 9.207622e-06  max resid 2.203081e-05 
    ## ... Similar to previous best
    ## Run 412 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744919  max resid 0.05105389 
    ## Run 413 stress 0.08233755 
    ## Run 414 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744151  max resid 0.05104063 
    ## Run 415 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 4.290845e-06  max resid 7.223673e-06 
    ## ... Similar to previous best
    ## Run 416 stress 0.07629235 
    ## Run 417 stress 0.07365787 
    ## ... Procrustes: rmse 9.500459e-05  max resid 0.0002257119 
    ## ... Similar to previous best
    ## Run 418 stress 0.07732935 
    ## Run 419 stress 0.08000471 
    ## Run 420 stress 0.08000476 
    ## Run 421 stress 0.08000467 
    ## Run 422 stress 0.07365784 
    ## ... Procrustes: rmse 2.424781e-05  max resid 5.310111e-05 
    ## ... Similar to previous best
    ## Run 423 stress 0.07629235 
    ## Run 424 stress 0.07365783 
    ## ... Procrustes: rmse 2.091604e-05  max resid 4.787206e-05 
    ## ... Similar to previous best
    ## Run 425 stress 0.08000468 
    ## Run 426 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744512  max resid 0.05105914 
    ## Run 427 stress 0.0828804 
    ## Run 428 stress 0.08000471 
    ## Run 429 stress 0.0773293 
    ## Run 430 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744105  max resid 0.05099807 
    ## Run 431 stress 0.08233769 
    ## Run 432 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744367  max resid 0.05102225 
    ## Run 433 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744851  max resid 0.05106825 
    ## Run 434 stress 0.08000478 
    ## Run 435 stress 0.07365784 
    ## ... Procrustes: rmse 5.390116e-05  max resid 0.00012368 
    ## ... Similar to previous best
    ## Run 436 stress 0.07365785 
    ## ... Procrustes: rmse 7.424289e-05  max resid 0.0001750388 
    ## ... Similar to previous best
    ## Run 437 stress 0.07365784 
    ## ... Procrustes: rmse 5.151653e-05  max resid 0.0001177016 
    ## ... Similar to previous best
    ## Run 438 stress 0.07629234 
    ## Run 439 stress 0.3493349 
    ## Run 440 stress 0.08000481 
    ## Run 441 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001461076  max resid 0.0003467042 
    ## ... Similar to previous best
    ## Run 442 stress 0.08000467 
    ## Run 443 stress 0.07629234 
    ## Run 444 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745121  max resid 0.05108716 
    ## Run 445 stress 0.08288054 
    ## Run 446 stress 0.07629243 
    ## Run 447 stress 0.07732935 
    ## Run 448 stress 0.07629248 
    ## Run 449 stress 0.07365783 
    ## ... Procrustes: rmse 3.473953e-06  max resid 5.876491e-06 
    ## ... Similar to previous best
    ## Run 450 stress 0.07365784 
    ## ... Procrustes: rmse 2.668022e-05  max resid 6.209969e-05 
    ## ... Similar to previous best
    ## Run 451 stress 0.07629242 
    ## Run 452 stress 0.07629234 
    ## Run 453 stress 0.08288052 
    ## Run 454 stress 0.0737823 
    ## ... Procrustes: rmse 0.01745429  max resid 0.05110718 
    ## Run 455 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744681  max resid 0.05105313 
    ## Run 456 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001178135  max resid 0.0002657216 
    ## ... Similar to previous best
    ## Run 457 stress 0.07629238 
    ## Run 458 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744537  max resid 0.05103674 
    ## Run 459 stress 0.07629234 
    ## Run 460 stress 0.0823376 
    ## Run 461 stress 0.07629234 
    ## Run 462 stress 0.07629234 
    ## Run 463 stress 0.08000468 
    ## Run 464 stress 0.08000468 
    ## Run 465 stress 0.07365785 
    ## ... Procrustes: rmse 6.333216e-05  max resid 0.0001488038 
    ## ... Similar to previous best
    ## Run 466 stress 0.07629241 
    ## Run 467 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001573585  max resid 0.0003698344 
    ## ... Similar to previous best
    ## Run 468 stress 0.07732935 
    ## Run 469 stress 0.07378229 
    ## ... Procrustes: rmse 0.01747246  max resid 0.05114139 
    ## Run 470 stress 0.07365785 
    ## ... Procrustes: rmse 7.03887e-05  max resid 0.0001661333 
    ## ... Similar to previous best
    ## Run 471 stress 0.07629241 
    ## Run 472 stress 0.07365784 
    ## ... Procrustes: rmse 4.254867e-05  max resid 0.0001008526 
    ## ... Similar to previous best
    ## Run 473 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001460491  max resid 0.0003471382 
    ## ... Similar to previous best
    ## Run 474 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745373  max resid 0.05108913 
    ## Run 475 stress 0.07629237 
    ## Run 476 stress 0.07629235 
    ## Run 477 stress 0.08000469 
    ## Run 478 stress 0.08000474 
    ## Run 479 stress 0.07378233 
    ## ... Procrustes: rmse 0.01743747  max resid 0.05097168 
    ## Run 480 stress 0.07365784 
    ## ... Procrustes: rmse 4.494502e-05  max resid 0.0001057181 
    ## ... Similar to previous best
    ## Run 481 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744686  max resid 0.05105492 
    ## Run 482 stress 0.07365795 
    ## ... Procrustes: rmse 0.0001861264  max resid 0.0004371024 
    ## ... Similar to previous best
    ## Run 483 stress 0.07365784 
    ## ... Procrustes: rmse 3.12242e-05  max resid 6.480044e-05 
    ## ... Similar to previous best
    ## Run 484 stress 0.07629242 
    ## Run 485 stress 0.07365783 
    ## ... Procrustes: rmse 2.363226e-05  max resid 5.59773e-05 
    ## ... Similar to previous best
    ## Run 486 stress 0.08000468 
    ## Run 487 stress 0.07365787 
    ## ... Procrustes: rmse 9.534861e-05  max resid 0.0002220014 
    ## ... Similar to previous best
    ## Run 488 stress 0.07629235 
    ## Run 489 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744303  max resid 0.05103213 
    ## Run 490 stress 0.07629235 
    ## Run 491 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745362  max resid 0.05109924 
    ## Run 492 stress 0.0762924 
    ## Run 493 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744303  max resid 0.05102853 
    ## Run 494 stress 0.08000471 
    ## Run 495 stress 0.08233762 
    ## Run 496 stress 0.08000472 
    ## Run 497 stress 0.07629234 
    ## Run 498 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744349  max resid 0.0510247 
    ## Run 499 stress 0.07629241 
    ## Run 500 stress 0.07365785 
    ## ... Procrustes: rmse 7.103789e-05  max resid 0.0001682031 
    ## ... Similar to previous best
    ## *** Best solution repeated 22 times

``` r
# Stratified lakes and ocean sites
SD_beta_geo_SO_NMDS <- metaMDS(SD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07970526 
    ## Run 1 stress 0.07970525 
    ## ... New best solution
    ## ... Procrustes: rmse 3.744714e-05  max resid 6.001898e-05 
    ## ... Similar to previous best
    ## Run 2 stress 0.2911511 
    ## Run 3 stress 0.07250812 
    ## ... New best solution
    ## ... Procrustes: rmse 0.08994375  max resid 0.2576108 
    ## Run 4 stress 0.07250812 
    ## ... New best solution
    ## ... Procrustes: rmse 3.637515e-05  max resid 0.0001026144 
    ## ... Similar to previous best
    ## Run 5 stress 0.06942777 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04030834  max resid 0.1229513 
    ## Run 6 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 3.032153e-05  max resid 7.806635e-05 
    ## ... Similar to previous best
    ## Run 7 stress 0.07970525 
    ## Run 8 stress 0.07428313 
    ## Run 9 stress 0.07970525 
    ## Run 10 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321712  max resid 0.03324989 
    ## Run 11 stress 0.08448457 
    ## Run 12 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325212  max resid 0.0333259 
    ## Run 13 stress 0.07250814 
    ## Run 14 stress 0.08340287 
    ## Run 15 stress 0.07844931 
    ## Run 16 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324352  max resid 0.0333059 
    ## Run 17 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 1.822204e-05  max resid 4.785371e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.07428313 
    ## Run 19 stress 0.07844918 
    ## Run 20 stress 0.07970525 
    ## Run 21 stress 0.08340286 
    ## Run 22 stress 0.06978194 
    ## ... Procrustes: rmse 0.0132641  max resid 0.03334894 
    ## Run 23 stress 0.06978201 
    ## ... Procrustes: rmse 0.01328547  max resid 0.03338984 
    ## Run 24 stress 0.07428315 
    ## Run 25 stress 0.08340294 
    ## Run 26 stress 0.08448439 
    ## Run 27 stress 0.07250816 
    ## Run 28 stress 0.07428316 
    ## Run 29 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326533  max resid 0.03335112 
    ## Run 30 stress 0.08340288 
    ## Run 31 stress 0.08448437 
    ## Run 32 stress 0.07428317 
    ## Run 33 stress 0.08340287 
    ## Run 34 stress 0.07970525 
    ## Run 35 stress 0.08340287 
    ## Run 36 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 5.761832e-06  max resid 1.619031e-05 
    ## ... Similar to previous best
    ## Run 37 stress 0.07428313 
    ## Run 38 stress 0.07250813 
    ## Run 39 stress 0.08340287 
    ## Run 40 stress 0.07970526 
    ## Run 41 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320877  max resid 0.03323028 
    ## Run 42 stress 0.07428313 
    ## Run 43 stress 0.07250812 
    ## Run 44 stress 0.07970525 
    ## Run 45 stress 0.08340287 
    ## Run 46 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326108  max resid 0.0333428 
    ## Run 47 stress 0.07250813 
    ## Run 48 stress 0.06978196 
    ## ... Procrustes: rmse 0.01327048  max resid 0.033362 
    ## Run 49 stress 0.07844914 
    ## Run 50 stress 0.07250819 
    ## Run 51 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324625  max resid 0.03331069 
    ## Run 52 stress 0.07970525 
    ## Run 53 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326486  max resid 0.03334873 
    ## Run 54 stress 0.06942776 
    ## ... Procrustes: rmse 2.876478e-06  max resid 6.647171e-06 
    ## ... Similar to previous best
    ## Run 55 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325109  max resid 0.03332305 
    ## Run 56 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322085  max resid 0.03326398 
    ## Run 57 stress 0.08340292 
    ## Run 58 stress 0.07428313 
    ## Run 59 stress 0.08340289 
    ## Run 60 stress 0.06942777 
    ## ... Procrustes: rmse 4.820539e-05  max resid 0.0001234661 
    ## ... Similar to previous best
    ## Run 61 stress 0.08340287 
    ## Run 62 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320799  max resid 0.03323119 
    ## Run 63 stress 0.06942776 
    ## ... Procrustes: rmse 1.741008e-05  max resid 4.327021e-05 
    ## ... Similar to previous best
    ## Run 64 stress 0.07428318 
    ## Run 65 stress 0.0834029 
    ## Run 66 stress 0.06942777 
    ## ... Procrustes: rmse 4.339948e-05  max resid 0.0001111063 
    ## ... Similar to previous best
    ## Run 67 stress 0.07428316 
    ## Run 68 stress 0.08340288 
    ## Run 69 stress 0.07250812 
    ## Run 70 stress 0.06978197 
    ## ... Procrustes: rmse 0.0131817  max resid 0.03317729 
    ## Run 71 stress 0.08448451 
    ## Run 72 stress 0.06942776 
    ## ... Procrustes: rmse 2.950368e-05  max resid 7.533743e-05 
    ## ... Similar to previous best
    ## Run 73 stress 0.08340287 
    ## Run 74 stress 0.0834029 
    ## Run 75 stress 0.07970526 
    ## Run 76 stress 0.06942777 
    ## ... Procrustes: rmse 4.892274e-05  max resid 0.0001252281 
    ## ... Similar to previous best
    ## Run 77 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.123429e-06  max resid 4.938975e-06 
    ## ... Similar to previous best
    ## Run 78 stress 0.07428315 
    ## Run 79 stress 0.06942776 
    ## ... Procrustes: rmse 5.144077e-06  max resid 1.351022e-05 
    ## ... Similar to previous best
    ## Run 80 stress 0.06942776 
    ## ... Procrustes: rmse 1.01453e-05  max resid 2.480991e-05 
    ## ... Similar to previous best
    ## Run 81 stress 0.08340287 
    ## Run 82 stress 0.07844958 
    ## Run 83 stress 0.07428313 
    ## Run 84 stress 0.08340286 
    ## Run 85 stress 0.08340288 
    ## Run 86 stress 0.06942777 
    ## ... Procrustes: rmse 3.15606e-05  max resid 8.422446e-05 
    ## ... Similar to previous best
    ## Run 87 stress 0.08448438 
    ## Run 88 stress 0.07250812 
    ## Run 89 stress 0.06942778 
    ## ... Procrustes: rmse 4.922372e-05  max resid 0.0001275284 
    ## ... Similar to previous best
    ## Run 90 stress 0.08448452 
    ## Run 91 stress 0.06942777 
    ## ... Procrustes: rmse 4.324357e-05  max resid 0.0001118728 
    ## ... Similar to previous best
    ## Run 92 stress 0.07428316 
    ## Run 93 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326112  max resid 0.03334514 
    ## Run 94 stress 0.07970526 
    ## Run 95 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324618  max resid 0.03331098 
    ## Run 96 stress 0.07250812 
    ## Run 97 stress 0.06942776 
    ## ... Procrustes: rmse 2.869918e-06  max resid 7.763553e-06 
    ## ... Similar to previous best
    ## Run 98 stress 0.06942776 
    ## ... Procrustes: rmse 2.910845e-05  max resid 7.502308e-05 
    ## ... Similar to previous best
    ## Run 99 stress 0.06942776 
    ## ... Procrustes: rmse 1.717486e-05  max resid 3.476715e-05 
    ## ... Similar to previous best
    ## Run 100 stress 0.08340293 
    ## Run 101 stress 0.07428315 
    ## Run 102 stress 0.07250813 
    ## Run 103 stress 0.06942776 
    ## ... Procrustes: rmse 4.372349e-06  max resid 1.043989e-05 
    ## ... Similar to previous best
    ## Run 104 stress 0.07428314 
    ## Run 105 stress 0.08340286 
    ## Run 106 stress 0.07428315 
    ## Run 107 stress 0.06978196 
    ## ... Procrustes: rmse 0.01327072  max resid 0.03336155 
    ## Run 108 stress 0.07428314 
    ## Run 109 stress 0.08340286 
    ## Run 110 stress 0.07250812 
    ## Run 111 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321484  max resid 0.03324565 
    ## Run 112 stress 0.07250813 
    ## Run 113 stress 0.07250814 
    ## Run 114 stress 0.06978192 
    ## ... Procrustes: rmse 0.0132535  max resid 0.03332569 
    ## Run 115 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326671  max resid 0.03335382 
    ## Run 116 stress 0.07250812 
    ## Run 117 stress 0.07844939 
    ## Run 118 stress 0.0834029 
    ## Run 119 stress 0.06978199 
    ## ... Procrustes: rmse 0.01327869  max resid 0.03337729 
    ## Run 120 stress 0.07250815 
    ## Run 121 stress 0.08340288 
    ## Run 122 stress 0.07428321 
    ## Run 123 stress 0.08340294 
    ## Run 124 stress 0.07428314 
    ## Run 125 stress 0.06942777 
    ## ... Procrustes: rmse 3.986908e-05  max resid 0.0001024342 
    ## ... Similar to previous best
    ## Run 126 stress 0.06942776 
    ## ... Procrustes: rmse 9.123769e-06  max resid 2.467238e-05 
    ## ... Similar to previous best
    ## Run 127 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326991  max resid 0.03335793 
    ## Run 128 stress 0.06942776 
    ## ... Procrustes: rmse 5.842482e-07  max resid 1.08859e-06 
    ## ... Similar to previous best
    ## Run 129 stress 0.07250812 
    ## Run 130 stress 0.08340287 
    ## Run 131 stress 0.07970525 
    ## Run 132 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323124  max resid 0.03328001 
    ## Run 133 stress 0.07970526 
    ## Run 134 stress 0.07428313 
    ## Run 135 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326312  max resid 0.03334908 
    ## Run 136 stress 0.08340292 
    ## Run 137 stress 0.07428315 
    ## Run 138 stress 0.06978192 
    ## ... Procrustes: rmse 0.0132482  max resid 0.0333162 
    ## Run 139 stress 0.08448435 
    ## Run 140 stress 0.07250812 
    ## Run 141 stress 0.07250812 
    ## Run 142 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326284  max resid 0.03334771 
    ## Run 143 stress 0.06978191 
    ## ... Procrustes: rmse 0.0132116  max resid 0.03324095 
    ## Run 144 stress 0.08340293 
    ## Run 145 stress 0.07970527 
    ## Run 146 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323895  max resid 0.03329767 
    ## Run 147 stress 0.07250815 
    ## Run 148 stress 0.08340289 
    ## Run 149 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322003  max resid 0.03325588 
    ## Run 150 stress 0.07428313 
    ## Run 151 stress 0.08340298 
    ## Run 152 stress 0.0742832 
    ## Run 153 stress 0.06942776 
    ## ... Procrustes: rmse 3.332215e-05  max resid 8.636408e-05 
    ## ... Similar to previous best
    ## Run 154 stress 0.08340287 
    ## Run 155 stress 0.07428313 
    ## Run 156 stress 0.07428314 
    ## Run 157 stress 0.07428313 
    ## Run 158 stress 0.06942777 
    ## ... Procrustes: rmse 1.816968e-05  max resid 5.388662e-05 
    ## ... Similar to previous best
    ## Run 159 stress 0.06942776 
    ## ... Procrustes: rmse 3.506376e-05  max resid 9.067705e-05 
    ## ... Similar to previous best
    ## Run 160 stress 0.08340292 
    ## Run 161 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320416  max resid 0.03322125 
    ## Run 162 stress 0.08340293 
    ## Run 163 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325207  max resid 0.0333222 
    ## Run 164 stress 0.07250812 
    ## Run 165 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325788  max resid 0.03333515 
    ## Run 166 stress 0.06942777 
    ## ... Procrustes: rmse 3.804751e-05  max resid 9.924835e-05 
    ## ... Similar to previous best
    ## Run 167 stress 0.07428313 
    ## Run 168 stress 0.08340287 
    ## Run 169 stress 0.07970525 
    ## Run 170 stress 0.06978189 
    ## ... Procrustes: rmse 0.0132234  max resid 0.0332653 
    ## Run 171 stress 0.06942776 
    ## ... Procrustes: rmse 1.09487e-05  max resid 2.895472e-05 
    ## ... Similar to previous best
    ## Run 172 stress 0.06942777 
    ## ... Procrustes: rmse 4.447951e-05  max resid 0.0001148821 
    ## ... Similar to previous best
    ## Run 173 stress 0.07428313 
    ## Run 174 stress 0.08340287 
    ## Run 175 stress 0.08340286 
    ## Run 176 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322698  max resid 0.03326776 
    ## Run 177 stress 0.08448441 
    ## Run 178 stress 0.06942776 
    ## ... Procrustes: rmse 3.322837e-05  max resid 8.610934e-05 
    ## ... Similar to previous best
    ## Run 179 stress 0.06942776 
    ## ... Procrustes: rmse 8.110181e-06  max resid 2.163991e-05 
    ## ... Similar to previous best
    ## Run 180 stress 0.07250813 
    ## Run 181 stress 0.07250814 
    ## Run 182 stress 0.07428313 
    ## Run 183 stress 0.08340287 
    ## Run 184 stress 0.06942776 
    ## ... Procrustes: rmse 2.373417e-05  max resid 6.149608e-05 
    ## ... Similar to previous best
    ## Run 185 stress 0.07250813 
    ## Run 186 stress 0.07428319 
    ## Run 187 stress 0.07250812 
    ## Run 188 stress 0.08448436 
    ## Run 189 stress 0.06942776 
    ## ... Procrustes: rmse 1.968706e-05  max resid 5.238705e-05 
    ## ... Similar to previous best
    ## Run 190 stress 0.07428313 
    ## Run 191 stress 0.08340287 
    ## Run 192 stress 0.07250814 
    ## Run 193 stress 0.06978191 
    ## ... Procrustes: rmse 0.01320083  max resid 0.03321614 
    ## Run 194 stress 0.0834029 
    ## Run 195 stress 0.069782 
    ## ... Procrustes: rmse 0.01327996  max resid 0.03338249 
    ## Run 196 stress 0.06942776 
    ## ... Procrustes: rmse 2.995869e-06  max resid 7.368419e-06 
    ## ... Similar to previous best
    ## Run 197 stress 0.07428313 
    ## Run 198 stress 0.07844934 
    ## Run 199 stress 0.08448454 
    ## Run 200 stress 0.08340289 
    ## Run 201 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132322  max resid 0.03328015 
    ## Run 202 stress 0.06942777 
    ## ... Procrustes: rmse 4.35653e-05  max resid 0.0001127068 
    ## ... Similar to previous best
    ## Run 203 stress 0.069782 
    ## ... Procrustes: rmse 0.01328377  max resid 0.03338804 
    ## Run 204 stress 0.07970525 
    ## Run 205 stress 0.07428315 
    ## Run 206 stress 0.0742832 
    ## Run 207 stress 0.07970526 
    ## Run 208 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325579  max resid 0.03333236 
    ## Run 209 stress 0.07428314 
    ## Run 210 stress 0.06978193 
    ## ... Procrustes: rmse 0.01318655  max resid 0.03318456 
    ## Run 211 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325742  max resid 0.03333578 
    ## Run 212 stress 0.08448446 
    ## Run 213 stress 0.07250812 
    ## Run 214 stress 0.07970525 
    ## Run 215 stress 0.07250815 
    ## Run 216 stress 0.08340291 
    ## Run 217 stress 0.07250815 
    ## Run 218 stress 0.06942776 
    ## ... Procrustes: rmse 2.1507e-05  max resid 5.595701e-05 
    ## ... Similar to previous best
    ## Run 219 stress 0.07428312 
    ## Run 220 stress 0.06978203 
    ## ... Procrustes: rmse 0.01315463  max resid 0.033112 
    ## Run 221 stress 0.07428321 
    ## Run 222 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 1.910182e-06  max resid 5.001519e-06 
    ## ... Similar to previous best
    ## Run 223 stress 0.07250813 
    ## Run 224 stress 0.06942777 
    ## ... Procrustes: rmse 3.159365e-05  max resid 8.334135e-05 
    ## ... Similar to previous best
    ## Run 225 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327082  max resid 0.03336329 
    ## Run 226 stress 0.08340286 
    ## Run 227 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323501  max resid 0.03328801 
    ## Run 228 stress 0.0834029 
    ## Run 229 stress 0.06942776 
    ## ... Procrustes: rmse 1.65485e-05  max resid 4.273327e-05 
    ## ... Similar to previous best
    ## Run 230 stress 0.07844965 
    ## Run 231 stress 0.07428316 
    ## Run 232 stress 0.08340287 
    ## Run 233 stress 0.07250813 
    ## Run 234 stress 0.07970526 
    ## Run 235 stress 0.07428312 
    ## Run 236 stress 0.06942776 
    ## ... Procrustes: rmse 2.431539e-05  max resid 6.309389e-05 
    ## ... Similar to previous best
    ## Run 237 stress 0.07428313 
    ## Run 238 stress 0.07250812 
    ## Run 239 stress 0.07428313 
    ## Run 240 stress 0.07970525 
    ## Run 241 stress 0.08340286 
    ## Run 242 stress 0.08340291 
    ## Run 243 stress 0.07970526 
    ## Run 244 stress 0.06942776 
    ## ... Procrustes: rmse 2.805569e-05  max resid 7.218853e-05 
    ## ... Similar to previous best
    ## Run 245 stress 0.07428317 
    ## Run 246 stress 0.06942777 
    ## ... Procrustes: rmse 3.620996e-05  max resid 9.343677e-05 
    ## ... Similar to previous best
    ## Run 247 stress 0.07250814 
    ## Run 248 stress 0.07970525 
    ## Run 249 stress 0.07428313 
    ## Run 250 stress 0.07970525 
    ## Run 251 stress 0.07970525 
    ## Run 252 stress 0.07428315 
    ## Run 253 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323421  max resid 0.03328588 
    ## Run 254 stress 0.06942776 
    ## ... Procrustes: rmse 1.093768e-05  max resid 2.83517e-05 
    ## ... Similar to previous best
    ## Run 255 stress 0.07844946 
    ## Run 256 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321418  max resid 0.03324327 
    ## Run 257 stress 0.06978196 
    ## ... Procrustes: rmse 0.01327009  max resid 0.03335902 
    ## Run 258 stress 0.07844919 
    ## Run 259 stress 0.07428316 
    ## Run 260 stress 0.07970525 
    ## Run 261 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319837  max resid 0.03321024 
    ## Run 262 stress 0.07970525 
    ## Run 263 stress 0.08340287 
    ## Run 264 stress 0.07970525 
    ## Run 265 stress 0.08340291 
    ## Run 266 stress 0.07970525 
    ## Run 267 stress 0.07970525 
    ## Run 268 stress 0.07428314 
    ## Run 269 stress 0.06942776 
    ## ... Procrustes: rmse 8.107614e-06  max resid 2.098706e-05 
    ## ... Similar to previous best
    ## Run 270 stress 0.08340288 
    ## Run 271 stress 0.08448452 
    ## Run 272 stress 0.07428313 
    ## Run 273 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325744  max resid 0.03333289 
    ## Run 274 stress 0.08340286 
    ## Run 275 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326087  max resid 0.03334225 
    ## Run 276 stress 0.06942776 
    ## ... Procrustes: rmse 3.197905e-05  max resid 8.246851e-05 
    ## ... Similar to previous best
    ## Run 277 stress 0.07428313 
    ## Run 278 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322125  max resid 0.03326268 
    ## Run 279 stress 0.07970525 
    ## Run 280 stress 0.07428313 
    ## Run 281 stress 0.07428317 
    ## Run 282 stress 0.07428312 
    ## Run 283 stress 0.08340286 
    ## Run 284 stress 0.07844934 
    ## Run 285 stress 0.07428313 
    ## Run 286 stress 0.06942776 
    ## ... Procrustes: rmse 3.680904e-06  max resid 9.661107e-06 
    ## ... Similar to previous best
    ## Run 287 stress 0.06978191 
    ## ... Procrustes: rmse 0.0132452  max resid 0.03330436 
    ## Run 288 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326876  max resid 0.03335721 
    ## Run 289 stress 0.08340289 
    ## Run 290 stress 0.07428314 
    ## Run 291 stress 0.07428315 
    ## Run 292 stress 0.08340287 
    ## Run 293 stress 0.07250817 
    ## Run 294 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319745  max resid 0.03320859 
    ## Run 295 stress 0.07844938 
    ## Run 296 stress 0.08448437 
    ## Run 297 stress 0.06942777 
    ## ... Procrustes: rmse 5.446249e-05  max resid 0.0001406597 
    ## ... Similar to previous best
    ## Run 298 stress 0.06942776 
    ## ... Procrustes: rmse 1.476503e-05  max resid 3.809671e-05 
    ## ... Similar to previous best
    ## Run 299 stress 0.08340299 
    ## Run 300 stress 0.07250812 
    ## Run 301 stress 0.08340288 
    ## Run 302 stress 0.08340287 
    ## Run 303 stress 0.07428313 
    ## Run 304 stress 0.08340287 
    ## Run 305 stress 0.08340289 
    ## Run 306 stress 0.07428312 
    ## Run 307 stress 0.06942776 
    ## ... Procrustes: rmse 2.07201e-05  max resid 5.354857e-05 
    ## ... Similar to previous best
    ## Run 308 stress 0.07428314 
    ## Run 309 stress 0.07428313 
    ## Run 310 stress 0.07970526 
    ## Run 311 stress 0.07970526 
    ## Run 312 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319417  max resid 0.03320075 
    ## Run 313 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326184  max resid 0.03334343 
    ## Run 314 stress 0.07250812 
    ## Run 315 stress 0.07250814 
    ## Run 316 stress 0.07250813 
    ## Run 317 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324289  max resid 0.03330682 
    ## Run 318 stress 0.07250812 
    ## Run 319 stress 0.07428314 
    ## Run 320 stress 0.07250812 
    ## Run 321 stress 0.08340287 
    ## Run 322 stress 0.07844929 
    ## Run 323 stress 0.07428315 
    ## Run 324 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323997  max resid 0.033298 
    ## Run 325 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322826  max resid 0.03327552 
    ## Run 326 stress 0.07428321 
    ## Run 327 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327496  max resid 0.03336985 
    ## Run 328 stress 0.08340291 
    ## Run 329 stress 0.07250814 
    ## Run 330 stress 0.07250812 
    ## Run 331 stress 0.07970526 
    ## Run 332 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322777  max resid 0.03327183 
    ## Run 333 stress 0.07250813 
    ## Run 334 stress 0.07970525 
    ## Run 335 stress 0.06942777 
    ## ... Procrustes: rmse 3.78583e-05  max resid 9.829912e-05 
    ## ... Similar to previous best
    ## Run 336 stress 0.07970525 
    ## Run 337 stress 0.07250813 
    ## Run 338 stress 0.06942778 
    ## ... Procrustes: rmse 4.666785e-05  max resid 0.0001207036 
    ## ... Similar to previous best
    ## Run 339 stress 0.06942778 
    ## ... Procrustes: rmse 1.606144e-05  max resid 4.454798e-05 
    ## ... Similar to previous best
    ## Run 340 stress 0.07428314 
    ## Run 341 stress 0.2585105 
    ## Run 342 stress 0.07428314 
    ## Run 343 stress 0.07250813 
    ## Run 344 stress 0.07970526 
    ## Run 345 stress 0.08340296 
    ## Run 346 stress 0.07250812 
    ## Run 347 stress 0.08340288 
    ## Run 348 stress 0.07970525 
    ## Run 349 stress 0.06978189 
    ## ... Procrustes: rmse 0.0132211  max resid 0.03325687 
    ## Run 350 stress 0.07428314 
    ## Run 351 stress 0.07428313 
    ## Run 352 stress 0.08340287 
    ## Run 353 stress 0.07428316 
    ## Run 354 stress 0.07250812 
    ## Run 355 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327205  max resid 0.03336408 
    ## Run 356 stress 0.07844918 
    ## Run 357 stress 0.07970525 
    ## Run 358 stress 0.06942776 
    ## ... Procrustes: rmse 3.980485e-06  max resid 8.253582e-06 
    ## ... Similar to previous best
    ## Run 359 stress 0.07428318 
    ## Run 360 stress 0.07970525 
    ## Run 361 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324317  max resid 0.03330293 
    ## Run 362 stress 0.07428313 
    ## Run 363 stress 0.07844937 
    ## Run 364 stress 0.07428321 
    ## Run 365 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322758  max resid 0.03327198 
    ## Run 366 stress 0.07428314 
    ## Run 367 stress 0.07970525 
    ## Run 368 stress 0.07428313 
    ## Run 369 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327478  max resid 0.03336834 
    ## Run 370 stress 0.07428312 
    ## Run 371 stress 0.06942778 
    ## ... Procrustes: rmse 6.039958e-05  max resid 0.0001558541 
    ## ... Similar to previous best
    ## Run 372 stress 0.07970527 
    ## Run 373 stress 0.07428313 
    ## Run 374 stress 0.08340293 
    ## Run 375 stress 0.07844961 
    ## Run 376 stress 0.07428316 
    ## Run 377 stress 0.07970526 
    ## Run 378 stress 0.08340286 
    ## Run 379 stress 0.06942776 
    ## ... Procrustes: rmse 1.611675e-05  max resid 4.375495e-05 
    ## ... Similar to previous best
    ## Run 380 stress 0.07970525 
    ## Run 381 stress 0.07250822 
    ## Run 382 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327182  max resid 0.03336099 
    ## Run 383 stress 0.07428313 
    ## Run 384 stress 0.06942776 
    ## ... Procrustes: rmse 1.934551e-05  max resid 5.007728e-05 
    ## ... Similar to previous best
    ## Run 385 stress 0.07428314 
    ## Run 386 stress 0.08448438 
    ## Run 387 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326889  max resid 0.03335609 
    ## Run 388 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321494  max resid 0.03324581 
    ## Run 389 stress 0.08340291 
    ## Run 390 stress 0.06942777 
    ## ... Procrustes: rmse 2.983165e-05  max resid 7.832943e-05 
    ## ... Similar to previous best
    ## Run 391 stress 0.07250815 
    ## Run 392 stress 0.07970527 
    ## Run 393 stress 0.07844939 
    ## Run 394 stress 0.07428314 
    ## Run 395 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325202  max resid 0.03332318 
    ## Run 396 stress 0.0742832 
    ## Run 397 stress 0.08340288 
    ## Run 398 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324711  max resid 0.0333112 
    ## Run 399 stress 0.07970525 
    ## Run 400 stress 0.07250813 
    ## Run 401 stress 0.07970525 
    ## Run 402 stress 0.06942776 
    ## ... Procrustes: rmse 2.121159e-05  max resid 5.494303e-05 
    ## ... Similar to previous best
    ## Run 403 stress 0.08448456 
    ## Run 404 stress 0.06942777 
    ## ... Procrustes: rmse 2.052193e-05  max resid 5.56816e-05 
    ## ... Similar to previous best
    ## Run 405 stress 0.06942776 
    ## ... Procrustes: rmse 2.839977e-06  max resid 6.917103e-06 
    ## ... Similar to previous best
    ## Run 406 stress 0.07428317 
    ## Run 407 stress 0.07250812 
    ## Run 408 stress 0.07844919 
    ## Run 409 stress 0.08448447 
    ## Run 410 stress 0.07970526 
    ## Run 411 stress 0.07428312 
    ## Run 412 stress 0.06942777 
    ## ... Procrustes: rmse 2.661866e-05  max resid 7.170321e-05 
    ## ... Similar to previous best
    ## Run 413 stress 0.07844945 
    ## Run 414 stress 0.06942776 
    ## ... Procrustes: rmse 7.370668e-06  max resid 1.955969e-05 
    ## ... Similar to previous best
    ## Run 415 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327346  max resid 0.03336797 
    ## Run 416 stress 0.07428313 
    ## Run 417 stress 0.07844943 
    ## Run 418 stress 0.07428317 
    ## Run 419 stress 0.06942776 
    ## ... Procrustes: rmse 8.484904e-06  max resid 2.242413e-05 
    ## ... Similar to previous best
    ## Run 420 stress 0.07970526 
    ## Run 421 stress 0.07428313 
    ## Run 422 stress 0.0742832 
    ## Run 423 stress 0.07970526 
    ## Run 424 stress 0.07250817 
    ## Run 425 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325536  max resid 0.03333191 
    ## Run 426 stress 0.06942776 
    ## ... Procrustes: rmse 2.861414e-06  max resid 7.639454e-06 
    ## ... Similar to previous best
    ## Run 427 stress 0.06942777 
    ## ... Procrustes: rmse 5.406118e-05  max resid 0.0001391915 
    ## ... Similar to previous best
    ## Run 428 stress 0.07428315 
    ## Run 429 stress 0.07844961 
    ## Run 430 stress 0.08340287 
    ## Run 431 stress 0.08340288 
    ## Run 432 stress 0.07428315 
    ## Run 433 stress 0.08340287 
    ## Run 434 stress 0.0834029 
    ## Run 435 stress 0.07428313 
    ## Run 436 stress 0.07428315 
    ## Run 437 stress 0.06942776 
    ## ... Procrustes: rmse 3.622743e-06  max resid 9.82945e-06 
    ## ... Similar to previous best
    ## Run 438 stress 0.08340296 
    ## Run 439 stress 0.06942776 
    ## ... Procrustes: rmse 1.21475e-05  max resid 3.341078e-05 
    ## ... Similar to previous best
    ## Run 440 stress 0.08340287 
    ## Run 441 stress 0.07428314 
    ## Run 442 stress 0.07250812 
    ## Run 443 stress 0.07250812 
    ## Run 444 stress 0.06942776 
    ## ... Procrustes: rmse 4.53663e-06  max resid 1.229903e-05 
    ## ... Similar to previous best
    ## Run 445 stress 0.07844955 
    ## Run 446 stress 0.07970525 
    ## Run 447 stress 0.07970525 
    ## Run 448 stress 0.07970525 
    ## Run 449 stress 0.07428314 
    ## Run 450 stress 0.07250812 
    ## Run 451 stress 0.08340287 
    ## Run 452 stress 0.07250816 
    ## Run 453 stress 0.06942776 
    ## ... Procrustes: rmse 2.308132e-05  max resid 5.948694e-05 
    ## ... Similar to previous best
    ## Run 454 stress 0.08448451 
    ## Run 455 stress 0.06942777 
    ## ... Procrustes: rmse 4.786758e-05  max resid 0.0001236485 
    ## ... Similar to previous best
    ## Run 456 stress 0.07250812 
    ## Run 457 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326655  max resid 0.03335192 
    ## Run 458 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322854  max resid 0.03327346 
    ## Run 459 stress 0.08340292 
    ## Run 460 stress 0.07428313 
    ## Run 461 stress 0.08340289 
    ## Run 462 stress 0.07250814 
    ## Run 463 stress 0.07428314 
    ## Run 464 stress 0.07428313 
    ## Run 465 stress 0.07428313 
    ## Run 466 stress 0.07970525 
    ## Run 467 stress 0.07970525 
    ## Run 468 stress 0.06978192 
    ## ... Procrustes: rmse 0.01323165  max resid 0.03328407 
    ## Run 469 stress 0.07250813 
    ## Run 470 stress 0.07250812 
    ## Run 471 stress 0.07250817 
    ## Run 472 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323403  max resid 0.0332839 
    ## Run 473 stress 0.07970525 
    ## Run 474 stress 0.06978196 
    ## ... Procrustes: rmse 0.01317374  max resid 0.03315595 
    ## Run 475 stress 0.06942776 
    ## ... Procrustes: rmse 1.393327e-05  max resid 3.601521e-05 
    ## ... Similar to previous best
    ## Run 476 stress 0.07250812 
    ## Run 477 stress 0.07970525 
    ## Run 478 stress 0.07428313 
    ## Run 479 stress 0.07970526 
    ## Run 480 stress 0.07428313 
    ## Run 481 stress 0.07250813 
    ## Run 482 stress 0.07970525 
    ## Run 483 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325891  max resid 0.03334202 
    ## Run 484 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324306  max resid 0.03330472 
    ## Run 485 stress 0.08340289 
    ## Run 486 stress 0.08340287 
    ## Run 487 stress 0.07250814 
    ## Run 488 stress 0.07428313 
    ## Run 489 stress 0.07428313 
    ## Run 490 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319956  max resid 0.03321391 
    ## Run 491 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326848  max resid 0.03335755 
    ## Run 492 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323367  max resid 0.03328601 
    ## Run 493 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322878  max resid 0.03327611 
    ## Run 494 stress 0.08340289 
    ## Run 495 stress 0.07970526 
    ## Run 496 stress 0.06942776 
    ## ... Procrustes: rmse 1.869674e-05  max resid 5.029285e-05 
    ## ... Similar to previous best
    ## Run 497 stress 0.06942778 
    ## ... Procrustes: rmse 5.879831e-05  max resid 0.0001512498 
    ## ... Similar to previous best
    ## Run 498 stress 0.08340288 
    ## Run 499 stress 0.07970525 
    ## Run 500 stress 0.08340287 
    ## *** Best solution repeated 37 times

``` r
# Mixed lakes
SD_beta_geo_M_NMDS <- metaMDS(SD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.001279543 
    ## Run 1 stress 8.09222e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0432457  max resid 0.05943944 
    ## Run 2 stress 0.000445898 
    ## ... Procrustes: rmse 0.01479051  max resid 0.02034606 
    ## Run 3 stress 0.001363604 
    ## Run 4 stress 0.1776015 
    ## Run 5 stress 0.001374625 
    ## Run 6 stress 0.001431988 
    ## Run 7 stress 0.001233888 
    ## Run 8 stress 9.022464e-05 
    ## ... Procrustes: rmse 0.0002244797  max resid 0.0004342569 
    ## ... Similar to previous best
    ## Run 9 stress 0.3120129 
    ## Run 10 stress 0.00124257 
    ## Run 11 stress 0.2528294 
    ## Run 12 stress 0.1491957 
    ## Run 13 stress 8.862028e-05 
    ## ... Procrustes: rmse 0.0002210536  max resid 0.0004384559 
    ## ... Similar to previous best
    ## Run 14 stress 0.0004698505 
    ## ... Procrustes: rmse 0.01578473  max resid 0.02171492 
    ## Run 15 stress 0.1491957 
    ## Run 16 stress 0.001181935 
    ## Run 17 stress 9.531391e-05 
    ## ... Procrustes: rmse 0.0001929004  max resid 0.0003074461 
    ## ... Similar to previous best
    ## Run 18 stress 0.0005189867 
    ## ... Procrustes: rmse 0.01621848  max resid 0.02231288 
    ## Run 19 stress 0.1990774 
    ## Run 20 stress 9.593137e-05 
    ## ... Procrustes: rmse 0.0001318363  max resid 0.0002178063 
    ## ... Similar to previous best
    ## Run 21 stress 0.0005176124 
    ## ... Procrustes: rmse 0.01657546  max resid 0.02280735 
    ## Run 22 stress 0.1491957 
    ## Run 23 stress 9.265785e-05 
    ## ... Procrustes: rmse 0.0001263198  max resid 0.0002413074 
    ## ... Similar to previous best
    ## Run 24 stress 8.741494e-05 
    ## ... Procrustes: rmse 0.000151866  max resid 0.0002785374 
    ## ... Similar to previous best
    ## Run 25 stress 8.784489e-05 
    ## ... Procrustes: rmse 0.0001654502  max resid 0.0002362748 
    ## ... Similar to previous best
    ## Run 26 stress 0.001453385 
    ## Run 27 stress 0.2830608 
    ## Run 28 stress 9.651372e-05 
    ## ... Procrustes: rmse 0.0001648844  max resid 0.0002840905 
    ## ... Similar to previous best
    ## Run 29 stress 0.2848513 
    ## Run 30 stress 8.995027e-05 
    ## ... Procrustes: rmse 0.0001310663  max resid 0.0002603384 
    ## ... Similar to previous best
    ## Run 31 stress 0.0004866703 
    ## ... Procrustes: rmse 0.01606733  max resid 0.02210693 
    ## Run 32 stress 0.1491957 
    ## Run 33 stress 8.876503e-05 
    ## ... Procrustes: rmse 0.0001024626  max resid 0.0001669046 
    ## ... Similar to previous best
    ## Run 34 stress 0.001227536 
    ## Run 35 stress 0.2520602 
    ## Run 36 stress 0.001418459 
    ## Run 37 stress 0.284852 
    ## Run 38 stress 0.001236818 
    ## Run 39 stress 0.001130797 
    ## Run 40 stress 0.284852 
    ## Run 41 stress 0.0004838882 
    ## ... Procrustes: rmse 0.01602289  max resid 0.02204369 
    ## Run 42 stress 0.001082231 
    ## Run 43 stress 8.24549e-05 
    ## ... Procrustes: rmse 0.0001032108  max resid 0.0001842527 
    ## ... Similar to previous best
    ## Run 44 stress 8.267403e-05 
    ## ... Procrustes: rmse 0.0002068157  max resid 0.0004218053 
    ## ... Similar to previous best
    ## Run 45 stress 8.550649e-05 
    ## ... Procrustes: rmse 0.001606959  max resid 0.002170024 
    ## ... Similar to previous best
    ## Run 46 stress 0.0004794593 
    ## ... Procrustes: rmse 0.01594867  max resid 0.02194115 
    ## Run 47 stress 8.856712e-05 
    ## ... Procrustes: rmse 0.0001536852  max resid 0.0002803056 
    ## ... Similar to previous best
    ## Run 48 stress 8.124103e-05 
    ## ... Procrustes: rmse 0.002232534  max resid 0.003014242 
    ## ... Similar to previous best
    ## Run 49 stress 0.001155983 
    ## Run 50 stress 0.000494428 
    ## ... Procrustes: rmse 0.01619725  max resid 0.02228394 
    ## Run 51 stress 7.951855e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001764271  max resid 0.0003873502 
    ## ... Similar to previous best
    ## Run 52 stress 9.057694e-05 
    ## ... Procrustes: rmse 0.0001856  max resid 0.0003766169 
    ## ... Similar to previous best
    ## Run 53 stress 0.2509018 
    ## Run 54 stress 8.928457e-05 
    ## ... Procrustes: rmse 2.252441e-05  max resid 4.209094e-05 
    ## ... Similar to previous best
    ## Run 55 stress 9.46288e-05 
    ## ... Procrustes: rmse 8.99725e-05  max resid 0.0002070687 
    ## ... Similar to previous best
    ## Run 56 stress 9.876664e-05 
    ## ... Procrustes: rmse 0.003524241  max resid 0.004794594 
    ## ... Similar to previous best
    ## Run 57 stress 9.468862e-05 
    ## ... Procrustes: rmse 0.0001223295  max resid 0.0001898877 
    ## ... Similar to previous best
    ## Run 58 stress 0.00136259 
    ## Run 59 stress 0.00124428 
    ## Run 60 stress 8.214931e-05 
    ## ... Procrustes: rmse 0.0001712393  max resid 0.0003478222 
    ## ... Similar to previous best
    ## Run 61 stress 0.001412869 
    ## Run 62 stress 0.001130658 
    ## Run 63 stress 0.001384918 
    ## Run 64 stress 0.001233798 
    ## Run 65 stress 0.1491957 
    ## Run 66 stress 0.001287008 
    ## Run 67 stress 9.100132e-05 
    ## ... Procrustes: rmse 8.135664e-05  max resid 0.0001713315 
    ## ... Similar to previous best
    ## Run 68 stress 9.618994e-05 
    ## ... Procrustes: rmse 0.0001912841  max resid 0.0003865456 
    ## ... Similar to previous best
    ## Run 69 stress 0.2520602 
    ## Run 70 stress 9.333808e-05 
    ## ... Procrustes: rmse 0.0001830835  max resid 0.0003720637 
    ## ... Similar to previous best
    ## Run 71 stress 0.001323997 
    ## Run 72 stress 9.986598e-05 
    ## ... Procrustes: rmse 9.467444e-05  max resid 0.0002094524 
    ## ... Similar to previous best
    ## Run 73 stress 0.0002907433 
    ## ... Procrustes: rmse 0.02053361  max resid 0.02819928 
    ## Run 74 stress 0.0004796799 
    ## ... Procrustes: rmse 0.01601012  max resid 0.02191986 
    ## Run 75 stress 9.520651e-05 
    ## ... Procrustes: rmse 0.0001657192  max resid 0.0003283917 
    ## ... Similar to previous best
    ## Run 76 stress 0.0004842756 
    ## ... Procrustes: rmse 0.01607713  max resid 0.02201163 
    ## Run 77 stress 0.30831 
    ## Run 78 stress 0.001311301 
    ## Run 79 stress 0.2528294 
    ## Run 80 stress 0.001148579 
    ## Run 81 stress 0.2570108 
    ## Run 82 stress 9.94852e-05 
    ## ... Procrustes: rmse 0.0001307434  max resid 0.0002241481 
    ## ... Similar to previous best
    ## Run 83 stress 0.001270194 
    ## Run 84 stress 0.0004421771 
    ## ... Procrustes: rmse 0.01536811  max resid 0.02103304 
    ## Run 85 stress 9.647061e-05 
    ## ... Procrustes: rmse 8.583766e-05  max resid 0.000195961 
    ## ... Similar to previous best
    ## Run 86 stress 0.1491957 
    ## Run 87 stress 0.001335488 
    ## Run 88 stress 9.172295e-05 
    ## ... Procrustes: rmse 8.575979e-05  max resid 0.0002015382 
    ## ... Similar to previous best
    ## Run 89 stress 0.1776012 
    ## Run 90 stress 8.664673e-05 
    ## ... Procrustes: rmse 0.002772243  max resid 0.003750312 
    ## ... Similar to previous best
    ## Run 91 stress 8.955971e-05 
    ## ... Procrustes: rmse 0.0001990451  max resid 0.0003983582 
    ## ... Similar to previous best
    ## Run 92 stress 0.0001292626 
    ## ... Procrustes: rmse 0.008279093  max resid 0.01125092 
    ## Run 93 stress 0.001278736 
    ## Run 94 stress 0.001276589 
    ## Run 95 stress 0.001388535 
    ## Run 96 stress 0.2848527 
    ## Run 97 stress 0.2440272 
    ## Run 98 stress 8.841183e-05 
    ## ... Procrustes: rmse 0.0001984372  max resid 0.0004181593 
    ## ... Similar to previous best
    ## Run 99 stress 0.1990774 
    ## Run 100 stress 9.69396e-05 
    ## ... Procrustes: rmse 0.0001595367  max resid 0.000319056 
    ## ... Similar to previous best
    ## Run 101 stress 9.915641e-05 
    ## ... Procrustes: rmse 0.0002073169  max resid 0.0004168546 
    ## ... Similar to previous best
    ## Run 102 stress 0.0009477758 
    ## Run 103 stress 9.811856e-05 
    ## ... Procrustes: rmse 0.0001169694  max resid 0.0002572259 
    ## ... Similar to previous best
    ## Run 104 stress 0.1491957 
    ## Run 105 stress 0.000918937 
    ## Run 106 stress 9.479012e-05 
    ## ... Procrustes: rmse 9.264671e-05  max resid 0.0001894291 
    ## ... Similar to previous best
    ## Run 107 stress 0.0008445624 
    ## Run 108 stress 0.001362849 
    ## Run 109 stress 0.001324929 
    ## Run 110 stress 0.001406977 
    ## Run 111 stress 9.82127e-05 
    ## ... Procrustes: rmse 0.0001856927  max resid 0.0003769927 
    ## ... Similar to previous best
    ## Run 112 stress 9.701033e-05 
    ## ... Procrustes: rmse 3.874074e-05  max resid 6.109445e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.1776015 
    ## Run 114 stress 0.000472094 
    ## ... Procrustes: rmse 0.01588144  max resid 0.02174281 
    ## Run 115 stress 0.001304462 
    ## Run 116 stress 0.1491957 
    ## Run 117 stress 9.975956e-05 
    ## ... Procrustes: rmse 0.0001302945  max resid 0.0002013586 
    ## ... Similar to previous best
    ## Run 118 stress 0.001158522 
    ## Run 119 stress 8.682922e-05 
    ## ... Procrustes: rmse 0.0001352988  max resid 0.0002184814 
    ## ... Similar to previous best
    ## Run 120 stress 0.001225356 
    ## Run 121 stress 8.294239e-05 
    ## ... Procrustes: rmse 7.738016e-05  max resid 0.0001604441 
    ## ... Similar to previous best
    ## Run 122 stress 0.001325032 
    ## Run 123 stress 0.001070399 
    ## Run 124 stress 0.1491957 
    ## Run 125 stress 0.001026955 
    ## Run 126 stress 9.564016e-05 
    ## ... Procrustes: rmse 0.0001891961  max resid 0.0003940147 
    ## ... Similar to previous best
    ## Run 127 stress 0.1491957 
    ## Run 128 stress 0.1491957 
    ## Run 129 stress 9.330876e-05 
    ## ... Procrustes: rmse 0.0002040595  max resid 0.0004207193 
    ## ... Similar to previous best
    ## Run 130 stress 0.000123703 
    ## ... Procrustes: rmse 0.008085679  max resid 0.0109838 
    ## Run 131 stress 9.842699e-05 
    ## ... Procrustes: rmse 9.099247e-05  max resid 0.0002118124 
    ## ... Similar to previous best
    ## Run 132 stress 0.1491957 
    ## Run 133 stress 0.1491957 
    ## Run 134 stress 8.301148e-05 
    ## ... Procrustes: rmse 0.004109225  max resid 0.005512422 
    ## ... Similar to previous best
    ## Run 135 stress 9.09616e-05 
    ## ... Procrustes: rmse 0.004287267  max resid 0.005762762 
    ## ... Similar to previous best
    ## Run 136 stress 9.330359e-05 
    ## ... Procrustes: rmse 0.000160135  max resid 0.0003478336 
    ## ... Similar to previous best
    ## Run 137 stress 0.0004657162 
    ## ... Procrustes: rmse 0.01577448  max resid 0.02159382 
    ## Run 138 stress 8.713681e-05 
    ## ... Procrustes: rmse 0.002263526  max resid 0.003039852 
    ## ... Similar to previous best
    ## Run 139 stress 0.0009473336 
    ## Run 140 stress 0.001336285 
    ## Run 141 stress 0.0004712602 
    ## ... Procrustes: rmse 0.01586845  max resid 0.02172348 
    ## Run 142 stress 8.588003e-05 
    ## ... Procrustes: rmse 0.0001129435  max resid 0.0001860533 
    ## ... Similar to previous best
    ## Run 143 stress 0.0009567321 
    ## Run 144 stress 0.0004137202 
    ## ... Procrustes: rmse 0.01479812  max resid 0.02024737 
    ## Run 145 stress 8.266422e-05 
    ## ... Procrustes: rmse 7.191697e-05  max resid 0.0001681386 
    ## ... Similar to previous best
    ## Run 146 stress 8.179137e-05 
    ## ... Procrustes: rmse 3.316996e-05  max resid 5.780088e-05 
    ## ... Similar to previous best
    ## Run 147 stress 9.984296e-05 
    ## ... Procrustes: rmse 0.000192253  max resid 0.0004157285 
    ## ... Similar to previous best
    ## Run 148 stress 9.998584e-05 
    ## ... Procrustes: rmse 0.0002097633  max resid 0.0004482204 
    ## ... Similar to previous best
    ## Run 149 stress 8.531615e-05 
    ## ... Procrustes: rmse 1.521832e-05  max resid 2.04632e-05 
    ## ... Similar to previous best
    ## Run 150 stress 9.737739e-05 
    ## ... Procrustes: rmse 9.09833e-05  max resid 0.0001876353 
    ## ... Similar to previous best
    ## Run 151 stress 0.001311196 
    ## Run 152 stress 8.87472e-05 
    ## ... Procrustes: rmse 0.0001124729  max resid 0.0001793073 
    ## ... Similar to previous best
    ## Run 153 stress 0.1491957 
    ## Run 154 stress 0.001159751 
    ## Run 155 stress 9.671197e-05 
    ## ... Procrustes: rmse 0.0002040647  max resid 0.0004072422 
    ## ... Similar to previous best
    ## Run 156 stress 9.155841e-05 
    ## ... Procrustes: rmse 0.0001961075  max resid 0.0003903134 
    ## ... Similar to previous best
    ## Run 157 stress 0.2848521 
    ## Run 158 stress 9.969815e-05 
    ## ... Procrustes: rmse 0.0001906475  max resid 0.0003698006 
    ## ... Similar to previous best
    ## Run 159 stress 0.2570108 
    ## Run 160 stress 0.279732 
    ## Run 161 stress 0.1990774 
    ## Run 162 stress 0.0004688197 
    ## ... Procrustes: rmse 0.01582749  max resid 0.02166765 
    ## Run 163 stress 8.57526e-05 
    ## ... Procrustes: rmse 0.0001449496  max resid 0.0002836385 
    ## ... Similar to previous best
    ## Run 164 stress 9.625154e-05 
    ## ... Procrustes: rmse 0.0001880898  max resid 0.0003656541 
    ## ... Similar to previous best
    ## Run 165 stress 0.1491957 
    ## Run 166 stress 0.2520602 
    ## Run 167 stress 0.2528294 
    ## Run 168 stress 0.0004836157 
    ## ... Procrustes: rmse 0.01606944  max resid 0.02200317 
    ## Run 169 stress 0.0004679329 
    ## ... Procrustes: rmse 0.01581227  max resid 0.02164611 
    ## Run 170 stress 0.0001830794 
    ## ... Procrustes: rmse 0.009865365  max resid 0.01343976 
    ## Run 171 stress 9.592957e-05 
    ## ... Procrustes: rmse 8.851081e-05  max resid 0.0001993932 
    ## ... Similar to previous best
    ## Run 172 stress 0.264697 
    ## Run 173 stress 8.741e-05 
    ## ... Procrustes: rmse 0.0001264832  max resid 0.0002050266 
    ## ... Similar to previous best
    ## Run 174 stress 0.0008844979 
    ## Run 175 stress 9.498442e-05 
    ## ... Procrustes: rmse 0.002740933  max resid 0.003648091 
    ## ... Similar to previous best
    ## Run 176 stress 7.087573e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001021318  max resid 0.0001601013 
    ## ... Similar to previous best
    ## Run 177 stress 0.0004385498 
    ## ... Procrustes: rmse 0.02521889  max resid 0.03478299 
    ## Run 178 stress 9.502123e-05 
    ## ... Procrustes: rmse 0.0001067593  max resid 0.0001615175 
    ## ... Similar to previous best
    ## Run 179 stress 0.001139538 
    ## Run 180 stress 0.2851063 
    ## Run 181 stress 0.2509018 
    ## Run 182 stress 0.001253805 
    ## Run 183 stress 9.915632e-05 
    ## ... Procrustes: rmse 0.0001073257  max resid 0.0001508909 
    ## ... Similar to previous best
    ## Run 184 stress 0.2848521 
    ## Run 185 stress 0.001294345 
    ## Run 186 stress 0.1491957 
    ## Run 187 stress 0.0001018741 
    ## ... Procrustes: rmse 0.007397283  max resid 0.01005785 
    ## Run 188 stress 0.0001060442 
    ## ... Procrustes: rmse 0.007511452  max resid 0.0102148 
    ## Run 189 stress 0.0004059952 
    ## ... Procrustes: rmse 0.01476197  max resid 0.02022083 
    ## Run 190 stress 0.001386984 
    ## Run 191 stress 0.00125891 
    ## Run 192 stress 0.0004662642 
    ## ... Procrustes: rmse 0.01582527  max resid 0.02169232 
    ## Run 193 stress 9.967913e-05 
    ## ... Procrustes: rmse 0.0001081442  max resid 0.0001610984 
    ## ... Similar to previous best
    ## Run 194 stress 8.034862e-05 
    ## ... Procrustes: rmse 0.001527373  max resid 0.002110684 
    ## ... Similar to previous best
    ## Run 195 stress 0.0001763967 
    ## ... Procrustes: rmse 0.01598776  max resid 0.02201167 
    ## Run 196 stress 0.0004513554 
    ## ... Procrustes: rmse 0.0155689  max resid 0.02133457 
    ## Run 197 stress 4.834141e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002060213  max resid 0.000344088 
    ## ... Similar to previous best
    ## Run 198 stress 9.749372e-05 
    ## ... Procrustes: rmse 0.0002473548  max resid 0.0004646056 
    ## ... Similar to previous best
    ## Run 199 stress 0.00131738 
    ## Run 200 stress 0.000499122 
    ## ... Procrustes: rmse 0.01634637  max resid 0.02304494 
    ## Run 201 stress 9.480718e-05 
    ## ... Procrustes: rmse 0.001518519  max resid 0.002564195 
    ## ... Similar to previous best
    ## Run 202 stress 9.163004e-05 
    ## ... Procrustes: rmse 0.0002639375  max resid 0.0005090002 
    ## ... Similar to previous best
    ## Run 203 stress 9.485624e-05 
    ## ... Procrustes: rmse 0.0003113627  max resid 0.0005562994 
    ## ... Similar to previous best
    ## Run 204 stress 9.842211e-05 
    ## ... Procrustes: rmse 0.0001611606  max resid 0.0003173704 
    ## ... Similar to previous best
    ## Run 205 stress 0.3083098 
    ## Run 206 stress 0.1491957 
    ## Run 207 stress 8.795908e-05 
    ## ... Procrustes: rmse 0.0006170673  max resid 0.001250537 
    ## ... Similar to previous best
    ## Run 208 stress 0.0004722904 
    ## ... Procrustes: rmse 0.01589711  max resid 0.02242453 
    ## Run 209 stress 0.00116593 
    ## Run 210 stress 9.981785e-05 
    ## ... Procrustes: rmse 0.0002033918  max resid 0.0004236436 
    ## ... Similar to previous best
    ## Run 211 stress 0.1491957 
    ## Run 212 stress 9.652982e-05 
    ## ... Procrustes: rmse 0.0003143546  max resid 0.0005517795 
    ## ... Similar to previous best
    ## Run 213 stress 0.1491957 
    ## Run 214 stress 0.2842805 
    ## Run 215 stress 0.001320953 
    ## Run 216 stress 0.001408 
    ## Run 217 stress 9.237919e-05 
    ## ... Procrustes: rmse 0.0002409643  max resid 0.0005359462 
    ## ... Similar to previous best
    ## Run 218 stress 0.001307471 
    ## Run 219 stress 0.0009675214 
    ## Run 220 stress 8.403093e-05 
    ## ... Procrustes: rmse 0.0002125917  max resid 0.0003427404 
    ## ... Similar to previous best
    ## Run 221 stress 8.85094e-05 
    ## ... Procrustes: rmse 0.0002200973  max resid 0.0003723639 
    ## ... Similar to previous best
    ## Run 222 stress 0.0004514882 
    ## ... Procrustes: rmse 0.01554247  max resid 0.02193485 
    ## Run 223 stress 9.929227e-05 
    ## ... Procrustes: rmse 0.0002360385  max resid 0.0005356736 
    ## ... Similar to previous best
    ## Run 224 stress 9.260742e-05 
    ## ... Procrustes: rmse 0.0001971151  max resid 0.000392601 
    ## ... Similar to previous best
    ## Run 225 stress 0.00122075 
    ## Run 226 stress 0.1776013 
    ## Run 227 stress 9.250847e-05 
    ## ... Procrustes: rmse 0.0002567841  max resid 0.0004887592 
    ## ... Similar to previous best
    ## Run 228 stress 0.2373687 
    ## Run 229 stress 0.001147057 
    ## Run 230 stress 0.001208949 
    ## Run 231 stress 9.060369e-05 
    ## ... Procrustes: rmse 0.0001222282  max resid 0.000191506 
    ## ... Similar to previous best
    ## Run 232 stress 0.2585028 
    ## Run 233 stress 9.530994e-05 
    ## ... Procrustes: rmse 0.0002372549  max resid 0.0005376295 
    ## ... Similar to previous best
    ## Run 234 stress 0.0005709272 
    ## Run 235 stress 0.0005081702 
    ## ... Procrustes: rmse 0.01649412  max resid 0.02324797 
    ## Run 236 stress 8.160546e-05 
    ## ... Procrustes: rmse 0.0002617344  max resid 0.0004232539 
    ## ... Similar to previous best
    ## Run 237 stress 0.0004857679 
    ## ... Procrustes: rmse 0.01612397  max resid 0.02273773 
    ## Run 238 stress 0.0008135208 
    ## Run 239 stress 0.1491957 
    ## Run 240 stress 0.001112123 
    ## Run 241 stress 0.2842805 
    ## Run 242 stress 0.1491957 
    ## Run 243 stress 8.691534e-05 
    ## ... Procrustes: rmse 0.0002929794  max resid 0.000441426 
    ## ... Similar to previous best
    ## Run 244 stress 0.1491957 
    ## Run 245 stress 9.267725e-05 
    ## ... Procrustes: rmse 0.0003188366  max resid 0.0004574404 
    ## ... Similar to previous best
    ## Run 246 stress 9.049024e-05 
    ## ... Procrustes: rmse 0.002005641  max resid 0.003250829 
    ## ... Similar to previous best
    ## Run 247 stress 0.0004761673 
    ## ... Procrustes: rmse 0.0159511  max resid 0.02249897 
    ## Run 248 stress 9.889953e-05 
    ## ... Procrustes: rmse 0.0002344602  max resid 0.0005366239 
    ## ... Similar to previous best
    ## Run 249 stress 0.1776016 
    ## Run 250 stress 0.0002653821 
    ## ... Procrustes: rmse 0.0119066  max resid 0.01691714 
    ## Run 251 stress 0.001210011 
    ## Run 252 stress 0.001284599 
    ## Run 253 stress 0.001206706 
    ## Run 254 stress 0.0008576869 
    ## Run 255 stress 0.000476012 
    ## ... Procrustes: rmse 0.01595368  max resid 0.02250183 
    ## Run 256 stress 8.811567e-05 
    ## ... Procrustes: rmse 0.001279589  max resid 0.002229474 
    ## ... Similar to previous best
    ## Run 257 stress 0.001257896 
    ## Run 258 stress 9.350085e-05 
    ## ... Procrustes: rmse 0.0002064138  max resid 0.000345003 
    ## ... Similar to previous best
    ## Run 259 stress 0.0003454285 
    ## ... Procrustes: rmse 0.01358428  max resid 0.01923251 
    ## Run 260 stress 0.001321489 
    ## Run 261 stress 9.451754e-05 
    ## ... Procrustes: rmse 0.0002607056  max resid 0.0004948821 
    ## ... Similar to previous best
    ## Run 262 stress 8.910362e-05 
    ## ... Procrustes: rmse 0.0001634473  max resid 0.0003164064 
    ## ... Similar to previous best
    ## Run 263 stress 0.1990774 
    ## Run 264 stress 9.746368e-05 
    ## ... Procrustes: rmse 0.0002431485  max resid 0.0005495099 
    ## ... Similar to previous best
    ## Run 265 stress 9.337924e-05 
    ## ... Procrustes: rmse 0.0002624815  max resid 0.0005011478 
    ## ... Similar to previous best
    ## Run 266 stress 0.001391377 
    ## Run 267 stress 0.001510802 
    ## Run 268 stress 9.845072e-05 
    ## ... Procrustes: rmse 0.0002451845  max resid 0.000539326 
    ## ... Similar to previous best
    ## Run 269 stress 0.001337187 
    ## Run 270 stress 0.1491957 
    ## Run 271 stress 8.795581e-05 
    ## ... Procrustes: rmse 0.0005437352  max resid 0.001128898 
    ## ... Similar to previous best
    ## Run 272 stress 9.131664e-05 
    ## ... Procrustes: rmse 0.0001958277  max resid 0.0003510193 
    ## ... Similar to previous best
    ## Run 273 stress 9.952094e-05 
    ## ... Procrustes: rmse 0.0002018879  max resid 0.0004001762 
    ## ... Similar to previous best
    ## Run 274 stress 0.1491957 
    ## Run 275 stress 8.590095e-05 
    ## ... Procrustes: rmse 0.0002112341  max resid 0.0003576846 
    ## ... Similar to previous best
    ## Run 276 stress 8.987499e-05 
    ## ... Procrustes: rmse 0.0002066627  max resid 0.0003691687 
    ## ... Similar to previous best
    ## Run 277 stress 0.3082762 
    ## Run 278 stress 0.1990774 
    ## Run 279 stress 0.00134688 
    ## Run 280 stress 0.001237406 
    ## Run 281 stress 0.001210706 
    ## Run 282 stress 0.3083098 
    ## Run 283 stress 9.500595e-05 
    ## ... Procrustes: rmse 0.0002590567  max resid 0.000399081 
    ## ... Similar to previous best
    ## Run 284 stress 9.987706e-05 
    ## ... Procrustes: rmse 0.01163391  max resid 0.01605199 
    ## Run 285 stress 0.001379516 
    ## Run 286 stress 9.321015e-05 
    ## ... Procrustes: rmse 0.0002307104  max resid 0.0005191831 
    ## ... Similar to previous best
    ## Run 287 stress 0.2851257 
    ## Run 288 stress 9.945096e-05 
    ## ... Procrustes: rmse 0.0002340966  max resid 0.000438127 
    ## ... Similar to previous best
    ## Run 289 stress 0.2852151 
    ## Run 290 stress 0.001332765 
    ## Run 291 stress 0.001339351 
    ## Run 292 stress 8.654136e-05 
    ## ... Procrustes: rmse 0.0004553673  max resid 0.0009771189 
    ## ... Similar to previous best
    ## Run 293 stress 0.1491957 
    ## Run 294 stress 0.001345581 
    ## Run 295 stress 8.345077e-05 
    ## ... Procrustes: rmse 0.0001672957  max resid 0.000313142 
    ## ... Similar to previous best
    ## Run 296 stress 7.646895e-05 
    ## ... Procrustes: rmse 0.0001701734  max resid 0.0003149973 
    ## ... Similar to previous best
    ## Run 297 stress 0.1491957 
    ## Run 298 stress 9.666194e-05 
    ## ... Procrustes: rmse 0.0003026437  max resid 0.0004601575 
    ## ... Similar to previous best
    ## Run 299 stress 0.0003672023 
    ## ... Procrustes: rmse 0.02290504  max resid 0.03159757 
    ## Run 300 stress 9.839502e-05 
    ## ... Procrustes: rmse 0.003934418  max resid 0.00549587 
    ## ... Similar to previous best
    ## Run 301 stress 0.001235905 
    ## Run 302 stress 8.92022e-05 
    ## ... Procrustes: rmse 0.0002398304  max resid 0.0004837607 
    ## ... Similar to previous best
    ## Run 303 stress 0.001354707 
    ## Run 304 stress 9.653136e-05 
    ## ... Procrustes: rmse 0.0001614758  max resid 0.0003202412 
    ## ... Similar to previous best
    ## Run 305 stress 9.760059e-05 
    ## ... Procrustes: rmse 0.000160533  max resid 0.0003166948 
    ## ... Similar to previous best
    ## Run 306 stress 9.039431e-05 
    ## ... Procrustes: rmse 0.0001596248  max resid 0.0002744713 
    ## ... Similar to previous best
    ## Run 307 stress 0.1491957 
    ## Run 308 stress 0.2795549 
    ## Run 309 stress 9.367287e-05 
    ## ... Procrustes: rmse 0.0002694346  max resid 0.0004124613 
    ## ... Similar to previous best
    ## Run 310 stress 0.1491957 
    ## Run 311 stress 0.001062524 
    ## Run 312 stress 0.001358088 
    ## Run 313 stress 0.1491957 
    ## Run 314 stress 9.805309e-05 
    ## ... Procrustes: rmse 0.0001611009  max resid 0.0003186424 
    ## ... Similar to previous best
    ## Run 315 stress 0.1776012 
    ## Run 316 stress 9.972675e-05 
    ## ... Procrustes: rmse 0.0002915528  max resid 0.0005601393 
    ## ... Similar to previous best
    ## Run 317 stress 0.001197248 
    ## Run 318 stress 0.1491957 
    ## Run 319 stress 0.001346496 
    ## Run 320 stress 0.001226909 
    ## Run 321 stress 0.001319721 
    ## Run 322 stress 8.507461e-05 
    ## ... Procrustes: rmse 0.000371528  max resid 0.0008284665 
    ## ... Similar to previous best
    ## Run 323 stress 9.436606e-05 
    ## ... Procrustes: rmse 0.001439066  max resid 0.002023971 
    ## ... Similar to previous best
    ## Run 324 stress 0.1491957 
    ## Run 325 stress 0.001361315 
    ## Run 326 stress 9.078493e-05 
    ## ... Procrustes: rmse 0.0002327688  max resid 0.0005040577 
    ## ... Similar to previous best
    ## Run 327 stress 0.1491957 
    ## Run 328 stress 8.441226e-05 
    ## ... Procrustes: rmse 0.0002280018  max resid 0.0004729616 
    ## ... Similar to previous best
    ## Run 329 stress 9.368217e-05 
    ## ... Procrustes: rmse 0.000231741  max resid 0.0005186503 
    ## ... Similar to previous best
    ## Run 330 stress 9.716884e-05 
    ## ... Procrustes: rmse 0.0002427014  max resid 0.0005040468 
    ## ... Similar to previous best
    ## Run 331 stress 0.2509018 
    ## Run 332 stress 0.0002262876 
    ## ... Procrustes: rmse 0.010978  max resid 0.01563634 
    ## Run 333 stress 8.454626e-05 
    ## ... Procrustes: rmse 0.0001654936  max resid 0.0003176011 
    ## ... Similar to previous best
    ## Run 334 stress 0.1491957 
    ## Run 335 stress 9.901428e-05 
    ## ... Procrustes: rmse 0.0002492858  max resid 0.0005482029 
    ## ... Similar to previous best
    ## Run 336 stress 9.83702e-05 
    ## ... Procrustes: rmse 0.0002324193  max resid 0.0005295517 
    ## ... Similar to previous best
    ## Run 337 stress 0.1990774 
    ## Run 338 stress 0.001313001 
    ## Run 339 stress 9.636337e-05 
    ## ... Procrustes: rmse 0.01111795  max resid 0.01532396 
    ## Run 340 stress 0.001141133 
    ## Run 341 stress 0.1990774 
    ## Run 342 stress 8.856707e-05 
    ## ... Procrustes: rmse 0.006872346  max resid 0.009471068 
    ## Run 343 stress 0.001249485 
    ## Run 344 stress 0.2646972 
    ## Run 345 stress 0.001329151 
    ## Run 346 stress 9.267139e-05 
    ## ... Procrustes: rmse 0.000308147  max resid 0.0005484738 
    ## ... Similar to previous best
    ## Run 347 stress 0.2696744 
    ## Run 348 stress 9.523791e-05 
    ## ... Procrustes: rmse 0.0002620495  max resid 0.0005039161 
    ## ... Similar to previous best
    ## Run 349 stress 0.001285238 
    ## Run 350 stress 9.053051e-05 
    ## ... Procrustes: rmse 0.000238032  max resid 0.0004760409 
    ## ... Similar to previous best
    ## Run 351 stress 0.2494706 
    ## Run 352 stress 9.123779e-05 
    ## ... Procrustes: rmse 0.0002855325  max resid 0.0005294102 
    ## ... Similar to previous best
    ## Run 353 stress 8.645751e-05 
    ## ... Procrustes: rmse 0.0002627451  max resid 0.0004977898 
    ## ... Similar to previous best
    ## Run 354 stress 0.000479679 
    ## ... Procrustes: rmse 0.01589963  max resid 0.02242921 
    ## Run 355 stress 0.001322845 
    ## Run 356 stress 9.808776e-05 
    ## ... Procrustes: rmse 0.0002002058  max resid 0.0003649371 
    ## ... Similar to previous best
    ## Run 357 stress 0.001266814 
    ## Run 358 stress 9.861378e-05 
    ## ... Procrustes: rmse 0.0002424394  max resid 0.0004988117 
    ## ... Similar to previous best
    ## Run 359 stress 8.444956e-05 
    ## ... Procrustes: rmse 0.001103066  max resid 0.00196747 
    ## ... Similar to previous best
    ## Run 360 stress 0.0004601349 
    ## ... Procrustes: rmse 0.01568828  max resid 0.02214124 
    ## Run 361 stress 0.2373687 
    ## Run 362 stress 0.1990774 
    ## Run 363 stress 0.1776012 
    ## Run 364 stress 0.0009512276 
    ## Run 365 stress 9.600972e-05 
    ## ... Procrustes: rmse 0.0002037784  max resid 0.0004212565 
    ## ... Similar to previous best
    ## Run 366 stress 9.090607e-05 
    ## ... Procrustes: rmse 0.0002049918  max resid 0.0004283307 
    ## ... Similar to previous best
    ## Run 367 stress 0.001255674 
    ## Run 368 stress 9.070645e-05 
    ## ... Procrustes: rmse 0.0002610707  max resid 0.0004904895 
    ## ... Similar to previous best
    ## Run 369 stress 0.1491957 
    ## Run 370 stress 0.0003865858 
    ## ... Procrustes: rmse 0.01438122  max resid 0.020332 
    ## Run 371 stress 0.1776018 
    ## Run 372 stress 0.1491957 
    ## Run 373 stress 0.0003023011 
    ## ... Procrustes: rmse 0.02075818  max resid 0.02862645 
    ## Run 374 stress 0.2848512 
    ## Run 375 stress 0.0004507831 
    ## ... Procrustes: rmse 0.0155327  max resid 0.02192088 
    ## Run 376 stress 0.1990774 
    ## Run 377 stress 0.001449652 
    ## Run 378 stress 9.919654e-05 
    ## ... Procrustes: rmse 0.007153585  max resid 0.01036082 
    ## Run 379 stress 0.0004879943 
    ## ... Procrustes: rmse 0.01616206  max resid 0.02279202 
    ## Run 380 stress 0.2528294 
    ## Run 381 stress 0.0001331138 
    ## ... Procrustes: rmse 0.008420913  max resid 0.01210545 
    ## Run 382 stress 9.402115e-05 
    ## ... Procrustes: rmse 0.0002316303  max resid 0.0005224187 
    ## ... Similar to previous best
    ## Run 383 stress 0.001330311 
    ## Run 384 stress 9.858885e-05 
    ## ... Procrustes: rmse 0.0001620471  max resid 0.0003234043 
    ## ... Similar to previous best
    ## Run 385 stress 0.000501958 
    ## ... Procrustes: rmse 0.01639303  max resid 0.02310854 
    ## Run 386 stress 0.0004616043 
    ## ... Procrustes: rmse 0.01571573  max resid 0.02217404 
    ## Run 387 stress 0.001343474 
    ## Run 388 stress 9.472957e-05 
    ## ... Procrustes: rmse 0.0002778617  max resid 0.0005184085 
    ## ... Similar to previous best
    ## Run 389 stress 0.001248123 
    ## Run 390 stress 0.2528294 
    ## Run 391 stress 9.238207e-05 
    ## ... Procrustes: rmse 0.0001622561  max resid 0.0003199927 
    ## ... Similar to previous best
    ## Run 392 stress 0.001392608 
    ## Run 393 stress 7.832425e-05 
    ## ... Procrustes: rmse 0.000257781  max resid 0.0004828969 
    ## ... Similar to previous best
    ## Run 394 stress 0.0005678459 
    ## Run 395 stress 0.001243537 
    ## Run 396 stress 9.80186e-05 
    ## ... Procrustes: rmse 0.0002655103  max resid 0.0005152459 
    ## ... Similar to previous best
    ## Run 397 stress 0.001296373 
    ## Run 398 stress 8.297729e-05 
    ## ... Procrustes: rmse 0.0004265601  max resid 0.0006281865 
    ## ... Similar to previous best
    ## Run 399 stress 8.787072e-05 
    ## ... Procrustes: rmse 0.0002067251  max resid 0.0003703582 
    ## ... Similar to previous best
    ## Run 400 stress 0.1491957 
    ## Run 401 stress 0.0004957942 
    ## ... Procrustes: rmse 0.01628994  max resid 0.02296663 
    ## Run 402 stress 0.0004206188 
    ## ... Procrustes: rmse 0.01500215  max resid 0.02119059 
    ## Run 403 stress 0.1491957 
    ## Run 404 stress 0.1491957 
    ## Run 405 stress 0.1776018 
    ## Run 406 stress 0.001234551 
    ## Run 407 stress 0.2848517 
    ## Run 408 stress 9.749518e-05 
    ## ... Procrustes: rmse 0.0002443048  max resid 0.0005288631 
    ## ... Similar to previous best
    ## Run 409 stress 9.181635e-05 
    ## ... Procrustes: rmse 0.0002862325  max resid 0.000543791 
    ## ... Similar to previous best
    ## Run 410 stress 0.0006683354 
    ## Run 411 stress 0.001292835 
    ## Run 412 stress 9.419983e-05 
    ## ... Procrustes: rmse 0.0002134514  max resid 0.0004031866 
    ## ... Similar to previous best
    ## Run 413 stress 0.001401271 
    ## Run 414 stress 9.251972e-05 
    ## ... Procrustes: rmse 0.0002070625  max resid 0.0003416164 
    ## ... Similar to previous best
    ## Run 415 stress 0.001327551 
    ## Run 416 stress 8.667864e-05 
    ## ... Procrustes: rmse 0.0002430063  max resid 0.000517418 
    ## ... Similar to previous best
    ## Run 417 stress 0.0006925893 
    ## Run 418 stress 8.998575e-05 
    ## ... Procrustes: rmse 0.0002010446  max resid 0.0003299899 
    ## ... Similar to previous best
    ## Run 419 stress 9.319275e-05 
    ## ... Procrustes: rmse 0.0001616915  max resid 0.0003166575 
    ## ... Similar to previous best
    ## Run 420 stress 9.340122e-05 
    ## ... Procrustes: rmse 0.006689125  max resid 0.009218441 
    ## Run 421 stress 0.0004613208 
    ## ... Procrustes: rmse 0.01569451  max resid 0.022147 
    ## Run 422 stress 8.510999e-05 
    ## ... Procrustes: rmse 0.0002806573  max resid 0.0005254173 
    ## ... Similar to previous best
    ## Run 423 stress 0.0004622469 
    ## ... Procrustes: rmse 0.01572932  max resid 0.02219264 
    ## Run 424 stress 9.009509e-05 
    ## ... Procrustes: rmse 0.0003593771  max resid 0.0007918946 
    ## ... Similar to previous best
    ## Run 425 stress 0.001384596 
    ## Run 426 stress 0.0006640716 
    ## Run 427 stress 9.817036e-05 
    ## ... Procrustes: rmse 0.0002456961  max resid 0.0005363753 
    ## ... Similar to previous best
    ## Run 428 stress 9.927047e-05 
    ## ... Procrustes: rmse 0.0004193082  max resid 0.0009199237 
    ## ... Similar to previous best
    ## Run 429 stress 0.001335401 
    ## Run 430 stress 0.001183935 
    ## Run 431 stress 0.1990774 
    ## Run 432 stress 9.081689e-05 
    ## ... Procrustes: rmse 0.0001759866  max resid 0.000319217 
    ## ... Similar to previous best
    ## Run 433 stress 0.1990774 
    ## Run 434 stress 9.725025e-05 
    ## ... Procrustes: rmse 0.0002437666  max resid 0.0005361476 
    ## ... Similar to previous best
    ## Run 435 stress 9.955835e-05 
    ## ... Procrustes: rmse 0.0002025519  max resid 0.0004096156 
    ## ... Similar to previous best
    ## Run 436 stress 9.817345e-05 
    ## ... Procrustes: rmse 0.0002008974  max resid 0.0003714487 
    ## ... Similar to previous best
    ## Run 437 stress 9.204813e-05 
    ## ... Procrustes: rmse 0.0002852736  max resid 0.0005332542 
    ## ... Similar to previous best
    ## Run 438 stress 0.0004583026 
    ## ... Procrustes: rmse 0.01564  max resid 0.0220713 
    ## Run 439 stress 0.0004938096 
    ## ... Procrustes: rmse 0.01625867  max resid 0.02292302 
    ## Run 440 stress 0.001260114 
    ## Run 441 stress 0.0004575565 
    ## ... Procrustes: rmse 0.01564853  max resid 0.02208118 
    ## Run 442 stress 0.2373687 
    ## Run 443 stress 0.3083098 
    ## Run 444 stress 0.0004869351 
    ## ... Procrustes: rmse 0.01614108  max resid 0.02276141 
    ## Run 445 stress 0.001261667 
    ## Run 446 stress 0.001231535 
    ## Run 447 stress 8.480829e-05 
    ## ... Procrustes: rmse 0.000262845  max resid 0.0004003956 
    ## ... Similar to previous best
    ## Run 448 stress 7.834532e-05 
    ## ... Procrustes: rmse 0.0002278761  max resid 0.0004827657 
    ## ... Similar to previous best
    ## Run 449 stress 9.872087e-05 
    ## ... Procrustes: rmse 0.0002029214  max resid 0.0004210833 
    ## ... Similar to previous best
    ## Run 450 stress 0.001386621 
    ## Run 451 stress 0.0004635202 
    ## ... Procrustes: rmse 0.01574517  max resid 0.02221798 
    ## Run 452 stress 8.264996e-05 
    ## ... Procrustes: rmse 0.0009994225  max resid 0.001826721 
    ## ... Similar to previous best
    ## Run 453 stress 0.1491957 
    ## Run 454 stress 9.246247e-05 
    ## ... Procrustes: rmse 0.0002081218  max resid 0.0003715335 
    ## ... Similar to previous best
    ## Run 455 stress 0.000477284 
    ## ... Procrustes: rmse 0.01598129  max resid 0.02254119 
    ## Run 456 stress 0.1491957 
    ## Run 457 stress 0.001384926 
    ## Run 458 stress 8.572162e-05 
    ## ... Procrustes: rmse 0.0002448195  max resid 0.0004980477 
    ## ... Similar to previous best
    ## Run 459 stress 0.0003993282 
    ## ... Procrustes: rmse 0.02382964  max resid 0.0328771 
    ## Run 460 stress 9.78636e-05 
    ## ... Procrustes: rmse 0.000203111  max resid 0.0004070457 
    ## ... Similar to previous best
    ## Run 461 stress 9.616443e-05 
    ## ... Procrustes: rmse 0.0001996782  max resid 0.0004176994 
    ## ... Similar to previous best
    ## Run 462 stress 0.001045741 
    ## Run 463 stress 9.633386e-05 
    ## ... Procrustes: rmse 0.0002067177  max resid 0.0003800902 
    ## ... Similar to previous best
    ## Run 464 stress 9.427673e-05 
    ## ... Procrustes: rmse 0.0002407054  max resid 0.0005398826 
    ## ... Similar to previous best
    ## Run 465 stress 0.2797322 
    ## Run 466 stress 8.823584e-05 
    ## ... Procrustes: rmse 0.0002288152  max resid 0.0004890884 
    ## ... Similar to previous best
    ## Run 467 stress 0.001427905 
    ## Run 468 stress 0.001249589 
    ## Run 469 stress 8.509493e-05 
    ## ... Procrustes: rmse 0.000280638  max resid 0.0005147691 
    ## ... Similar to previous best
    ## Run 470 stress 0.1491957 
    ## Run 471 stress 0.0002653543 
    ## ... Procrustes: rmse 0.01190668  max resid 0.01691723 
    ## Run 472 stress 0.2854403 
    ## Run 473 stress 0.2829397 
    ## Run 474 stress 8.744124e-05 
    ## ... Procrustes: rmse 0.0006522993  max resid 0.001307365 
    ## ... Similar to previous best
    ## Run 475 stress 0.2842598 
    ## Run 476 stress 8.94863e-05 
    ## ... Procrustes: rmse 0.0002757285  max resid 0.00046126 
    ## ... Similar to previous best
    ## Run 477 stress 9.20764e-05 
    ## ... Procrustes: rmse 0.0002078593  max resid 0.0004151649 
    ## ... Similar to previous best
    ## Run 478 stress 0.0004820633 
    ## ... Procrustes: rmse 0.01606343  max resid 0.02265414 
    ## Run 479 stress 9.707033e-05 
    ## ... Procrustes: rmse 0.002905296  max resid 0.003992917 
    ## ... Similar to previous best
    ## Run 480 stress 0.0004980965 
    ## ... Procrustes: rmse 0.01632753  max resid 0.02301808 
    ## Run 481 stress 0.3075141 
    ## Run 482 stress 8.769123e-05 
    ## ... Procrustes: rmse 0.003668825  max resid 0.005095378 
    ## ... Similar to previous best
    ## Run 483 stress 0.001060742 
    ## Run 484 stress 0.0008415447 
    ## Run 485 stress 0.001171535 
    ## Run 486 stress 0.0003769183 
    ## ... Procrustes: rmse 0.02320955  max resid 0.03201917 
    ## Run 487 stress 0.001178799 
    ## Run 488 stress 0.001218038 
    ## Run 489 stress 0.00132105 
    ## Run 490 stress 0.0004631797 
    ## ... Procrustes: rmse 0.01574545  max resid 0.0222149 
    ## Run 491 stress 0.1776014 
    ## Run 492 stress 0.001306067 
    ## Run 493 stress 0.0004512547 
    ## ... Procrustes: rmse 0.01553925  max resid 0.02193027 
    ## Run 494 stress 9.401507e-05 
    ## ... Procrustes: rmse 0.0002863145  max resid 0.0005388817 
    ## ... Similar to previous best
    ## Run 495 stress 8.791597e-05 
    ## ... Procrustes: rmse 0.0002743326  max resid 0.0004636581 
    ## ... Similar to previous best
    ## Run 496 stress 9.46161e-05 
    ## ... Procrustes: rmse 0.0006370896  max resid 0.0008988681 
    ## ... Similar to previous best
    ## Run 497 stress 9.576855e-05 
    ## ... Procrustes: rmse 0.0002897689  max resid 0.0004811209 
    ## ... Similar to previous best
    ## Run 498 stress 0.2830608 
    ## Run 499 stress 0.001209403 
    ## Run 500 stress 0.001443626 
    ## *** Best solution repeated 114 times

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
    ##                                    NMDS1     NMDS2     r2 Pr(>r)  
    ## env[surveyed_sites_env, c(34)]  0.997820 -0.065952 0.7388  0.031 *
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
    ## env[surveyed_sites_env, c(34)]  0.997820 -0.065952 0.7388  0.031 *
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
    ## temperature_median -0.48783 -0.87294 0.0586  0.821  
    ## salinity_median     0.99782 -0.06595 0.7388  0.029 *
    ## oxygen_median       0.79316  0.60902 0.5362  0.828  
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
    ## temperature_median -0.48783 -0.87294 0.0586  1.000  
    ## salinity_median     0.99782 -0.06595 0.7388  0.087 .
    ## oxygen_median       0.79316  0.60902 0.5362  1.000  
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
    ## temperature_median -0.56529 -0.82489 0.0778  0.802  
    ## salinity_median     0.99173  0.12834 0.7458  0.040 *
    ## oxygen_median       0.94271  0.33362 0.3552  0.903  
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
    ## temperature_median -0.56529 -0.82489 0.0778   1.00
    ## salinity_median     0.99173  0.12834 0.7458   0.12
    ## oxygen_median       0.94271  0.33362 0.3552   1.00
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
    ## temperature_median -0.10614  0.99435 0.0705  0.680   
    ## salinity_median    -0.78732 -0.61655 0.7842  0.008 **
    ## oxygen_median      -0.96712  0.25433 0.7188  0.043 * 
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
    ## salinity_median    -0.78732 -0.61655 0.7842  0.024 *
    ## oxygen_median      -0.96712  0.25433 0.7188  0.129  
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
    ## temperature_median  0.0001598 -1.0000000 0.0466  0.705
    ## salinity_median    -0.0003527  1.0000000 0.5383  0.376
    ## oxygen_median      -0.0033738 -0.9999900 0.7122  0.889
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
    ## temperature_median  0.0001598 -1.0000000 0.0466      1
    ## salinity_median    -0.0003527  1.0000000 0.5383      1
    ## oxygen_median      -0.0033738 -0.9999900 0.7122      1
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
    ##                         NMDS1      NMDS2     r2 Pr(>r)  
    ## temperature_median -0.0021693  1.0000000 0.0350  0.878  
    ## salinity_median     0.0024111 -1.0000000 0.6084  0.067 .
    ## oxygen_median       0.0078681 -0.9999700 0.6140  0.086 .
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
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## temperature_median -0.0021693  1.0000000 0.0350  1.000
    ## salinity_median     0.0024111 -1.0000000 0.6084  0.201
    ## oxygen_median       0.0078681 -0.9999700 0.6140  0.258
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
    ## temperature_median       -0.35881 -0.93341 0.1964  0.607
    ## salinity_median          -0.70036  0.71379 0.5991  0.127
    ## oxygen_median             0.98543 -0.17006 0.5185  0.181
    ## distance_to_ocean_mean_m  0.01326 -0.99991 0.1726  0.575
    ## max_depth                -0.49998 -0.86604 0.0383  0.888
    ## logArea                  -0.26513 -0.96421 0.2144  0.532
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
    ## salinity_median          -0.70036  0.71379 0.5991  0.762
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
    ## distance_to_ocean_min_m -0.99979  0.02043 0.6139  0.151
    ## max_depth               -0.17728 -0.98416 0.1185  0.776
    ## logArea                  0.26569 -0.96406 0.2080  0.118
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
    ## distance_to_ocean_min_m -0.99979  0.02043 0.6139  0.453
    ## max_depth               -0.17728 -0.98416 0.1185  1.000
    ## logArea                  0.26569 -0.96406 0.2080  0.354
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
    ## distance_to_ocean_min_m -0.99979  0.02043 0.6139  0.122
    ## max_depth               -0.17728 -0.98416 0.1185  0.767
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
    ## distance_to_ocean_min_m -0.99979  0.02043 0.6139  0.366
    ## max_depth               -0.17728 -0.98416 0.1185  1.000
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
    ## distance_to_ocean_min_m -0.94032  0.34028 0.5170  0.151
    ## max_depth               -0.21482 -0.97665 0.1851  0.654
    ## logArea                 -0.06539  0.99786 0.0118  0.951
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
    ## distance_to_ocean_min_m -0.94032  0.34028 0.5170  0.453
    ## max_depth               -0.21482 -0.97665 0.1851  1.000
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
    ## distance_to_ocean_min_m  0.86176 -0.50732 0.3198  0.076 .
    ## max_depth               -0.68923 -0.72454 0.5689  0.021 *
    ## logArea                 -0.84908  0.52826 0.2993  0.235  
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
    ## distance_to_ocean_min_m  0.86176 -0.50732 0.3198  0.228  
    ## max_depth               -0.68923 -0.72454 0.5689  0.063 .
    ## logArea                 -0.84908  0.52826 0.2993  0.705  
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
    ## distance_to_ocean_min_m  0.86963  0.49371 0.6793  0.199
    ## max_depth                0.68534  0.72822 0.0510  0.968
    ## logArea                 -0.52042  0.85391 0.2187  0.553
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
    ## distance_to_ocean_min_m  0.86963  0.49371 0.6793  0.597
    ## max_depth                0.68534  0.72822 0.0510  1.000
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
    ## distance_to_ocean_min_m -0.0007238  1.0000000 0.6650  0.065 . 
    ## max_depth                0.0018532 -1.0000000 0.8251  0.014 * 
    ## logArea                  0.0015105  1.0000000 0.8209  0.005 **
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
    ## distance_to_ocean_min_m -0.0007238  1.0000000 0.6650  0.195  
    ## max_depth                0.0018532 -1.0000000 0.8251  0.042 *
    ## logArea                  0.0015105  1.0000000 0.8209  0.015 *
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
    ##       Significance: 0.252 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.290 0.318 0.335 0.351 
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
    ##       Significance: 0.002 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.574 0.602 0.618 0.640 
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
    ##       Significance: 0.493 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.585 0.617 0.649 0.690 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_mant_pv <- rbind(SD_beta_env_mant_t$signif, SD_beta_env_mant_s$signif, SD_beta_env_mant_o$signif)
SD_beta_env_mant_pv <- SD_beta_env_mant_pv[,1]
(SD_beta_env_mant_pv <- p.adjust(SD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.756 0.006 1.000

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
    ##       Significance: 0.318 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.316 0.353 0.387 0.438 
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
    ##       Significance: 0.011 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.556 0.601 0.647 0.668 
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
    ##       Significance: 0.596 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.486 0.532 0.563 0.599 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_MS_mant_pv <- rbind(SD_beta_env_MS_mant_t$signif, SD_beta_env_MS_mant_s$signif, SD_beta_env_MS_mant_o$signif)
SD_beta_env_MS_mant_pv <- SD_beta_env_MS_mant_pv[,1]
(SD_beta_env_MS_mant_pv <- p.adjust(SD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 0.954 0.033 1.000

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
    ##       Significance: 0.772 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.137 0.187 0.223 0.370 
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
    ##       Significance: 0.044 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.352 0.410 0.477 0.545 
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
    ##       Significance: 0.042 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.368 0.503 0.587 0.678 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_OM_mant_pv <- rbind(SD_beta_env_OM_mant_t$signif, SD_beta_env_OM_mant_s$signif, SD_beta_env_OM_mant_o$signif)
SD_beta_env_OM_mant_pv <- SD_beta_env_OM_mant_pv[,1]
(SD_beta_env_OM_mant_pv <- p.adjust(SD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.132 0.126

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
    ##       Significance: 0.209 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0262 0.0507 0.0708 0.0824 
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
    ##       Significance: 0.009 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.456 0.493 0.524 0.547 
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
    ##       Significance: 0.398 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.709 0.736 0.750 0.767 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_SO_mant_pv <- rbind(SD_beta_env_SO_mant_t$signif, SD_beta_env_SO_mant_s$signif, SD_beta_env_SO_mant_o$signif)
SD_beta_env_SO_mant_pv <- SD_beta_env_SO_mant_pv[,1]
(SD_beta_env_SO_mant_pv <- p.adjust(SD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 0.627 0.027 1.000

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
    ##       Significance: 0.75 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.266 0.357 0.447 0.717 
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
    ##       Significance: 0.168 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.270 0.333 0.370 0.433 
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
    ##       Significance: 0.093 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.249 0.400 0.498 0.567 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_M_mant_pv <- rbind(SD_beta_env_M_mant_t$signif, SD_beta_env_M_mant_s$signif, SD_beta_env_M_mant_o$signif)
SD_beta_env_M_mant_pv <- SD_beta_env_M_mant_pv[,1]
(SD_beta_env_M_mant_pv <- p.adjust(SD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.504 0.279

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
    ##       Significance: 0.245 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.230 0.286 0.320 0.385 
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
    ##       Significance: 0.481 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.541 0.576 0.603 0.624 
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
    ##       Significance: 0.622 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.189 0.223 0.250 0.280 
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
    ##       Significance: 0.385 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0800 0.0923 0.1038 0.1083 
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
    ##       Significance: 0.399 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.435 0.480 0.522 0.559 
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
    ##       Significance: 0.759 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.227 0.260 0.297 0.328 
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
    ##       Significance: 0.762 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0651 0.0945 0.1152 0.1378 
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
    ##       Significance: 0.591 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.244 0.271 0.318 0.351 
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
    ##       Significance: 0.062 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.184 0.254 0.310 0.375 
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
    ##       Significance: 0.096 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.199 0.234 0.262 0.291 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_OM_mant_pv <- rbind( SD_beta_geo_OM_mant_dmean$signif, SD_beta_geo_OM_mant_md$signif, SD_beta_geo_OM_mant_la$signif)
SD_beta_geo_OM_mant_pv <- SD_beta_geo_OM_mant_pv[,1]
(SD_beta_geo_OM_mant_pv <- p.adjust(SD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.186 0.288

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
    ##       Significance: 0.497 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.677 0.702 0.732 0.762 
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
    ##       Significance: 0.61 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.130 0.163 0.192 0.226 
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
    ##       Significance: 0.22 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.276 0.293 0.310 0.324 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_SO_mant_pv <- rbind( SD_beta_geo_SO_mant_dmean$signif, SD_beta_geo_SO_mant_md$signif, SD_beta_geo_SO_mant_la$signif)
SD_beta_geo_SO_mant_pv <- SD_beta_geo_SO_mant_pv[,1]
(SD_beta_geo_SO_mant_pv <- p.adjust(SD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.00 1.00 0.66

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
    ##       Significance: 0.465 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.251 0.386 0.508 0.727 
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
    ## 0.303 0.449 0.545 0.601 
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
    ##       Significance: 0.037 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.221 0.374 0.473 0.542 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_M_mant_pv <- rbind( SD_beta_geo_M_mant_dmean$signif, SD_beta_geo_M_mant_md$signif, SD_beta_geo_M_mant_la$signif)
SD_beta_geo_M_mant_pv <- SD_beta_geo_M_mant_pv[,1]
(SD_beta_geo_M_mant_pv <- p.adjust(SD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.006 0.111

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
        legend.position = "right", 
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
