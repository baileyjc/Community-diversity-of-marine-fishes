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
    ## 3 Stratified vs Reference  1 0.6261283 1.909357 0.2143091   0.102      0.612
    ## 4          Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.041      0.246
    ## 5      Mixed vs Reference  1 0.5979203 1.966332 0.2193017   0.116      0.696
    ## 6      Ocean vs Reference  1 0.5599217 1.628213 0.2456489   0.147      0.882
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
    ## 2 Stratified vs Ocean  1 1.3093089 3.913330 0.2459152   0.002      0.006   *
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.063      0.189

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
    ## 3      Mixed vs Ocean  1 0.6396174 2.135322 0.1625633   0.005      0.015   .

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
    ## 2 Stratified vs Ocean  1 1.3260066 4.147587 0.2931657   0.005      0.015   .
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.061      0.183

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
    ## Run 1 stress 0.1079381 
    ## Run 2 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214737  max resid 0.07898057 
    ## Run 3 stress 0.1135598 
    ## Run 4 stress 0.1135597 
    ## Run 5 stress 0.1160544 
    ## Run 6 stress 0.1135598 
    ## Run 7 stress 0.1097989 
    ## Run 8 stress 0.1088547 
    ## Run 9 stress 0.1050338 
    ## ... Procrustes: rmse 0.008884858  max resid 0.02488359 
    ## Run 10 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307313  max resid 0.07848374 
    ## Run 11 stress 0.1083993 
    ## Run 12 stress 0.1063128 
    ## Run 13 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214813  max resid 0.07898466 
    ## Run 14 stress 0.1083998 
    ## Run 15 stress 0.1138986 
    ## Run 16 stress 0.1050144 
    ## ... Procrustes: rmse 0.006662754  max resid 0.02374762 
    ## Run 17 stress 0.113594 
    ## Run 18 stress 0.1083501 
    ## Run 19 stress 0.1079381 
    ## Run 20 stress 0.1050338 
    ## ... Procrustes: rmse 0.008865351  max resid 0.02481863 
    ## Run 21 stress 0.1174537 
    ## Run 22 stress 0.1063128 
    ## Run 23 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214871  max resid 0.07899245 
    ## Run 24 stress 0.1097988 
    ## Run 25 stress 0.1083499 
    ## Run 26 stress 0.1079382 
    ## Run 27 stress 0.1050144 
    ## ... Procrustes: rmse 0.006666977  max resid 0.02377699 
    ## Run 28 stress 0.1097987 
    ## Run 29 stress 0.1083996 
    ## Run 30 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214874  max resid 0.0789905 
    ## Run 31 stress 0.1063128 
    ## Run 32 stress 0.1088548 
    ## Run 33 stress 0.1138395 
    ## Run 34 stress 0.113594 
    ## Run 35 stress 0.113594 
    ## Run 36 stress 0.1083999 
    ## Run 37 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214983  max resid 0.07899561 
    ## Run 38 stress 0.1049903 
    ## ... New best solution
    ## ... Procrustes: rmse 7.386947e-05  max resid 0.0001990407 
    ## ... Similar to previous best
    ## Run 39 stress 0.1063128 
    ## Run 40 stress 0.1136646 
    ## Run 41 stress 0.1135597 
    ## Run 42 stress 0.1088547 
    ## Run 43 stress 0.1051403 
    ## ... Procrustes: rmse 0.02307539  max resid 0.0784787 
    ## Run 44 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307834  max resid 0.0784895 
    ## Run 45 stress 0.1083991 
    ## Run 46 stress 0.1049903 
    ## ... Procrustes: rmse 2.682219e-05  max resid 8.714319e-05 
    ## ... Similar to previous best
    ## Run 47 stress 0.1309385 
    ## Run 48 stress 0.1083996 
    ## Run 49 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307936  max resid 0.07849315 
    ## Run 50 stress 0.1051923 
    ## ... Procrustes: rmse 0.0221462  max resid 0.07900739 
    ## Run 51 stress 0.1082722 
    ## Run 52 stress 0.1097988 
    ## Run 53 stress 0.1058423 
    ## Run 54 stress 0.1160548 
    ## Run 55 stress 0.1051923 
    ## ... Procrustes: rmse 0.02215487  max resid 0.07904889 
    ## Run 56 stress 0.1084965 
    ## Run 57 stress 0.1097988 
    ## Run 58 stress 0.1050338 
    ## ... Procrustes: rmse 0.008897233  max resid 0.02503245 
    ## Run 59 stress 0.1058705 
    ## Run 60 stress 0.1135597 
    ## Run 61 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307896  max resid 0.07849014 
    ## Run 62 stress 0.1083997 
    ## Run 63 stress 0.1050144 
    ## ... Procrustes: rmse 0.006655759  max resid 0.02375386 
    ## Run 64 stress 0.1084964 
    ## Run 65 stress 0.1135597 
    ## Run 66 stress 0.1063128 
    ## Run 67 stress 0.1051402 
    ## ... Procrustes: rmse 0.0230792  max resid 0.0784921 
    ## Run 68 stress 0.1051402 
    ## ... Procrustes: rmse 0.0230785  max resid 0.0784892 
    ## Run 69 stress 0.1136646 
    ## Run 70 stress 0.1160537 
    ## Run 71 stress 0.1102721 
    ## Run 72 stress 0.1063128 
    ## Run 73 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214555  max resid 0.07900125 
    ## Run 74 stress 0.1051923 
    ## ... Procrustes: rmse 0.0221502  max resid 0.07903365 
    ## Run 75 stress 0.1063128 
    ## Run 76 stress 0.1079382 
    ## Run 77 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214465  max resid 0.07899135 
    ## Run 78 stress 0.1063128 
    ## Run 79 stress 0.1136647 
    ## Run 80 stress 0.1051402 
    ## ... Procrustes: rmse 0.02308363  max resid 0.07850891 
    ## Run 81 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214228  max resid 0.0789781 
    ## Run 82 stress 0.1135597 
    ## Run 83 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214499  max resid 0.07899893 
    ## Run 84 stress 0.1051924 
    ## ... Procrustes: rmse 0.02213976  max resid 0.07895841 
    ## Run 85 stress 0.1079383 
    ## Run 86 stress 0.1088547 
    ## Run 87 stress 0.1050144 
    ## ... Procrustes: rmse 0.006656779  max resid 0.02375754 
    ## Run 88 stress 0.1135598 
    ## Run 89 stress 0.1051402 
    ## ... Procrustes: rmse 0.02308021  max resid 0.07847848 
    ## Run 90 stress 0.1058423 
    ## Run 91 stress 0.1049903 
    ## ... Procrustes: rmse 7.010274e-05  max resid 0.0001774912 
    ## ... Similar to previous best
    ## Run 92 stress 0.1136646 
    ## Run 93 stress 0.1135598 
    ## Run 94 stress 0.1135597 
    ## Run 95 stress 0.1063128 
    ## Run 96 stress 0.1136646 
    ## Run 97 stress 0.1063128 
    ## Run 98 stress 0.1084965 
    ## Run 99 stress 0.1084964 
    ## Run 100 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307978  max resid 0.07845895 
    ## Run 101 stress 0.1079382 
    ## Run 102 stress 0.1135597 
    ## Run 103 stress 0.113594 
    ## Run 104 stress 0.1084966 
    ## Run 105 stress 0.1058708 
    ## Run 106 stress 0.1050338 
    ## ... Procrustes: rmse 0.00889812  max resid 0.02503241 
    ## Run 107 stress 0.1088547 
    ## Run 108 stress 0.1083501 
    ## Run 109 stress 0.113594 
    ## Run 110 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307861  max resid 0.07848997 
    ## Run 111 stress 0.1084965 
    ## Run 112 stress 0.1082722 
    ## Run 113 stress 0.1051923 
    ## ... Procrustes: rmse 0.02215052  max resid 0.07901401 
    ## Run 114 stress 0.1079381 
    ## Run 115 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214413  max resid 0.0789916 
    ## Run 116 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307854  max resid 0.07848981 
    ## Run 117 stress 0.110272 
    ## Run 118 stress 0.1082722 
    ## Run 119 stress 0.1083501 
    ## Run 120 stress 0.1097988 
    ## Run 121 stress 0.1058423 
    ## Run 122 stress 0.1136646 
    ## Run 123 stress 0.1084966 
    ## Run 124 stress 0.1136646 
    ## Run 125 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214294  max resid 0.07898291 
    ## Run 126 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214603  max resid 0.07900562 
    ## Run 127 stress 0.1058423 
    ## Run 128 stress 0.1135596 
    ## Run 129 stress 0.1063128 
    ## Run 130 stress 0.1083501 
    ## Run 131 stress 0.1049903 
    ## ... New best solution
    ## ... Procrustes: rmse 2.636634e-05  max resid 6.271474e-05 
    ## ... Similar to previous best
    ## Run 132 stress 0.1049903 
    ## ... Procrustes: rmse 6.832608e-05  max resid 0.0001700157 
    ## ... Similar to previous best
    ## Run 133 stress 0.1063128 
    ## Run 134 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214559  max resid 0.07898551 
    ## Run 135 stress 0.1079382 
    ## Run 136 stress 0.107938 
    ## Run 137 stress 0.1083995 
    ## Run 138 stress 0.1050338 
    ## ... Procrustes: rmse 0.008897508  max resid 0.02499437 
    ## Run 139 stress 0.107938 
    ## Run 140 stress 0.10835 
    ## Run 141 stress 0.1160548 
    ## Run 142 stress 0.1063128 
    ## Run 143 stress 0.1160538 
    ## Run 144 stress 0.1063128 
    ## Run 145 stress 0.1084964 
    ## Run 146 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307436  max resid 0.07847651 
    ## Run 147 stress 0.1084965 
    ## Run 148 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307661  max resid 0.0784844 
    ## Run 149 stress 0.1051402 
    ## ... Procrustes: rmse 0.0230819  max resid 0.07849964 
    ## Run 150 stress 0.1160549 
    ## Run 151 stress 0.10835 
    ## Run 152 stress 0.1097989 
    ## Run 153 stress 0.1079381 
    ## Run 154 stress 0.1063128 
    ## Run 155 stress 0.1049903 
    ## ... Procrustes: rmse 9.627325e-05  max resid 0.0002657657 
    ## ... Similar to previous best
    ## Run 156 stress 0.1050144 
    ## ... Procrustes: rmse 0.006655951  max resid 0.02374772 
    ## Run 157 stress 0.1051402 
    ## ... Procrustes: rmse 0.02308081  max resid 0.07849854 
    ## Run 158 stress 0.1135598 
    ## Run 159 stress 0.1063128 
    ## Run 160 stress 0.1135596 
    ## Run 161 stress 0.1102721 
    ## Run 162 stress 0.1088547 
    ## Run 163 stress 0.1051923 
    ## ... Procrustes: rmse 0.0221448  max resid 0.07899303 
    ## Run 164 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214789  max resid 0.07900157 
    ## Run 165 stress 0.1063128 
    ## Run 166 stress 0.1050144 
    ## ... Procrustes: rmse 0.006657127  max resid 0.02375323 
    ## Run 167 stress 0.1174537 
    ## Run 168 stress 0.1051924 
    ## ... Procrustes: rmse 0.0221404  max resid 0.07894882 
    ## Run 169 stress 0.113594 
    ## Run 170 stress 0.1097988 
    ## Run 171 stress 0.1136647 
    ## Run 172 stress 0.1049903 
    ## ... Procrustes: rmse 6.302789e-05  max resid 0.0001677583 
    ## ... Similar to previous best
    ## Run 173 stress 0.1084965 
    ## Run 174 stress 0.1083502 
    ## Run 175 stress 0.1063129 
    ## Run 176 stress 0.1063128 
    ## Run 177 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307776  max resid 0.07849386 
    ## Run 178 stress 0.1050144 
    ## ... Procrustes: rmse 0.006660039  max resid 0.02374059 
    ## Run 179 stress 0.1097987 
    ## Run 180 stress 0.1097988 
    ## Run 181 stress 0.1049903 
    ## ... New best solution
    ## ... Procrustes: rmse 4.900617e-06  max resid 9.187916e-06 
    ## ... Similar to previous best
    ## Run 182 stress 0.1160538 
    ## Run 183 stress 0.1079381 
    ## Run 184 stress 0.1084964 
    ## Run 185 stress 0.1050144 
    ## ... Procrustes: rmse 0.006658564  max resid 0.02375508 
    ## Run 186 stress 0.1136647 
    ## Run 187 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214289  max resid 0.07896903 
    ## Run 188 stress 0.1084964 
    ## Run 189 stress 0.1135599 
    ## Run 190 stress 0.1160572 
    ## Run 191 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214603  max resid 0.07900138 
    ## Run 192 stress 0.1049903 
    ## ... Procrustes: rmse 5.200901e-05  max resid 0.0001352612 
    ## ... Similar to previous best
    ## Run 193 stress 0.1063128 
    ## Run 194 stress 0.116055 
    ## Run 195 stress 0.1063128 
    ## Run 196 stress 0.1084965 
    ## Run 197 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307451  max resid 0.07847927 
    ## Run 198 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307779  max resid 0.07849393 
    ## Run 199 stress 0.1136646 
    ## Run 200 stress 0.1050144 
    ## ... Procrustes: rmse 0.006662577  max resid 0.02378868 
    ## Run 201 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307963  max resid 0.07849658 
    ## Run 202 stress 0.1083499 
    ## Run 203 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307573  max resid 0.07848776 
    ## Run 204 stress 0.1050338 
    ## ... Procrustes: rmse 0.008898668  max resid 0.02499038 
    ## Run 205 stress 0.10835 
    ## Run 206 stress 0.1135598 
    ## Run 207 stress 0.1049903 
    ## ... Procrustes: rmse 5.024045e-06  max resid 1.416704e-05 
    ## ... Similar to previous best
    ## Run 208 stress 0.1051923 
    ## ... Procrustes: rmse 0.0221437  max resid 0.07897765 
    ## Run 209 stress 0.1135598 
    ## Run 210 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214628  max resid 0.07899391 
    ## Run 211 stress 0.1138396 
    ## Run 212 stress 0.1097987 
    ## Run 213 stress 0.1063129 
    ## Run 214 stress 0.1051923 
    ## ... Procrustes: rmse 0.02215474  max resid 0.07904277 
    ## Run 215 stress 0.1135597 
    ## Run 216 stress 0.107938 
    ## Run 217 stress 0.1058423 
    ## Run 218 stress 0.1063128 
    ## Run 219 stress 0.1050338 
    ## ... Procrustes: rmse 0.008904911  max resid 0.02501951 
    ## Run 220 stress 0.1135596 
    ## Run 221 stress 0.1050144 
    ## ... Procrustes: rmse 0.006658134  max resid 0.02375605 
    ## Run 222 stress 0.1102722 
    ## Run 223 stress 0.1138394 
    ## Run 224 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307642  max resid 0.07849002 
    ## Run 225 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307685  max resid 0.07847707 
    ## Run 226 stress 0.1160552 
    ## Run 227 stress 0.1097987 
    ## Run 228 stress 0.1082722 
    ## Run 229 stress 0.1084964 
    ## Run 230 stress 0.1063129 
    ## Run 231 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214447  max resid 0.07898564 
    ## Run 232 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307483  max resid 0.07848099 
    ## Run 233 stress 0.1136647 
    ## Run 234 stress 0.1136646 
    ## Run 235 stress 0.1079383 
    ## Run 236 stress 0.1174538 
    ## Run 237 stress 0.113594 
    ## Run 238 stress 0.1138393 
    ## Run 239 stress 0.1135599 
    ## Run 240 stress 0.1051402 
    ## ... Procrustes: rmse 0.02307958  max resid 0.07849778 
    ## Run 241 stress 0.1079382 
    ## Run 242 stress 0.1135599 
    ## Run 243 stress 0.1083501 
    ## Run 244 stress 0.1049903 
    ## ... Procrustes: rmse 1.019846e-05  max resid 2.366678e-05 
    ## ... Similar to previous best
    ## Run 245 stress 0.1136646 
    ## Run 246 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 0.006714249  max resid 0.02571064 
    ## Run 247 stress 0.1160552 
    ## Run 248 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181858  max resid 0.07691689 
    ## Run 249 stress 0.1050338 
    ## ... Procrustes: rmse 0.006252536  max resid 0.02184587 
    ## Run 250 stress 0.1049792 
    ## ... Procrustes: rmse 3.280448e-05  max resid 9.526798e-05 
    ## ... Similar to previous best
    ## Run 251 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181147  max resid 0.07689661 
    ## Run 252 stress 0.1063128 
    ## Run 253 stress 0.1097987 
    ## Run 254 stress 0.1135941 
    ## Run 255 stress 0.1079382 
    ## Run 256 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284225  max resid 0.077325 
    ## Run 257 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284404  max resid 0.07733298 
    ## Run 258 stress 0.1165301 
    ## Run 259 stress 0.1084964 
    ## Run 260 stress 0.1088547 
    ## Run 261 stress 0.1097987 
    ## Run 262 stress 0.1135597 
    ## Run 263 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284428  max resid 0.07732813 
    ## Run 264 stress 0.1049903 
    ## ... Procrustes: rmse 0.006705405  max resid 0.02572104 
    ## Run 265 stress 0.1063128 
    ## Run 266 stress 0.1135598 
    ## Run 267 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179984  max resid 0.07685496 
    ## Run 268 stress 0.1050338 
    ## ... Procrustes: rmse 0.006245012  max resid 0.02181713 
    ## Run 269 stress 0.1102722 
    ## Run 270 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180988  max resid 0.07688924 
    ## Run 271 stress 0.1063128 
    ## Run 272 stress 0.1135597 
    ## Run 273 stress 0.1049903 
    ## ... Procrustes: rmse 0.006717111  max resid 0.02578397 
    ## Run 274 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284394  max resid 0.07732314 
    ## Run 275 stress 0.1135598 
    ## Run 276 stress 0.1135598 
    ## Run 277 stress 0.1084964 
    ## Run 278 stress 0.10835 
    ## Run 279 stress 0.113594 
    ## Run 280 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181121  max resid 0.07690159 
    ## Run 281 stress 0.1079381 
    ## Run 282 stress 0.1063128 
    ## Run 283 stress 0.1079383 
    ## Run 284 stress 0.1309391 
    ## Run 285 stress 0.1063128 
    ## Run 286 stress 0.1088547 
    ## Run 287 stress 0.1138394 
    ## Run 288 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284471  max resid 0.07731921 
    ## Run 289 stress 0.1079382 
    ## Run 290 stress 0.1063128 
    ## Run 291 stress 0.1136646 
    ## Run 292 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284509  max resid 0.0773369 
    ## Run 293 stress 0.113594 
    ## Run 294 stress 0.1084965 
    ## Run 295 stress 0.1049903 
    ## ... Procrustes: rmse 0.006710046  max resid 0.02574653 
    ## Run 296 stress 0.1136646 
    ## Run 297 stress 0.1082722 
    ## Run 298 stress 0.1135597 
    ## Run 299 stress 0.1097987 
    ## Run 300 stress 0.1050338 
    ## ... Procrustes: rmse 0.006261064  max resid 0.02189049 
    ## Run 301 stress 0.1135597 
    ## Run 302 stress 0.1160552 
    ## Run 303 stress 0.1050338 
    ## ... Procrustes: rmse 0.006268194  max resid 0.02194417 
    ## Run 304 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181168  max resid 0.07690323 
    ## Run 305 stress 0.1160544 
    ## Run 306 stress 0.1082722 
    ## Run 307 stress 0.1088547 
    ## Run 308 stress 0.113594 
    ## Run 309 stress 0.1135597 
    ## Run 310 stress 0.1063128 
    ## Run 311 stress 0.1063128 
    ## Run 312 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284445  max resid 0.07733361 
    ## Run 313 stress 0.1084965 
    ## Run 314 stress 0.1063128 
    ## Run 315 stress 0.1050144 
    ## ... Procrustes: rmse 0.009655743  max resid 0.0260778 
    ## Run 316 stress 0.1088548 
    ## Run 317 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284701  max resid 0.07735316 
    ## Run 318 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181785  max resid 0.07691505 
    ## Run 319 stress 0.1084965 
    ## Run 320 stress 0.1102722 
    ## Run 321 stress 0.1050144 
    ## ... Procrustes: rmse 0.009656751  max resid 0.02608393 
    ## Run 322 stress 0.1084966 
    ## Run 323 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284415  max resid 0.07733212 
    ## Run 324 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181423  max resid 0.0769022 
    ## Run 325 stress 0.1097988 
    ## Run 326 stress 0.1049903 
    ## ... Procrustes: rmse 0.006723395  max resid 0.02580554 
    ## Run 327 stress 0.1050144 
    ## ... Procrustes: rmse 0.009680435  max resid 0.02621906 
    ## Run 328 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284445  max resid 0.07728529 
    ## Run 329 stress 0.1050338 
    ## ... Procrustes: rmse 0.006257418  max resid 0.021871 
    ## Run 330 stress 0.10835 
    ## Run 331 stress 0.1058423 
    ## Run 332 stress 0.1084964 
    ## Run 333 stress 0.1051924 
    ## ... Procrustes: rmse 0.02283525  max resid 0.07736303 
    ## Run 334 stress 0.113594 
    ## Run 335 stress 0.107938 
    ## Run 336 stress 0.1136646 
    ## Run 337 stress 0.1082722 
    ## Run 338 stress 0.1135598 
    ## Run 339 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284763  max resid 0.07731068 
    ## Run 340 stress 0.1102722 
    ## Run 341 stress 0.1063128 
    ## Run 342 stress 0.1135596 
    ## Run 343 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284291  max resid 0.07734942 
    ## Run 344 stress 0.1063128 
    ## Run 345 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180894  max resid 0.0768719 
    ## Run 346 stress 0.1097987 
    ## Run 347 stress 0.1097989 
    ## Run 348 stress 0.1050144 
    ## ... Procrustes: rmse 0.009648992  max resid 0.02602351 
    ## Run 349 stress 0.1063128 
    ## Run 350 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181577  max resid 0.07690896 
    ## Run 351 stress 0.1079381 
    ## Run 352 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284429  max resid 0.07734295 
    ## Run 353 stress 0.1082722 
    ## Run 354 stress 0.1083502 
    ## Run 355 stress 0.1063128 
    ## Run 356 stress 0.1135597 
    ## Run 357 stress 0.1083996 
    ## Run 358 stress 0.1084965 
    ## Run 359 stress 0.1136646 
    ## Run 360 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181243  max resid 0.0768973 
    ## Run 361 stress 0.1136646 
    ## Run 362 stress 0.1135941 
    ## Run 363 stress 0.1063128 
    ## Run 364 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228411  max resid 0.07729381 
    ## Run 365 stress 0.10835 
    ## Run 366 stress 0.1135598 
    ## Run 367 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218103  max resid 0.07690007 
    ## Run 368 stress 0.1083999 
    ## Run 369 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180849  max resid 0.076866 
    ## Run 370 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284334  max resid 0.0773181 
    ## Run 371 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284376  max resid 0.07734795 
    ## Run 372 stress 0.1082724 
    ## Run 373 stress 0.1082722 
    ## Run 374 stress 0.1138394 
    ## Run 375 stress 0.1049903 
    ## ... Procrustes: rmse 0.006714854  max resid 0.02576861 
    ## Run 376 stress 0.1082722 
    ## Run 377 stress 0.1136646 
    ## Run 378 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180848  max resid 0.07686118 
    ## Run 379 stress 0.1084965 
    ## Run 380 stress 0.1079381 
    ## Run 381 stress 0.1135598 
    ## Run 382 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181057  max resid 0.07689575 
    ## Run 383 stress 0.1135598 
    ## Run 384 stress 0.1063128 
    ## Run 385 stress 0.1160563 
    ## Run 386 stress 0.1138985 
    ## Run 387 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284488  max resid 0.07733583 
    ## Run 388 stress 0.1309397 
    ## Run 389 stress 0.1135597 
    ## Run 390 stress 0.1050338 
    ## ... Procrustes: rmse 0.006273162  max resid 0.0219472 
    ## Run 391 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 3.818753e-05  max resid 0.0001088474 
    ## ... Similar to previous best
    ## Run 392 stress 0.1050144 
    ## ... Procrustes: rmse 0.009669176  max resid 0.0261008 
    ## Run 393 stress 0.1138396 
    ## Run 394 stress 0.1063128 
    ## Run 395 stress 0.1051924 
    ## ... Procrustes: rmse 0.02283861  max resid 0.07725556 
    ## Run 396 stress 0.1063128 
    ## Run 397 stress 0.1136646 
    ## Run 398 stress 0.130938 
    ## Run 399 stress 0.1050338 
    ## ... Procrustes: rmse 0.006254733  max resid 0.02186412 
    ## Run 400 stress 0.1102723 
    ## Run 401 stress 0.1051402 
    ## ... Procrustes: rmse 0.0217964  max resid 0.07683464 
    ## Run 402 stress 0.1050144 
    ## ... Procrustes: rmse 0.009667281  max resid 0.02609845 
    ## Run 403 stress 0.1082722 
    ## Run 404 stress 0.1102721 
    ## Run 405 stress 0.1097987 
    ## Run 406 stress 0.1084002 
    ## Run 407 stress 0.1050338 
    ## ... Procrustes: rmse 0.00625485  max resid 0.02185505 
    ## Run 408 stress 0.1051402 
    ## ... Procrustes: rmse 0.0217955  max resid 0.07682591 
    ## Run 409 stress 0.1084965 
    ## Run 410 stress 0.1135598 
    ## Run 411 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283458  max resid 0.07727911 
    ## Run 412 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283294  max resid 0.07726728 
    ## Run 413 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179486  max resid 0.0768302 
    ## Run 414 stress 0.1083501 
    ## Run 415 stress 0.1088547 
    ## Run 416 stress 0.1082722 
    ## Run 417 stress 0.1135597 
    ## Run 418 stress 0.107938 
    ## Run 419 stress 0.1135597 
    ## Run 420 stress 0.1102724 
    ## Run 421 stress 0.1102724 
    ## Run 422 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179479  max resid 0.07682935 
    ## Run 423 stress 0.1050338 
    ## ... Procrustes: rmse 0.006255184  max resid 0.02185324 
    ## Run 424 stress 0.1136646 
    ## Run 425 stress 0.1136646 
    ## Run 426 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179506  max resid 0.07681541 
    ## Run 427 stress 0.1082722 
    ## Run 428 stress 0.1135597 
    ## Run 429 stress 0.1082722 
    ## Run 430 stress 0.1079381 
    ## Run 431 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 2.198787e-05  max resid 6.015956e-05 
    ## ... Similar to previous best
    ## Run 432 stress 0.1049903 
    ## ... Procrustes: rmse 0.006778088  max resid 0.02603412 
    ## Run 433 stress 0.1050144 
    ## ... Procrustes: rmse 0.009681854  max resid 0.02620916 
    ## Run 434 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180358  max resid 0.07686101 
    ## Run 435 stress 0.11356 
    ## Run 436 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284595  max resid 0.07728933 
    ## Run 437 stress 0.1063129 
    ## Run 438 stress 0.1050144 
    ## ... Procrustes: rmse 0.009661601  max resid 0.02608445 
    ## Run 439 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284201  max resid 0.07732055 
    ## Run 440 stress 0.1084964 
    ## Run 441 stress 0.1138984 
    ## Run 442 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180349  max resid 0.07684719 
    ## Run 443 stress 0.1050338 
    ## ... Procrustes: rmse 0.006261448  max resid 0.02187662 
    ## Run 444 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180367  max resid 0.07686125 
    ## Run 445 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284004  max resid 0.07730368 
    ## Run 446 stress 0.1088547 
    ## Run 447 stress 0.1049903 
    ## ... Procrustes: rmse 0.006730427  max resid 0.02581077 
    ## Run 448 stress 0.1174542 
    ## Run 449 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283968  max resid 0.07729432 
    ## Run 450 stress 0.1082722 
    ## Run 451 stress 0.1135596 
    ## Run 452 stress 0.1063128 
    ## Run 453 stress 0.10835 
    ## Run 454 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284086  max resid 0.07730802 
    ## Run 455 stress 0.1063128 
    ## Run 456 stress 0.1097987 
    ## Run 457 stress 0.1084966 
    ## Run 458 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181326  max resid 0.07689284 
    ## Run 459 stress 0.1135598 
    ## Run 460 stress 0.1050144 
    ## ... Procrustes: rmse 0.009670693  max resid 0.02613625 
    ## Run 461 stress 0.1084965 
    ## Run 462 stress 0.1083499 
    ## Run 463 stress 0.1135598 
    ## Run 464 stress 0.1135941 
    ## Run 465 stress 0.1063128 
    ## Run 466 stress 0.1084966 
    ## Run 467 stress 0.1083499 
    ## Run 468 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180517  max resid 0.0768715 
    ## Run 469 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180344  max resid 0.07685972 
    ## Run 470 stress 0.1083501 
    ## Run 471 stress 0.1084964 
    ## Run 472 stress 0.1135598 
    ## Run 473 stress 0.1084965 
    ## Run 474 stress 0.1050144 
    ## ... Procrustes: rmse 0.00965654  max resid 0.0260424 
    ## Run 475 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284105  max resid 0.0772956 
    ## Run 476 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181243  max resid 0.0768902 
    ## Run 477 stress 0.1058423 
    ## Run 478 stress 0.1097988 
    ## Run 479 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180374  max resid 0.07686219 
    ## Run 480 stress 0.1050144 
    ## ... Procrustes: rmse 0.00967671  max resid 0.02616432 
    ## Run 481 stress 0.1160541 
    ## Run 482 stress 0.1084965 
    ## Run 483 stress 0.1088547 
    ## Run 484 stress 0.1049903 
    ## ... Procrustes: rmse 0.006731412  max resid 0.0258172 
    ## Run 485 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284327  max resid 0.07731401 
    ## Run 486 stress 0.1083501 
    ## Run 487 stress 0.1160544 
    ## Run 488 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181009  max resid 0.07688127 
    ## Run 489 stress 0.1063128 
    ## Run 490 stress 0.1097988 
    ## Run 491 stress 0.1079382 
    ## Run 492 stress 0.1063128 
    ## Run 493 stress 0.1088547 
    ## Run 494 stress 0.113594 
    ## Run 495 stress 0.1050144 
    ## ... Procrustes: rmse 0.009662716  max resid 0.02609013 
    ## Run 496 stress 0.1083501 
    ## Run 497 stress 0.1063128 
    ## Run 498 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180377  max resid 0.07685487 
    ## Run 499 stress 0.1083996 
    ## Run 500 stress 0.1135597 
    ## *** Best solution repeated 1 times

``` r
# Surveyed sites 
SD_beta_NMDS <- metaMDS(SD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09021173 
    ## Run 1 stress 0.09044623 
    ## ... Procrustes: rmse 0.01003653  max resid 0.03411346 
    ## Run 2 stress 0.08938548 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04159085  max resid 0.1192769 
    ## Run 3 stress 0.08938551 
    ## ... Procrustes: rmse 0.0003786196  max resid 0.001119865 
    ## ... Similar to previous best
    ## Run 4 stress 0.08926078 
    ## ... New best solution
    ## ... Procrustes: rmse 0.011765  max resid 0.04006784 
    ## Run 5 stress 0.1092126 
    ## Run 6 stress 0.09088606 
    ## Run 7 stress 0.09044611 
    ## Run 8 stress 0.1076302 
    ## Run 9 stress 0.0901871 
    ## Run 10 stress 0.1056895 
    ## Run 11 stress 0.1092131 
    ## Run 12 stress 0.08946331 
    ## ... Procrustes: rmse 0.03754047  max resid 0.1178074 
    ## Run 13 stress 0.1079015 
    ## Run 14 stress 0.08938961 
    ## ... Procrustes: rmse 0.03598851  max resid 0.1183379 
    ## Run 15 stress 0.08946331 
    ## ... Procrustes: rmse 0.0375356  max resid 0.1177983 
    ## Run 16 stress 0.1071322 
    ## Run 17 stress 0.09130101 
    ## Run 18 stress 0.1087862 
    ## Run 19 stress 0.1118951 
    ## Run 20 stress 0.09130094 
    ## Run 21 stress 0.08938542 
    ## ... Procrustes: rmse 0.01179564  max resid 0.03976933 
    ## Run 22 stress 0.1092131 
    ## Run 23 stress 0.09109108 
    ## Run 24 stress 0.08938538 
    ## ... Procrustes: rmse 0.01181866  max resid 0.03976598 
    ## Run 25 stress 0.08938969 
    ## ... Procrustes: rmse 0.03595905  max resid 0.118291 
    ## Run 26 stress 0.08926106 
    ## ... Procrustes: rmse 0.0001891038  max resid 0.0006468045 
    ## ... Similar to previous best
    ## Run 27 stress 0.1056902 
    ## Run 28 stress 0.09503445 
    ## Run 29 stress 0.09594335 
    ## Run 30 stress 0.1071321 
    ## Run 31 stress 0.08938962 
    ## ... Procrustes: rmse 0.0359791  max resid 0.118324 
    ## Run 32 stress 0.1074434 
    ## Run 33 stress 0.109145 
    ## Run 34 stress 0.1074434 
    ## Run 35 stress 0.08926074 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000358288  max resid 0.0008659373 
    ## ... Similar to previous best
    ## Run 36 stress 0.08938545 
    ## ... Procrustes: rmse 0.01177668  max resid 0.03947498 
    ## Run 37 stress 0.09592144 
    ## Run 38 stress 0.0894633 
    ## ... Procrustes: rmse 0.03744767  max resid 0.1177244 
    ## Run 39 stress 0.08938538 
    ## ... Procrustes: rmse 0.01175254  max resid 0.03942597 
    ## Run 40 stress 0.1061297 
    ## Run 41 stress 0.1060433 
    ## Run 42 stress 0.09039131 
    ## Run 43 stress 0.1060428 
    ## Run 44 stress 0.1091878 
    ## Run 45 stress 0.09130102 
    ## Run 46 stress 0.08926099 
    ## ... Procrustes: rmse 0.000484042  max resid 0.001345761 
    ## ... Similar to previous best
    ## Run 47 stress 0.08946652 
    ## ... Procrustes: rmse 0.03315273  max resid 0.1162814 
    ## Run 48 stress 0.08946659 
    ## ... Procrustes: rmse 0.03313751  max resid 0.1162556 
    ## Run 49 stress 0.09039129 
    ## Run 50 stress 0.1104184 
    ## Run 51 stress 0.08938975 
    ## ... Procrustes: rmse 0.03586724  max resid 0.1181892 
    ## Run 52 stress 0.0893854 
    ## ... Procrustes: rmse 0.01173269  max resid 0.03943687 
    ## Run 53 stress 0.1080637 
    ## Run 54 stress 0.1071322 
    ## Run 55 stress 0.08938963 
    ## ... Procrustes: rmse 0.03588667  max resid 0.1182185 
    ## Run 56 stress 0.1075422 
    ## Run 57 stress 0.1102603 
    ## Run 58 stress 0.1074431 
    ## Run 59 stress 0.1060088 
    ## Run 60 stress 0.09021171 
    ## Run 61 stress 0.09503422 
    ## Run 62 stress 0.09087339 
    ## Run 63 stress 0.09087342 
    ## Run 64 stress 0.1052648 
    ## Run 65 stress 0.08938961 
    ## ... Procrustes: rmse 0.0359029  max resid 0.118252 
    ## Run 66 stress 0.1052651 
    ## Run 67 stress 0.1108316 
    ## Run 68 stress 0.08926078 
    ## ... Procrustes: rmse 0.0003426732  max resid 0.0008285216 
    ## ... Similar to previous best
    ## Run 69 stress 0.1060093 
    ## Run 70 stress 0.09503413 
    ## Run 71 stress 0.105265 
    ## Run 72 stress 0.09503413 
    ## Run 73 stress 0.1061306 
    ## Run 74 stress 0.08938964 
    ## ... Procrustes: rmse 0.03588411  max resid 0.1182153 
    ## Run 75 stress 0.08938539 
    ## ... Procrustes: rmse 0.01177455  max resid 0.03942181 
    ## Run 76 stress 0.1071322 
    ## Run 77 stress 0.08926093 
    ## ... Procrustes: rmse 0.0002951199  max resid 0.0009507564 
    ## ... Similar to previous best
    ## Run 78 stress 0.1103693 
    ## Run 79 stress 0.1088392 
    ## Run 80 stress 0.09018712 
    ## Run 81 stress 0.1068605 
    ## Run 82 stress 0.109256 
    ## Run 83 stress 0.08938546 
    ## ... Procrustes: rmse 0.01182211  max resid 0.03945986 
    ## Run 84 stress 0.09018708 
    ## Run 85 stress 0.08938963 
    ## ... Procrustes: rmse 0.03588891  max resid 0.1182226 
    ## Run 86 stress 0.08938543 
    ## ... Procrustes: rmse 0.01178563  max resid 0.03938409 
    ## Run 87 stress 0.1092094 
    ## Run 88 stress 0.08951717 
    ## ... Procrustes: rmse 0.03500033  max resid 0.1156255 
    ## Run 89 stress 0.09099565 
    ## Run 90 stress 0.08946653 
    ## ... Procrustes: rmse 0.03315188  max resid 0.1162795 
    ## Run 91 stress 0.1061301 
    ## Run 92 stress 0.08938962 
    ## ... Procrustes: rmse 0.0359138  max resid 0.1182629 
    ## Run 93 stress 0.1075421 
    ## Run 94 stress 0.0950344 
    ## Run 95 stress 0.08926066 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001156859  max resid 0.0003375256 
    ## ... Similar to previous best
    ## Run 96 stress 0.08946335 
    ## ... Procrustes: rmse 0.03746335  max resid 0.1176998 
    ## Run 97 stress 0.08946657 
    ## ... Procrustes: rmse 0.03317005  max resid 0.1162756 
    ## Run 98 stress 0.08946652 
    ## ... Procrustes: rmse 0.03318127  max resid 0.1162953 
    ## Run 99 stress 0.08926081 
    ## ... Procrustes: rmse 0.000342796  max resid 0.001018701 
    ## ... Similar to previous best
    ## Run 100 stress 0.1056912 
    ## Run 101 stress 0.08938553 
    ## ... Procrustes: rmse 0.01173378  max resid 0.03953068 
    ## Run 102 stress 0.0892607 
    ## ... Procrustes: rmse 0.000262098  max resid 0.0007621065 
    ## ... Similar to previous best
    ## Run 103 stress 0.09039129 
    ## Run 104 stress 0.09775549 
    ## Run 105 stress 0.0892609 
    ## ... Procrustes: rmse 0.0004258884  max resid 0.001316865 
    ## ... Similar to previous best
    ## Run 106 stress 0.08926065 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001967393  max resid 0.0005490675 
    ## ... Similar to previous best
    ## Run 107 stress 0.08938553 
    ## ... Procrustes: rmse 0.01182858  max resid 0.03961228 
    ## Run 108 stress 0.09592126 
    ## Run 109 stress 0.08938971 
    ## ... Procrustes: rmse 0.03593571  max resid 0.1182618 
    ## Run 110 stress 0.1061302 
    ## Run 111 stress 0.08938539 
    ## ... Procrustes: rmse 0.01182496  max resid 0.03973866 
    ## Run 112 stress 0.09180897 
    ## Run 113 stress 0.10569 
    ## Run 114 stress 0.1065376 
    ## Run 115 stress 0.09087339 
    ## Run 116 stress 0.09503414 
    ## Run 117 stress 0.08938541 
    ## ... Procrustes: rmse 0.01183564  max resid 0.0397355 
    ## Run 118 stress 0.08938963 
    ## ... Procrustes: rmse 0.03595478  max resid 0.1182913 
    ## Run 119 stress 0.1075421 
    ## Run 120 stress 0.1118783 
    ## Run 121 stress 0.1074433 
    ## Run 122 stress 0.08938972 
    ## ... Procrustes: rmse 0.03593409  max resid 0.1182568 
    ## Run 123 stress 0.09018708 
    ## Run 124 stress 0.1079015 
    ## Run 125 stress 0.1067487 
    ## Run 126 stress 0.08926087 
    ## ... Procrustes: rmse 0.000200195  max resid 0.000647489 
    ## ... Similar to previous best
    ## Run 127 stress 0.08926067 
    ## ... Procrustes: rmse 0.0002452171  max resid 0.0008112988 
    ## ... Similar to previous best
    ## Run 128 stress 0.08938544 
    ## ... Procrustes: rmse 0.0118567  max resid 0.03966933 
    ## Run 129 stress 0.1075421 
    ## Run 130 stress 0.08938538 
    ## ... Procrustes: rmse 0.01182952  max resid 0.03972154 
    ## Run 131 stress 0.1075423 
    ## Run 132 stress 0.1075422 
    ## Run 133 stress 0.1085128 
    ## Run 134 stress 0.08926089 
    ## ... Procrustes: rmse 0.0002211612  max resid 0.0006935838 
    ## ... Similar to previous best
    ## Run 135 stress 0.08938539 
    ## ... Procrustes: rmse 0.01181825  max resid 0.03973644 
    ## Run 136 stress 0.1074434 
    ## Run 137 stress 0.1091448 
    ## Run 138 stress 0.1071321 
    ## Run 139 stress 0.08946329 
    ## ... Procrustes: rmse 0.03752249  max resid 0.1177794 
    ## Run 140 stress 0.0893854 
    ## ... Procrustes: rmse 0.01180089  max resid 0.03966226 
    ## Run 141 stress 0.08938967 
    ## ... Procrustes: rmse 0.03594332  max resid 0.1182725 
    ## Run 142 stress 0.1075421 
    ## Run 143 stress 0.1074434 
    ## Run 144 stress 0.1076307 
    ## Run 145 stress 0.1065368 
    ## Run 146 stress 0.09039136 
    ## Run 147 stress 0.1085124 
    ## Run 148 stress 0.09021169 
    ## Run 149 stress 0.1056895 
    ## Run 150 stress 0.1056906 
    ## Run 151 stress 0.1052652 
    ## Run 152 stress 0.09087343 
    ## Run 153 stress 0.08938548 
    ## ... Procrustes: rmse 0.01187609  max resid 0.03967181 
    ## Run 154 stress 0.09021165 
    ## Run 155 stress 0.09503415 
    ## Run 156 stress 0.0894634 
    ## ... Procrustes: rmse 0.03749216  max resid 0.1177389 
    ## Run 157 stress 0.1075423 
    ## Run 158 stress 0.09178286 
    ## Run 159 stress 0.1080637 
    ## Run 160 stress 0.09592134 
    ## Run 161 stress 0.0950342 
    ## Run 162 stress 0.08938962 
    ## ... Procrustes: rmse 0.03598217  max resid 0.1183364 
    ## Run 163 stress 0.1067896 
    ## Run 164 stress 0.1056904 
    ## Run 165 stress 0.08938545 
    ## ... Procrustes: rmse 0.01178263  max resid 0.03974044 
    ## Run 166 stress 0.09018714 
    ## Run 167 stress 0.09503425 
    ## Run 168 stress 0.1067853 
    ## Run 169 stress 0.1052647 
    ## Run 170 stress 0.0910911 
    ## Run 171 stress 0.08938539 
    ## ... Procrustes: rmse 0.01182604  max resid 0.03973654 
    ## Run 172 stress 0.1071323 
    ## Run 173 stress 0.1071321 
    ## Run 174 stress 0.1065358 
    ## Run 175 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001577155  max resid 0.0005083975 
    ## ... Similar to previous best
    ## Run 176 stress 0.08926081 
    ## ... Procrustes: rmse 0.0003277879  max resid 0.001041936 
    ## ... Similar to previous best
    ## Run 177 stress 0.09099571 
    ## Run 178 stress 0.09503427 
    ## Run 179 stress 0.106043 
    ## Run 180 stress 0.0950343 
    ## Run 181 stress 0.09018708 
    ## Run 182 stress 0.08938971 
    ## ... Procrustes: rmse 0.03591355  max resid 0.1182312 
    ## Run 183 stress 0.1092094 
    ## Run 184 stress 0.1091878 
    ## Run 185 stress 0.08946339 
    ## ... Procrustes: rmse 0.03747195  max resid 0.1177012 
    ## Run 186 stress 0.08926081 
    ## ... Procrustes: rmse 0.0003215041  max resid 0.00103542 
    ## ... Similar to previous best
    ## Run 187 stress 0.08938962 
    ## ... Procrustes: rmse 0.03595253  max resid 0.1182911 
    ## Run 188 stress 0.08926101 
    ## ... Procrustes: rmse 0.0003826453  max resid 0.001213238 
    ## ... Similar to previous best
    ## Run 189 stress 0.09612772 
    ## Run 190 stress 0.09109127 
    ## Run 191 stress 0.09130095 
    ## Run 192 stress 0.08938537 
    ## ... Procrustes: rmse 0.01184778  max resid 0.03967609 
    ## Run 193 stress 0.08926111 
    ## ... Procrustes: rmse 0.0004840047  max resid 0.001518853 
    ## ... Similar to previous best
    ## Run 194 stress 0.08926072 
    ## ... Procrustes: rmse 0.0002408886  max resid 0.0007733648 
    ## ... Similar to previous best
    ## Run 195 stress 0.1076304 
    ## Run 196 stress 0.0950345 
    ## Run 197 stress 0.09021166 
    ## Run 198 stress 0.08938544 
    ## ... Procrustes: rmse 0.01189596  max resid 0.03966003 
    ## Run 199 stress 0.1074433 
    ## Run 200 stress 0.09087339 
    ## Run 201 stress 0.1075421 
    ## Run 202 stress 0.08938979 
    ## ... Procrustes: rmse 0.03589741  max resid 0.1182057 
    ## Run 203 stress 0.08926106 
    ## ... Procrustes: rmse 0.0004760266  max resid 0.001552049 
    ## ... Similar to previous best
    ## Run 204 stress 0.09021165 
    ## Run 205 stress 0.1116432 
    ## Run 206 stress 0.0895172 
    ## ... Procrustes: rmse 0.03504178  max resid 0.1156436 
    ## Run 207 stress 0.1052647 
    ## Run 208 stress 0.1076306 
    ## Run 209 stress 0.09130086 
    ## Run 210 stress 0.08926065 
    ## ... Procrustes: rmse 6.951925e-05  max resid 0.0001717767 
    ## ... Similar to previous best
    ## Run 211 stress 0.1091825 
    ## Run 212 stress 0.08938967 
    ## ... Procrustes: rmse 0.03591686  max resid 0.1182326 
    ## Run 213 stress 0.1065327 
    ## Run 214 stress 0.1092127 
    ## Run 215 stress 0.08926074 
    ## ... Procrustes: rmse 0.0001237648  max resid 0.0003503 
    ## ... Similar to previous best
    ## Run 216 stress 0.08926063 
    ## ... Procrustes: rmse 2.777749e-05  max resid 7.357442e-05 
    ## ... Similar to previous best
    ## Run 217 stress 0.1071322 
    ## Run 218 stress 0.1071321 
    ## Run 219 stress 0.1101445 
    ## Run 220 stress 0.08938542 
    ## ... Procrustes: rmse 0.0118158  max resid 0.03970368 
    ## Run 221 stress 0.1063192 
    ## Run 222 stress 0.1063192 
    ## Run 223 stress 0.08938966 
    ## ... Procrustes: rmse 0.03591768  max resid 0.1182344 
    ## Run 224 stress 0.1074432 
    ## Run 225 stress 0.0893855 
    ## ... Procrustes: rmse 0.01183492  max resid 0.03976029 
    ## Run 226 stress 0.08946659 
    ## ... Procrustes: rmse 0.03324041  max resid 0.116393 
    ## Run 227 stress 0.09088627 
    ## Run 228 stress 0.08946653 
    ## ... Procrustes: rmse 0.03319023  max resid 0.116297 
    ## Run 229 stress 0.08926064 
    ## ... Procrustes: rmse 0.0001208478  max resid 0.0003432832 
    ## ... Similar to previous best
    ## Run 230 stress 0.09088614 
    ## Run 231 stress 0.1075422 
    ## Run 232 stress 0.09130085 
    ## Run 233 stress 0.1068391 
    ## Run 234 stress 0.1060428 
    ## Run 235 stress 0.08938552 
    ## ... Procrustes: rmse 0.0117684  max resid 0.03958373 
    ## Run 236 stress 0.09018707 
    ## Run 237 stress 0.1120101 
    ## Run 238 stress 0.1071321 
    ## Run 239 stress 0.09503414 
    ## Run 240 stress 0.08938962 
    ## ... Procrustes: rmse 0.03593214  max resid 0.1182625 
    ## Run 241 stress 0.1075422 
    ## Run 242 stress 0.09109114 
    ## Run 243 stress 0.1063186 
    ## Run 244 stress 0.08926115 
    ## ... Procrustes: rmse 0.0005200582  max resid 0.001649595 
    ## ... Similar to previous best
    ## Run 245 stress 0.09503425 
    ## Run 246 stress 0.1076304 
    ## Run 247 stress 0.1076303 
    ## Run 248 stress 0.08946653 
    ## ... Procrustes: rmse 0.03319203  max resid 0.1163069 
    ## Run 249 stress 0.09099519 
    ## Run 250 stress 0.09018711 
    ## Run 251 stress 0.08926079 
    ## ... Procrustes: rmse 0.0003108613  max resid 0.001006169 
    ## ... Similar to previous best
    ## Run 252 stress 0.1075421 
    ## Run 253 stress 0.1080635 
    ## Run 254 stress 0.08946656 
    ## ... Procrustes: rmse 0.03322923  max resid 0.1163832 
    ## Run 255 stress 0.1056902 
    ## Run 256 stress 0.090886 
    ## Run 257 stress 0.1075423 
    ## Run 258 stress 0.1060086 
    ## Run 259 stress 0.1052647 
    ## Run 260 stress 0.1056905 
    ## Run 261 stress 0.1071322 
    ## Run 262 stress 0.09088606 
    ## Run 263 stress 0.09592134 
    ## Run 264 stress 0.09503439 
    ## Run 265 stress 0.1071323 
    ## Run 266 stress 0.09130088 
    ## Run 267 stress 0.08946657 
    ## ... Procrustes: rmse 0.03318342  max resid 0.1162875 
    ## Run 268 stress 0.09775569 
    ## Run 269 stress 0.08926095 
    ## ... Procrustes: rmse 0.0004154112  max resid 0.001375877 
    ## ... Similar to previous best
    ## Run 270 stress 0.09503431 
    ## Run 271 stress 0.09099581 
    ## Run 272 stress 0.09775568 
    ## Run 273 stress 0.1086561 
    ## Run 274 stress 0.1063193 
    ## Run 275 stress 0.08938961 
    ## ... Procrustes: rmse 0.03594211  max resid 0.1182713 
    ## Run 276 stress 0.09039131 
    ## Run 277 stress 0.1065365 
    ## Run 278 stress 0.09109121 
    ## Run 279 stress 0.1056899 
    ## Run 280 stress 0.09044611 
    ## Run 281 stress 0.0908736 
    ## Run 282 stress 0.08938538 
    ## ... Procrustes: rmse 0.01182718  max resid 0.03966086 
    ## Run 283 stress 0.09087341 
    ## Run 284 stress 0.1089907 
    ## Run 285 stress 0.09088603 
    ## Run 286 stress 0.1064195 
    ## Run 287 stress 0.08926064 
    ## ... Procrustes: rmse 4.021033e-05  max resid 0.0001059956 
    ## ... Similar to previous best
    ## Run 288 stress 0.1105342 
    ## Run 289 stress 0.08938964 
    ## ... Procrustes: rmse 0.03592447  max resid 0.1182458 
    ## Run 290 stress 0.1052651 
    ## Run 291 stress 0.0894665 
    ## ... Procrustes: rmse 0.0332079  max resid 0.116331 
    ## Run 292 stress 0.09178289 
    ## Run 293 stress 0.09109108 
    ## Run 294 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 4.31369e-05  max resid 0.0001405368 
    ## ... Similar to previous best
    ## Run 295 stress 0.08938962 
    ## ... Procrustes: rmse 0.0359658  max resid 0.1183084 
    ## Run 296 stress 0.08926078 
    ## ... Procrustes: rmse 0.0001869521  max resid 0.0004974722 
    ## ... Similar to previous best
    ## Run 297 stress 0.08926082 
    ## ... Procrustes: rmse 0.0002827882  max resid 0.0008375323 
    ## ... Similar to previous best
    ## Run 298 stress 0.08938968 
    ## ... Procrustes: rmse 0.03592406  max resid 0.1182418 
    ## Run 299 stress 0.08938964 
    ## ... Procrustes: rmse 0.03595292  max resid 0.1182882 
    ## Run 300 stress 0.0893854 
    ## ... Procrustes: rmse 0.01180573  max resid 0.03968058 
    ## Run 301 stress 0.1061306 
    ## Run 302 stress 0.08946653 
    ## ... Procrustes: rmse 0.03320315  max resid 0.1163182 
    ## Run 303 stress 0.1074431 
    ## Run 304 stress 0.08938964 
    ## ... Procrustes: rmse 0.0359713  max resid 0.1183236 
    ## Run 305 stress 0.09109112 
    ## Run 306 stress 0.08926064 
    ## ... Procrustes: rmse 6.85493e-05  max resid 0.000190986 
    ## ... Similar to previous best
    ## Run 307 stress 0.1075423 
    ## Run 308 stress 0.1075421 
    ## Run 309 stress 0.1060094 
    ## Run 310 stress 0.09503413 
    ## Run 311 stress 0.1071321 
    ## Run 312 stress 0.1074432 
    ## Run 313 stress 0.08926068 
    ## ... Procrustes: rmse 0.0001171138  max resid 0.0003905042 
    ## ... Similar to previous best
    ## Run 314 stress 0.1075421 
    ## Run 315 stress 0.1060432 
    ## Run 316 stress 0.08938964 
    ## ... Procrustes: rmse 0.03597202  max resid 0.118317 
    ## Run 317 stress 0.1108232 
    ## Run 318 stress 0.09503445 
    ## Run 319 stress 0.112832 
    ## Run 320 stress 0.08926085 
    ## ... Procrustes: rmse 0.0003037443  max resid 0.0009875914 
    ## ... Similar to previous best
    ## Run 321 stress 0.0893897 
    ## ... Procrustes: rmse 0.03592051  max resid 0.1182372 
    ## Run 322 stress 0.08926063 
    ## ... Procrustes: rmse 4.597475e-05  max resid 0.0001047026 
    ## ... Similar to previous best
    ## Run 323 stress 0.08938962 
    ## ... Procrustes: rmse 0.03594918  max resid 0.118284 
    ## Run 324 stress 0.0893855 
    ## ... Procrustes: rmse 0.01176193  max resid 0.03960625 
    ## Run 325 stress 0.1079019 
    ## Run 326 stress 0.08938539 
    ## ... Procrustes: rmse 0.01183898  max resid 0.03964695 
    ## Run 327 stress 0.09503428 
    ## Run 328 stress 0.1066195 
    ## Run 329 stress 0.08926064 
    ## ... Procrustes: rmse 7.497113e-05  max resid 0.0001667893 
    ## ... Similar to previous best
    ## Run 330 stress 0.08926104 
    ## ... Procrustes: rmse 0.000418373  max resid 0.001375624 
    ## ... Similar to previous best
    ## Run 331 stress 0.08926079 
    ## ... Procrustes: rmse 0.0002041054  max resid 0.0005784325 
    ## ... Similar to previous best
    ## Run 332 stress 0.09087353 
    ## Run 333 stress 0.1064188 
    ## Run 334 stress 0.08926098 
    ## ... Procrustes: rmse 0.0003949753  max resid 0.001292184 
    ## ... Similar to previous best
    ## Run 335 stress 0.08946331 
    ## ... Procrustes: rmse 0.03749622  max resid 0.1177407 
    ## Run 336 stress 0.0892609 
    ## ... Procrustes: rmse 0.0003514611  max resid 0.00111958 
    ## ... Similar to previous best
    ## Run 337 stress 0.1071322 
    ## Run 338 stress 0.08938541 
    ## ... Procrustes: rmse 0.01180968  max resid 0.03970441 
    ## Run 339 stress 0.1108634 
    ## Run 340 stress 0.1076304 
    ## Run 341 stress 0.08938545 
    ## ... Procrustes: rmse 0.01180244  max resid 0.03960146 
    ## Run 342 stress 0.08946339 
    ## ... Procrustes: rmse 0.03755648  max resid 0.1178364 
    ## Run 343 stress 0.08926072 
    ## ... Procrustes: rmse 0.0002036003  max resid 0.0006543191 
    ## ... Similar to previous best
    ## Run 344 stress 0.08938539 
    ## ... Procrustes: rmse 0.01185122  max resid 0.03970989 
    ## Run 345 stress 0.08938539 
    ## ... Procrustes: rmse 0.01184718  max resid 0.03971239 
    ## Run 346 stress 0.09039133 
    ## Run 347 stress 0.1074433 
    ## Run 348 stress 0.1101442 
    ## Run 349 stress 0.09262372 
    ## Run 350 stress 0.08938544 
    ## ... Procrustes: rmse 0.01183857  max resid 0.03970527 
    ## Run 351 stress 0.1056897 
    ## Run 352 stress 0.08938978 
    ## ... Procrustes: rmse 0.03590949  max resid 0.1182226 
    ## Run 353 stress 0.1052647 
    ## Run 354 stress 0.08926064 
    ## ... Procrustes: rmse 6.557912e-05  max resid 0.000210851 
    ## ... Similar to previous best
    ## Run 355 stress 0.1108326 
    ## Run 356 stress 0.09592128 
    ## Run 357 stress 0.09088605 
    ## Run 358 stress 0.1075421 
    ## Run 359 stress 0.09087345 
    ## Run 360 stress 0.1052648 
    ## Run 361 stress 0.08946653 
    ## ... Procrustes: rmse 0.03320423  max resid 0.1163206 
    ## Run 362 stress 0.1056899 
    ## Run 363 stress 0.08938545 
    ## ... Procrustes: rmse 0.01188014  max resid 0.03964841 
    ## Run 364 stress 0.09503444 
    ## Run 365 stress 0.09099566 
    ## Run 366 stress 0.1108231 
    ## Run 367 stress 0.08938562 
    ## ... Procrustes: rmse 0.01188012  max resid 0.03957847 
    ## Run 368 stress 0.09775566 
    ## Run 369 stress 0.09592153 
    ## Run 370 stress 0.1087557 
    ## Run 371 stress 0.1062567 
    ## Run 372 stress 0.1056895 
    ## Run 373 stress 0.1074432 
    ## Run 374 stress 0.1093483 
    ## Run 375 stress 0.1056893 
    ## Run 376 stress 0.1075421 
    ## Run 377 stress 0.1091824 
    ## Run 378 stress 0.09592152 
    ## Run 379 stress 0.08926065 
    ## ... Procrustes: rmse 9.673889e-05  max resid 0.0003097775 
    ## ... Similar to previous best
    ## Run 380 stress 0.08938538 
    ## ... Procrustes: rmse 0.0118424  max resid 0.03969395 
    ## Run 381 stress 0.0917829 
    ## Run 382 stress 0.08938963 
    ## ... Procrustes: rmse 0.03596677  max resid 0.1183156 
    ## Run 383 stress 0.09087345 
    ## Run 384 stress 0.09087339 
    ## Run 385 stress 0.10613 
    ## Run 386 stress 0.09039129 
    ## Run 387 stress 0.1062569 
    ## Run 388 stress 0.1074433 
    ## Run 389 stress 0.08926074 
    ## ... Procrustes: rmse 0.000220211  max resid 0.0006984725 
    ## ... Similar to previous best
    ## Run 390 stress 0.08938966 
    ## ... Procrustes: rmse 0.03597366  max resid 0.1183293 
    ## Run 391 stress 0.1092096 
    ## Run 392 stress 0.09088604 
    ## Run 393 stress 0.08926065 
    ## ... Procrustes: rmse 6.999475e-05  max resid 0.0001600155 
    ## ... Similar to previous best
    ## Run 394 stress 0.08926066 
    ## ... Procrustes: rmse 0.0001189136  max resid 0.0003779813 
    ## ... Similar to previous best
    ## Run 395 stress 0.08938962 
    ## ... Procrustes: rmse 0.03594243  max resid 0.1182722 
    ## Run 396 stress 0.1061302 
    ## Run 397 stress 0.08938544 
    ## ... Procrustes: rmse 0.0118813  max resid 0.03973237 
    ## Run 398 stress 0.1101448 
    ## Run 399 stress 0.09044611 
    ## Run 400 stress 0.1061301 
    ## Run 401 stress 0.1076305 
    ## Run 402 stress 0.08938543 
    ## ... Procrustes: rmse 0.01188051  max resid 0.03973104 
    ## Run 403 stress 0.1056896 
    ## Run 404 stress 0.1085127 
    ## Run 405 stress 0.09039132 
    ## Run 406 stress 0.09088611 
    ## Run 407 stress 0.1056893 
    ## Run 408 stress 0.09503441 
    ## Run 409 stress 0.110145 
    ## Run 410 stress 0.1052652 
    ## Run 411 stress 0.08946349 
    ## ... Procrustes: rmse 0.03746463  max resid 0.117678 
    ## Run 412 stress 0.0904462 
    ## Run 413 stress 0.1088372 
    ## Run 414 stress 0.08938964 
    ## ... Procrustes: rmse 0.03593344  max resid 0.1182577 
    ## Run 415 stress 0.1086563 
    ## Run 416 stress 0.1104189 
    ## Run 417 stress 0.1075421 
    ## Run 418 stress 0.09087339 
    ## Run 419 stress 0.09087342 
    ## Run 420 stress 0.09592142 
    ## Run 421 stress 0.09021174 
    ## Run 422 stress 0.09109108 
    ## Run 423 stress 0.1068388 
    ## Run 424 stress 0.09130087 
    ## Run 425 stress 0.08946651 
    ## ... Procrustes: rmse 0.03321356  max resid 0.1163383 
    ## Run 426 stress 0.08926085 
    ## ... Procrustes: rmse 0.0003173455  max resid 0.001035553 
    ## ... Similar to previous best
    ## Run 427 stress 0.1071321 
    ## Run 428 stress 0.09503421 
    ## Run 429 stress 0.08938541 
    ## ... Procrustes: rmse 0.01182928  max resid 0.03971461 
    ## Run 430 stress 0.1071321 
    ## Run 431 stress 0.0950344 
    ## Run 432 stress 0.08946331 
    ## ... Procrustes: rmse 0.03749844  max resid 0.1177435 
    ## Run 433 stress 0.1095629 
    ## Run 434 stress 0.1096137 
    ## Run 435 stress 0.1056906 
    ## Run 436 stress 0.1056897 
    ## Run 437 stress 0.1075422 
    ## Run 438 stress 0.1056897 
    ## Run 439 stress 0.09178294 
    ## Run 440 stress 0.08926063 
    ## ... Procrustes: rmse 1.356311e-05  max resid 3.489998e-05 
    ## ... Similar to previous best
    ## Run 441 stress 0.1061297 
    ## Run 442 stress 0.1056897 
    ## Run 443 stress 0.1065365 
    ## Run 444 stress 0.08946654 
    ## ... Procrustes: rmse 0.03319874  max resid 0.1163115 
    ## Run 445 stress 0.0959213 
    ## Run 446 stress 0.1065373 
    ## Run 447 stress 0.09044611 
    ## Run 448 stress 0.1052651 
    ## Run 449 stress 0.09178295 
    ## Run 450 stress 0.0903914 
    ## Run 451 stress 0.1118949 
    ## Run 452 stress 0.1105341 
    ## Run 453 stress 0.08926078 
    ## ... Procrustes: rmse 0.0002563673  max resid 0.0008143886 
    ## ... Similar to previous best
    ## Run 454 stress 0.08938962 
    ## ... Procrustes: rmse 0.03593965  max resid 0.1182697 
    ## Run 455 stress 0.1071321 
    ## Run 456 stress 0.1067587 
    ## Run 457 stress 0.08938544 
    ## ... Procrustes: rmse 0.01180308  max resid 0.03971733 
    ## Run 458 stress 0.09044616 
    ## Run 459 stress 0.08938968 
    ## ... Procrustes: rmse 0.0359762  max resid 0.1183367 
    ## Run 460 stress 0.08938963 
    ## ... Procrustes: rmse 0.03596294  max resid 0.1182982 
    ## Run 461 stress 0.1074433 
    ## Run 462 stress 0.08926082 
    ## ... Procrustes: rmse 0.0002830702  max resid 0.0008842011 
    ## ... Similar to previous best
    ## Run 463 stress 0.1074432 
    ## Run 464 stress 0.09044612 
    ## Run 465 stress 0.1052647 
    ## Run 466 stress 0.08946651 
    ## ... Procrustes: rmse 0.03321455  max resid 0.1163394 
    ## Run 467 stress 0.1056893 
    ## Run 468 stress 0.0908861 
    ## Run 469 stress 0.0908734 
    ## Run 470 stress 0.08951716 
    ## ... Procrustes: rmse 0.03505966  max resid 0.115641 
    ## Run 471 stress 0.0894665 
    ## ... Procrustes: rmse 0.03321737  max resid 0.1163423 
    ## Run 472 stress 0.08938964 
    ## ... Procrustes: rmse 0.03597114  max resid 0.1183224 
    ## Run 473 stress 0.1067876 
    ## Run 474 stress 0.09228506 
    ## Run 475 stress 0.1066205 
    ## Run 476 stress 0.1056904 
    ## Run 477 stress 0.110372 
    ## Run 478 stress 0.1075422 
    ## Run 479 stress 0.09039129 
    ## Run 480 stress 0.09099521 
    ## Run 481 stress 0.09503422 
    ## Run 482 stress 0.08938961 
    ## ... Procrustes: rmse 0.03596249  max resid 0.118305 
    ## Run 483 stress 0.09503423 
    ## Run 484 stress 0.1056903 
    ## Run 485 stress 0.09130088 
    ## Run 486 stress 0.1065371 
    ## Run 487 stress 0.08946653 
    ## ... Procrustes: rmse 0.03320827  max resid 0.1163297 
    ## Run 488 stress 0.0922851 
    ## Run 489 stress 0.1092128 
    ## Run 490 stress 0.1063191 
    ## Run 491 stress 0.0910911 
    ## Run 492 stress 0.09088606 
    ## Run 493 stress 0.09503441 
    ## Run 494 stress 0.09087343 
    ## Run 495 stress 0.09109117 
    ## Run 496 stress 0.08946329 
    ## ... Procrustes: rmse 0.03750618  max resid 0.1177578 
    ## Run 497 stress 0.1061303 
    ## Run 498 stress 0.1101447 
    ## Run 499 stress 0.1060429 
    ## Run 500 stress 0.08938962 
    ## ... Procrustes: rmse 0.03594195  max resid 0.1182718 
    ## *** Best solution repeated 22 times

``` r
round(SD_beta_NMDS$stress, digits = 2)
```

    ## [1] 0.09

``` r
SD_beta_rep_NMDS <- metaMDS(SD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.3323292 
    ## Run 1 stress 0.3366013 
    ## Run 2 stress 0.3330036 
    ## Run 3 stress 0.3401924 
    ## Run 4 stress 0.3259695 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1929728  max resid 0.2913449 
    ## Run 5 stress 0.3415134 
    ## Run 6 stress 0.3441287 
    ## Run 7 stress 0.3371296 
    ## Run 8 stress 0.3377831 
    ## Run 9 stress 0.3286362 
    ## Run 10 stress 0.3338046 
    ## Run 11 stress 0.3281392 
    ## Run 12 stress 0.3378244 
    ## Run 13 stress 0.3350498 
    ## Run 14 stress 0.3351606 
    ## Run 15 stress 0.3320623 
    ## Run 16 stress 0.337876 
    ## Run 17 stress 0.3445183 
    ## Run 18 stress 0.3807945 
    ## Run 19 stress 0.3397402 
    ## Run 20 stress 0.3423883 
    ## Run 21 stress 0.3426539 
    ## Run 22 stress 0.338709 
    ## Run 23 stress 0.3489196 
    ## Run 24 stress 0.334207 
    ## Run 25 stress 0.342263 
    ## Run 26 stress 0.3354895 
    ## Run 27 stress 0.3363671 
    ## Run 28 stress 0.3395042 
    ## Run 29 stress 0.3338637 
    ## Run 30 stress 0.3458779 
    ## Run 31 stress 0.3406533 
    ## Run 32 stress 0.3469102 
    ## Run 33 stress 0.3497669 
    ## Run 34 stress 0.3307556 
    ## Run 35 stress 0.3339723 
    ## Run 36 stress 0.3371437 
    ## Run 37 stress 0.3476278 
    ## Run 38 stress 0.3298718 
    ## Run 39 stress 0.3282959 
    ## Run 40 stress 0.3345534 
    ## Run 41 stress 0.3364141 
    ## Run 42 stress 0.3339999 
    ## Run 43 stress 0.3313691 
    ## Run 44 stress 0.3305348 
    ## Run 45 stress 0.3514874 
    ## Run 46 stress 0.3296635 
    ## Run 47 stress 0.3341564 
    ## Run 48 stress 0.3330522 
    ## Run 49 stress 0.3319781 
    ## Run 50 stress 0.3324543 
    ## Run 51 stress 0.339769 
    ## Run 52 stress 0.3428678 
    ## Run 53 stress 0.3446944 
    ## Run 54 stress 0.3302087 
    ## Run 55 stress 0.3331162 
    ## Run 56 stress 0.3306197 
    ## Run 57 stress 0.3401496 
    ## Run 58 stress 0.3317256 
    ## Run 59 stress 0.3478028 
    ## Run 60 stress 0.3530935 
    ## Run 61 stress 0.3426533 
    ## Run 62 stress 0.3446687 
    ## Run 63 stress 0.3341751 
    ## Run 64 stress 0.3331313 
    ## Run 65 stress 0.3300056 
    ## Run 66 stress 0.3335969 
    ## Run 67 stress 0.334226 
    ## Run 68 stress 0.3318269 
    ## Run 69 stress 0.3426377 
    ## Run 70 stress 0.3480072 
    ## Run 71 stress 0.3327994 
    ## Run 72 stress 0.3372658 
    ## Run 73 stress 0.3444469 
    ## Run 74 stress 0.3475163 
    ## Run 75 stress 0.3298565 
    ## Run 76 stress 0.3349956 
    ## Run 77 stress 0.3488153 
    ## Run 78 stress 0.3420633 
    ## Run 79 stress 0.3362836 
    ## Run 80 stress 0.3242725 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1630896  max resid 0.4003051 
    ## Run 81 stress 0.3305084 
    ## Run 82 stress 0.3300636 
    ## Run 83 stress 0.334803 
    ## Run 84 stress 0.3476129 
    ## Run 85 stress 0.3331105 
    ## Run 86 stress 0.340324 
    ## Run 87 stress 0.3378001 
    ## Run 88 stress 0.3302552 
    ## Run 89 stress 0.3292484 
    ## Run 90 stress 0.3457915 
    ## Run 91 stress 0.3397205 
    ## Run 92 stress 0.3308172 
    ## Run 93 stress 0.3293426 
    ## Run 94 stress 0.338028 
    ## Run 95 stress 0.3431053 
    ## Run 96 stress 0.3439072 
    ## Run 97 stress 0.3408141 
    ## Run 98 stress 0.3311996 
    ## Run 99 stress 0.3301505 
    ## Run 100 stress 0.3467236 
    ## Run 101 stress 0.343118 
    ## Run 102 stress 0.3327624 
    ## Run 103 stress 0.3401676 
    ## Run 104 stress 0.3392667 
    ## Run 105 stress 0.3401816 
    ## Run 106 stress 0.3435891 
    ## Run 107 stress 0.3435186 
    ## Run 108 stress 0.3349809 
    ## Run 109 stress 0.3276352 
    ## Run 110 stress 0.3463447 
    ## Run 111 stress 0.3349836 
    ## Run 112 stress 0.3361434 
    ## Run 113 stress 0.3501682 
    ## Run 114 stress 0.3329783 
    ## Run 115 stress 0.3299186 
    ## Run 116 stress 0.3377503 
    ## Run 117 stress 0.3481308 
    ## Run 118 stress 0.3405619 
    ## Run 119 stress 0.3454584 
    ## Run 120 stress 0.346219 
    ## Run 121 stress 0.3408677 
    ## Run 122 stress 0.330966 
    ## Run 123 stress 0.3279737 
    ## Run 124 stress 0.337744 
    ## Run 125 stress 0.3294626 
    ## Run 126 stress 0.3499885 
    ## Run 127 stress 0.3438669 
    ## Run 128 stress 0.3490856 
    ## Run 129 stress 0.3400765 
    ## Run 130 stress 0.34037 
    ## Run 131 stress 0.3295809 
    ## Run 132 stress 0.3413308 
    ## Run 133 stress 0.3368286 
    ## Run 134 stress 0.3354691 
    ## Run 135 stress 0.3330295 
    ## Run 136 stress 0.3350107 
    ## Run 137 stress 0.3511715 
    ## Run 138 stress 0.3418651 
    ## Run 139 stress 0.3286158 
    ## Run 140 stress 0.342672 
    ## Run 141 stress 0.3369088 
    ## Run 142 stress 0.3311596 
    ## Run 143 stress 0.3348051 
    ## Run 144 stress 0.3330868 
    ## Run 145 stress 0.3362356 
    ## Run 146 stress 0.3279387 
    ## Run 147 stress 0.3330004 
    ## Run 148 stress 0.3379817 
    ## Run 149 stress 0.3378388 
    ## Run 150 stress 0.3336651 
    ## Run 151 stress 0.3407871 
    ## Run 152 stress 0.3467941 
    ## Run 153 stress 0.3328672 
    ## Run 154 stress 0.3513575 
    ## Run 155 stress 0.3349228 
    ## Run 156 stress 0.3348674 
    ## Run 157 stress 0.3347449 
    ## Run 158 stress 0.3348201 
    ## Run 159 stress 0.3300074 
    ## Run 160 stress 0.335586 
    ## Run 161 stress 0.3302882 
    ## Run 162 stress 0.33905 
    ## Run 163 stress 0.3421028 
    ## Run 164 stress 0.3256435 
    ## Run 165 stress 0.3377422 
    ## Run 166 stress 0.3436659 
    ## Run 167 stress 0.338451 
    ## Run 168 stress 0.3303848 
    ## Run 169 stress 0.3357844 
    ## Run 170 stress 0.3331018 
    ## Run 171 stress 0.3329004 
    ## Run 172 stress 0.3482075 
    ## Run 173 stress 0.3351272 
    ## Run 174 stress 0.3358819 
    ## Run 175 stress 0.3299073 
    ## Run 176 stress 0.3469689 
    ## Run 177 stress 0.3350024 
    ## Run 178 stress 0.3401302 
    ## Run 179 stress 0.3367642 
    ## Run 180 stress 0.3382324 
    ## Run 181 stress 0.3337049 
    ## Run 182 stress 0.3507843 
    ## Run 183 stress 0.3276266 
    ## Run 184 stress 0.3418333 
    ## Run 185 stress 0.3365034 
    ## Run 186 stress 0.3460759 
    ## Run 187 stress 0.3286024 
    ## Run 188 stress 0.3484406 
    ## Run 189 stress 0.3356566 
    ## Run 190 stress 0.3424944 
    ## Run 191 stress 0.3415377 
    ## Run 192 stress 0.3348487 
    ## Run 193 stress 0.3462856 
    ## Run 194 stress 0.341576 
    ## Run 195 stress 0.3351117 
    ## Run 196 stress 0.3340822 
    ## Run 197 stress 0.346515 
    ## Run 198 stress 0.3282024 
    ## Run 199 stress 0.3443044 
    ## Run 200 stress 0.3377703 
    ## Run 201 stress 0.3382187 
    ## Run 202 stress 0.3350482 
    ## Run 203 stress 0.3440934 
    ## Run 204 stress 0.333167 
    ## Run 205 stress 0.3384951 
    ## Run 206 stress 0.3334074 
    ## Run 207 stress 0.3343282 
    ## Run 208 stress 0.3408813 
    ## Run 209 stress 0.3337501 
    ## Run 210 stress 0.3601396 
    ## Run 211 stress 0.3346532 
    ## Run 212 stress 0.3402336 
    ## Run 213 stress 0.3319643 
    ## Run 214 stress 0.325835 
    ## Run 215 stress 0.3397265 
    ## Run 216 stress 0.3491611 
    ## Run 217 stress 0.3343224 
    ## Run 218 stress 0.3407249 
    ## Run 219 stress 0.3290578 
    ## Run 220 stress 0.3457416 
    ## Run 221 stress 0.3306089 
    ## Run 222 stress 0.3374909 
    ## Run 223 stress 0.3394412 
    ## Run 224 stress 0.3375215 
    ## Run 225 stress 0.3478782 
    ## Run 226 stress 0.3330329 
    ## Run 227 stress 0.3399684 
    ## Run 228 stress 0.3351612 
    ## Run 229 stress 0.3345432 
    ## Run 230 stress 0.3334264 
    ## Run 231 stress 0.3308862 
    ## Run 232 stress 0.3437942 
    ## Run 233 stress 0.3423222 
    ## Run 234 stress 0.3332828 
    ## Run 235 stress 0.3467102 
    ## Run 236 stress 0.3445503 
    ## Run 237 stress 0.3294559 
    ## Run 238 stress 0.3490918 
    ## Run 239 stress 0.3457776 
    ## Run 240 stress 0.3366743 
    ## Run 241 stress 0.3452902 
    ## Run 242 stress 0.330522 
    ## Run 243 stress 0.3366632 
    ## Run 244 stress 0.3354439 
    ## Run 245 stress 0.3381743 
    ## Run 246 stress 0.3336896 
    ## Run 247 stress 0.3450978 
    ## Run 248 stress 0.3257467 
    ## Run 249 stress 0.3367893 
    ## Run 250 stress 0.3475592 
    ## Run 251 stress 0.3316413 
    ## Run 252 stress 0.3475075 
    ## Run 253 stress 0.3338607 
    ## Run 254 stress 0.333751 
    ## Run 255 stress 0.3393451 
    ## Run 256 stress 0.33116 
    ## Run 257 stress 0.3347735 
    ## Run 258 stress 0.3330827 
    ## Run 259 stress 0.3426567 
    ## Run 260 stress 0.344244 
    ## Run 261 stress 0.3338938 
    ## Run 262 stress 0.3395195 
    ## Run 263 stress 0.3273307 
    ## Run 264 stress 0.3384735 
    ## Run 265 stress 0.3318162 
    ## Run 266 stress 0.3339676 
    ## Run 267 stress 0.3369831 
    ## Run 268 stress 0.3444913 
    ## Run 269 stress 0.3437785 
    ## Run 270 stress 0.3376294 
    ## Run 271 stress 0.3349802 
    ## Run 272 stress 0.342933 
    ## Run 273 stress 0.3361324 
    ## Run 274 stress 0.337583 
    ## Run 275 stress 0.3306609 
    ## Run 276 stress 0.34092 
    ## Run 277 stress 0.3307253 
    ## Run 278 stress 0.3415979 
    ## Run 279 stress 0.332396 
    ## Run 280 stress 0.338995 
    ## Run 281 stress 0.3353495 
    ## Run 282 stress 0.3455324 
    ## Run 283 stress 0.3362439 
    ## Run 284 stress 0.3400287 
    ## Run 285 stress 0.3331781 
    ## Run 286 stress 0.3435624 
    ## Run 287 stress 0.328814 
    ## Run 288 stress 0.3349079 
    ## Run 289 stress 0.3313388 
    ## Run 290 stress 0.3331892 
    ## Run 291 stress 0.3458261 
    ## Run 292 stress 0.3486638 
    ## Run 293 stress 0.3328943 
    ## Run 294 stress 0.3423202 
    ## Run 295 stress 0.3337862 
    ## Run 296 stress 0.3288997 
    ## Run 297 stress 0.3401038 
    ## Run 298 stress 0.3445629 
    ## Run 299 stress 0.3383937 
    ## Run 300 stress 0.3346804 
    ## Run 301 stress 0.3363224 
    ## Run 302 stress 0.333736 
    ## Run 303 stress 0.3347723 
    ## Run 304 stress 0.3388679 
    ## Run 305 stress 0.3410684 
    ## Run 306 stress 0.3468819 
    ## Run 307 stress 0.3446423 
    ## Run 308 stress 0.3366023 
    ## Run 309 stress 0.3293851 
    ## Run 310 stress 0.3330499 
    ## Run 311 stress 0.335524 
    ## Run 312 stress 0.3457041 
    ## Run 313 stress 0.3287222 
    ## Run 314 stress 0.3336231 
    ## Run 315 stress 0.3325573 
    ## Run 316 stress 0.3499293 
    ## Run 317 stress 0.3482621 
    ## Run 318 stress 0.3349713 
    ## Run 319 stress 0.3351398 
    ## Run 320 stress 0.3411715 
    ## Run 321 stress 0.3351298 
    ## Run 322 stress 0.3269736 
    ## Run 323 stress 0.3409906 
    ## Run 324 stress 0.3563198 
    ## Run 325 stress 0.3433565 
    ## Run 326 stress 0.3456568 
    ## Run 327 stress 0.3302223 
    ## Run 328 stress 0.3310237 
    ## Run 329 stress 0.3297976 
    ## Run 330 stress 0.3414431 
    ## Run 331 stress 0.3363171 
    ## Run 332 stress 0.3435569 
    ## Run 333 stress 0.3312042 
    ## Run 334 stress 0.3663518 
    ## Run 335 stress 0.3497303 
    ## Run 336 stress 0.3301354 
    ## Run 337 stress 0.3312331 
    ## Run 338 stress 0.3422868 
    ## Run 339 stress 0.3418134 
    ## Run 340 stress 0.3302446 
    ## Run 341 stress 0.3345645 
    ## Run 342 stress 0.3348981 
    ## Run 343 stress 0.3374553 
    ## Run 344 stress 0.3274056 
    ## Run 345 stress 0.3403228 
    ## Run 346 stress 0.3349865 
    ## Run 347 stress 0.3400151 
    ## Run 348 stress 0.3475384 
    ## Run 349 stress 0.3338678 
    ## Run 350 stress 0.3267563 
    ## Run 351 stress 0.3461117 
    ## Run 352 stress 0.345572 
    ## Run 353 stress 0.3466465 
    ## Run 354 stress 0.330855 
    ## Run 355 stress 0.3327571 
    ## Run 356 stress 0.3315928 
    ## Run 357 stress 0.3252342 
    ## Run 358 stress 0.3350718 
    ## Run 359 stress 0.3302704 
    ## Run 360 stress 0.3303052 
    ## Run 361 stress 0.3324332 
    ## Run 362 stress 0.3390498 
    ## Run 363 stress 0.3421688 
    ## Run 364 stress 0.3340562 
    ## Run 365 stress 0.3327665 
    ## Run 366 stress 0.3455103 
    ## Run 367 stress 0.3773331 
    ## Run 368 stress 0.3427755 
    ## Run 369 stress 0.3493954 
    ## Run 370 stress 0.3357085 
    ## Run 371 stress 0.3500098 
    ## Run 372 stress 0.3454542 
    ## Run 373 stress 0.3317517 
    ## Run 374 stress 0.3365048 
    ## Run 375 stress 0.332401 
    ## Run 376 stress 0.3488461 
    ## Run 377 stress 0.3290485 
    ## Run 378 stress 0.3309745 
    ## Run 379 stress 0.3502059 
    ## Run 380 stress 0.3390544 
    ## Run 381 stress 0.3436483 
    ## Run 382 stress 0.3350292 
    ## Run 383 stress 0.3458465 
    ## Run 384 stress 0.3453806 
    ## Run 385 stress 0.3313799 
    ## Run 386 stress 0.3320206 
    ## Run 387 stress 0.3381721 
    ## Run 388 stress 0.3418935 
    ## Run 389 stress 0.3417603 
    ## Run 390 stress 0.3428683 
    ## Run 391 stress 0.3435401 
    ## Run 392 stress 0.3436807 
    ## Run 393 stress 0.3346925 
    ## Run 394 stress 0.3413961 
    ## Run 395 stress 0.3423716 
    ## Run 396 stress 0.3399107 
    ## Run 397 stress 0.3432802 
    ## Run 398 stress 0.343045 
    ## Run 399 stress 0.3270133 
    ## Run 400 stress 0.3262088 
    ## Run 401 stress 0.3297483 
    ## Run 402 stress 0.3362434 
    ## Run 403 stress 0.3335461 
    ## Run 404 stress 0.3362154 
    ## Run 405 stress 0.3427537 
    ## Run 406 stress 0.3287421 
    ## Run 407 stress 0.3290789 
    ## Run 408 stress 0.3350235 
    ## Run 409 stress 0.3441227 
    ## Run 410 stress 0.3308241 
    ## Run 411 stress 0.3299899 
    ## Run 412 stress 0.3454151 
    ## Run 413 stress 0.3379834 
    ## Run 414 stress 0.3319055 
    ## Run 415 stress 0.3372705 
    ## Run 416 stress 0.3378383 
    ## Run 417 stress 0.3343525 
    ## Run 418 stress 0.3282043 
    ## Run 419 stress 0.338415 
    ## Run 420 stress 0.3322439 
    ## Run 421 stress 0.3407747 
    ## Run 422 stress 0.3317464 
    ## Run 423 stress 0.3270908 
    ## Run 424 stress 0.3421173 
    ## Run 425 stress 0.338523 
    ## Run 426 stress 0.3275845 
    ## Run 427 stress 0.3271147 
    ## Run 428 stress 0.3334286 
    ## Run 429 stress 0.3347997 
    ## Run 430 stress 0.3326347 
    ## Run 431 stress 0.3374434 
    ## Run 432 stress 0.3289921 
    ## Run 433 stress 0.3471921 
    ## Run 434 stress 0.3491603 
    ## Run 435 stress 0.3416722 
    ## Run 436 stress 0.332399 
    ## Run 437 stress 0.3326869 
    ## Run 438 stress 0.3392309 
    ## Run 439 stress 0.3311284 
    ## Run 440 stress 0.3428153 
    ## Run 441 stress 0.3438833 
    ## Run 442 stress 0.3385699 
    ## Run 443 stress 0.345881 
    ## Run 444 stress 0.3355273 
    ## Run 445 stress 0.3364274 
    ## Run 446 stress 0.3451392 
    ## Run 447 stress 0.3331102 
    ## Run 448 stress 0.3354258 
    ## Run 449 stress 0.3312324 
    ## Run 450 stress 0.3326956 
    ## Run 451 stress 0.3396123 
    ## Run 452 stress 0.3411919 
    ## Run 453 stress 0.3348609 
    ## Run 454 stress 0.3387228 
    ## Run 455 stress 0.3280932 
    ## Run 456 stress 0.3400811 
    ## Run 457 stress 0.3807935 
    ## Run 458 stress 0.3397543 
    ## Run 459 stress 0.3341873 
    ## Run 460 stress 0.3502161 
    ## Run 461 stress 0.3464482 
    ## Run 462 stress 0.354423 
    ## Run 463 stress 0.3306984 
    ## Run 464 stress 0.3428242 
    ## Run 465 stress 0.3306921 
    ## Run 466 stress 0.3432298 
    ## Run 467 stress 0.3334312 
    ## Run 468 stress 0.3338607 
    ## Run 469 stress 0.3482835 
    ## Run 470 stress 0.3335328 
    ## Run 471 stress 0.3392754 
    ## Run 472 stress 0.3459787 
    ## Run 473 stress 0.341619 
    ## Run 474 stress 0.3338528 
    ## Run 475 stress 0.3298204 
    ## Run 476 stress 0.3454611 
    ## Run 477 stress 0.3388363 
    ## Run 478 stress 0.3291826 
    ## Run 479 stress 0.3463214 
    ## Run 480 stress 0.3408457 
    ## Run 481 stress 0.3383454 
    ## Run 482 stress 0.3377149 
    ## Run 483 stress 0.3399514 
    ## Run 484 stress 0.3432207 
    ## Run 485 stress 0.3376051 
    ## Run 486 stress 0.34571 
    ## Run 487 stress 0.330549 
    ## Run 488 stress 0.3343805 
    ## Run 489 stress 0.3429129 
    ## Run 490 stress 0.3440306 
    ## Run 491 stress 0.3288823 
    ## Run 492 stress 0.3324104 
    ## Run 493 stress 0.3338089 
    ## Run 494 stress 0.3350625 
    ## Run 495 stress 0.3308107 
    ## Run 496 stress 0.3342716 
    ## Run 497 stress 0.3345112 
    ## Run 498 stress 0.3308616 
    ## Run 499 stress 0.3276451 
    ## Run 500 stress 0.3411421 
    ## *** Best solution was not repeated -- monoMDS stopping criteria:
    ##      1: no. of iterations >= maxit
    ##    499: stress ratio > sratmax

``` r
SD_beta_ric_NMDS <- metaMDS(SD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.01879552 
    ## Run 1 stress 0.02486143 
    ## Run 2 stress 0.01879595 
    ## ... Procrustes: rmse 0.001313348  max resid 0.002704535 
    ## ... Similar to previous best
    ## Run 3 stress 0.02486108 
    ## Run 4 stress 0.01879582 
    ## ... Procrustes: rmse 0.001274553  max resid 0.002624894 
    ## ... Similar to previous best
    ## Run 5 stress 0.01879582 
    ## ... Procrustes: rmse 0.0001383358  max resid 0.0002837959 
    ## ... Similar to previous best
    ## Run 6 stress 0.02509456 
    ## Run 7 stress 0.01879537 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001047346  max resid 0.002157328 
    ## ... Similar to previous best
    ## Run 8 stress 0.0250945 
    ## Run 9 stress 0.01910389 
    ## ... Procrustes: rmse 0.01333472  max resid 0.02729765 
    ## Run 10 stress 0.02486148 
    ## Run 11 stress 0.01953612 
    ## Run 12 stress 0.018866 
    ## ... Procrustes: rmse 0.00515277  max resid 0.01028918 
    ## Run 13 stress 0.01879593 
    ## ... Procrustes: rmse 0.0002573139  max resid 0.000529616 
    ## ... Similar to previous best
    ## Run 14 stress 0.018832 
    ## ... Procrustes: rmse 0.002251981  max resid 0.005044072 
    ## ... Similar to previous best
    ## Run 15 stress 0.01879571 
    ## ... Procrustes: rmse 0.0001771111  max resid 0.0003648178 
    ## ... Similar to previous best
    ## Run 16 stress 0.02476741 
    ## Run 17 stress 0.01879583 
    ## ... Procrustes: rmse 0.0002198909  max resid 0.0004526537 
    ## ... Similar to previous best
    ## Run 18 stress 0.01879564 
    ## ... Procrustes: rmse 0.001106971  max resid 0.002281947 
    ## ... Similar to previous best
    ## Run 19 stress 0.02495789 
    ## Run 20 stress 0.01879581 
    ## ... Procrustes: rmse 0.0002142722  max resid 0.0004410716 
    ## ... Similar to previous best
    ## Run 21 stress 0.02492192 
    ## Run 22 stress 0.01879596 
    ## ... Procrustes: rmse 0.001229483  max resid 0.002533603 
    ## ... Similar to previous best
    ## Run 23 stress 0.01879554 
    ## ... Procrustes: rmse 9.432754e-05  max resid 0.0001942769 
    ## ... Similar to previous best
    ## Run 24 stress 0.02492195 
    ## Run 25 stress 0.02509483 
    ## Run 26 stress 0.02509457 
    ## Run 27 stress 0.01879558 
    ## ... Procrustes: rmse 0.0001141523  max resid 0.0002350156 
    ## ... Similar to previous best
    ## Run 28 stress 0.02004735 
    ## Run 29 stress 0.02509478 
    ## Run 30 stress 0.02509467 
    ## Run 31 stress 0.01879566 
    ## ... Procrustes: rmse 0.001116882  max resid 0.002302341 
    ## ... Similar to previous best
    ## Run 32 stress 0.01883352 
    ## ... Procrustes: rmse 0.002388353  max resid 0.005136792 
    ## ... Similar to previous best
    ## Run 33 stress 0.0187957 
    ## ... Procrustes: rmse 0.0001738022  max resid 0.0003580179 
    ## ... Similar to previous best
    ## Run 34 stress 0.01883508 
    ## ... Procrustes: rmse 0.002685489  max resid 0.005375422 
    ## ... Similar to previous best
    ## Run 35 stress 0.02509479 
    ## Run 36 stress 0.01879563 
    ## ... Procrustes: rmse 0.001101233  max resid 0.002270182 
    ## ... Similar to previous best
    ## Run 37 stress 0.01879558 
    ## ... Procrustes: rmse 0.001076855  max resid 0.002220066 
    ## ... Similar to previous best
    ## Run 38 stress 0.02492164 
    ## Run 39 stress 0.02520894 
    ## Run 40 stress 0.02492183 
    ## Run 41 stress 0.02492177 
    ## Run 42 stress 0.2038604 
    ## Run 43 stress 0.01883627 
    ## ... Procrustes: rmse 0.002815764  max resid 0.005463227 
    ## ... Similar to previous best
    ## Run 44 stress 0.02492183 
    ## Run 45 stress 0.01879663 
    ## ... Procrustes: rmse 0.0005086965  max resid 0.001047782 
    ## ... Similar to previous best
    ## Run 46 stress 0.01879577 
    ## ... Procrustes: rmse 0.0002046203  max resid 0.000421505 
    ## ... Similar to previous best
    ## Run 47 stress 0.01892345 
    ## ... Procrustes: rmse 0.008107008  max resid 0.0166385 
    ## Run 48 stress 0.01888234 
    ## ... Procrustes: rmse 0.006065956  max resid 0.01220041 
    ## Run 49 stress 0.01879601 
    ## ... Procrustes: rmse 0.0002879677  max resid 0.0005927772 
    ## ... Similar to previous best
    ## Run 50 stress 0.01884211 
    ## ... Procrustes: rmse 0.003247945  max resid 0.006238313 
    ## ... Similar to previous best
    ## Run 51 stress 0.02492191 
    ## Run 52 stress 0.02492204 
    ## Run 53 stress 0.02509456 
    ## Run 54 stress 0.02492434 
    ## Run 55 stress 0.02509639 
    ## Run 56 stress 0.01886377 
    ## ... Procrustes: rmse 0.00501243  max resid 0.009991576 
    ## Run 57 stress 0.02492186 
    ## Run 58 stress 0.02520909 
    ## Run 59 stress 0.2241378 
    ## Run 60 stress 0.02486134 
    ## Run 61 stress 0.01879558 
    ## ... Procrustes: rmse 0.0001169616  max resid 0.0002409908 
    ## ... Similar to previous best
    ## Run 62 stress 0.02486128 
    ## Run 63 stress 0.01879733 
    ## ... Procrustes: rmse 0.000704237  max resid 0.001450683 
    ## ... Similar to previous best
    ## Run 64 stress 0.01883379 
    ## ... Procrustes: rmse 0.00253042  max resid 0.005254018 
    ## ... Similar to previous best
    ## Run 65 stress 0.02492349 
    ## Run 66 stress 0.02520891 
    ## Run 67 stress 0.01879578 
    ## ... Procrustes: rmse 0.000209354  max resid 0.000431215 
    ## ... Similar to previous best
    ## Run 68 stress 0.02509442 
    ## Run 69 stress 0.01879556 
    ## ... Procrustes: rmse 0.0001046583  max resid 0.0002155614 
    ## ... Similar to previous best
    ## Run 70 stress 0.02486134 
    ## Run 71 stress 0.0250945 
    ## Run 72 stress 0.01891054 
    ## ... Procrustes: rmse 0.007713847  max resid 0.01583634 
    ## Run 73 stress 0.02486129 
    ## Run 74 stress 0.02520899 
    ## Run 75 stress 0.01879818 
    ## ... Procrustes: rmse 0.0009036017  max resid 0.00186156 
    ## ... Similar to previous best
    ## Run 76 stress 0.01879593 
    ## ... Procrustes: rmse 0.00123483  max resid 0.002544697 
    ## ... Similar to previous best
    ## Run 77 stress 0.01882684 
    ## ... Procrustes: rmse 0.00391796  max resid 0.008062698 
    ## ... Similar to previous best
    ## Run 78 stress 0.01913827 
    ## ... Procrustes: rmse 0.01341698  max resid 0.02745591 
    ## Run 79 stress 0.01884116 
    ## ... Procrustes: rmse 0.004820686  max resid 0.009917308 
    ## ... Similar to previous best
    ## Run 80 stress 0.02486136 
    ## Run 81 stress 0.01880538 
    ## ... Procrustes: rmse 0.002020112  max resid 0.004162566 
    ## ... Similar to previous best
    ## Run 82 stress 0.02486134 
    ## Run 83 stress 0.02509448 
    ## Run 84 stress 0.02509454 
    ## Run 85 stress 0.01879544 
    ## ... Procrustes: rmse 4.261626e-05  max resid 8.756341e-05 
    ## ... Similar to previous best
    ## Run 86 stress 0.0248613 
    ## Run 87 stress 0.01879611 
    ## ... Procrustes: rmse 0.0003224413  max resid 0.0006637586 
    ## ... Similar to previous best
    ## Run 88 stress 0.2039098 
    ## Run 89 stress 0.02492161 
    ## Run 90 stress 0.01885544 
    ## ... Procrustes: rmse 0.00559119  max resid 0.01149728 
    ## Run 91 stress 0.01879502 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0005325658  max resid 0.001084297 
    ## ... Similar to previous best
    ## Run 92 stress 0.02492163 
    ## Run 93 stress 0.01879565 
    ## ... Procrustes: rmse 0.0006745165  max resid 0.00137556 
    ## ... Similar to previous best
    ## Run 94 stress 0.02509461 
    ## Run 95 stress 0.02509452 
    ## Run 96 stress 0.02509435 
    ## Run 97 stress 0.02486134 
    ## Run 98 stress 0.018796 
    ## ... Procrustes: rmse 0.0007307048  max resid 0.001512917 
    ## ... Similar to previous best
    ## Run 99 stress 0.2039098 
    ## Run 100 stress 0.02509469 
    ## Run 101 stress 0.02509483 
    ## Run 102 stress 0.0248613 
    ## Run 103 stress 0.0250944 
    ## Run 104 stress 0.02520881 
    ## Run 105 stress 0.347268 
    ## Run 106 stress 0.01879582 
    ## ... Procrustes: rmse 0.0007541243  max resid 0.00154074 
    ## ... Similar to previous best
    ## Run 107 stress 0.01883644 
    ## ... Procrustes: rmse 0.003327658  max resid 0.006427727 
    ## ... Similar to previous best
    ## Run 108 stress 0.02509442 
    ## Run 109 stress 0.02492192 
    ## Run 110 stress 0.01879616 
    ## ... Procrustes: rmse 0.0008717438  max resid 0.001783348 
    ## ... Similar to previous best
    ## Run 111 stress 0.01879586 
    ## ... Procrustes: rmse 0.0006745774  max resid 0.001397439 
    ## ... Similar to previous best
    ## Run 112 stress 0.02509435 
    ## Run 113 stress 0.01879567 
    ## ... Procrustes: rmse 0.000693341  max resid 0.001415386 
    ## ... Similar to previous best
    ## Run 114 stress 0.02492167 
    ## Run 115 stress 0.02486132 
    ## Run 116 stress 0.02520881 
    ## Run 117 stress 0.01879575 
    ## ... Procrustes: rmse 0.0007269312  max resid 0.001484753 
    ## ... Similar to previous best
    ## Run 118 stress 0.02492175 
    ## Run 119 stress 0.01879576 
    ## ... Procrustes: rmse 0.0007345259  max resid 0.00150043 
    ## ... Similar to previous best
    ## Run 120 stress 0.02492164 
    ## Run 121 stress 0.02509444 
    ## Run 122 stress 0.01882333 
    ## ... Procrustes: rmse 0.004194015  max resid 0.00861996 
    ## ... Similar to previous best
    ## Run 123 stress 0.02492194 
    ## Run 124 stress 0.02509463 
    ## Run 125 stress 0.01879593 
    ## ... Procrustes: rmse 0.0008021901  max resid 0.001639892 
    ## ... Similar to previous best
    ## Run 126 stress 0.02509474 
    ## Run 127 stress 0.02492211 
    ## Run 128 stress 0.02509452 
    ## Run 129 stress 0.02486114 
    ## Run 130 stress 0.02045619 
    ## Run 131 stress 0.01879546 
    ## ... Procrustes: rmse 0.0004848302  max resid 0.001006595 
    ## ... Similar to previous best
    ## Run 132 stress 0.0248612 
    ## Run 133 stress 0.01879599 
    ## ... Procrustes: rmse 0.0008117954  max resid 0.00165978 
    ## ... Similar to previous best
    ## Run 134 stress 0.01879805 
    ## ... Procrustes: rmse 0.001375339  max resid 0.002822447 
    ## ... Similar to previous best
    ## Run 135 stress 0.01883687 
    ## ... Procrustes: rmse 0.003375703  max resid 0.006530586 
    ## ... Similar to previous best
    ## Run 136 stress 0.01879578 
    ## ... Procrustes: rmse 0.0007416833  max resid 0.001515225 
    ## ... Similar to previous best
    ## Run 137 stress 0.01879596 
    ## ... Procrustes: rmse 0.0008158448  max resid 0.001668067 
    ## ... Similar to previous best
    ## Run 138 stress 0.02509455 
    ## Run 139 stress 0.01879551 
    ## ... Procrustes: rmse 0.0005128767  max resid 0.001064595 
    ## ... Similar to previous best
    ## Run 140 stress 0.02509473 
    ## Run 141 stress 0.01879546 
    ## ... Procrustes: rmse 0.0005810771  max resid 0.00118379 
    ## ... Similar to previous best
    ## Run 142 stress 0.01879545 
    ## ... Procrustes: rmse 0.0004811944  max resid 0.0009992193 
    ## ... Similar to previous best
    ## Run 143 stress 0.01884165 
    ## ... Procrustes: rmse 0.003862465  max resid 0.007569469 
    ## ... Similar to previous best
    ## Run 144 stress 0.0248612 
    ## Run 145 stress 0.02509472 
    ## Run 146 stress 0.0188686 
    ## ... Procrustes: rmse 0.006338391  max resid 0.01300709 
    ## Run 147 stress 0.01879562 
    ## ... Procrustes: rmse 0.0006628821  max resid 0.001352578 
    ## ... Similar to previous best
    ## Run 148 stress 0.02509467 
    ## Run 149 stress 0.01879561 
    ## ... Procrustes: rmse 0.0006539023  max resid 0.001334027 
    ## ... Similar to previous best
    ## Run 150 stress 0.02520904 
    ## Run 151 stress 0.01879583 
    ## ... Procrustes: rmse 0.0007530987  max resid 0.001538762 
    ## ... Similar to previous best
    ## Run 152 stress 0.0188507 
    ## ... Procrustes: rmse 0.005050325  max resid 0.01038204 
    ## Run 153 stress 0.02509487 
    ## Run 154 stress 0.02509444 
    ## Run 155 stress 0.01879563 
    ## ... Procrustes: rmse 0.0006719801  max resid 0.001371428 
    ## ... Similar to previous best
    ## Run 156 stress 0.02492168 
    ## Run 157 stress 0.0187999 
    ## ... Procrustes: rmse 0.00171983  max resid 0.003531786 
    ## ... Similar to previous best
    ## Run 158 stress 0.01881147 
    ## ... Procrustes: rmse 0.00319222  max resid 0.006562032 
    ## ... Similar to previous best
    ## Run 159 stress 0.01879582 
    ## ... Procrustes: rmse 0.0006572191  max resid 0.001361735 
    ## ... Similar to previous best
    ## Run 160 stress 0.0250944 
    ## Run 161 stress 0.01879582 
    ## ... Procrustes: rmse 0.0007575763  max resid 0.001547937 
    ## ... Similar to previous best
    ## Run 162 stress 0.0250944 
    ## Run 163 stress 0.02509482 
    ## Run 164 stress 0.02509493 
    ## Run 165 stress 0.01885348 
    ## ... Procrustes: rmse 0.004768128  max resid 0.009483427 
    ## ... Similar to previous best
    ## Run 166 stress 0.01879574 
    ## ... Procrustes: rmse 0.0007211953  max resid 0.00147296 
    ## ... Similar to previous best
    ## Run 167 stress 0.01879582 
    ## ... Procrustes: rmse 0.0006483113  max resid 0.001343444 
    ## ... Similar to previous best
    ## Run 168 stress 0.01879593 
    ## ... Procrustes: rmse 0.0008021688  max resid 0.001639897 
    ## ... Similar to previous best
    ## Run 169 stress 0.02520908 
    ## Run 170 stress 0.0248611 
    ## Run 171 stress 0.01879562 
    ## ... Procrustes: rmse 0.0006653942  max resid 0.0013578 
    ## ... Similar to previous best
    ## Run 172 stress 0.02509464 
    ## Run 173 stress 0.0187958 
    ## ... Procrustes: rmse 0.0007507751  max resid 0.001533921 
    ## ... Similar to previous best
    ## Run 174 stress 0.01879843 
    ## ... Procrustes: rmse 0.001482246  max resid 0.003041657 
    ## ... Similar to previous best
    ## Run 175 stress 0.02509452 
    ## Run 176 stress 0.02486126 
    ## Run 177 stress 0.01879695 
    ## ... Procrustes: rmse 0.001134311  max resid 0.002324616 
    ## ... Similar to previous best
    ## Run 178 stress 0.02492169 
    ## Run 179 stress 0.0250945 
    ## Run 180 stress 0.0188439 
    ## ... Procrustes: rmse 0.003993202  max resid 0.007845272 
    ## ... Similar to previous best
    ## Run 181 stress 0.02492187 
    ## Run 182 stress 0.01879584 
    ## ... Procrustes: rmse 0.000766423  max resid 0.00156615 
    ## ... Similar to previous best
    ## Run 183 stress 0.01879594 
    ## ... Procrustes: rmse 0.0008052973  max resid 0.001646249 
    ## ... Similar to previous best
    ## Run 184 stress 0.02486131 
    ## Run 185 stress 0.02476845 
    ## Run 186 stress 0.02492164 
    ## Run 187 stress 0.01884943 
    ## ... Procrustes: rmse 0.004476347  max resid 0.00886984 
    ## ... Similar to previous best
    ## Run 188 stress 0.01880333 
    ## ... Procrustes: rmse 0.002264347  max resid 0.004652792 
    ## ... Similar to previous best
    ## Run 189 stress 0.01883075 
    ## ... Procrustes: rmse 0.002556893  max resid 0.00495186 
    ## ... Similar to previous best
    ## Run 190 stress 0.02509476 
    ## Run 191 stress 0.02492173 
    ## Run 192 stress 0.02509475 
    ## Run 193 stress 0.01896147 
    ## ... Procrustes: rmse 0.009768229  max resid 0.01987577 
    ## Run 194 stress 0.02492183 
    ## Run 195 stress 0.02486113 
    ## Run 196 stress 0.01883565 
    ## ... Procrustes: rmse 0.003238439  max resid 0.006237478 
    ## ... Similar to previous best
    ## Run 197 stress 0.02492189 
    ## Run 198 stress 0.02486137 
    ## Run 199 stress 0.01879617 
    ## ... Procrustes: rmse 0.000779167  max resid 0.001612565 
    ## ... Similar to previous best
    ## Run 200 stress 0.01901579 
    ## ... Procrustes: rmse 0.0114539  max resid 0.0234719 
    ## Run 201 stress 0.01879599 
    ## ... Procrustes: rmse 0.0008149892  max resid 0.001666375 
    ## ... Similar to previous best
    ## Run 202 stress 0.01879608 
    ## ... Procrustes: rmse 0.0008430493  max resid 0.001723907 
    ## ... Similar to previous best
    ## Run 203 stress 0.0248612 
    ## Run 204 stress 0.01879565 
    ## ... Procrustes: rmse 0.0006805834  max resid 0.001389169 
    ## ... Similar to previous best
    ## Run 205 stress 0.02486132 
    ## Run 206 stress 0.02486109 
    ## Run 207 stress 0.01882202 
    ## ... Procrustes: rmse 0.004098836  max resid 0.008424843 
    ## ... Similar to previous best
    ## Run 208 stress 0.01879833 
    ## ... Procrustes: rmse 0.001465713  max resid 0.00300748 
    ## ... Similar to previous best
    ## Run 209 stress 0.01879597 
    ## ... Procrustes: rmse 0.0008182417  max resid 0.001672982 
    ## ... Similar to previous best
    ## Run 210 stress 0.02509435 
    ## Run 211 stress 0.01879888 
    ## ... Procrustes: rmse 0.001577775  max resid 0.003238655 
    ## ... Similar to previous best
    ## Run 212 stress 0.01884862 
    ## ... Procrustes: rmse 0.004464495  max resid 0.008846076 
    ## ... Similar to previous best
    ## Run 213 stress 0.01884888 
    ## ... Procrustes: rmse 0.004469179  max resid 0.00885482 
    ## ... Similar to previous best
    ## Run 214 stress 0.01906282 
    ## ... Procrustes: rmse 0.01281056  max resid 0.02613508 
    ## Run 215 stress 0.02486125 
    ## Run 216 stress 0.02486134 
    ## Run 217 stress 0.0249217 
    ## Run 218 stress 0.02509461 
    ## Run 219 stress 0.01879587 
    ## ... Procrustes: rmse 0.0007712401  max resid 0.001575998 
    ## ... Similar to previous best
    ## Run 220 stress 0.02509471 
    ## Run 221 stress 0.0187957 
    ## ... Procrustes: rmse 0.0006072135  max resid 0.001258823 
    ## ... Similar to previous best
    ## Run 222 stress 0.0187957 
    ## ... Procrustes: rmse 0.0007031566  max resid 0.001435649 
    ## ... Similar to previous best
    ## Run 223 stress 0.02509464 
    ## Run 224 stress 0.01881056 
    ## ... Procrustes: rmse 0.003101321  max resid 0.006375541 
    ## ... Similar to previous best
    ## Run 225 stress 0.01879826 
    ## ... Procrustes: rmse 0.001452355  max resid 0.00297981 
    ## ... Similar to previous best
    ## Run 226 stress 0.01879573 
    ## ... Procrustes: rmse 0.0006215373  max resid 0.001288321 
    ## ... Similar to previous best
    ## Run 227 stress 0.01880537 
    ## ... Procrustes: rmse 0.002533995  max resid 0.005208001 
    ## ... Similar to previous best
    ## Run 228 stress 0.02509451 
    ## Run 229 stress 0.02520917 
    ## Run 230 stress 0.0249216 
    ## Run 231 stress 0.01881049 
    ## ... Procrustes: rmse 0.003111065  max resid 0.006395804 
    ## ... Similar to previous best
    ## Run 232 stress 0.0187956 
    ## ... Procrustes: rmse 0.0005583719  max resid 0.001158305 
    ## ... Similar to previous best
    ## Run 233 stress 0.01879882 
    ## ... Procrustes: rmse 0.001533528  max resid 0.003146567 
    ## ... Similar to previous best
    ## Run 234 stress 0.02509454 
    ## Run 235 stress 0.01879516 
    ## ... Procrustes: rmse 0.0003689761  max resid 0.0007452324 
    ## ... Similar to previous best
    ## Run 236 stress 0.02492349 
    ## Run 237 stress 0.02509442 
    ## Run 238 stress 0.02492181 
    ## Run 239 stress 0.02509441 
    ## Run 240 stress 0.0187957 
    ## ... Procrustes: rmse 0.0007075849  max resid 0.001444872 
    ## ... Similar to previous best
    ## Run 241 stress 0.02509431 
    ## Run 242 stress 0.01880704 
    ## ... Procrustes: rmse 0.002741943  max resid 0.00563631 
    ## ... Similar to previous best
    ## Run 243 stress 0.01884265 
    ## ... Procrustes: rmse 0.003937027  max resid 0.007730145 
    ## ... Similar to previous best
    ## Run 244 stress 0.01879559 
    ## ... Procrustes: rmse 0.0006555856  max resid 0.001337575 
    ## ... Similar to previous best
    ## Run 245 stress 0.02509439 
    ## Run 246 stress 0.02486136 
    ## Run 247 stress 0.2999566 
    ## Run 248 stress 0.02486139 
    ## Run 249 stress 0.02492173 
    ## Run 250 stress 0.02486111 
    ## Run 251 stress 0.02509471 
    ## Run 252 stress 0.02509427 
    ## Run 253 stress 0.02492176 
    ## Run 254 stress 0.02492164 
    ## Run 255 stress 0.02492186 
    ## Run 256 stress 0.02486103 
    ## Run 257 stress 0.01886499 
    ## ... Procrustes: rmse 0.00658288  max resid 0.01352046 
    ## Run 258 stress 0.01879574 
    ## ... Procrustes: rmse 0.0006226353  max resid 0.001290557 
    ## ... Similar to previous best
    ## Run 259 stress 0.02492182 
    ## Run 260 stress 0.02486129 
    ## Run 261 stress 0.2039098 
    ## Run 262 stress 0.01920782 
    ## ... Procrustes: rmse 0.01376091  max resid 0.02808342 
    ## Run 263 stress 0.02492177 
    ## Run 264 stress 0.01879579 
    ## ... Procrustes: rmse 0.0007458237  max resid 0.001523711 
    ## ... Similar to previous best
    ## Run 265 stress 0.2038604 
    ## Run 266 stress 0.02486132 
    ## Run 267 stress 0.02492192 
    ## Run 268 stress 0.0190212 
    ## ... Procrustes: rmse 0.01190528  max resid 0.02439134 
    ## Run 269 stress 0.0190049 
    ## ... Procrustes: rmse 0.01116778  max resid 0.02275721 
    ## Run 270 stress 0.02486103 
    ## Run 271 stress 0.01879575 
    ## ... Procrustes: rmse 0.00072707  max resid 0.001485072 
    ## ... Similar to previous best
    ## Run 272 stress 0.0188498 
    ## ... Procrustes: rmse 0.005834349  max resid 0.01198522 
    ## Run 273 stress 0.02492194 
    ## Run 274 stress 0.01887561 
    ## ... Procrustes: rmse 0.007076037  max resid 0.01453079 
    ## Run 275 stress 0.02492198 
    ## Run 276 stress 0.0250946 
    ## Run 277 stress 0.01879556 
    ## ... Procrustes: rmse 0.0006344804  max resid 0.001294085 
    ## ... Similar to previous best
    ## Run 278 stress 0.01879565 
    ## ... Procrustes: rmse 0.000681918  max resid 0.001391923 
    ## ... Similar to previous best
    ## Run 279 stress 0.02492188 
    ## Run 280 stress 0.01879596 
    ## ... Procrustes: rmse 0.0008035146  max resid 0.001642464 
    ## ... Similar to previous best
    ## Run 281 stress 0.02486096 
    ## Run 282 stress 0.01879595 
    ## ... Procrustes: rmse 0.0008124553  max resid 0.001661063 
    ## ... Similar to previous best
    ## Run 283 stress 0.02509447 
    ## Run 284 stress 0.01881244 
    ## ... Procrustes: rmse 0.00329745  max resid 0.006779371 
    ## ... Similar to previous best
    ## Run 285 stress 0.02509465 
    ## Run 286 stress 0.02509474 
    ## Run 287 stress 0.02493106 
    ## Run 288 stress 0.01880614 
    ## ... Procrustes: rmse 0.002644478  max resid 0.005436178 
    ## ... Similar to previous best
    ## Run 289 stress 0.01879552 
    ## ... Procrustes: rmse 0.0006186197  max resid 0.001261322 
    ## ... Similar to previous best
    ## Run 290 stress 0.02499223 
    ## Run 291 stress 0.01879543 
    ## ... Procrustes: rmse 0.0005664312  max resid 0.001153611 
    ## ... Similar to previous best
    ## Run 292 stress 0.0187956 
    ## ... Procrustes: rmse 0.0006499758  max resid 0.001325841 
    ## ... Similar to previous best
    ## Run 293 stress 0.01879579 
    ## ... Procrustes: rmse 0.0007325485  max resid 0.001496384 
    ## ... Similar to previous best
    ## Run 294 stress 0.02486135 
    ## Run 295 stress 0.2262188 
    ## Run 296 stress 0.01879589 
    ## ... Procrustes: rmse 0.0006869083  max resid 0.001422809 
    ## ... Similar to previous best
    ## Run 297 stress 0.0188127 
    ## ... Procrustes: rmse 0.003323684  max resid 0.006833099 
    ## ... Similar to previous best
    ## Run 298 stress 0.02486125 
    ## Run 299 stress 0.02509448 
    ## Run 300 stress 0.268335 
    ## Run 301 stress 0.01879566 
    ## ... Procrustes: rmse 0.0006886786  max resid 0.001405883 
    ## ... Similar to previous best
    ## Run 302 stress 0.01879597 
    ## ... Procrustes: rmse 0.0008043344  max resid 0.001644144 
    ## ... Similar to previous best
    ## Run 303 stress 0.02509478 
    ## Run 304 stress 0.02509428 
    ## Run 305 stress 0.0188337 
    ## ... Procrustes: rmse 0.0030005  max resid 0.005727635 
    ## ... Similar to previous best
    ## Run 306 stress 0.02509446 
    ## Run 307 stress 0.01879581 
    ## ... Procrustes: rmse 0.000645006  max resid 0.001336475 
    ## ... Similar to previous best
    ## Run 308 stress 0.02492183 
    ## Run 309 stress 0.01879598 
    ## ... Procrustes: rmse 0.0008048205  max resid 0.001645398 
    ## ... Similar to previous best
    ## Run 310 stress 0.02486113 
    ## Run 311 stress 0.0248614 
    ## Run 312 stress 0.02486122 
    ## Run 313 stress 0.02509466 
    ## Run 314 stress 0.02492195 
    ## Run 315 stress 0.0187994 
    ## ... Procrustes: rmse 0.001666178  max resid 0.003431261 
    ## ... Similar to previous best
    ## Run 316 stress 0.01879579 
    ## ... Procrustes: rmse 0.0006291073  max resid 0.001304733 
    ## ... Similar to previous best
    ## Run 317 stress 0.02509457 
    ## Run 318 stress 0.3747865 
    ## Run 319 stress 0.01879549 
    ## ... Procrustes: rmse 0.0006000488  max resid 0.001222987 
    ## ... Similar to previous best
    ## Run 320 stress 0.02492205 
    ## Run 321 stress 0.02509482 
    ## Run 322 stress 0.01879596 
    ## ... Procrustes: rmse 0.0007950785  max resid 0.001625326 
    ## ... Similar to previous best
    ## Run 323 stress 0.01879574 
    ## ... Procrustes: rmse 0.0007200393  max resid 0.00147058 
    ## ... Similar to previous best
    ## Run 324 stress 0.01879584 
    ## ... Procrustes: rmse 0.0007453032  max resid 0.001522405 
    ## ... Similar to previous best
    ## Run 325 stress 0.02509474 
    ## Run 326 stress 0.01880197 
    ## ... Procrustes: rmse 0.00126177  max resid 0.002579144 
    ## ... Similar to previous best
    ## Run 327 stress 0.01879566 
    ## ... Procrustes: rmse 0.0006891129  max resid 0.001406773 
    ## ... Similar to previous best
    ## Run 328 stress 0.02509447 
    ## Run 329 stress 0.01881156 
    ## ... Procrustes: rmse 0.003167473  max resid 0.006512864 
    ## ... Similar to previous best
    ## Run 330 stress 0.01879564 
    ## ... Procrustes: rmse 0.0006783358  max resid 0.001384565 
    ## ... Similar to previous best
    ## Run 331 stress 0.02492194 
    ## Run 332 stress 0.01879585 
    ## ... Procrustes: rmse 0.0007623287  max resid 0.001557786 
    ## ... Similar to previous best
    ## Run 333 stress 0.02492985 
    ## Run 334 stress 0.02520878 
    ## Run 335 stress 0.01879574 
    ## ... Procrustes: rmse 0.000725627  max resid 0.001482089 
    ## ... Similar to previous best
    ## Run 336 stress 0.01880184 
    ## ... Procrustes: rmse 0.002081225  max resid 0.004275958 
    ## ... Similar to previous best
    ## Run 337 stress 0.02509454 
    ## Run 338 stress 0.01879586 
    ## ... Procrustes: rmse 0.0007744459  max resid 0.001582734 
    ## ... Similar to previous best
    ## Run 339 stress 0.02509462 
    ## Run 340 stress 0.0188589 
    ## ... Procrustes: rmse 0.00498381  max resid 0.009933672 
    ## ... Similar to previous best
    ## Run 341 stress 0.3747789 
    ## Run 342 stress 0.01888489 
    ## ... Procrustes: rmse 0.007472293  max resid 0.01534058 
    ## Run 343 stress 0.02509463 
    ## Run 344 stress 0.01883085 
    ## ... Procrustes: rmse 0.002586029  max resid 0.004973756 
    ## ... Similar to previous best
    ## Run 345 stress 0.01894093 
    ## ... Procrustes: rmse 0.008886668  max resid 0.01805153 
    ## Run 346 stress 0.02509444 
    ## Run 347 stress 0.02520872 
    ## Run 348 stress 0.01879582 
    ## ... Procrustes: rmse 0.0007560482  max resid 0.001544825 
    ## ... Similar to previous best
    ## Run 349 stress 0.02509469 
    ## Run 350 stress 0.02492854 
    ## Run 351 stress 0.01887997 
    ## ... Procrustes: rmse 0.00726797  max resid 0.01492356 
    ## Run 352 stress 0.01880527 
    ## ... Procrustes: rmse 0.002537197  max resid 0.005216173 
    ## ... Similar to previous best
    ## Run 353 stress 0.01879975 
    ## ... Procrustes: rmse 0.001742517  max resid 0.003577909 
    ## ... Similar to previous best
    ## Run 354 stress 0.02492199 
    ## Run 355 stress 0.02509469 
    ## Run 356 stress 0.02492173 
    ## Run 357 stress 0.02509436 
    ## Run 358 stress 0.3581175 
    ## Run 359 stress 0.02486113 
    ## Run 360 stress 0.01879548 
    ## ... Procrustes: rmse 0.0005979236  max resid 0.001218597 
    ## ... Similar to previous best
    ## Run 361 stress 0.02509457 
    ## Run 362 stress 0.01885022 
    ## ... Procrustes: rmse 0.004587254  max resid 0.009105459 
    ## ... Similar to previous best
    ## Run 363 stress 0.01881503 
    ## ... Procrustes: rmse 0.003534618  max resid 0.007265888 
    ## ... Similar to previous best
    ## Run 364 stress 0.02509471 
    ## Run 365 stress 0.01879561 
    ## ... Procrustes: rmse 0.0005590798  max resid 0.001159815 
    ## ... Similar to previous best
    ## Run 366 stress 0.01879913 
    ## ... Procrustes: rmse 0.00162432  max resid 0.003334629 
    ## ... Similar to previous best
    ## Run 367 stress 0.02509458 
    ## Run 368 stress 0.02509455 
    ## Run 369 stress 0.01882211 
    ## ... Procrustes: rmse 0.003098048  max resid 0.006367297 
    ## ... Similar to previous best
    ## Run 370 stress 0.01879604 
    ## ... Procrustes: rmse 0.000832191  max resid 0.001701559 
    ## ... Similar to previous best
    ## Run 371 stress 0.01879553 
    ## ... Procrustes: rmse 0.0006253367  max resid 0.001275174 
    ## ... Similar to previous best
    ## Run 372 stress 0.01879563 
    ## ... Procrustes: rmse 0.0006733357  max resid 0.001374232 
    ## ... Similar to previous best
    ## Run 373 stress 0.0187959 
    ## ... Procrustes: rmse 0.0007797082  max resid 0.001593399 
    ## ... Similar to previous best
    ## Run 374 stress 0.01881885 
    ## ... Procrustes: rmse 0.003852978  max resid 0.007919829 
    ## ... Similar to previous best
    ## Run 375 stress 0.02509443 
    ## Run 376 stress 0.02492194 
    ## Run 377 stress 0.0187959 
    ## ... Procrustes: rmse 0.000775691  max resid 0.001585085 
    ## ... Similar to previous best
    ## Run 378 stress 0.02509424 
    ## Run 379 stress 0.01879596 
    ## ... Procrustes: rmse 0.0008079678  max resid 0.001651889 
    ## ... Similar to previous best
    ## Run 380 stress 0.02486124 
    ## Run 381 stress 0.02492183 
    ## Run 382 stress 0.0187958 
    ## ... Procrustes: rmse 0.0007494944  max resid 0.001531325 
    ## ... Similar to previous best
    ## Run 383 stress 0.01879558 
    ## ... Procrustes: rmse 0.0006463915  max resid 0.001318647 
    ## ... Similar to previous best
    ## Run 384 stress 0.02520885 
    ## Run 385 stress 0.02509484 
    ## Run 386 stress 0.02492166 
    ## Run 387 stress 0.02650418 
    ## Run 388 stress 0.0250945 
    ## Run 389 stress 0.02509451 
    ## Run 390 stress 0.01879595 
    ## ... Procrustes: rmse 0.0008105117  max resid 0.00165705 
    ## ... Similar to previous best
    ## Run 391 stress 0.0250947 
    ## Run 392 stress 0.01879721 
    ## ... Procrustes: rmse 0.001204912  max resid 0.0024699 
    ## ... Similar to previous best
    ## Run 393 stress 0.02509448 
    ## Run 394 stress 0.02509468 
    ## Run 395 stress 0.0249221 
    ## Run 396 stress 0.01885442 
    ## ... Procrustes: rmse 0.004900421  max resid 0.009763128 
    ## ... Similar to previous best
    ## Run 397 stress 0.02509474 
    ## Run 398 stress 0.01879541 
    ## ... Procrustes: rmse 0.0005582846  max resid 0.001136779 
    ## ... Similar to previous best
    ## Run 399 stress 0.01883119 
    ## ... Procrustes: rmse 0.002652646  max resid 0.005048237 
    ## ... Similar to previous best
    ## Run 400 stress 0.01917924 
    ## ... Procrustes: rmse 0.01555444  max resid 0.03193873 
    ## Run 401 stress 0.01880282 
    ## ... Procrustes: rmse 0.002221473  max resid 0.004565262 
    ## ... Similar to previous best
    ## Run 402 stress 0.01879559 
    ## ... Procrustes: rmse 0.0006556677  max resid 0.001337754 
    ## ... Similar to previous best
    ## Run 403 stress 0.0204953 
    ## Run 404 stress 0.01901494 
    ## ... Procrustes: rmse 0.01124326  max resid 0.02290734 
    ## Run 405 stress 0.01879596 
    ## ... Procrustes: rmse 0.0007882957  max resid 0.00161099 
    ## ... Similar to previous best
    ## Run 406 stress 0.02509478 
    ## Run 407 stress 0.02486489 
    ## Run 408 stress 0.01941028 
    ## Run 409 stress 0.01883813 
    ## ... Procrustes: rmse 0.003512588  max resid 0.006823826 
    ## ... Similar to previous best
    ## Run 410 stress 0.01887021 
    ## ... Procrustes: rmse 0.005399653  max resid 0.01080921 
    ## Run 411 stress 0.01879588 
    ## ... Procrustes: rmse 0.0007853462  max resid 0.001605189 
    ## ... Similar to previous best
    ## Run 412 stress 0.01884011 
    ## ... Procrustes: rmse 0.003711131  max resid 0.007247994 
    ## ... Similar to previous best
    ## Run 413 stress 0.01879585 
    ## ... Procrustes: rmse 0.0007727155  max resid 0.001579164 
    ## ... Similar to previous best
    ## Run 414 stress 0.02509471 
    ## Run 415 stress 0.02492182 
    ## Run 416 stress 0.01919185 
    ## ... Procrustes: rmse 0.01581262  max resid 0.0324361 
    ## Run 417 stress 0.01879548 
    ## ... Procrustes: rmse 0.0004981442  max resid 0.001034223 
    ## ... Similar to previous best
    ## Run 418 stress 0.02509472 
    ## Run 419 stress 0.01879548 
    ## ... Procrustes: rmse 0.000593779  max resid 0.001210031 
    ## ... Similar to previous best
    ## Run 420 stress 0.0187959 
    ## ... Procrustes: rmse 0.0007914177  max resid 0.001617731 
    ## ... Similar to previous best
    ## Run 421 stress 0.01879564 
    ## ... Procrustes: rmse 0.0006766665  max resid 0.001381033 
    ## ... Similar to previous best
    ## Run 422 stress 0.02509456 
    ## Run 423 stress 0.02492195 
    ## Run 424 stress 0.01895344 
    ## ... Procrustes: rmse 0.009857089  max resid 0.02022137 
    ## Run 425 stress 0.01879592 
    ## ... Procrustes: rmse 0.0007822797  max resid 0.001598936 
    ## ... Similar to previous best
    ## Run 426 stress 0.01879589 
    ## ... Procrustes: rmse 0.0007857838  max resid 0.001606081 
    ## ... Similar to previous best
    ## Run 427 stress 0.02509475 
    ## Run 428 stress 0.3513261 
    ## Run 429 stress 0.02492158 
    ## Run 430 stress 0.02486138 
    ## Run 431 stress 0.0249268 
    ## Run 432 stress 0.02492183 
    ## Run 433 stress 0.02492182 
    ## Run 434 stress 0.01879735 
    ## ... Procrustes: rmse 0.001240105  max resid 0.002542533 
    ## ... Similar to previous best
    ## Run 435 stress 0.0250948 
    ## Run 436 stress 0.02492194 
    ## Run 437 stress 0.01889387 
    ## ... Procrustes: rmse 0.007840857  max resid 0.01609616 
    ## Run 438 stress 0.0187961 
    ## ... Procrustes: rmse 0.0008568781  max resid 0.001752451 
    ## ... Similar to previous best
    ## Run 439 stress 0.02492175 
    ## Run 440 stress 0.02509482 
    ## Run 441 stress 0.01879779 
    ## ... Procrustes: rmse 0.001348821  max resid 0.002766308 
    ## ... Similar to previous best
    ## Run 442 stress 0.01888867 
    ## ... Procrustes: rmse 0.007631881  max resid 0.01566764 
    ## Run 443 stress 0.02509443 
    ## Run 444 stress 0.01879562 
    ## ... Procrustes: rmse 0.0006694562  max resid 0.001366225 
    ## ... Similar to previous best
    ## Run 445 stress 0.02509462 
    ## Run 446 stress 0.02492197 
    ## Run 447 stress 0.01883705 
    ## ... Procrustes: rmse 0.003393639  max resid 0.006568722 
    ## ... Similar to previous best
    ## Run 448 stress 0.01879588 
    ## ... Procrustes: rmse 0.0007717087  max resid 0.001577085 
    ## ... Similar to previous best
    ## Run 449 stress 0.01879574 
    ## ... Procrustes: rmse 0.0007199219  max resid 0.001470331 
    ## ... Similar to previous best
    ## Run 450 stress 0.02509479 
    ## Run 451 stress 0.01883849 
    ## ... Procrustes: rmse 0.005197366  max resid 0.01067921 
    ## Run 452 stress 0.02497163 
    ## Run 453 stress 0.0250947 
    ## Run 454 stress 0.02509452 
    ## Run 455 stress 0.02492159 
    ## Run 456 stress 0.02509485 
    ## Run 457 stress 0.01880851 
    ## ... Procrustes: rmse 0.002894662  max resid 0.00595145 
    ## ... Similar to previous best
    ## Run 458 stress 0.02509437 
    ## Run 459 stress 0.02509467 
    ## Run 460 stress 0.02509479 
    ## Run 461 stress 0.02486111 
    ## Run 462 stress 0.02509453 
    ## Run 463 stress 0.01879574 
    ## ... Procrustes: rmse 0.000725272  max resid 0.001481314 
    ## ... Similar to previous best
    ## Run 464 stress 0.02509483 
    ## Run 465 stress 0.01879544 
    ## ... Procrustes: rmse 0.0005743591  max resid 0.001169728 
    ## ... Similar to previous best
    ## Run 466 stress 0.01890559 
    ## ... Procrustes: rmse 0.008293498  max resid 0.01702016 
    ## Run 467 stress 0.01879746 
    ## ... Procrustes: rmse 0.001269241  max resid 0.002603046 
    ## ... Similar to previous best
    ## Run 468 stress 0.01888479 
    ## ... Procrustes: rmse 0.006692184  max resid 0.01350591 
    ## Run 469 stress 0.02520895 
    ## Run 470 stress 0.02509462 
    ## Run 471 stress 0.02486115 
    ## Run 472 stress 0.01886685 
    ## ... Procrustes: rmse 0.005701838  max resid 0.01144212 
    ## Run 473 stress 0.01879578 
    ## ... Procrustes: rmse 0.0007401311  max resid 0.001511888 
    ## ... Similar to previous best
    ## Run 474 stress 0.01879965 
    ## ... Procrustes: rmse 0.001700403  max resid 0.003503388 
    ## ... Similar to previous best
    ## Run 475 stress 0.02492177 
    ## Run 476 stress 0.02509476 
    ## Run 477 stress 0.01890431 
    ## ... Procrustes: rmse 0.007496934  max resid 0.01518165 
    ## Run 478 stress 0.02492181 
    ## Run 479 stress 0.01879556 
    ## ... Procrustes: rmse 0.0006265122  max resid 0.001277633 
    ## ... Similar to previous best
    ## Run 480 stress 0.02509435 
    ## Run 481 stress 0.01879613 
    ## ... Procrustes: rmse 0.0008803146  max resid 0.001800926 
    ## ... Similar to previous best
    ## Run 482 stress 0.02509489 
    ## Run 483 stress 0.01879595 
    ## ... Procrustes: rmse 0.0007121931  max resid 0.001474857 
    ## ... Similar to previous best
    ## Run 484 stress 0.0250946 
    ## Run 485 stress 0.01879592 
    ## ... Procrustes: rmse 0.0008002188  max resid 0.001635861 
    ## ... Similar to previous best
    ## Run 486 stress 0.0250943 
    ## Run 487 stress 0.02509451 
    ## Run 488 stress 0.01880759 
    ## ... Procrustes: rmse 0.002805578  max resid 0.005767601 
    ## ... Similar to previous best
    ## Run 489 stress 0.02486128 
    ## Run 490 stress 0.02492212 
    ## Run 491 stress 0.0188319 
    ## ... Procrustes: rmse 0.002748428  max resid 0.005186162 
    ## ... Similar to previous best
    ## Run 492 stress 0.01879772 
    ## ... Procrustes: rmse 0.001316904  max resid 0.002701186 
    ## ... Similar to previous best
    ## Run 493 stress 0.01879593 
    ## ... Procrustes: rmse 0.0007021486  max resid 0.0014542 
    ## ... Similar to previous best
    ## Run 494 stress 0.01879571 
    ## ... Procrustes: rmse 0.0007121539  max resid 0.001454282 
    ## ... Similar to previous best
    ## Run 495 stress 0.01879577 
    ## ... Procrustes: rmse 0.0007377906  max resid 0.001507148 
    ## ... Similar to previous best
    ## Run 496 stress 0.02486118 
    ## Run 497 stress 0.3773891 
    ## Run 498 stress 0.02509443 
    ## Run 499 stress 0.02509457 
    ## Run 500 stress 0.01879602 
    ## ... Procrustes: rmse 0.000825148  max resid 0.001687058 
    ## ... Similar to previous best
    ## *** Best solution repeated 166 times

``` r
# Mixed and stratified lakes
SD_beta_MS_NMDS <- metaMDS(SD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09977534 
    ## Run 2 stress 0.09308948 
    ## Run 3 stress 0.09030405 
    ## Run 4 stress 0.09286084 
    ## Run 5 stress 0.09465907 
    ## Run 6 stress 0.09400517 
    ## Run 7 stress 0.09721206 
    ## Run 8 stress 0.08503462 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001220524  max resid 0.003402932 
    ## ... Similar to previous best
    ## Run 9 stress 0.08973862 
    ## Run 10 stress 0.08503455 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000150532  max resid 0.0004011701 
    ## ... Similar to previous best
    ## Run 11 stress 0.09337247 
    ## Run 12 stress 0.09030408 
    ## Run 13 stress 0.09286083 
    ## Run 14 stress 0.09535508 
    ## Run 15 stress 0.08503489 
    ## ... Procrustes: rmse 0.0005246271  max resid 0.001098889 
    ## ... Similar to previous best
    ## Run 16 stress 0.1054548 
    ## Run 17 stress 0.09268326 
    ## Run 18 stress 0.08973879 
    ## Run 19 stress 0.09465917 
    ## Run 20 stress 0.09380577 
    ## Run 21 stress 0.08773479 
    ## Run 22 stress 0.0933729 
    ## Run 23 stress 0.09590504 
    ## Run 24 stress 0.09145332 
    ## Run 25 stress 0.09286104 
    ## Run 26 stress 0.3126992 
    ## Run 27 stress 0.09168928 
    ## Run 28 stress 0.08973881 
    ## Run 29 stress 0.08973876 
    ## Run 30 stress 0.09030395 
    ## Run 31 stress 0.08973879 
    ## Run 32 stress 0.09721202 
    ## Run 33 stress 0.09407432 
    ## Run 34 stress 0.08503492 
    ## ... Procrustes: rmse 0.000553735  max resid 0.001160636 
    ## ... Similar to previous best
    ## Run 35 stress 0.09030396 
    ## Run 36 stress 0.09030404 
    ## Run 37 stress 0.09969977 
    ## Run 38 stress 0.09407985 
    ## Run 39 stress 0.09465904 
    ## Run 40 stress 0.09337267 
    ## Run 41 stress 0.09145324 
    ## Run 42 stress 0.09464437 
    ## Run 43 stress 0.09145334 
    ## Run 44 stress 0.09407972 
    ## Run 45 stress 0.09969975 
    ## Run 46 stress 0.09374237 
    ## Run 47 stress 0.09381846 
    ## Run 48 stress 0.09159088 
    ## Run 49 stress 0.09159086 
    ## Run 50 stress 0.09159102 
    ## Run 51 stress 0.09268329 
    ## Run 52 stress 0.09407965 
    ## Run 53 stress 0.1038047 
    ## Run 54 stress 0.09465905 
    ## Run 55 stress 0.09268333 
    ## Run 56 stress 0.08503566 
    ## ... Procrustes: rmse 0.001020316  max resid 0.00262216 
    ## ... Similar to previous best
    ## Run 57 stress 0.09407991 
    ## Run 58 stress 0.1038045 
    ## Run 59 stress 0.09168938 
    ## Run 60 stress 0.08503683 
    ## ... Procrustes: rmse 0.001574822  max resid 0.004018256 
    ## ... Similar to previous best
    ## Run 61 stress 0.08503509 
    ## ... Procrustes: rmse 0.0006508559  max resid 0.001335771 
    ## ... Similar to previous best
    ## Run 62 stress 0.09465917 
    ## Run 63 stress 0.08503501 
    ## ... Procrustes: rmse 0.0005978248  max resid 0.001254528 
    ## ... Similar to previous best
    ## Run 64 stress 0.09337266 
    ## Run 65 stress 0.09268354 
    ## Run 66 stress 0.09407112 
    ## Run 67 stress 0.08503591 
    ## ... Procrustes: rmse 0.001156364  max resid 0.002979781 
    ## ... Similar to previous best
    ## Run 68 stress 0.08973865 
    ## Run 69 stress 0.0940799 
    ## Run 70 stress 0.09159084 
    ## Run 71 stress 0.0844027 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01417179  max resid 0.04243614 
    ## Run 72 stress 0.08440263 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004674801  max resid 0.0009251534 
    ## ... Similar to previous best
    ## Run 73 stress 0.09308957 
    ## Run 74 stress 0.09407967 
    ## Run 75 stress 0.08503469 
    ## Run 76 stress 0.09760835 
    ## Run 77 stress 0.09408005 
    ## Run 78 stress 0.09308972 
    ## Run 79 stress 0.09539198 
    ## Run 80 stress 0.09721196 
    ## Run 81 stress 0.09400519 
    ## Run 82 stress 0.08773485 
    ## Run 83 stress 0.2574234 
    ## Run 84 stress 0.09159095 
    ## Run 85 stress 0.09030401 
    ## Run 86 stress 0.09286085 
    ## Run 87 stress 0.08503482 
    ## Run 88 stress 0.08973881 
    ## Run 89 stress 0.09445511 
    ## Run 90 stress 0.08440261 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003553523  max resid 0.0006803682 
    ## ... Similar to previous best
    ## Run 91 stress 0.08503459 
    ## Run 92 stress 0.09407988 
    ## Run 93 stress 0.09374287 
    ## Run 94 stress 0.09374352 
    ## Run 95 stress 0.09407383 
    ## Run 96 stress 0.09969981 
    ## Run 97 stress 0.1052851 
    ## Run 98 stress 0.09030394 
    ## Run 99 stress 0.09337249 
    ## Run 100 stress 0.0877347 
    ## Run 101 stress 0.09308954 
    ## Run 102 stress 0.1084716 
    ## Run 103 stress 0.09030407 
    ## Run 104 stress 0.3084172 
    ## Run 105 stress 0.09030396 
    ## Run 106 stress 0.08503536 
    ## Run 107 stress 0.09407966 
    ## Run 108 stress 0.09407971 
    ## Run 109 stress 0.0850347 
    ## Run 110 stress 0.08773466 
    ## Run 111 stress 0.08773471 
    ## Run 112 stress 0.09416144 
    ## Run 113 stress 0.09445496 
    ## Run 114 stress 0.09337282 
    ## Run 115 stress 0.09337242 
    ## Run 116 stress 0.09168957 
    ## Run 117 stress 0.08973863 
    ## Run 118 stress 0.08773481 
    ## Run 119 stress 0.08440261 
    ## ... Procrustes: rmse 0.000337307  max resid 0.0006608736 
    ## ... Similar to previous best
    ## Run 120 stress 0.09539192 
    ## Run 121 stress 0.0877347 
    ## Run 122 stress 0.09465906 
    ## Run 123 stress 0.09464463 
    ## Run 124 stress 0.1026215 
    ## Run 125 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001762076  max resid 0.0003495732 
    ## ... Similar to previous best
    ## Run 126 stress 0.08440256 
    ## ... Procrustes: rmse 0.000105317  max resid 0.0001906879 
    ## ... Similar to previous best
    ## Run 127 stress 0.09760832 
    ## Run 128 stress 0.1052853 
    ## Run 129 stress 0.1053425 
    ## Run 130 stress 0.09030402 
    ## Run 131 stress 0.09159086 
    ## Run 132 stress 0.09407986 
    ## Run 133 stress 0.09407106 
    ## Run 134 stress 0.09145324 
    ## Run 135 stress 0.09286083 
    ## Run 136 stress 0.09416133 
    ## Run 137 stress 0.09374266 
    ## Run 138 stress 0.0877347 
    ## Run 139 stress 0.09308952 
    ## Run 140 stress 0.08773485 
    ## Run 141 stress 0.0844026 
    ## ... Procrustes: rmse 0.0001558249  max resid 0.0002875161 
    ## ... Similar to previous best
    ## Run 142 stress 0.09286084 
    ## Run 143 stress 0.09159085 
    ## Run 144 stress 0.09465905 
    ## Run 145 stress 0.09030394 
    ## Run 146 stress 0.09337235 
    ## Run 147 stress 0.08973885 
    ## Run 148 stress 0.09407967 
    ## Run 149 stress 0.09416135 
    ## Run 150 stress 0.09030413 
    ## Run 151 stress 0.08773473 
    ## Run 152 stress 0.09286084 
    ## Run 153 stress 0.08773473 
    ## Run 154 stress 0.08773482 
    ## Run 155 stress 0.09168926 
    ## Run 156 stress 0.09030402 
    ## Run 157 stress 0.1038054 
    ## Run 158 stress 0.0916893 
    ## Run 159 stress 0.09535417 
    ## Run 160 stress 0.09407969 
    ## Run 161 stress 0.09760835 
    ## Run 162 stress 0.08973868 
    ## Run 163 stress 0.09337277 
    ## Run 164 stress 0.09168942 
    ## Run 165 stress 0.09286083 
    ## Run 166 stress 0.08773487 
    ## Run 167 stress 0.08503514 
    ## Run 168 stress 0.0850352 
    ## Run 169 stress 0.09030393 
    ## Run 170 stress 0.1038046 
    ## Run 171 stress 0.09969994 
    ## Run 172 stress 0.09030409 
    ## Run 173 stress 0.09669864 
    ## Run 174 stress 0.0940437 
    ## Run 175 stress 0.09145322 
    ## Run 176 stress 0.08973869 
    ## Run 177 stress 0.08503507 
    ## Run 178 stress 0.2758279 
    ## Run 179 stress 0.08973881 
    ## Run 180 stress 0.09465908 
    ## Run 181 stress 0.08440254 
    ## ... Procrustes: rmse 8.068998e-05  max resid 0.0001547508 
    ## ... Similar to previous best
    ## Run 182 stress 0.09030402 
    ## Run 183 stress 0.103804 
    ## Run 184 stress 0.08773475 
    ## Run 185 stress 0.1038039 
    ## Run 186 stress 0.08503461 
    ## Run 187 stress 0.09268336 
    ## Run 188 stress 0.09145326 
    ## Run 189 stress 0.09403433 
    ## Run 190 stress 0.09721199 
    ## Run 191 stress 0.08440262 
    ## ... Procrustes: rmse 0.0001719711  max resid 0.0003328634 
    ## ... Similar to previous best
    ## Run 192 stress 0.09268396 
    ## Run 193 stress 0.08503615 
    ## Run 194 stress 0.08503459 
    ## Run 195 stress 0.08773468 
    ## Run 196 stress 0.09445483 
    ## Run 197 stress 0.1026216 
    ## Run 198 stress 0.08503485 
    ## Run 199 stress 0.09337298 
    ## Run 200 stress 0.09268319 
    ## Run 201 stress 0.09407986 
    ## Run 202 stress 0.0844026 
    ## ... Procrustes: rmse 0.0001535355  max resid 0.0002896026 
    ## ... Similar to previous best
    ## Run 203 stress 0.08440258 
    ## ... Procrustes: rmse 0.0001229527  max resid 0.0002550961 
    ## ... Similar to previous best
    ## Run 204 stress 0.09159088 
    ## Run 205 stress 0.3309708 
    ## Run 206 stress 0.0940799 
    ## Run 207 stress 0.09286084 
    ## Run 208 stress 0.09465907 
    ## Run 209 stress 0.09407472 
    ## Run 210 stress 0.09159088 
    ## Run 211 stress 0.09159083 
    ## Run 212 stress 0.08503461 
    ## Run 213 stress 0.09145329 
    ## Run 214 stress 0.09464374 
    ## Run 215 stress 0.09407127 
    ## Run 216 stress 0.0933725 
    ## Run 217 stress 0.08773471 
    ## Run 218 stress 0.0940798 
    ## Run 219 stress 0.09030407 
    ## Run 220 stress 0.09268388 
    ## Run 221 stress 0.1038046 
    ## Run 222 stress 0.09374276 
    ## Run 223 stress 0.08440254 
    ## ... Procrustes: rmse 8.069051e-05  max resid 0.0001505872 
    ## ... Similar to previous best
    ## Run 224 stress 0.09712929 
    ## Run 225 stress 0.09908192 
    ## Run 226 stress 0.09374224 
    ## Run 227 stress 0.09286085 
    ## Run 228 stress 0.08503504 
    ## Run 229 stress 0.09445475 
    ## Run 230 stress 0.09268358 
    ## Run 231 stress 0.08503464 
    ## Run 232 stress 0.08973877 
    ## Run 233 stress 0.2707116 
    ## Run 234 stress 0.1053428 
    ## Run 235 stress 0.09416141 
    ## Run 236 stress 0.09030402 
    ## Run 237 stress 0.08503455 
    ## Run 238 stress 0.09408017 
    ## Run 239 stress 0.09590922 
    ## Run 240 stress 0.08503479 
    ## Run 241 stress 0.09337236 
    ## Run 242 stress 0.08503505 
    ## Run 243 stress 0.09447345 
    ## Run 244 stress 0.09416132 
    ## Run 245 stress 0.08503503 
    ## Run 246 stress 0.08503473 
    ## Run 247 stress 0.09400528 
    ## Run 248 stress 0.09268348 
    ## Run 249 stress 0.08973864 
    ## Run 250 stress 0.08973862 
    ## Run 251 stress 0.09030399 
    ## Run 252 stress 0.09969978 
    ## Run 253 stress 0.1026218 
    ## Run 254 stress 0.09535455 
    ## Run 255 stress 0.09407987 
    ## Run 256 stress 0.09337271 
    ## Run 257 stress 0.1038089 
    ## Run 258 stress 0.09407105 
    ## Run 259 stress 0.0877348 
    ## Run 260 stress 0.09286084 
    ## Run 261 stress 0.09407968 
    ## Run 262 stress 0.09539201 
    ## Run 263 stress 0.08973875 
    ## Run 264 stress 0.09760829 
    ## Run 265 stress 0.09337237 
    ## Run 266 stress 0.09337238 
    ## Run 267 stress 0.09374227 
    ## Run 268 stress 0.09407123 
    ## Run 269 stress 0.08440263 
    ## ... Procrustes: rmse 0.0002347002  max resid 0.000467476 
    ## ... Similar to previous best
    ## Run 270 stress 0.09969965 
    ## Run 271 stress 0.08773483 
    ## Run 272 stress 0.09445463 
    ## Run 273 stress 0.08440254 
    ## ... Procrustes: rmse 7.072788e-05  max resid 0.0001338563 
    ## ... Similar to previous best
    ## Run 274 stress 0.0938058 
    ## Run 275 stress 0.09445487 
    ## Run 276 stress 0.09286084 
    ## Run 277 stress 0.08440255 
    ## ... Procrustes: rmse 8.849063e-05  max resid 0.0001859701 
    ## ... Similar to previous best
    ## Run 278 stress 0.0940798 
    ## Run 279 stress 0.09760819 
    ## Run 280 stress 0.09407148 
    ## Run 281 stress 0.09268396 
    ## Run 282 stress 0.09970001 
    ## Run 283 stress 0.09159084 
    ## Run 284 stress 0.09030393 
    ## Run 285 stress 0.08503542 
    ## Run 286 stress 0.103804 
    ## Run 287 stress 0.08503496 
    ## Run 288 stress 0.09168953 
    ## Run 289 stress 0.08503456 
    ## Run 290 stress 0.09407124 
    ## Run 291 stress 0.09145324 
    ## Run 292 stress 0.08773466 
    ## Run 293 stress 0.09030411 
    ## Run 294 stress 0.09610799 
    ## Run 295 stress 0.09030413 
    ## Run 296 stress 0.08440254 
    ## ... Procrustes: rmse 7.510591e-05  max resid 0.0001605382 
    ## ... Similar to previous best
    ## Run 297 stress 0.09380586 
    ## Run 298 stress 0.09030406 
    ## Run 299 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001984542  max resid 0.0003881941 
    ## ... Similar to previous best
    ## Run 300 stress 0.09337277 
    ## Run 301 stress 0.09381846 
    ## Run 302 stress 0.08503594 
    ## Run 303 stress 0.105285 
    ## Run 304 stress 0.09030394 
    ## Run 305 stress 0.09168963 
    ## Run 306 stress 0.09337265 
    ## Run 307 stress 0.0996997 
    ## Run 308 stress 0.2835514 
    ## Run 309 stress 0.08503507 
    ## Run 310 stress 0.08503541 
    ## Run 311 stress 0.0938058 
    ## Run 312 stress 0.08503455 
    ## Run 313 stress 0.09407973 
    ## Run 314 stress 0.1038053 
    ## Run 315 stress 0.09969982 
    ## Run 316 stress 0.09030395 
    ## Run 317 stress 0.1038049 
    ## Run 318 stress 0.0933729 
    ## Run 319 stress 0.08973863 
    ## Run 320 stress 0.08503502 
    ## Run 321 stress 0.09590792 
    ## Run 322 stress 0.09465908 
    ## Run 323 stress 0.08973892 
    ## Run 324 stress 0.09030404 
    ## Run 325 stress 0.09337237 
    ## Run 326 stress 0.09268343 
    ## Run 327 stress 0.1026215 
    ## Run 328 stress 0.09159092 
    ## Run 329 stress 0.09380577 
    ## Run 330 stress 0.09374263 
    ## Run 331 stress 0.09465907 
    ## Run 332 stress 0.09159086 
    ## Run 333 stress 0.09268323 
    ## Run 334 stress 0.09407965 
    ## Run 335 stress 0.08503638 
    ## Run 336 stress 0.09407335 
    ## Run 337 stress 0.09030397 
    ## Run 338 stress 0.09969972 
    ## Run 339 stress 0.09447333 
    ## Run 340 stress 0.09030394 
    ## Run 341 stress 0.09374255 
    ## Run 342 stress 0.08973883 
    ## Run 343 stress 0.09168947 
    ## Run 344 stress 0.095909 
    ## Run 345 stress 0.09268401 
    ## Run 346 stress 0.09407982 
    ## Run 347 stress 0.09465905 
    ## Run 348 stress 0.08503717 
    ## Run 349 stress 0.1038047 
    ## Run 350 stress 0.09145338 
    ## Run 351 stress 0.0996999 
    ## Run 352 stress 0.08503456 
    ## Run 353 stress 0.1006357 
    ## Run 354 stress 0.08440264 
    ## ... Procrustes: rmse 0.0001855055  max resid 0.000340724 
    ## ... Similar to previous best
    ## Run 355 stress 0.09969995 
    ## Run 356 stress 0.3347414 
    ## Run 357 stress 0.09400518 
    ## Run 358 stress 0.09268385 
    ## Run 359 stress 0.09159084 
    ## Run 360 stress 0.08773474 
    ## Run 361 stress 0.0916894 
    ## Run 362 stress 0.08773479 
    ## Run 363 stress 0.08773471 
    ## Run 364 stress 0.09159091 
    ## Run 365 stress 0.0933723 
    ## Run 366 stress 0.09417778 
    ## Run 367 stress 0.09145337 
    ## Run 368 stress 0.09308997 
    ## Run 369 stress 0.09168942 
    ## Run 370 stress 0.09030404 
    ## Run 371 stress 0.1038042 
    ## Run 372 stress 0.08503661 
    ## Run 373 stress 0.08773471 
    ## Run 374 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001006205  max resid 0.0001732253 
    ## ... Similar to previous best
    ## Run 375 stress 0.09145323 
    ## Run 376 stress 0.090304 
    ## Run 377 stress 0.09159083 
    ## Run 378 stress 0.08773468 
    ## Run 379 stress 0.2739367 
    ## Run 380 stress 0.09308951 
    ## Run 381 stress 0.09268349 
    ## Run 382 stress 0.08973884 
    ## Run 383 stress 0.09337239 
    ## Run 384 stress 0.0850349 
    ## Run 385 stress 0.0940749 
    ## Run 386 stress 0.09374345 
    ## Run 387 stress 0.09407968 
    ## Run 388 stress 0.08503518 
    ## Run 389 stress 0.09308973 
    ## Run 390 stress 0.08440257 
    ## ... Procrustes: rmse 0.0001171347  max resid 0.0002460099 
    ## ... Similar to previous best
    ## Run 391 stress 0.08503516 
    ## Run 392 stress 0.09030399 
    ## Run 393 stress 0.09030397 
    ## Run 394 stress 0.09286088 
    ## Run 395 stress 0.09337247 
    ## Run 396 stress 0.1072866 
    ## Run 397 stress 0.09030394 
    ## Run 398 stress 0.0996997 
    ## Run 399 stress 0.090304 
    ## Run 400 stress 0.09286083 
    ## Run 401 stress 0.08440269 
    ## ... Procrustes: rmse 0.0002264631  max resid 0.0004105863 
    ## ... Similar to previous best
    ## Run 402 stress 0.1053423 
    ## Run 403 stress 0.0915909 
    ## Run 404 stress 0.09308951 
    ## Run 405 stress 0.09145323 
    ## Run 406 stress 0.09337275 
    ## Run 407 stress 0.08503499 
    ## Run 408 stress 0.1103048 
    ## Run 409 stress 0.09030413 
    ## Run 410 stress 0.0850349 
    ## Run 411 stress 0.09286085 
    ## Run 412 stress 0.09308976 
    ## Run 413 stress 0.08503543 
    ## Run 414 stress 0.09286083 
    ## Run 415 stress 0.09337266 
    ## Run 416 stress 0.09337269 
    ## Run 417 stress 0.1052849 
    ## Run 418 stress 0.08503491 
    ## Run 419 stress 0.1026216 
    ## Run 420 stress 0.09721201 
    ## Run 421 stress 0.08440258 
    ## ... Procrustes: rmse 9.206574e-05  max resid 0.0001744903 
    ## ... Similar to previous best
    ## Run 422 stress 0.09407972 
    ## Run 423 stress 0.09464479 
    ## Run 424 stress 0.09168941 
    ## Run 425 stress 0.09030404 
    ## Run 426 stress 0.09308963 
    ## Run 427 stress 0.08773475 
    ## Run 428 stress 0.09416144 
    ## Run 429 stress 0.08503501 
    ## Run 430 stress 0.09168945 
    ## Run 431 stress 0.09159083 
    ## Run 432 stress 0.09159084 
    ## Run 433 stress 0.09374295 
    ## Run 434 stress 0.09760835 
    ## Run 435 stress 0.09535449 
    ## Run 436 stress 0.09159089 
    ## Run 437 stress 0.09030413 
    ## Run 438 stress 0.09407998 
    ## Run 439 stress 0.09374179 
    ## Run 440 stress 0.09407986 
    ## Run 441 stress 0.09168941 
    ## Run 442 stress 0.0877347 
    ## Run 443 stress 0.08973876 
    ## Run 444 stress 0.08503474 
    ## Run 445 stress 0.09407136 
    ## Run 446 stress 0.09416137 
    ## Run 447 stress 0.09308976 
    ## Run 448 stress 0.08440274 
    ## ... Procrustes: rmse 0.0002250472  max resid 0.00041022 
    ## ... Similar to previous best
    ## Run 449 stress 0.08973866 
    ## Run 450 stress 0.09145331 
    ## Run 451 stress 0.1038046 
    ## Run 452 stress 0.09407399 
    ## Run 453 stress 0.09308947 
    ## Run 454 stress 0.09535511 
    ## Run 455 stress 0.09445498 
    ## Run 456 stress 0.09286085 
    ## Run 457 stress 0.09286084 
    ## Run 458 stress 0.08440275 
    ## ... Procrustes: rmse 0.0002664119  max resid 0.0004829787 
    ## ... Similar to previous best
    ## Run 459 stress 0.09159084 
    ## Run 460 stress 0.09407965 
    ## Run 461 stress 0.09407985 
    ## Run 462 stress 0.08440262 
    ## ... Procrustes: rmse 0.0001743251  max resid 0.0003225216 
    ## ... Similar to previous best
    ## Run 463 stress 0.08773474 
    ## Run 464 stress 0.105285 
    ## Run 465 stress 0.09286087 
    ## Run 466 stress 0.1038044 
    ## Run 467 stress 0.09380578 
    ## Run 468 stress 0.08773483 
    ## Run 469 stress 0.09159088 
    ## Run 470 stress 0.09416141 
    ## Run 471 stress 0.09308955 
    ## Run 472 stress 0.09168926 
    ## Run 473 stress 0.09168938 
    ## Run 474 stress 0.08440252 
    ## ... Procrustes: rmse 9.695134e-06  max resid 2.16309e-05 
    ## ... Similar to previous best
    ## Run 475 stress 0.09337244 
    ## Run 476 stress 0.08503465 
    ## Run 477 stress 0.09590666 
    ## Run 478 stress 0.08503491 
    ## Run 479 stress 0.09268339 
    ## Run 480 stress 0.09268329 
    ## Run 481 stress 0.08503674 
    ## Run 482 stress 0.08773487 
    ## Run 483 stress 0.09407988 
    ## Run 484 stress 0.09268376 
    ## Run 485 stress 0.3056852 
    ## Run 486 stress 0.09321725 
    ## Run 487 stress 0.09030415 
    ## Run 488 stress 0.09404372 
    ## Run 489 stress 0.09286092 
    ## Run 490 stress 0.09159102 
    ## Run 491 stress 0.09286091 
    ## Run 492 stress 0.1038047 
    ## Run 493 stress 0.09969988 
    ## Run 494 stress 0.0897387 
    ## Run 495 stress 0.09030398 
    ## Run 496 stress 0.09969989 
    ## Run 497 stress 0.08503494 
    ## Run 498 stress 0.08973864 
    ## Run 499 stress 0.08973879 
    ## Run 500 stress 0.09145322 
    ## *** Best solution repeated 22 times

``` r
# Ocean sites and mixed lakes
SD_beta_OM_NMDS <- metaMDS(SD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07365783 
    ## Run 1 stress 0.08000468 
    ## Run 2 stress 0.07629237 
    ## Run 3 stress 0.07365785 
    ## ... Procrustes: rmse 6.722695e-05  max resid 0.0001566428 
    ## ... Similar to previous best
    ## Run 4 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743037  max resid 0.05098635 
    ## Run 5 stress 0.07629234 
    ## Run 6 stress 0.07629245 
    ## Run 7 stress 0.07629237 
    ## Run 8 stress 0.08000471 
    ## Run 9 stress 0.07732936 
    ## Run 10 stress 0.08000468 
    ## Run 11 stress 0.07629241 
    ## Run 12 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174472  max resid 0.0510635 
    ## Run 13 stress 0.08288042 
    ## Run 14 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743761  max resid 0.05101547 
    ## Run 15 stress 0.08000471 
    ## Run 16 stress 0.2950735 
    ## Run 17 stress 0.07732925 
    ## Run 18 stress 0.07365783 
    ## ... Procrustes: rmse 1.029755e-05  max resid 2.26477e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745124  max resid 0.05107305 
    ## Run 20 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744758  max resid 0.05105709 
    ## Run 21 stress 0.07629235 
    ## Run 22 stress 0.08288054 
    ## Run 23 stress 0.07629234 
    ## Run 24 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001434063  max resid 0.0003420031 
    ## ... Similar to previous best
    ## Run 25 stress 0.07629242 
    ## Run 26 stress 0.07629239 
    ## Run 27 stress 0.07378231 
    ## ... Procrustes: rmse 0.01743458  max resid 0.05103473 
    ## Run 28 stress 0.07365785 
    ## ... Procrustes: rmse 6.259445e-05  max resid 0.0001436909 
    ## ... Similar to previous best
    ## Run 29 stress 0.07629237 
    ## Run 30 stress 0.08000467 
    ## Run 31 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001235401  max resid 0.0002922211 
    ## ... Similar to previous best
    ## Run 32 stress 0.07365784 
    ## ... Procrustes: rmse 4.057575e-05  max resid 9.66128e-05 
    ## ... Similar to previous best
    ## Run 33 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001440943  max resid 0.0003331943 
    ## ... Similar to previous best
    ## Run 34 stress 0.08000471 
    ## Run 35 stress 0.08288034 
    ## Run 36 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001109052  max resid 0.0002521377 
    ## ... Similar to previous best
    ## Run 37 stress 0.07629242 
    ## Run 38 stress 0.08000472 
    ## Run 39 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744414  max resid 0.05105257 
    ## Run 40 stress 0.07365784 
    ## ... Procrustes: rmse 4.94525e-05  max resid 0.0001161669 
    ## ... Similar to previous best
    ## Run 41 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001296227  max resid 0.0003095151 
    ## ... Similar to previous best
    ## Run 42 stress 0.07629234 
    ## Run 43 stress 0.0800047 
    ## Run 44 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744524  max resid 0.05104936 
    ## Run 45 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744532  max resid 0.05104468 
    ## Run 46 stress 0.07732928 
    ## Run 47 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001226371  max resid 0.0002883977 
    ## ... Similar to previous best
    ## Run 48 stress 0.07629238 
    ## Run 49 stress 0.07629243 
    ## Run 50 stress 0.07365784 
    ## ... Procrustes: rmse 3.190599e-05  max resid 7.502414e-05 
    ## ... Similar to previous best
    ## Run 51 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744675  max resid 0.05104319 
    ## Run 52 stress 0.07629233 
    ## Run 53 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744121  max resid 0.051025 
    ## Run 54 stress 0.08000472 
    ## Run 55 stress 0.07629244 
    ## Run 56 stress 0.07629236 
    ## Run 57 stress 0.07629234 
    ## Run 58 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744675  max resid 0.05104539 
    ## Run 59 stress 0.07629244 
    ## Run 60 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744217  max resid 0.05102333 
    ## Run 61 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744652  max resid 0.05104066 
    ## Run 62 stress 0.08288037 
    ## Run 63 stress 0.07629247 
    ## Run 64 stress 0.07629242 
    ## Run 65 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001081712  max resid 0.0002577757 
    ## ... Similar to previous best
    ## Run 66 stress 0.07629234 
    ## Run 67 stress 0.07629243 
    ## Run 68 stress 0.08000478 
    ## Run 69 stress 0.07365784 
    ## ... Procrustes: rmse 3.972001e-05  max resid 9.479151e-05 
    ## ... Similar to previous best
    ## Run 70 stress 0.07629237 
    ## Run 71 stress 0.08233757 
    ## Run 72 stress 0.07629243 
    ## Run 73 stress 0.07378227 
    ## ... Procrustes: rmse 0.01742697  max resid 0.05096331 
    ## Run 74 stress 0.07629235 
    ## Run 75 stress 0.07629233 
    ## Run 76 stress 0.07629233 
    ## Run 77 stress 0.0762924 
    ## Run 78 stress 0.08000467 
    ## Run 79 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001513494  max resid 0.000358114 
    ## ... Similar to previous best
    ## Run 80 stress 0.07629236 
    ## Run 81 stress 0.07378229 
    ## ... Procrustes: rmse 0.01743619  max resid 0.05100612 
    ## Run 82 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001234864  max resid 0.0002779923 
    ## ... Similar to previous best
    ## Run 83 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 4.67412e-06  max resid 8.824454e-06 
    ## ... Similar to previous best
    ## Run 84 stress 0.08000469 
    ## Run 85 stress 0.2493376 
    ## Run 86 stress 0.07629234 
    ## Run 87 stress 0.07629235 
    ## Run 88 stress 0.08000473 
    ## Run 89 stress 0.07365783 
    ## ... Procrustes: rmse 3.862005e-06  max resid 6.582629e-06 
    ## ... Similar to previous best
    ## Run 90 stress 0.07629238 
    ## Run 91 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001065172  max resid 0.0002493834 
    ## ... Similar to previous best
    ## Run 92 stress 0.07365792 
    ## ... Procrustes: rmse 0.000151719  max resid 0.0003561027 
    ## ... Similar to previous best
    ## Run 93 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744737  max resid 0.05104899 
    ## Run 94 stress 0.07732942 
    ## Run 95 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001170001  max resid 0.0002790555 
    ## ... Similar to previous best
    ## Run 96 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744967  max resid 0.0510643 
    ## Run 97 stress 0.07365783 
    ## ... Procrustes: rmse 4.664507e-06  max resid 8.773073e-06 
    ## ... Similar to previous best
    ## Run 98 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174433  max resid 0.0510138 
    ## Run 99 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001200524  max resid 0.0002832676 
    ## ... Similar to previous best
    ## Run 100 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744635  max resid 0.05104259 
    ## Run 101 stress 0.07629245 
    ## Run 102 stress 0.0828804 
    ## Run 103 stress 0.0800047 
    ## Run 104 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174537  max resid 0.05108813 
    ## Run 105 stress 0.07629236 
    ## Run 106 stress 0.08000474 
    ## Run 107 stress 0.0773293 
    ## Run 108 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744472  max resid 0.05105747 
    ## Run 109 stress 0.07365786 
    ## ... Procrustes: rmse 9.392391e-05  max resid 0.0002159365 
    ## ... Similar to previous best
    ## Run 110 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174498  max resid 0.05106622 
    ## Run 111 stress 0.351322 
    ## Run 112 stress 0.07378232 
    ## ... Procrustes: rmse 0.01744192  max resid 0.05098995 
    ## Run 113 stress 0.08000477 
    ## Run 114 stress 0.0737823 
    ## ... Procrustes: rmse 0.0174419  max resid 0.05099747 
    ## Run 115 stress 0.07629234 
    ## Run 116 stress 0.07629238 
    ## Run 117 stress 0.08000478 
    ## Run 118 stress 0.07365787 
    ## ... Procrustes: rmse 9.947501e-05  max resid 0.0002381218 
    ## ... Similar to previous best
    ## Run 119 stress 0.07629243 
    ## Run 120 stress 0.07365794 
    ## ... Procrustes: rmse 0.0001573244  max resid 0.0003642315 
    ## ... Similar to previous best
    ## Run 121 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744589  max resid 0.05103761 
    ## Run 122 stress 0.07365787 
    ## ... Procrustes: rmse 8.734477e-05  max resid 0.0001982423 
    ## ... Similar to previous best
    ## Run 123 stress 0.08288042 
    ## Run 124 stress 0.0762924 
    ## Run 125 stress 0.08000468 
    ## Run 126 stress 0.07365785 
    ## ... Procrustes: rmse 8.115055e-05  max resid 0.0001898588 
    ## ... Similar to previous best
    ## Run 127 stress 0.07629234 
    ## Run 128 stress 0.258405 
    ## Run 129 stress 0.07629239 
    ## Run 130 stress 0.07365783 
    ## ... Procrustes: rmse 5.16896e-06  max resid 9.177983e-06 
    ## ... Similar to previous best
    ## Run 131 stress 0.07629235 
    ## Run 132 stress 0.07365783 
    ## ... Procrustes: rmse 1.144598e-05  max resid 2.576632e-05 
    ## ... Similar to previous best
    ## Run 133 stress 0.2584052 
    ## Run 134 stress 0.08288037 
    ## Run 135 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001465653  max resid 0.0003485824 
    ## ... Similar to previous best
    ## Run 136 stress 0.07365784 
    ## ... Procrustes: rmse 5.87516e-05  max resid 0.0001389066 
    ## ... Similar to previous best
    ## Run 137 stress 0.07629237 
    ## Run 138 stress 0.07365785 
    ## ... Procrustes: rmse 6.90215e-05  max resid 0.0001624432 
    ## ... Similar to previous best
    ## Run 139 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001236167  max resid 0.0002854281 
    ## ... Similar to previous best
    ## Run 140 stress 0.08000474 
    ## Run 141 stress 0.0762924 
    ## Run 142 stress 0.07365785 
    ## ... Procrustes: rmse 2.551867e-05  max resid 5.571832e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.07365783 
    ## ... Procrustes: rmse 6.517597e-06  max resid 1.448527e-05 
    ## ... Similar to previous best
    ## Run 144 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745154  max resid 0.05106879 
    ## Run 145 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 4.969329e-06  max resid 1.185221e-05 
    ## ... Similar to previous best
    ## Run 146 stress 0.07365783 
    ## ... Procrustes: rmse 7.542183e-06  max resid 1.237828e-05 
    ## ... Similar to previous best
    ## Run 147 stress 0.07365784 
    ## ... Procrustes: rmse 5.463769e-05  max resid 0.0001273842 
    ## ... Similar to previous best
    ## Run 148 stress 0.07732945 
    ## Run 149 stress 0.07629241 
    ## Run 150 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744494  max resid 0.05103384 
    ## Run 151 stress 0.07365783 
    ## ... Procrustes: rmse 1.166532e-05  max resid 2.430768e-05 
    ## ... Similar to previous best
    ## Run 152 stress 0.0800047 
    ## Run 153 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744353  max resid 0.05105905 
    ## Run 154 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744657  max resid 0.05105838 
    ## Run 155 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001434949  max resid 0.0003405436 
    ## ... Similar to previous best
    ## Run 156 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001418216  max resid 0.0003397049 
    ## ... Similar to previous best
    ## Run 157 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744692  max resid 0.05104671 
    ## Run 158 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001185916  max resid 0.0002795242 
    ## ... Similar to previous best
    ## Run 159 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174468  max resid 0.05104655 
    ## Run 160 stress 0.07732935 
    ## Run 161 stress 0.07365784 
    ## ... Procrustes: rmse 4.76709e-05  max resid 0.0001120283 
    ## ... Similar to previous best
    ## Run 162 stress 0.07365784 
    ## ... Procrustes: rmse 2.444494e-05  max resid 4.998375e-05 
    ## ... Similar to previous best
    ## Run 163 stress 0.07365784 
    ## ... Procrustes: rmse 6.156706e-05  max resid 0.000146211 
    ## ... Similar to previous best
    ## Run 164 stress 0.08000476 
    ## Run 165 stress 0.07629234 
    ## Run 166 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744405  max resid 0.05102607 
    ## Run 167 stress 0.07365784 
    ## ... Procrustes: rmse 2.399851e-05  max resid 4.714438e-05 
    ## ... Similar to previous best
    ## Run 168 stress 0.07629233 
    ## Run 169 stress 0.07629237 
    ## Run 170 stress 0.08000475 
    ## Run 171 stress 0.07629237 
    ## Run 172 stress 0.07365783 
    ## ... Procrustes: rmse 1.687006e-05  max resid 4.163124e-05 
    ## ... Similar to previous best
    ## Run 173 stress 0.07629237 
    ## Run 174 stress 0.07732928 
    ## Run 175 stress 0.07365783 
    ## ... Procrustes: rmse 8.615376e-06  max resid 2.004634e-05 
    ## ... Similar to previous best
    ## Run 176 stress 0.07629241 
    ## Run 177 stress 0.08000468 
    ## Run 178 stress 0.07629236 
    ## Run 179 stress 0.07629246 
    ## Run 180 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745395  max resid 0.05108658 
    ## Run 181 stress 0.07629233 
    ## Run 182 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 6.919037e-06  max resid 1.619072e-05 
    ## ... Similar to previous best
    ## Run 183 stress 0.07732929 
    ## Run 184 stress 0.07629239 
    ## Run 185 stress 0.07629236 
    ## Run 186 stress 0.07629244 
    ## Run 187 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001094562  max resid 0.0002579038 
    ## ... Similar to previous best
    ## Run 188 stress 0.07365783 
    ## ... Procrustes: rmse 1.129175e-05  max resid 2.689761e-05 
    ## ... Similar to previous best
    ## Run 189 stress 0.0736579 
    ## ... Procrustes: rmse 0.000121927  max resid 0.0002911861 
    ## ... Similar to previous best
    ## Run 190 stress 0.08000471 
    ## Run 191 stress 0.0737823 
    ## ... Procrustes: rmse 0.01743991  max resid 0.05099382 
    ## Run 192 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745116  max resid 0.05108039 
    ## Run 193 stress 0.08288054 
    ## Run 194 stress 0.07732932 
    ## Run 195 stress 0.08233753 
    ## Run 196 stress 0.07629236 
    ## Run 197 stress 0.07378234 
    ## ... Procrustes: rmse 0.01744671  max resid 0.05108219 
    ## Run 198 stress 0.07365783 
    ## ... Procrustes: rmse 1.049367e-05  max resid 2.11104e-05 
    ## ... Similar to previous best
    ## Run 199 stress 0.07732929 
    ## Run 200 stress 0.07629243 
    ## Run 201 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744876  max resid 0.05106664 
    ## Run 202 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744811  max resid 0.0510562 
    ## Run 203 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744401  max resid 0.05103622 
    ## Run 204 stress 0.07629239 
    ## Run 205 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744222  max resid 0.05105643 
    ## Run 206 stress 0.08000467 
    ## Run 207 stress 0.08288065 
    ## Run 208 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174453  max resid 0.05104399 
    ## Run 209 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744409  max resid 0.05103575 
    ## Run 210 stress 0.07629251 
    ## Run 211 stress 0.07365794 
    ## ... Procrustes: rmse 0.0001795279  max resid 0.0004236844 
    ## ... Similar to previous best
    ## Run 212 stress 0.07629237 
    ## Run 213 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745585  max resid 0.0510853 
    ## Run 214 stress 0.0800047 
    ## Run 215 stress 0.08000471 
    ## Run 216 stress 0.07378235 
    ## ... Procrustes: rmse 0.01743661  max resid 0.05095888 
    ## Run 217 stress 0.08233755 
    ## Run 218 stress 0.07365783 
    ## ... Procrustes: rmse 1.881757e-05  max resid 4.434849e-05 
    ## ... Similar to previous best
    ## Run 219 stress 0.08000468 
    ## Run 220 stress 0.07629246 
    ## Run 221 stress 0.08233757 
    ## Run 222 stress 0.07629233 
    ## Run 223 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744454  max resid 0.0510541 
    ## Run 224 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744285  max resid 0.05105004 
    ## Run 225 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744246  max resid 0.05102872 
    ## Run 226 stress 0.07365786 
    ## ... Procrustes: rmse 5.488597e-05  max resid 0.0001172607 
    ## ... Similar to previous best
    ## Run 227 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001010307  max resid 0.0002339844 
    ## ... Similar to previous best
    ## Run 228 stress 0.08233754 
    ## Run 229 stress 0.07365785 
    ## ... Procrustes: rmse 6.366315e-05  max resid 0.0001520878 
    ## ... Similar to previous best
    ## Run 230 stress 0.07365786 
    ## ... Procrustes: rmse 5.918556e-05  max resid 0.0001365198 
    ## ... Similar to previous best
    ## Run 231 stress 0.08233756 
    ## Run 232 stress 0.07629233 
    ## Run 233 stress 0.07378226 
    ## ... Procrustes: rmse 0.017441  max resid 0.05101944 
    ## Run 234 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744609  max resid 0.05106305 
    ## Run 235 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744464  max resid 0.05104109 
    ## Run 236 stress 0.0800047 
    ## Run 237 stress 0.07365783 
    ## ... Procrustes: rmse 1.192292e-05  max resid 2.809664e-05 
    ## ... Similar to previous best
    ## Run 238 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744253  max resid 0.05103148 
    ## Run 239 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744952  max resid 0.05107425 
    ## Run 240 stress 0.07629235 
    ## Run 241 stress 0.08000476 
    ## Run 242 stress 0.0762924 
    ## Run 243 stress 0.07629233 
    ## Run 244 stress 0.08000482 
    ## Run 245 stress 0.0762924 
    ## Run 246 stress 0.07365783 
    ## ... Procrustes: rmse 7.126942e-06  max resid 1.516083e-05 
    ## ... Similar to previous best
    ## Run 247 stress 0.07365785 
    ## ... Procrustes: rmse 6.249203e-05  max resid 0.0001489225 
    ## ... Similar to previous best
    ## Run 248 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001533985  max resid 0.000361344 
    ## ... Similar to previous best
    ## Run 249 stress 0.08233758 
    ## Run 250 stress 0.07365784 
    ## ... Procrustes: rmse 3.512044e-05  max resid 7.988952e-05 
    ## ... Similar to previous best
    ## Run 251 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001466594  max resid 0.0003465691 
    ## ... Similar to previous best
    ## Run 252 stress 0.07365783 
    ## ... Procrustes: rmse 5.884505e-06  max resid 1.054323e-05 
    ## ... Similar to previous best
    ## Run 253 stress 0.3495749 
    ## Run 254 stress 0.08233753 
    ## Run 255 stress 0.07629233 
    ## Run 256 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744261  max resid 0.05101965 
    ## Run 257 stress 0.08000467 
    ## Run 258 stress 0.08000467 
    ## Run 259 stress 0.07365784 
    ## ... Procrustes: rmse 3.501874e-05  max resid 8.076834e-05 
    ## ... Similar to previous best
    ## Run 260 stress 0.07365783 
    ## ... Procrustes: rmse 8.202387e-06  max resid 1.829017e-05 
    ## ... Similar to previous best
    ## Run 261 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744219  max resid 0.05101649 
    ## Run 262 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 6.624335e-06  max resid 1.529863e-05 
    ## ... Similar to previous best
    ## Run 263 stress 0.07629245 
    ## Run 264 stress 0.07365785 
    ## ... Procrustes: rmse 7.025975e-05  max resid 0.0001648109 
    ## ... Similar to previous best
    ## Run 265 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001220365  max resid 0.0002859986 
    ## ... Similar to previous best
    ## Run 266 stress 0.07629237 
    ## Run 267 stress 0.08000468 
    ## Run 268 stress 0.07365783 
    ## ... Procrustes: rmse 2.795703e-05  max resid 6.499829e-05 
    ## ... Similar to previous best
    ## Run 269 stress 0.07365786 
    ## ... Procrustes: rmse 8.943957e-05  max resid 0.0002098903 
    ## ... Similar to previous best
    ## Run 270 stress 0.08233754 
    ## Run 271 stress 0.07365786 
    ## ... Procrustes: rmse 6.469491e-05  max resid 0.0001412916 
    ## ... Similar to previous best
    ## Run 272 stress 0.07365784 
    ## ... Procrustes: rmse 5.8288e-05  max resid 0.0001371165 
    ## ... Similar to previous best
    ## Run 273 stress 0.07732925 
    ## Run 274 stress 0.08000475 
    ## Run 275 stress 0.07629236 
    ## Run 276 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001051732  max resid 0.000245224 
    ## ... Similar to previous best
    ## Run 277 stress 0.07629245 
    ## Run 278 stress 0.0828804 
    ## Run 279 stress 0.07629236 
    ## Run 280 stress 0.07629233 
    ## Run 281 stress 0.07365785 
    ## ... Procrustes: rmse 6.856648e-05  max resid 0.0001613611 
    ## ... Similar to previous best
    ## Run 282 stress 0.08000474 
    ## Run 283 stress 0.07629245 
    ## Run 284 stress 0.07629234 
    ## Run 285 stress 0.08000468 
    ## Run 286 stress 0.07365785 
    ## ... Procrustes: rmse 6.543843e-05  max resid 0.0001518186 
    ## ... Similar to previous best
    ## Run 287 stress 0.07378234 
    ## ... Procrustes: rmse 0.01744125  max resid 0.0509879 
    ## Run 288 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001128214  max resid 0.0002636034 
    ## ... Similar to previous best
    ## Run 289 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174491  max resid 0.05106317 
    ## Run 290 stress 0.08000475 
    ## Run 291 stress 0.0737823 
    ## ... Procrustes: rmse 0.01745237  max resid 0.05104797 
    ## Run 292 stress 0.07365783 
    ## ... Procrustes: rmse 8.947317e-06  max resid 1.927692e-05 
    ## ... Similar to previous best
    ## Run 293 stress 0.07629238 
    ## Run 294 stress 0.07629235 
    ## Run 295 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001188341  max resid 0.0002801979 
    ## ... Similar to previous best
    ## Run 296 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001170612  max resid 0.0002735255 
    ## ... Similar to previous best
    ## Run 297 stress 0.07629239 
    ## Run 298 stress 0.07365784 
    ## ... Procrustes: rmse 4.248384e-05  max resid 9.943453e-05 
    ## ... Similar to previous best
    ## Run 299 stress 0.07629243 
    ## Run 300 stress 0.08288048 
    ## Run 301 stress 0.0800047 
    ## Run 302 stress 0.07732935 
    ## Run 303 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744704  max resid 0.05105558 
    ## Run 304 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744065  max resid 0.05105179 
    ## Run 305 stress 0.07365786 
    ## ... Procrustes: rmse 8.629392e-05  max resid 0.0002045469 
    ## ... Similar to previous best
    ## Run 306 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001632996  max resid 0.000388695 
    ## ... Similar to previous best
    ## Run 307 stress 0.07378227 
    ## ... Procrustes: rmse 0.017443  max resid 0.0510203 
    ## Run 308 stress 0.07365784 
    ## ... Procrustes: rmse 3.882611e-05  max resid 9.255797e-05 
    ## ... Similar to previous best
    ## Run 309 stress 0.07732942 
    ## Run 310 stress 0.08233767 
    ## Run 311 stress 0.07629237 
    ## Run 312 stress 0.08000469 
    ## Run 313 stress 0.07365785 
    ## ... Procrustes: rmse 8.208772e-05  max resid 0.0001946875 
    ## ... Similar to previous best
    ## Run 314 stress 0.0800047 
    ## Run 315 stress 0.0737823 
    ## ... Procrustes: rmse 0.01745446  max resid 0.0510544 
    ## Run 316 stress 0.07365783 
    ## ... Procrustes: rmse 7.532085e-06  max resid 1.739761e-05 
    ## ... Similar to previous best
    ## Run 317 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744944  max resid 0.0510675 
    ## Run 318 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744746  max resid 0.05105716 
    ## Run 319 stress 0.07365786 
    ## ... Procrustes: rmse 8.733835e-05  max resid 0.0002080219 
    ## ... Similar to previous best
    ## Run 320 stress 0.08000471 
    ## Run 321 stress 0.07629242 
    ## Run 322 stress 0.07629237 
    ## Run 323 stress 0.07629235 
    ## Run 324 stress 0.07378229 
    ## ... Procrustes: rmse 0.01746549  max resid 0.05114849 
    ## Run 325 stress 0.07629241 
    ## Run 326 stress 0.08233751 
    ## Run 327 stress 0.07732936 
    ## Run 328 stress 0.07629234 
    ## Run 329 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744646  max resid 0.05104493 
    ## Run 330 stress 0.07365784 
    ## ... Procrustes: rmse 3.775477e-05  max resid 8.881604e-05 
    ## ... Similar to previous best
    ## Run 331 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744548  max resid 0.05106191 
    ## Run 332 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744295  max resid 0.05105664 
    ## Run 333 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744552  max resid 0.05103476 
    ## Run 334 stress 0.07732931 
    ## Run 335 stress 0.07365784 
    ## ... Procrustes: rmse 1.600149e-05  max resid 2.634551e-05 
    ## ... Similar to previous best
    ## Run 336 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744421  max resid 0.05103257 
    ## Run 337 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744916  max resid 0.05106608 
    ## Run 338 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744616  max resid 0.05105006 
    ## Run 339 stress 0.07629255 
    ## Run 340 stress 0.07365794 
    ## ... Procrustes: rmse 0.0001518544  max resid 0.0003627576 
    ## ... Similar to previous best
    ## Run 341 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744161  max resid 0.05100288 
    ## Run 342 stress 0.07629241 
    ## Run 343 stress 0.07378231 
    ## ... Procrustes: rmse 0.01743985  max resid 0.05098962 
    ## Run 344 stress 0.07629246 
    ## Run 345 stress 0.07365783 
    ## ... Procrustes: rmse 2.10077e-05  max resid 4.896169e-05 
    ## ... Similar to previous best
    ## Run 346 stress 0.07629235 
    ## Run 347 stress 0.08000481 
    ## Run 348 stress 0.08233777 
    ## Run 349 stress 0.07629246 
    ## Run 350 stress 0.07365783 
    ## ... Procrustes: rmse 2.804197e-06  max resid 6.147494e-06 
    ## ... Similar to previous best
    ## Run 351 stress 0.07365783 
    ## ... Procrustes: rmse 2.60009e-05  max resid 5.821547e-05 
    ## ... Similar to previous best
    ## Run 352 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744553  max resid 0.05104622 
    ## Run 353 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744536  max resid 0.05103845 
    ## Run 354 stress 0.0828804 
    ## Run 355 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744474  max resid 0.05103584 
    ## Run 356 stress 0.08000474 
    ## Run 357 stress 0.07732929 
    ## Run 358 stress 0.07629233 
    ## Run 359 stress 0.08000482 
    ## Run 360 stress 0.08233757 
    ## Run 361 stress 0.08000469 
    ## Run 362 stress 0.07629234 
    ## Run 363 stress 0.07629233 
    ## Run 364 stress 0.08000468 
    ## Run 365 stress 0.07629239 
    ## Run 366 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744696  max resid 0.05106173 
    ## Run 367 stress 0.08233759 
    ## Run 368 stress 0.07629249 
    ## Run 369 stress 0.08233761 
    ## Run 370 stress 0.07629237 
    ## Run 371 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744719  max resid 0.0510658 
    ## Run 372 stress 0.07629237 
    ## Run 373 stress 0.0762924 
    ## Run 374 stress 0.07365783 
    ## ... Procrustes: rmse 3.450246e-05  max resid 8.126585e-05 
    ## ... Similar to previous best
    ## Run 375 stress 0.07629238 
    ## Run 376 stress 0.08000468 
    ## Run 377 stress 0.07732927 
    ## Run 378 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744884  max resid 0.05107368 
    ## Run 379 stress 0.07629237 
    ## Run 380 stress 0.07629244 
    ## Run 381 stress 0.07629234 
    ## Run 382 stress 0.08233768 
    ## Run 383 stress 0.08288042 
    ## Run 384 stress 0.08000474 
    ## Run 385 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744795  max resid 0.0510611 
    ## Run 386 stress 0.08000468 
    ## Run 387 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744743  max resid 0.05106017 
    ## Run 388 stress 0.08000473 
    ## Run 389 stress 0.0800047 
    ## Run 390 stress 0.07365783 
    ## ... Procrustes: rmse 2.751334e-05  max resid 6.579832e-05 
    ## ... Similar to previous best
    ## Run 391 stress 0.08000467 
    ## Run 392 stress 0.07629238 
    ## Run 393 stress 0.07378228 
    ## ... Procrustes: rmse 0.01742485  max resid 0.05096855 
    ## Run 394 stress 0.07629236 
    ## Run 395 stress 0.07629235 
    ## Run 396 stress 0.08000472 
    ## Run 397 stress 0.08288039 
    ## Run 398 stress 0.07365783 
    ## ... Procrustes: rmse 9.018111e-06  max resid 2.065983e-05 
    ## ... Similar to previous best
    ## Run 399 stress 0.08288051 
    ## Run 400 stress 0.07629234 
    ## Run 401 stress 0.07732945 
    ## Run 402 stress 0.07629234 
    ## Run 403 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744503  max resid 0.05105662 
    ## Run 404 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001035811  max resid 0.0002457973 
    ## ... Similar to previous best
    ## Run 405 stress 0.08000468 
    ## Run 406 stress 0.07629234 
    ## Run 407 stress 0.08000467 
    ## Run 408 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001141895  max resid 0.0002680424 
    ## ... Similar to previous best
    ## Run 409 stress 0.07378232 
    ## ... Procrustes: rmse 0.01744616  max resid 0.05100994 
    ## Run 410 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745511  max resid 0.05106876 
    ## Run 411 stress 0.07629237 
    ## Run 412 stress 0.07365784 
    ## ... Procrustes: rmse 3.010038e-05  max resid 6.987851e-05 
    ## ... Similar to previous best
    ## Run 413 stress 0.07732928 
    ## Run 414 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744469  max resid 0.05103532 
    ## Run 415 stress 0.07365785 
    ## ... Procrustes: rmse 6.357752e-05  max resid 0.000150665 
    ## ... Similar to previous best
    ## Run 416 stress 0.07365783 
    ## ... Procrustes: rmse 2.879129e-05  max resid 6.398377e-05 
    ## ... Similar to previous best
    ## Run 417 stress 0.0828804 
    ## Run 418 stress 0.07629233 
    ## Run 419 stress 0.07629236 
    ## Run 420 stress 0.07629236 
    ## Run 421 stress 0.07629247 
    ## Run 422 stress 0.08000468 
    ## Run 423 stress 0.07365784 
    ## ... Procrustes: rmse 2.49897e-05  max resid 5.898303e-05 
    ## ... Similar to previous best
    ## Run 424 stress 0.07629242 
    ## Run 425 stress 0.07378232 
    ## ... Procrustes: rmse 0.01743098  max resid 0.0509469 
    ## Run 426 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744414  max resid 0.05103914 
    ## Run 427 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001057087  max resid 0.0002485971 
    ## ... Similar to previous best
    ## Run 428 stress 0.07732925 
    ## Run 429 stress 0.0800047 
    ## Run 430 stress 0.08288067 
    ## Run 431 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001182914  max resid 0.0002802283 
    ## ... Similar to previous best
    ## Run 432 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174479  max resid 0.05105929 
    ## Run 433 stress 0.08000468 
    ## Run 434 stress 0.08000467 
    ## Run 435 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744798  max resid 0.05106244 
    ## Run 436 stress 0.08233758 
    ## Run 437 stress 0.07629236 
    ## Run 438 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001300921  max resid 0.0003092484 
    ## ... Similar to previous best
    ## Run 439 stress 0.08000469 
    ## Run 440 stress 0.07732927 
    ## Run 441 stress 0.08000475 
    ## Run 442 stress 0.0736579 
    ## ... Procrustes: rmse 6.836231e-05  max resid 0.0001591254 
    ## ... Similar to previous best
    ## Run 443 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001006771  max resid 0.0002363569 
    ## ... Similar to previous best
    ## Run 444 stress 0.07629235 
    ## Run 445 stress 0.07732926 
    ## Run 446 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744702  max resid 0.05102475 
    ## Run 447 stress 0.08288048 
    ## Run 448 stress 0.08000477 
    ## Run 449 stress 0.07629235 
    ## Run 450 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174523  max resid 0.05108519 
    ## Run 451 stress 0.07365783 
    ## ... Procrustes: rmse 3.040084e-05  max resid 6.964405e-05 
    ## ... Similar to previous best
    ## Run 452 stress 0.07629251 
    ## Run 453 stress 0.07365783 
    ## ... Procrustes: rmse 7.977659e-06  max resid 1.601505e-05 
    ## ... Similar to previous best
    ## Run 454 stress 0.07365783 
    ## ... Procrustes: rmse 1.796168e-05  max resid 4.261042e-05 
    ## ... Similar to previous best
    ## Run 455 stress 0.08000467 
    ## Run 456 stress 0.07629235 
    ## Run 457 stress 0.07365784 
    ## ... Procrustes: rmse 3.940359e-05  max resid 9.247485e-05 
    ## ... Similar to previous best
    ## Run 458 stress 0.0828804 
    ## Run 459 stress 0.08000467 
    ## Run 460 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744553  max resid 0.0510437 
    ## Run 461 stress 0.07365783 
    ## ... Procrustes: rmse 2.425957e-05  max resid 5.710679e-05 
    ## ... Similar to previous best
    ## Run 462 stress 0.07629236 
    ## Run 463 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174817  max resid 0.05119858 
    ## Run 464 stress 0.07629246 
    ## Run 465 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001231698  max resid 0.0002934641 
    ## ... Similar to previous best
    ## Run 466 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744493  max resid 0.05103629 
    ## Run 467 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744395  max resid 0.05103199 
    ## Run 468 stress 0.07365783 
    ## ... Procrustes: rmse 9.196244e-06  max resid 2.289771e-05 
    ## ... Similar to previous best
    ## Run 469 stress 0.08233755 
    ## Run 470 stress 0.08233772 
    ## Run 471 stress 0.08000471 
    ## Run 472 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001087571  max resid 0.0002540996 
    ## ... Similar to previous best
    ## Run 473 stress 0.0737823 
    ## ... Procrustes: rmse 0.01747107  max resid 0.0511302 
    ## Run 474 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745049  max resid 0.05107401 
    ## Run 475 stress 0.07629235 
    ## Run 476 stress 0.08233761 
    ## Run 477 stress 0.08000468 
    ## Run 478 stress 0.07629243 
    ## Run 479 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744698  max resid 0.05105507 
    ## Run 480 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001176551  max resid 0.0002811654 
    ## ... Similar to previous best
    ## Run 481 stress 0.07365784 
    ## ... Procrustes: rmse 2.073915e-05  max resid 4.822799e-05 
    ## ... Similar to previous best
    ## Run 482 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744245  max resid 0.05100717 
    ## Run 483 stress 0.08000469 
    ## Run 484 stress 0.07365784 
    ## ... Procrustes: rmse 2.803284e-05  max resid 5.828017e-05 
    ## ... Similar to previous best
    ## Run 485 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001672013  max resid 0.000395846 
    ## ... Similar to previous best
    ## Run 486 stress 0.07365785 
    ## ... Procrustes: rmse 6.384515e-05  max resid 0.0001484342 
    ## ... Similar to previous best
    ## Run 487 stress 0.07629243 
    ## Run 488 stress 0.08233762 
    ## Run 489 stress 0.07629242 
    ## Run 490 stress 0.07629243 
    ## Run 491 stress 0.08000469 
    ## Run 492 stress 0.07732932 
    ## Run 493 stress 0.08233774 
    ## Run 494 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174377  max resid 0.05100941 
    ## Run 495 stress 0.07629242 
    ## Run 496 stress 0.2667134 
    ## Run 497 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744926  max resid 0.05106919 
    ## Run 498 stress 0.07365783 
    ## ... Procrustes: rmse 2.238644e-05  max resid 5.292541e-05 
    ## ... Similar to previous best
    ## Run 499 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001311601  max resid 0.0003049335 
    ## ... Similar to previous best
    ## Run 500 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744492  max resid 0.05103936 
    ## *** Best solution repeated 56 times

``` r
# Stratified lakes and ocean sites
SD_beta_SO_NMDS <- metaMDS(SD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07970526 
    ## Run 1 stress 0.06942777 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1019668  max resid 0.263953 
    ## Run 2 stress 0.07428316 
    ## Run 3 stress 0.08340287 
    ## Run 4 stress 0.08448452 
    ## Run 5 stress 0.07428318 
    ## Run 6 stress 0.07844912 
    ## Run 7 stress 0.08340288 
    ## Run 8 stress 0.07250816 
    ## Run 9 stress 0.07250812 
    ## Run 10 stress 0.08340288 
    ## Run 11 stress 0.07970528 
    ## Run 12 stress 0.07428313 
    ## Run 13 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323924  max resid 0.03330029 
    ## Run 14 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324408  max resid 0.03330975 
    ## Run 15 stress 0.07970525 
    ## Run 16 stress 0.07970525 
    ## Run 17 stress 0.0834029 
    ## Run 18 stress 0.07428315 
    ## Run 19 stress 0.07428313 
    ## Run 20 stress 0.08340295 
    ## Run 21 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326797  max resid 0.03335856 
    ## Run 22 stress 0.07428316 
    ## Run 23 stress 0.08448435 
    ## Run 24 stress 0.07428313 
    ## Run 25 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.624082e-05  max resid 7.801582e-05 
    ## ... Similar to previous best
    ## Run 26 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326499  max resid 0.03334847 
    ## Run 27 stress 0.08448442 
    ## Run 28 stress 0.08340288 
    ## Run 29 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324847  max resid 0.03331409 
    ## Run 30 stress 0.07428313 
    ## Run 31 stress 0.07970525 
    ## Run 32 stress 0.08340287 
    ## Run 33 stress 0.07970525 
    ## Run 34 stress 0.06978195 
    ## ... Procrustes: rmse 0.01325717  max resid 0.03333395 
    ## Run 35 stress 0.07428314 
    ## Run 36 stress 0.07428318 
    ## Run 37 stress 0.08340286 
    ## Run 38 stress 0.07428316 
    ## Run 39 stress 0.08340286 
    ## Run 40 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325443  max resid 0.03332119 
    ## Run 41 stress 0.08340287 
    ## Run 42 stress 0.08340287 
    ## Run 43 stress 0.07428318 
    ## Run 44 stress 0.07428313 
    ## Run 45 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323214  max resid 0.03327984 
    ## Run 46 stress 0.07250812 
    ## Run 47 stress 0.07428316 
    ## Run 48 stress 0.07428315 
    ## Run 49 stress 0.07428316 
    ## Run 50 stress 0.07428315 
    ## Run 51 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325662  max resid 0.03333088 
    ## Run 52 stress 0.07844949 
    ## Run 53 stress 0.07428313 
    ## Run 54 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 8.587129e-06  max resid 2.517885e-05 
    ## ... Similar to previous best
    ## Run 55 stress 0.06942776 
    ## ... Procrustes: rmse 7.517031e-06  max resid 1.983512e-05 
    ## ... Similar to previous best
    ## Run 56 stress 0.07970525 
    ## Run 57 stress 0.07428314 
    ## Run 58 stress 0.08340287 
    ## Run 59 stress 0.07428317 
    ## Run 60 stress 0.07844915 
    ## Run 61 stress 0.06942776 
    ## ... Procrustes: rmse 2.622287e-05  max resid 6.786384e-05 
    ## ... Similar to previous best
    ## Run 62 stress 0.07428312 
    ## Run 63 stress 0.08340291 
    ## Run 64 stress 0.07428314 
    ## Run 65 stress 0.07970525 
    ## Run 66 stress 0.06942779 
    ## ... Procrustes: rmse 2.067544e-05  max resid 6.625564e-05 
    ## ... Similar to previous best
    ## Run 67 stress 0.07428313 
    ## Run 68 stress 0.08340287 
    ## Run 69 stress 0.07970525 
    ## Run 70 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327457  max resid 0.0333681 
    ## Run 71 stress 0.0697819 
    ## ... Procrustes: rmse 0.013214  max resid 0.03324267 
    ## Run 72 stress 0.07970526 
    ## Run 73 stress 0.07250812 
    ## Run 74 stress 0.06978193 
    ## ... Procrustes: rmse 0.01318834  max resid 0.03318716 
    ## Run 75 stress 0.07970526 
    ## Run 76 stress 0.08340287 
    ## Run 77 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321216  max resid 0.03323928 
    ## Run 78 stress 0.07428313 
    ## Run 79 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325798  max resid 0.03333415 
    ## Run 80 stress 0.0834029 
    ## Run 81 stress 0.07970525 
    ## Run 82 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323806  max resid 0.03329471 
    ## Run 83 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326658  max resid 0.0333531 
    ## Run 84 stress 0.07250813 
    ## Run 85 stress 0.08448451 
    ## Run 86 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325181  max resid 0.03332212 
    ## Run 87 stress 0.07250812 
    ## Run 88 stress 0.07970526 
    ## Run 89 stress 0.06942776 
    ## ... Procrustes: rmse 1.608885e-05  max resid 4.18703e-05 
    ## ... Similar to previous best
    ## Run 90 stress 0.08340287 
    ## Run 91 stress 0.07428316 
    ## Run 92 stress 0.07428313 
    ## Run 93 stress 0.07970525 
    ## Run 94 stress 0.08340288 
    ## Run 95 stress 0.07428313 
    ## Run 96 stress 0.07970526 
    ## Run 97 stress 0.07970525 
    ## Run 98 stress 0.08340286 
    ## Run 99 stress 0.06942776 
    ## ... Procrustes: rmse 2.713788e-05  max resid 7.029711e-05 
    ## ... Similar to previous best
    ## Run 100 stress 0.07428313 
    ## Run 101 stress 0.08340291 
    ## Run 102 stress 0.07428313 
    ## Run 103 stress 0.07970526 
    ## Run 104 stress 0.08340294 
    ## Run 105 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 5.559886e-06  max resid 1.337174e-05 
    ## ... Similar to previous best
    ## Run 106 stress 0.07428314 
    ## Run 107 stress 0.08448446 
    ## Run 108 stress 0.07844959 
    ## Run 109 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327718  max resid 0.03337777 
    ## Run 110 stress 0.07970525 
    ## Run 111 stress 0.08340287 
    ## Run 112 stress 0.07250813 
    ## Run 113 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324536  max resid 0.03331434 
    ## Run 114 stress 0.07970526 
    ## Run 115 stress 0.07250813 
    ## Run 116 stress 0.0834029 
    ## Run 117 stress 0.07250812 
    ## Run 118 stress 0.07970526 
    ## Run 119 stress 0.08340286 
    ## Run 120 stress 0.07250813 
    ## Run 121 stress 0.07970526 
    ## Run 122 stress 0.06942776 
    ## ... Procrustes: rmse 7.349461e-06  max resid 2.043454e-05 
    ## ... Similar to previous best
    ## Run 123 stress 0.07844921 
    ## Run 124 stress 0.08340287 
    ## Run 125 stress 0.06942778 
    ## ... Procrustes: rmse 6.599325e-05  max resid 0.0001689805 
    ## ... Similar to previous best
    ## Run 126 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326243  max resid 0.03334877 
    ## Run 127 stress 0.06942776 
    ## ... Procrustes: rmse 2.550443e-05  max resid 6.452845e-05 
    ## ... Similar to previous best
    ## Run 128 stress 0.06942776 
    ## ... Procrustes: rmse 8.16872e-06  max resid 1.823796e-05 
    ## ... Similar to previous best
    ## Run 129 stress 0.06978197 
    ## ... Procrustes: rmse 0.01317329  max resid 0.03315938 
    ## Run 130 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327841  max resid 0.03338151 
    ## Run 131 stress 0.06942776 
    ## ... Procrustes: rmse 9.139576e-06  max resid 2.337092e-05 
    ## ... Similar to previous best
    ## Run 132 stress 0.06942776 
    ## ... Procrustes: rmse 5.058782e-06  max resid 1.025097e-05 
    ## ... Similar to previous best
    ## Run 133 stress 0.06978199 
    ## ... Procrustes: rmse 0.01317283  max resid 0.03315488 
    ## Run 134 stress 0.07428315 
    ## Run 135 stress 0.08340292 
    ## Run 136 stress 0.08340296 
    ## Run 137 stress 0.07428316 
    ## Run 138 stress 0.06942776 
    ## ... Procrustes: rmse 7.229809e-06  max resid 1.621379e-05 
    ## ... Similar to previous best
    ## Run 139 stress 0.08340292 
    ## Run 140 stress 0.07428313 
    ## Run 141 stress 0.07250812 
    ## Run 142 stress 0.07428317 
    ## Run 143 stress 0.07428313 
    ## Run 144 stress 0.06942776 
    ## ... Procrustes: rmse 1.300713e-05  max resid 3.357381e-05 
    ## ... Similar to previous best
    ## Run 145 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 4.053859e-06  max resid 8.803975e-06 
    ## ... Similar to previous best
    ## Run 146 stress 0.06942777 
    ## ... Procrustes: rmse 4.381778e-05  max resid 0.00011277 
    ## ... Similar to previous best
    ## Run 147 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321995  max resid 0.03325536 
    ## Run 148 stress 0.0834029 
    ## Run 149 stress 0.07428316 
    ## Run 150 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321343  max resid 0.03324134 
    ## Run 151 stress 0.07250812 
    ## Run 152 stress 0.07428313 
    ## Run 153 stress 0.06942776 
    ## ... Procrustes: rmse 2.963867e-05  max resid 7.636985e-05 
    ## ... Similar to previous best
    ## Run 154 stress 0.07970525 
    ## Run 155 stress 0.07428314 
    ## Run 156 stress 0.08340289 
    ## Run 157 stress 0.07970526 
    ## Run 158 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327019  max resid 0.03336302 
    ## Run 159 stress 0.07970525 
    ## Run 160 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324796  max resid 0.03331304 
    ## Run 161 stress 0.07250812 
    ## Run 162 stress 0.06978193 
    ## ... Procrustes: rmse 0.01319165  max resid 0.03319362 
    ## Run 163 stress 0.08340287 
    ## Run 164 stress 0.06942776 
    ## ... Procrustes: rmse 1.699968e-05  max resid 4.51217e-05 
    ## ... Similar to previous best
    ## Run 165 stress 0.07970526 
    ## Run 166 stress 0.08448443 
    ## Run 167 stress 0.07428317 
    ## Run 168 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325857  max resid 0.03333804 
    ## Run 169 stress 0.08448441 
    ## Run 170 stress 0.06942776 
    ## ... Procrustes: rmse 1.863897e-05  max resid 5.015341e-05 
    ## ... Similar to previous best
    ## Run 171 stress 0.06978203 
    ## ... Procrustes: rmse 0.01328635  max resid 0.03339705 
    ## Run 172 stress 0.06978202 
    ## ... Procrustes: rmse 0.01324468  max resid 0.0333013 
    ## Run 173 stress 0.06942776 
    ## ... Procrustes: rmse 2.888498e-05  max resid 7.467683e-05 
    ## ... Similar to previous best
    ## Run 174 stress 0.07250814 
    ## Run 175 stress 0.08340287 
    ## Run 176 stress 0.08448452 
    ## Run 177 stress 0.07250812 
    ## Run 178 stress 0.08448446 
    ## Run 179 stress 0.08340287 
    ## Run 180 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326888  max resid 0.03335681 
    ## Run 181 stress 0.07250814 
    ## Run 182 stress 0.06942777 
    ## ... Procrustes: rmse 4.80182e-05  max resid 0.0001235303 
    ## ... Similar to previous best
    ## Run 183 stress 0.07428315 
    ## Run 184 stress 0.07250815 
    ## Run 185 stress 0.07250814 
    ## Run 186 stress 0.07844925 
    ## Run 187 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325537  max resid 0.03332907 
    ## Run 188 stress 0.08340286 
    ## Run 189 stress 0.08340286 
    ## Run 190 stress 0.07428313 
    ## Run 191 stress 0.07250813 
    ## Run 192 stress 0.07428315 
    ## Run 193 stress 0.07428315 
    ## Run 194 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 5.536109e-07  max resid 1.507749e-06 
    ## ... Similar to previous best
    ## Run 195 stress 0.06942776 
    ## ... Procrustes: rmse 2.461821e-06  max resid 6.835765e-06 
    ## ... Similar to previous best
    ## Run 196 stress 0.07428315 
    ## Run 197 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320904  max resid 0.03323304 
    ## Run 198 stress 0.07428313 
    ## Run 199 stress 0.07970525 
    ## Run 200 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321151  max resid 0.03323637 
    ## Run 201 stress 0.07428312 
    ## Run 202 stress 0.06942777 
    ## ... Procrustes: rmse 3.626514e-05  max resid 9.388701e-05 
    ## ... Similar to previous best
    ## Run 203 stress 0.06978193 
    ## ... Procrustes: rmse 0.01318537  max resid 0.03318211 
    ## Run 204 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322417  max resid 0.03326787 
    ## Run 205 stress 0.06978192 
    ## ... Procrustes: rmse 0.0132491  max resid 0.03331842 
    ## Run 206 stress 0.07970525 
    ## Run 207 stress 0.08448441 
    ## Run 208 stress 0.07428313 
    ## Run 209 stress 0.07428313 
    ## Run 210 stress 0.08340287 
    ## Run 211 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324989  max resid 0.03331835 
    ## Run 212 stress 0.08340288 
    ## Run 213 stress 0.07250812 
    ## Run 214 stress 0.06942777 
    ## ... Procrustes: rmse 4.126318e-05  max resid 0.0001065618 
    ## ... Similar to previous best
    ## Run 215 stress 0.08340295 
    ## Run 216 stress 0.07250813 
    ## Run 217 stress 0.08340287 
    ## Run 218 stress 0.07970526 
    ## Run 219 stress 0.08340288 
    ## Run 220 stress 0.08448439 
    ## Run 221 stress 0.08340287 
    ## Run 222 stress 0.07250812 
    ## Run 223 stress 0.07844948 
    ## Run 224 stress 0.07428314 
    ## Run 225 stress 0.08340286 
    ## Run 226 stress 0.07428317 
    ## Run 227 stress 0.06942776 
    ## ... Procrustes: rmse 4.216241e-06  max resid 8.868934e-06 
    ## ... Similar to previous best
    ## Run 228 stress 0.07970526 
    ## Run 229 stress 0.07428318 
    ## Run 230 stress 0.07250814 
    ## Run 231 stress 0.08340289 
    ## Run 232 stress 0.08340287 
    ## Run 233 stress 0.08340287 
    ## Run 234 stress 0.07970525 
    ## Run 235 stress 0.07428313 
    ## Run 236 stress 0.07970526 
    ## Run 237 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321357  max resid 0.03324298 
    ## Run 238 stress 0.08340299 
    ## Run 239 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319413  max resid 0.03319946 
    ## Run 240 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321745  max resid 0.03325066 
    ## Run 241 stress 0.07970526 
    ## Run 242 stress 0.08448446 
    ## Run 243 stress 0.07428313 
    ## Run 244 stress 0.08448436 
    ## Run 245 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326836  max resid 0.03335496 
    ## Run 246 stress 0.07970526 
    ## Run 247 stress 0.07428313 
    ## Run 248 stress 0.08340286 
    ## Run 249 stress 0.06942776 
    ## ... Procrustes: rmse 5.599354e-06  max resid 1.414811e-05 
    ## ... Similar to previous best
    ## Run 250 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320416  max resid 0.0332221 
    ## Run 251 stress 0.07250814 
    ## Run 252 stress 0.08448445 
    ## Run 253 stress 0.07970526 
    ## Run 254 stress 0.07428313 
    ## Run 255 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325773  max resid 0.03333444 
    ## Run 256 stress 0.06942777 
    ## ... Procrustes: rmse 3.672757e-05  max resid 9.493083e-05 
    ## ... Similar to previous best
    ## Run 257 stress 0.07428314 
    ## Run 258 stress 0.07428316 
    ## Run 259 stress 0.07428313 
    ## Run 260 stress 0.07970526 
    ## Run 261 stress 0.07250812 
    ## Run 262 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324184  max resid 0.0333009 
    ## Run 263 stress 0.08340291 
    ## Run 264 stress 0.06942777 
    ## ... Procrustes: rmse 3.58043e-05  max resid 9.201142e-05 
    ## ... Similar to previous best
    ## Run 265 stress 0.07970525 
    ## Run 266 stress 0.07428318 
    ## Run 267 stress 0.07250814 
    ## Run 268 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323994  max resid 0.03329913 
    ## Run 269 stress 0.07428313 
    ## Run 270 stress 0.08340297 
    ## Run 271 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319809  max resid 0.03320893 
    ## Run 272 stress 0.07428313 
    ## Run 273 stress 0.08448446 
    ## Run 274 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326567  max resid 0.03335272 
    ## Run 275 stress 0.08340287 
    ## Run 276 stress 0.07428313 
    ## Run 277 stress 0.07970525 
    ## Run 278 stress 0.07250814 
    ## Run 279 stress 0.06942776 
    ## ... Procrustes: rmse 1.786494e-05  max resid 4.681258e-05 
    ## ... Similar to previous best
    ## Run 280 stress 0.07428318 
    ## Run 281 stress 0.08448445 
    ## Run 282 stress 0.06942777 
    ## ... Procrustes: rmse 4.803137e-05  max resid 0.0001232774 
    ## ... Similar to previous best
    ## Run 283 stress 0.07970526 
    ## Run 284 stress 0.07428314 
    ## Run 285 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132402  max resid 0.03329942 
    ## Run 286 stress 0.08340294 
    ## Run 287 stress 0.07428314 
    ## Run 288 stress 0.08340287 
    ## Run 289 stress 0.06978193 
    ## ... Procrustes: rmse 0.01324617  max resid 0.03331437 
    ## Run 290 stress 0.08340287 
    ## Run 291 stress 0.06978199 
    ## ... Procrustes: rmse 0.01324326  max resid 0.03329983 
    ## Run 292 stress 0.08340289 
    ## Run 293 stress 0.06978199 
    ## ... Procrustes: rmse 0.01327699  max resid 0.0333735 
    ## Run 294 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322162  max resid 0.03326101 
    ## Run 295 stress 0.07428318 
    ## Run 296 stress 0.083403 
    ## Run 297 stress 0.07428314 
    ## Run 298 stress 0.07844919 
    ## Run 299 stress 0.08340287 
    ## Run 300 stress 0.07970525 
    ## Run 301 stress 0.0834029 
    ## Run 302 stress 0.08448449 
    ## Run 303 stress 0.06942777 
    ## ... Procrustes: rmse 4.381243e-05  max resid 0.0001124761 
    ## ... Similar to previous best
    ## Run 304 stress 0.08448439 
    ## Run 305 stress 0.07428313 
    ## Run 306 stress 0.069782 
    ## ... Procrustes: rmse 0.01328302  max resid 0.03338542 
    ## Run 307 stress 0.07844949 
    ## Run 308 stress 0.08340287 
    ## Run 309 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321649  max resid 0.03324588 
    ## Run 310 stress 0.07970525 
    ## Run 311 stress 0.06942776 
    ## ... Procrustes: rmse 2.547472e-05  max resid 6.590207e-05 
    ## ... Similar to previous best
    ## Run 312 stress 0.07428313 
    ## Run 313 stress 0.07250813 
    ## Run 314 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324785  max resid 0.03331918 
    ## Run 315 stress 0.07428314 
    ## Run 316 stress 0.07428313 
    ## Run 317 stress 0.07250812 
    ## Run 318 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325776  max resid 0.03333616 
    ## Run 319 stress 0.08340288 
    ## Run 320 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324344  max resid 0.03330505 
    ## Run 321 stress 0.07250817 
    ## Run 322 stress 0.08448446 
    ## Run 323 stress 0.07250814 
    ## Run 324 stress 0.07970525 
    ## Run 325 stress 0.08340288 
    ## Run 326 stress 0.0834029 
    ## Run 327 stress 0.07428316 
    ## Run 328 stress 0.06942777 
    ## ... Procrustes: rmse 1.412769e-05  max resid 2.491279e-05 
    ## ... Similar to previous best
    ## Run 329 stress 0.06942776 
    ## ... Procrustes: rmse 1.003291e-05  max resid 2.521507e-05 
    ## ... Similar to previous best
    ## Run 330 stress 0.08340287 
    ## Run 331 stress 0.06942776 
    ## ... Procrustes: rmse 2.005675e-05  max resid 5.153888e-05 
    ## ... Similar to previous best
    ## Run 332 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324516  max resid 0.03330825 
    ## Run 333 stress 0.06942776 
    ## ... Procrustes: rmse 8.689986e-06  max resid 2.211812e-05 
    ## ... Similar to previous best
    ## Run 334 stress 0.08448447 
    ## Run 335 stress 0.06942776 
    ## ... Procrustes: rmse 1.24918e-05  max resid 3.103361e-05 
    ## ... Similar to previous best
    ## Run 336 stress 0.08448435 
    ## Run 337 stress 0.07844926 
    ## Run 338 stress 0.07250813 
    ## Run 339 stress 0.07428313 
    ## Run 340 stress 0.07250814 
    ## Run 341 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323017  max resid 0.03327953 
    ## Run 342 stress 0.07844914 
    ## Run 343 stress 0.07428317 
    ## Run 344 stress 0.07970525 
    ## Run 345 stress 0.07250812 
    ## Run 346 stress 0.07970526 
    ## Run 347 stress 0.06942776 
    ## ... Procrustes: rmse 1.134177e-05  max resid 2.946696e-05 
    ## ... Similar to previous best
    ## Run 348 stress 0.08340291 
    ## Run 349 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323504  max resid 0.03328725 
    ## Run 350 stress 0.06942776 
    ## ... Procrustes: rmse 2.226025e-05  max resid 5.95737e-05 
    ## ... Similar to previous best
    ## Run 351 stress 0.07428312 
    ## Run 352 stress 0.07970526 
    ## Run 353 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323374  max resid 0.03328504 
    ## Run 354 stress 0.08340293 
    ## Run 355 stress 0.07428314 
    ## Run 356 stress 0.06942776 
    ## ... Procrustes: rmse 7.012771e-06  max resid 1.753824e-05 
    ## ... Similar to previous best
    ## Run 357 stress 0.07250814 
    ## Run 358 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322835  max resid 0.03327324 
    ## Run 359 stress 0.08340287 
    ## Run 360 stress 0.06942776 
    ## ... Procrustes: rmse 4.841631e-06  max resid 1.274427e-05 
    ## ... Similar to previous best
    ## Run 361 stress 0.07428314 
    ## Run 362 stress 0.07970525 
    ## Run 363 stress 0.08340286 
    ## Run 364 stress 0.08340293 
    ## Run 365 stress 0.08340286 
    ## Run 366 stress 0.07428312 
    ## Run 367 stress 0.08340286 
    ## Run 368 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324762  max resid 0.03331343 
    ## Run 369 stress 0.07970525 
    ## Run 370 stress 0.06942776 
    ## ... Procrustes: rmse 6.336427e-06  max resid 1.623071e-05 
    ## ... Similar to previous best
    ## Run 371 stress 0.08340293 
    ## Run 372 stress 0.07428316 
    ## Run 373 stress 0.07970525 
    ## Run 374 stress 0.08340287 
    ## Run 375 stress 0.07428318 
    ## Run 376 stress 0.08340287 
    ## Run 377 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323658  max resid 0.03328902 
    ## Run 378 stress 0.08340287 
    ## Run 379 stress 0.07970525 
    ## Run 380 stress 0.06942779 
    ## ... Procrustes: rmse 3.753751e-05  max resid 0.000103687 
    ## ... Similar to previous best
    ## Run 381 stress 0.07428313 
    ## Run 382 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324707  max resid 0.03331276 
    ## Run 383 stress 0.08340287 
    ## Run 384 stress 0.06942776 
    ## ... Procrustes: rmse 1.07855e-05  max resid 2.802407e-05 
    ## ... Similar to previous best
    ## Run 385 stress 0.08448435 
    ## Run 386 stress 0.08340288 
    ## Run 387 stress 0.07250813 
    ## Run 388 stress 0.07250815 
    ## Run 389 stress 0.07428314 
    ## Run 390 stress 0.07428315 
    ## Run 391 stress 0.08340287 
    ## Run 392 stress 0.07428315 
    ## Run 393 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325555  max resid 0.03332996 
    ## Run 394 stress 0.07970526 
    ## Run 395 stress 0.08340286 
    ## Run 396 stress 0.07250813 
    ## Run 397 stress 0.08448439 
    ## Run 398 stress 0.07428319 
    ## Run 399 stress 0.07250812 
    ## Run 400 stress 0.07844927 
    ## Run 401 stress 0.07428317 
    ## Run 402 stress 0.07250812 
    ## Run 403 stress 0.07970525 
    ## Run 404 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319971  max resid 0.03321328 
    ## Run 405 stress 0.06978192 
    ## ... Procrustes: rmse 0.01321475  max resid 0.03324305 
    ## Run 406 stress 0.0784494 
    ## Run 407 stress 0.06978194 
    ## ... Procrustes: rmse 0.01318332  max resid 0.03317795 
    ## Run 408 stress 0.08340295 
    ## Run 409 stress 0.07428313 
    ## Run 410 stress 0.07250812 
    ## Run 411 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323672  max resid 0.03329139 
    ## Run 412 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326096  max resid 0.03334181 
    ## Run 413 stress 0.08340287 
    ## Run 414 stress 0.07428315 
    ## Run 415 stress 0.06942776 
    ## ... Procrustes: rmse 1.287835e-05  max resid 3.292356e-05 
    ## ... Similar to previous best
    ## Run 416 stress 0.07428318 
    ## Run 417 stress 0.08340287 
    ## Run 418 stress 0.07428316 
    ## Run 419 stress 0.07250812 
    ## Run 420 stress 0.07250814 
    ## Run 421 stress 0.07428314 
    ## Run 422 stress 0.07970526 
    ## Run 423 stress 0.06942776 
    ## ... Procrustes: rmse 3.789332e-06  max resid 9.480183e-06 
    ## ... Similar to previous best
    ## Run 424 stress 0.07970525 
    ## Run 425 stress 0.07970525 
    ## Run 426 stress 0.06978196 
    ## ... Procrustes: rmse 0.01327038  max resid 0.03335909 
    ## Run 427 stress 0.07970525 
    ## Run 428 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323015  max resid 0.03327689 
    ## Run 429 stress 0.06978191 
    ## ... Procrustes: rmse 0.0131977  max resid 0.03320763 
    ## Run 430 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324753  max resid 0.03331556 
    ## Run 431 stress 0.06978193 
    ## ... Procrustes: rmse 0.01324446  max resid 0.03331076 
    ## Run 432 stress 0.07844918 
    ## Run 433 stress 0.0834029 
    ## Run 434 stress 0.08448452 
    ## Run 435 stress 0.07428313 
    ## Run 436 stress 0.07250813 
    ## Run 437 stress 0.07250812 
    ## Run 438 stress 0.07250813 
    ## Run 439 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326703  max resid 0.03335423 
    ## Run 440 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325751  max resid 0.03333592 
    ## Run 441 stress 0.07428313 
    ## Run 442 stress 0.06978191 
    ## ... Procrustes: rmse 0.01321954  max resid 0.03325236 
    ## Run 443 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324224  max resid 0.03330371 
    ## Run 444 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325275  max resid 0.03332575 
    ## Run 445 stress 0.07844919 
    ## Run 446 stress 0.07844943 
    ## Run 447 stress 0.08340287 
    ## Run 448 stress 0.07250815 
    ## Run 449 stress 0.06978199 
    ## ... Procrustes: rmse 0.01316701  max resid 0.03314269 
    ## Run 450 stress 0.07428316 
    ## Run 451 stress 0.06942776 
    ## ... Procrustes: rmse 4.409149e-06  max resid 7.823426e-06 
    ## ... Similar to previous best
    ## Run 452 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321564  max resid 0.03324633 
    ## Run 453 stress 0.07428313 
    ## Run 454 stress 0.07428319 
    ## Run 455 stress 0.07428314 
    ## Run 456 stress 0.06942776 
    ## ... Procrustes: rmse 1.163498e-06  max resid 2.662396e-06 
    ## ... Similar to previous best
    ## Run 457 stress 0.07970526 
    ## Run 458 stress 0.07428318 
    ## Run 459 stress 0.07428313 
    ## Run 460 stress 0.06942776 
    ## ... Procrustes: rmse 1.600598e-05  max resid 4.280505e-05 
    ## ... Similar to previous best
    ## Run 461 stress 0.07970525 
    ## Run 462 stress 0.07250812 
    ## Run 463 stress 0.08340287 
    ## Run 464 stress 0.08340288 
    ## Run 465 stress 0.06942776 
    ## ... Procrustes: rmse 2.202412e-05  max resid 5.655757e-05 
    ## ... Similar to previous best
    ## Run 466 stress 0.07250815 
    ## Run 467 stress 0.07250812 
    ## Run 468 stress 0.06978193 
    ## ... Procrustes: rmse 0.01319573  max resid 0.03320072 
    ## Run 469 stress 0.06978194 
    ## ... Procrustes: rmse 0.01318878  max resid 0.03319184 
    ## Run 470 stress 0.08340287 
    ## Run 471 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322493  max resid 0.03326654 
    ## Run 472 stress 0.06942776 
    ## ... Procrustes: rmse 5.086941e-06  max resid 1.331977e-05 
    ## ... Similar to previous best
    ## Run 473 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326054  max resid 0.03334038 
    ## Run 474 stress 0.08448435 
    ## Run 475 stress 0.07250812 
    ## Run 476 stress 0.08340288 
    ## Run 477 stress 0.07970525 
    ## Run 478 stress 0.07970525 
    ## Run 479 stress 0.07970526 
    ## Run 480 stress 0.07428314 
    ## Run 481 stress 0.07428313 
    ## Run 482 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325503  max resid 0.03332968 
    ## Run 483 stress 0.08448457 
    ## Run 484 stress 0.07250812 
    ## Run 485 stress 0.07250813 
    ## Run 486 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325366  max resid 0.03332848 
    ## Run 487 stress 0.07250812 
    ## Run 488 stress 0.07428313 
    ## Run 489 stress 0.07970526 
    ## Run 490 stress 0.07428313 
    ## Run 491 stress 0.08340289 
    ## Run 492 stress 0.07970525 
    ## Run 493 stress 0.07428316 
    ## Run 494 stress 0.08340287 
    ## Run 495 stress 0.06942776 
    ## ... Procrustes: rmse 3.502174e-06  max resid 5.616962e-06 
    ## ... Similar to previous best
    ## Run 496 stress 0.06942776 
    ## ... Procrustes: rmse 5.761444e-06  max resid 1.445205e-05 
    ## ... Similar to previous best
    ## Run 497 stress 0.08340289 
    ## Run 498 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324001  max resid 0.033298 
    ## Run 499 stress 0.07428313 
    ## Run 500 stress 0.08448439 
    ## *** Best solution repeated 33 times

``` r
# Stratified lakes
SD_beta_S_NMDS <- metaMDS(SD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1410449 
    ## Run 1 stress 0.1415299 
    ## ... Procrustes: rmse 0.1479841  max resid 0.2321199 
    ## Run 2 stress 0.1415299 
    ## ... Procrustes: rmse 0.1479839  max resid 0.23212 
    ## Run 3 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 0.2673125  max resid 0.5425959 
    ## Run 4 stress 0.1320258 
    ## ... Procrustes: rmse 2.327379e-06  max resid 4.631148e-06 
    ## ... Similar to previous best
    ## Run 5 stress 0.1598154 
    ## Run 6 stress 0.1415299 
    ## Run 7 stress 0.1998955 
    ## Run 8 stress 0.1407298 
    ## Run 9 stress 0.2008162 
    ## Run 10 stress 0.2842805 
    ## Run 11 stress 0.2385029 
    ## Run 12 stress 0.1407207 
    ## Run 13 stress 0.2074935 
    ## Run 14 stress 0.1415299 
    ## Run 15 stress 0.2008162 
    ## Run 16 stress 0.1410449 
    ## Run 17 stress 0.1693717 
    ## Run 18 stress 0.2520602 
    ## Run 19 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 4.420276e-06  max resid 9.016586e-06 
    ## ... Similar to previous best
    ## Run 20 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 1.796548e-06  max resid 3.617845e-06 
    ## ... Similar to previous best
    ## Run 21 stress 0.2015531 
    ## Run 22 stress 0.1320258 
    ## ... Procrustes: rmse 6.093099e-07  max resid 8.873662e-07 
    ## ... Similar to previous best
    ## Run 23 stress 0.1663544 
    ## Run 24 stress 0.1445811 
    ## Run 25 stress 0.1962476 
    ## Run 26 stress 0.1320258 
    ## ... Procrustes: rmse 1.67212e-06  max resid 3.514294e-06 
    ## ... Similar to previous best
    ## Run 27 stress 0.1551863 
    ## Run 28 stress 0.1415299 
    ## Run 29 stress 0.1410448 
    ## Run 30 stress 0.1551863 
    ## Run 31 stress 0.2550043 
    ## Run 32 stress 0.1383682 
    ## Run 33 stress 0.1434197 
    ## Run 34 stress 0.2227526 
    ## Run 35 stress 0.1583016 
    ## Run 36 stress 0.1410448 
    ## Run 37 stress 0.141045 
    ## Run 38 stress 0.1415299 
    ## Run 39 stress 0.1693717 
    ## Run 40 stress 0.1830332 
    ## Run 41 stress 0.2384387 
    ## Run 42 stress 0.2074935 
    ## Run 43 stress 0.1998955 
    ## Run 44 stress 0.1383681 
    ## Run 45 stress 0.2008159 
    ## Run 46 stress 0.1320258 
    ## ... Procrustes: rmse 2.579456e-06  max resid 5.314518e-06 
    ## ... Similar to previous best
    ## Run 47 stress 0.3083098 
    ## Run 48 stress 0.1410452 
    ## Run 49 stress 0.2008158 
    ## Run 50 stress 0.1383681 
    ## Run 51 stress 0.1320258 
    ## ... Procrustes: rmse 8.84629e-07  max resid 1.63239e-06 
    ## ... Similar to previous best
    ## Run 52 stress 0.1407209 
    ## Run 53 stress 0.141045 
    ## Run 54 stress 0.2550043 
    ## Run 55 stress 0.1383681 
    ## Run 56 stress 0.1445811 
    ## Run 57 stress 0.2008162 
    ## Run 58 stress 0.1551863 
    ## Run 59 stress 0.1320258 
    ## ... Procrustes: rmse 7.159151e-07  max resid 1.337737e-06 
    ## ... Similar to previous best
    ## Run 60 stress 0.1410448 
    ## Run 61 stress 0.1693717 
    ## Run 62 stress 0.1410449 
    ## Run 63 stress 0.2005461 
    ## Run 64 stress 0.1407298 
    ## Run 65 stress 0.1320258 
    ## ... Procrustes: rmse 1.165022e-06  max resid 2.434336e-06 
    ## ... Similar to previous best
    ## Run 66 stress 0.1320258 
    ## ... Procrustes: rmse 5.80787e-07  max resid 8.55661e-07 
    ## ... Similar to previous best
    ## Run 67 stress 0.1320258 
    ## ... Procrustes: rmse 1.7804e-06  max resid 3.572965e-06 
    ## ... Similar to previous best
    ## Run 68 stress 0.2842805 
    ## Run 69 stress 0.1587623 
    ## Run 70 stress 0.1407298 
    ## Run 71 stress 0.1407298 
    ## Run 72 stress 0.1410448 
    ## Run 73 stress 0.1663544 
    ## Run 74 stress 0.1572668 
    ## Run 75 stress 0.1771966 
    ## Run 76 stress 0.1415299 
    ## Run 77 stress 0.1383682 
    ## Run 78 stress 0.1320258 
    ## ... Procrustes: rmse 1.337571e-06  max resid 2.775041e-06 
    ## ... Similar to previous best
    ## Run 79 stress 0.1434197 
    ## Run 80 stress 0.1970303 
    ## Run 81 stress 0.1551863 
    ## Run 82 stress 0.1693717 
    ## Run 83 stress 0.1407206 
    ## Run 84 stress 0.2520602 
    ## Run 85 stress 0.1771969 
    ## Run 86 stress 0.1320258 
    ## ... Procrustes: rmse 1.115462e-06  max resid 2.253063e-06 
    ## ... Similar to previous best
    ## Run 87 stress 0.1320258 
    ## ... Procrustes: rmse 2.62732e-07  max resid 4.623101e-07 
    ## ... Similar to previous best
    ## Run 88 stress 0.1320258 
    ## ... Procrustes: rmse 1.269083e-06  max resid 2.609839e-06 
    ## ... Similar to previous best
    ## Run 89 stress 0.1320258 
    ## ... Procrustes: rmse 8.61066e-07  max resid 1.653009e-06 
    ## ... Similar to previous best
    ## Run 90 stress 0.1407298 
    ## Run 91 stress 0.1572668 
    ## Run 92 stress 0.2538018 
    ## Run 93 stress 0.1383682 
    ## Run 94 stress 0.1410451 
    ## Run 95 stress 0.1410448 
    ## Run 96 stress 0.1572668 
    ## Run 97 stress 0.1415299 
    ## Run 98 stress 0.1410451 
    ## Run 99 stress 0.1410452 
    ## Run 100 stress 0.1693717 
    ## Run 101 stress 0.1551863 
    ## Run 102 stress 0.1551863 
    ## Run 103 stress 0.1693717 
    ## Run 104 stress 0.2852159 
    ## Run 105 stress 0.1320258 
    ## ... Procrustes: rmse 6.901546e-07  max resid 1.080012e-06 
    ## ... Similar to previous best
    ## Run 106 stress 0.1407298 
    ## Run 107 stress 0.1693717 
    ## Run 108 stress 0.2373121 
    ## Run 109 stress 0.1415299 
    ## Run 110 stress 0.1572668 
    ## Run 111 stress 0.1407298 
    ## Run 112 stress 0.1407298 
    ## Run 113 stress 0.1776793 
    ## Run 114 stress 0.1410453 
    ## Run 115 stress 0.2838148 
    ## Run 116 stress 0.1383681 
    ## Run 117 stress 0.2008158 
    ## Run 118 stress 0.1415299 
    ## Run 119 stress 0.1410449 
    ## Run 120 stress 0.1587623 
    ## Run 121 stress 0.1410451 
    ## Run 122 stress 0.2008166 
    ## Run 123 stress 0.1407207 
    ## Run 124 stress 0.1407298 
    ## Run 125 stress 0.1598154 
    ## Run 126 stress 0.1583016 
    ## Run 127 stress 0.1410448 
    ## Run 128 stress 0.200816 
    ## Run 129 stress 0.1407207 
    ## Run 130 stress 0.1410454 
    ## Run 131 stress 0.1383681 
    ## Run 132 stress 0.1410454 
    ## Run 133 stress 0.1445811 
    ## Run 134 stress 0.1551863 
    ## Run 135 stress 0.2612474 
    ## Run 136 stress 0.1407298 
    ## Run 137 stress 0.1320258 
    ## ... Procrustes: rmse 9.870403e-07  max resid 2.05762e-06 
    ## ... Similar to previous best
    ## Run 138 stress 0.2550043 
    ## Run 139 stress 0.2496383 
    ## Run 140 stress 0.3083098 
    ## Run 141 stress 0.1320258 
    ## ... Procrustes: rmse 3.47911e-07  max resid 4.883861e-07 
    ## ... Similar to previous best
    ## Run 142 stress 0.1572668 
    ## Run 143 stress 0.2416997 
    ## Run 144 stress 0.1407209 
    ## Run 145 stress 0.1383682 
    ## Run 146 stress 0.1383681 
    ## Run 147 stress 0.1415299 
    ## Run 148 stress 0.1771966 
    ## Run 149 stress 0.1410448 
    ## Run 150 stress 0.2224299 
    ## Run 151 stress 0.1407298 
    ## Run 152 stress 0.1320258 
    ## ... Procrustes: rmse 2.973706e-07  max resid 5.838966e-07 
    ## ... Similar to previous best
    ## Run 153 stress 0.1572668 
    ## Run 154 stress 0.1320258 
    ## ... Procrustes: rmse 5.812647e-07  max resid 1.17753e-06 
    ## ... Similar to previous best
    ## Run 155 stress 0.1320258 
    ## ... Procrustes: rmse 1.610524e-06  max resid 3.257827e-06 
    ## ... Similar to previous best
    ## Run 156 stress 0.2227526 
    ## Run 157 stress 0.1771966 
    ## Run 158 stress 0.1383681 
    ## Run 159 stress 0.2118436 
    ## Run 160 stress 0.1320258 
    ## ... Procrustes: rmse 1.845076e-06  max resid 3.570474e-06 
    ## ... Similar to previous best
    ## Run 161 stress 0.1551863 
    ## Run 162 stress 0.1320258 
    ## ... Procrustes: rmse 1.245518e-06  max resid 2.530888e-06 
    ## ... Similar to previous best
    ## Run 163 stress 0.1410452 
    ## Run 164 stress 0.1383681 
    ## Run 165 stress 0.1407298 
    ## Run 166 stress 0.1320258 
    ## ... Procrustes: rmse 2.499473e-06  max resid 5.058791e-06 
    ## ... Similar to previous best
    ## Run 167 stress 0.1551863 
    ## Run 168 stress 0.1320258 
    ## ... Procrustes: rmse 2.910179e-06  max resid 6.040116e-06 
    ## ... Similar to previous best
    ## Run 169 stress 0.1407298 
    ## Run 170 stress 0.1572668 
    ## Run 171 stress 0.1434197 
    ## Run 172 stress 0.1415299 
    ## Run 173 stress 0.2224298 
    ## Run 174 stress 0.3083098 
    ## Run 175 stress 0.1796865 
    ## Run 176 stress 0.1587623 
    ## Run 177 stress 0.1320258 
    ## ... Procrustes: rmse 1.348792e-06  max resid 2.720601e-06 
    ## ... Similar to previous best
    ## Run 178 stress 0.2501066 
    ## Run 179 stress 0.1410451 
    ## Run 180 stress 0.1970303 
    ## Run 181 stress 0.1320258 
    ## ... Procrustes: rmse 3.756452e-07  max resid 7.047354e-07 
    ## ... Similar to previous best
    ## Run 182 stress 0.1830332 
    ## Run 183 stress 0.2015531 
    ## Run 184 stress 0.1434197 
    ## Run 185 stress 0.1320258 
    ## ... Procrustes: rmse 1.130685e-06  max resid 2.106901e-06 
    ## ... Similar to previous best
    ## Run 186 stress 0.1776793 
    ## Run 187 stress 0.1383681 
    ## Run 188 stress 0.1383681 
    ## Run 189 stress 0.1551863 
    ## Run 190 stress 0.2950733 
    ## Run 191 stress 0.1929519 
    ## Run 192 stress 0.1929519 
    ## Run 193 stress 0.2384387 
    ## Run 194 stress 0.1693717 
    ## Run 195 stress 0.1572668 
    ## Run 196 stress 0.1415299 
    ## Run 197 stress 0.1771969 
    ## Run 198 stress 0.1572668 
    ## Run 199 stress 0.2842805 
    ## Run 200 stress 0.2510079 
    ## Run 201 stress 0.1445811 
    ## Run 202 stress 0.1551863 
    ## Run 203 stress 0.1830332 
    ## Run 204 stress 0.1587623 
    ## Run 205 stress 0.1320258 
    ## ... Procrustes: rmse 1.493363e-06  max resid 3.028815e-06 
    ## ... Similar to previous best
    ## Run 206 stress 0.1693717 
    ## Run 207 stress 0.1776793 
    ## Run 208 stress 0.1410452 
    ## Run 209 stress 0.2501066 
    ## Run 210 stress 0.1587623 
    ## Run 211 stress 0.1771969 
    ## Run 212 stress 0.1415299 
    ## Run 213 stress 0.1407298 
    ## Run 214 stress 0.1383681 
    ## Run 215 stress 0.1572668 
    ## Run 216 stress 0.1407298 
    ## Run 217 stress 0.1320258 
    ## ... Procrustes: rmse 7.558944e-07  max resid 9.705691e-07 
    ## ... Similar to previous best
    ## Run 218 stress 0.1551863 
    ## Run 219 stress 0.2385029 
    ## Run 220 stress 0.1415299 
    ## Run 221 stress 0.1551863 
    ## Run 222 stress 0.1583016 
    ## Run 223 stress 0.1320258 
    ## ... Procrustes: rmse 5.92122e-07  max resid 9.511032e-07 
    ## ... Similar to previous best
    ## Run 224 stress 0.1410449 
    ## Run 225 stress 0.1320258 
    ## ... Procrustes: rmse 1.968075e-06  max resid 3.966462e-06 
    ## ... Similar to previous best
    ## Run 226 stress 0.1771969 
    ## Run 227 stress 0.1410451 
    ## Run 228 stress 0.1998955 
    ## Run 229 stress 0.1551863 
    ## Run 230 stress 0.1693719 
    ## Run 231 stress 0.1551863 
    ## Run 232 stress 0.1320258 
    ## ... Procrustes: rmse 7.72295e-07  max resid 1.567399e-06 
    ## ... Similar to previous best
    ## Run 233 stress 0.1410449 
    ## Run 234 stress 0.2085297 
    ## Run 235 stress 0.1628606 
    ## Run 236 stress 0.1407298 
    ## Run 237 stress 0.141045 
    ## Run 238 stress 0.1771966 
    ## Run 239 stress 0.1415299 
    ## Run 240 stress 0.1415299 
    ## Run 241 stress 0.1383682 
    ## Run 242 stress 0.1796864 
    ## Run 243 stress 0.1407209 
    ## Run 244 stress 0.2332576 
    ## Run 245 stress 0.1929519 
    ## Run 246 stress 0.2005461 
    ## Run 247 stress 0.2501063 
    ## Run 248 stress 0.1410454 
    ## Run 249 stress 0.1415299 
    ## Run 250 stress 0.1383681 
    ## Run 251 stress 0.1693717 
    ## Run 252 stress 0.1415299 
    ## Run 253 stress 0.2227527 
    ## Run 254 stress 0.1771969 
    ## Run 255 stress 0.3083099 
    ## Run 256 stress 0.1693717 
    ## Run 257 stress 0.3083098 
    ## Run 258 stress 0.1415299 
    ## Run 259 stress 0.1383682 
    ## Run 260 stress 0.1415299 
    ## Run 261 stress 0.141045 
    ## Run 262 stress 0.2615852 
    ## Run 263 stress 0.1415299 
    ## Run 264 stress 0.1320258 
    ## ... Procrustes: rmse 9.47616e-07  max resid 1.79281e-06 
    ## ... Similar to previous best
    ## Run 265 stress 0.2848514 
    ## Run 266 stress 0.1407298 
    ## Run 267 stress 0.1383681 
    ## Run 268 stress 0.2711712 
    ## Run 269 stress 0.1771969 
    ## Run 270 stress 0.1434197 
    ## Run 271 stress 0.1802752 
    ## Run 272 stress 0.1551863 
    ## Run 273 stress 0.2008158 
    ## Run 274 stress 0.2538018 
    ## Run 275 stress 0.1771969 
    ## Run 276 stress 0.1320258 
    ## ... Procrustes: rmse 2.425453e-06  max resid 5.044081e-06 
    ## ... Similar to previous best
    ## Run 277 stress 0.1407298 
    ## Run 278 stress 0.1434197 
    ## Run 279 stress 0.1628606 
    ## Run 280 stress 0.1407298 
    ## Run 281 stress 0.2015531 
    ## Run 282 stress 0.2456775 
    ## Run 283 stress 0.1320258 
    ## ... Procrustes: rmse 5.551808e-07  max resid 9.402352e-07 
    ## ... Similar to previous best
    ## Run 284 stress 0.1663544 
    ## Run 285 stress 0.1998955 
    ## Run 286 stress 0.1583016 
    ## Run 287 stress 0.1572668 
    ## Run 288 stress 0.1407209 
    ## Run 289 stress 0.1407298 
    ## Run 290 stress 0.1320258 
    ## ... Procrustes: rmse 5.023758e-07  max resid 1.014614e-06 
    ## ... Similar to previous best
    ## Run 291 stress 0.1929519 
    ## Run 292 stress 0.1434197 
    ## Run 293 stress 0.2848515 
    ## Run 294 stress 0.1320258 
    ## ... Procrustes: rmse 2.931403e-07  max resid 4.986293e-07 
    ## ... Similar to previous best
    ## Run 295 stress 0.1320258 
    ## ... Procrustes: rmse 1.407992e-06  max resid 2.89699e-06 
    ## ... Similar to previous best
    ## Run 296 stress 0.2227527 
    ## Run 297 stress 0.1663544 
    ## Run 298 stress 0.1407207 
    ## Run 299 stress 0.1830332 
    ## Run 300 stress 0.1383681 
    ## Run 301 stress 0.1320258 
    ## ... Procrustes: rmse 1.897201e-06  max resid 3.890226e-06 
    ## ... Similar to previous best
    ## Run 302 stress 0.1410454 
    ## Run 303 stress 0.1410448 
    ## Run 304 stress 0.3002107 
    ## Run 305 stress 0.1320258 
    ## ... Procrustes: rmse 1.390537e-06  max resid 2.677467e-06 
    ## ... Similar to previous best
    ## Run 306 stress 0.1383681 
    ## Run 307 stress 0.1693717 
    ## Run 308 stress 0.1776793 
    ## Run 309 stress 0.1693717 
    ## Run 310 stress 0.1572668 
    ## Run 311 stress 0.250851 
    ## Run 312 stress 0.3083099 
    ## Run 313 stress 0.1415299 
    ## Run 314 stress 0.1320258 
    ## ... Procrustes: rmse 1.0814e-06  max resid 2.127894e-06 
    ## ... Similar to previous best
    ## Run 315 stress 0.1320258 
    ## ... Procrustes: rmse 5.881005e-07  max resid 1.035756e-06 
    ## ... Similar to previous best
    ## Run 316 stress 0.1415299 
    ## Run 317 stress 0.1407208 
    ## Run 318 stress 0.1415299 
    ## Run 319 stress 0.1383682 
    ## Run 320 stress 0.1628606 
    ## Run 321 stress 0.1415299 
    ## Run 322 stress 0.1410454 
    ## Run 323 stress 0.1445811 
    ## Run 324 stress 0.1383681 
    ## Run 325 stress 0.1415299 
    ## Run 326 stress 0.1663544 
    ## Run 327 stress 0.1415299 
    ## Run 328 stress 0.1320258 
    ## ... Procrustes: rmse 9.10603e-07  max resid 1.790876e-06 
    ## ... Similar to previous best
    ## Run 329 stress 0.1320258 
    ## ... Procrustes: rmse 1.761767e-06  max resid 3.617065e-06 
    ## ... Similar to previous best
    ## Run 330 stress 0.1415299 
    ## Run 331 stress 0.2008164 
    ## Run 332 stress 0.1410451 
    ## Run 333 stress 0.1830332 
    ## Run 334 stress 0.200816 
    ## Run 335 stress 0.1415299 
    ## Run 336 stress 0.1663544 
    ## Run 337 stress 0.1572668 
    ## Run 338 stress 0.1628606 
    ## Run 339 stress 0.3047186 
    ## Run 340 stress 0.1410448 
    ## Run 341 stress 0.1830332 
    ## Run 342 stress 0.1383681 
    ## Run 343 stress 0.2373121 
    ## Run 344 stress 0.1970303 
    ## Run 345 stress 0.1583016 
    ## Run 346 stress 0.1407207 
    ## Run 347 stress 0.1320258 
    ## ... Procrustes: rmse 8.291094e-07  max resid 1.561111e-06 
    ## ... Similar to previous best
    ## Run 348 stress 0.1383681 
    ## Run 349 stress 0.1320258 
    ## ... Procrustes: rmse 2.909072e-06  max resid 6.025762e-06 
    ## ... Similar to previous best
    ## Run 350 stress 0.2384387 
    ## Run 351 stress 0.1410452 
    ## Run 352 stress 0.1415299 
    ## Run 353 stress 0.1572668 
    ## Run 354 stress 0.1407207 
    ## Run 355 stress 0.1415299 
    ## Run 356 stress 0.1320258 
    ## ... Procrustes: rmse 4.370513e-07  max resid 6.78698e-07 
    ## ... Similar to previous best
    ## Run 357 stress 0.1383682 
    ## Run 358 stress 0.1598154 
    ## Run 359 stress 0.1572668 
    ## Run 360 stress 0.1830332 
    ## Run 361 stress 0.1320258 
    ## ... Procrustes: rmse 3.57153e-06  max resid 7.246137e-06 
    ## ... Similar to previous best
    ## Run 362 stress 0.141045 
    ## Run 363 stress 0.2074935 
    ## Run 364 stress 0.1320258 
    ## ... Procrustes: rmse 1.492266e-06  max resid 3.085343e-06 
    ## ... Similar to previous best
    ## Run 365 stress 0.2416997 
    ## Run 366 stress 0.1970303 
    ## Run 367 stress 0.1970303 
    ## Run 368 stress 0.1320258 
    ## ... Procrustes: rmse 1.614001e-06  max resid 3.131216e-06 
    ## ... Similar to previous best
    ## Run 369 stress 0.141045 
    ## Run 370 stress 0.1407298 
    ## Run 371 stress 0.1415299 
    ## Run 372 stress 0.1929519 
    ## Run 373 stress 0.2514078 
    ## Run 374 stress 0.1415299 
    ## Run 375 stress 0.1407298 
    ## Run 376 stress 0.1415299 
    ## Run 377 stress 0.1320258 
    ## ... Procrustes: rmse 1.019301e-06  max resid 1.977555e-06 
    ## ... Similar to previous best
    ## Run 378 stress 0.1407298 
    ## Run 379 stress 0.1415299 
    ## Run 380 stress 0.1583016 
    ## Run 381 stress 0.1320258 
    ## ... Procrustes: rmse 7.10369e-07  max resid 1.44193e-06 
    ## ... Similar to previous best
    ## Run 382 stress 0.141045 
    ## Run 383 stress 0.1320258 
    ## ... Procrustes: rmse 9.551429e-07  max resid 1.966545e-06 
    ## ... Similar to previous best
    ## Run 384 stress 0.1410448 
    ## Run 385 stress 0.2227526 
    ## Run 386 stress 0.1663544 
    ## Run 387 stress 0.1693717 
    ## Run 388 stress 0.1771966 
    ## Run 389 stress 0.1598154 
    ## Run 390 stress 0.1320258 
    ## ... Procrustes: rmse 9.159215e-07  max resid 1.829767e-06 
    ## ... Similar to previous best
    ## Run 391 stress 0.2962402 
    ## Run 392 stress 0.2008165 
    ## Run 393 stress 0.1415299 
    ## Run 394 stress 0.1551863 
    ## Run 395 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 2.042465e-07  max resid 2.922831e-07 
    ## ... Similar to previous best
    ## Run 396 stress 0.1320258 
    ## ... Procrustes: rmse 1.347865e-06  max resid 2.71665e-06 
    ## ... Similar to previous best
    ## Run 397 stress 0.1410448 
    ## Run 398 stress 0.1572668 
    ## Run 399 stress 0.1410451 
    ## Run 400 stress 0.1320258 
    ## ... Procrustes: rmse 1.093536e-06  max resid 2.297199e-06 
    ## ... Similar to previous best
    ## Run 401 stress 0.1407209 
    ## Run 402 stress 0.1383681 
    ## Run 403 stress 0.1583016 
    ## Run 404 stress 0.1407298 
    ## Run 405 stress 0.1320258 
    ## ... Procrustes: rmse 1.156749e-06  max resid 1.643936e-06 
    ## ... Similar to previous best
    ## Run 406 stress 0.1383682 
    ## Run 407 stress 0.1551863 
    ## Run 408 stress 0.1693717 
    ## Run 409 stress 0.1410452 
    ## Run 410 stress 0.1434197 
    ## Run 411 stress 0.1407298 
    ## Run 412 stress 0.1572668 
    ## Run 413 stress 0.1771966 
    ## Run 414 stress 0.1693717 
    ## Run 415 stress 0.2415237 
    ## Run 416 stress 0.1320258 
    ## ... Procrustes: rmse 4.427927e-07  max resid 8.459903e-07 
    ## ... Similar to previous best
    ## Run 417 stress 0.2852173 
    ## Run 418 stress 0.1572668 
    ## Run 419 stress 0.1320258 
    ## ... Procrustes: rmse 2.290117e-06  max resid 4.717939e-06 
    ## ... Similar to previous best
    ## Run 420 stress 0.141045 
    ## Run 421 stress 0.1415299 
    ## Run 422 stress 0.1383682 
    ## Run 423 stress 0.1320258 
    ## ... Procrustes: rmse 8.968593e-07  max resid 1.476846e-06 
    ## ... Similar to previous best
    ## Run 424 stress 0.1410449 
    ## Run 425 stress 0.1320258 
    ## ... Procrustes: rmse 9.171026e-07  max resid 1.770719e-06 
    ## ... Similar to previous best
    ## Run 426 stress 0.1415299 
    ## Run 427 stress 0.1415299 
    ## Run 428 stress 0.1407298 
    ## Run 429 stress 0.1445811 
    ## Run 430 stress 0.1410449 
    ## Run 431 stress 0.1410449 
    ## Run 432 stress 0.1415299 
    ## Run 433 stress 0.1551863 
    ## Run 434 stress 0.1410451 
    ## Run 435 stress 0.1572668 
    ## Run 436 stress 0.1572668 
    ## Run 437 stress 0.2385029 
    ## Run 438 stress 0.1598154 
    ## Run 439 stress 0.1628606 
    ## Run 440 stress 0.1410451 
    ## Run 441 stress 0.1407298 
    ## Run 442 stress 0.1693717 
    ## Run 443 stress 0.1583016 
    ## Run 444 stress 0.3083098 
    ## Run 445 stress 0.2085297 
    ## Run 446 stress 0.2612821 
    ## Run 447 stress 0.2510079 
    ## Run 448 stress 0.1796865 
    ## Run 449 stress 0.1434197 
    ## Run 450 stress 0.1407298 
    ## Run 451 stress 0.1583016 
    ## Run 452 stress 0.1407298 
    ## Run 453 stress 0.1434197 
    ## Run 454 stress 0.1415299 
    ## Run 455 stress 0.1320258 
    ## ... Procrustes: rmse 1.073698e-06  max resid 2.138459e-06 
    ## ... Similar to previous best
    ## Run 456 stress 0.1320258 
    ## ... Procrustes: rmse 2.081814e-06  max resid 4.301976e-06 
    ## ... Similar to previous best
    ## Run 457 stress 0.1415299 
    ## Run 458 stress 0.1407206 
    ## Run 459 stress 0.2600627 
    ## Run 460 stress 0.1415299 
    ## Run 461 stress 0.1407207 
    ## Run 462 stress 0.1383681 
    ## Run 463 stress 0.1693717 
    ## Run 464 stress 0.1407298 
    ## Run 465 stress 0.1771966 
    ## Run 466 stress 0.1415299 
    ## Run 467 stress 0.1383681 
    ## Run 468 stress 0.1776793 
    ## Run 469 stress 0.1410449 
    ## Run 470 stress 0.1383682 
    ## Run 471 stress 0.1415299 
    ## Run 472 stress 0.1598154 
    ## Run 473 stress 0.2550043 
    ## Run 474 stress 0.1572668 
    ## Run 475 stress 0.1415299 
    ## Run 476 stress 0.1572668 
    ## Run 477 stress 0.1383681 
    ## Run 478 stress 0.1320258 
    ## ... Procrustes: rmse 2.081207e-06  max resid 4.322561e-06 
    ## ... Similar to previous best
    ## Run 479 stress 0.1407208 
    ## Run 480 stress 0.1320258 
    ## ... Procrustes: rmse 7.583322e-07  max resid 1.603087e-06 
    ## ... Similar to previous best
    ## Run 481 stress 0.1551863 
    ## Run 482 stress 0.1572668 
    ## Run 483 stress 0.1320258 
    ## ... Procrustes: rmse 1.110799e-06  max resid 2.174606e-06 
    ## ... Similar to previous best
    ## Run 484 stress 0.1407298 
    ## Run 485 stress 0.1572668 
    ## Run 486 stress 0.1415299 
    ## Run 487 stress 0.1663544 
    ## Run 488 stress 0.1970303 
    ## Run 489 stress 0.1320258 
    ## ... Procrustes: rmse 2.029254e-06  max resid 3.049746e-06 
    ## ... Similar to previous best
    ## Run 490 stress 0.1415299 
    ## Run 491 stress 0.1383682 
    ## Run 492 stress 0.1572668 
    ## Run 493 stress 0.1410454 
    ## Run 494 stress 0.1407298 
    ## Run 495 stress 0.1693717 
    ## Run 496 stress 0.1771966 
    ## Run 497 stress 0.1410449 
    ## Run 498 stress 0.1572668 
    ## Run 499 stress 0.1407207 
    ## Run 500 stress 0.1383681 
    ## *** Best solution repeated 14 times

``` r
### Environmental
# Surveyed sites 
SD_beta_env_NMDS <- metaMDS(SD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07669178 
    ## Run 1 stress 0.08317339 
    ## Run 2 stress 0.08317306 
    ## Run 3 stress 0.08476081 
    ## Run 4 stress 0.08274915 
    ## Run 5 stress 0.08252117 
    ## Run 6 stress 0.08097762 
    ## Run 7 stress 0.08182118 
    ## Run 8 stress 0.08182117 
    ## Run 9 stress 0.07669151 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004027874  max resid 0.0006934588 
    ## ... Similar to previous best
    ## Run 10 stress 0.08097755 
    ## Run 11 stress 0.08355841 
    ## Run 12 stress 0.08175703 
    ## Run 13 stress 0.08927503 
    ## Run 14 stress 0.07731516 
    ## Run 15 stress 0.08249622 
    ## Run 16 stress 0.08104917 
    ## Run 17 stress 0.0801295 
    ## Run 18 stress 0.07669161 
    ## ... Procrustes: rmse 0.0001206556  max resid 0.0002075508 
    ## ... Similar to previous best
    ## Run 19 stress 0.07943324 
    ## Run 20 stress 0.08835156 
    ## Run 21 stress 0.08012974 
    ## Run 22 stress 0.0817571 
    ## Run 23 stress 0.08144922 
    ## Run 24 stress 0.08433704 
    ## Run 25 stress 0.08217341 
    ## Run 26 stress 0.08433712 
    ## Run 27 stress 0.08175701 
    ## Run 28 stress 0.08182124 
    ## Run 29 stress 0.08144937 
    ## Run 30 stress 0.08144861 
    ## Run 31 stress 0.08329767 
    ## Run 32 stress 0.08104896 
    ## Run 33 stress 0.08003197 
    ## Run 34 stress 0.08104924 
    ## Run 35 stress 0.0851415 
    ## Run 36 stress 0.07669155 
    ## ... Procrustes: rmse 0.0001030362  max resid 0.0001932596 
    ## ... Similar to previous best
    ## Run 37 stress 0.084337 
    ## Run 38 stress 0.08104929 
    ## Run 39 stress 0.09022702 
    ## Run 40 stress 0.07951472 
    ## Run 41 stress 0.08217345 
    ## Run 42 stress 0.08182121 
    ## Run 43 stress 0.08252116 
    ## Run 44 stress 0.08252135 
    ## Run 45 stress 0.08097739 
    ## Run 46 stress 0.08144807 
    ## Run 47 stress 0.08335426 
    ## Run 48 stress 0.08421041 
    ## Run 49 stress 0.08252125 
    ## Run 50 stress 0.08241491 
    ## Run 51 stress 0.0800319 
    ## Run 52 stress 0.08684138 
    ## Run 53 stress 0.07669171 
    ## ... Procrustes: rmse 0.0003615467  max resid 0.0006235203 
    ## ... Similar to previous best
    ## Run 54 stress 0.07871761 
    ## Run 55 stress 0.08443304 
    ## Run 56 stress 0.08265407 
    ## Run 57 stress 0.08104907 
    ## Run 58 stress 0.07960658 
    ## Run 59 stress 0.07669151 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001443071  max resid 0.0002511276 
    ## ... Similar to previous best
    ## Run 60 stress 0.08182123 
    ## Run 61 stress 0.08329751 
    ## Run 62 stress 0.08097739 
    ## Run 63 stress 0.07871757 
    ## Run 64 stress 0.08265316 
    ## Run 65 stress 0.08968626 
    ## Run 66 stress 0.08493093 
    ## Run 67 stress 0.08405285 
    ## Run 68 stress 0.08305496 
    ## Run 69 stress 0.0794896 
    ## Run 70 stress 0.08337145 
    ## Run 71 stress 0.07669156 
    ## ... Procrustes: rmse 6.458515e-05  max resid 0.0001222342 
    ## ... Similar to previous best
    ## Run 72 stress 0.0825011 
    ## Run 73 stress 0.07868382 
    ## Run 74 stress 0.08305414 
    ## Run 75 stress 0.0836026 
    ## Run 76 stress 0.07868334 
    ## Run 77 stress 0.07669174 
    ## ... Procrustes: rmse 0.0003574661  max resid 0.0006212382 
    ## ... Similar to previous best
    ## Run 78 stress 0.08283543 
    ## Run 79 stress 0.08175712 
    ## Run 80 stress 0.08329763 
    ## Run 81 stress 0.08360261 
    ## Run 82 stress 0.08252128 
    ## Run 83 stress 0.08898794 
    ## Run 84 stress 0.07868342 
    ## Run 85 stress 0.07871754 
    ## Run 86 stress 0.08305468 
    ## Run 87 stress 0.08957564 
    ## Run 88 stress 0.07871744 
    ## Run 89 stress 0.07871727 
    ## Run 90 stress 0.08337141 
    ## Run 91 stress 0.0795147 
    ## Run 92 stress 0.08657631 
    ## Run 93 stress 0.09053918 
    ## Run 94 stress 0.08175707 
    ## Run 95 stress 0.08182108 
    ## Run 96 stress 0.0848512 
    ## Run 97 stress 0.08584975 
    ## Run 98 stress 0.08013491 
    ## Run 99 stress 0.0821734 
    ## Run 100 stress 0.08553088 
    ## Run 101 stress 0.08218538 
    ## Run 102 stress 0.07960695 
    ## Run 103 stress 0.08217339 
    ## Run 104 stress 0.08433684 
    ## Run 105 stress 0.08684121 
    ## Run 106 stress 0.0826538 
    ## Run 107 stress 0.08250157 
    ## Run 108 stress 0.07949 
    ## Run 109 stress 0.08003186 
    ## Run 110 stress 0.07669161 
    ## ... Procrustes: rmse 0.0002445533  max resid 0.0004911949 
    ## ... Similar to previous best
    ## Run 111 stress 0.08673729 
    ## Run 112 stress 0.07871744 
    ## Run 113 stress 0.07868345 
    ## Run 114 stress 0.08250137 
    ## Run 115 stress 0.08443317 
    ## Run 116 stress 0.08443376 
    ## Run 117 stress 0.0837572 
    ## Run 118 stress 0.07868331 
    ## Run 119 stress 0.08972617 
    ## Run 120 stress 0.08957525 
    ## Run 121 stress 0.08175712 
    ## Run 122 stress 0.08097749 
    ## Run 123 stress 0.08283584 
    ## Run 124 stress 0.0831441 
    ## Run 125 stress 0.08375744 
    ## Run 126 stress 0.08422715 
    ## Run 127 stress 0.08182117 
    ## Run 128 stress 0.08175688 
    ## Run 129 stress 0.08485135 
    ## Run 130 stress 0.08433698 
    ## Run 131 stress 0.08493091 
    ## Run 132 stress 0.08144897 
    ## Run 133 stress 0.08104917 
    ## Run 134 stress 0.08616045 
    ## Run 135 stress 0.07669154 
    ## ... Procrustes: rmse 0.0001802149  max resid 0.0003138054 
    ## ... Similar to previous best
    ## Run 136 stress 0.0786834 
    ## Run 137 stress 0.0849312 
    ## Run 138 stress 0.08362765 
    ## Run 139 stress 0.08240533 
    ## Run 140 stress 0.08305368 
    ## Run 141 stress 0.08684172 
    ## Run 142 stress 0.07731524 
    ## Run 143 stress 0.08305434 
    ## Run 144 stress 0.07871757 
    ## Run 145 stress 0.08182142 
    ## Run 146 stress 0.07731519 
    ## Run 147 stress 0.08175693 
    ## Run 148 stress 0.09196164 
    ## Run 149 stress 0.08012997 
    ## Run 150 stress 0.0897262 
    ## Run 151 stress 0.08175689 
    ## Run 152 stress 0.08556629 
    ## Run 153 stress 0.08097781 
    ## Run 154 stress 0.08104914 
    ## Run 155 stress 0.08547316 
    ## Run 156 stress 0.08422717 
    ## Run 157 stress 0.08553061 
    ## Run 158 stress 0.07871739 
    ## Run 159 stress 0.08362766 
    ## Run 160 stress 0.08175695 
    ## Run 161 stress 0.08421025 
    ## Run 162 stress 0.08003194 
    ## Run 163 stress 0.08375705 
    ## Run 164 stress 0.08420999 
    ## Run 165 stress 0.08252119 
    ## Run 166 stress 0.08443263 
    ## Run 167 stress 0.08317318 
    ## Run 168 stress 0.08305435 
    ## Run 169 stress 0.08189368 
    ## Run 170 stress 0.08265392 
    ## Run 171 stress 0.08317315 
    ## Run 172 stress 0.08360259 
    ## Run 173 stress 0.07871733 
    ## Run 174 stress 0.08421032 
    ## Run 175 stress 0.08250172 
    ## Run 176 stress 0.07669176 
    ## ... Procrustes: rmse 0.0003377546  max resid 0.000578107 
    ## ... Similar to previous best
    ## Run 177 stress 0.08283539 
    ## Run 178 stress 0.08957608 
    ## Run 179 stress 0.08370302 
    ## Run 180 stress 0.08972602 
    ## Run 181 stress 0.08175719 
    ## Run 182 stress 0.08433712 
    ## Run 183 stress 0.08012953 
    ## Run 184 stress 0.07871762 
    ## Run 185 stress 0.08175705 
    ## Run 186 stress 0.07943392 
    ## Run 187 stress 0.08175701 
    ## Run 188 stress 0.08175693 
    ## Run 189 stress 0.08252129 
    ## Run 190 stress 0.07669168 
    ## ... Procrustes: rmse 0.0001939547  max resid 0.0003400559 
    ## ... Similar to previous best
    ## Run 191 stress 0.08499342 
    ## Run 192 stress 0.08337178 
    ## Run 193 stress 0.08252152 
    ## Run 194 stress 0.08684154 
    ## Run 195 stress 0.08175693 
    ## Run 196 stress 0.08433684 
    ## Run 197 stress 0.0817569 
    ## Run 198 stress 0.0870486 
    ## Run 199 stress 0.08097742 
    ## Run 200 stress 0.08217341 
    ## Run 201 stress 0.07868381 
    ## Run 202 stress 0.08471609 
    ## Run 203 stress 0.08175691 
    ## Run 204 stress 0.07669196 
    ## ... Procrustes: rmse 0.0004275134  max resid 0.0007415365 
    ## ... Similar to previous best
    ## Run 205 stress 0.0796066 
    ## Run 206 stress 0.08337214 
    ## Run 207 stress 0.08305476 
    ## Run 208 stress 0.08362758 
    ## Run 209 stress 0.08144833 
    ## Run 210 stress 0.08972611 
    ## Run 211 stress 0.08387167 
    ## Run 212 stress 0.08189375 
    ## Run 213 stress 0.08421004 
    ## Run 214 stress 0.0794337 
    ## Run 215 stress 0.08217344 
    ## Run 216 stress 0.08104887 
    ## Run 217 stress 0.08189368 
    ## Run 218 stress 0.07871734 
    ## Run 219 stress 0.08097759 
    ## Run 220 stress 0.08243402 
    ## Run 221 stress 0.08835149 
    ## Run 222 stress 0.08433689 
    ## Run 223 stress 0.08305472 
    ## Run 224 stress 0.08265452 
    ## Run 225 stress 0.08175708 
    ## Run 226 stress 0.0853018 
    ## Run 227 stress 0.0787174 
    ## Run 228 stress 0.083144 
    ## Run 229 stress 0.08249606 
    ## Run 230 stress 0.0826539 
    ## Run 231 stress 0.08104885 
    ## Run 232 stress 0.08218514 
    ## Run 233 stress 0.08684074 
    ## Run 234 stress 0.08305322 
    ## Run 235 stress 0.08684189 
    ## Run 236 stress 0.08283568 
    ## Run 237 stress 0.08144872 
    ## Run 238 stress 0.08355843 
    ## Run 239 stress 0.07669158 
    ## ... Procrustes: rmse 0.000164122  max resid 0.0003665765 
    ## ... Similar to previous best
    ## Run 240 stress 0.08217346 
    ## Run 241 stress 0.08003205 
    ## Run 242 stress 0.07669158 
    ## ... Procrustes: rmse 8.074443e-05  max resid 0.0001474962 
    ## ... Similar to previous best
    ## Run 243 stress 0.08175704 
    ## Run 244 stress 0.08422684 
    ## Run 245 stress 0.08013006 
    ## Run 246 stress 0.08493154 
    ## Run 247 stress 0.0795149 
    ## Run 248 stress 0.08241487 
    ## Run 249 stress 0.08360276 
    ## Run 250 stress 0.08097768 
    ## Run 251 stress 0.07669164 
    ## ... Procrustes: rmse 0.0002876605  max resid 0.000499452 
    ## ... Similar to previous best
    ## Run 252 stress 0.08175695 
    ## Run 253 stress 0.08317367 
    ## Run 254 stress 0.08501477 
    ## Run 255 stress 0.08097768 
    ## Run 256 stress 0.0766916 
    ## ... Procrustes: rmse 0.0001395007  max resid 0.0002506097 
    ## ... Similar to previous best
    ## Run 257 stress 0.07669175 
    ## ... Procrustes: rmse 0.0002040733  max resid 0.0003511031 
    ## ... Similar to previous best
    ## Run 258 stress 0.08144885 
    ## Run 259 stress 0.08233258 
    ## Run 260 stress 0.08249618 
    ## Run 261 stress 0.0766915 
    ## ... New best solution
    ## ... Procrustes: rmse 6.213926e-05  max resid 0.0001116996 
    ## ... Similar to previous best
    ## Run 262 stress 0.08144912 
    ## Run 263 stress 0.08305492 
    ## Run 264 stress 0.08252142 
    ## Run 265 stress 0.08527016 
    ## Run 266 stress 0.08433703 
    ## Run 267 stress 0.09012168 
    ## Run 268 stress 0.07868366 
    ## Run 269 stress 0.08283569 
    ## Run 270 stress 0.08249607 
    ## Run 271 stress 0.08877279 
    ## Run 272 stress 0.08252131 
    ## Run 273 stress 0.08500332 
    ## Run 274 stress 0.08405274 
    ## Run 275 stress 0.0787177 
    ## Run 276 stress 0.08241494 
    ## Run 277 stress 0.08314384 
    ## Run 278 stress 0.08175691 
    ## Run 279 stress 0.08175708 
    ## Run 280 stress 0.08097754 
    ## Run 281 stress 0.07731517 
    ## Run 282 stress 0.08337178 
    ## Run 283 stress 0.08374713 
    ## Run 284 stress 0.08003188 
    ## Run 285 stress 0.08657522 
    ## Run 286 stress 0.08189395 
    ## Run 287 stress 0.08470879 
    ## Run 288 stress 0.08421013 
    ## Run 289 stress 0.0801296 
    ## Run 290 stress 0.08375707 
    ## Run 291 stress 0.08507049 
    ## Run 292 stress 0.0833716 
    ## Run 293 stress 0.08104886 
    ## Run 294 stress 0.08265453 
    ## Run 295 stress 0.07871737 
    ## Run 296 stress 0.08527028 
    ## Run 297 stress 0.0809774 
    ## Run 298 stress 0.0833715 
    ## Run 299 stress 0.0931669 
    ## Run 300 stress 0.08218519 
    ## Run 301 stress 0.09160893 
    ## Run 302 stress 0.08063928 
    ## Run 303 stress 0.08485164 
    ## Run 304 stress 0.08471601 
    ## Run 305 stress 0.08362754 
    ## Run 306 stress 0.08175698 
    ## Run 307 stress 0.0868408 
    ## Run 308 stress 0.08405264 
    ## Run 309 stress 0.08175699 
    ## Run 310 stress 0.07868365 
    ## Run 311 stress 0.08360264 
    ## Run 312 stress 0.08314371 
    ## Run 313 stress 0.081449 
    ## Run 314 stress 0.08499325 
    ## Run 315 stress 0.08360273 
    ## Run 316 stress 0.08175688 
    ## Run 317 stress 0.08175691 
    ## Run 318 stress 0.08360262 
    ## Run 319 stress 0.07960665 
    ## Run 320 stress 0.08097764 
    ## Run 321 stress 0.08252131 
    ## Run 322 stress 0.08527038 
    ## Run 323 stress 0.08175687 
    ## Run 324 stress 0.08527022 
    ## Run 325 stress 0.08265377 
    ## Run 326 stress 0.08012952 
    ## Run 327 stress 0.07669183 
    ## ... Procrustes: rmse 0.0003501851  max resid 0.0006449663 
    ## ... Similar to previous best
    ## Run 328 stress 0.08355833 
    ## Run 329 stress 0.07669154 
    ## ... Procrustes: rmse 8.670026e-05  max resid 0.0001470198 
    ## ... Similar to previous best
    ## Run 330 stress 0.07960666 
    ## Run 331 stress 0.08175705 
    ## Run 332 stress 0.0817571 
    ## Run 333 stress 0.08265429 
    ## Run 334 stress 0.07669157 
    ## ... Procrustes: rmse 0.0001560247  max resid 0.000265402 
    ## ... Similar to previous best
    ## Run 335 stress 0.08182119 
    ## Run 336 stress 0.08252124 
    ## Run 337 stress 0.07868376 
    ## Run 338 stress 0.08684171 
    ## Run 339 stress 0.08104908 
    ## Run 340 stress 0.08218527 
    ## Run 341 stress 0.08337163 
    ## Run 342 stress 0.08013445 
    ## Run 343 stress 0.08835181 
    ## Run 344 stress 0.0830553 
    ## Run 345 stress 0.07943363 
    ## Run 346 stress 0.08003211 
    ## Run 347 stress 0.08684104 
    ## Run 348 stress 0.07868357 
    ## Run 349 stress 0.08012954 
    ## Run 350 stress 0.08012986 
    ## Run 351 stress 0.08012972 
    ## Run 352 stress 0.08144863 
    ## Run 353 stress 0.07868365 
    ## Run 354 stress 0.07871748 
    ## Run 355 stress 0.07669156 
    ## ... Procrustes: rmse 0.0001462236  max resid 0.0002538633 
    ## ... Similar to previous best
    ## Run 356 stress 0.07871755 
    ## Run 357 stress 0.08314374 
    ## Run 358 stress 0.08493125 
    ## Run 359 stress 0.08835188 
    ## Run 360 stress 0.08012976 
    ## Run 361 stress 0.08144887 
    ## Run 362 stress 0.08683243 
    ## Run 363 stress 0.08189389 
    ## Run 364 stress 0.08104915 
    ## Run 365 stress 0.08175701 
    ## Run 366 stress 0.08477437 
    ## Run 367 stress 0.08968622 
    ## Run 368 stress 0.08317317 
    ## Run 369 stress 0.08339627 
    ## Run 370 stress 0.07669164 
    ## ... Procrustes: rmse 0.0002113418  max resid 0.0003628722 
    ## ... Similar to previous best
    ## Run 371 stress 0.07960672 
    ## Run 372 stress 0.08527043 
    ## Run 373 stress 0.08355806 
    ## Run 374 stress 0.08182119 
    ## Run 375 stress 0.08527026 
    ## Run 376 stress 0.08835132 
    ## Run 377 stress 0.08428748 
    ## Run 378 stress 0.0766915 
    ## ... Procrustes: rmse 5.640349e-05  max resid 0.0001164193 
    ## ... Similar to previous best
    ## Run 379 stress 0.08682638 
    ## Run 380 stress 0.08012978 
    ## Run 381 stress 0.08097763 
    ## Run 382 stress 0.08233242 
    ## Run 383 stress 0.0909274 
    ## Run 384 stress 0.08265376 
    ## Run 385 stress 0.08305055 
    ## Run 386 stress 0.08265365 
    ## Run 387 stress 0.08012943 
    ## Run 388 stress 0.07871757 
    ## Run 389 stress 0.08175692 
    ## Run 390 stress 0.08249611 
    ## Run 391 stress 0.08548277 
    ## Run 392 stress 0.08418575 
    ## Run 393 stress 0.07868346 
    ## Run 394 stress 0.08003223 
    ## Run 395 stress 0.08182119 
    ## Run 396 stress 0.07731526 
    ## Run 397 stress 0.08249604 
    ## Run 398 stress 0.09022508 
    ## Run 399 stress 0.08337137 
    ## Run 400 stress 0.07669166 
    ## ... Procrustes: rmse 0.0002374375  max resid 0.0004266415 
    ## ... Similar to previous best
    ## Run 401 stress 0.08470942 
    ## Run 402 stress 0.08370318 
    ## Run 403 stress 0.08704837 
    ## Run 404 stress 0.08144957 
    ## Run 405 stress 0.08265475 
    ## Run 406 stress 0.08063951 
    ## Run 407 stress 0.08175702 
    ## Run 408 stress 0.07868359 
    ## Run 409 stress 0.08514224 
    ## Run 410 stress 0.07669153 
    ## ... Procrustes: rmse 0.0001175802  max resid 0.0002654371 
    ## ... Similar to previous best
    ## Run 411 stress 0.07868336 
    ## Run 412 stress 0.08003195 
    ## Run 413 stress 0.08464819 
    ## Run 414 stress 0.08233211 
    ## Run 415 stress 0.08013717 
    ## Run 416 stress 0.08144884 
    ## Run 417 stress 0.08274927 
    ## Run 418 stress 0.08421037 
    ## Run 419 stress 0.08250149 
    ## Run 420 stress 0.08433684 
    ## Run 421 stress 0.07951478 
    ## Run 422 stress 0.08241491 
    ## Run 423 stress 0.08182111 
    ## Run 424 stress 0.08252123 
    ## Run 425 stress 0.08175705 
    ## Run 426 stress 0.08317307 
    ## Run 427 stress 0.0830532 
    ## Run 428 stress 0.08182108 
    ## Run 429 stress 0.08265433 
    ## Run 430 stress 0.08530176 
    ## Run 431 stress 0.08182131 
    ## Run 432 stress 0.08453229 
    ## Run 433 stress 0.0801295 
    ## Run 434 stress 0.08175693 
    ## Run 435 stress 0.08097759 
    ## Run 436 stress 0.08243367 
    ## Run 437 stress 0.08337175 
    ## Run 438 stress 0.08673716 
    ## Run 439 stress 0.08250169 
    ## Run 440 stress 0.08175722 
    ## Run 441 stress 0.07669167 
    ## ... Procrustes: rmse 0.0002432831  max resid 0.0004231191 
    ## ... Similar to previous best
    ## Run 442 stress 0.08175691 
    ## Run 443 stress 0.08003182 
    ## Run 444 stress 0.08242187 
    ## Run 445 stress 0.09196159 
    ## Run 446 stress 0.08527026 
    ## Run 447 stress 0.08422734 
    ## Run 448 stress 0.08249595 
    ## Run 449 stress 0.0787173 
    ## Run 450 stress 0.08339811 
    ## Run 451 stress 0.07669179 
    ## ... Procrustes: rmse 0.0002887291  max resid 0.0004999453 
    ## ... Similar to previous best
    ## Run 452 stress 0.08433699 
    ## Run 453 stress 0.08337178 
    ## Run 454 stress 0.08175703 
    ## Run 455 stress 0.0818212 
    ## Run 456 stress 0.08252126 
    ## Run 457 stress 0.08012946 
    ## Run 458 stress 0.08304899 
    ## Run 459 stress 0.08104891 
    ## Run 460 stress 0.08421028 
    ## Run 461 stress 0.08175693 
    ## Run 462 stress 0.08527038 
    ## Run 463 stress 0.08175696 
    ## Run 464 stress 0.08265371 
    ## Run 465 stress 0.07669171 
    ## ... Procrustes: rmse 0.0002764547  max resid 0.0004801277 
    ## ... Similar to previous best
    ## Run 466 stress 0.08443195 
    ## Run 467 stress 0.08003202 
    ## Run 468 stress 0.07868355 
    ## Run 469 stress 0.08835136 
    ## Run 470 stress 0.08104891 
    ## Run 471 stress 0.08003193 
    ## Run 472 stress 0.0766917 
    ## ... Procrustes: rmse 0.0002771543  max resid 0.0005225684 
    ## ... Similar to previous best
    ## Run 473 stress 0.0850037 
    ## Run 474 stress 0.08252127 
    ## Run 475 stress 0.08360268 
    ## Run 476 stress 0.07949011 
    ## Run 477 stress 0.08470925 
    ## Run 478 stress 0.08233244 
    ## Run 479 stress 0.08485157 
    ## Run 480 stress 0.08175709 
    ## Run 481 stress 0.08097748 
    ## Run 482 stress 0.0801297 
    ## Run 483 stress 0.07871731 
    ## Run 484 stress 0.08217341 
    ## Run 485 stress 0.08527047 
    ## Run 486 stress 0.08175699 
    ## Run 487 stress 0.08217375 
    ## Run 488 stress 0.08499348 
    ## Run 489 stress 0.08003197 
    ## Run 490 stress 0.08104889 
    ## Run 491 stress 0.0817569 
    ## Run 492 stress 0.07669155 
    ## ... Procrustes: rmse 0.0001332612  max resid 0.0002313726 
    ## ... Similar to previous best
    ## Run 493 stress 0.08274935 
    ## Run 494 stress 0.07871726 
    ## Run 495 stress 0.08470908 
    ## Run 496 stress 0.08314392 
    ## Run 497 stress 0.08568814 
    ## Run 498 stress 0.08485138 
    ## Run 499 stress 0.07943352 
    ## Run 500 stress 0.08584081 
    ## *** Best solution repeated 14 times

``` r
# Mixed and stratified lakes
SD_beta_env_MS_NMDS <- metaMDS(SD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09145329 
    ## Run 2 stress 0.09030404 
    ## Run 3 stress 0.09465905 
    ## Run 4 stress 0.08973885 
    ## Run 5 stress 0.2936875 
    ## Run 6 stress 0.09407977 
    ## Run 7 stress 0.09416149 
    ## Run 8 stress 0.09374285 
    ## Run 9 stress 0.08440279 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01443912  max resid 0.04316214 
    ## Run 10 stress 0.09159091 
    ## Run 11 stress 0.09374232 
    ## Run 12 stress 0.090304 
    ## Run 13 stress 0.09030393 
    ## Run 14 stress 0.09381846 
    ## Run 15 stress 0.08773497 
    ## Run 16 stress 0.09168943 
    ## Run 17 stress 0.08440273 
    ## ... New best solution
    ## ... Procrustes: rmse 3.76487e-05  max resid 9.190879e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.09465906 
    ## Run 19 stress 0.09321886 
    ## Run 20 stress 0.08503487 
    ## Run 21 stress 0.09374218 
    ## Run 22 stress 0.0933729 
    ## Run 23 stress 0.09403428 
    ## Run 24 stress 0.0915909 
    ## Run 25 stress 0.08503521 
    ## Run 26 stress 0.09969967 
    ## Run 27 stress 0.09407991 
    ## Run 28 stress 0.08503528 
    ## Run 29 stress 0.09760817 
    ## Run 30 stress 0.0850363 
    ## Run 31 stress 0.09145334 
    ## Run 32 stress 0.2612629 
    ## Run 33 stress 0.08503478 
    ## Run 34 stress 0.08440262 
    ## ... New best solution
    ## ... Procrustes: rmse 9.989827e-05  max resid 0.000196528 
    ## ... Similar to previous best
    ## Run 35 stress 0.090304 
    ## Run 36 stress 0.09168945 
    ## Run 37 stress 0.09030399 
    ## Run 38 stress 0.08773484 
    ## Run 39 stress 0.09380599 
    ## Run 40 stress 0.09337235 
    ## Run 41 stress 0.09374191 
    ## Run 42 stress 0.1038048 
    ## Run 43 stress 0.09145321 
    ## Run 44 stress 0.09374238 
    ## Run 45 stress 0.09416137 
    ## Run 46 stress 0.08440254 
    ## ... New best solution
    ## ... Procrustes: rmse 8.555776e-05  max resid 0.0001748874 
    ## ... Similar to previous best
    ## Run 47 stress 0.08440265 
    ## ... Procrustes: rmse 0.000112795  max resid 0.000218171 
    ## ... Similar to previous best
    ## Run 48 stress 0.08503641 
    ## Run 49 stress 0.09590623 
    ## Run 50 stress 0.09416134 
    ## Run 51 stress 0.09464449 
    ## Run 52 stress 0.09416141 
    ## Run 53 stress 0.09712953 
    ## Run 54 stress 0.09268337 
    ## Run 55 stress 0.0937424 
    ## Run 56 stress 0.08973877 
    ## Run 57 stress 0.09308946 
    ## Run 58 stress 0.08503479 
    ## Run 59 stress 0.09539197 
    ## Run 60 stress 0.0926836 
    ## Run 61 stress 0.08440256 
    ## ... Procrustes: rmse 5.269422e-05  max resid 9.127662e-05 
    ## ... Similar to previous best
    ## Run 62 stress 0.09268368 
    ## Run 63 stress 0.08503634 
    ## Run 64 stress 0.09407345 
    ## Run 65 stress 0.09400527 
    ## Run 66 stress 0.09381865 
    ## Run 67 stress 0.09308947 
    ## Run 68 stress 0.1026216 
    ## Run 69 stress 0.08773485 
    ## Run 70 stress 0.08440261 
    ## ... Procrustes: rmse 8.847985e-05  max resid 0.000157327 
    ## ... Similar to previous best
    ## Run 71 stress 0.09407128 
    ## Run 72 stress 0.08503516 
    ## Run 73 stress 0.09268365 
    ## Run 74 stress 0.09030397 
    ## Run 75 stress 0.08503642 
    ## Run 76 stress 0.09030417 
    ## Run 77 stress 0.103805 
    ## Run 78 stress 0.08773474 
    ## Run 79 stress 0.09286084 
    ## Run 80 stress 0.1038083 
    ## Run 81 stress 0.09535443 
    ## Run 82 stress 0.09145324 
    ## Run 83 stress 0.09030398 
    ## Run 84 stress 0.09374281 
    ## Run 85 stress 0.09030402 
    ## Run 86 stress 0.09381844 
    ## Run 87 stress 0.09030396 
    ## Run 88 stress 0.09465907 
    ## Run 89 stress 0.09308968 
    ## Run 90 stress 0.08773468 
    ## Run 91 stress 0.1052849 
    ## Run 92 stress 0.09030404 
    ## Run 93 stress 0.09286084 
    ## Run 94 stress 0.09400519 
    ## Run 95 stress 0.08973867 
    ## Run 96 stress 0.09286085 
    ## Run 97 stress 0.09268347 
    ## Run 98 stress 0.09145321 
    ## Run 99 stress 0.08503501 
    ## Run 100 stress 0.09286084 
    ## Run 101 stress 0.09465905 
    ## Run 102 stress 0.2461312 
    ## Run 103 stress 0.09408004 
    ## Run 104 stress 0.09268365 
    ## Run 105 stress 0.08503533 
    ## Run 106 stress 0.09539199 
    ## Run 107 stress 0.0850357 
    ## Run 108 stress 0.09145322 
    ## Run 109 stress 0.08440264 
    ## ... Procrustes: rmse 0.000105501  max resid 0.0001892881 
    ## ... Similar to previous best
    ## Run 110 stress 0.105285 
    ## Run 111 stress 0.09168926 
    ## Run 112 stress 0.09337259 
    ## Run 113 stress 0.09337276 
    ## Run 114 stress 0.08503469 
    ## Run 115 stress 0.08973872 
    ## Run 116 stress 0.09337231 
    ## Run 117 stress 0.08773468 
    ## Run 118 stress 0.09337278 
    ## Run 119 stress 0.09465914 
    ## Run 120 stress 0.09374223 
    ## Run 121 stress 0.09969972 
    ## Run 122 stress 0.09268397 
    ## Run 123 stress 0.09590895 
    ## Run 124 stress 0.09030396 
    ## Run 125 stress 0.103804 
    ## Run 126 stress 0.09168927 
    ## Run 127 stress 0.09168926 
    ## Run 128 stress 0.09030393 
    ## Run 129 stress 0.08503476 
    ## Run 130 stress 0.09308957 
    ## Run 131 stress 0.08503552 
    ## Run 132 stress 0.09321542 
    ## Run 133 stress 0.0916894 
    ## Run 134 stress 0.08503601 
    ## Run 135 stress 0.09407104 
    ## Run 136 stress 0.08973867 
    ## Run 137 stress 0.3411762 
    ## Run 138 stress 0.09145321 
    ## Run 139 stress 0.0844026 
    ## ... Procrustes: rmse 7.749024e-05  max resid 0.0001398993 
    ## ... Similar to previous best
    ## Run 140 stress 0.08503465 
    ## Run 141 stress 0.0926836 
    ## Run 142 stress 0.08773473 
    ## Run 143 stress 0.09168954 
    ## Run 144 stress 0.09030394 
    ## Run 145 stress 0.1005996 
    ## Run 146 stress 0.1052849 
    ## Run 147 stress 0.09721197 
    ## Run 148 stress 0.08773466 
    ## Run 149 stress 0.10196 
    ## Run 150 stress 0.09145325 
    ## Run 151 stress 0.09760831 
    ## Run 152 stress 0.09407103 
    ## Run 153 stress 0.09308971 
    ## Run 154 stress 0.09416138 
    ## Run 155 stress 0.09268323 
    ## Run 156 stress 0.09407983 
    ## Run 157 stress 0.08503474 
    ## Run 158 stress 0.09969994 
    ## Run 159 stress 0.08773482 
    ## Run 160 stress 0.09145332 
    ## Run 161 stress 0.09159085 
    ## Run 162 stress 0.08773468 
    ## Run 163 stress 0.09590519 
    ## Run 164 stress 0.09030418 
    ## Run 165 stress 0.0916894 
    ## Run 166 stress 0.09308965 
    ## Run 167 stress 0.09407316 
    ## Run 168 stress 0.09969985 
    ## Run 169 stress 0.09168971 
    ## Run 170 stress 0.08773479 
    ## Run 171 stress 0.09337284 
    ## Run 172 stress 0.08503544 
    ## Run 173 stress 0.08440261 
    ## ... Procrustes: rmse 8.856374e-05  max resid 0.0001720174 
    ## ... Similar to previous best
    ## Run 174 stress 0.09381852 
    ## Run 175 stress 0.09030398 
    ## Run 176 stress 0.09407993 
    ## Run 177 stress 0.1038081 
    ## Run 178 stress 0.08440267 
    ## ... Procrustes: rmse 0.0001361463  max resid 0.0002639477 
    ## ... Similar to previous best
    ## Run 179 stress 0.09268387 
    ## Run 180 stress 0.0877347 
    ## Run 181 stress 0.09407992 
    ## Run 182 stress 0.09464465 
    ## Run 183 stress 0.09308982 
    ## Run 184 stress 0.08503471 
    ## Run 185 stress 0.09380576 
    ## Run 186 stress 0.09760822 
    ## Run 187 stress 0.08773487 
    ## Run 188 stress 0.09145335 
    ## Run 189 stress 0.09030412 
    ## Run 190 stress 0.0850346 
    ## Run 191 stress 0.09268326 
    ## Run 192 stress 0.09268352 
    ## Run 193 stress 0.09374288 
    ## Run 194 stress 0.09145325 
    ## Run 195 stress 0.09030401 
    ## Run 196 stress 0.09535426 
    ## Run 197 stress 0.08973866 
    ## Run 198 stress 0.08503489 
    ## Run 199 stress 0.08773467 
    ## Run 200 stress 0.09145321 
    ## Run 201 stress 0.08503486 
    ## Run 202 stress 0.09030396 
    ## Run 203 stress 0.09407122 
    ## Run 204 stress 0.08503551 
    ## Run 205 stress 0.09445489 
    ## Run 206 stress 0.09159086 
    ## Run 207 stress 0.09168946 
    ## Run 208 stress 0.09159087 
    ## Run 209 stress 0.09535429 
    ## Run 210 stress 0.09374334 
    ## Run 211 stress 0.1004612 
    ## Run 212 stress 0.08773466 
    ## Run 213 stress 0.09969982 
    ## Run 214 stress 0.09145327 
    ## Run 215 stress 0.09416157 
    ## Run 216 stress 0.08503474 
    ## Run 217 stress 0.09465916 
    ## Run 218 stress 0.08773468 
    ## Run 219 stress 0.09417787 
    ## Run 220 stress 0.09969985 
    ## Run 221 stress 0.09374271 
    ## Run 222 stress 0.09464507 
    ## Run 223 stress 0.0930896 
    ## Run 224 stress 0.09337258 
    ## Run 225 stress 0.09416136 
    ## Run 226 stress 0.09145328 
    ## Run 227 stress 0.08973875 
    ## Run 228 stress 0.08973881 
    ## Run 229 stress 0.09337287 
    ## Run 230 stress 0.09590932 
    ## Run 231 stress 0.08503606 
    ## Run 232 stress 0.09268372 
    ## Run 233 stress 0.1038045 
    ## Run 234 stress 0.08503476 
    ## Run 235 stress 0.09308968 
    ## Run 236 stress 0.2688247 
    ## Run 237 stress 0.09408002 
    ## Run 238 stress 0.0926836 
    ## Run 239 stress 0.08503484 
    ## Run 240 stress 0.09539198 
    ## Run 241 stress 0.1095796 
    ## Run 242 stress 0.09308994 
    ## Run 243 stress 0.09159089 
    ## Run 244 stress 0.1026216 
    ## Run 245 stress 0.08503584 
    ## Run 246 stress 0.09145321 
    ## Run 247 stress 0.08773467 
    ## Run 248 stress 0.09308953 
    ## Run 249 stress 0.09286088 
    ## Run 250 stress 0.08440258 
    ## ... Procrustes: rmse 0.0001338428  max resid 0.0004086489 
    ## ... Similar to previous best
    ## Run 251 stress 0.09337293 
    ## Run 252 stress 0.09168941 
    ## Run 253 stress 0.09268334 
    ## Run 254 stress 0.2492117 
    ## Run 255 stress 0.09308961 
    ## Run 256 stress 0.09030394 
    ## Run 257 stress 0.09445467 
    ## Run 258 stress 0.09286098 
    ## Run 259 stress 0.09030402 
    ## Run 260 stress 0.09337229 
    ## Run 261 stress 0.09407997 
    ## Run 262 stress 0.09445479 
    ## Run 263 stress 0.08503581 
    ## Run 264 stress 0.09445501 
    ## Run 265 stress 0.09268371 
    ## Run 266 stress 0.09380576 
    ## Run 267 stress 0.09030399 
    ## Run 268 stress 0.09465915 
    ## Run 269 stress 0.08440265 
    ## ... Procrustes: rmse 0.0001281886  max resid 0.0002177531 
    ## ... Similar to previous best
    ## Run 270 stress 0.09168936 
    ## Run 271 stress 0.08973868 
    ## Run 272 stress 0.08440257 
    ## ... Procrustes: rmse 0.0002492867  max resid 0.0004822009 
    ## ... Similar to previous best
    ## Run 273 stress 0.09145323 
    ## Run 274 stress 0.0940797 
    ## Run 275 stress 0.09030404 
    ## Run 276 stress 0.08440259 
    ## ... Procrustes: rmse 6.76795e-05  max resid 0.0001230113 
    ## ... Similar to previous best
    ## Run 277 stress 0.09030393 
    ## Run 278 stress 0.08503499 
    ## Run 279 stress 0.08503529 
    ## Run 280 stress 0.09308952 
    ## Run 281 stress 0.0941616 
    ## Run 282 stress 0.09539193 
    ## Run 283 stress 0.09400529 
    ## Run 284 stress 0.09416132 
    ## Run 285 stress 0.09969982 
    ## Run 286 stress 0.09308946 
    ## Run 287 stress 0.08503699 
    ## Run 288 stress 0.08440253 
    ## ... New best solution
    ## ... Procrustes: rmse 6.178377e-05  max resid 0.0001227662 
    ## ... Similar to previous best
    ## Run 289 stress 0.09464453 
    ## Run 290 stress 0.09380586 
    ## Run 291 stress 0.09416132 
    ## Run 292 stress 0.09168946 
    ## Run 293 stress 0.2574362 
    ## Run 294 stress 0.103804 
    ## Run 295 stress 0.09308947 
    ## Run 296 stress 0.09407972 
    ## Run 297 stress 0.09337266 
    ## Run 298 stress 0.2643459 
    ## Run 299 stress 0.09030401 
    ## Run 300 stress 0.08973875 
    ## Run 301 stress 0.08773467 
    ## Run 302 stress 0.09416138 
    ## Run 303 stress 0.09159083 
    ## Run 304 stress 0.08773494 
    ## Run 305 stress 0.09465904 
    ## Run 306 stress 0.09337285 
    ## Run 307 stress 0.09030401 
    ## Run 308 stress 0.08773465 
    ## Run 309 stress 0.09400525 
    ## Run 310 stress 0.08773468 
    ## Run 311 stress 0.09286085 
    ## Run 312 stress 0.09286083 
    ## Run 313 stress 0.1052851 
    ## Run 314 stress 0.0844028 
    ## ... Procrustes: rmse 0.0002712752  max resid 0.0004859289 
    ## ... Similar to previous best
    ## Run 315 stress 0.1026215 
    ## Run 316 stress 0.0926833 
    ## Run 317 stress 0.08503488 
    ## Run 318 stress 0.1038043 
    ## Run 319 stress 0.09145324 
    ## Run 320 stress 0.09407971 
    ## Run 321 stress 0.09969964 
    ## Run 322 stress 0.08773471 
    ## Run 323 stress 0.09590905 
    ## Run 324 stress 0.09416135 
    ## Run 325 stress 0.09590507 
    ## Run 326 stress 0.09268346 
    ## Run 327 stress 0.1038047 
    ## Run 328 stress 0.09447309 
    ## Run 329 stress 0.08973863 
    ## Run 330 stress 0.08503468 
    ## Run 331 stress 0.09407998 
    ## Run 332 stress 0.08503499 
    ## Run 333 stress 0.09403406 
    ## Run 334 stress 0.09145323 
    ## Run 335 stress 0.09374274 
    ## Run 336 stress 0.08973871 
    ## Run 337 stress 0.09030405 
    ## Run 338 stress 0.08440255 
    ## ... Procrustes: rmse 7.506251e-05  max resid 0.0001621803 
    ## ... Similar to previous best
    ## Run 339 stress 0.09407979 
    ## Run 340 stress 0.09465904 
    ## Run 341 stress 0.09381849 
    ## Run 342 stress 0.09416132 
    ## Run 343 stress 0.0940797 
    ## Run 344 stress 0.09286085 
    ## Run 345 stress 0.09030403 
    ## Run 346 stress 0.09416134 
    ## Run 347 stress 0.09030397 
    ## Run 348 stress 0.09030415 
    ## Run 349 stress 0.1038035 
    ## Run 350 stress 0.09464447 
    ## Run 351 stress 0.09416148 
    ## Run 352 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 5.332323e-05  max resid 0.0001045925 
    ## ... Similar to previous best
    ## Run 353 stress 0.09445481 
    ## Run 354 stress 0.09030408 
    ## Run 355 stress 0.09539205 
    ## Run 356 stress 0.09145324 
    ## Run 357 stress 0.08773467 
    ## Run 358 stress 0.1052852 
    ## Run 359 stress 0.248884 
    ## Run 360 stress 0.09168948 
    ## Run 361 stress 0.2914661 
    ## Run 362 stress 0.08440258 
    ## ... Procrustes: rmse 0.0001615498  max resid 0.0003401048 
    ## ... Similar to previous best
    ## Run 363 stress 0.09464293 
    ## Run 364 stress 0.08440254 
    ## ... Procrustes: rmse 0.0001126811  max resid 0.0002182841 
    ## ... Similar to previous best
    ## Run 365 stress 0.09407969 
    ## Run 366 stress 0.09590942 
    ## Run 367 stress 0.09321574 
    ## Run 368 stress 0.09669829 
    ## Run 369 stress 0.09030406 
    ## Run 370 stress 0.08973863 
    ## Run 371 stress 0.08503598 
    ## Run 372 stress 0.0877348 
    ## Run 373 stress 0.09407981 
    ## Run 374 stress 0.09168933 
    ## Run 375 stress 0.08773471 
    ## Run 376 stress 0.09268369 
    ## Run 377 stress 0.08973872 
    ## Run 378 stress 0.09159095 
    ## Run 379 stress 0.09268382 
    ## Run 380 stress 0.08973866 
    ## Run 381 stress 0.0850348 
    ## Run 382 stress 0.1054536 
    ## Run 383 stress 0.09159087 
    ## Run 384 stress 0.09168936 
    ## Run 385 stress 0.08773476 
    ## Run 386 stress 0.09030393 
    ## Run 387 stress 0.1038039 
    ## Run 388 stress 0.09145321 
    ## Run 389 stress 0.08773469 
    ## Run 390 stress 0.09321796 
    ## Run 391 stress 0.09760815 
    ## Run 392 stress 0.09969985 
    ## Run 393 stress 0.090304 
    ## Run 394 stress 0.08973862 
    ## Run 395 stress 0.09374248 
    ## Run 396 stress 0.09268324 
    ## Run 397 stress 0.09969973 
    ## Run 398 stress 0.09145322 
    ## Run 399 stress 0.08503468 
    ## Run 400 stress 0.09464505 
    ## Run 401 stress 0.09308964 
    ## Run 402 stress 0.09374314 
    ## Run 403 stress 0.09145348 
    ## Run 404 stress 0.09337277 
    ## Run 405 stress 0.09030412 
    ## Run 406 stress 0.09380585 
    ## Run 407 stress 0.09590936 
    ## Run 408 stress 0.08973877 
    ## Run 409 stress 0.0850347 
    ## Run 410 stress 0.09337267 
    ## Run 411 stress 0.08503465 
    ## Run 412 stress 0.08973862 
    ## Run 413 stress 0.08440262 
    ## ... Procrustes: rmse 0.0001953667  max resid 0.0003642871 
    ## ... Similar to previous best
    ## Run 414 stress 0.08503565 
    ## Run 415 stress 0.09168953 
    ## Run 416 stress 0.08973871 
    ## Run 417 stress 0.08503488 
    ## Run 418 stress 0.09407972 
    ## Run 419 stress 0.10958 
    ## Run 420 stress 0.09030396 
    ## Run 421 stress 0.090304 
    ## Run 422 stress 0.09721198 
    ## Run 423 stress 0.08503593 
    ## Run 424 stress 0.08503571 
    ## Run 425 stress 0.3164639 
    ## Run 426 stress 0.09030406 
    ## Run 427 stress 0.09030408 
    ## Run 428 stress 0.09374282 
    ## Run 429 stress 0.08503486 
    ## Run 430 stress 0.08773491 
    ## Run 431 stress 0.08440275 
    ## ... Procrustes: rmse 0.0002886511  max resid 0.0005741943 
    ## ... Similar to previous best
    ## Run 432 stress 0.09308946 
    ## Run 433 stress 0.09145322 
    ## Run 434 stress 0.09381844 
    ## Run 435 stress 0.08503534 
    ## Run 436 stress 0.09407972 
    ## Run 437 stress 0.09268369 
    ## Run 438 stress 0.09969975 
    ## Run 439 stress 0.09969964 
    ## Run 440 stress 0.09374305 
    ## Run 441 stress 0.08503472 
    ## Run 442 stress 0.09168929 
    ## Run 443 stress 0.111379 
    ## Run 444 stress 0.09721199 
    ## Run 445 stress 0.08503493 
    ## Run 446 stress 0.09268376 
    ## Run 447 stress 0.09168962 
    ## Run 448 stress 0.09969982 
    ## Run 449 stress 0.1052848 
    ## Run 450 stress 0.09465905 
    ## Run 451 stress 0.09308964 
    ## Run 452 stress 0.09145324 
    ## Run 453 stress 0.08773468 
    ## Run 454 stress 0.08973881 
    ## Run 455 stress 0.3003119 
    ## Run 456 stress 0.09030396 
    ## Run 457 stress 0.09465905 
    ## Run 458 stress 0.09374245 
    ## Run 459 stress 0.09535501 
    ## Run 460 stress 0.09168933 
    ## Run 461 stress 0.08503478 
    ## Run 462 stress 0.08773469 
    ## Run 463 stress 0.08973867 
    ## Run 464 stress 0.09408013 
    ## Run 465 stress 0.09380588 
    ## Run 466 stress 0.09030401 
    ## Run 467 stress 0.08503481 
    ## Run 468 stress 0.0941614 
    ## Run 469 stress 0.09268414 
    ## Run 470 stress 0.09445496 
    ## Run 471 stress 0.08973888 
    ## Run 472 stress 0.09030399 
    ## Run 473 stress 0.09286097 
    ## Run 474 stress 0.08503486 
    ## Run 475 stress 0.09969966 
    ## Run 476 stress 0.08973883 
    ## Run 477 stress 0.09374303 
    ## Run 478 stress 0.08440255 
    ## ... Procrustes: rmse 8.596412e-05  max resid 0.000168639 
    ## ... Similar to previous best
    ## Run 479 stress 0.08503507 
    ## Run 480 stress 0.08503482 
    ## Run 481 stress 0.09760831 
    ## Run 482 stress 0.09407966 
    ## Run 483 stress 0.09404356 
    ## Run 484 stress 0.08503508 
    ## Run 485 stress 0.09407984 
    ## Run 486 stress 0.08773465 
    ## Run 487 stress 0.0915909 
    ## Run 488 stress 0.09416131 
    ## Run 489 stress 0.09159085 
    ## Run 490 stress 0.09268403 
    ## Run 491 stress 0.0926835 
    ## Run 492 stress 0.08440251 
    ## ... Procrustes: rmse 3.707619e-05  max resid 6.915048e-05 
    ## ... Similar to previous best
    ## Run 493 stress 0.09337261 
    ## Run 494 stress 0.08773481 
    ## Run 495 stress 0.09286084 
    ## Run 496 stress 0.09308975 
    ## Run 497 stress 0.09030398 
    ## Run 498 stress 0.09159098 
    ## Run 499 stress 0.08973868 
    ## Run 500 stress 0.09610776 
    ## *** Best solution repeated 7 times

``` r
# Ocean sites and mixed lakes
SD_beta_env_OM_NMDS <- metaMDS(SD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.06623458 
    ## Run 1 stress 0.07411945 
    ## Run 2 stress 0.06477902 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1532118  max resid 0.265361 
    ## Run 3 stress 0.06623459 
    ## Run 4 stress 0.07166973 
    ## Run 5 stress 0.06623462 
    ## Run 6 stress 0.06623459 
    ## Run 7 stress 0.07411941 
    ## Run 8 stress 0.0612452 
    ## ... New best solution
    ## ... Procrustes: rmse 0.09078582  max resid 0.2170445 
    ## Run 9 stress 0.06623457 
    ## Run 10 stress 0.06123362 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01521342  max resid 0.04230428 
    ## Run 11 stress 0.07411943 
    ## Run 12 stress 0.0647795 
    ## Run 13 stress 0.06478001 
    ## Run 14 stress 0.06124508 
    ## ... Procrustes: rmse 0.01514247  max resid 0.04192778 
    ## Run 15 stress 0.06623459 
    ## Run 16 stress 0.06623467 
    ## Run 17 stress 0.07166983 
    ## Run 18 stress 0.07411941 
    ## Run 19 stress 0.06124503 
    ## ... Procrustes: rmse 0.01509306  max resid 0.04179647 
    ## Run 20 stress 0.2362366 
    ## Run 21 stress 0.0741194 
    ## Run 22 stress 0.2802425 
    ## Run 23 stress 0.06123361 
    ## ... New best solution
    ## ... Procrustes: rmse 5.078527e-05  max resid 0.0001207847 
    ## ... Similar to previous best
    ## Run 24 stress 0.07411938 
    ## Run 25 stress 0.0662346 
    ## Run 26 stress 0.06623458 
    ## Run 27 stress 0.06124489 
    ## ... Procrustes: rmse 0.01487631  max resid 0.04121813 
    ## Run 28 stress 0.2345379 
    ## Run 29 stress 0.06124489 
    ## ... Procrustes: rmse 0.01488054  max resid 0.04122955 
    ## Run 30 stress 0.08226158 
    ## Run 31 stress 0.0741194 
    ## Run 32 stress 0.06477969 
    ## Run 33 stress 0.07411945 
    ## Run 34 stress 0.06124501 
    ## ... Procrustes: rmse 0.01504492  max resid 0.04166997 
    ## Run 35 stress 0.0716697 
    ## Run 36 stress 0.08226171 
    ## Run 37 stress 0.06477857 
    ## Run 38 stress 0.07166976 
    ## Run 39 stress 0.0741195 
    ## Run 40 stress 0.06623457 
    ## Run 41 stress 0.06623463 
    ## Run 42 stress 0.06623465 
    ## Run 43 stress 0.06477958 
    ## Run 44 stress 0.2362366 
    ## Run 45 stress 0.06124497 
    ## ... Procrustes: rmse 0.01495915  max resid 0.04143803 
    ## Run 46 stress 0.06623458 
    ## Run 47 stress 0.06477908 
    ## Run 48 stress 0.07166973 
    ## Run 49 stress 0.06623458 
    ## Run 50 stress 0.06124494 
    ## ... Procrustes: rmse 0.01496833  max resid 0.04146496 
    ## Run 51 stress 0.08226158 
    ## Run 52 stress 0.06477959 
    ## Run 53 stress 0.0741194 
    ## Run 54 stress 0.0741194 
    ## Run 55 stress 0.08226169 
    ## Run 56 stress 0.07411939 
    ## Run 57 stress 0.06477952 
    ## Run 58 stress 0.06623461 
    ## Run 59 stress 0.06623457 
    ## Run 60 stress 0.08226185 
    ## Run 61 stress 0.06623458 
    ## Run 62 stress 0.06623461 
    ## Run 63 stress 0.07166983 
    ## Run 64 stress 0.3266727 
    ## Run 65 stress 0.07166988 
    ## Run 66 stress 0.06123361 
    ## ... New best solution
    ## ... Procrustes: rmse 3.900723e-05  max resid 8.816599e-05 
    ## ... Similar to previous best
    ## Run 67 stress 0.0662346 
    ## Run 68 stress 0.08226165 
    ## Run 69 stress 0.08226157 
    ## Run 70 stress 0.06477968 
    ## Run 71 stress 0.06623458 
    ## Run 72 stress 0.06477893 
    ## Run 73 stress 0.0647787 
    ## Run 74 stress 0.06477852 
    ## Run 75 stress 0.07166977 
    ## Run 76 stress 0.06124522 
    ## ... Procrustes: rmse 0.0152335  max resid 0.042173 
    ## Run 77 stress 0.0662346 
    ## Run 78 stress 0.06623467 
    ## Run 79 stress 0.07411945 
    ## Run 80 stress 0.06623459 
    ## Run 81 stress 0.06623465 
    ## Run 82 stress 0.06477829 
    ## Run 83 stress 0.08226165 
    ## Run 84 stress 0.0741194 
    ## Run 85 stress 0.06477916 
    ## Run 86 stress 0.07411945 
    ## Run 87 stress 0.06477993 
    ## Run 88 stress 0.07411945 
    ## Run 89 stress 0.06623457 
    ## Run 90 stress 0.08226169 
    ## Run 91 stress 0.06124512 
    ## ... Procrustes: rmse 0.01517064  max resid 0.04200262 
    ## Run 92 stress 0.06124498 
    ## ... Procrustes: rmse 0.01503888  max resid 0.04165016 
    ## Run 93 stress 0.07411943 
    ## Run 94 stress 0.07411948 
    ## Run 95 stress 0.06477875 
    ## Run 96 stress 0.07166973 
    ## Run 97 stress 0.07411946 
    ## Run 98 stress 0.06124505 
    ## ... Procrustes: rmse 0.01440265  max resid 0.03994391 
    ## Run 99 stress 0.07411939 
    ## Run 100 stress 0.06477868 
    ## Run 101 stress 0.06123361 
    ## ... Procrustes: rmse 8.072661e-05  max resid 0.000103023 
    ## ... Similar to previous best
    ## Run 102 stress 0.06477899 
    ## Run 103 stress 0.06477953 
    ## Run 104 stress 0.06623458 
    ## Run 105 stress 0.3114145 
    ## Run 106 stress 0.06623457 
    ## Run 107 stress 0.2345379 
    ## Run 108 stress 0.06477877 
    ## Run 109 stress 0.07411948 
    ## Run 110 stress 0.07411952 
    ## Run 111 stress 0.06124493 
    ## ... Procrustes: rmse 0.01498492  max resid 0.04150604 
    ## Run 112 stress 0.09469136 
    ## Run 113 stress 0.07411938 
    ## Run 114 stress 0.06623458 
    ## Run 115 stress 0.06477831 
    ## Run 116 stress 0.06124504 
    ## ... Procrustes: rmse 0.01510467  max resid 0.04182608 
    ## Run 117 stress 0.06124494 
    ## ... Procrustes: rmse 0.01499038  max resid 0.04152149 
    ## Run 118 stress 0.2362366 
    ## Run 119 stress 0.07411949 
    ## Run 120 stress 0.06124513 
    ## ... Procrustes: rmse 0.01517552  max resid 0.04201745 
    ## Run 121 stress 0.07166969 
    ## Run 122 stress 0.07166978 
    ## Run 123 stress 0.06623459 
    ## Run 124 stress 0.06477876 
    ## Run 125 stress 0.2965857 
    ## Run 126 stress 0.06623468 
    ## Run 127 stress 0.06623457 
    ## Run 128 stress 0.07411938 
    ## Run 129 stress 0.07166973 
    ## Run 130 stress 0.07166964 
    ## Run 131 stress 0.06123361 
    ## ... Procrustes: rmse 9.223424e-05  max resid 0.0001203989 
    ## ... Similar to previous best
    ## Run 132 stress 0.07411939 
    ## Run 133 stress 0.08226161 
    ## Run 134 stress 0.06623462 
    ## Run 135 stress 0.06623458 
    ## Run 136 stress 0.07166969 
    ## Run 137 stress 0.07166974 
    ## Run 138 stress 0.06623458 
    ## Run 139 stress 0.06623459 
    ## Run 140 stress 0.06124514 
    ## ... Procrustes: rmse 0.01518599  max resid 0.04204481 
    ## Run 141 stress 0.06623459 
    ## Run 142 stress 0.2520061 
    ## Run 143 stress 0.06477917 
    ## Run 144 stress 0.3181213 
    ## Run 145 stress 0.06124498 
    ## ... Procrustes: rmse 0.01448024  max resid 0.04015218 
    ## Run 146 stress 0.08226158 
    ## Run 147 stress 0.06623465 
    ## Run 148 stress 0.06124501 
    ## ... Procrustes: rmse 0.01507333  max resid 0.04174315 
    ## Run 149 stress 0.06477904 
    ## Run 150 stress 0.07411938 
    ## Run 151 stress 0.07411942 
    ## Run 152 stress 0.07411939 
    ## Run 153 stress 0.2324715 
    ## Run 154 stress 0.06477832 
    ## Run 155 stress 0.06623459 
    ## Run 156 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 5.1024e-05  max resid 6.180689e-05 
    ## ... Similar to previous best
    ## Run 157 stress 0.06623459 
    ## Run 158 stress 0.06623457 
    ## Run 159 stress 0.07411941 
    ## Run 160 stress 0.06124507 
    ## ... Procrustes: rmse 0.01515229  max resid 0.04195498 
    ## Run 161 stress 0.07166966 
    ## Run 162 stress 0.2362366 
    ## Run 163 stress 0.07411947 
    ## Run 164 stress 0.06477939 
    ## Run 165 stress 0.06623457 
    ## Run 166 stress 0.08226172 
    ## Run 167 stress 0.06123361 
    ## ... Procrustes: rmse 2.997082e-05  max resid 4.557891e-05 
    ## ... Similar to previous best
    ## Run 168 stress 0.06477901 
    ## Run 169 stress 0.06477878 
    ## Run 170 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 1.025509e-05  max resid 1.245947e-05 
    ## ... Similar to previous best
    ## Run 171 stress 0.06124504 
    ## ... Procrustes: rmse 0.01511729  max resid 0.04186104 
    ## Run 172 stress 0.07411938 
    ## Run 173 stress 0.06477875 
    ## Run 174 stress 0.06623458 
    ## Run 175 stress 0.07411948 
    ## Run 176 stress 0.3151248 
    ## Run 177 stress 0.0612336 
    ## ... Procrustes: rmse 1.142758e-05  max resid 2.456059e-05 
    ## ... Similar to previous best
    ## Run 178 stress 0.06477979 
    ## Run 179 stress 0.06623463 
    ## Run 180 stress 0.07166967 
    ## Run 181 stress 0.08226156 
    ## Run 182 stress 0.0662346 
    ## Run 183 stress 0.06477856 
    ## Run 184 stress 0.07411946 
    ## Run 185 stress 0.07411945 
    ## Run 186 stress 0.07166967 
    ## Run 187 stress 0.06477863 
    ## Run 188 stress 0.0647788 
    ## Run 189 stress 0.07411944 
    ## Run 190 stress 0.06123361 
    ## ... Procrustes: rmse 2.912609e-05  max resid 4.033121e-05 
    ## ... Similar to previous best
    ## Run 191 stress 0.06623457 
    ## Run 192 stress 0.07166967 
    ## Run 193 stress 0.0716697 
    ## Run 194 stress 0.06123361 
    ## ... Procrustes: rmse 1.811996e-05  max resid 2.66049e-05 
    ## ... Similar to previous best
    ## Run 195 stress 0.07411945 
    ## Run 196 stress 0.06477881 
    ## Run 197 stress 0.07411941 
    ## Run 198 stress 0.061245 
    ## ... Procrustes: rmse 0.01508317  max resid 0.04177006 
    ## Run 199 stress 0.0662346 
    ## Run 200 stress 0.0647787 
    ## Run 201 stress 0.340734 
    ## Run 202 stress 0.0612451 
    ## ... Procrustes: rmse 0.01438228  max resid 0.03988914 
    ## Run 203 stress 0.07166965 
    ## Run 204 stress 0.06623466 
    ## Run 205 stress 0.06623459 
    ## Run 206 stress 0.08226176 
    ## Run 207 stress 0.06623459 
    ## Run 208 stress 0.07411948 
    ## Run 209 stress 0.06124508 
    ## ... Procrustes: rmse 0.01515082  max resid 0.04195056 
    ## Run 210 stress 0.07411939 
    ## Run 211 stress 0.07411955 
    ## Run 212 stress 0.0662347 
    ## Run 213 stress 0.06477943 
    ## Run 214 stress 0.06477891 
    ## Run 215 stress 0.06477859 
    ## Run 216 stress 0.08226169 
    ## Run 217 stress 0.3157917 
    ## Run 218 stress 0.07166987 
    ## Run 219 stress 0.07411939 
    ## Run 220 stress 0.2381016 
    ## Run 221 stress 0.06477889 
    ## Run 222 stress 0.07411944 
    ## Run 223 stress 0.07411938 
    ## Run 224 stress 0.07411954 
    ## Run 225 stress 0.2362366 
    ## Run 226 stress 0.06123361 
    ## ... Procrustes: rmse 2.253497e-05  max resid 3.967582e-05 
    ## ... Similar to previous best
    ## Run 227 stress 0.07411946 
    ## Run 228 stress 0.06124505 
    ## ... Procrustes: rmse 0.01513002  max resid 0.0418951 
    ## Run 229 stress 0.07166974 
    ## Run 230 stress 0.0716697 
    ## Run 231 stress 0.2520056 
    ## Run 232 stress 0.06623464 
    ## Run 233 stress 0.061245 
    ## ... Procrustes: rmse 0.01508435  max resid 0.04177378 
    ## Run 234 stress 0.06623459 
    ## Run 235 stress 0.06477951 
    ## Run 236 stress 0.06123364 
    ## ... Procrustes: rmse 4.608394e-05  max resid 7.027788e-05 
    ## ... Similar to previous best
    ## Run 237 stress 0.06623461 
    ## Run 238 stress 0.06123361 
    ## ... Procrustes: rmse 9.694182e-06  max resid 1.799225e-05 
    ## ... Similar to previous best
    ## Run 239 stress 0.07166964 
    ## Run 240 stress 0.07411938 
    ## Run 241 stress 0.06477897 
    ## Run 242 stress 0.06623458 
    ## Run 243 stress 0.06623459 
    ## Run 244 stress 0.06124523 
    ## ... Procrustes: rmse 0.01526377  max resid 0.04225331 
    ## Run 245 stress 0.07411945 
    ## Run 246 stress 0.06124505 
    ## ... Procrustes: rmse 0.01513156  max resid 0.0418993 
    ## Run 247 stress 0.08226167 
    ## Run 248 stress 0.0741194 
    ## Run 249 stress 0.07166977 
    ## Run 250 stress 0.06124504 
    ## ... Procrustes: rmse 0.01511876  max resid 0.04186489 
    ## Run 251 stress 0.0662346 
    ## Run 252 stress 0.06477896 
    ## Run 253 stress 0.06124504 
    ## ... Procrustes: rmse 0.01512173  max resid 0.04187303 
    ## Run 254 stress 0.06124499 
    ## ... Procrustes: rmse 0.0150667  max resid 0.04172597 
    ## Run 255 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 1.982557e-06  max resid 3.490409e-06 
    ## ... Similar to previous best
    ## Run 256 stress 0.06477881 
    ## Run 257 stress 0.07411941 
    ## Run 258 stress 0.08226156 
    ## Run 259 stress 0.08226158 
    ## Run 260 stress 0.06623457 
    ## Run 261 stress 0.06124492 
    ## ... Procrustes: rmse 0.01498324  max resid 0.04150291 
    ## Run 262 stress 0.07411938 
    ## Run 263 stress 0.06124549 
    ## ... Procrustes: rmse 0.01417696  max resid 0.03934134 
    ## Run 264 stress 0.08226156 
    ## Run 265 stress 0.0647798 
    ## Run 266 stress 0.06477839 
    ## Run 267 stress 0.06477993 
    ## Run 268 stress 0.07411939 
    ## Run 269 stress 0.0822616 
    ## Run 270 stress 0.06124511 
    ## ... Procrustes: rmse 0.01436756  max resid 0.03985091 
    ## Run 271 stress 0.06623461 
    ## Run 272 stress 0.08226161 
    ## Run 273 stress 0.06477858 
    ## Run 274 stress 0.06124496 
    ## ... Procrustes: rmse 0.01503859  max resid 0.04165151 
    ## Run 275 stress 0.0741195 
    ## Run 276 stress 0.06123361 
    ## ... Procrustes: rmse 3.138202e-05  max resid 4.123497e-05 
    ## ... Similar to previous best
    ## Run 277 stress 0.06477831 
    ## Run 278 stress 0.0612451 
    ## ... Procrustes: rmse 0.01517415  max resid 0.04201329 
    ## Run 279 stress 0.06477962 
    ## Run 280 stress 0.0741194 
    ## Run 281 stress 0.06623458 
    ## Run 282 stress 0.07166976 
    ## Run 283 stress 0.06477971 
    ## Run 284 stress 0.06623462 
    ## Run 285 stress 0.07411947 
    ## Run 286 stress 0.06477934 
    ## Run 287 stress 0.06124508 
    ## ... Procrustes: rmse 0.01515455  max resid 0.04196041 
    ## Run 288 stress 0.06124503 
    ## ... Procrustes: rmse 0.01511379  max resid 0.04185166 
    ## Run 289 stress 0.06477885 
    ## Run 290 stress 0.07411946 
    ## Run 291 stress 0.07166982 
    ## Run 292 stress 0.0930082 
    ## Run 293 stress 0.06123361 
    ## ... Procrustes: rmse 1.429227e-05  max resid 2.608923e-05 
    ## ... Similar to previous best
    ## Run 294 stress 0.06124512 
    ## ... Procrustes: rmse 0.01436084  max resid 0.03983335 
    ## Run 295 stress 0.07166973 
    ## Run 296 stress 0.06124506 
    ## ... Procrustes: rmse 0.01513848  max resid 0.04191868 
    ## Run 297 stress 0.07411946 
    ## Run 298 stress 0.07411944 
    ## Run 299 stress 0.06123361 
    ## ... Procrustes: rmse 4.742072e-05  max resid 6.37204e-05 
    ## ... Similar to previous best
    ## Run 300 stress 0.07411944 
    ## Run 301 stress 0.06623459 
    ## Run 302 stress 0.0612452 
    ## ... Procrustes: rmse 0.01524591  max resid 0.04220535 
    ## Run 303 stress 0.08226157 
    ## Run 304 stress 0.07411942 
    ## Run 305 stress 0.07411954 
    ## Run 306 stress 0.06623457 
    ## Run 307 stress 0.06623458 
    ## Run 308 stress 0.0741194 
    ## Run 309 stress 0.06123361 
    ## ... Procrustes: rmse 5.438101e-05  max resid 7.257848e-05 
    ## ... Similar to previous best
    ## Run 310 stress 0.07166971 
    ## Run 311 stress 0.06623457 
    ## Run 312 stress 0.07166964 
    ## Run 313 stress 0.07411942 
    ## Run 314 stress 0.06477897 
    ## Run 315 stress 0.08226164 
    ## Run 316 stress 0.06623458 
    ## Run 317 stress 0.06477864 
    ## Run 318 stress 0.07411946 
    ## Run 319 stress 0.06623458 
    ## Run 320 stress 0.0662346 
    ## Run 321 stress 0.06623459 
    ## Run 322 stress 0.07166981 
    ## Run 323 stress 0.0647795 
    ## Run 324 stress 0.2361811 
    ## Run 325 stress 0.07166965 
    ## Run 326 stress 0.06623462 
    ## Run 327 stress 0.07166972 
    ## Run 328 stress 0.07411938 
    ## Run 329 stress 0.06477834 
    ## Run 330 stress 0.06623458 
    ## Run 331 stress 0.06477876 
    ## Run 332 stress 0.0612336 
    ## ... Procrustes: rmse 2.758574e-06  max resid 4.722824e-06 
    ## ... Similar to previous best
    ## Run 333 stress 0.07411946 
    ## Run 334 stress 0.0741195 
    ## Run 335 stress 0.07411945 
    ## Run 336 stress 0.06623458 
    ## Run 337 stress 0.07411949 
    ## Run 338 stress 0.0612336 
    ## ... Procrustes: rmse 7.297646e-06  max resid 9.871721e-06 
    ## ... Similar to previous best
    ## Run 339 stress 0.06623457 
    ## Run 340 stress 0.06623457 
    ## Run 341 stress 0.07166974 
    ## Run 342 stress 0.07411939 
    ## Run 343 stress 0.08226159 
    ## Run 344 stress 0.07411942 
    ## Run 345 stress 0.06124491 
    ## ... Procrustes: rmse 0.01496124  max resid 0.04144448 
    ## Run 346 stress 0.0612449 
    ## ... Procrustes: rmse 0.01495053  max resid 0.04141446 
    ## Run 347 stress 0.3019511 
    ## Run 348 stress 0.06623459 
    ## Run 349 stress 0.06623463 
    ## Run 350 stress 0.2324715 
    ## Run 351 stress 0.06124501 
    ## ... Procrustes: rmse 0.015095  max resid 0.04180201 
    ## Run 352 stress 0.07166968 
    ## Run 353 stress 0.2548129 
    ## Run 354 stress 0.06623457 
    ## Run 355 stress 0.06623458 
    ## Run 356 stress 0.06123361 
    ## ... Procrustes: rmse 2.028732e-05  max resid 2.953833e-05 
    ## ... Similar to previous best
    ## Run 357 stress 0.06477844 
    ## Run 358 stress 0.06123361 
    ## ... Procrustes: rmse 7.110764e-06  max resid 1.546362e-05 
    ## ... Similar to previous best
    ## Run 359 stress 0.0741195 
    ## Run 360 stress 0.06623461 
    ## Run 361 stress 0.06477827 
    ## Run 362 stress 0.07411952 
    ## Run 363 stress 0.06477983 
    ## Run 364 stress 0.06123362 
    ## ... Procrustes: rmse 5.326946e-05  max resid 7.494647e-05 
    ## ... Similar to previous best
    ## Run 365 stress 0.06623461 
    ## Run 366 stress 0.07166964 
    ## Run 367 stress 0.0662346 
    ## Run 368 stress 0.07411938 
    ## Run 369 stress 0.06477908 
    ## Run 370 stress 0.06124513 
    ## ... Procrustes: rmse 0.01519347  max resid 0.04206492 
    ## Run 371 stress 0.3407337 
    ## Run 372 stress 0.061245 
    ## ... Procrustes: rmse 0.01507597  max resid 0.04175092 
    ## Run 373 stress 0.07411949 
    ## Run 374 stress 0.07411941 
    ## Run 375 stress 0.2362366 
    ## Run 376 stress 0.08226167 
    ## Run 377 stress 0.07166978 
    ## Run 378 stress 0.07166975 
    ## Run 379 stress 0.06124497 
    ## ... Procrustes: rmse 0.01504033  max resid 0.04165568 
    ## Run 380 stress 0.0662346 
    ## Run 381 stress 0.06623465 
    ## Run 382 stress 0.06623462 
    ## Run 383 stress 0.06623461 
    ## Run 384 stress 0.06477898 
    ## Run 385 stress 0.0741194 
    ## Run 386 stress 0.06124526 
    ## ... Procrustes: rmse 0.01528185  max resid 0.04230112 
    ## Run 387 stress 0.07411943 
    ## Run 388 stress 0.06477941 
    ## Run 389 stress 0.06623463 
    ## Run 390 stress 0.0662347 
    ## Run 391 stress 0.0662346 
    ## Run 392 stress 0.06124513 
    ## ... Procrustes: rmse 0.01512809  max resid 0.04189372 
    ## Run 393 stress 0.07166978 
    ## Run 394 stress 0.06478012 
    ## Run 395 stress 0.09300823 
    ## Run 396 stress 0.07166987 
    ## Run 397 stress 0.07166978 
    ## Run 398 stress 0.07411945 
    ## Run 399 stress 0.07166971 
    ## Run 400 stress 0.07411952 
    ## Run 401 stress 0.06124534 
    ## ... Procrustes: rmse 0.01423419  max resid 0.03949447 
    ## Run 402 stress 0.06477907 
    ## Run 403 stress 0.06623457 
    ## Run 404 stress 0.06623457 
    ## Run 405 stress 0.08226167 
    ## Run 406 stress 0.06623458 
    ## Run 407 stress 0.06477906 
    ## Run 408 stress 0.06623461 
    ## Run 409 stress 0.0612451 
    ## ... Procrustes: rmse 0.01515177  max resid 0.04195207 
    ## Run 410 stress 0.06123361 
    ## ... Procrustes: rmse 4.586765e-05  max resid 6.014403e-05 
    ## ... Similar to previous best
    ## Run 411 stress 0.06124532 
    ## ... Procrustes: rmse 0.01527118  max resid 0.04227642 
    ## Run 412 stress 0.06124522 
    ## ... Procrustes: rmse 0.01525544  max resid 0.04223076 
    ## Run 413 stress 0.06124486 
    ## ... Procrustes: rmse 0.01473962  max resid 0.04084984 
    ## Run 414 stress 0.06123361 
    ## ... Procrustes: rmse 4.578314e-05  max resid 8.184885e-05 
    ## ... Similar to previous best
    ## Run 415 stress 0.08226158 
    ## Run 416 stress 0.2361811 
    ## Run 417 stress 0.06477884 
    ## Run 418 stress 0.06124514 
    ## ... Procrustes: rmse 0.01519841  max resid 0.04207807 
    ## Run 419 stress 0.06623459 
    ## Run 420 stress 0.06623457 
    ## Run 421 stress 0.0741194 
    ## Run 422 stress 0.06123362 
    ## ... Procrustes: rmse 6.705351e-05  max resid 9.044328e-05 
    ## ... Similar to previous best
    ## Run 423 stress 0.08226168 
    ## Run 424 stress 0.07411939 
    ## Run 425 stress 0.0822616 
    ## Run 426 stress 0.06623457 
    ## Run 427 stress 0.06123364 
    ## ... Procrustes: rmse 9.150789e-05  max resid 0.0001139568 
    ## ... Similar to previous best
    ## Run 428 stress 0.0741194 
    ## Run 429 stress 0.06124489 
    ## ... Procrustes: rmse 0.01492794  max resid 0.04135473 
    ## Run 430 stress 0.06124517 
    ## ... Procrustes: rmse 0.01521555  max resid 0.04212536 
    ## Run 431 stress 0.06477897 
    ## Run 432 stress 0.0716697 
    ## Run 433 stress 0.06623458 
    ## Run 434 stress 0.07411941 
    ## Run 435 stress 0.06623462 
    ## Run 436 stress 0.06124512 
    ## ... Procrustes: rmse 0.01518272  max resid 0.04203604 
    ## Run 437 stress 0.06623457 
    ## Run 438 stress 0.06124505 
    ## ... Procrustes: rmse 0.01512933  max resid 0.04189415 
    ## Run 439 stress 0.0741194 
    ## Run 440 stress 0.06124527 
    ## ... Procrustes: rmse 0.01523055  max resid 0.04216792 
    ## Run 441 stress 0.07411947 
    ## Run 442 stress 0.0647787 
    ## Run 443 stress 0.07411941 
    ## Run 444 stress 0.07411944 
    ## Run 445 stress 0.06124517 
    ## ... Procrustes: rmse 0.0152233  max resid 0.04214495 
    ## Run 446 stress 0.07411946 
    ## Run 447 stress 0.06477892 
    ## Run 448 stress 0.06623458 
    ## Run 449 stress 0.08226159 
    ## Run 450 stress 0.06123362 
    ## ... Procrustes: rmse 4.038577e-05  max resid 5.440988e-05 
    ## ... Similar to previous best
    ## Run 451 stress 0.06477969 
    ## Run 452 stress 0.2381016 
    ## Run 453 stress 0.06124502 
    ## ... Procrustes: rmse 0.01509954  max resid 0.04181428 
    ## Run 454 stress 0.07411942 
    ## Run 455 stress 0.08226158 
    ## Run 456 stress 0.07166985 
    ## Run 457 stress 0.07411943 
    ## Run 458 stress 0.06477901 
    ## Run 459 stress 0.06477865 
    ## Run 460 stress 0.07166986 
    ## Run 461 stress 0.07166967 
    ## Run 462 stress 0.06124491 
    ## ... Procrustes: rmse 0.01459828  max resid 0.04047179 
    ## Run 463 stress 0.07166968 
    ## Run 464 stress 0.07166981 
    ## Run 465 stress 0.06124504 
    ## ... Procrustes: rmse 0.01446316  max resid 0.04010532 
    ## Run 466 stress 0.07411945 
    ## Run 467 stress 0.06123362 
    ## ... Procrustes: rmse 6.120052e-05  max resid 8.242623e-05 
    ## ... Similar to previous best
    ## Run 468 stress 0.06478006 
    ## Run 469 stress 0.3407347 
    ## Run 470 stress 0.08226178 
    ## Run 471 stress 0.07166964 
    ## Run 472 stress 0.06124493 
    ## ... Procrustes: rmse 0.01498813  max resid 0.04151457 
    ## Run 473 stress 0.0741195 
    ## Run 474 stress 0.06477843 
    ## Run 475 stress 0.06477957 
    ## Run 476 stress 0.07411942 
    ## Run 477 stress 0.06623469 
    ## Run 478 stress 0.09469042 
    ## Run 479 stress 0.07411945 
    ## Run 480 stress 0.08226165 
    ## Run 481 stress 0.07411947 
    ## Run 482 stress 0.07166972 
    ## Run 483 stress 0.06623463 
    ## Run 484 stress 0.06477966 
    ## Run 485 stress 0.07166971 
    ## Run 486 stress 0.07411951 
    ## Run 487 stress 0.06623458 
    ## Run 488 stress 0.07166979 
    ## Run 489 stress 0.3407336 
    ## Run 490 stress 0.06124507 
    ## ... Procrustes: rmse 0.01514376  max resid 0.04193204 
    ## Run 491 stress 0.07411948 
    ## Run 492 stress 0.06124507 
    ## ... Procrustes: rmse 0.01513996  max resid 0.04192309 
    ## Run 493 stress 0.06477832 
    ## Run 494 stress 0.06123361 
    ## ... Procrustes: rmse 4.131689e-05  max resid 6.191789e-05 
    ## ... Similar to previous best
    ## Run 495 stress 0.06477952 
    ## Run 496 stress 0.08226178 
    ## Run 497 stress 0.06477904 
    ## Run 498 stress 0.0612454 
    ## ... Procrustes: rmse 0.01421223  max resid 0.03942972 
    ## Run 499 stress 0.2518588 
    ## Run 500 stress 0.06124497 
    ## ... Procrustes: rmse 0.01504685  max resid 0.04167306 
    ## *** Best solution repeated 17 times

``` r
# Stratified lakes and ocean sites
SD_beta_env_SO_NMDS <- metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.59483e-05 
    ## Run 1 stress 9.739026e-05 
    ## ... Procrustes: rmse 0.0001838397  max resid 0.0004639504 
    ## ... Similar to previous best
    ## Run 2 stress 9.491934e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001782711  max resid 0.000446651 
    ## ... Similar to previous best
    ## Run 3 stress 9.177963e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001684101  max resid 0.0002856731 
    ## ... Similar to previous best
    ## Run 4 stress 9.843344e-05 
    ## ... Procrustes: rmse 0.0002180475  max resid 0.0003611731 
    ## ... Similar to previous best
    ## Run 5 stress 9.219632e-05 
    ## ... Procrustes: rmse 0.0001435806  max resid 0.0002760377 
    ## ... Similar to previous best
    ## Run 6 stress 9.683031e-05 
    ## ... Procrustes: rmse 1.264943e-05  max resid 2.54755e-05 
    ## ... Similar to previous best
    ## Run 7 stress 9.938598e-05 
    ## ... Procrustes: rmse 1.969842e-05  max resid 3.455737e-05 
    ## ... Similar to previous best
    ## Run 8 stress 8.833905e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000165862  max resid 0.0002729302 
    ## ... Similar to previous best
    ## Run 9 stress 9.747773e-05 
    ## ... Procrustes: rmse 0.0002229855  max resid 0.0004799202 
    ## ... Similar to previous best
    ## Run 10 stress 9.585684e-05 
    ## ... Procrustes: rmse 9.503068e-05  max resid 0.0002151741 
    ## ... Similar to previous best
    ## Run 11 stress 9.552015e-05 
    ## ... Procrustes: rmse 0.0001666956  max resid 0.0003460856 
    ## ... Similar to previous best
    ## Run 12 stress 9.975023e-05 
    ## ... Procrustes: rmse 0.0001431305  max resid 0.000351508 
    ## ... Similar to previous best
    ## Run 13 stress 9.801228e-05 
    ## ... Procrustes: rmse 0.0002232427  max resid 0.0004799616 
    ## ... Similar to previous best
    ## Run 14 stress 9.99643e-05 
    ## ... Procrustes: rmse 9.434606e-05  max resid 0.0002267314 
    ## ... Similar to previous best
    ## Run 15 stress 9.765205e-05 
    ## ... Procrustes: rmse 1.781092e-05  max resid 3.304114e-05 
    ## ... Similar to previous best
    ## Run 16 stress 9.595076e-05 
    ## ... Procrustes: rmse 0.0001668825  max resid 0.0003223807 
    ## ... Similar to previous best
    ## Run 17 stress 9.738654e-05 
    ## ... Procrustes: rmse 0.0001805719  max resid 0.0004629331 
    ## ... Similar to previous best
    ## Run 18 stress 9.778487e-05 
    ## ... Procrustes: rmse 0.0001733177  max resid 0.0002795738 
    ## ... Similar to previous best
    ## Run 19 stress 9.733011e-05 
    ## ... Procrustes: rmse 0.0001485777  max resid 0.0002292773 
    ## ... Similar to previous best
    ## Run 20 stress 9.701714e-05 
    ## ... Procrustes: rmse 0.0001795636  max resid 0.0004491849 
    ## ... Similar to previous best
    ## Run 21 stress 9.886768e-05 
    ## ... Procrustes: rmse 0.0002073127  max resid 0.0004340273 
    ## ... Similar to previous best
    ## Run 22 stress 9.269463e-05 
    ## ... Procrustes: rmse 0.0001730235  max resid 0.0003177731 
    ## ... Similar to previous best
    ## Run 23 stress 9.166485e-05 
    ## ... Procrustes: rmse 0.0001359578  max resid 0.0003366578 
    ## ... Similar to previous best
    ## Run 24 stress 9.243048e-05 
    ## ... Procrustes: rmse 1.06868e-05  max resid 1.447145e-05 
    ## ... Similar to previous best
    ## Run 25 stress 9.063609e-05 
    ## ... Procrustes: rmse 0.0001631409  max resid 0.0002691165 
    ## ... Similar to previous best
    ## Run 26 stress 8.665191e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001640908  max resid 0.0003257745 
    ## ... Similar to previous best
    ## Run 27 stress 9.991848e-05 
    ## ... Procrustes: rmse 9.497916e-05  max resid 0.0002330046 
    ## ... Similar to previous best
    ## Run 28 stress 9.645129e-05 
    ## ... Procrustes: rmse 0.0001511298  max resid 0.0003531831 
    ## ... Similar to previous best
    ## Run 29 stress 9.363412e-05 
    ## ... Procrustes: rmse 0.0001519582  max resid 0.0003445903 
    ## ... Similar to previous best
    ## Run 30 stress 9.784739e-05 
    ## ... Procrustes: rmse 0.0001332297  max resid 0.0003173071 
    ## ... Similar to previous best
    ## Run 31 stress 9.544579e-05 
    ## ... Procrustes: rmse 0.000201447  max resid 0.0004366977 
    ## ... Similar to previous best
    ## Run 32 stress 9.410021e-05 
    ## ... Procrustes: rmse 0.0001979444  max resid 0.0004279843 
    ## ... Similar to previous best
    ## Run 33 stress 9.502593e-05 
    ## ... Procrustes: rmse 0.0001870246  max resid 0.0003507502 
    ## ... Similar to previous best
    ## Run 34 stress 9.231733e-05 
    ## ... Procrustes: rmse 0.0001613153  max resid 0.0003239504 
    ## ... Similar to previous best
    ## Run 35 stress 9.496605e-05 
    ## ... Procrustes: rmse 0.0001577375  max resid 0.0003623737 
    ## ... Similar to previous best
    ## Run 36 stress 9.181635e-05 
    ## ... Procrustes: rmse 0.0001664047  max resid 0.0003344401 
    ## ... Similar to previous best
    ## Run 37 stress 9.127964e-05 
    ## ... Procrustes: rmse 2.121403e-05  max resid 3.315408e-05 
    ## ... Similar to previous best
    ## Run 38 stress 9.857421e-05 
    ## ... Procrustes: rmse 0.0001670431  max resid 0.0003561348 
    ## ... Similar to previous best
    ## Run 39 stress 9.541818e-05 
    ## ... Procrustes: rmse 0.0001590093  max resid 0.0003494481 
    ## ... Similar to previous best
    ## Run 40 stress 9.563577e-05 
    ## ... Procrustes: rmse 6.522085e-05  max resid 0.0001543281 
    ## ... Similar to previous best
    ## Run 41 stress 9.791148e-05 
    ## ... Procrustes: rmse 0.0001504513  max resid 0.0003317234 
    ## ... Similar to previous best
    ## Run 42 stress 9.398705e-05 
    ## ... Procrustes: rmse 0.000169759  max resid 0.0003394197 
    ## ... Similar to previous best
    ## Run 43 stress 9.729983e-05 
    ## ... Procrustes: rmse 0.00019438  max resid 0.0003621783 
    ## ... Similar to previous best
    ## Run 44 stress 9.731753e-05 
    ## ... Procrustes: rmse 0.0001342373  max resid 0.0003184167 
    ## ... Similar to previous best
    ## Run 45 stress 9.886034e-05 
    ## ... Procrustes: rmse 0.0001979603  max resid 0.0004422729 
    ## ... Similar to previous best
    ## Run 46 stress 9.774925e-05 
    ## ... Procrustes: rmse 0.0001798986  max resid 0.0003582273 
    ## ... Similar to previous best
    ## Run 47 stress 9.615516e-05 
    ## ... Procrustes: rmse 0.0001858358  max resid 0.0003886043 
    ## ... Similar to previous best
    ## Run 48 stress 9.592124e-05 
    ## ... Procrustes: rmse 0.0001718216  max resid 0.0003403636 
    ## ... Similar to previous best
    ## Run 49 stress 9.410932e-05 
    ## ... Procrustes: rmse 0.0001150837  max resid 0.0002099223 
    ## ... Similar to previous best
    ## Run 50 stress 9.972659e-05 
    ## ... Procrustes: rmse 0.0001231174  max resid 0.0002203627 
    ## ... Similar to previous best
    ## Run 51 stress 9.325319e-05 
    ## ... Procrustes: rmse 0.0001994549  max resid 0.0004311353 
    ## ... Similar to previous best
    ## Run 52 stress 9.700932e-05 
    ## ... Procrustes: rmse 0.0001627245  max resid 0.0003050688 
    ## ... Similar to previous best
    ## Run 53 stress 9.77126e-05 
    ## ... Procrustes: rmse 0.0001734948  max resid 0.0003438108 
    ## ... Similar to previous best
    ## Run 54 stress 9.811747e-05 
    ## ... Procrustes: rmse 5.917882e-05  max resid 9.611429e-05 
    ## ... Similar to previous best
    ## Run 55 stress 9.917248e-05 
    ## ... Procrustes: rmse 0.0001962721  max resid 0.0003647478 
    ## ... Similar to previous best
    ## Run 56 stress 9.262398e-05 
    ## ... Procrustes: rmse 0.0001375946  max resid 0.0003171432 
    ## ... Similar to previous best
    ## Run 57 stress 9.210221e-05 
    ## ... Procrustes: rmse 0.0001521007  max resid 0.0003429077 
    ## ... Similar to previous best
    ## Run 58 stress 9.863138e-05 
    ## ... Procrustes: rmse 0.0001614811  max resid 0.0003211915 
    ## ... Similar to previous best
    ## Run 59 stress 9.912773e-05 
    ## ... Procrustes: rmse 0.0001872029  max resid 0.0003261298 
    ## ... Similar to previous best
    ## Run 60 stress 9.650337e-05 
    ## ... Procrustes: rmse 0.000116708  max resid 0.0002105603 
    ## ... Similar to previous best
    ## Run 61 stress 9.902066e-05 
    ## ... Procrustes: rmse 0.0001827826  max resid 0.0003428565 
    ## ... Similar to previous best
    ## Run 62 stress 9.999497e-05 
    ## ... Procrustes: rmse 0.0001641188  max resid 0.0003724111 
    ## ... Similar to previous best
    ## Run 63 stress 9.529874e-05 
    ## ... Procrustes: rmse 0.0001688661  max resid 0.000334107 
    ## ... Similar to previous best
    ## Run 64 stress 7.655096e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001573493  max resid 0.0003103832 
    ## ... Similar to previous best
    ## Run 65 stress 9.299191e-05 
    ## ... Procrustes: rmse 0.0002135278  max resid 0.0003859188 
    ## ... Similar to previous best
    ## Run 66 stress 9.589414e-05 
    ## ... Procrustes: rmse 0.0001841312  max resid 0.0002682963 
    ## ... Similar to previous best
    ## Run 67 stress 8.447991e-05 
    ## ... Procrustes: rmse 0.0001556448  max resid 0.0002318021 
    ## ... Similar to previous best
    ## Run 68 stress 9.366654e-05 
    ## ... Procrustes: rmse 0.0001417477  max resid 0.0002638739 
    ## ... Similar to previous best
    ## Run 69 stress 8.893431e-05 
    ## ... Procrustes: rmse 0.0001574352  max resid 0.0002320691 
    ## ... Similar to previous best
    ## Run 70 stress 9.610284e-05 
    ## ... Procrustes: rmse 0.0001925662  max resid 0.0003637764 
    ## ... Similar to previous best
    ## Run 71 stress 9.714978e-05 
    ## ... Procrustes: rmse 0.0001625573  max resid 0.000344638 
    ## ... Similar to previous best
    ## Run 72 stress 9.468055e-05 
    ## ... Procrustes: rmse 0.0001690334  max resid 0.0002743705 
    ## ... Similar to previous best
    ## Run 73 stress 9.823147e-05 
    ## ... Procrustes: rmse 0.0001953741  max resid 0.0003592921 
    ## ... Similar to previous best
    ## Run 74 stress 9.75185e-05 
    ## ... Procrustes: rmse 0.0002051702  max resid 0.0003548338 
    ## ... Similar to previous best
    ## Run 75 stress 9.795479e-05 
    ## ... Procrustes: rmse 0.0001965865  max resid 0.0003533852 
    ## ... Similar to previous best
    ## Run 76 stress 9.976507e-05 
    ## ... Procrustes: rmse 0.000199856  max resid 0.0003669482 
    ## ... Similar to previous best
    ## Run 77 stress 9.275031e-05 
    ## ... Procrustes: rmse 0.0001845266  max resid 0.0003525302 
    ## ... Similar to previous best
    ## Run 78 stress 9.420133e-05 
    ## ... Procrustes: rmse 0.0001643885  max resid 0.0002458808 
    ## ... Similar to previous best
    ## Run 79 stress 9.453037e-05 
    ## ... Procrustes: rmse 0.0001989596  max resid 0.000348788 
    ## ... Similar to previous best
    ## Run 80 stress 9.344885e-05 
    ## ... Procrustes: rmse 0.0001695935  max resid 0.0002556783 
    ## ... Similar to previous best
    ## Run 81 stress 9.756035e-05 
    ## ... Procrustes: rmse 0.0001780277  max resid 0.0002661691 
    ## ... Similar to previous best
    ## Run 82 stress 8.013286e-05 
    ## ... Procrustes: rmse 0.0001488158  max resid 0.0002738947 
    ## ... Similar to previous best
    ## Run 83 stress 9.268481e-05 
    ## ... Procrustes: rmse 0.000161144  max resid 0.0003314959 
    ## ... Similar to previous best
    ## Run 84 stress 8.120909e-05 
    ## ... Procrustes: rmse 0.000161802  max resid 0.0002836121 
    ## ... Similar to previous best
    ## Run 85 stress 9.511412e-05 
    ## ... Procrustes: rmse 0.000190948  max resid 0.0003291796 
    ## ... Similar to previous best
    ## Run 86 stress 8.904636e-05 
    ## ... Procrustes: rmse 0.0001564952  max resid 0.0002628181 
    ## ... Similar to previous best
    ## Run 87 stress 9.931034e-05 
    ## ... Procrustes: rmse 0.0001138201  max resid 0.0002457857 
    ## ... Similar to previous best
    ## Run 88 stress 9.430804e-05 
    ## ... Procrustes: rmse 0.0001894767  max resid 0.0003554709 
    ## ... Similar to previous best
    ## Run 89 stress 9.67453e-05 
    ## ... Procrustes: rmse 0.0001714846  max resid 0.0002572793 
    ## ... Similar to previous best
    ## Run 90 stress 9.74004e-05 
    ## ... Procrustes: rmse 0.0002184451  max resid 0.0003062336 
    ## ... Similar to previous best
    ## Run 91 stress 8.716354e-05 
    ## ... Procrustes: rmse 0.0001532016  max resid 0.0003260298 
    ## ... Similar to previous best
    ## Run 92 stress 9.540324e-05 
    ## ... Procrustes: rmse 0.0001843587  max resid 0.0002768679 
    ## ... Similar to previous best
    ## Run 93 stress 9.259824e-05 
    ## ... Procrustes: rmse 0.0001103284  max resid 0.0002124861 
    ## ... Similar to previous best
    ## Run 94 stress 9.794885e-05 
    ## ... Procrustes: rmse 0.0001853953  max resid 0.0002878856 
    ## ... Similar to previous best
    ## Run 95 stress 9.664387e-05 
    ## ... Procrustes: rmse 0.0001862892  max resid 0.0002707078 
    ## ... Similar to previous best
    ## Run 96 stress 9.693974e-05 
    ## ... Procrustes: rmse 0.0002191833  max resid 0.0002992543 
    ## ... Similar to previous best
    ## Run 97 stress 9.591798e-05 
    ## ... Procrustes: rmse 0.0001704484  max resid 0.0002584645 
    ## ... Similar to previous best
    ## Run 98 stress 9.900332e-05 
    ## ... Procrustes: rmse 0.0001959276  max resid 0.0003656059 
    ## ... Similar to previous best
    ## Run 99 stress 8.864277e-05 
    ## ... Procrustes: rmse 0.0001550857  max resid 0.0002549343 
    ## ... Similar to previous best
    ## Run 100 stress 8.894376e-05 
    ## ... Procrustes: rmse 0.0001373709  max resid 0.0002829759 
    ## ... Similar to previous best
    ## Run 101 stress 8.529538e-05 
    ## ... Procrustes: rmse 0.000161102  max resid 0.0002632623 
    ## ... Similar to previous best
    ## Run 102 stress 9.591838e-05 
    ## ... Procrustes: rmse 0.0001997019  max resid 0.0003498918 
    ## ... Similar to previous best
    ## Run 103 stress 9.751503e-05 
    ## ... Procrustes: rmse 0.0001739718  max resid 0.0002580829 
    ## ... Similar to previous best
    ## Run 104 stress 8.771165e-05 
    ## ... Procrustes: rmse 0.000142854  max resid 0.0002350219 
    ## ... Similar to previous best
    ## Run 105 stress 8.835331e-05 
    ## ... Procrustes: rmse 0.0001864536  max resid 0.0003330335 
    ## ... Similar to previous best
    ## Run 106 stress 7.615689e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001106308  max resid 0.0002238833 
    ## ... Similar to previous best
    ## Run 107 stress 9.707664e-05 
    ## ... Procrustes: rmse 0.0001453075  max resid 0.0002757907 
    ## ... Similar to previous best
    ## Run 108 stress 9.845097e-05 
    ## ... Procrustes: rmse 0.0001747324  max resid 0.0003072297 
    ## ... Similar to previous best
    ## Run 109 stress 8.494643e-05 
    ## ... Procrustes: rmse 0.0001855695  max resid 0.0003991701 
    ## ... Similar to previous best
    ## Run 110 stress 9.125344e-05 
    ## ... Procrustes: rmse 0.0002069267  max resid 0.000430303 
    ## ... Similar to previous best
    ## Run 111 stress 0.2280257 
    ## Run 112 stress 9.560985e-05 
    ## ... Procrustes: rmse 0.0001661243  max resid 0.0003188283 
    ## ... Similar to previous best
    ## Run 113 stress 9.828423e-05 
    ## ... Procrustes: rmse 0.0001757448  max resid 0.0004043598 
    ## ... Similar to previous best
    ## Run 114 stress 8.886209e-05 
    ## ... Procrustes: rmse 0.0001910114  max resid 0.0003942758 
    ## ... Similar to previous best
    ## Run 115 stress 9.439334e-05 
    ## ... Procrustes: rmse 0.0001690503  max resid 0.0003029515 
    ## ... Similar to previous best
    ## Run 116 stress 9.838369e-05 
    ## ... Procrustes: rmse 0.0002213791  max resid 0.0004374289 
    ## ... Similar to previous best
    ## Run 117 stress 9.900301e-05 
    ## ... Procrustes: rmse 0.0002111013  max resid 0.000424729 
    ## ... Similar to previous best
    ## Run 118 stress 9.960392e-05 
    ## ... Procrustes: rmse 0.0002116621  max resid 0.0004409891 
    ## ... Similar to previous best
    ## Run 119 stress 9.176499e-05 
    ## ... Procrustes: rmse 3.844174e-05  max resid 7.002736e-05 
    ## ... Similar to previous best
    ## Run 120 stress 9.645965e-05 
    ## ... Procrustes: rmse 0.0001220523  max resid 0.0002622057 
    ## ... Similar to previous best
    ## Run 121 stress 9.762052e-05 
    ## ... Procrustes: rmse 0.000104965  max resid 0.0002029153 
    ## ... Similar to previous best
    ## Run 122 stress 9.154212e-05 
    ## ... Procrustes: rmse 0.0002044179  max resid 0.0003175228 
    ## ... Similar to previous best
    ## Run 123 stress 9.716935e-05 
    ## ... Procrustes: rmse 0.0001749839  max resid 0.0003579652 
    ## ... Similar to previous best
    ## Run 124 stress 9.764521e-05 
    ## ... Procrustes: rmse 0.0002382007  max resid 0.0004664901 
    ## ... Similar to previous best
    ## Run 125 stress 9.792074e-05 
    ## ... Procrustes: rmse 0.0001034915  max resid 0.0001981549 
    ## ... Similar to previous best
    ## Run 126 stress 9.28386e-05 
    ## ... Procrustes: rmse 0.0001871767  max resid 0.0003799618 
    ## ... Similar to previous best
    ## Run 127 stress 8.880668e-05 
    ## ... Procrustes: rmse 0.0001882884  max resid 0.0003847727 
    ## ... Similar to previous best
    ## Run 128 stress 9.777062e-05 
    ## ... Procrustes: rmse 0.0001679502  max resid 0.0003238596 
    ## ... Similar to previous best
    ## Run 129 stress 9.92972e-05 
    ## ... Procrustes: rmse 5.206814e-05  max resid 7.533311e-05 
    ## ... Similar to previous best
    ## Run 130 stress 9.5245e-05 
    ## ... Procrustes: rmse 2.776272e-05  max resid 5.270509e-05 
    ## ... Similar to previous best
    ## Run 131 stress 9.879832e-05 
    ## ... Procrustes: rmse 0.0001866672  max resid 0.0004233817 
    ## ... Similar to previous best
    ## Run 132 stress 9.914959e-05 
    ## ... Procrustes: rmse 0.0002237243  max resid 0.0004551426 
    ## ... Similar to previous best
    ## Run 133 stress 9.595847e-05 
    ## ... Procrustes: rmse 4.594198e-05  max resid 6.996799e-05 
    ## ... Similar to previous best
    ## Run 134 stress 9.642341e-05 
    ## ... Procrustes: rmse 0.0001632014  max resid 0.0003118426 
    ## ... Similar to previous best
    ## Run 135 stress 9.425083e-05 
    ## ... Procrustes: rmse 0.000135669  max resid 0.0003179059 
    ## ... Similar to previous best
    ## Run 136 stress 9.308536e-05 
    ## ... Procrustes: rmse 9.672089e-05  max resid 0.0001960697 
    ## ... Similar to previous best
    ## Run 137 stress 9.735603e-05 
    ## ... Procrustes: rmse 0.0001410657  max resid 0.0003140554 
    ## ... Similar to previous best
    ## Run 138 stress 9.485879e-05 
    ## ... Procrustes: rmse 0.0001659353  max resid 0.0002839718 
    ## ... Similar to previous best
    ## Run 139 stress 9.437594e-05 
    ## ... Procrustes: rmse 0.0002029042  max resid 0.000413601 
    ## ... Similar to previous best
    ## Run 140 stress 9.626446e-05 
    ## ... Procrustes: rmse 0.0001892311  max resid 0.0003894865 
    ## ... Similar to previous best
    ## Run 141 stress 8.826663e-05 
    ## ... Procrustes: rmse 0.0001522576  max resid 0.0003054467 
    ## ... Similar to previous best
    ## Run 142 stress 9.731849e-05 
    ## ... Procrustes: rmse 0.0002091613  max resid 0.0004261847 
    ## ... Similar to previous best
    ## Run 143 stress 9.339227e-05 
    ## ... Procrustes: rmse 5.295427e-05  max resid 0.0001270588 
    ## ... Similar to previous best
    ## Run 144 stress 9.851824e-05 
    ## ... Procrustes: rmse 0.0002102928  max resid 0.0003420866 
    ## ... Similar to previous best
    ## Run 145 stress 9.660109e-05 
    ## ... Procrustes: rmse 0.0001633369  max resid 0.0003203182 
    ## ... Similar to previous best
    ## Run 146 stress 9.463415e-05 
    ## ... Procrustes: rmse 0.0002211242  max resid 0.0003523789 
    ## ... Similar to previous best
    ## Run 147 stress 8.876347e-05 
    ## ... Procrustes: rmse 8.427431e-05  max resid 0.0001836989 
    ## ... Similar to previous best
    ## Run 148 stress 8.524144e-05 
    ## ... Procrustes: rmse 0.0001635404  max resid 0.0003632661 
    ## ... Similar to previous best
    ## Run 149 stress 9.714494e-05 
    ## ... Procrustes: rmse 0.0001801182  max resid 0.0004010164 
    ## ... Similar to previous best
    ## Run 150 stress 9.368642e-05 
    ## ... Procrustes: rmse 0.0001774514  max resid 0.0004075682 
    ## ... Similar to previous best
    ## Run 151 stress 9.974434e-05 
    ## ... Procrustes: rmse 0.0001712186  max resid 0.0003285526 
    ## ... Similar to previous best
    ## Run 152 stress 9.145683e-05 
    ## ... Procrustes: rmse 0.0001100809  max resid 0.0002843087 
    ## ... Similar to previous best
    ## Run 153 stress 9.670257e-05 
    ## ... Procrustes: rmse 0.0001789535  max resid 0.0004238297 
    ## ... Similar to previous best
    ## Run 154 stress 9.078669e-05 
    ## ... Procrustes: rmse 9.052902e-05  max resid 0.0001888065 
    ## ... Similar to previous best
    ## Run 155 stress 0.206342 
    ## Run 156 stress 9.472004e-05 
    ## ... Procrustes: rmse 0.0001092832  max resid 0.0002148451 
    ## ... Similar to previous best
    ## Run 157 stress 9.744778e-05 
    ## ... Procrustes: rmse 0.0001686641  max resid 0.0003259288 
    ## ... Similar to previous best
    ## Run 158 stress 8.688148e-05 
    ## ... Procrustes: rmse 0.0001491863  max resid 0.0003017982 
    ## ... Similar to previous best
    ## Run 159 stress 9.25509e-05 
    ## ... Procrustes: rmse 0.0001593875  max resid 0.0003147253 
    ## ... Similar to previous best
    ## Run 160 stress 9.8363e-05 
    ## ... Procrustes: rmse 0.000215275  max resid 0.0004451072 
    ## ... Similar to previous best
    ## Run 161 stress 9.575794e-05 
    ## ... Procrustes: rmse 0.0001008772  max resid 0.0001993372 
    ## ... Similar to previous best
    ## Run 162 stress 7.119588e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 7.235209e-05  max resid 0.0001778423 
    ## ... Similar to previous best
    ## Run 163 stress 9.046942e-05 
    ## ... Procrustes: rmse 0.0001327306  max resid 0.0002578637 
    ## ... Similar to previous best
    ## Run 164 stress 9.543049e-05 
    ## ... Procrustes: rmse 0.0001406106  max resid 0.0002295944 
    ## ... Similar to previous best
    ## Run 165 stress 9.495635e-05 
    ## ... Procrustes: rmse 0.000157437  max resid 0.0003014262 
    ## ... Similar to previous best
    ## Run 166 stress 0.2595628 
    ## Run 167 stress 8.338543e-05 
    ## ... Procrustes: rmse 3.415938e-05  max resid 7.273277e-05 
    ## ... Similar to previous best
    ## Run 168 stress 9.056839e-05 
    ## ... Procrustes: rmse 0.0001577277  max resid 0.0003665734 
    ## ... Similar to previous best
    ## Run 169 stress 9.831864e-05 
    ## ... Procrustes: rmse 7.609076e-05  max resid 0.0001283565 
    ## ... Similar to previous best
    ## Run 170 stress 8.536567e-05 
    ## ... Procrustes: rmse 9.679598e-05  max resid 0.0001602558 
    ## ... Similar to previous best
    ## Run 171 stress 9.868084e-05 
    ## ... Procrustes: rmse 0.0001336041  max resid 0.0002945868 
    ## ... Similar to previous best
    ## Run 172 stress 9.683843e-05 
    ## ... Procrustes: rmse 0.0001496991  max resid 0.000285226 
    ## ... Similar to previous best
    ## Run 173 stress 9.599305e-05 
    ## ... Procrustes: rmse 0.0001842917  max resid 0.000366218 
    ## ... Similar to previous best
    ## Run 174 stress 8.192809e-05 
    ## ... Procrustes: rmse 0.0001636844  max resid 0.0003340757 
    ## ... Similar to previous best
    ## Run 175 stress 7.87627e-05 
    ## ... Procrustes: rmse 0.0001187103  max resid 0.0001800095 
    ## ... Similar to previous best
    ## Run 176 stress 9.797562e-05 
    ## ... Procrustes: rmse 0.0001597851  max resid 0.0003036208 
    ## ... Similar to previous best
    ## Run 177 stress 8.76152e-05 
    ## ... Procrustes: rmse 0.0001902059  max resid 0.0003888222 
    ## ... Similar to previous best
    ## Run 178 stress 9.811669e-05 
    ## ... Procrustes: rmse 0.0001743175  max resid 0.0002619449 
    ## ... Similar to previous best
    ## Run 179 stress 9.442858e-05 
    ## ... Procrustes: rmse 0.0001889038  max resid 0.0003867171 
    ## ... Similar to previous best
    ## Run 180 stress 9.561395e-05 
    ## ... Procrustes: rmse 0.0001676311  max resid 0.0003771431 
    ## ... Similar to previous best
    ## Run 181 stress 9.892796e-05 
    ## ... Procrustes: rmse 0.0001523527  max resid 0.0002742626 
    ## ... Similar to previous best
    ## Run 182 stress 9.837209e-05 
    ## ... Procrustes: rmse 0.0001736189  max resid 0.000401771 
    ## ... Similar to previous best
    ## Run 183 stress 9.599192e-05 
    ## ... Procrustes: rmse 0.0002074565  max resid 0.0004130584 
    ## ... Similar to previous best
    ## Run 184 stress 9.430105e-05 
    ## ... Procrustes: rmse 0.0001984074  max resid 0.0004005124 
    ## ... Similar to previous best
    ## Run 185 stress 7.344955e-05 
    ## ... Procrustes: rmse 9.9606e-05  max resid 0.0002400966 
    ## ... Similar to previous best
    ## Run 186 stress 9.774353e-05 
    ## ... Procrustes: rmse 0.000210595  max resid 0.0004179826 
    ## ... Similar to previous best
    ## Run 187 stress 9.409843e-05 
    ## ... Procrustes: rmse 0.0001856015  max resid 0.0003107761 
    ## ... Similar to previous best
    ## Run 188 stress 9.354183e-05 
    ## ... Procrustes: rmse 0.0001627905  max resid 0.0002476157 
    ## ... Similar to previous best
    ## Run 189 stress 8.474862e-05 
    ## ... Procrustes: rmse 9.853702e-05  max resid 0.0001945575 
    ## ... Similar to previous best
    ## Run 190 stress 9.891172e-05 
    ## ... Procrustes: rmse 0.0001705899  max resid 0.0003790864 
    ## ... Similar to previous best
    ## Run 191 stress 9.159864e-05 
    ## ... Procrustes: rmse 0.0001784462  max resid 0.0003042766 
    ## ... Similar to previous best
    ## Run 192 stress 9.604682e-05 
    ## ... Procrustes: rmse 0.0001582765  max resid 0.0003016717 
    ## ... Similar to previous best
    ## Run 193 stress 9.841044e-05 
    ## ... Procrustes: rmse 0.0001929111  max resid 0.0003185415 
    ## ... Similar to previous best
    ## Run 194 stress 9.41661e-05 
    ## ... Procrustes: rmse 0.0001491917  max resid 0.0002833765 
    ## ... Similar to previous best
    ## Run 195 stress 9.521126e-05 
    ## ... Procrustes: rmse 6.813594e-05  max resid 0.0001168662 
    ## ... Similar to previous best
    ## Run 196 stress 9.308981e-05 
    ## ... Procrustes: rmse 0.0001665157  max resid 0.0003128142 
    ## ... Similar to previous best
    ## Run 197 stress 9.54579e-05 
    ## ... Procrustes: rmse 0.0001916136  max resid 0.0003664189 
    ## ... Similar to previous best
    ## Run 198 stress 9.899251e-05 
    ## ... Procrustes: rmse 5.691021e-05  max resid 9.062552e-05 
    ## ... Similar to previous best
    ## Run 199 stress 7.427132e-05 
    ## ... Procrustes: rmse 0.0001327116  max resid 0.0002450704 
    ## ... Similar to previous best
    ## Run 200 stress 7.769209e-05 
    ## ... Procrustes: rmse 0.0001430982  max resid 0.0002917381 
    ## ... Similar to previous best
    ## Run 201 stress 9.736128e-05 
    ## ... Procrustes: rmse 0.0001941816  max resid 0.0003212616 
    ## ... Similar to previous best
    ## Run 202 stress 9.886332e-05 
    ## ... Procrustes: rmse 0.0001432628  max resid 0.0002826745 
    ## ... Similar to previous best
    ## Run 203 stress 9.742617e-05 
    ## ... Procrustes: rmse 0.0001535859  max resid 0.0002831164 
    ## ... Similar to previous best
    ## Run 204 stress 9.954556e-05 
    ## ... Procrustes: rmse 0.0001455104  max resid 0.000307275 
    ## ... Similar to previous best
    ## Run 205 stress 9.666039e-05 
    ## ... Procrustes: rmse 0.0001335507  max resid 0.0002928673 
    ## ... Similar to previous best
    ## Run 206 stress 9.089449e-05 
    ## ... Procrustes: rmse 0.000160488  max resid 0.0003814566 
    ## ... Similar to previous best
    ## Run 207 stress 9.874132e-05 
    ## ... Procrustes: rmse 0.0001457007  max resid 0.0002665798 
    ## ... Similar to previous best
    ## Run 208 stress 8.924996e-05 
    ## ... Procrustes: rmse 0.0001577671  max resid 0.0003798888 
    ## ... Similar to previous best
    ## Run 209 stress 9.752717e-05 
    ## ... Procrustes: rmse 0.0001728927  max resid 0.0002598061 
    ## ... Similar to previous best
    ## Run 210 stress 9.678719e-05 
    ## ... Procrustes: rmse 0.0002005215  max resid 0.0003166153 
    ## ... Similar to previous best
    ## Run 211 stress 9.543908e-05 
    ## ... Procrustes: rmse 0.0001677748  max resid 0.0003755885 
    ## ... Similar to previous best
    ## Run 212 stress 8.881013e-05 
    ## ... Procrustes: rmse 0.0001196328  max resid 0.0001831475 
    ## ... Similar to previous best
    ## Run 213 stress 9.573495e-05 
    ## ... Procrustes: rmse 0.0001915108  max resid 0.0003165165 
    ## ... Similar to previous best
    ## Run 214 stress 9.633603e-05 
    ## ... Procrustes: rmse 0.000195787  max resid 0.0003201169 
    ## ... Similar to previous best
    ## Run 215 stress 0.2480349 
    ## Run 216 stress 9.495838e-05 
    ## ... Procrustes: rmse 0.0001634169  max resid 0.0003123138 
    ## ... Similar to previous best
    ## Run 217 stress 9.707791e-05 
    ## ... Procrustes: rmse 0.0001749643  max resid 0.0003627243 
    ## ... Similar to previous best
    ## Run 218 stress 7.253361e-05 
    ## ... Procrustes: rmse 0.0001002155  max resid 0.0002370225 
    ## ... Similar to previous best
    ## Run 219 stress 9.480185e-05 
    ## ... Procrustes: rmse 0.0001478639  max resid 0.0002730058 
    ## ... Similar to previous best
    ## Run 220 stress 9.818383e-05 
    ## ... Procrustes: rmse 0.0002120333  max resid 0.000418632 
    ## ... Similar to previous best
    ## Run 221 stress 9.734893e-05 
    ## ... Procrustes: rmse 0.0001281831  max resid 0.0002411293 
    ## ... Similar to previous best
    ## Run 222 stress 9.214942e-05 
    ## ... Procrustes: rmse 0.0001281268  max resid 0.0003013493 
    ## ... Similar to previous best
    ## Run 223 stress 9.739229e-05 
    ## ... Procrustes: rmse 0.0001705918  max resid 0.0003823909 
    ## ... Similar to previous best
    ## Run 224 stress 7.776825e-05 
    ## ... Procrustes: rmse 0.0001229917  max resid 0.0002540912 
    ## ... Similar to previous best
    ## Run 225 stress 9.76314e-05 
    ## ... Procrustes: rmse 0.0002012996  max resid 0.0003170535 
    ## ... Similar to previous best
    ## Run 226 stress 9.290046e-05 
    ## ... Procrustes: rmse 0.0001403335  max resid 0.0002738166 
    ## ... Similar to previous best
    ## Run 227 stress 9.762542e-05 
    ## ... Procrustes: rmse 0.0002050238  max resid 0.0003205228 
    ## ... Similar to previous best
    ## Run 228 stress 9.404982e-05 
    ## ... Procrustes: rmse 0.0001890533  max resid 0.0003842609 
    ## ... Similar to previous best
    ## Run 229 stress 9.554569e-05 
    ## ... Procrustes: rmse 0.0001883513  max resid 0.000313722 
    ## ... Similar to previous best
    ## Run 230 stress 9.363216e-05 
    ## ... Procrustes: rmse 0.0001659533  max resid 0.000376609 
    ## ... Similar to previous best
    ## Run 231 stress 9.542216e-05 
    ## ... Procrustes: rmse 0.0001906966  max resid 0.0003155143 
    ## ... Similar to previous best
    ## Run 232 stress 9.88375e-05 
    ## ... Procrustes: rmse 0.0001538978  max resid 0.0002811011 
    ## ... Similar to previous best
    ## Run 233 stress 9.480699e-05 
    ## ... Procrustes: rmse 0.0002007195  max resid 0.0004037051 
    ## ... Similar to previous best
    ## Run 234 stress 9.132987e-05 
    ## ... Procrustes: rmse 0.0001601366  max resid 0.0003268907 
    ## ... Similar to previous best
    ## Run 235 stress 9.36903e-05 
    ## ... Procrustes: rmse 9.848801e-05  max resid 0.0002081917 
    ## ... Similar to previous best
    ## Run 236 stress 9.504178e-05 
    ## ... Procrustes: rmse 0.0001648893  max resid 0.0003106746 
    ## ... Similar to previous best
    ## Run 237 stress 9.973741e-05 
    ## ... Procrustes: rmse 0.0001729563  max resid 0.0003881302 
    ## ... Similar to previous best
    ## Run 238 stress 9.976704e-05 
    ## ... Procrustes: rmse 0.0001980561  max resid 0.0003236562 
    ## ... Similar to previous best
    ## Run 239 stress 8.910409e-05 
    ## ... Procrustes: rmse 0.0001506396  max resid 0.0003221832 
    ## ... Similar to previous best
    ## Run 240 stress 8.460454e-05 
    ## ... Procrustes: rmse 0.0001383577  max resid 0.0003028564 
    ## ... Similar to previous best
    ## Run 241 stress 9.707959e-05 
    ## ... Procrustes: rmse 0.0001503553  max resid 0.0002862227 
    ## ... Similar to previous best
    ## Run 242 stress 9.330251e-05 
    ## ... Procrustes: rmse 0.0001497369  max resid 0.0002757393 
    ## ... Similar to previous best
    ## Run 243 stress 9.839121e-05 
    ## ... Procrustes: rmse 0.0001974366  max resid 0.0003220222 
    ## ... Similar to previous best
    ## Run 244 stress 9.696445e-05 
    ## ... Procrustes: rmse 7.063902e-05  max resid 0.0001211196 
    ## ... Similar to previous best
    ## Run 245 stress 9.434289e-05 
    ## ... Procrustes: rmse 0.0001377084  max resid 0.0002272577 
    ## ... Similar to previous best
    ## Run 246 stress 9.018066e-05 
    ## ... Procrustes: rmse 0.0001624871  max resid 0.0002852687 
    ## ... Similar to previous best
    ## Run 247 stress 9.124406e-05 
    ## ... Procrustes: rmse 0.000183192  max resid 0.0003743716 
    ## ... Similar to previous best
    ## Run 248 stress 9.623505e-05 
    ## ... Procrustes: rmse 0.0001423219  max resid 0.000230569 
    ## ... Similar to previous best
    ## Run 249 stress 9.605323e-05 
    ## ... Procrustes: rmse 0.0001924651  max resid 0.0003896626 
    ## ... Similar to previous best
    ## Run 250 stress 9.325017e-05 
    ## ... Procrustes: rmse 0.0001587244  max resid 0.0003208662 
    ## ... Similar to previous best
    ## Run 251 stress 9.411841e-05 
    ## ... Procrustes: rmse 6.691071e-05  max resid 0.0001133962 
    ## ... Similar to previous best
    ## Run 252 stress 9.293682e-05 
    ## ... Procrustes: rmse 0.0001364965  max resid 0.000238713 
    ## ... Similar to previous best
    ## Run 253 stress 9.476729e-05 
    ## ... Procrustes: rmse 0.000144543  max resid 0.0002792666 
    ## ... Similar to previous best
    ## Run 254 stress 9.967563e-05 
    ## ... Procrustes: rmse 0.0002097898  max resid 0.0003257617 
    ## ... Similar to previous best
    ## Run 255 stress 9.747423e-05 
    ## ... Procrustes: rmse 0.0001945379  max resid 0.000379738 
    ## ... Similar to previous best
    ## Run 256 stress 9.343972e-05 
    ## ... Procrustes: rmse 8.953703e-05  max resid 0.000179619 
    ## ... Similar to previous best
    ## Run 257 stress 9.938983e-05 
    ## ... Procrustes: rmse 7.777385e-05  max resid 0.0001318873 
    ## ... Similar to previous best
    ## Run 258 stress 9.73241e-05 
    ## ... Procrustes: rmse 0.0001826645  max resid 0.0003624929 
    ## ... Similar to previous best
    ## Run 259 stress 9.95538e-05 
    ## ... Procrustes: rmse 7.027875e-05  max resid 0.0001199997 
    ## ... Similar to previous best
    ## Run 260 stress 8.830601e-05 
    ## ... Procrustes: rmse 0.0001850244  max resid 0.0003008834 
    ## ... Similar to previous best
    ## Run 261 stress 9.894167e-05 
    ## ... Procrustes: rmse 0.0001758535  max resid 0.0002643512 
    ## ... Similar to previous best
    ## Run 262 stress 9.063414e-05 
    ## ... Procrustes: rmse 8.021532e-05  max resid 0.0001884621 
    ## ... Similar to previous best
    ## Run 263 stress 9.413282e-05 
    ## ... Procrustes: rmse 0.0001343415  max resid 0.0002939022 
    ## ... Similar to previous best
    ## Run 264 stress 9.333193e-05 
    ## ... Procrustes: rmse 0.0001556026  max resid 0.0003078765 
    ## ... Similar to previous best
    ## Run 265 stress 9.786231e-05 
    ## ... Procrustes: rmse 0.0001612242  max resid 0.0002405308 
    ## ... Similar to previous best
    ## Run 266 stress 9.798507e-05 
    ## ... Procrustes: rmse 0.000154883  max resid 0.0002855084 
    ## ... Similar to previous best
    ## Run 267 stress 8.905742e-05 
    ## ... Procrustes: rmse 0.0001920095  max resid 0.0003911468 
    ## ... Similar to previous best
    ## Run 268 stress 9.722541e-05 
    ## ... Procrustes: rmse 9.504397e-05  max resid 0.0002106685 
    ## ... Similar to previous best
    ## Run 269 stress 9.504325e-05 
    ## ... Procrustes: rmse 0.0001505746  max resid 0.0002784484 
    ## ... Similar to previous best
    ## Run 270 stress 9.482687e-05 
    ## ... Procrustes: rmse 0.000144338  max resid 0.0002794796 
    ## ... Similar to previous best
    ## Run 271 stress 9.391858e-05 
    ## ... Procrustes: rmse 0.0001350695  max resid 0.0002952593 
    ## ... Similar to previous best
    ## Run 272 stress 9.831633e-05 
    ## ... Procrustes: rmse 0.0001301349  max resid 0.0002888859 
    ## ... Similar to previous best
    ## Run 273 stress 9.385174e-05 
    ## ... Procrustes: rmse 9.772881e-05  max resid 0.0002173692 
    ## ... Similar to previous best
    ## Run 274 stress 9.539168e-05 
    ## ... Procrustes: rmse 0.0001803728  max resid 0.000357362 
    ## ... Similar to previous best
    ## Run 275 stress 9.87925e-05 
    ## ... Procrustes: rmse 0.0001690373  max resid 0.0003635096 
    ## ... Similar to previous best
    ## Run 276 stress 9.807974e-05 
    ## ... Procrustes: rmse 0.0001560419  max resid 0.0002974278 
    ## ... Similar to previous best
    ## Run 277 stress 9.644425e-05 
    ## ... Procrustes: rmse 0.0001712151  max resid 0.0003542461 
    ## ... Similar to previous best
    ## Run 278 stress 9.913011e-05 
    ## ... Procrustes: rmse 0.0002049317  max resid 0.0004095846 
    ## ... Similar to previous best
    ## Run 279 stress 9.482757e-05 
    ## ... Procrustes: rmse 0.0001369914  max resid 0.0002979548 
    ## ... Similar to previous best
    ## Run 280 stress 8.740587e-05 
    ## ... Procrustes: rmse 4.55905e-05  max resid 7.487133e-05 
    ## ... Similar to previous best
    ## Run 281 stress 8.990284e-05 
    ## ... Procrustes: rmse 0.0001445125  max resid 0.0002295388 
    ## ... Similar to previous best
    ## Run 282 stress 9.543185e-05 
    ## ... Procrustes: rmse 0.0001969904  max resid 0.0003122933 
    ## ... Similar to previous best
    ## Run 283 stress 9.899606e-05 
    ## ... Procrustes: rmse 0.0001449139  max resid 0.00030633 
    ## ... Similar to previous best
    ## Run 284 stress 8.661975e-05 
    ## ... Procrustes: rmse 0.0001512603  max resid 0.0003309038 
    ## ... Similar to previous best
    ## Run 285 stress 9.815577e-05 
    ## ... Procrustes: rmse 0.0001422865  max resid 0.0003027165 
    ## ... Similar to previous best
    ## Run 286 stress 9.212892e-05 
    ## ... Procrustes: rmse 6.41515e-05  max resid 0.0001193567 
    ## ... Similar to previous best
    ## Run 287 stress 9.873402e-05 
    ## ... Procrustes: rmse 0.0001616785  max resid 0.0003060736 
    ## ... Similar to previous best
    ## Run 288 stress 9.902018e-05 
    ## ... Procrustes: rmse 0.0001445125  max resid 0.0002818721 
    ## ... Similar to previous best
    ## Run 289 stress 9.846878e-05 
    ## ... Procrustes: rmse 0.0001152169  max resid 0.0001973053 
    ## ... Similar to previous best
    ## Run 290 stress 9.311028e-05 
    ## ... Procrustes: rmse 0.0001615597  max resid 0.0003713612 
    ## ... Similar to previous best
    ## Run 291 stress 9.606438e-05 
    ## ... Procrustes: rmse 0.0001798496  max resid 0.0003601826 
    ## ... Similar to previous best
    ## Run 292 stress 9.978007e-05 
    ## ... Procrustes: rmse 0.000168661  max resid 0.0003574107 
    ## ... Similar to previous best
    ## Run 293 stress 9.844056e-05 
    ## ... Procrustes: rmse 0.0001425312  max resid 0.0003037945 
    ## ... Similar to previous best
    ## Run 294 stress 9.480519e-05 
    ## ... Procrustes: rmse 0.0001682775  max resid 0.0003617631 
    ## ... Similar to previous best
    ## Run 295 stress 9.729902e-05 
    ## ... Procrustes: rmse 0.0002097776  max resid 0.0004162611 
    ## ... Similar to previous best
    ## Run 296 stress 9.576715e-05 
    ## ... Procrustes: rmse 0.0001694401  max resid 0.0002554903 
    ## ... Similar to previous best
    ## Run 297 stress 9.867991e-05 
    ## ... Procrustes: rmse 0.0001974395  max resid 0.000395795 
    ## ... Similar to previous best
    ## Run 298 stress 9.856681e-05 
    ## ... Procrustes: rmse 0.0001890675  max resid 0.0003138155 
    ## ... Similar to previous best
    ## Run 299 stress 9.506311e-05 
    ## ... Procrustes: rmse 0.0001452169  max resid 0.0002759788 
    ## ... Similar to previous best
    ## Run 300 stress 8.580189e-05 
    ## ... Procrustes: rmse 0.0001321989  max resid 0.0002609464 
    ## ... Similar to previous best
    ## Run 301 stress 9.40507e-05 
    ## ... Procrustes: rmse 9.769372e-05  max resid 0.0002163575 
    ## ... Similar to previous best
    ## Run 302 stress 9.635889e-05 
    ## ... Procrustes: rmse 0.0001684079  max resid 0.0003799526 
    ## ... Similar to previous best
    ## Run 303 stress 9.776789e-05 
    ## ... Procrustes: rmse 5.856343e-05  max resid 0.0001174397 
    ## ... Similar to previous best
    ## Run 304 stress 8.830141e-05 
    ## ... Procrustes: rmse 0.0001412849  max resid 0.0002972612 
    ## ... Similar to previous best
    ## Run 305 stress 9.547571e-05 
    ## ... Procrustes: rmse 0.0001524975  max resid 0.0002883846 
    ## ... Similar to previous best
    ## Run 306 stress 8.663233e-05 
    ## ... Procrustes: rmse 0.0001234229  max resid 0.000284239 
    ## ... Similar to previous best
    ## Run 307 stress 9.777532e-05 
    ## ... Procrustes: rmse 0.0001731098  max resid 0.0003862348 
    ## ... Similar to previous best
    ## Run 308 stress 9.625346e-05 
    ## ... Procrustes: rmse 0.0001801274  max resid 0.0003653015 
    ## ... Similar to previous best
    ## Run 309 stress 9.968287e-05 
    ## ... Procrustes: rmse 0.0001990764  max resid 0.0003482193 
    ## ... Similar to previous best
    ## Run 310 stress 9.442487e-05 
    ## ... Procrustes: rmse 0.0001677074  max resid 0.0002943492 
    ## ... Similar to previous best
    ## Run 311 stress 5.60075e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001097997  max resid 0.0001433972 
    ## ... Similar to previous best
    ## Run 312 stress 9.506214e-05 
    ## ... Procrustes: rmse 0.0001843722  max resid 0.0002622783 
    ## ... Similar to previous best
    ## Run 313 stress 9.571512e-05 
    ## ... Procrustes: rmse 0.0001854266  max resid 0.0002737561 
    ## ... Similar to previous best
    ## Run 314 stress 9.47352e-05 
    ## ... Procrustes: rmse 0.0001936215  max resid 0.0002818678 
    ## ... Similar to previous best
    ## Run 315 stress 9.221629e-05 
    ## ... Procrustes: rmse 0.0001716189  max resid 0.0002689059 
    ## ... Similar to previous best
    ## Run 316 stress 9.671329e-05 
    ## ... Procrustes: rmse 0.0001514126  max resid 0.0002263258 
    ## ... Similar to previous best
    ## Run 317 stress 9.122795e-05 
    ## ... Procrustes: rmse 0.0001639891  max resid 0.0002343082 
    ## ... Similar to previous best
    ## Run 318 stress 9.763074e-05 
    ## ... Procrustes: rmse 0.0001916525  max resid 0.0002680807 
    ## ... Similar to previous best
    ## Run 319 stress 9.856003e-05 
    ## ... Procrustes: rmse 0.0001626777  max resid 0.0002435124 
    ## ... Similar to previous best
    ## Run 320 stress 9.072042e-05 
    ## ... Procrustes: rmse 0.0001182115  max resid 0.000193593 
    ## ... Similar to previous best
    ## Run 321 stress 7.994264e-05 
    ## ... Procrustes: rmse 0.0001416794  max resid 0.0001966403 
    ## ... Similar to previous best
    ## Run 322 stress 9.591657e-05 
    ## ... Procrustes: rmse 0.0001785072  max resid 0.0002620044 
    ## ... Similar to previous best
    ## Run 323 stress 9.625937e-05 
    ## ... Procrustes: rmse 0.0001753633  max resid 0.0002644883 
    ## ... Similar to previous best
    ## Run 324 stress 9.713384e-05 
    ## ... Procrustes: rmse 0.0001952794  max resid 0.0002633587 
    ## ... Similar to previous best
    ## Run 325 stress 9.404321e-05 
    ## ... Procrustes: rmse 0.0001794411  max resid 0.0002837894 
    ## ... Similar to previous best
    ## Run 326 stress 9.173332e-05 
    ## ... Procrustes: rmse 0.0001204616  max resid 0.0001997937 
    ## ... Similar to previous best
    ## Run 327 stress 9.383729e-05 
    ## ... Procrustes: rmse 0.0001786412  max resid 0.000278834 
    ## ... Similar to previous best
    ## Run 328 stress 8.485755e-05 
    ## ... Procrustes: rmse 0.0001294848  max resid 0.0001824835 
    ## ... Similar to previous best
    ## Run 329 stress 9.732841e-05 
    ## ... Procrustes: rmse 0.0002050839  max resid 0.0002917365 
    ## ... Similar to previous best
    ## Run 330 stress 9.530438e-05 
    ## ... Procrustes: rmse 0.0001834174  max resid 0.0002715398 
    ## ... Similar to previous best
    ## Run 331 stress 9.723656e-05 
    ## ... Procrustes: rmse 0.0002012548  max resid 0.0002946064 
    ## ... Similar to previous best
    ## Run 332 stress 9.372322e-05 
    ## ... Procrustes: rmse 0.0001235521  max resid 0.0001796174 
    ## ... Similar to previous best
    ## Run 333 stress 9.747957e-05 
    ## ... Procrustes: rmse 0.0001998374  max resid 0.0002819388 
    ## ... Similar to previous best
    ## Run 334 stress 9.218401e-05 
    ## ... Procrustes: rmse 0.0001484245  max resid 0.0002259451 
    ## ... Similar to previous best
    ## Run 335 stress 9.349547e-05 
    ## ... Procrustes: rmse 0.0002065347  max resid 0.0002751859 
    ## ... Similar to previous best
    ## Run 336 stress 9.530456e-05 
    ## ... Procrustes: rmse 0.0001496234  max resid 0.0002314894 
    ## ... Similar to previous best
    ## Run 337 stress 8.512083e-05 
    ## ... Procrustes: rmse 0.000154607  max resid 0.0002225906 
    ## ... Similar to previous best
    ## Run 338 stress 9.499455e-05 
    ## ... Procrustes: rmse 0.0001850511  max resid 0.0002725423 
    ## ... Similar to previous best
    ## Run 339 stress 9.99702e-05 
    ## ... Procrustes: rmse 0.0002201028  max resid 0.0003000303 
    ## ... Similar to previous best
    ## Run 340 stress 9.321704e-05 
    ## ... Procrustes: rmse 0.0002127789  max resid 0.0002970236 
    ## ... Similar to previous best
    ## Run 341 stress 9.746945e-05 
    ## ... Procrustes: rmse 0.0001914167  max resid 0.0002625057 
    ## ... Similar to previous best
    ## Run 342 stress 9.933702e-05 
    ## ... Procrustes: rmse 0.0001851555  max resid 0.0002669695 
    ## ... Similar to previous best
    ## Run 343 stress 9.549809e-05 
    ## ... Procrustes: rmse 0.0001524966  max resid 0.0002649857 
    ## ... Similar to previous best
    ## Run 344 stress 7.384904e-05 
    ## ... Procrustes: rmse 0.0001581126  max resid 0.0002092641 
    ## ... Similar to previous best
    ## Run 345 stress 7.681017e-05 
    ## ... Procrustes: rmse 0.0001277894  max resid 0.0001741058 
    ## ... Similar to previous best
    ## Run 346 stress 9.613963e-05 
    ## ... Procrustes: rmse 0.000101678  max resid 0.0001890514 
    ## ... Similar to previous best
    ## Run 347 stress 9.849452e-05 
    ## ... Procrustes: rmse 0.0002133589  max resid 0.0002985388 
    ## ... Similar to previous best
    ## Run 348 stress 9.932325e-05 
    ## ... Procrustes: rmse 0.0001932779  max resid 0.0002783541 
    ## ... Similar to previous best
    ## Run 349 stress 8.477847e-05 
    ## ... Procrustes: rmse 0.0001649897  max resid 0.0002200456 
    ## ... Similar to previous best
    ## Run 350 stress 9.404319e-05 
    ## ... Procrustes: rmse 0.0001668884  max resid 0.0002414598 
    ## ... Similar to previous best
    ## Run 351 stress 9.782338e-05 
    ## ... Procrustes: rmse 0.0001979332  max resid 0.0002994905 
    ## ... Similar to previous best
    ## Run 352 stress 9.360468e-05 
    ## ... Procrustes: rmse 0.0002037629  max resid 0.000300852 
    ## ... Similar to previous best
    ## Run 353 stress 9.894678e-05 
    ## ... Procrustes: rmse 0.0001979138  max resid 0.0002908522 
    ## ... Similar to previous best
    ## Run 354 stress 8.984839e-05 
    ## ... Procrustes: rmse 0.0001579115  max resid 0.0002331669 
    ## ... Similar to previous best
    ## Run 355 stress 8.684821e-05 
    ## ... Procrustes: rmse 0.0001715363  max resid 0.0002488544 
    ## ... Similar to previous best
    ## Run 356 stress 7.06486e-05 
    ## ... Procrustes: rmse 0.000120405  max resid 0.0001639435 
    ## ... Similar to previous best
    ## Run 357 stress 9.85584e-05 
    ## ... Procrustes: rmse 0.0001337549  max resid 0.0001745603 
    ## ... Similar to previous best
    ## Run 358 stress 6.844047e-05 
    ## ... Procrustes: rmse 0.0001141698  max resid 0.0001524452 
    ## ... Similar to previous best
    ## Run 359 stress 9.077445e-05 
    ## ... Procrustes: rmse 0.0001654973  max resid 0.0002239847 
    ## ... Similar to previous best
    ## Run 360 stress 9.57764e-05 
    ## ... Procrustes: rmse 0.0002123912  max resid 0.0003055365 
    ## ... Similar to previous best
    ## Run 361 stress 9.221343e-05 
    ## ... Procrustes: rmse 0.0001685015  max resid 0.0002506759 
    ## ... Similar to previous best
    ## Run 362 stress 9.437484e-05 
    ## ... Procrustes: rmse 0.0002079385  max resid 0.0002966194 
    ## ... Similar to previous best
    ## Run 363 stress 9.701149e-05 
    ## ... Procrustes: rmse 0.0001810894  max resid 0.0002744283 
    ## ... Similar to previous best
    ## Run 364 stress 9.712523e-05 
    ## ... Procrustes: rmse 0.0001512156  max resid 0.0002594657 
    ## ... Similar to previous best
    ## Run 365 stress 9.324724e-05 
    ## ... Procrustes: rmse 0.0001798749  max resid 0.0002557182 
    ## ... Similar to previous best
    ## Run 366 stress 8.95674e-05 
    ## ... Procrustes: rmse 0.0001933992  max resid 0.0002552248 
    ## ... Similar to previous best
    ## Run 367 stress 9.274032e-05 
    ## ... Procrustes: rmse 0.0001841871  max resid 0.0002422143 
    ## ... Similar to previous best
    ## Run 368 stress 9.046149e-05 
    ## ... Procrustes: rmse 0.0001434508  max resid 0.0002144307 
    ## ... Similar to previous best
    ## Run 369 stress 9.59434e-05 
    ## ... Procrustes: rmse 0.0002044759  max resid 0.0002935066 
    ## ... Similar to previous best
    ## Run 370 stress 8.475045e-05 
    ## ... Procrustes: rmse 0.0001108084  max resid 0.0001653786 
    ## ... Similar to previous best
    ## Run 371 stress 8.909251e-05 
    ## ... Procrustes: rmse 0.000121843  max resid 0.0002108708 
    ## ... Similar to previous best
    ## Run 372 stress 9.058076e-05 
    ## ... Procrustes: rmse 0.0001740181  max resid 0.0002412338 
    ## ... Similar to previous best
    ## Run 373 stress 9.940376e-05 
    ## ... Procrustes: rmse 0.0001948169  max resid 0.0002634384 
    ## ... Similar to previous best
    ## Run 374 stress 9.430602e-05 
    ## ... Procrustes: rmse 0.000188619  max resid 0.000270715 
    ## ... Similar to previous best
    ## Run 375 stress 9.770061e-05 
    ## ... Procrustes: rmse 0.0001596461  max resid 0.0002226682 
    ## ... Similar to previous best
    ## Run 376 stress 9.284031e-05 
    ## ... Procrustes: rmse 0.0001785425  max resid 0.0002428229 
    ## ... Similar to previous best
    ## Run 377 stress 9.502988e-05 
    ## ... Procrustes: rmse 0.0001730632  max resid 0.0002397931 
    ## ... Similar to previous best
    ## Run 378 stress 9.006957e-05 
    ## ... Procrustes: rmse 0.0001480136  max resid 0.0002338535 
    ## ... Similar to previous best
    ## Run 379 stress 9.852786e-05 
    ## ... Procrustes: rmse 0.0002116437  max resid 0.0003138276 
    ## ... Similar to previous best
    ## Run 380 stress 9.225941e-05 
    ## ... Procrustes: rmse 0.0001891032  max resid 0.0002602023 
    ## ... Similar to previous best
    ## Run 381 stress 9.398595e-05 
    ## ... Procrustes: rmse 0.0002024505  max resid 0.0002855187 
    ## ... Similar to previous best
    ## Run 382 stress 9.965003e-05 
    ## ... Procrustes: rmse 0.0001836096  max resid 0.0002639007 
    ## ... Similar to previous best
    ## Run 383 stress 9.207914e-05 
    ## ... Procrustes: rmse 0.0001419801  max resid 0.0002023112 
    ## ... Similar to previous best
    ## Run 384 stress 9.407984e-05 
    ## ... Procrustes: rmse 0.0001963999  max resid 0.0002712718 
    ## ... Similar to previous best
    ## Run 385 stress 9.244611e-05 
    ## ... Procrustes: rmse 0.0001668419  max resid 0.0002504651 
    ## ... Similar to previous best
    ## Run 386 stress 9.022962e-05 
    ## ... Procrustes: rmse 0.0001725474  max resid 0.0002286244 
    ## ... Similar to previous best
    ## Run 387 stress 9.365812e-05 
    ## ... Procrustes: rmse 0.0002042322  max resid 0.0002820159 
    ## ... Similar to previous best
    ## Run 388 stress 9.218999e-05 
    ## ... Procrustes: rmse 0.0001955288  max resid 0.0002771861 
    ## ... Similar to previous best
    ## Run 389 stress 9.372988e-05 
    ## ... Procrustes: rmse 0.0002009728  max resid 0.0002846607 
    ## ... Similar to previous best
    ## Run 390 stress 9.43541e-05 
    ## ... Procrustes: rmse 0.0001786681  max resid 0.0002444113 
    ## ... Similar to previous best
    ## Run 391 stress 9.895429e-05 
    ## ... Procrustes: rmse 0.0001799046  max resid 0.0002595103 
    ## ... Similar to previous best
    ## Run 392 stress 8.930387e-05 
    ## ... Procrustes: rmse 0.0001196484  max resid 0.0001926319 
    ## ... Similar to previous best
    ## Run 393 stress 9.89658e-05 
    ## ... Procrustes: rmse 0.0001831867  max resid 0.00023198 
    ## ... Similar to previous best
    ## Run 394 stress 9.339474e-05 
    ## ... Procrustes: rmse 0.0001921338  max resid 0.0002800944 
    ## ... Similar to previous best
    ## Run 395 stress 9.741504e-05 
    ## ... Procrustes: rmse 0.0001910215  max resid 0.000272499 
    ## ... Similar to previous best
    ## Run 396 stress 8.720506e-05 
    ## ... Procrustes: rmse 7.577157e-05  max resid 0.0001375195 
    ## ... Similar to previous best
    ## Run 397 stress 9.919469e-05 
    ## ... Procrustes: rmse 0.0001921338  max resid 0.0002686878 
    ## ... Similar to previous best
    ## Run 398 stress 9.5823e-05 
    ## ... Procrustes: rmse 0.0002071283  max resid 0.0002889715 
    ## ... Similar to previous best
    ## Run 399 stress 9.384919e-05 
    ## ... Procrustes: rmse 0.0001918895  max resid 0.0002535042 
    ## ... Similar to previous best
    ## Run 400 stress 9.151949e-05 
    ## ... Procrustes: rmse 0.0001994854  max resid 0.0002721311 
    ## ... Similar to previous best
    ## Run 401 stress 9.661367e-05 
    ## ... Procrustes: rmse 0.0001791707  max resid 0.0002517075 
    ## ... Similar to previous best
    ## Run 402 stress 9.564646e-05 
    ## ... Procrustes: rmse 0.0001992839  max resid 0.000299708 
    ## ... Similar to previous best
    ## Run 403 stress 9.457039e-05 
    ## ... Procrustes: rmse 0.0001680901  max resid 0.0002578914 
    ## ... Similar to previous best
    ## Run 404 stress 8.706326e-05 
    ## ... Procrustes: rmse 0.0001635532  max resid 0.000222394 
    ## ... Similar to previous best
    ## Run 405 stress 9.731333e-05 
    ## ... Procrustes: rmse 0.0002043747  max resid 0.0002881831 
    ## ... Similar to previous best
    ## Run 406 stress 9.82914e-05 
    ## ... Procrustes: rmse 0.0001913377  max resid 0.000269926 
    ## ... Similar to previous best
    ## Run 407 stress 9.466285e-05 
    ## ... Procrustes: rmse 0.0001984498  max resid 0.0002807202 
    ## ... Similar to previous best
    ## Run 408 stress 9.263509e-05 
    ## ... Procrustes: rmse 0.0001720962  max resid 0.0002661094 
    ## ... Similar to previous best
    ## Run 409 stress 9.755255e-05 
    ## ... Procrustes: rmse 0.000152224  max resid 0.0001903737 
    ## ... Similar to previous best
    ## Run 410 stress 8.901378e-05 
    ## ... Procrustes: rmse 0.0001891562  max resid 0.0002714223 
    ## ... Similar to previous best
    ## Run 411 stress 9.725509e-05 
    ## ... Procrustes: rmse 0.0002012049  max resid 0.0002952448 
    ## ... Similar to previous best
    ## Run 412 stress 9.461419e-05 
    ## ... Procrustes: rmse 0.0001825032  max resid 0.0002681078 
    ## ... Similar to previous best
    ## Run 413 stress 8.985605e-05 
    ## ... Procrustes: rmse 0.0001330925  max resid 0.0001899779 
    ## ... Similar to previous best
    ## Run 414 stress 9.552193e-05 
    ## ... Procrustes: rmse 0.0002021704  max resid 0.0002925026 
    ## ... Similar to previous best
    ## Run 415 stress 9.458174e-05 
    ## ... Procrustes: rmse 0.000204834  max resid 0.0003060267 
    ## ... Similar to previous best
    ## Run 416 stress 9.419649e-05 
    ## ... Procrustes: rmse 0.0001946648  max resid 0.0002697563 
    ## ... Similar to previous best
    ## Run 417 stress 9.627719e-05 
    ## ... Procrustes: rmse 0.0001858785  max resid 0.0002610249 
    ## ... Similar to previous best
    ## Run 418 stress 9.792851e-05 
    ## ... Procrustes: rmse 0.0001876121  max resid 0.0002817454 
    ## ... Similar to previous best
    ## Run 419 stress 9.898135e-05 
    ## ... Procrustes: rmse 0.0002203217  max resid 0.0002975514 
    ## ... Similar to previous best
    ## Run 420 stress 8.910241e-05 
    ## ... Procrustes: rmse 0.0002011012  max resid 0.0002752906 
    ## ... Similar to previous best
    ## Run 421 stress 9.712538e-05 
    ## ... Procrustes: rmse 0.0001890987  max resid 0.0002675767 
    ## ... Similar to previous best
    ## Run 422 stress 9.524322e-05 
    ## ... Procrustes: rmse 0.0001906487  max resid 0.0002910147 
    ## ... Similar to previous best
    ## Run 423 stress 9.84054e-05 
    ## ... Procrustes: rmse 0.0002176538  max resid 0.0002953019 
    ## ... Similar to previous best
    ## Run 424 stress 9.729634e-05 
    ## ... Procrustes: rmse 0.0001953541  max resid 0.0002634386 
    ## ... Similar to previous best
    ## Run 425 stress 9.834635e-05 
    ## ... Procrustes: rmse 0.0001777834  max resid 0.0002439494 
    ## ... Similar to previous best
    ## Run 426 stress 9.961889e-05 
    ## ... Procrustes: rmse 0.000185572  max resid 0.000273547 
    ## ... Similar to previous best
    ## Run 427 stress 9.882312e-05 
    ## ... Procrustes: rmse 0.0001946997  max resid 0.0002854714 
    ## ... Similar to previous best
    ## Run 428 stress 9.698757e-05 
    ## ... Procrustes: rmse 0.0001842832  max resid 0.0002811593 
    ## ... Similar to previous best
    ## Run 429 stress 9.789326e-05 
    ## ... Procrustes: rmse 0.000201085  max resid 0.0002986733 
    ## ... Similar to previous best
    ## Run 430 stress 9.647378e-05 
    ## ... Procrustes: rmse 0.0001658486  max resid 0.0002436461 
    ## ... Similar to previous best
    ## Run 431 stress 9.982075e-05 
    ## ... Procrustes: rmse 0.0002122582  max resid 0.0002991579 
    ## ... Similar to previous best
    ## Run 432 stress 9.613254e-05 
    ## ... Procrustes: rmse 0.0001245724  max resid 0.0002206808 
    ## ... Similar to previous best
    ## Run 433 stress 8.099868e-05 
    ## ... Procrustes: rmse 0.0001233667  max resid 0.000175073 
    ## ... Similar to previous best
    ## Run 434 stress 9.332343e-05 
    ## ... Procrustes: rmse 0.0002005294  max resid 0.0002872433 
    ## ... Similar to previous best
    ## Run 435 stress 9.250208e-05 
    ## ... Procrustes: rmse 0.0001196733  max resid 0.0001639649 
    ## ... Similar to previous best
    ## Run 436 stress 9.763504e-05 
    ## ... Procrustes: rmse 0.0001763648  max resid 0.0002188442 
    ## ... Similar to previous best
    ## Run 437 stress 9.890107e-05 
    ## ... Procrustes: rmse 0.0001982984  max resid 0.0002848916 
    ## ... Similar to previous best
    ## Run 438 stress 9.966858e-05 
    ## ... Procrustes: rmse 0.0002182758  max resid 0.0003064655 
    ## ... Similar to previous best
    ## Run 439 stress 9.201986e-05 
    ## ... Procrustes: rmse 0.0001998471  max resid 0.0002668097 
    ## ... Similar to previous best
    ## Run 440 stress 8.894955e-05 
    ## ... Procrustes: rmse 0.000175622  max resid 0.0002405939 
    ## ... Similar to previous best
    ## Run 441 stress 9.667623e-05 
    ## ... Procrustes: rmse 0.0002183643  max resid 0.0003044257 
    ## ... Similar to previous best
    ## Run 442 stress 8.945749e-05 
    ## ... Procrustes: rmse 0.0001820719  max resid 0.0002524061 
    ## ... Similar to previous best
    ## Run 443 stress 8.650407e-05 
    ## ... Procrustes: rmse 0.0001547003  max resid 0.0002601191 
    ## ... Similar to previous best
    ## Run 444 stress 9.949428e-05 
    ## ... Procrustes: rmse 0.0001809721  max resid 0.0002726276 
    ## ... Similar to previous best
    ## Run 445 stress 8.694075e-05 
    ## ... Procrustes: rmse 0.0001889579  max resid 0.0002501759 
    ## ... Similar to previous best
    ## Run 446 stress 8.550467e-05 
    ## ... Procrustes: rmse 0.0001474242  max resid 0.0001912465 
    ## ... Similar to previous best
    ## Run 447 stress 9.26902e-05 
    ## ... Procrustes: rmse 0.000178802  max resid 0.0002557592 
    ## ... Similar to previous best
    ## Run 448 stress 8.861193e-05 
    ## ... Procrustes: rmse 0.0001706151  max resid 0.0002374708 
    ## ... Similar to previous best
    ## Run 449 stress 8.479499e-05 
    ## ... Procrustes: rmse 0.0001551897  max resid 0.0002772997 
    ## ... Similar to previous best
    ## Run 450 stress 9.657097e-05 
    ## ... Procrustes: rmse 0.0001849491  max resid 0.0002480105 
    ## ... Similar to previous best
    ## Run 451 stress 9.176947e-05 
    ## ... Procrustes: rmse 0.0001821005  max resid 0.0002560385 
    ## ... Similar to previous best
    ## Run 452 stress 9.189141e-05 
    ## ... Procrustes: rmse 0.0001740347  max resid 0.0002472046 
    ## ... Similar to previous best
    ## Run 453 stress 7.158975e-05 
    ## ... Procrustes: rmse 0.0001253106  max resid 0.0001665083 
    ## ... Similar to previous best
    ## Run 454 stress 9.801326e-05 
    ## ... Procrustes: rmse 0.0001832298  max resid 0.0002692605 
    ## ... Similar to previous best
    ## Run 455 stress 9.075884e-05 
    ## ... Procrustes: rmse 0.0001640457  max resid 0.0002113855 
    ## ... Similar to previous best
    ## Run 456 stress 9.010297e-05 
    ## ... Procrustes: rmse 0.0001956492  max resid 0.0002867289 
    ## ... Similar to previous best
    ## Run 457 stress 9.490251e-05 
    ## ... Procrustes: rmse 0.0001712695  max resid 0.0002471926 
    ## ... Similar to previous best
    ## Run 458 stress 8.694659e-05 
    ## ... Procrustes: rmse 0.0001647657  max resid 0.0002245726 
    ## ... Similar to previous best
    ## Run 459 stress 9.843643e-05 
    ## ... Procrustes: rmse 0.000193504  max resid 0.0002563868 
    ## ... Similar to previous best
    ## Run 460 stress 9.783428e-05 
    ## ... Procrustes: rmse 0.0001895652  max resid 0.0002638052 
    ## ... Similar to previous best
    ## Run 461 stress 9.766207e-05 
    ## ... Procrustes: rmse 0.0002072741  max resid 0.0002879727 
    ## ... Similar to previous best
    ## Run 462 stress 9.497352e-05 
    ## ... Procrustes: rmse 0.0001825456  max resid 0.0002661979 
    ## ... Similar to previous best
    ## Run 463 stress 9.674709e-05 
    ## ... Procrustes: rmse 0.0001910103  max resid 0.0002665543 
    ## ... Similar to previous best
    ## Run 464 stress 9.595492e-05 
    ## ... Procrustes: rmse 0.0002062394  max resid 0.0002899984 
    ## ... Similar to previous best
    ## Run 465 stress 9.461507e-05 
    ## ... Procrustes: rmse 0.0001710273  max resid 0.0002195049 
    ## ... Similar to previous best
    ## Run 466 stress 9.75061e-05 
    ## ... Procrustes: rmse 0.0002086231  max resid 0.0002936962 
    ## ... Similar to previous best
    ## Run 467 stress 8.633633e-05 
    ## ... Procrustes: rmse 0.0001644892  max resid 0.0002277783 
    ## ... Similar to previous best
    ## Run 468 stress 9.01857e-05 
    ## ... Procrustes: rmse 0.00016552  max resid 0.000243693 
    ## ... Similar to previous best
    ## Run 469 stress 9.276585e-05 
    ## ... Procrustes: rmse 0.0001877771  max resid 0.0002557896 
    ## ... Similar to previous best
    ## Run 470 stress 9.823292e-05 
    ## ... Procrustes: rmse 0.0001748591  max resid 0.0002659065 
    ## ... Similar to previous best
    ## Run 471 stress 9.432517e-05 
    ## ... Procrustes: rmse 0.0002104558  max resid 0.0003039897 
    ## ... Similar to previous best
    ## Run 472 stress 9.93815e-05 
    ## ... Procrustes: rmse 0.0002144403  max resid 0.0003014987 
    ## ... Similar to previous best
    ## Run 473 stress 9.257889e-05 
    ## ... Procrustes: rmse 0.0001698117  max resid 0.000251542 
    ## ... Similar to previous best
    ## Run 474 stress 9.465733e-05 
    ## ... Procrustes: rmse 0.0001942847  max resid 0.0002846641 
    ## ... Similar to previous best
    ## Run 475 stress 9.78027e-05 
    ## ... Procrustes: rmse 0.0001334904  max resid 0.0002039734 
    ## ... Similar to previous best
    ## Run 476 stress 9.98633e-05 
    ## ... Procrustes: rmse 0.0001925398  max resid 0.000298844 
    ## ... Similar to previous best
    ## Run 477 stress 9.164999e-05 
    ## ... Procrustes: rmse 0.0001506803  max resid 0.0002015096 
    ## ... Similar to previous best
    ## Run 478 stress 9.54282e-05 
    ## ... Procrustes: rmse 0.0002060205  max resid 0.0002880103 
    ## ... Similar to previous best
    ## Run 479 stress 9.582581e-05 
    ## ... Procrustes: rmse 0.000186316  max resid 0.0002721201 
    ## ... Similar to previous best
    ## Run 480 stress 8.636738e-05 
    ## ... Procrustes: rmse 0.0001640461  max resid 0.0002630798 
    ## ... Similar to previous best
    ## Run 481 stress 9.737609e-05 
    ## ... Procrustes: rmse 0.0001964321  max resid 0.0002745636 
    ## ... Similar to previous best
    ## Run 482 stress 9.9757e-05 
    ## ... Procrustes: rmse 0.0001934946  max resid 0.000284919 
    ## ... Similar to previous best
    ## Run 483 stress 9.951113e-05 
    ## ... Procrustes: rmse 0.0001863472  max resid 0.0002862825 
    ## ... Similar to previous best
    ## Run 484 stress 9.862883e-05 
    ## ... Procrustes: rmse 0.0001897598  max resid 0.0002659639 
    ## ... Similar to previous best
    ## Run 485 stress 8.959475e-05 
    ## ... Procrustes: rmse 0.0001710828  max resid 0.000271106 
    ## ... Similar to previous best
    ## Run 486 stress 9.616518e-05 
    ## ... Procrustes: rmse 0.0001680848  max resid 0.0002629279 
    ## ... Similar to previous best
    ## Run 487 stress 9.659022e-05 
    ## ... Procrustes: rmse 0.0002004649  max resid 0.0002839858 
    ## ... Similar to previous best
    ## Run 488 stress 9.876017e-05 
    ## ... Procrustes: rmse 0.0001597908  max resid 0.0002348329 
    ## ... Similar to previous best
    ## Run 489 stress 9.200757e-05 
    ## ... Procrustes: rmse 0.0001491462  max resid 0.0002550818 
    ## ... Similar to previous best
    ## Run 490 stress 9.136952e-05 
    ## ... Procrustes: rmse 0.0001980139  max resid 0.0002908018 
    ## ... Similar to previous best
    ## Run 491 stress 9.282702e-05 
    ## ... Procrustes: rmse 0.0001853178  max resid 0.0002711061 
    ## ... Similar to previous best
    ## Run 492 stress 9.257951e-05 
    ## ... Procrustes: rmse 0.0001418262  max resid 0.000201741 
    ## ... Similar to previous best
    ## Run 493 stress 9.355237e-05 
    ## ... Procrustes: rmse 0.0001577479  max resid 0.0002370497 
    ## ... Similar to previous best
    ## Run 494 stress 9.678551e-05 
    ## ... Procrustes: rmse 0.00021231  max resid 0.0003066905 
    ## ... Similar to previous best
    ## Run 495 stress 9.259836e-05 
    ## ... Procrustes: rmse 0.000152019  max resid 0.0002182697 
    ## ... Similar to previous best
    ## Run 496 stress 9.676915e-05 
    ## ... Procrustes: rmse 0.0001864068  max resid 0.0002574405 
    ## ... Similar to previous best
    ## Run 497 stress 9.806084e-05 
    ## ... Procrustes: rmse 0.0001880511  max resid 0.0002822359 
    ## ... Similar to previous best
    ## Run 498 stress 9.927199e-05 
    ## ... Procrustes: rmse 0.0001998301  max resid 0.000268116 
    ## ... Similar to previous best
    ## Run 499 stress 9.212719e-05 
    ## ... Procrustes: rmse 0.0001891696  max resid 0.0002774385 
    ## ... Similar to previous best
    ## Run 500 stress 9.859685e-05 
    ## ... Procrustes: rmse 0.000213394  max resid 0.0002989235 
    ## ... Similar to previous best
    ## *** Best solution repeated 190 times

    ## Warning in metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Mixed lakes
SD_beta_env_M_NMDS <- metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.001279543 
    ## Run 1 stress 0.2509018 
    ## Run 2 stress 9.530432e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04314259  max resid 0.05937547 
    ## Run 3 stress 9.063325e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001988673  max resid 0.0004186471 
    ## ... Similar to previous best
    ## Run 4 stress 0.1491957 
    ## Run 5 stress 0.1776019 
    ## Run 6 stress 0.001200212 
    ## Run 7 stress 0.001068991 
    ## Run 8 stress 9.493862e-05 
    ## ... Procrustes: rmse 8.127995e-05  max resid 0.0001587538 
    ## ... Similar to previous best
    ## Run 9 stress 0.0001304869 
    ## ... Procrustes: rmse 0.008293833  max resid 0.01124753 
    ## Run 10 stress 9.407984e-05 
    ## ... Procrustes: rmse 0.000207473  max resid 0.0004362096 
    ## ... Similar to previous best
    ## Run 11 stress 9.761272e-05 
    ## ... Procrustes: rmse 0.0002173631  max resid 0.0004611981 
    ## ... Similar to previous best
    ## Run 12 stress 0.001180805 
    ## Run 13 stress 0.0004668667 
    ## ... Procrustes: rmse 0.01577982  max resid 0.02157834 
    ## Run 14 stress 8.83867e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.002410583  max resid 0.003309862 
    ## ... Similar to previous best
    ## Run 15 stress 0.0007729858 
    ## Run 16 stress 0.001299033 
    ## Run 17 stress 0.1491957 
    ## Run 18 stress 0.1491957 
    ## Run 19 stress 0.0004964594 
    ## ... Procrustes: rmse 0.01548599  max resid 0.02553646 
    ## Run 20 stress 8.464373e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.002483049  max resid 0.003496284 
    ## ... Similar to previous best
    ## Run 21 stress 0.001326822 
    ## Run 22 stress 0.001215203 
    ## Run 23 stress 9.644184e-05 
    ## ... Procrustes: rmse 0.000211095  max resid 0.0003097452 
    ## ... Similar to previous best
    ## Run 24 stress 0.3083098 
    ## Run 25 stress 0.0004872458 
    ## ... Procrustes: rmse 0.01615376  max resid 0.02229977 
    ## Run 26 stress 9.129198e-05 
    ## ... Procrustes: rmse 0.0001349572  max resid 0.0002222999 
    ## ... Similar to previous best
    ## Run 27 stress 0.001413505 
    ## Run 28 stress 9.892487e-05 
    ## ... Procrustes: rmse 0.0001713482  max resid 0.0003064807 
    ## ... Similar to previous best
    ## Run 29 stress 0.001381536 
    ## Run 30 stress 9.11671e-05 
    ## ... Procrustes: rmse 0.0001506694  max resid 0.0002821092 
    ## ... Similar to previous best
    ## Run 31 stress 9.253848e-05 
    ## ... Procrustes: rmse 0.0002045923  max resid 0.0002978738 
    ## ... Similar to previous best
    ## Run 32 stress 7.912299e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001204893  max resid 0.0001703096 
    ## ... Similar to previous best
    ## Run 33 stress 0.0004159855 
    ## ... Procrustes: rmse 0.01493609  max resid 0.02057098 
    ## Run 34 stress 0.001270473 
    ## Run 35 stress 0.0004869478 
    ## ... Procrustes: rmse 0.01616196  max resid 0.02226245 
    ## Run 36 stress 0.001548143 
    ## Run 37 stress 0.2354287 
    ## Run 38 stress 0.0005917903 
    ## Run 39 stress 0.1491957 
    ## Run 40 stress 9.758333e-05 
    ## ... Procrustes: rmse 0.0001575026  max resid 0.0002160905 
    ## ... Similar to previous best
    ## Run 41 stress 0.001218603 
    ## Run 42 stress 9.493099e-05 
    ## ... Procrustes: rmse 0.0001839273  max resid 0.0002989321 
    ## ... Similar to previous best
    ## Run 43 stress 0.001268309 
    ## Run 44 stress 0.2354286 
    ## Run 45 stress 8.300046e-05 
    ## ... Procrustes: rmse 0.001296883  max resid 0.001816209 
    ## ... Similar to previous best
    ## Run 46 stress 8.508663e-05 
    ## ... Procrustes: rmse 0.0004210696  max resid 0.0006141555 
    ## ... Similar to previous best
    ## Run 47 stress 9.789789e-05 
    ## ... Procrustes: rmse 0.0001660535  max resid 0.0002418037 
    ## ... Similar to previous best
    ## Run 48 stress 9.500784e-05 
    ## ... Procrustes: rmse 0.0001755901  max resid 0.0002742212 
    ## ... Similar to previous best
    ## Run 49 stress 0.001369619 
    ## Run 50 stress 0.1491957 
    ## Run 51 stress 0.0009448833 
    ## Run 52 stress 9.271142e-05 
    ## ... Procrustes: rmse 0.000136139  max resid 0.0002009126 
    ## ... Similar to previous best
    ## Run 53 stress 0.001255534 
    ## Run 54 stress 9.052101e-05 
    ## ... Procrustes: rmse 0.0001884312  max resid 0.0002326466 
    ## ... Similar to previous best
    ## Run 55 stress 0.001342064 
    ## Run 56 stress 8.539223e-05 
    ## ... Procrustes: rmse 0.0004356918  max resid 0.0006745438 
    ## ... Similar to previous best
    ## Run 57 stress 0.001296311 
    ## Run 58 stress 0.001311088 
    ## Run 59 stress 8.72552e-05 
    ## ... Procrustes: rmse 0.00662055  max resid 0.009097757 
    ## Run 60 stress 9.146264e-05 
    ## ... Procrustes: rmse 0.0001693365  max resid 0.0002883291 
    ## ... Similar to previous best
    ## Run 61 stress 0.3083098 
    ## Run 62 stress 0.001427286 
    ## Run 63 stress 6.977308e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.002241025  max resid 0.003103624 
    ## ... Similar to previous best
    ## Run 64 stress 0.001402296 
    ## Run 65 stress 0.1491957 
    ## Run 66 stress 0.001424546 
    ## Run 67 stress 0.0005168998 
    ## ... Procrustes: rmse 0.01441261  max resid 0.01983516 
    ## Run 68 stress 0.177601 
    ## Run 69 stress 0.001357951 
    ## Run 70 stress 0.1491957 
    ## Run 71 stress 9.654091e-05 
    ## ... Procrustes: rmse 0.00223105  max resid 0.002954945 
    ## ... Similar to previous best
    ## Run 72 stress 0.0006811409 
    ## Run 73 stress 0.001079414 
    ## Run 74 stress 9.918103e-05 
    ## ... Procrustes: rmse 0.002202602  max resid 0.002975566 
    ## ... Similar to previous best
    ## Run 75 stress 9.267136e-05 
    ## ... Procrustes: rmse 0.002266967  max resid 0.003227469 
    ## ... Similar to previous best
    ## Run 76 stress 9.126559e-05 
    ## ... Procrustes: rmse 0.002194052  max resid 0.002920148 
    ## ... Similar to previous best
    ## Run 77 stress 7.60579e-05 
    ## ... Procrustes: rmse 0.001931972  max resid 0.002679441 
    ## ... Similar to previous best
    ## Run 78 stress 9.845037e-05 
    ## ... Procrustes: rmse 0.002241298  max resid 0.003056843 
    ## ... Similar to previous best
    ## Run 79 stress 0.001163623 
    ## Run 80 stress 9.526548e-05 
    ## ... Procrustes: rmse 0.002202167  max resid 0.003014284 
    ## ... Similar to previous best
    ## Run 81 stress 9.415717e-05 
    ## ... Procrustes: rmse 0.002187265  max resid 0.002922506 
    ## ... Similar to previous best
    ## Run 82 stress 0.0004717504 
    ## ... Procrustes: rmse 0.01366553  max resid 0.01880416 
    ## Run 83 stress 9.325669e-05 
    ## ... Procrustes: rmse 0.002206417  max resid 0.003059405 
    ## ... Similar to previous best
    ## Run 84 stress 9.173202e-05 
    ## ... Procrustes: rmse 0.002193279  max resid 0.002961326 
    ## ... Similar to previous best
    ## Run 85 stress 9.646512e-05 
    ## ... Procrustes: rmse 0.009444016  max resid 0.01675277 
    ## Run 86 stress 9.948606e-05 
    ## ... Procrustes: rmse 0.002226708  max resid 0.002951121 
    ## ... Similar to previous best
    ## Run 87 stress 9.452727e-05 
    ## ... Procrustes: rmse 0.002225993  max resid 0.002930715 
    ## ... Similar to previous best
    ## Run 88 stress 9.932497e-05 
    ## ... Procrustes: rmse 0.002265521  max resid 0.003098095 
    ## ... Similar to previous best
    ## Run 89 stress 7.19206e-05 
    ## ... Procrustes: rmse 0.002423139  max resid 0.005448526 
    ## ... Similar to previous best
    ## Run 90 stress 0.001001375 
    ## Run 91 stress 0.000474705 
    ## ... Procrustes: rmse 0.01371946  max resid 0.0188781 
    ## Run 92 stress 9.741823e-05 
    ## ... Procrustes: rmse 0.002265766  max resid 0.003102562 
    ## ... Similar to previous best
    ## Run 93 stress 0.1491957 
    ## Run 94 stress 0.000472169 
    ## ... Procrustes: rmse 0.01367672  max resid 0.01881915 
    ## Run 95 stress 9.195069e-05 
    ## ... Procrustes: rmse 0.00225322  max resid 0.004423952 
    ## ... Similar to previous best
    ## Run 96 stress 0.177602 
    ## Run 97 stress 0.0002538465 
    ## ... Procrustes: rmse 0.01876966  max resid 0.02978772 
    ## Run 98 stress 0.0004844369 
    ## ... Procrustes: rmse 0.01388041  max resid 0.01910041 
    ## Run 99 stress 0.001485578 
    ## Run 100 stress 8.410374e-05 
    ## ... Procrustes: rmse 0.0021566  max resid 0.002936167 
    ## ... Similar to previous best
    ## Run 101 stress 9.864e-05 
    ## ... Procrustes: rmse 0.002264648  max resid 0.003097792 
    ## ... Similar to previous best
    ## Run 102 stress 0.2373687 
    ## Run 103 stress 0.2854408 
    ## Run 104 stress 0.00140012 
    ## Run 105 stress 0.0003819878 
    ## ... Procrustes: rmse 0.0120725  max resid 0.01660559 
    ## Run 106 stress 0.1990774 
    ## Run 107 stress 9.361818e-05 
    ## ... Procrustes: rmse 0.002163195  max resid 0.002899011 
    ## ... Similar to previous best
    ## Run 108 stress 0.1491957 
    ## Run 109 stress 0.001344925 
    ## Run 110 stress 0.001358304 
    ## Run 111 stress 0.001346405 
    ## Run 112 stress 0.001335148 
    ## Run 113 stress 9.539468e-05 
    ## ... Procrustes: rmse 0.002241648  max resid 0.003087599 
    ## ... Similar to previous best
    ## Run 114 stress 0.0004733569 
    ## ... Procrustes: rmse 0.01369572  max resid 0.01884374 
    ## Run 115 stress 9.261009e-05 
    ## ... Procrustes: rmse 0.00222247  max resid 0.002897336 
    ## ... Similar to previous best
    ## Run 116 stress 0.0007597521 
    ## Run 117 stress 9.598544e-05 
    ## ... Procrustes: rmse 0.002210591  max resid 0.002974558 
    ## ... Similar to previous best
    ## Run 118 stress 0.0004555946 
    ## ... Procrustes: rmse 0.01339445  max resid 0.01842944 
    ## Run 119 stress 0.0004744353 
    ## ... Procrustes: rmse 0.01371253  max resid 0.01886939 
    ## Run 120 stress 0.1491957 
    ## Run 121 stress 9.635551e-05 
    ## ... Procrustes: rmse 0.002203374  max resid 0.003069848 
    ## ... Similar to previous best
    ## Run 122 stress 8.236057e-05 
    ## ... Procrustes: rmse 0.002267633  max resid 0.003495425 
    ## ... Similar to previous best
    ## Run 123 stress 9.571675e-05 
    ## ... Procrustes: rmse 0.002192853  max resid 0.002975427 
    ## ... Similar to previous best
    ## Run 124 stress 0.0005124424 
    ## ... Procrustes: rmse 0.01396752  max resid 0.01922781 
    ## Run 125 stress 0.2570113 
    ## Run 126 stress 9.188418e-05 
    ## ... Procrustes: rmse 0.002265108  max resid 0.00307969 
    ## ... Similar to previous best
    ## Run 127 stress 9.226318e-05 
    ## ... Procrustes: rmse 0.002168905  max resid 0.002905512 
    ## ... Similar to previous best
    ## Run 128 stress 0.177601 
    ## Run 129 stress 9.775149e-05 
    ## ... Procrustes: rmse 0.002203617  max resid 0.003072695 
    ## ... Similar to previous best
    ## Run 130 stress 9.983833e-05 
    ## ... Procrustes: rmse 0.002183541  max resid 0.002963655 
    ## ... Similar to previous best
    ## Run 131 stress 0.001203444 
    ## Run 132 stress 9.141664e-05 
    ## ... Procrustes: rmse 0.002171011  max resid 0.002908277 
    ## ... Similar to previous best
    ## Run 133 stress 9.698115e-05 
    ## ... Procrustes: rmse 0.002331039  max resid 0.005148726 
    ## ... Similar to previous best
    ## Run 134 stress 9.69137e-05 
    ## ... Procrustes: rmse 0.006866473  max resid 0.01305559 
    ## Run 135 stress 0.0004990816 
    ## ... Procrustes: rmse 0.01412489  max resid 0.01943704 
    ## Run 136 stress 0.0003316397 
    ## ... Procrustes: rmse 0.01109449  max resid 0.0152561 
    ## Run 137 stress 0.1990774 
    ## Run 138 stress 0.0005084947 
    ## ... Procrustes: rmse 0.01427736  max resid 0.01965143 
    ## Run 139 stress 0.1491957 
    ## Run 140 stress 9.71043e-05 
    ## ... Procrustes: rmse 0.002265323  max resid 0.003098016 
    ## ... Similar to previous best
    ## Run 141 stress 0.0005017171 
    ## ... Procrustes: rmse 0.01416783  max resid 0.01949802 
    ## Run 142 stress 9.070161e-05 
    ## ... Procrustes: rmse 0.00226443  max resid 0.003076233 
    ## ... Similar to previous best
    ## Run 143 stress 0.1491957 
    ## Run 144 stress 0.001088156 
    ## Run 145 stress 9.429354e-05 
    ## ... Procrustes: rmse 0.002264863  max resid 0.003085037 
    ## ... Similar to previous best
    ## Run 146 stress 0.1491957 
    ## Run 147 stress 0.1491957 
    ## Run 148 stress 0.0004665641 
    ## ... Procrustes: rmse 0.01358172  max resid 0.01868812 
    ## Run 149 stress 9.716729e-05 
    ## ... Procrustes: rmse 0.002205809  max resid 0.003070458 
    ## ... Similar to previous best
    ## Run 150 stress 0.001334032 
    ## Run 151 stress 0.0013929 
    ## Run 152 stress 0.001142486 
    ## Run 153 stress 0.1776012 
    ## Run 154 stress 0.1990774 
    ## Run 155 stress 9.877699e-05 
    ## ... Procrustes: rmse 0.002222632  max resid 0.003051886 
    ## ... Similar to previous best
    ## Run 156 stress 0.0005949587 
    ## Run 157 stress 0.001277775 
    ## Run 158 stress 9.625139e-05 
    ## ... Procrustes: rmse 0.002199971  max resid 0.003062051 
    ## ... Similar to previous best
    ## Run 159 stress 0.257864 
    ## Run 160 stress 0.001225981 
    ## Run 161 stress 0.1990774 
    ## Run 162 stress 9.501795e-05 
    ## ... Procrustes: rmse 0.002195098  max resid 0.003029841 
    ## ... Similar to previous best
    ## Run 163 stress 0.000223153 
    ## ... Procrustes: rmse 0.0086928  max resid 0.01194186 
    ## Run 164 stress 0.0006026796 
    ## Run 165 stress 0.001254118 
    ## Run 166 stress 0.0004156248 
    ## ... Procrustes: rmse 0.01268625  max resid 0.01745249 
    ## Run 167 stress 0.0003439734 
    ## ... Procrustes: rmse 0.01131661  max resid 0.01556249 
    ## Run 168 stress 0.0004443407 
    ## ... Procrustes: rmse 0.01319493  max resid 0.01815366 
    ## Run 169 stress 0.1776019 
    ## Run 170 stress 0.0004822693 
    ## ... Procrustes: rmse 0.01384461  max resid 0.01905061 
    ## Run 171 stress 0.001370146 
    ## Run 172 stress 7.590612e-05 
    ## ... Procrustes: rmse 0.002970758  max resid 0.004020977 
    ## ... Similar to previous best
    ## Run 173 stress 9.901219e-05 
    ## ... Procrustes: rmse 0.002265096  max resid 0.003082623 
    ## ... Similar to previous best
    ## Run 174 stress 0.2509018 
    ## Run 175 stress 0.1491957 
    ## Run 176 stress 0.0001413361 
    ## ... Procrustes: rmse 0.006454326  max resid 0.008852771 
    ## Run 177 stress 8.927221e-05 
    ## ... Procrustes: rmse 0.0009533568  max resid 0.001290797 
    ## ... Similar to previous best
    ## Run 178 stress 0.1491957 
    ## Run 179 stress 8.714251e-05 
    ## ... Procrustes: rmse 0.002267291  max resid 0.002990594 
    ## ... Similar to previous best
    ## Run 180 stress 0.1990774 
    ## Run 181 stress 0.0001267234 
    ## ... Procrustes: rmse 0.005956991  max resid 0.008166357 
    ## Run 182 stress 0.0004583842 
    ## ... Procrustes: rmse 0.01344024  max resid 0.01849297 
    ## Run 183 stress 0.0004736721 
    ## ... Procrustes: rmse 0.0136396  max resid 0.01877096 
    ## Run 184 stress 0.001279288 
    ## Run 185 stress 9.558809e-05 
    ## ... Procrustes: rmse 0.002248543  max resid 0.003190543 
    ## ... Similar to previous best
    ## Run 186 stress 0.1776015 
    ## Run 187 stress 0.2528294 
    ## Run 188 stress 0.28294 
    ## Run 189 stress 8.008201e-05 
    ## ... Procrustes: rmse 0.002229931  max resid 0.003816347 
    ## ... Similar to previous best
    ## Run 190 stress 8.51633e-05 
    ## ... Procrustes: rmse 0.002273407  max resid 0.003184735 
    ## ... Similar to previous best
    ## Run 191 stress 9.266523e-05 
    ## ... Procrustes: rmse 0.00224296  max resid 0.003083991 
    ## ... Similar to previous best
    ## Run 192 stress 0.000428353 
    ## ... Procrustes: rmse 0.01290539  max resid 0.01775478 
    ## Run 193 stress 8.324597e-05 
    ## ... Procrustes: rmse 0.002209911  max resid 0.003034137 
    ## ... Similar to previous best
    ## Run 194 stress 0.001405565 
    ## Run 195 stress 0.001221917 
    ## Run 196 stress 0.2440272 
    ## Run 197 stress 0.1776012 
    ## Run 198 stress 0.2848525 
    ## Run 199 stress 9.601934e-05 
    ## ... Procrustes: rmse 0.002204087  max resid 0.00307012 
    ## ... Similar to previous best
    ## Run 200 stress 0.2568281 
    ## Run 201 stress 9.627768e-05 
    ## ... Procrustes: rmse 0.002265437  max resid 0.003079024 
    ## ... Similar to previous best
    ## Run 202 stress 0.001158328 
    ## Run 203 stress 9.1965e-05 
    ## ... Procrustes: rmse 0.002240516  max resid 0.003069743 
    ## ... Similar to previous best
    ## Run 204 stress 0.001093962 
    ## Run 205 stress 0.0009618888 
    ## Run 206 stress 9.102322e-05 
    ## ... Procrustes: rmse 0.002191544  max resid 0.002943518 
    ## ... Similar to previous best
    ## Run 207 stress 0.0007547942 
    ## Run 208 stress 0.1990774 
    ## Run 209 stress 9.791252e-05 
    ## ... Procrustes: rmse 0.002208699  max resid 0.003068291 
    ## ... Similar to previous best
    ## Run 210 stress 0.001492876 
    ## Run 211 stress 7.615435e-05 
    ## ... Procrustes: rmse 0.002185849  max resid 0.003721302 
    ## ... Similar to previous best
    ## Run 212 stress 9.243126e-05 
    ## ... Procrustes: rmse 0.002420826  max resid 0.005411609 
    ## ... Similar to previous best
    ## Run 213 stress 0.1491957 
    ## Run 214 stress 9.962876e-05 
    ## ... Procrustes: rmse 0.005037731  max resid 0.00689988 
    ## Run 215 stress 7.778873e-05 
    ## ... Procrustes: rmse 0.002203612  max resid 0.002931325 
    ## ... Similar to previous best
    ## Run 216 stress 0.001349774 
    ## Run 217 stress 0.2520602 
    ## Run 218 stress 0.001434129 
    ## Run 219 stress 9.797493e-05 
    ## ... Procrustes: rmse 0.002163944  max resid 0.002911786 
    ## ... Similar to previous best
    ## Run 220 stress 7.732758e-05 
    ## ... Procrustes: rmse 0.002215326  max resid 0.003019195 
    ## ... Similar to previous best
    ## Run 221 stress 9.018535e-05 
    ## ... Procrustes: rmse 0.002379968  max resid 0.003252538 
    ## ... Similar to previous best
    ## Run 222 stress 0.0004858726 
    ## ... Procrustes: rmse 0.01390617  max resid 0.01913558 
    ## Run 223 stress 0.0002820039 
    ## ... Procrustes: rmse 0.01005319  max resid 0.01381932 
    ## Run 224 stress 9.302189e-05 
    ## ... Procrustes: rmse 0.002223936  max resid 0.002897943 
    ## ... Similar to previous best
    ## Run 225 stress 0.0004619506 
    ## ... Procrustes: rmse 0.01327259  max resid 0.01826213 
    ## Run 226 stress 0.0012735 
    ## Run 227 stress 0.2520602 
    ## Run 228 stress 0.1990774 
    ## Run 229 stress 9.621578e-05 
    ## ... Procrustes: rmse 0.002241302  max resid 0.003091963 
    ## ... Similar to previous best
    ## Run 230 stress 5.717998e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.002267551  max resid 0.00302102 
    ## ... Similar to previous best
    ## Run 231 stress 0.1491957 
    ## Run 232 stress 0.2520602 
    ## Run 233 stress 8.377553e-05 
    ## ... Procrustes: rmse 9.115543e-05  max resid 0.0001872082 
    ## ... Similar to previous best
    ## Run 234 stress 0.001216004 
    ## Run 235 stress 9.246446e-05 
    ## ... Procrustes: rmse 9.546361e-05  max resid 0.0001406451 
    ## ... Similar to previous best
    ## Run 236 stress 0.0004661462 
    ## ... Procrustes: rmse 0.01583395  max resid 0.02169026 
    ## Run 237 stress 9.255893e-05 
    ## ... Procrustes: rmse 0.0001075027  max resid 0.0001709877 
    ## ... Similar to previous best
    ## Run 238 stress 0.0008054491 
    ## Run 239 stress 0.2440272 
    ## Run 240 stress 9.592323e-05 
    ## ... Procrustes: rmse 0.0001017567  max resid 0.0001527294 
    ## ... Similar to previous best
    ## Run 241 stress 8.65476e-05 
    ## ... Procrustes: rmse 0.002475241  max resid 0.0033167 
    ## ... Similar to previous best
    ## Run 242 stress 0.001321772 
    ## Run 243 stress 9.50884e-05 
    ## ... Procrustes: rmse 0.0001817971  max resid 0.0003409671 
    ## ... Similar to previous best
    ## Run 244 stress 8.758248e-05 
    ## ... Procrustes: rmse 0.0001587458  max resid 0.0003122386 
    ## ... Similar to previous best
    ## Run 245 stress 0.001400166 
    ## Run 246 stress 9.059874e-05 
    ## ... Procrustes: rmse 0.0001721947  max resid 0.000257548 
    ## ... Similar to previous best
    ## Run 247 stress 0.2842805 
    ## Run 248 stress 0.001334951 
    ## Run 249 stress 0.001354235 
    ## Run 250 stress 9.22993e-05 
    ## ... Procrustes: rmse 9.851993e-05  max resid 0.0001439995 
    ## ... Similar to previous best
    ## Run 251 stress 0.1990774 
    ## Run 252 stress 0.0004638497 
    ## ... Procrustes: rmse 0.0157815  max resid 0.02161764 
    ## Run 253 stress 9.442918e-05 
    ## ... Procrustes: rmse 0.0001053027  max resid 0.0001492773 
    ## ... Similar to previous best
    ## Run 254 stress 0.001342427 
    ## Run 255 stress 0.001226648 
    ## Run 256 stress 8.855607e-05 
    ## ... Procrustes: rmse 0.0002026198  max resid 0.0003496558 
    ## ... Similar to previous best
    ## Run 257 stress 8.486184e-05 
    ## ... Procrustes: rmse 0.008108418  max resid 0.01120922 
    ## Run 258 stress 9.288783e-05 
    ## ... Procrustes: rmse 0.0009652284  max resid 0.001365641 
    ## ... Similar to previous best
    ## Run 259 stress 0.1776015 
    ## Run 260 stress 0.0005238304 
    ## ... Procrustes: rmse 0.02759358  max resid 0.03807776 
    ## Run 261 stress 0.001369737 
    ## Run 262 stress 9.625835e-05 
    ## ... Procrustes: rmse 0.0001286475  max resid 0.0002246848 
    ## ... Similar to previous best
    ## Run 263 stress 9.119769e-05 
    ## ... Procrustes: rmse 0.003873107  max resid 0.005327581 
    ## ... Similar to previous best
    ## Run 264 stress 8.077834e-05 
    ## ... Procrustes: rmse 0.0004898816  max resid 0.0007008552 
    ## ... Similar to previous best
    ## Run 265 stress 0.001206811 
    ## Run 266 stress 0.001264456 
    ## Run 267 stress 0.0004859824 
    ## ... Procrustes: rmse 0.01615998  max resid 0.02214558 
    ## Run 268 stress 9.688109e-05 
    ## ... Procrustes: rmse 0.0007398641  max resid 0.0009811719 
    ## ... Similar to previous best
    ## Run 269 stress 9.267661e-05 
    ## ... Procrustes: rmse 9.71824e-05  max resid 0.0001428675 
    ## ... Similar to previous best
    ## Run 270 stress 0.0001840086 
    ## ... Procrustes: rmse 0.009899739  max resid 0.01350248 
    ## Run 271 stress 0.0004618125 
    ## ... Procrustes: rmse 0.01576455  max resid 0.02159645 
    ## Run 272 stress 0.1491957 
    ## Run 273 stress 0.1491957 
    ## Run 274 stress 0.264697 
    ## Run 275 stress 9.999801e-05 
    ## ... Procrustes: rmse 0.00733368  max resid 0.009960886 
    ## Run 276 stress 0.1776017 
    ## Run 277 stress 7.360489e-05 
    ## ... Procrustes: rmse 0.0008781781  max resid 0.001162605 
    ## ... Similar to previous best
    ## Run 278 stress 0.2520602 
    ## Run 279 stress 0.001058604 
    ## Run 280 stress 0.1491957 
    ## Run 281 stress 0.1990774 
    ## Run 282 stress 9.551039e-05 
    ## ... Procrustes: rmse 0.000159377  max resid 0.000311047 
    ## ... Similar to previous best
    ## Run 283 stress 7.840901e-05 
    ## ... Procrustes: rmse 0.0001568266  max resid 0.0002991742 
    ## ... Similar to previous best
    ## Run 284 stress 0.0004579526 
    ## ... Procrustes: rmse 0.01569872  max resid 0.02150541 
    ## Run 285 stress 0.0003024855 
    ## ... Procrustes: rmse 0.01275814  max resid 0.01744828 
    ## Run 286 stress 8.733291e-05 
    ## ... Procrustes: rmse 0.0001061454  max resid 0.0001811801 
    ## ... Similar to previous best
    ## Run 287 stress 9.638459e-05 
    ## ... Procrustes: rmse 0.0001064097  max resid 0.0001492005 
    ## ... Similar to previous best
    ## Run 288 stress 9.57424e-05 
    ## ... Procrustes: rmse 0.000168569  max resid 0.0002595014 
    ## ... Similar to previous best
    ## Run 289 stress 0.0009504169 
    ## Run 290 stress 0.0004657116 
    ## ... Procrustes: rmse 0.01583079  max resid 0.02168719 
    ## Run 291 stress 9.836511e-05 
    ## ... Procrustes: rmse 0.004814149  max resid 0.006538825 
    ## ... Similar to previous best
    ## Run 292 stress 9.448037e-05 
    ## ... Procrustes: rmse 0.000107366  max resid 0.0001556866 
    ## ... Similar to previous best
    ## Run 293 stress 0.001246695 
    ## Run 294 stress 0.1990774 
    ## Run 295 stress 0.1491957 
    ## Run 296 stress 0.0004438368 
    ## ... Procrustes: rmse 0.02538215  max resid 0.03501327 
    ## Run 297 stress 0.1491957 
    ## Run 298 stress 0.0009507348 
    ## Run 299 stress 9.376714e-05 
    ## ... Procrustes: rmse 0.0002141264  max resid 0.0003771302 
    ## ... Similar to previous best
    ## Run 300 stress 9.30481e-05 
    ## ... Procrustes: rmse 0.0001064212  max resid 0.0001334248 
    ## ... Similar to previous best
    ## Run 301 stress 9.493976e-05 
    ## ... Procrustes: rmse 0.0001312488  max resid 0.0002429125 
    ## ... Similar to previous best
    ## Run 302 stress 0.001325024 
    ## Run 303 stress 0.1491957 
    ## Run 304 stress 0.1491957 
    ## Run 305 stress 0.0004624871 
    ## ... Procrustes: rmse 0.01577572  max resid 0.02161134 
    ## Run 306 stress 0.000440051 
    ## ... Procrustes: rmse 0.0153811  max resid 0.02106958 
    ## Run 307 stress 0.2842805 
    ## Run 308 stress 0.2852122 
    ## Run 309 stress 0.001167116 
    ## Run 310 stress 0.001373775 
    ## Run 311 stress 9.997717e-05 
    ## ... Procrustes: rmse 0.0001901507  max resid 0.0002753302 
    ## ... Similar to previous best
    ## Run 312 stress 0.0001319958 
    ## ... Procrustes: rmse 0.008424851  max resid 0.01146809 
    ## Run 313 stress 0.0001452838 
    ## ... Procrustes: rmse 0.008803177  max resid 0.01198919 
    ## Run 314 stress 0.0004833454 
    ## ... Procrustes: rmse 0.0158089  max resid 0.02165774 
    ## Run 315 stress 0.00136773 
    ## Run 316 stress 9.133053e-05 
    ## ... Procrustes: rmse 0.003129211  max resid 0.004178925 
    ## ... Similar to previous best
    ## Run 317 stress 0.00110609 
    ## Run 318 stress 0.0009196712 
    ## Run 319 stress 0.0004448481 
    ## ... Procrustes: rmse 0.01547252  max resid 0.02119317 
    ## Run 320 stress 0.0004502054 
    ## ... Procrustes: rmse 0.0155647  max resid 0.02132008 
    ## Run 321 stress 9.010939e-05 
    ## ... Procrustes: rmse 0.0001729605  max resid 0.000251102 
    ## ... Similar to previous best
    ## Run 322 stress 0.001323205 
    ## Run 323 stress 9.383348e-05 
    ## ... Procrustes: rmse 0.0001277409  max resid 0.0002387249 
    ## ... Similar to previous best
    ## Run 324 stress 0.001371995 
    ## Run 325 stress 0.001334303 
    ## Run 326 stress 8.917901e-05 
    ## ... Procrustes: rmse 9.416738e-05  max resid 0.0001361946 
    ## ... Similar to previous best
    ## Run 327 stress 0.177601 
    ## Run 328 stress 0.2842805 
    ## Run 329 stress 0.1491957 
    ## Run 330 stress 9.937747e-05 
    ## ... Procrustes: rmse 0.0001187421  max resid 0.00015508 
    ## ... Similar to previous best
    ## Run 331 stress 0.3083098 
    ## Run 332 stress 0.001180309 
    ## Run 333 stress 0.3083098 
    ## Run 334 stress 0.0004736915 
    ## ... Procrustes: rmse 0.01596477  max resid 0.02187212 
    ## Run 335 stress 8.534148e-05 
    ## ... Procrustes: rmse 0.006999345  max resid 0.009631641 
    ## Run 336 stress 0.3083098 
    ## Run 337 stress 9.631185e-05 
    ## ... Procrustes: rmse 0.0002187378  max resid 0.0003856627 
    ## ... Similar to previous best
    ## Run 338 stress 0.001293189 
    ## Run 339 stress 0.001280546 
    ## Run 340 stress 3.65717e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0009484543  max resid 0.001281629 
    ## ... Similar to previous best
    ## Run 341 stress 9.990154e-05 
    ## ... Procrustes: rmse 0.0009662824  max resid 0.001516666 
    ## ... Similar to previous best
    ## Run 342 stress 9.124449e-05 
    ## ... Procrustes: rmse 0.0009603983  max resid 0.001489385 
    ## ... Similar to previous best
    ## Run 343 stress 0.1491957 
    ## Run 344 stress 0.1491957 
    ## Run 345 stress 0.0002211654 
    ## ... Procrustes: rmse 0.016959  max resid 0.02338545 
    ## Run 346 stress 0.0006788763 
    ## Run 347 stress 9.803191e-05 
    ## ... Procrustes: rmse 0.000913578  max resid 0.00133475 
    ## ... Similar to previous best
    ## Run 348 stress 9.235193e-05 
    ## ... Procrustes: rmse 0.0009497087  max resid 0.001435518 
    ## ... Similar to previous best
    ## Run 349 stress 0.2373687 
    ## Run 350 stress 0.2797323 
    ## Run 351 stress 0.001189954 
    ## Run 352 stress 0.001230063 
    ## Run 353 stress 0.1491957 
    ## Run 354 stress 0.2373687 
    ## Run 355 stress 0.000468573 
    ## ... Procrustes: rmse 0.01564314  max resid 0.02337279 
    ## Run 356 stress 9.782775e-05 
    ## ... Procrustes: rmse 0.0009568267  max resid 0.001489664 
    ## ... Similar to previous best
    ## Run 357 stress 0.001289753 
    ## Run 358 stress 0.0009767609 
    ## Run 359 stress 7.851782e-05 
    ## ... Procrustes: rmse 0.001018096  max resid 0.001973162 
    ## ... Similar to previous best
    ## Run 360 stress 0.001407448 
    ## Run 361 stress 0.001306767 
    ## Run 362 stress 0.001181291 
    ## Run 363 stress 0.001256097 
    ## Run 364 stress 9.515438e-05 
    ## ... Procrustes: rmse 0.003020943  max resid 0.005786146 
    ## ... Similar to previous best
    ## Run 365 stress 0.3083098 
    ## Run 366 stress 0.001304436 
    ## Run 367 stress 0.001282933 
    ## Run 368 stress 0.2528294 
    ## Run 369 stress 0.1990774 
    ## Run 370 stress 9.624409e-05 
    ## ... Procrustes: rmse 0.0009903066  max resid 0.001433903 
    ## ... Similar to previous best
    ## Run 371 stress 0.1491957 
    ## Run 372 stress 9.841377e-05 
    ## ... Procrustes: rmse 0.0009979584  max resid 0.001320414 
    ## ... Similar to previous best
    ## Run 373 stress 0.0004429638 
    ## ... Procrustes: rmse 0.01520287  max resid 0.02276381 
    ## Run 374 stress 0.0004799965 
    ## ... Procrustes: rmse 0.01583281  max resid 0.02363535 
    ## Run 375 stress 8.732308e-05 
    ## ... Procrustes: rmse 0.00101588  max resid 0.001817718 
    ## ... Similar to previous best
    ## Run 376 stress 0.001188016 
    ## Run 377 stress 0.0001134847 
    ## ... Procrustes: rmse 0.007609185  max resid 0.0122422 
    ## Run 378 stress 8.640207e-05 
    ## ... Procrustes: rmse 0.0009955837  max resid 0.00145783 
    ## ... Similar to previous best
    ## Run 379 stress 9.1496e-05 
    ## ... Procrustes: rmse 0.007059937  max resid 0.009739961 
    ## Run 380 stress 9.725057e-05 
    ## ... Procrustes: rmse 0.0009876326  max resid 0.001427528 
    ## ... Similar to previous best
    ## Run 381 stress 9.175468e-05 
    ## ... Procrustes: rmse 0.0009544899  max resid 0.001460677 
    ## ... Similar to previous best
    ## Run 382 stress 4.967095e-05 
    ## ... Procrustes: rmse 0.001021033  max resid 0.00134186 
    ## ... Similar to previous best
    ## Run 383 stress 9.310908e-05 
    ## ... Procrustes: rmse 0.0008890002  max resid 0.001219428 
    ## ... Similar to previous best
    ## Run 384 stress 0.0004188561 
    ## ... Procrustes: rmse 0.01477751  max resid 0.02217329 
    ## Run 385 stress 9.354656e-05 
    ## ... Procrustes: rmse 0.0009677818  max resid 0.001469162 
    ## ... Similar to previous best
    ## Run 386 stress 8.825787e-05 
    ## ... Procrustes: rmse 0.0009610662  max resid 0.001458854 
    ## ... Similar to previous best
    ## Run 387 stress 0.001330753 
    ## Run 388 stress 0.0001047481 
    ## ... Procrustes: rmse 0.007252137  max resid 0.01174524 
    ## Run 389 stress 0.1776018 
    ## Run 390 stress 0.001171173 
    ## Run 391 stress 0.001342903 
    ## Run 392 stress 9.561116e-05 
    ## ... Procrustes: rmse 0.0009543892  max resid 0.001479319 
    ## ... Similar to previous best
    ## Run 393 stress 9.818636e-05 
    ## ... Procrustes: rmse 0.0009425255  max resid 0.001406259 
    ## ... Similar to previous best
    ## Run 394 stress 8.906811e-05 
    ## ... Procrustes: rmse 0.0009315237  max resid 0.001349143 
    ## ... Similar to previous best
    ## Run 395 stress 9.34324e-05 
    ## ... Procrustes: rmse 0.001066724  max resid 0.002314447 
    ## ... Similar to previous best
    ## Run 396 stress 0.0012895 
    ## Run 397 stress 0.001394707 
    ## Run 398 stress 8.090693e-05 
    ## ... Procrustes: rmse 0.001025816  max resid 0.001420564 
    ## ... Similar to previous best
    ## Run 399 stress 0.2842805 
    ## Run 400 stress 0.0013912 
    ## Run 401 stress 0.0006459562 
    ## Run 402 stress 0.1776011 
    ## Run 403 stress 8.992203e-05 
    ## ... Procrustes: rmse 0.0009398537  max resid 0.001245859 
    ## ... Similar to previous best
    ## Run 404 stress 0.001280193 
    ## Run 405 stress 0.001245747 
    ## Run 406 stress 0.2373687 
    ## Run 407 stress 9.72116e-05 
    ## ... Procrustes: rmse 0.0008808525  max resid 0.001230039 
    ## ... Similar to previous best
    ## Run 408 stress 3.662008e-05 
    ## ... Procrustes: rmse 0.0009938826  max resid 0.001384777 
    ## ... Similar to previous best
    ## Run 409 stress 0.001264618 
    ## Run 410 stress 0.001367782 
    ## Run 411 stress 0.2528294 
    ## Run 412 stress 0.0004608084 
    ## ... Procrustes: rmse 0.01547895  max resid 0.02314668 
    ## Run 413 stress 8.372799e-05 
    ## ... Procrustes: rmse 0.001027331  max resid 0.001423708 
    ## ... Similar to previous best
    ## Run 414 stress 0.1990774 
    ## Run 415 stress 0.001348643 
    ## Run 416 stress 8.360943e-05 
    ## ... Procrustes: rmse 0.0009037253  max resid 0.001234489 
    ## ... Similar to previous best
    ## Run 417 stress 0.001210713 
    ## Run 418 stress 9.25879e-05 
    ## ... Procrustes: rmse 0.000894447  max resid 0.001235958 
    ## ... Similar to previous best
    ## Run 419 stress 7.494501e-05 
    ## ... Procrustes: rmse 0.0009299159  max resid 0.001265119 
    ## ... Similar to previous best
    ## Run 420 stress 0.001066733 
    ## Run 421 stress 0.00041954 
    ## ... Procrustes: rmse 0.0147907  max resid 0.02219508 
    ## Run 422 stress 0.1776011 
    ## Run 423 stress 0.2795549 
    ## Run 424 stress 0.001192435 
    ## Run 425 stress 9.663644e-05 
    ## ... Procrustes: rmse 0.0008811375  max resid 0.001228527 
    ## ... Similar to previous best
    ## Run 426 stress 0.001337806 
    ## Run 427 stress 0.1776017 
    ## Run 428 stress 8.488986e-05 
    ## ... Procrustes: rmse 0.001013944  max resid 0.001437646 
    ## ... Similar to previous best
    ## Run 429 stress 0.001290384 
    ## Run 430 stress 9.310159e-05 
    ## ... Procrustes: rmse 0.0009570865  max resid 0.001471984 
    ## ... Similar to previous best
    ## Run 431 stress 0.001122174 
    ## Run 432 stress 8.004252e-05 
    ## ... Procrustes: rmse 0.000905144  max resid 0.001224294 
    ## ... Similar to previous best
    ## Run 433 stress 9.114959e-05 
    ## ... Procrustes: rmse 0.0009259515  max resid 0.001314359 
    ## ... Similar to previous best
    ## Run 434 stress 0.177602 
    ## Run 435 stress 8.675821e-05 
    ## ... Procrustes: rmse 0.001013249  max resid 0.001436828 
    ## ... Similar to previous best
    ## Run 436 stress 0.1990774 
    ## Run 437 stress 0.00121198 
    ## Run 438 stress 0.3082974 
    ## Run 439 stress 0.0002356606 
    ## ... Procrustes: rmse 0.0110371  max resid 0.01699993 
    ## Run 440 stress 0.1990774 
    ## Run 441 stress 9.10924e-05 
    ## ... Procrustes: rmse 0.0009576007  max resid 0.001446231 
    ## ... Similar to previous best
    ## Run 442 stress 0.1990774 
    ## Run 443 stress 0.001246123 
    ## Run 444 stress 9.784201e-05 
    ## ... Procrustes: rmse 0.0009148536  max resid 0.001315957 
    ## ... Similar to previous best
    ## Run 445 stress 0.1776015 
    ## Run 446 stress 0.1776011 
    ## Run 447 stress 0.001373285 
    ## Run 448 stress 0.1990774 
    ## Run 449 stress 0.001301336 
    ## Run 450 stress 9.543928e-05 
    ## ... Procrustes: rmse 0.0004639475  max resid 0.0006928428 
    ## ... Similar to previous best
    ## Run 451 stress 0.1491957 
    ## Run 452 stress 0.001272573 
    ## Run 453 stress 0.0002462836 
    ## ... Procrustes: rmse 0.01128723  max resid 0.01734661 
    ## Run 454 stress 0.1491957 
    ## Run 455 stress 9.839414e-05 
    ## ... Procrustes: rmse 0.001007706  max resid 0.00153431 
    ## ... Similar to previous best
    ## Run 456 stress 0.0004522296 
    ## ... Procrustes: rmse 0.01536223  max resid 0.02298437 
    ## Run 457 stress 0.001173416 
    ## Run 458 stress 9.467428e-05 
    ## ... Procrustes: rmse 0.0009482216  max resid 0.001440455 
    ## ... Similar to previous best
    ## Run 459 stress 9.927208e-05 
    ## ... Procrustes: rmse 0.0009475238  max resid 0.001461277 
    ## ... Similar to previous best
    ## Run 460 stress 8.529012e-05 
    ## ... Procrustes: rmse 0.0009848191  max resid 0.001356076 
    ## ... Similar to previous best
    ## Run 461 stress 0.001359554 
    ## Run 462 stress 8.531617e-05 
    ## ... Procrustes: rmse 0.001024179  max resid 0.001472517 
    ## ... Similar to previous best
    ## Run 463 stress 9.203266e-05 
    ## ... Procrustes: rmse 0.0009276643  max resid 0.001336198 
    ## ... Similar to previous best
    ## Run 464 stress 0.001163842 
    ## Run 465 stress 0.001418008 
    ## Run 466 stress 0.2520602 
    ## Run 467 stress 0.1491957 
    ## Run 468 stress 9.460975e-05 
    ## ... Procrustes: rmse 0.0009515834  max resid 0.00144861 
    ## ... Similar to previous best
    ## Run 469 stress 0.0005057782 
    ## ... Procrustes: rmse 0.01624896  max resid 0.02420965 
    ## Run 470 stress 0.1776015 
    ## Run 471 stress 0.0004856345 
    ## ... Procrustes: rmse 0.01592796  max resid 0.02376637 
    ## Run 472 stress 0.1990774 
    ## Run 473 stress 0.0008290189 
    ## Run 474 stress 0.2373687 
    ## Run 475 stress 0.0004820412 
    ## ... Procrustes: rmse 0.01576805  max resid 0.02354526 
    ## Run 476 stress 9.072836e-05 
    ## ... Procrustes: rmse 0.0009208021  max resid 0.001312791 
    ## ... Similar to previous best
    ## Run 477 stress 0.177602 
    ## Run 478 stress 0.3083099 
    ## Run 479 stress 9.7295e-05 
    ## ... Procrustes: rmse 0.0009698077  max resid 0.001469606 
    ## ... Similar to previous best
    ## Run 480 stress 0.001280373 
    ## Run 481 stress 0.1990774 
    ## Run 482 stress 0.2570108 
    ## Run 483 stress 0.2440272 
    ## Run 484 stress 8.031833e-05 
    ## ... Procrustes: rmse 0.0009621337  max resid 0.001252844 
    ## ... Similar to previous best
    ## Run 485 stress 0.0003934756 
    ## ... Procrustes: rmse 0.02295533  max resid 0.03168173 
    ## Run 486 stress 0.2354286 
    ## Run 487 stress 0.00133147 
    ## Run 488 stress 9.941261e-05 
    ## ... Procrustes: rmse 0.0009584555  max resid 0.001468881 
    ## ... Similar to previous best
    ## Run 489 stress 0.001277124 
    ## Run 490 stress 0.1491957 
    ## Run 491 stress 0.0004676987 
    ## ... Procrustes: rmse 0.01561949  max resid 0.02333961 
    ## Run 492 stress 8.075345e-05 
    ## ... Procrustes: rmse 0.001205083  max resid 0.002761448 
    ## ... Similar to previous best
    ## Run 493 stress 9.564289e-05 
    ## ... Procrustes: rmse 0.0009385003  max resid 0.00125081 
    ## ... Similar to previous best
    ## Run 494 stress 8.335101e-05 
    ## ... Procrustes: rmse 0.0009996918  max resid 0.001475235 
    ## ... Similar to previous best
    ## Run 495 stress 0.1990774 
    ## Run 496 stress 9.115754e-05 
    ## ... Procrustes: rmse 0.0009694692  max resid 0.001496653 
    ## ... Similar to previous best
    ## Run 497 stress 0.001403976 
    ## Run 498 stress 9.00475e-05 
    ## ... Procrustes: rmse 0.0009682823  max resid 0.001448264 
    ## ... Similar to previous best
    ## Run 499 stress 0.000482297 
    ## ... Procrustes: rmse 0.0158669  max resid 0.02368201 
    ## Run 500 stress 0.00147094 
    ## *** Best solution repeated 55 times

    ## Warning in metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
SD_beta_geo_NMDS <- metaMDS(SD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09021173 
    ## Run 1 stress 0.09503459 
    ## Run 2 stress 0.09503418 
    ## Run 3 stress 0.1074433 
    ## Run 4 stress 0.1062596 
    ## Run 5 stress 0.08946652 
    ## ... New best solution
    ## ... Procrustes: rmse 0.02559492  max resid 0.07781279 
    ## Run 6 stress 0.08926098 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03327807  max resid 0.1157532 
    ## Run 7 stress 0.1071038 
    ## Run 8 stress 0.08938539 
    ## ... Procrustes: rmse 0.01181979  max resid 0.03982073 
    ## Run 9 stress 0.08926071 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000391196  max resid 0.0009876583 
    ## ... Similar to previous best
    ## Run 10 stress 0.1111383 
    ## Run 11 stress 0.109188 
    ## Run 12 stress 0.1074432 
    ## Run 13 stress 0.1061302 
    ## Run 14 stress 0.08951728 
    ## ... Procrustes: rmse 0.03498179  max resid 0.1155898 
    ## Run 15 stress 0.08926109 
    ## ... Procrustes: rmse 0.0004464214  max resid 0.001156221 
    ## ... Similar to previous best
    ## Run 16 stress 0.08926076 
    ## ... Procrustes: rmse 0.0002026478  max resid 0.0004685252 
    ## ... Similar to previous best
    ## Run 17 stress 0.09018709 
    ## Run 18 stress 0.1075421 
    ## Run 19 stress 0.1108229 
    ## Run 20 stress 0.1067595 
    ## Run 21 stress 0.09503449 
    ## Run 22 stress 0.1067397 
    ## Run 23 stress 0.1075421 
    ## Run 24 stress 0.1064189 
    ## Run 25 stress 0.2559422 
    ## Run 26 stress 0.1056897 
    ## Run 27 stress 0.1101456 
    ## Run 28 stress 0.1076304 
    ## Run 29 stress 0.09099532 
    ## Run 30 stress 0.09503446 
    ## Run 31 stress 0.09109109 
    ## Run 32 stress 0.0904461 
    ## Run 33 stress 0.08946328 
    ## ... Procrustes: rmse 0.03747952  max resid 0.1177654 
    ## Run 34 stress 0.1087584 
    ## Run 35 stress 0.0894634 
    ## ... Procrustes: rmse 0.03743955  max resid 0.1177012 
    ## Run 36 stress 0.106043 
    ## Run 37 stress 0.08938963 
    ## ... Procrustes: rmse 0.03594633  max resid 0.1183095 
    ## Run 38 stress 0.09044609 
    ## Run 39 stress 0.08926097 
    ## ... Procrustes: rmse 0.0003595372  max resid 0.0009649763 
    ## ... Similar to previous best
    ## Run 40 stress 0.1066193 
    ## Run 41 stress 0.09109109 
    ## Run 42 stress 0.08938539 
    ## ... Procrustes: rmse 0.0117527  max resid 0.0394942 
    ## Run 43 stress 0.08938972 
    ## ... Procrustes: rmse 0.0358939  max resid 0.1182231 
    ## Run 44 stress 0.09178285 
    ## Run 45 stress 0.1076304 
    ## Run 46 stress 0.1060434 
    ## Run 47 stress 0.08938539 
    ## ... Procrustes: rmse 0.01175268  max resid 0.039492 
    ## Run 48 stress 0.1074433 
    ## Run 49 stress 0.08938968 
    ## ... Procrustes: rmse 0.03590525  max resid 0.1182381 
    ## Run 50 stress 0.09503416 
    ## Run 51 stress 0.1088391 
    ## Run 52 stress 0.1092209 
    ## Run 53 stress 0.09503463 
    ## Run 54 stress 0.1052647 
    ## Run 55 stress 0.08938963 
    ## ... Procrustes: rmse 0.03594523  max resid 0.1183085 
    ## Run 56 stress 0.08938963 
    ## ... Procrustes: rmse 0.03594429  max resid 0.1183032 
    ## Run 57 stress 0.08938544 
    ## ... Procrustes: rmse 0.01168233  max resid 0.03943337 
    ## Run 58 stress 0.09021165 
    ## Run 59 stress 0.08938964 
    ## ... Procrustes: rmse 0.03591439  max resid 0.1182402 
    ## Run 60 stress 0.1063193 
    ## Run 61 stress 0.1105335 
    ## Run 62 stress 0.1075421 
    ## Run 63 stress 0.0893898 
    ## ... Procrustes: rmse 0.03588421  max resid 0.1182081 
    ## Run 64 stress 0.08926073 
    ## ... Procrustes: rmse 7.85346e-05  max resid 0.000240812 
    ## ... Similar to previous best
    ## Run 65 stress 0.08951726 
    ## ... Procrustes: rmse 0.03499251  max resid 0.1156282 
    ## Run 66 stress 0.105265 
    ## Run 67 stress 0.08946332 
    ## ... Procrustes: rmse 0.03745465  max resid 0.1177261 
    ## Run 68 stress 0.08926067 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001388569  max resid 0.0003392165 
    ## ... Similar to previous best
    ## Run 69 stress 0.08946655 
    ## ... Procrustes: rmse 0.03323839  max resid 0.1164061 
    ## Run 70 stress 0.09262351 
    ## Run 71 stress 0.08938963 
    ## ... Procrustes: rmse 0.0359377  max resid 0.1182703 
    ## Run 72 stress 0.08938967 
    ## ... Procrustes: rmse 0.03596985  max resid 0.1183228 
    ## Run 73 stress 0.09018707 
    ## Run 74 stress 0.09180881 
    ## Run 75 stress 0.109213 
    ## Run 76 stress 0.09039142 
    ## Run 77 stress 0.09503418 
    ## Run 78 stress 0.1088372 
    ## Run 79 stress 0.09612885 
    ## Run 80 stress 0.09180892 
    ## Run 81 stress 0.08938549 
    ## ... Procrustes: rmse 0.01188978  max resid 0.03966387 
    ## Run 82 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595506  max resid 0.1182988 
    ## Run 83 stress 0.1067894 
    ## Run 84 stress 0.1052647 
    ## Run 85 stress 0.1056897 
    ## Run 86 stress 0.1060432 
    ## Run 87 stress 0.1065358 
    ## Run 88 stress 0.09109111 
    ## Run 89 stress 0.09130094 
    ## Run 90 stress 0.1056904 
    ## Run 91 stress 0.09044607 
    ## Run 92 stress 0.1056901 
    ## Run 93 stress 0.1092129 
    ## Run 94 stress 0.09130089 
    ## Run 95 stress 0.09021171 
    ## Run 96 stress 0.1056901 
    ## Run 97 stress 0.08938546 
    ## ... Procrustes: rmse 0.01180211  max resid 0.03963573 
    ## Run 98 stress 0.08946653 
    ## ... Procrustes: rmse 0.03320687  max resid 0.1163273 
    ## Run 99 stress 0.08926077 
    ## ... Procrustes: rmse 0.0002177642  max resid 0.0006831119 
    ## ... Similar to previous best
    ## Run 100 stress 0.08938966 
    ## ... Procrustes: rmse 0.03592807  max resid 0.1182542 
    ## Run 101 stress 0.09021165 
    ## Run 102 stress 0.1091825 
    ## Run 103 stress 0.1105336 
    ## Run 104 stress 0.08938547 
    ## ... Procrustes: rmse 0.01189946  max resid 0.03972683 
    ## Run 105 stress 0.1074432 
    ## Run 106 stress 0.1074433 
    ## Run 107 stress 0.1075421 
    ## Run 108 stress 0.08946662 
    ## ... Procrustes: rmse 0.03318794  max resid 0.1162977 
    ## Run 109 stress 0.09228484 
    ## Run 110 stress 0.09109109 
    ## Run 111 stress 0.0893856 
    ## ... Procrustes: rmse 0.01188226  max resid 0.03963028 
    ## Run 112 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183251  max resid 0.03971445 
    ## Run 113 stress 0.09503415 
    ## Run 114 stress 0.108839 
    ## Run 115 stress 0.0892607 
    ## ... Procrustes: rmse 0.0001498523  max resid 0.0003668195 
    ## ... Similar to previous best
    ## Run 116 stress 0.09775551 
    ## Run 117 stress 0.08946653 
    ## ... Procrustes: rmse 0.03320794  max resid 0.1163337 
    ## Run 118 stress 0.1101443 
    ## Run 119 stress 0.1066196 
    ## Run 120 stress 0.1065351 
    ## Run 121 stress 0.09018707 
    ## Run 122 stress 0.1095631 
    ## Run 123 stress 0.09021167 
    ## Run 124 stress 0.1087575 
    ## Run 125 stress 0.08938542 
    ## ... Procrustes: rmse 0.01188614  max resid 0.03971318 
    ## Run 126 stress 0.0893854 
    ## ... Procrustes: rmse 0.01184278  max resid 0.0396649 
    ## Run 127 stress 0.1075421 
    ## Run 128 stress 0.1087863 
    ## Run 129 stress 0.09021169 
    ## Run 130 stress 0.08926081 
    ## ... Procrustes: rmse 0.0002307107  max resid 0.000695656 
    ## ... Similar to previous best
    ## Run 131 stress 0.1076302 
    ## Run 132 stress 0.08938976 
    ## ... Procrustes: rmse 0.03591511  max resid 0.1182359 
    ## Run 133 stress 0.09018708 
    ## Run 134 stress 0.09018707 
    ## Run 135 stress 0.09039129 
    ## Run 136 stress 0.1056896 
    ## Run 137 stress 0.1061303 
    ## Run 138 stress 0.09503433 
    ## Run 139 stress 0.1092095 
    ## Run 140 stress 0.1092946 
    ## Run 141 stress 0.08938961 
    ## ... Procrustes: rmse 0.03594588  max resid 0.1182857 
    ## Run 142 stress 0.08938539 
    ## ... Procrustes: rmse 0.01185213  max resid 0.03973374 
    ## Run 143 stress 0.0895173 
    ## ... Procrustes: rmse 0.03510585  max resid 0.1157018 
    ## Run 144 stress 0.09612803 
    ## Run 145 stress 0.08938968 
    ## ... Procrustes: rmse 0.03592549  max resid 0.1182485 
    ## Run 146 stress 0.0922849 
    ## Run 147 stress 0.0892607 
    ## ... Procrustes: rmse 0.000161049  max resid 0.0005471734 
    ## ... Similar to previous best
    ## Run 148 stress 0.1067588 
    ## Run 149 stress 0.08938963 
    ## ... Procrustes: rmse 0.03595806  max resid 0.1182901 
    ## Run 150 stress 0.08938539 
    ## ... Procrustes: rmse 0.01184579  max resid 0.03972384 
    ## Run 151 stress 0.08926065 
    ## ... New best solution
    ## ... Procrustes: rmse 8.255546e-05  max resid 0.000242902 
    ## ... Similar to previous best
    ## Run 152 stress 0.08926068 
    ## ... Procrustes: rmse 0.0002232357  max resid 0.0005605565 
    ## ... Similar to previous best
    ## Run 153 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183296  max resid 0.03971891 
    ## Run 154 stress 0.08938539 
    ## ... Procrustes: rmse 0.01184258  max resid 0.03973276 
    ## Run 155 stress 0.0903913 
    ## Run 156 stress 0.1052648 
    ## Run 157 stress 0.08926064 
    ## ... New best solution
    ## ... Procrustes: rmse 2.089444e-05  max resid 6.073691e-05 
    ## ... Similar to previous best
    ## Run 158 stress 0.09039131 
    ## Run 159 stress 0.1120099 
    ## Run 160 stress 0.1108634 
    ## Run 161 stress 0.1060433 
    ## Run 162 stress 0.0894665 
    ## ... Procrustes: rmse 0.0332282  max resid 0.1163586 
    ## Run 163 stress 0.1060429 
    ## Run 164 stress 0.08946658 
    ## ... Procrustes: rmse 0.03320401  max resid 0.1163168 
    ## Run 165 stress 0.1075421 
    ## Run 166 stress 0.08926067 
    ## ... Procrustes: rmse 2.895569e-05  max resid 7.607549e-05 
    ## ... Similar to previous best
    ## Run 167 stress 0.09039136 
    ## Run 168 stress 0.08938537 
    ## ... Procrustes: rmse 0.0118376  max resid 0.03969419 
    ## Run 169 stress 0.08938964 
    ## ... Procrustes: rmse 0.03594635  max resid 0.1182808 
    ## Run 170 stress 0.105265 
    ## Run 171 stress 0.08946339 
    ## ... Procrustes: rmse 0.03748723  max resid 0.1177271 
    ## Run 172 stress 0.08938968 
    ## ... Procrustes: rmse 0.03598458  max resid 0.1183489 
    ## Run 173 stress 0.1096134 
    ## Run 174 stress 0.09503456 
    ## Run 175 stress 0.09503432 
    ## Run 176 stress 0.08926091 
    ## ... Procrustes: rmse 0.0002766261  max resid 0.0008809682 
    ## ... Similar to previous best
    ## Run 177 stress 0.1063191 
    ## Run 178 stress 0.09592141 
    ## Run 179 stress 0.08926074 
    ## ... Procrustes: rmse 0.0001495768  max resid 0.0005013345 
    ## ... Similar to previous best
    ## Run 180 stress 0.09039129 
    ## Run 181 stress 0.1075421 
    ## Run 182 stress 0.1075423 
    ## Run 183 stress 0.09775546 
    ## Run 184 stress 0.1091448 
    ## Run 185 stress 0.08926065 
    ## ... Procrustes: rmse 0.0001363414  max resid 0.0003370252 
    ## ... Similar to previous best
    ## Run 186 stress 0.09044608 
    ## Run 187 stress 0.1060094 
    ## Run 188 stress 0.1067485 
    ## Run 189 stress 0.09612798 
    ## Run 190 stress 0.09099587 
    ## Run 191 stress 0.08938973 
    ## ... Procrustes: rmse 0.03592639  max resid 0.1182506 
    ## Run 192 stress 0.1071322 
    ## Run 193 stress 0.08938964 
    ## ... Procrustes: rmse 0.03594277  max resid 0.1182727 
    ## Run 194 stress 0.08938537 
    ## ... Procrustes: rmse 0.01182665  max resid 0.03968874 
    ## Run 195 stress 0.09018707 
    ## Run 196 stress 0.1111382 
    ## Run 197 stress 0.0893854 
    ## ... Procrustes: rmse 0.01181217  max resid 0.03969938 
    ## Run 198 stress 0.08926068 
    ## ... Procrustes: rmse 0.0002231783  max resid 0.0006846782 
    ## ... Similar to previous best
    ## Run 199 stress 0.1104472 
    ## Run 200 stress 0.1056904 
    ## Run 201 stress 0.08951716 
    ## ... Procrustes: rmse 0.0350755  max resid 0.1156723 
    ## Run 202 stress 0.09228483 
    ## Run 203 stress 0.08946653 
    ## ... Procrustes: rmse 0.03321748  max resid 0.1163392 
    ## Run 204 stress 0.09099545 
    ## Run 205 stress 0.08946652 
    ## ... Procrustes: rmse 0.0332443  max resid 0.1163979 
    ## Run 206 stress 0.08946332 
    ## ... Procrustes: rmse 0.03750218  max resid 0.1177519 
    ## Run 207 stress 0.109213 
    ## Run 208 stress 0.08926087 
    ## ... Procrustes: rmse 0.000253636  max resid 0.0008085258 
    ## ... Similar to previous best
    ## Run 209 stress 0.1089903 
    ## Run 210 stress 0.09021168 
    ## Run 211 stress 0.1056893 
    ## Run 212 stress 0.1060436 
    ## Run 213 stress 0.09088605 
    ## Run 214 stress 0.1062563 
    ## Run 215 stress 0.1052651 
    ## Run 216 stress 0.08938962 
    ## ... Procrustes: rmse 0.03594801  max resid 0.1182807 
    ## Run 217 stress 0.08938963 
    ## ... Procrustes: rmse 0.03595182  max resid 0.1182844 
    ## Run 218 stress 0.1071029 
    ## Run 219 stress 0.08951717 
    ## ... Procrustes: rmse 0.03506547  max resid 0.1156541 
    ## Run 220 stress 0.1056893 
    ## Run 221 stress 0.08926082 
    ## ... Procrustes: rmse 0.0002187328  max resid 0.0007023457 
    ## ... Similar to previous best
    ## Run 222 stress 0.1091828 
    ## Run 223 stress 0.08926084 
    ## ... Procrustes: rmse 0.0002163559  max resid 0.0004640879 
    ## ... Similar to previous best
    ## Run 224 stress 0.1056896 
    ## Run 225 stress 0.1108713 
    ## Run 226 stress 0.1092128 
    ## Run 227 stress 0.08926079 
    ## ... Procrustes: rmse 0.000331693  max resid 0.001020788 
    ## ... Similar to previous best
    ## Run 228 stress 0.09039137 
    ## Run 229 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001046582  max resid 0.0002640455 
    ## ... Similar to previous best
    ## Run 230 stress 0.08926091 
    ## ... Procrustes: rmse 0.0003764073  max resid 0.001178186 
    ## ... Similar to previous best
    ## Run 231 stress 0.08946328 
    ## ... Procrustes: rmse 0.0375095  max resid 0.1177751 
    ## Run 232 stress 0.09592125 
    ## Run 233 stress 0.08938976 
    ## ... Procrustes: rmse 0.03590162  max resid 0.1182174 
    ## Run 234 stress 0.1075421 
    ## Run 235 stress 0.09099542 
    ## Run 236 stress 0.09109113 
    ## Run 237 stress 0.1075421 
    ## Run 238 stress 0.09503418 
    ## Run 239 stress 0.1071321 
    ## Run 240 stress 0.09130092 
    ## Run 241 stress 0.08938552 
    ## ... Procrustes: rmse 0.01177326  max resid 0.03949919 
    ## Run 242 stress 0.08938962 
    ## ... Procrustes: rmse 0.03594759  max resid 0.118287 
    ## Run 243 stress 0.09109108 
    ## Run 244 stress 0.1089198 
    ## Run 245 stress 0.1056903 
    ## Run 246 stress 0.08946654 
    ## ... Procrustes: rmse 0.03322409  max resid 0.1163852 
    ## Run 247 stress 0.08946341 
    ## ... Procrustes: rmse 0.03746135  max resid 0.1177025 
    ## Run 248 stress 0.09088623 
    ## Run 249 stress 0.09109108 
    ## Run 250 stress 0.08938966 
    ## ... Procrustes: rmse 0.03591612  max resid 0.1182374 
    ## Run 251 stress 0.1105344 
    ## Run 252 stress 0.09503436 
    ## Run 253 stress 0.1096135 
    ## Run 254 stress 0.08938554 
    ## ... Procrustes: rmse 0.01175881  max resid 0.03966183 
    ## Run 255 stress 0.09503441 
    ## Run 256 stress 0.09039129 
    ## Run 257 stress 0.08938539 
    ## ... Procrustes: rmse 0.01180783  max resid 0.0396277 
    ## Run 258 stress 0.09503453 
    ## Run 259 stress 0.08926069 
    ## ... Procrustes: rmse 0.000169942  max resid 0.0005478728 
    ## ... Similar to previous best
    ## Run 260 stress 0.08946653 
    ## ... Procrustes: rmse 0.03322276  max resid 0.1163747 
    ## Run 261 stress 0.08938968 
    ## ... Procrustes: rmse 0.03591209  max resid 0.1182308 
    ## Run 262 stress 0.08926095 
    ## ... Procrustes: rmse 0.0003976539  max resid 0.001247836 
    ## ... Similar to previous best
    ## Run 263 stress 0.08946331 
    ## ... Procrustes: rmse 0.03748412  max resid 0.1177327 
    ## Run 264 stress 0.08926066 
    ## ... Procrustes: rmse 0.0001128408  max resid 0.0002757166 
    ## ... Similar to previous best
    ## Run 265 stress 0.1092094 
    ## Run 266 stress 0.1075423 
    ## Run 267 stress 0.1080636 
    ## Run 268 stress 0.1064197 
    ## Run 269 stress 0.08946331 
    ## ... Procrustes: rmse 0.03750769  max resid 0.1177625 
    ## Run 270 stress 0.1074432 
    ## Run 271 stress 0.1075423 
    ## Run 272 stress 0.1061307 
    ## Run 273 stress 0.1071023 
    ## Run 274 stress 0.1071322 
    ## Run 275 stress 0.1110405 
    ## Run 276 stress 0.08926078 
    ## ... Procrustes: rmse 0.0002259196  max resid 0.0007195448 
    ## ... Similar to previous best
    ## Run 277 stress 0.08946333 
    ## ... Procrustes: rmse 0.03747717  max resid 0.1177218 
    ## Run 278 stress 0.1067507 
    ## Run 279 stress 0.09109113 
    ## Run 280 stress 0.0950342 
    ## Run 281 stress 0.1075421 
    ## Run 282 stress 0.08938544 
    ## ... Procrustes: rmse 0.01184987  max resid 0.03966462 
    ## Run 283 stress 0.1075423 
    ## Run 284 stress 0.09503441 
    ## Run 285 stress 0.1064203 
    ## Run 286 stress 0.1065375 
    ## Run 287 stress 0.08926088 
    ## ... Procrustes: rmse 0.0003612616  max resid 0.001113848 
    ## ... Similar to previous best
    ## Run 288 stress 0.09039147 
    ## Run 289 stress 0.09018708 
    ## Run 290 stress 0.0892611 
    ## ... Procrustes: rmse 0.0004808362  max resid 0.001449913 
    ## ... Similar to previous best
    ## Run 291 stress 0.1065365 
    ## Run 292 stress 0.1052649 
    ## Run 293 stress 0.08946652 
    ## ... Procrustes: rmse 0.03319454  max resid 0.1163153 
    ## Run 294 stress 0.08926081 
    ## ... Procrustes: rmse 0.000300644  max resid 0.0009704602 
    ## ... Similar to previous best
    ## Run 295 stress 0.09228503 
    ## Run 296 stress 0.1074433 
    ## Run 297 stress 0.08946657 
    ## ... Procrustes: rmse 0.03318245  max resid 0.1162934 
    ## Run 298 stress 0.09612762 
    ## Run 299 stress 0.08926069 
    ## ... Procrustes: rmse 7.993463e-05  max resid 0.0002216281 
    ## ... Similar to previous best
    ## Run 300 stress 0.1063193 
    ## Run 301 stress 0.1074434 
    ## Run 302 stress 0.1092129 
    ## Run 303 stress 0.1103709 
    ## Run 304 stress 0.09503428 
    ## Run 305 stress 0.1052647 
    ## Run 306 stress 0.10569 
    ## Run 307 stress 0.09039153 
    ## Run 308 stress 0.0959215 
    ## Run 309 stress 0.1074432 
    ## Run 310 stress 0.1052647 
    ## Run 311 stress 0.08938541 
    ## ... Procrustes: rmse 0.01178812  max resid 0.039634 
    ## Run 312 stress 0.110845 
    ## Run 313 stress 0.1071322 
    ## Run 314 stress 0.1091829 
    ## Run 315 stress 0.0959214 
    ## Run 316 stress 0.08938962 
    ## ... Procrustes: rmse 0.03595251  max resid 0.118296 
    ## Run 317 stress 0.09099569 
    ## Run 318 stress 0.1052648 
    ## Run 319 stress 0.08951716 
    ## ... Procrustes: rmse 0.03503886  max resid 0.1156247 
    ## Run 320 stress 0.09088601 
    ## Run 321 stress 0.1092092 
    ## Run 322 stress 0.1062569 
    ## Run 323 stress 0.08938962 
    ## ... Procrustes: rmse 0.03593056  max resid 0.1182594 
    ## Run 324 stress 0.1087868 
    ## Run 325 stress 0.08938965 
    ## ... Procrustes: rmse 0.03595487  max resid 0.118301 
    ## Run 326 stress 0.09021179 
    ## Run 327 stress 0.09503422 
    ## Run 328 stress 0.09503418 
    ## Run 329 stress 0.1074432 
    ## Run 330 stress 0.1076307 
    ## Run 331 stress 0.08938979 
    ## ... Procrustes: rmse 0.03589625  max resid 0.1182066 
    ## Run 332 stress 0.08926083 
    ## ... Procrustes: rmse 0.0002692667  max resid 0.0007914578 
    ## ... Similar to previous best
    ## Run 333 stress 0.09039137 
    ## Run 334 stress 0.1052649 
    ## Run 335 stress 0.1060094 
    ## Run 336 stress 0.109221 
    ## Run 337 stress 0.1112204 
    ## Run 338 stress 0.09039131 
    ## Run 339 stress 0.1052647 
    ## Run 340 stress 0.09109108 
    ## Run 341 stress 0.108656 
    ## Run 342 stress 0.1085125 
    ## Run 343 stress 0.08938538 
    ## ... Procrustes: rmse 0.01180323  max resid 0.0396134 
    ## Run 344 stress 0.08938963 
    ## ... Procrustes: rmse 0.03592613  max resid 0.1182548 
    ## Run 345 stress 0.08926072 
    ## ... Procrustes: rmse 0.0002358095  max resid 0.0006970573 
    ## ... Similar to previous best
    ## Run 346 stress 0.09087356 
    ## Run 347 stress 0.1064427 
    ## Run 348 stress 0.09109111 
    ## Run 349 stress 0.08938975 
    ## ... Procrustes: rmse 0.03590157  max resid 0.1182132 
    ## Run 350 stress 0.09503417 
    ## Run 351 stress 0.1092095 
    ## Run 352 stress 0.08926077 
    ## ... Procrustes: rmse 0.0002556923  max resid 0.000615078 
    ## ... Similar to previous best
    ## Run 353 stress 0.1068614 
    ## Run 354 stress 0.106737 
    ## Run 355 stress 0.1075421 
    ## Run 356 stress 0.09044631 
    ## Run 357 stress 0.1060435 
    ## Run 358 stress 0.09503436 
    ## Run 359 stress 0.09592169 
    ## Run 360 stress 0.1065377 
    ## Run 361 stress 0.08946656 
    ## ... Procrustes: rmse 0.03318185  max resid 0.1162942 
    ## Run 362 stress 0.1091505 
    ## Run 363 stress 0.10569 
    ## Run 364 stress 0.09044609 
    ## Run 365 stress 0.08938539 
    ## ... Procrustes: rmse 0.01182901  max resid 0.0396377 
    ## Run 366 stress 0.09178289 
    ## Run 367 stress 0.089261 
    ## ... Procrustes: rmse 0.0004304929  max resid 0.001294112 
    ## ... Similar to previous best
    ## Run 368 stress 0.1076303 
    ## Run 369 stress 0.1074432 
    ## Run 370 stress 0.1052648 
    ## Run 371 stress 0.09180879 
    ## Run 372 stress 0.1076303 
    ## Run 373 stress 0.1085126 
    ## Run 374 stress 0.0893898 
    ## ... Procrustes: rmse 0.03598551  max resid 0.1183439 
    ## Run 375 stress 0.1092361 
    ## Run 376 stress 0.1061296 
    ## Run 377 stress 0.08926065 
    ## ... Procrustes: rmse 0.0001153803  max resid 0.000373329 
    ## ... Similar to previous best
    ## Run 378 stress 0.1071322 
    ## Run 379 stress 0.08946329 
    ## ... Procrustes: rmse 0.03749206  max resid 0.1177455 
    ## Run 380 stress 0.08951717 
    ## ... Procrustes: rmse 0.03504895  max resid 0.115629 
    ## Run 381 stress 0.08938964 
    ## ... Procrustes: rmse 0.03595905  max resid 0.1183076 
    ## Run 382 stress 0.08938538 
    ## ... Procrustes: rmse 0.01181401  max resid 0.0396331 
    ## Run 383 stress 0.08938961 
    ## ... Procrustes: rmse 0.03594809  max resid 0.1182858 
    ## Run 384 stress 0.08926089 
    ## ... Procrustes: rmse 0.0003203626  max resid 0.0009479839 
    ## ... Similar to previous best
    ## Run 385 stress 0.08946333 
    ## ... Procrustes: rmse 0.03748091  max resid 0.1177261 
    ## Run 386 stress 0.10632 
    ## Run 387 stress 0.08926093 
    ## ... Procrustes: rmse 0.0003963006  max resid 0.001148825 
    ## ... Similar to previous best
    ## Run 388 stress 0.08946651 
    ## ... Procrustes: rmse 0.03321653  max resid 0.1163595 
    ## Run 389 stress 0.08926074 
    ## ... Procrustes: rmse 0.0001950034  max resid 0.0006517973 
    ## ... Similar to previous best
    ## Run 390 stress 0.1071322 
    ## Run 391 stress 0.08938542 
    ## ... Procrustes: rmse 0.01186717  max resid 0.03962141 
    ## Run 392 stress 0.08926091 
    ## ... Procrustes: rmse 0.0003799003  max resid 0.001183368 
    ## ... Similar to previous best
    ## Run 393 stress 0.1076304 
    ## Run 394 stress 0.08938965 
    ## ... Procrustes: rmse 0.03595994  max resid 0.1183123 
    ## Run 395 stress 0.1071322 
    ## Run 396 stress 0.08926077 
    ## ... Procrustes: rmse 0.0002452343  max resid 0.0007937176 
    ## ... Similar to previous best
    ## Run 397 stress 0.08946669 
    ## ... Procrustes: rmse 0.03316047  max resid 0.1162463 
    ## Run 398 stress 0.09099557 
    ## Run 399 stress 0.0913009 
    ## Run 400 stress 0.1111369 
    ## Run 401 stress 0.1118954 
    ## Run 402 stress 0.1108457 
    ## Run 403 stress 0.1060433 
    ## Run 404 stress 0.106043 
    ## Run 405 stress 0.08938961 
    ## ... Procrustes: rmse 0.03594378  max resid 0.1182774 
    ## Run 406 stress 0.09503427 
    ## Run 407 stress 0.109561 
    ## Run 408 stress 0.08938561 
    ## ... Procrustes: rmse 0.01175137  max resid 0.0396755 
    ## Run 409 stress 0.08926072 
    ## ... Procrustes: rmse 0.0002270771  max resid 0.0006500895 
    ## ... Similar to previous best
    ## Run 410 stress 0.1101451 
    ## Run 411 stress 0.1092208 
    ## Run 412 stress 0.105265 
    ## Run 413 stress 0.08938983 
    ## ... Procrustes: rmse 0.03589189  max resid 0.1181974 
    ## Run 414 stress 0.0892608 
    ## ... Procrustes: rmse 0.0002952849  max resid 0.0008725557 
    ## ... Similar to previous best
    ## Run 415 stress 0.08938538 
    ## ... Procrustes: rmse 0.01181466  max resid 0.0396123 
    ## Run 416 stress 0.1074432 
    ## Run 417 stress 0.08938963 
    ## ... Procrustes: rmse 0.03595553  max resid 0.1183013 
    ## Run 418 stress 0.08946328 
    ## ... Procrustes: rmse 0.03749966  max resid 0.1177586 
    ## Run 419 stress 0.08938546 
    ## ... Procrustes: rmse 0.01185087  max resid 0.03967762 
    ## Run 420 stress 0.09178289 
    ## Run 421 stress 0.09088607 
    ## Run 422 stress 0.08938975 
    ## ... Procrustes: rmse 0.03590236  max resid 0.1182156 
    ## Run 423 stress 0.08926072 
    ## ... Procrustes: rmse 0.0002345028  max resid 0.0006880819 
    ## ... Similar to previous best
    ## Run 424 stress 0.1103732 
    ## Run 425 stress 0.1104217 
    ## Run 426 stress 0.1068683 
    ## Run 427 stress 0.1071321 
    ## Run 428 stress 0.08938967 
    ## ... Procrustes: rmse 0.03596373  max resid 0.118322 
    ## Run 429 stress 0.08926065 
    ## ... Procrustes: rmse 8.705286e-05  max resid 0.0002269863 
    ## ... Similar to previous best
    ## Run 430 stress 0.1060431 
    ## Run 431 stress 0.09594359 
    ## Run 432 stress 0.08926082 
    ## ... Procrustes: rmse 0.0003140343  max resid 0.0009556614 
    ## ... Similar to previous best
    ## Run 433 stress 0.08938537 
    ## ... Procrustes: rmse 0.01181738  max resid 0.03961347 
    ## Run 434 stress 0.08926078 
    ## ... Procrustes: rmse 0.0002222761  max resid 0.0007129351 
    ## ... Similar to previous best
    ## Run 435 stress 0.1091826 
    ## Run 436 stress 0.1092128 
    ## Run 437 stress 0.1096133 
    ## Run 438 stress 0.1092562 
    ## Run 439 stress 0.1071322 
    ## Run 440 stress 0.1088392 
    ## Run 441 stress 0.09099549 
    ## Run 442 stress 0.1056904 
    ## Run 443 stress 0.08938976 
    ## ... Procrustes: rmse 0.03598164  max resid 0.1183408 
    ## Run 444 stress 0.08938965 
    ## ... Procrustes: rmse 0.03592446  max resid 0.1182532 
    ## Run 445 stress 0.09099544 
    ## Run 446 stress 0.08938969 
    ## ... Procrustes: rmse 0.03591082  max resid 0.1182303 
    ## Run 447 stress 0.09592146 
    ## Run 448 stress 0.1111368 
    ## Run 449 stress 0.09178282 
    ## Run 450 stress 0.08938553 
    ## ... Procrustes: rmse 0.0118775  max resid 0.03957892 
    ## Run 451 stress 0.08946332 
    ## ... Procrustes: rmse 0.0374767  max resid 0.1177191 
    ## Run 452 stress 0.1101448 
    ## Run 453 stress 0.09594347 
    ## Run 454 stress 0.08951722 
    ## ... Procrustes: rmse 0.03501488  max resid 0.1155987 
    ## Run 455 stress 0.08926085 
    ## ... Procrustes: rmse 0.0003401581  max resid 0.001050075 
    ## ... Similar to previous best
    ## Run 456 stress 0.08938974 
    ## ... Procrustes: rmse 0.03590304  max resid 0.1182166 
    ## Run 457 stress 0.08946653 
    ## ... Procrustes: rmse 0.03319281  max resid 0.1163111 
    ## Run 458 stress 0.08946663 
    ## ... Procrustes: rmse 0.03317032  max resid 0.1162779 
    ## Run 459 stress 0.1052647 
    ## Run 460 stress 0.1060429 
    ## Run 461 stress 0.1074432 
    ## Run 462 stress 0.09088601 
    ## Run 463 stress 0.09044615 
    ## Run 464 stress 0.1062594 
    ## Run 465 stress 0.08946654 
    ## ... Procrustes: rmse 0.03318669  max resid 0.1163006 
    ## Run 466 stress 0.09592149 
    ## Run 467 stress 0.09018707 
    ## Run 468 stress 0.08926089 
    ## ... Procrustes: rmse 0.0003661138  max resid 0.001133328 
    ## ... Similar to previous best
    ## Run 469 stress 0.1105347 
    ## Run 470 stress 0.09178291 
    ## Run 471 stress 0.08946328 
    ## ... Procrustes: rmse 0.03750243  max resid 0.1177642 
    ## Run 472 stress 0.08938966 
    ## ... Procrustes: rmse 0.03596096  max resid 0.1183126 
    ## Run 473 stress 0.1091824 
    ## Run 474 stress 0.0893854 
    ## ... Procrustes: rmse 0.01182898  max resid 0.03958066 
    ## Run 475 stress 0.08926117 
    ## ... Procrustes: rmse 0.0004626283  max resid 0.001469668 
    ## ... Similar to previous best
    ## Run 476 stress 0.1092092 
    ## Run 477 stress 0.09503423 
    ## Run 478 stress 0.08938542 
    ## ... Procrustes: rmse 0.01176817  max resid 0.03957287 
    ## Run 479 stress 0.08926075 
    ## ... Procrustes: rmse 0.0002647407  max resid 0.0007913387 
    ## ... Similar to previous best
    ## Run 480 stress 0.1112208 
    ## Run 481 stress 0.08938548 
    ## ... Procrustes: rmse 0.01177152  max resid 0.03965359 
    ## Run 482 stress 0.1074434 
    ## Run 483 stress 0.08938538 
    ## ... Procrustes: rmse 0.01181795  max resid 0.03963572 
    ## Run 484 stress 0.08938976 
    ## ... Procrustes: rmse 0.03590056  max resid 0.118216 
    ## Run 485 stress 0.08938967 
    ## ... Procrustes: rmse 0.03596303  max resid 0.1183213 
    ## Run 486 stress 0.1060429 
    ## Run 487 stress 0.09109111 
    ## Run 488 stress 0.1052652 
    ## Run 489 stress 0.09503452 
    ## Run 490 stress 0.1101446 
    ## Run 491 stress 0.1128304 
    ## Run 492 stress 0.09109113 
    ## Run 493 stress 0.09109115 
    ## Run 494 stress 0.105265 
    ## Run 495 stress 0.08938962 
    ## ... Procrustes: rmse 0.03593015  max resid 0.1182581 
    ## Run 496 stress 0.0892607 
    ## ... Procrustes: rmse 0.000206596  max resid 0.000559921 
    ## ... Similar to previous best
    ## Run 497 stress 0.08946329 
    ## ... Procrustes: rmse 0.03748781  max resid 0.1177395 
    ## Run 498 stress 0.08946331 
    ## ... Procrustes: rmse 0.03748088  max resid 0.117731 
    ## Run 499 stress 0.08926096 
    ## ... Procrustes: rmse 0.000411346  max resid 0.001236581 
    ## ... Similar to previous best
    ## Run 500 stress 0.1071322 
    ## *** Best solution repeated 32 times

``` r
# Mixed and stratified lakes
SD_beta_geo_MS_NMDS <- metaMDS(SD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09407982 
    ## Run 2 stress 0.0914534 
    ## Run 3 stress 0.09168936 
    ## Run 4 stress 0.09337277 
    ## Run 5 stress 0.09268397 
    ## Run 6 stress 0.09286083 
    ## Run 7 stress 0.09407983 
    ## Run 8 stress 0.09408013 
    ## Run 9 stress 0.09610753 
    ## Run 10 stress 0.09969987 
    ## Run 11 stress 0.09969971 
    ## Run 12 stress 0.3254553 
    ## Run 13 stress 0.09159087 
    ## Run 14 stress 0.08440255 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01443197  max resid 0.04297608 
    ## Run 15 stress 0.09337233 
    ## Run 16 stress 0.09407979 
    ## Run 17 stress 0.08440254 
    ## ... New best solution
    ## ... Procrustes: rmse 2.642263e-05  max resid 5.019977e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.09445499 
    ## Run 19 stress 0.09159102 
    ## Run 20 stress 0.09159093 
    ## Run 21 stress 0.09408004 
    ## Run 22 stress 0.08503497 
    ## Run 23 stress 0.09969991 
    ## Run 24 stress 0.08973867 
    ## Run 25 stress 0.09308951 
    ## Run 26 stress 0.08503627 
    ## Run 27 stress 0.09168934 
    ## Run 28 stress 0.09416141 
    ## Run 29 stress 0.09159084 
    ## Run 30 stress 0.08440266 
    ## ... Procrustes: rmse 0.0001470759  max resid 0.0003161458 
    ## ... Similar to previous best
    ## Run 31 stress 0.09407986 
    ## Run 32 stress 0.08503488 
    ## Run 33 stress 0.09760825 
    ## Run 34 stress 0.08503509 
    ## Run 35 stress 0.08973885 
    ## Run 36 stress 0.08503538 
    ## Run 37 stress 0.09381844 
    ## Run 38 stress 0.09465905 
    ## Run 39 stress 0.1113498 
    ## Run 40 stress 0.09969964 
    ## Run 41 stress 0.08503627 
    ## Run 42 stress 0.09308961 
    ## Run 43 stress 0.09308955 
    ## Run 44 stress 0.09268343 
    ## Run 45 stress 0.09381852 
    ## Run 46 stress 0.2643453 
    ## Run 47 stress 0.09465905 
    ## Run 48 stress 0.2511969 
    ## Run 49 stress 0.08773476 
    ## Run 50 stress 0.09159083 
    ## Run 51 stress 0.09407992 
    ## Run 52 stress 0.09465909 
    ## Run 53 stress 0.08503455 
    ## Run 54 stress 0.08973876 
    ## Run 55 stress 0.09590946 
    ## Run 56 stress 0.09286086 
    ## Run 57 stress 0.09159094 
    ## Run 58 stress 0.09374352 
    ## Run 59 stress 0.09374278 
    ## Run 60 stress 0.09168955 
    ## Run 61 stress 0.09407968 
    ## Run 62 stress 0.0850347 
    ## Run 63 stress 0.0938059 
    ## Run 64 stress 0.08440255 
    ## ... Procrustes: rmse 4.824514e-05  max resid 7.95086e-05 
    ## ... Similar to previous best
    ## Run 65 stress 0.1038036 
    ## Run 66 stress 0.08973871 
    ## Run 67 stress 0.08503653 
    ## Run 68 stress 0.09374368 
    ## Run 69 stress 0.090304 
    ## Run 70 stress 0.09416153 
    ## Run 71 stress 0.08773467 
    ## Run 72 stress 0.08503573 
    ## Run 73 stress 0.09145328 
    ## Run 74 stress 0.09286092 
    ## Run 75 stress 0.09969963 
    ## Run 76 stress 0.08503672 
    ## Run 77 stress 0.09159087 
    ## Run 78 stress 0.09408008 
    ## Run 79 stress 0.08503496 
    ## Run 80 stress 0.08503569 
    ## Run 81 stress 0.09268362 
    ## Run 82 stress 0.09760831 
    ## Run 83 stress 0.09337247 
    ## Run 84 stress 0.09268372 
    ## Run 85 stress 0.0996997 
    ## Run 86 stress 0.09465914 
    ## Run 87 stress 0.09286089 
    ## Run 88 stress 0.0940799 
    ## Run 89 stress 0.09407104 
    ## Run 90 stress 0.08973863 
    ## Run 91 stress 0.09268377 
    ## Run 92 stress 0.08440258 
    ## ... Procrustes: rmse 7.010151e-05  max resid 0.0001699302 
    ## ... Similar to previous best
    ## Run 93 stress 0.09030402 
    ## Run 94 stress 0.08773482 
    ## Run 95 stress 0.0959086 
    ## Run 96 stress 0.09159083 
    ## Run 97 stress 0.08503511 
    ## Run 98 stress 0.09286093 
    ## Run 99 stress 0.09760844 
    ## Run 100 stress 0.09030402 
    ## Run 101 stress 0.09539199 
    ## Run 102 stress 0.08440277 
    ## ... Procrustes: rmse 0.0001888545  max resid 0.0003581724 
    ## ... Similar to previous best
    ## Run 103 stress 0.09407986 
    ## Run 104 stress 0.2465478 
    ## Run 105 stress 0.09286085 
    ## Run 106 stress 0.09535522 
    ## Run 107 stress 0.09445465 
    ## Run 108 stress 0.09407974 
    ## Run 109 stress 0.09712971 
    ## Run 110 stress 0.09407986 
    ## Run 111 stress 0.09268348 
    ## Run 112 stress 0.09447319 
    ## Run 113 stress 0.08440255 
    ## ... Procrustes: rmse 1.805666e-05  max resid 3.529711e-05 
    ## ... Similar to previous best
    ## Run 114 stress 0.09760813 
    ## Run 115 stress 0.08503614 
    ## Run 116 stress 0.09374248 
    ## Run 117 stress 0.08503463 
    ## Run 118 stress 0.33558 
    ## Run 119 stress 0.0940798 
    ## Run 120 stress 0.09400534 
    ## Run 121 stress 0.09030401 
    ## Run 122 stress 0.29859 
    ## Run 123 stress 0.09030418 
    ## Run 124 stress 0.08503654 
    ## Run 125 stress 0.09407989 
    ## Run 126 stress 0.0926836 
    ## Run 127 stress 0.09400521 
    ## Run 128 stress 0.1038036 
    ## Run 129 stress 0.09308963 
    ## Run 130 stress 0.08773483 
    ## Run 131 stress 0.0941615 
    ## Run 132 stress 0.09168928 
    ## Run 133 stress 0.08440254 
    ## ... New best solution
    ## ... Procrustes: rmse 4.477575e-05  max resid 7.961903e-05 
    ## ... Similar to previous best
    ## Run 134 stress 0.09030393 
    ## Run 135 stress 0.08503508 
    ## Run 136 stress 0.08773475 
    ## Run 137 stress 0.09447333 
    ## Run 138 stress 0.089739 
    ## Run 139 stress 0.09337253 
    ## Run 140 stress 0.09159083 
    ## Run 141 stress 0.09337265 
    ## Run 142 stress 0.1052849 
    ## Run 143 stress 0.09465905 
    ## Run 144 stress 0.09308972 
    ## Run 145 stress 0.09030408 
    ## Run 146 stress 0.0903041 
    ## Run 147 stress 0.09168952 
    ## Run 148 stress 0.09407996 
    ## Run 149 stress 0.08440258 
    ## ... Procrustes: rmse 7.907265e-05  max resid 0.0001469634 
    ## ... Similar to previous best
    ## Run 150 stress 0.09268368 
    ## Run 151 stress 0.09590908 
    ## Run 152 stress 0.09145321 
    ## Run 153 stress 0.09337263 
    ## Run 154 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 8.261646e-05  max resid 0.0001582433 
    ## ... Similar to previous best
    ## Run 155 stress 0.08773467 
    ## Run 156 stress 0.0996997 
    ## Run 157 stress 0.08773479 
    ## Run 158 stress 0.08440265 
    ## ... Procrustes: rmse 0.0002207107  max resid 0.0004876567 
    ## ... Similar to previous best
    ## Run 159 stress 0.09308948 
    ## Run 160 stress 0.08503464 
    ## Run 161 stress 0.08440254 
    ## ... Procrustes: rmse 0.0001141048  max resid 0.0002196762 
    ## ... Similar to previous best
    ## Run 162 stress 0.09400532 
    ## Run 163 stress 0.09760831 
    ## Run 164 stress 0.09030394 
    ## Run 165 stress 0.08440254 
    ## ... Procrustes: rmse 7.696911e-05  max resid 0.0001609979 
    ## ... Similar to previous best
    ## Run 166 stress 0.09408002 
    ## Run 167 stress 0.09030399 
    ## Run 168 stress 0.09969981 
    ## Run 169 stress 0.09030401 
    ## Run 170 stress 0.08503488 
    ## Run 171 stress 0.1038042 
    ## Run 172 stress 0.09407976 
    ## Run 173 stress 0.3138878 
    ## Run 174 stress 0.09030393 
    ## Run 175 stress 0.08773466 
    ## Run 176 stress 0.08973876 
    ## Run 177 stress 0.09268407 
    ## Run 178 stress 0.1052851 
    ## Run 179 stress 0.09030393 
    ## Run 180 stress 0.09159084 
    ## Run 181 stress 0.08503549 
    ## Run 182 stress 0.09760823 
    ## Run 183 stress 0.08440267 
    ## ... Procrustes: rmse 0.0002315468  max resid 0.0003999384 
    ## ... Similar to previous best
    ## Run 184 stress 0.09407984 
    ## Run 185 stress 0.1097016 
    ## Run 186 stress 0.08503512 
    ## Run 187 stress 0.08773471 
    ## Run 188 stress 0.09030413 
    ## Run 189 stress 0.09465912 
    ## Run 190 stress 0.08503499 
    ## Run 191 stress 0.09447306 
    ## Run 192 stress 0.09321708 
    ## Run 193 stress 0.09268346 
    ## Run 194 stress 0.08440252 
    ## ... Procrustes: rmse 6.408984e-05  max resid 0.0001252988 
    ## ... Similar to previous best
    ## Run 195 stress 0.09447343 
    ## Run 196 stress 0.09374386 
    ## Run 197 stress 0.0916893 
    ## Run 198 stress 0.09465921 
    ## Run 199 stress 0.09969979 
    ## Run 200 stress 0.09030402 
    ## Run 201 stress 0.09145323 
    ## Run 202 stress 0.09286087 
    ## Run 203 stress 0.09407981 
    ## Run 204 stress 0.0850354 
    ## Run 205 stress 0.09464431 
    ## Run 206 stress 0.09969996 
    ## Run 207 stress 0.09030406 
    ## Run 208 stress 0.09590947 
    ## Run 209 stress 0.09308949 
    ## Run 210 stress 0.09416136 
    ## Run 211 stress 0.09407118 
    ## Run 212 stress 0.09286093 
    ## Run 213 stress 0.09268359 
    ## Run 214 stress 0.08973875 
    ## Run 215 stress 0.09337299 
    ## Run 216 stress 0.09268368 
    ## Run 217 stress 0.08440262 
    ## ... Procrustes: rmse 0.0002002558  max resid 0.0003738022 
    ## ... Similar to previous best
    ## Run 218 stress 0.09268334 
    ## Run 219 stress 0.090304 
    ## Run 220 stress 0.0844028 
    ## ... Procrustes: rmse 0.0002849879  max resid 0.0005224125 
    ## ... Similar to previous best
    ## Run 221 stress 0.09721199 
    ## Run 222 stress 0.09268421 
    ## Run 223 stress 0.08503482 
    ## Run 224 stress 0.08503676 
    ## Run 225 stress 0.09416156 
    ## Run 226 stress 0.09308953 
    ## Run 227 stress 0.1038086 
    ## Run 228 stress 0.09417793 
    ## Run 229 stress 0.09721201 
    ## Run 230 stress 0.09145323 
    ## Run 231 stress 0.08773486 
    ## Run 232 stress 0.08503462 
    ## Run 233 stress 0.08440253 
    ## ... Procrustes: rmse 9.798481e-05  max resid 0.0001879543 
    ## ... Similar to previous best
    ## Run 234 stress 0.08973889 
    ## Run 235 stress 0.09168942 
    ## Run 236 stress 0.09977533 
    ## Run 237 stress 0.08973882 
    ## Run 238 stress 0.0850353 
    ## Run 239 stress 0.09168947 
    ## Run 240 stress 0.09159084 
    ## Run 241 stress 0.08503564 
    ## Run 242 stress 0.09721199 
    ## Run 243 stress 0.09145321 
    ## Run 244 stress 0.09159091 
    ## Run 245 stress 0.09030407 
    ## Run 246 stress 0.08503556 
    ## Run 247 stress 0.08440254 
    ## ... Procrustes: rmse 0.0001098645  max resid 0.0002050647 
    ## ... Similar to previous best
    ## Run 248 stress 0.0940799 
    ## Run 249 stress 0.08503463 
    ## Run 250 stress 0.08503513 
    ## Run 251 stress 0.09539199 
    ## Run 252 stress 0.0877347 
    ## Run 253 stress 0.08773476 
    ## Run 254 stress 0.09030397 
    ## Run 255 stress 0.1038041 
    ## Run 256 stress 0.09168949 
    ## Run 257 stress 0.08503511 
    ## Run 258 stress 0.08503464 
    ## Run 259 stress 0.09030413 
    ## Run 260 stress 0.09286087 
    ## Run 261 stress 0.08773467 
    ## Run 262 stress 0.09465905 
    ## Run 263 stress 0.09030396 
    ## Run 264 stress 0.08440253 
    ## ... Procrustes: rmse 9.824065e-05  max resid 0.0001999386 
    ## ... Similar to previous best
    ## Run 265 stress 0.09337249 
    ## Run 266 stress 0.09030403 
    ## Run 267 stress 0.08503676 
    ## Run 268 stress 0.329594 
    ## Run 269 stress 0.0850347 
    ## Run 270 stress 0.09381849 
    ## Run 271 stress 0.09417776 
    ## Run 272 stress 0.09268362 
    ## Run 273 stress 0.09407992 
    ## Run 274 stress 0.0928609 
    ## Run 275 stress 0.09535445 
    ## Run 276 stress 0.08503483 
    ## Run 277 stress 0.09159087 
    ## Run 278 stress 0.08773473 
    ## Run 279 stress 0.09030396 
    ## Run 280 stress 0.09030402 
    ## Run 281 stress 0.09380606 
    ## Run 282 stress 0.09321353 
    ## Run 283 stress 0.09416151 
    ## Run 284 stress 0.0941614 
    ## Run 285 stress 0.09286087 
    ## Run 286 stress 0.08503478 
    ## Run 287 stress 0.09286084 
    ## Run 288 stress 0.09408014 
    ## Run 289 stress 0.09268386 
    ## Run 290 stress 0.09760818 
    ## Run 291 stress 0.09407981 
    ## Run 292 stress 0.09408011 
    ## Run 293 stress 0.09404357 
    ## Run 294 stress 0.08440262 
    ## ... Procrustes: rmse 0.0001952383  max resid 0.0004159481 
    ## ... Similar to previous best
    ## Run 295 stress 0.08503482 
    ## Run 296 stress 0.08973883 
    ## Run 297 stress 0.08503506 
    ## Run 298 stress 0.09030393 
    ## Run 299 stress 0.09168949 
    ## Run 300 stress 0.09337238 
    ## Run 301 stress 0.09159084 
    ## Run 302 stress 0.08503511 
    ## Run 303 stress 0.09374252 
    ## Run 304 stress 0.09337286 
    ## Run 305 stress 0.09381848 
    ## Run 306 stress 0.09381859 
    ## Run 307 stress 0.08773479 
    ## Run 308 stress 0.09030398 
    ## Run 309 stress 0.08440266 
    ## ... Procrustes: rmse 0.0002262075  max resid 0.0004764004 
    ## ... Similar to previous best
    ## Run 310 stress 0.1052849 
    ## Run 311 stress 0.09381846 
    ## Run 312 stress 0.100461 
    ## Run 313 stress 0.09030405 
    ## Run 314 stress 0.09308962 
    ## Run 315 stress 0.09445501 
    ## Run 316 stress 0.09337248 
    ## Run 317 stress 0.09308952 
    ## Run 318 stress 0.09337279 
    ## Run 319 stress 0.09337287 
    ## Run 320 stress 0.09407992 
    ## Run 321 stress 0.09374307 
    ## Run 322 stress 0.08773469 
    ## Run 323 stress 0.09407368 
    ## Run 324 stress 0.09159083 
    ## Run 325 stress 0.09417784 
    ## Run 326 stress 0.09030404 
    ## Run 327 stress 0.09760827 
    ## Run 328 stress 0.08773479 
    ## Run 329 stress 0.09030397 
    ## Run 330 stress 0.08440252 
    ## ... Procrustes: rmse 6.386535e-05  max resid 0.0001242377 
    ## ... Similar to previous best
    ## Run 331 stress 0.09030393 
    ## Run 332 stress 0.09539201 
    ## Run 333 stress 0.09465904 
    ## Run 334 stress 0.09168944 
    ## Run 335 stress 0.09969962 
    ## Run 336 stress 0.08503465 
    ## Run 337 stress 0.08973875 
    ## Run 338 stress 0.09159083 
    ## Run 339 stress 0.09721207 
    ## Run 340 stress 0.08503516 
    ## Run 341 stress 0.08973876 
    ## Run 342 stress 0.09268375 
    ## Run 343 stress 0.09380582 
    ## Run 344 stress 0.09381872 
    ## Run 345 stress 0.09308952 
    ## Run 346 stress 0.09464429 
    ## Run 347 stress 0.1038045 
    ## Run 348 stress 0.09374337 
    ## Run 349 stress 0.09159084 
    ## Run 350 stress 0.1038042 
    ## Run 351 stress 0.09380575 
    ## Run 352 stress 0.08503488 
    ## Run 353 stress 0.09374291 
    ## Run 354 stress 0.09535509 
    ## Run 355 stress 0.08503652 
    ## Run 356 stress 0.09407969 
    ## Run 357 stress 0.1005996 
    ## Run 358 stress 0.09969962 
    ## Run 359 stress 0.09030405 
    ## Run 360 stress 0.09030396 
    ## Run 361 stress 0.09159094 
    ## Run 362 stress 0.09159085 
    ## Run 363 stress 0.08973887 
    ## Run 364 stress 0.09337271 
    ## Run 365 stress 0.08503487 
    ## Run 366 stress 0.09374266 
    ## Run 367 stress 0.08503499 
    ## Run 368 stress 0.09374305 
    ## Run 369 stress 0.09030398 
    ## Run 370 stress 0.0914533 
    ## Run 371 stress 0.0941614 
    ## Run 372 stress 0.09465904 
    ## Run 373 stress 0.3399 
    ## Run 374 stress 0.09400525 
    ## Run 375 stress 0.09168952 
    ## Run 376 stress 0.09286087 
    ## Run 377 stress 0.090304 
    ## Run 378 stress 0.09145321 
    ## Run 379 stress 0.09159083 
    ## Run 380 stress 0.09337261 
    ## Run 381 stress 0.0928609 
    ## Run 382 stress 0.09030394 
    ## Run 383 stress 0.09145321 
    ## Run 384 stress 0.09380582 
    ## Run 385 stress 0.09308985 
    ## Run 386 stress 0.09337257 
    ## Run 387 stress 0.09286087 
    ## Run 388 stress 0.09465906 
    ## Run 389 stress 0.09465908 
    ## Run 390 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001628818  max resid 0.0003018428 
    ## ... Similar to previous best
    ## Run 391 stress 0.0932178 
    ## Run 392 stress 0.08503603 
    ## Run 393 stress 0.09969987 
    ## Run 394 stress 0.0897388 
    ## Run 395 stress 0.09337229 
    ## Run 396 stress 0.0937429 
    ## Run 397 stress 0.2780284 
    ## Run 398 stress 0.09465906 
    ## Run 399 stress 0.09374208 
    ## Run 400 stress 0.09465904 
    ## Run 401 stress 0.09168926 
    ## Run 402 stress 0.0850361 
    ## Run 403 stress 0.09030406 
    ## Run 404 stress 0.08503545 
    ## Run 405 stress 0.09721196 
    ## Run 406 stress 0.09030401 
    ## Run 407 stress 0.1026216 
    ## Run 408 stress 0.08503477 
    ## Run 409 stress 0.09400535 
    ## Run 410 stress 0.1026216 
    ## Run 411 stress 0.08973868 
    ## Run 412 stress 0.09159085 
    ## Run 413 stress 0.09407967 
    ## Run 414 stress 0.09721202 
    ## Run 415 stress 0.0940799 
    ## Run 416 stress 0.09268404 
    ## Run 417 stress 0.09416152 
    ## Run 418 stress 0.09400534 
    ## Run 419 stress 0.08503467 
    ## Run 420 stress 0.090304 
    ## Run 421 stress 0.09030394 
    ## Run 422 stress 0.0996998 
    ## Run 423 stress 0.09539208 
    ## Run 424 stress 0.09416133 
    ## Run 425 stress 0.09535513 
    ## Run 426 stress 0.08440253 
    ## ... Procrustes: rmse 8.256353e-05  max resid 0.000178969 
    ## ... Similar to previous best
    ## Run 427 stress 0.09145334 
    ## Run 428 stress 0.09539196 
    ## Run 429 stress 0.09030399 
    ## Run 430 stress 0.09447302 
    ## Run 431 stress 0.09268366 
    ## Run 432 stress 0.08503573 
    ## Run 433 stress 0.08973891 
    ## Run 434 stress 0.09465915 
    ## Run 435 stress 0.09337299 
    ## Run 436 stress 0.09760839 
    ## Run 437 stress 0.09539198 
    ## Run 438 stress 0.09268354 
    ## Run 439 stress 0.09969972 
    ## Run 440 stress 0.09969985 
    ## Run 441 stress 0.094034 
    ## Run 442 stress 0.08503586 
    ## Run 443 stress 0.1052849 
    ## Run 444 stress 0.09407118 
    ## Run 445 stress 0.09407979 
    ## Run 446 stress 0.09374239 
    ## Run 447 stress 0.08973863 
    ## Run 448 stress 0.09407983 
    ## Run 449 stress 0.09374314 
    ## Run 450 stress 0.09286084 
    ## Run 451 stress 0.09969966 
    ## Run 452 stress 0.08973892 
    ## Run 453 stress 0.08773484 
    ## Run 454 stress 0.09168943 
    ## Run 455 stress 0.09030418 
    ## Run 456 stress 0.09416146 
    ## Run 457 stress 0.0953547 
    ## Run 458 stress 0.09465904 
    ## Run 459 stress 0.09465909 
    ## Run 460 stress 0.0903042 
    ## Run 461 stress 0.09465904 
    ## Run 462 stress 0.09969975 
    ## Run 463 stress 0.08503654 
    ## Run 464 stress 0.09159085 
    ## Run 465 stress 0.09030402 
    ## Run 466 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001768122  max resid 0.0003328502 
    ## ... Similar to previous best
    ## Run 467 stress 0.1038044 
    ## Run 468 stress 0.08503465 
    ## Run 469 stress 0.1052849 
    ## Run 470 stress 0.09286083 
    ## Run 471 stress 0.0928609 
    ## Run 472 stress 0.09760816 
    ## Run 473 stress 0.09030408 
    ## Run 474 stress 0.09465911 
    ## Run 475 stress 0.08440254 
    ## ... Procrustes: rmse 8.727573e-05  max resid 0.0001843867 
    ## ... Similar to previous best
    ## Run 476 stress 0.0850351 
    ## Run 477 stress 0.0914533 
    ## Run 478 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001764237  max resid 0.0003324958 
    ## ... Similar to previous best
    ## Run 479 stress 0.08503658 
    ## Run 480 stress 0.09168937 
    ## Run 481 stress 0.09145322 
    ## Run 482 stress 0.08973864 
    ## Run 483 stress 0.09434653 
    ## Run 484 stress 0.09969966 
    ## Run 485 stress 0.08503463 
    ## Run 486 stress 0.0850346 
    ## Run 487 stress 0.09286083 
    ## Run 488 stress 0.09308946 
    ## Run 489 stress 0.09416135 
    ## Run 490 stress 0.1038042 
    ## Run 491 stress 0.09407991 
    ## Run 492 stress 0.09168937 
    ## Run 493 stress 0.09168944 
    ## Run 494 stress 0.09145326 
    ## Run 495 stress 0.08773472 
    ## Run 496 stress 0.1038045 
    ## Run 497 stress 0.1072864 
    ## Run 498 stress 0.1038043 
    ## Run 499 stress 0.09447311 
    ## Run 500 stress 0.08503514 
    ## *** Best solution repeated 19 times

``` r
# Ocean sites and mixed lakes
SD_beta_geo_OM_NMDS <- metaMDS(SD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07365783 
    ## Run 1 stress 0.07732931 
    ## Run 2 stress 0.07365785 
    ## ... Procrustes: rmse 6.204959e-05  max resid 0.0001461593 
    ## ... Similar to previous best
    ## Run 3 stress 0.08000468 
    ## Run 4 stress 0.08288056 
    ## Run 5 stress 0.07365783 
    ## ... Procrustes: rmse 2.020975e-05  max resid 4.52553e-05 
    ## ... Similar to previous best
    ## Run 6 stress 0.0828804 
    ## Run 7 stress 0.07629238 
    ## Run 8 stress 0.07732927 
    ## Run 9 stress 0.07378226 
    ## ... Procrustes: rmse 0.017445  max resid 0.05104887 
    ## Run 10 stress 0.0800048 
    ## Run 11 stress 0.3153601 
    ## Run 12 stress 0.07629234 
    ## Run 13 stress 0.07629247 
    ## Run 14 stress 0.07732944 
    ## Run 15 stress 0.07365783 
    ## ... Procrustes: rmse 1.701281e-05  max resid 3.519201e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.08000467 
    ## Run 17 stress 0.07365785 
    ## ... Procrustes: rmse 7.4254e-05  max resid 0.0001736535 
    ## ... Similar to previous best
    ## Run 18 stress 0.07629235 
    ## Run 19 stress 0.08233764 
    ## Run 20 stress 0.08233775 
    ## Run 21 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001049018  max resid 0.0002477259 
    ## ... Similar to previous best
    ## Run 22 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744782  max resid 0.05106717 
    ## Run 23 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744375  max resid 0.05103684 
    ## Run 24 stress 0.07629235 
    ## Run 25 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744157  max resid 0.05102082 
    ## Run 26 stress 0.07629234 
    ## Run 27 stress 0.07365788 
    ## ... Procrustes: rmse 7.615404e-05  max resid 0.0001813403 
    ## ... Similar to previous best
    ## Run 28 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744215  max resid 0.05103549 
    ## Run 29 stress 0.07365783 
    ## ... Procrustes: rmse 1.258845e-05  max resid 2.533803e-05 
    ## ... Similar to previous best
    ## Run 30 stress 0.08233778 
    ## Run 31 stress 0.07378232 
    ## ... Procrustes: rmse 0.01744561  max resid 0.05101733 
    ## Run 32 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744209  max resid 0.05101708 
    ## Run 33 stress 0.08233776 
    ## Run 34 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174341  max resid 0.05099946 
    ## Run 35 stress 0.08000468 
    ## Run 36 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001184222  max resid 0.000280812 
    ## ... Similar to previous best
    ## Run 37 stress 0.07365784 
    ## ... Procrustes: rmse 3.688509e-05  max resid 8.703697e-05 
    ## ... Similar to previous best
    ## Run 38 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001059511  max resid 0.0002512776 
    ## ... Similar to previous best
    ## Run 39 stress 0.08233775 
    ## Run 40 stress 0.07629235 
    ## Run 41 stress 0.07365786 
    ## ... Procrustes: rmse 7.85903e-05  max resid 0.0001878699 
    ## ... Similar to previous best
    ## Run 42 stress 0.07629238 
    ## Run 43 stress 0.07629248 
    ## Run 44 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001107494  max resid 0.0002588266 
    ## ... Similar to previous best
    ## Run 45 stress 0.07629238 
    ## Run 46 stress 0.07629243 
    ## Run 47 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 9.467625e-06  max resid 2.010992e-05 
    ## ... Similar to previous best
    ## Run 48 stress 0.08000469 
    ## Run 49 stress 0.07629235 
    ## Run 50 stress 0.08000475 
    ## Run 51 stress 0.07629233 
    ## Run 52 stress 0.07732935 
    ## Run 53 stress 0.07378232 
    ## ... Procrustes: rmse 0.01744335  max resid 0.05099833 
    ## Run 54 stress 0.07365788 
    ## ... Procrustes: rmse 0.000121034  max resid 0.0002864197 
    ## ... Similar to previous best
    ## Run 55 stress 0.07732944 
    ## Run 56 stress 0.07365788 
    ## ... Procrustes: rmse 5.582778e-05  max resid 0.000109162 
    ## ... Similar to previous best
    ## Run 57 stress 0.07732946 
    ## Run 58 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744692  max resid 0.05104112 
    ## Run 59 stress 0.07629238 
    ## Run 60 stress 0.07365785 
    ## ... Procrustes: rmse 7.6058e-05  max resid 0.0001756826 
    ## ... Similar to previous best
    ## Run 61 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744713  max resid 0.05104551 
    ## Run 62 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745559  max resid 0.05105864 
    ## Run 63 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744716  max resid 0.05104793 
    ## Run 64 stress 0.08288039 
    ## Run 65 stress 0.08000467 
    ## Run 66 stress 0.08000474 
    ## Run 67 stress 0.07732939 
    ## Run 68 stress 0.07365785 
    ## ... Procrustes: rmse 6.483033e-05  max resid 0.0001544832 
    ## ... Similar to previous best
    ## Run 69 stress 0.07365786 
    ## ... Procrustes: rmse 8.79984e-05  max resid 0.0002078128 
    ## ... Similar to previous best
    ## Run 70 stress 0.07365784 
    ## ... Procrustes: rmse 3.290448e-05  max resid 7.44211e-05 
    ## ... Similar to previous best
    ## Run 71 stress 0.07365786 
    ## ... Procrustes: rmse 9.028417e-05  max resid 0.0002132733 
    ## ... Similar to previous best
    ## Run 72 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744706  max resid 0.05102473 
    ## Run 73 stress 0.07378227 
    ## ... Procrustes: rmse 0.01746843  max resid 0.05114131 
    ## Run 74 stress 0.07365783 
    ## ... Procrustes: rmse 2.770216e-05  max resid 6.518167e-05 
    ## ... Similar to previous best
    ## Run 75 stress 0.08233757 
    ## Run 76 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001384572  max resid 0.0003254488 
    ## ... Similar to previous best
    ## Run 77 stress 0.08000478 
    ## Run 78 stress 0.07732932 
    ## Run 79 stress 0.07365785 
    ## ... Procrustes: rmse 2.359666e-05  max resid 4.032013e-05 
    ## ... Similar to previous best
    ## Run 80 stress 0.07365783 
    ## ... Procrustes: rmse 9.978038e-06  max resid 2.69417e-05 
    ## ... Similar to previous best
    ## Run 81 stress 0.08000474 
    ## Run 82 stress 0.07365785 
    ## ... Procrustes: rmse 2.619936e-05  max resid 5.631387e-05 
    ## ... Similar to previous best
    ## Run 83 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743487  max resid 0.05100352 
    ## Run 84 stress 0.07378233 
    ## ... Procrustes: rmse 0.01745021  max resid 0.05101793 
    ## Run 85 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744813  max resid 0.05105212 
    ## Run 86 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744837  max resid 0.05105229 
    ## Run 87 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744662  max resid 0.05103422 
    ## Run 88 stress 0.07629236 
    ## Run 89 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745349  max resid 0.05107346 
    ## Run 90 stress 0.07629234 
    ## Run 91 stress 0.0773293 
    ## Run 92 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174444  max resid 0.05102412 
    ## Run 93 stress 0.07365783 
    ## ... Procrustes: rmse 1.096767e-05  max resid 2.564858e-05 
    ## ... Similar to previous best
    ## Run 94 stress 0.07629244 
    ## Run 95 stress 0.07365796 
    ## ... Procrustes: rmse 0.0001887716  max resid 0.0004473142 
    ## ... Similar to previous best
    ## Run 96 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001068524  max resid 0.0002451657 
    ## ... Similar to previous best
    ## Run 97 stress 0.08000471 
    ## Run 98 stress 0.08233753 
    ## Run 99 stress 0.07732933 
    ## Run 100 stress 0.07629234 
    ## Run 101 stress 0.08000469 
    ## Run 102 stress 0.07732927 
    ## Run 103 stress 0.07629233 
    ## Run 104 stress 0.07629234 
    ## Run 105 stress 0.08000474 
    ## Run 106 stress 0.08000473 
    ## Run 107 stress 0.07629241 
    ## Run 108 stress 0.07365783 
    ## ... Procrustes: rmse 1.410935e-05  max resid 2.69497e-05 
    ## ... Similar to previous best
    ## Run 109 stress 0.07365783 
    ## ... Procrustes: rmse 1.197092e-05  max resid 2.722393e-05 
    ## ... Similar to previous best
    ## Run 110 stress 0.07365784 
    ## ... Procrustes: rmse 5.192491e-05  max resid 0.0001194942 
    ## ... Similar to previous best
    ## Run 111 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745074  max resid 0.05106507 
    ## Run 112 stress 0.07365794 
    ## ... Procrustes: rmse 0.000174836  max resid 0.0004144106 
    ## ... Similar to previous best
    ## Run 113 stress 0.08000469 
    ## Run 114 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744334  max resid 0.05100995 
    ## Run 115 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744745  max resid 0.05105031 
    ## Run 116 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744675  max resid 0.05103457 
    ## Run 117 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001549584  max resid 0.0003620861 
    ## ... Similar to previous best
    ## Run 118 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744883  max resid 0.0510642 
    ## Run 119 stress 0.07629234 
    ## Run 120 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744406  max resid 0.05106098 
    ## Run 121 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001271819  max resid 0.0002978627 
    ## ... Similar to previous best
    ## Run 122 stress 0.07629237 
    ## Run 123 stress 0.07365785 
    ## ... Procrustes: rmse 5.947733e-05  max resid 0.000138016 
    ## ... Similar to previous best
    ## Run 124 stress 0.07629234 
    ## Run 125 stress 0.07629233 
    ## Run 126 stress 0.07365785 
    ## ... Procrustes: rmse 7.085782e-05  max resid 0.0001629572 
    ## ... Similar to previous best
    ## Run 127 stress 0.07629245 
    ## Run 128 stress 0.08000472 
    ## Run 129 stress 0.07629243 
    ## Run 130 stress 0.07365783 
    ## ... Procrustes: rmse 2.131235e-05  max resid 5.014268e-05 
    ## ... Similar to previous best
    ## Run 131 stress 0.2574761 
    ## Run 132 stress 0.07365783 
    ## ... Procrustes: rmse 2.294535e-05  max resid 5.360244e-05 
    ## ... Similar to previous best
    ## Run 133 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743734  max resid 0.0510059 
    ## Run 134 stress 0.07629233 
    ## Run 135 stress 0.0800047 
    ## Run 136 stress 0.07365785 
    ## ... Procrustes: rmse 5.594945e-05  max resid 0.000125182 
    ## ... Similar to previous best
    ## Run 137 stress 0.0828804 
    ## Run 138 stress 0.07365784 
    ## ... Procrustes: rmse 3.422418e-05  max resid 7.339991e-05 
    ## ... Similar to previous best
    ## Run 139 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744704  max resid 0.05102088 
    ## Run 140 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745209  max resid 0.05107878 
    ## Run 141 stress 0.07629246 
    ## Run 142 stress 0.07365783 
    ## ... Procrustes: rmse 1.106029e-05  max resid 2.502924e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.07365784 
    ## ... Procrustes: rmse 4.744513e-05  max resid 0.0001099525 
    ## ... Similar to previous best
    ## Run 144 stress 0.08233756 
    ## Run 145 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001115996  max resid 0.0002629971 
    ## ... Similar to previous best
    ## Run 146 stress 0.07732942 
    ## Run 147 stress 0.07365783 
    ## ... Procrustes: rmse 2.019677e-05  max resid 4.884062e-05 
    ## ... Similar to previous best
    ## Run 148 stress 0.08000471 
    ## Run 149 stress 0.08233778 
    ## Run 150 stress 0.07629241 
    ## Run 151 stress 0.07732937 
    ## Run 152 stress 0.07378233 
    ## ... Procrustes: rmse 0.0174393  max resid 0.05105401 
    ## Run 153 stress 0.07365784 
    ## ... Procrustes: rmse 1.345893e-05  max resid 2.523695e-05 
    ## ... Similar to previous best
    ## Run 154 stress 0.07365787 
    ## ... Procrustes: rmse 9.556758e-05  max resid 0.0002236221 
    ## ... Similar to previous best
    ## Run 155 stress 0.07365783 
    ## ... Procrustes: rmse 2.881023e-05  max resid 6.789113e-05 
    ## ... Similar to previous best
    ## Run 156 stress 0.07629236 
    ## Run 157 stress 0.08288044 
    ## Run 158 stress 0.07629244 
    ## Run 159 stress 0.07365784 
    ## ... Procrustes: rmse 4.952174e-05  max resid 0.0001162173 
    ## ... Similar to previous best
    ## Run 160 stress 0.07629235 
    ## Run 161 stress 0.08000468 
    ## Run 162 stress 0.288024 
    ## Run 163 stress 0.08233753 
    ## Run 164 stress 0.08000468 
    ## Run 165 stress 0.08000467 
    ## Run 166 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745669  max resid 0.0511036 
    ## Run 167 stress 0.07365784 
    ## ... Procrustes: rmse 3.78395e-05  max resid 8.8922e-05 
    ## ... Similar to previous best
    ## Run 168 stress 0.08000468 
    ## Run 169 stress 0.08233759 
    ## Run 170 stress 0.2437517 
    ## Run 171 stress 0.07365783 
    ## ... Procrustes: rmse 3.394385e-06  max resid 7.015768e-06 
    ## ... Similar to previous best
    ## Run 172 stress 0.07365785 
    ## ... Procrustes: rmse 3.182916e-05  max resid 6.630058e-05 
    ## ... Similar to previous best
    ## Run 173 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745031  max resid 0.05106807 
    ## Run 174 stress 0.08233754 
    ## Run 175 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744981  max resid 0.05105122 
    ## Run 176 stress 0.07378232 
    ## ... Procrustes: rmse 0.01746162  max resid 0.05113878 
    ## Run 177 stress 0.07629233 
    ## Run 178 stress 0.0762924 
    ## Run 179 stress 0.07365784 
    ## ... Procrustes: rmse 4.377861e-05  max resid 0.0001025381 
    ## ... Similar to previous best
    ## Run 180 stress 0.07732934 
    ## Run 181 stress 0.08000478 
    ## Run 182 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744997  max resid 0.05106693 
    ## Run 183 stress 0.08000479 
    ## Run 184 stress 0.08000473 
    ## Run 185 stress 0.07365784 
    ## ... Procrustes: rmse 4.542858e-05  max resid 0.0001083059 
    ## ... Similar to previous best
    ## Run 186 stress 0.07629234 
    ## Run 187 stress 0.07629235 
    ## Run 188 stress 0.07365788 
    ## ... Procrustes: rmse 7.938653e-05  max resid 0.0001887484 
    ## ... Similar to previous best
    ## Run 189 stress 0.07365797 
    ## ... Procrustes: rmse 0.0002010513  max resid 0.0004736406 
    ## ... Similar to previous best
    ## Run 190 stress 0.08233758 
    ## Run 191 stress 0.07365785 
    ## ... Procrustes: rmse 7.124832e-05  max resid 0.0001692758 
    ## ... Similar to previous best
    ## Run 192 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744691  max resid 0.05106382 
    ## Run 193 stress 0.08000471 
    ## Run 194 stress 0.0800047 
    ## Run 195 stress 0.07732925 
    ## Run 196 stress 0.07732932 
    ## Run 197 stress 0.3265265 
    ## Run 198 stress 0.07732944 
    ## Run 199 stress 0.07365784 
    ## ... Procrustes: rmse 4.779915e-05  max resid 0.0001119865 
    ## ... Similar to previous best
    ## Run 200 stress 0.07365783 
    ## ... Procrustes: rmse 2.459672e-05  max resid 5.800887e-05 
    ## ... Similar to previous best
    ## Run 201 stress 0.0762924 
    ## Run 202 stress 0.07365784 
    ## ... Procrustes: rmse 4.840918e-05  max resid 0.0001135425 
    ## ... Similar to previous best
    ## Run 203 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744817  max resid 0.05105274 
    ## Run 204 stress 0.07365783 
    ## ... Procrustes: rmse 3.539847e-05  max resid 8.325951e-05 
    ## ... Similar to previous best
    ## Run 205 stress 0.08000475 
    ## Run 206 stress 0.07365786 
    ## ... Procrustes: rmse 9.095299e-05  max resid 0.0002153809 
    ## ... Similar to previous best
    ## Run 207 stress 0.08000477 
    ## Run 208 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745254  max resid 0.05108083 
    ## Run 209 stress 0.08000467 
    ## Run 210 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744825  max resid 0.05106661 
    ## Run 211 stress 0.07732935 
    ## Run 212 stress 0.0762924 
    ## Run 213 stress 0.08233774 
    ## Run 214 stress 0.07629238 
    ## Run 215 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744878  max resid 0.05106252 
    ## Run 216 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174522  max resid 0.05107708 
    ## Run 217 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745371  max resid 0.05109089 
    ## Run 218 stress 0.07365783 
    ## ... Procrustes: rmse 2.928767e-05  max resid 6.832477e-05 
    ## ... Similar to previous best
    ## Run 219 stress 0.07629235 
    ## Run 220 stress 0.07732928 
    ## Run 221 stress 0.08233755 
    ## Run 222 stress 0.08288038 
    ## Run 223 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001685967  max resid 0.0003957156 
    ## ... Similar to previous best
    ## Run 224 stress 0.07629234 
    ## Run 225 stress 0.0773293 
    ## Run 226 stress 0.08000468 
    ## Run 227 stress 0.07365787 
    ## ... Procrustes: rmse 9.240305e-05  max resid 0.0002204957 
    ## ... Similar to previous best
    ## Run 228 stress 0.07365785 
    ## ... Procrustes: rmse 6.64909e-05  max resid 0.0001543886 
    ## ... Similar to previous best
    ## Run 229 stress 0.0762924 
    ## Run 230 stress 0.07629236 
    ## Run 231 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745074  max resid 0.05105889 
    ## Run 232 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744532  max resid 0.05102907 
    ## Run 233 stress 0.07629234 
    ## Run 234 stress 0.07629233 
    ## Run 235 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174475  max resid 0.05104301 
    ## Run 236 stress 0.07365785 
    ## ... Procrustes: rmse 5.780326e-05  max resid 0.0001375348 
    ## ... Similar to previous best
    ## Run 237 stress 0.07629234 
    ## Run 238 stress 0.07365785 
    ## ... Procrustes: rmse 6.580156e-05  max resid 0.0001526546 
    ## ... Similar to previous best
    ## Run 239 stress 0.07378229 
    ## ... Procrustes: rmse 0.0174456  max resid 0.05101224 
    ## Run 240 stress 0.08233764 
    ## Run 241 stress 0.07365783 
    ## ... Procrustes: rmse 7.566416e-06  max resid 1.561115e-05 
    ## ... Similar to previous best
    ## Run 242 stress 0.07629237 
    ## Run 243 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745381  max resid 0.05107464 
    ## Run 244 stress 0.07732932 
    ## Run 245 stress 0.08233759 
    ## Run 246 stress 0.07365783 
    ## ... Procrustes: rmse 1.202914e-05  max resid 2.344716e-05 
    ## ... Similar to previous best
    ## Run 247 stress 0.07732942 
    ## Run 248 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745406  max resid 0.05108811 
    ## Run 249 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745228  max resid 0.05107834 
    ## Run 250 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744546  max resid 0.05104167 
    ## Run 251 stress 0.07365783 
    ## ... Procrustes: rmse 1.598445e-05  max resid 3.871395e-05 
    ## ... Similar to previous best
    ## Run 252 stress 0.07629241 
    ## Run 253 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744372  max resid 0.05102298 
    ## Run 254 stress 0.07629236 
    ## Run 255 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745143  max resid 0.05107506 
    ## Run 256 stress 0.08000479 
    ## Run 257 stress 0.07629233 
    ## Run 258 stress 0.07629237 
    ## Run 259 stress 0.07629238 
    ## Run 260 stress 0.07365788 
    ## ... Procrustes: rmse 9.838774e-05  max resid 0.0002348392 
    ## ... Similar to previous best
    ## Run 261 stress 0.08000467 
    ## Run 262 stress 0.07365784 
    ## ... Procrustes: rmse 4.804581e-05  max resid 0.0001130177 
    ## ... Similar to previous best
    ## Run 263 stress 0.08000468 
    ## Run 264 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743721  max resid 0.05100496 
    ## Run 265 stress 0.08233762 
    ## Run 266 stress 0.07365783 
    ## ... Procrustes: rmse 2.74935e-05  max resid 6.508716e-05 
    ## ... Similar to previous best
    ## Run 267 stress 0.0828805 
    ## Run 268 stress 0.07365783 
    ## ... Procrustes: rmse 7.953986e-06  max resid 1.567625e-05 
    ## ... Similar to previous best
    ## Run 269 stress 0.08000469 
    ## Run 270 stress 0.07378232 
    ## ... Procrustes: rmse 0.01747253  max resid 0.0511771 
    ## Run 271 stress 0.07629233 
    ## Run 272 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745706  max resid 0.05109902 
    ## Run 273 stress 0.07629237 
    ## Run 274 stress 0.0762924 
    ## Run 275 stress 0.07629234 
    ## Run 276 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174439  max resid 0.05104637 
    ## Run 277 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174526  max resid 0.05106698 
    ## Run 278 stress 0.08000467 
    ## Run 279 stress 0.3512886 
    ## Run 280 stress 0.07365784 
    ## ... Procrustes: rmse 6.235804e-05  max resid 0.0001468386 
    ## ... Similar to previous best
    ## Run 281 stress 0.0737823 
    ## ... Procrustes: rmse 0.0174783  max resid 0.05119902 
    ## Run 282 stress 0.07629235 
    ## Run 283 stress 0.07378228 
    ## ... Procrustes: rmse 0.017444  max resid 0.05102376 
    ## Run 284 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001451526  max resid 0.0003424651 
    ## ... Similar to previous best
    ## Run 285 stress 0.08000475 
    ## Run 286 stress 0.07629249 
    ## Run 287 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744679  max resid 0.0510423 
    ## Run 288 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001061478  max resid 0.0002469217 
    ## ... Similar to previous best
    ## Run 289 stress 0.07629243 
    ## Run 290 stress 0.07732933 
    ## Run 291 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744774  max resid 0.05104923 
    ## Run 292 stress 0.07629234 
    ## Run 293 stress 0.07629243 
    ## Run 294 stress 0.2790898 
    ## Run 295 stress 0.07629243 
    ## Run 296 stress 0.07732932 
    ## Run 297 stress 0.08288051 
    ## Run 298 stress 0.07365783 
    ## ... Procrustes: rmse 2.105585e-05  max resid 4.978838e-05 
    ## ... Similar to previous best
    ## Run 299 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743156  max resid 0.0509757 
    ## Run 300 stress 0.08000471 
    ## Run 301 stress 0.07732934 
    ## Run 302 stress 0.08288075 
    ## Run 303 stress 0.08000469 
    ## Run 304 stress 0.07365785 
    ## ... Procrustes: rmse 5.689567e-05  max resid 0.0001239429 
    ## ... Similar to previous best
    ## Run 305 stress 0.07365783 
    ## ... Procrustes: rmse 7.446726e-06  max resid 1.57879e-05 
    ## ... Similar to previous best
    ## Run 306 stress 0.07629238 
    ## Run 307 stress 0.07629238 
    ## Run 308 stress 0.08288058 
    ## Run 309 stress 0.07365784 
    ## ... Procrustes: rmse 5.066334e-05  max resid 0.0001160408 
    ## ... Similar to previous best
    ## Run 310 stress 0.07732933 
    ## Run 311 stress 0.07629236 
    ## Run 312 stress 0.07365784 
    ## ... Procrustes: rmse 4.552431e-05  max resid 0.0001065852 
    ## ... Similar to previous best
    ## Run 313 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001119303  max resid 0.0002678229 
    ## ... Similar to previous best
    ## Run 314 stress 0.07732927 
    ## Run 315 stress 0.07365783 
    ## ... Procrustes: rmse 6.363805e-06  max resid 1.421018e-05 
    ## ... Similar to previous best
    ## Run 316 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744793  max resid 0.05106236 
    ## Run 317 stress 0.07365783 
    ## ... Procrustes: rmse 2.54075e-05  max resid 6.009238e-05 
    ## ... Similar to previous best
    ## Run 318 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001296174  max resid 0.0003033491 
    ## ... Similar to previous best
    ## Run 319 stress 0.0800047 
    ## Run 320 stress 0.08000468 
    ## Run 321 stress 0.08288041 
    ## Run 322 stress 0.07629249 
    ## Run 323 stress 0.07365785 
    ## ... Procrustes: rmse 7.108469e-05  max resid 0.0001695516 
    ## ... Similar to previous best
    ## Run 324 stress 0.08000474 
    ## Run 325 stress 0.07732925 
    ## Run 326 stress 0.07365784 
    ## ... Procrustes: rmse 4.905073e-05  max resid 0.0001160793 
    ## ... Similar to previous best
    ## Run 327 stress 0.08000467 
    ## Run 328 stress 0.08000468 
    ## Run 329 stress 0.07732934 
    ## Run 330 stress 0.08000468 
    ## Run 331 stress 0.0800047 
    ## Run 332 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001324457  max resid 0.0003123224 
    ## ... Similar to previous best
    ## Run 333 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001287422  max resid 0.000308474 
    ## ... Similar to previous best
    ## Run 334 stress 0.08288038 
    ## Run 335 stress 0.07365783 
    ## ... Procrustes: rmse 1.253074e-05  max resid 2.934939e-05 
    ## ... Similar to previous best
    ## Run 336 stress 0.07732926 
    ## Run 337 stress 0.07629248 
    ## Run 338 stress 0.0800047 
    ## Run 339 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174492  max resid 0.05104366 
    ## Run 340 stress 0.07629235 
    ## Run 341 stress 0.0800047 
    ## Run 342 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745179  max resid 0.05107717 
    ## Run 343 stress 0.07365786 
    ## ... Procrustes: rmse 7.932014e-05  max resid 0.00018807 
    ## ... Similar to previous best
    ## Run 344 stress 0.07365785 
    ## ... Procrustes: rmse 2.282769e-05  max resid 3.290023e-05 
    ## ... Similar to previous best
    ## Run 345 stress 0.0800047 
    ## Run 346 stress 0.07365784 
    ## ... Procrustes: rmse 2.609241e-05  max resid 4.946043e-05 
    ## ... Similar to previous best
    ## Run 347 stress 0.07365783 
    ## ... Procrustes: rmse 1.436449e-05  max resid 3.025264e-05 
    ## ... Similar to previous best
    ## Run 348 stress 0.0737823 
    ## ... Procrustes: rmse 0.01745573  max resid 0.05110161 
    ## Run 349 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174595  max resid 0.05108541 
    ## Run 350 stress 0.07365784 
    ## ... Procrustes: rmse 4.858735e-05  max resid 0.0001131938 
    ## ... Similar to previous best
    ## Run 351 stress 0.08233761 
    ## Run 352 stress 0.07629234 
    ## Run 353 stress 0.07629247 
    ## Run 354 stress 0.07378228 
    ## ... Procrustes: rmse 0.01746634  max resid 0.05111182 
    ## Run 355 stress 0.0762924 
    ## Run 356 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744651  max resid 0.05104815 
    ## Run 357 stress 0.07365784 
    ## ... Procrustes: rmse 2.738758e-05  max resid 6.09753e-05 
    ## ... Similar to previous best
    ## Run 358 stress 0.0737823 
    ## ... Procrustes: rmse 0.01745821  max resid 0.05111147 
    ## Run 359 stress 0.0737823 
    ## ... Procrustes: rmse 0.01746197  max resid 0.05108163 
    ## Run 360 stress 0.08000469 
    ## Run 361 stress 0.07629233 
    ## Run 362 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744749  max resid 0.05104903 
    ## Run 363 stress 0.08233759 
    ## Run 364 stress 0.07629244 
    ## Run 365 stress 0.07365787 
    ## ... Procrustes: rmse 9.061307e-05  max resid 0.0002053213 
    ## ... Similar to previous best
    ## Run 366 stress 0.08233755 
    ## Run 367 stress 0.07378233 
    ## ... Procrustes: rmse 0.01743976  max resid 0.05105408 
    ## Run 368 stress 0.07365783 
    ## ... Procrustes: rmse 2.000034e-05  max resid 4.621273e-05 
    ## ... Similar to previous best
    ## Run 369 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174352  max resid 0.05100512 
    ## Run 370 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744876  max resid 0.05105558 
    ## Run 371 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745144  max resid 0.05106289 
    ## Run 372 stress 0.07365783 
    ## ... Procrustes: rmse 2.582115e-05  max resid 5.843088e-05 
    ## ... Similar to previous best
    ## Run 373 stress 0.08000468 
    ## Run 374 stress 0.08000468 
    ## Run 375 stress 0.07629237 
    ## Run 376 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001508289  max resid 0.0003468877 
    ## ... Similar to previous best
    ## Run 377 stress 0.07365785 
    ## ... Procrustes: rmse 5.568409e-05  max resid 0.0001246383 
    ## ... Similar to previous best
    ## Run 378 stress 0.08000468 
    ## Run 379 stress 0.08233755 
    ## Run 380 stress 0.07629243 
    ## Run 381 stress 0.08233754 
    ## Run 382 stress 0.07732929 
    ## Run 383 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744637  max resid 0.05103693 
    ## Run 384 stress 0.07629242 
    ## Run 385 stress 0.08288049 
    ## Run 386 stress 0.07629245 
    ## Run 387 stress 0.08000468 
    ## Run 388 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001202325  max resid 0.0002876028 
    ## ... Similar to previous best
    ## Run 389 stress 0.08000478 
    ## Run 390 stress 0.07365786 
    ## ... Procrustes: rmse 9.135686e-05  max resid 0.0002158019 
    ## ... Similar to previous best
    ## Run 391 stress 0.07378228 
    ## ... Procrustes: rmse 0.01747394  max resid 0.05117133 
    ## Run 392 stress 0.07365785 
    ## ... Procrustes: rmse 7.454279e-05  max resid 0.0001751091 
    ## ... Similar to previous best
    ## Run 393 stress 0.3582285 
    ## Run 394 stress 0.2756221 
    ## Run 395 stress 0.07365786 
    ## ... Procrustes: rmse 2.395404e-05  max resid 4.951112e-05 
    ## ... Similar to previous best
    ## Run 396 stress 0.08000467 
    ## Run 397 stress 0.07365783 
    ## ... Procrustes: rmse 1.666868e-05  max resid 3.919737e-05 
    ## ... Similar to previous best
    ## Run 398 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744304  max resid 0.05100464 
    ## Run 399 stress 0.07629235 
    ## Run 400 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001308711  max resid 0.0002997636 
    ## ... Similar to previous best
    ## Run 401 stress 0.07629233 
    ## Run 402 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744314  max resid 0.05105906 
    ## Run 403 stress 0.07365784 
    ## ... Procrustes: rmse 4.599382e-05  max resid 0.0001069979 
    ## ... Similar to previous best
    ## Run 404 stress 0.07365785 
    ## ... Procrustes: rmse 7.653447e-05  max resid 0.000172737 
    ## ... Similar to previous best
    ## Run 405 stress 0.07365784 
    ## ... Procrustes: rmse 3.875556e-05  max resid 8.437171e-05 
    ## ... Similar to previous best
    ## Run 406 stress 0.08000471 
    ## Run 407 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744971  max resid 0.05106522 
    ## Run 408 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744806  max resid 0.05104766 
    ## Run 409 stress 0.07365784 
    ## ... Procrustes: rmse 2.632134e-05  max resid 5.302519e-05 
    ## ... Similar to previous best
    ## Run 410 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001576404  max resid 0.0003768656 
    ## ... Similar to previous best
    ## Run 411 stress 0.07629235 
    ## Run 412 stress 0.07629234 
    ## Run 413 stress 0.08000469 
    ## Run 414 stress 0.07629234 
    ## Run 415 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744806  max resid 0.05105967 
    ## Run 416 stress 0.07629239 
    ## Run 417 stress 0.0762924 
    ## Run 418 stress 0.07365784 
    ## ... Procrustes: rmse 5.578233e-05  max resid 0.0001323787 
    ## ... Similar to previous best
    ## Run 419 stress 0.08288045 
    ## Run 420 stress 0.07365784 
    ## ... Procrustes: rmse 4.81099e-05  max resid 0.0001137392 
    ## ... Similar to previous best
    ## Run 421 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745076  max resid 0.05105195 
    ## Run 422 stress 0.07365783 
    ## ... Procrustes: rmse 8.385373e-06  max resid 1.770222e-05 
    ## ... Similar to previous best
    ## Run 423 stress 0.3495556 
    ## Run 424 stress 0.08233758 
    ## Run 425 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744773  max resid 0.05102457 
    ## Run 426 stress 0.07365786 
    ## ... Procrustes: rmse 8.492913e-05  max resid 0.0002027279 
    ## ... Similar to previous best
    ## Run 427 stress 0.07365785 
    ## ... Procrustes: rmse 7.337378e-05  max resid 0.0001722688 
    ## ... Similar to previous best
    ## Run 428 stress 0.07365783 
    ## ... Procrustes: rmse 2.238124e-05  max resid 4.948543e-05 
    ## ... Similar to previous best
    ## Run 429 stress 0.07629247 
    ## Run 430 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744839  max resid 0.05104194 
    ## Run 431 stress 0.08233757 
    ## Run 432 stress 0.07732936 
    ## Run 433 stress 0.08000468 
    ## Run 434 stress 0.07629241 
    ## Run 435 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001546291  max resid 0.0003649815 
    ## ... Similar to previous best
    ## Run 436 stress 0.2870318 
    ## Run 437 stress 0.08288039 
    ## Run 438 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744965  max resid 0.05105625 
    ## Run 439 stress 0.08000476 
    ## Run 440 stress 0.0773293 
    ## Run 441 stress 0.07629233 
    ## Run 442 stress 0.07629238 
    ## Run 443 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001391669  max resid 0.0003317019 
    ## ... Similar to previous best
    ## Run 444 stress 0.07629237 
    ## Run 445 stress 0.07732941 
    ## Run 446 stress 0.07365785 
    ## ... Procrustes: rmse 6.923698e-05  max resid 0.0001635321 
    ## ... Similar to previous best
    ## Run 447 stress 0.0800047 
    ## Run 448 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744678  max resid 0.05106313 
    ## Run 449 stress 0.07365784 
    ## ... Procrustes: rmse 1.665422e-05  max resid 3.155378e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744475  max resid 0.05103752 
    ## Run 451 stress 0.08000469 
    ## Run 452 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745164  max resid 0.05107852 
    ## Run 453 stress 0.08000471 
    ## Run 454 stress 0.07365783 
    ## ... Procrustes: rmse 8.040265e-06  max resid 1.716642e-05 
    ## ... Similar to previous best
    ## Run 455 stress 0.08000469 
    ## Run 456 stress 0.08000467 
    ## Run 457 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745387  max resid 0.05109047 
    ## Run 458 stress 0.07365785 
    ## ... Procrustes: rmse 6.502698e-05  max resid 0.0001510705 
    ## ... Similar to previous best
    ## Run 459 stress 0.08000471 
    ## Run 460 stress 0.08000471 
    ## Run 461 stress 0.0800048 
    ## Run 462 stress 0.07629235 
    ## Run 463 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744711  max resid 0.05104425 
    ## Run 464 stress 0.07629243 
    ## Run 465 stress 0.08233766 
    ## Run 466 stress 0.07629239 
    ## Run 467 stress 0.07629234 
    ## Run 468 stress 0.08000473 
    ## Run 469 stress 0.07365788 
    ## ... Procrustes: rmse 9.860586e-05  max resid 0.0002291107 
    ## ... Similar to previous best
    ## Run 470 stress 0.07629237 
    ## Run 471 stress 0.07629235 
    ## Run 472 stress 0.08233775 
    ## Run 473 stress 0.07378235 
    ## ... Procrustes: rmse 0.01743815  max resid 0.05096103 
    ## Run 474 stress 0.07629236 
    ## Run 475 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 4.387327e-06  max resid 9.694821e-06 
    ## ... Similar to previous best
    ## Run 476 stress 0.07629238 
    ## Run 477 stress 0.07378229 
    ## ... Procrustes: rmse 0.01748355  max resid 0.05119218 
    ## Run 478 stress 0.07629234 
    ## Run 479 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744865  max resid 0.05106993 
    ## Run 480 stress 0.08000469 
    ## Run 481 stress 0.08288045 
    ## Run 482 stress 0.08000468 
    ## Run 483 stress 0.07378232 
    ## ... Procrustes: rmse 0.01743851  max resid 0.05104962 
    ## Run 484 stress 0.07732929 
    ## Run 485 stress 0.07365784 
    ## ... Procrustes: rmse 3.516826e-05  max resid 8.407204e-05 
    ## ... Similar to previous best
    ## Run 486 stress 0.07629237 
    ## Run 487 stress 0.07365784 
    ## ... Procrustes: rmse 5.105424e-05  max resid 0.000121252 
    ## ... Similar to previous best
    ## Run 488 stress 0.08000474 
    ## Run 489 stress 0.07365783 
    ## ... Procrustes: rmse 2.335446e-05  max resid 5.443626e-05 
    ## ... Similar to previous best
    ## Run 490 stress 0.07365792 
    ## ... Procrustes: rmse 0.000145223  max resid 0.0003469768 
    ## ... Similar to previous best
    ## Run 491 stress 0.07365786 
    ## ... Procrustes: rmse 8.694026e-05  max resid 0.0002060242 
    ## ... Similar to previous best
    ## Run 492 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174477  max resid 0.05105057 
    ## Run 493 stress 0.08000471 
    ## Run 494 stress 0.07629236 
    ## Run 495 stress 0.08288035 
    ## Run 496 stress 0.08288057 
    ## Run 497 stress 0.08000468 
    ## Run 498 stress 0.07629233 
    ## Run 499 stress 0.07365784 
    ## ... Procrustes: rmse 3.946797e-05  max resid 9.070167e-05 
    ## ... Similar to previous best
    ## Run 500 stress 0.08000471 
    ## *** Best solution repeated 7 times

``` r
# Stratified lakes and ocean sites
SD_beta_geo_SO_NMDS <- metaMDS(SD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07970526 
    ## Run 1 stress 0.0697819 
    ## ... New best solution
    ## ... Procrustes: rmse 0.09977989  max resid 0.2618694 
    ## Run 2 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01322897  max resid 0.03351544 
    ## Run 3 stress 0.07970525 
    ## Run 4 stress 0.08340289 
    ## Run 5 stress 0.07428312 
    ## Run 6 stress 0.07970527 
    ## Run 7 stress 0.07250812 
    ## Run 8 stress 0.08340293 
    ## Run 9 stress 0.07428314 
    ## Run 10 stress 0.07428313 
    ## Run 11 stress 0.07428314 
    ## Run 12 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319722  max resid 0.03320939 
    ## Run 13 stress 0.06978195 
    ## ... Procrustes: rmse 0.01317881  max resid 0.03316798 
    ## Run 14 stress 0.08340288 
    ## Run 15 stress 0.07970525 
    ## Run 16 stress 0.07428316 
    ## Run 17 stress 0.06942777 
    ## ... Procrustes: rmse 4.776959e-05  max resid 0.000123694 
    ## ... Similar to previous best
    ## Run 18 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326075  max resid 0.03334088 
    ## Run 19 stress 0.07428316 
    ## Run 20 stress 0.08340288 
    ## Run 21 stress 0.07428319 
    ## Run 22 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326572  max resid 0.03335282 
    ## Run 23 stress 0.07428317 
    ## Run 24 stress 0.07250812 
    ## Run 25 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325125  max resid 0.03332144 
    ## Run 26 stress 0.07970528 
    ## Run 27 stress 0.07428321 
    ## Run 28 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325818  max resid 0.03333875 
    ## Run 29 stress 0.08340286 
    ## Run 30 stress 0.07970526 
    ## Run 31 stress 0.06942776 
    ## ... Procrustes: rmse 2.054973e-05  max resid 5.449244e-05 
    ## ... Similar to previous best
    ## Run 32 stress 0.07428313 
    ## Run 33 stress 0.07428313 
    ## Run 34 stress 0.07428317 
    ## Run 35 stress 0.08448442 
    ## Run 36 stress 0.07428313 
    ## Run 37 stress 0.07250814 
    ## Run 38 stress 0.07970526 
    ## Run 39 stress 0.07970526 
    ## Run 40 stress 0.06978192 
    ## ... Procrustes: rmse 0.013251  max resid 0.03332212 
    ## Run 41 stress 0.08340286 
    ## Run 42 stress 0.07428315 
    ## Run 43 stress 0.07428314 
    ## Run 44 stress 0.07250813 
    ## Run 45 stress 0.08340293 
    ## Run 46 stress 0.07428316 
    ## Run 47 stress 0.07428322 
    ## Run 48 stress 0.07428313 
    ## Run 49 stress 0.07428314 
    ## Run 50 stress 0.07428319 
    ## Run 51 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 5.781181e-06  max resid 1.631511e-05 
    ## ... Similar to previous best
    ## Run 52 stress 0.07428313 
    ## Run 53 stress 0.07250816 
    ## Run 54 stress 0.07250812 
    ## Run 55 stress 0.08340288 
    ## Run 56 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326746  max resid 0.03335336 
    ## Run 57 stress 0.06942776 
    ## ... Procrustes: rmse 4.797683e-06  max resid 1.505233e-05 
    ## ... Similar to previous best
    ## Run 58 stress 0.08340287 
    ## Run 59 stress 0.07250812 
    ## Run 60 stress 0.07250812 
    ## Run 61 stress 0.08340288 
    ## Run 62 stress 0.0844844 
    ## Run 63 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326866  max resid 0.0333557 
    ## Run 64 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323422  max resid 0.03328688 
    ## Run 65 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323751  max resid 0.03329607 
    ## Run 66 stress 0.07250813 
    ## Run 67 stress 0.07844944 
    ## Run 68 stress 0.08448437 
    ## Run 69 stress 0.07970525 
    ## Run 70 stress 0.07428317 
    ## Run 71 stress 0.07970525 
    ## Run 72 stress 0.08340288 
    ## Run 73 stress 0.07970525 
    ## Run 74 stress 0.07970526 
    ## Run 75 stress 0.07250812 
    ## Run 76 stress 0.07844926 
    ## Run 77 stress 0.06978191 
    ## ... Procrustes: rmse 0.0132057  max resid 0.03322907 
    ## Run 78 stress 0.07970526 
    ## Run 79 stress 0.08340293 
    ## Run 80 stress 0.06942776 
    ## ... Procrustes: rmse 2.130086e-05  max resid 5.469436e-05 
    ## ... Similar to previous best
    ## Run 81 stress 0.08448438 
    ## Run 82 stress 0.06942776 
    ## ... Procrustes: rmse 7.023888e-06  max resid 1.876676e-05 
    ## ... Similar to previous best
    ## Run 83 stress 0.08340294 
    ## Run 84 stress 0.08448436 
    ## Run 85 stress 0.07428313 
    ## Run 86 stress 0.07844926 
    ## Run 87 stress 0.07428316 
    ## Run 88 stress 0.08340295 
    ## Run 89 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 8.261078e-07  max resid 2.660525e-06 
    ## ... Similar to previous best
    ## Run 90 stress 0.07428314 
    ## Run 91 stress 0.06978198 
    ## ... Procrustes: rmse 0.01316932  max resid 0.03314905 
    ## Run 92 stress 0.07844916 
    ## Run 93 stress 0.08340289 
    ## Run 94 stress 0.07428321 
    ## Run 95 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326729  max resid 0.03335212 
    ## Run 96 stress 0.07970525 
    ## Run 97 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323948  max resid 0.03329828 
    ## Run 98 stress 0.08340288 
    ## Run 99 stress 0.0834029 
    ## Run 100 stress 0.08340287 
    ## Run 101 stress 0.06942777 
    ## ... Procrustes: rmse 3.496827e-05  max resid 9.217361e-05 
    ## ... Similar to previous best
    ## Run 102 stress 0.06942777 
    ## ... Procrustes: rmse 1.569715e-05  max resid 4.453831e-05 
    ## ... Similar to previous best
    ## Run 103 stress 0.08340294 
    ## Run 104 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323366  max resid 0.03328463 
    ## Run 105 stress 0.07428313 
    ## Run 106 stress 0.07852811 
    ## Run 107 stress 0.08340287 
    ## Run 108 stress 0.06942777 
    ## ... Procrustes: rmse 2.259445e-05  max resid 6.17338e-05 
    ## ... Similar to previous best
    ## Run 109 stress 0.07428313 
    ## Run 110 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324814  max resid 0.03331493 
    ## Run 111 stress 0.08340287 
    ## Run 112 stress 0.0742832 
    ## Run 113 stress 0.08340287 
    ## Run 114 stress 0.07428314 
    ## Run 115 stress 0.07428314 
    ## Run 116 stress 0.06942776 
    ## ... Procrustes: rmse 1.065852e-05  max resid 2.201927e-05 
    ## ... Similar to previous best
    ## Run 117 stress 0.0834029 
    ## Run 118 stress 0.07970525 
    ## Run 119 stress 0.08340287 
    ## Run 120 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326277  max resid 0.033346 
    ## Run 121 stress 0.07428313 
    ## Run 122 stress 0.07970525 
    ## Run 123 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321936  max resid 0.03325733 
    ## Run 124 stress 0.08340287 
    ## Run 125 stress 0.07428319 
    ## Run 126 stress 0.07970525 
    ## Run 127 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326409  max resid 0.03334914 
    ## Run 128 stress 0.06978193 
    ## ... Procrustes: rmse 0.01324992  max resid 0.03332149 
    ## Run 129 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326718  max resid 0.03335523 
    ## Run 130 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325049  max resid 0.0333206 
    ## Run 131 stress 0.07428313 
    ## Run 132 stress 0.07428318 
    ## Run 133 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326057  max resid 0.03334102 
    ## Run 134 stress 0.07970525 
    ## Run 135 stress 0.07428313 
    ## Run 136 stress 0.07428319 
    ## Run 137 stress 0.07428313 
    ## Run 138 stress 0.08340286 
    ## Run 139 stress 0.06942776 
    ## ... Procrustes: rmse 3.852501e-06  max resid 1.038145e-05 
    ## ... Similar to previous best
    ## Run 140 stress 0.07428313 
    ## Run 141 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327362  max resid 0.03336993 
    ## Run 142 stress 0.07428316 
    ## Run 143 stress 0.08448446 
    ## Run 144 stress 0.06942776 
    ## ... Procrustes: rmse 1.84881e-05  max resid 4.764955e-05 
    ## ... Similar to previous best
    ## Run 145 stress 0.07970526 
    ## Run 146 stress 0.07428314 
    ## Run 147 stress 0.08340287 
    ## Run 148 stress 0.06942777 
    ## ... Procrustes: rmse 4.90752e-05  max resid 0.0001261711 
    ## ... Similar to previous best
    ## Run 149 stress 0.06942776 
    ## ... Procrustes: rmse 1.725068e-05  max resid 4.417829e-05 
    ## ... Similar to previous best
    ## Run 150 stress 0.07844946 
    ## Run 151 stress 0.08340287 
    ## Run 152 stress 0.08340286 
    ## Run 153 stress 0.08340289 
    ## Run 154 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319676  max resid 0.03320608 
    ## Run 155 stress 0.08340287 
    ## Run 156 stress 0.07428313 
    ## Run 157 stress 0.07250812 
    ## Run 158 stress 0.06942777 
    ## ... Procrustes: rmse 1.837714e-05  max resid 5.296916e-05 
    ## ... Similar to previous best
    ## Run 159 stress 0.07250813 
    ## Run 160 stress 0.07250813 
    ## Run 161 stress 0.07844954 
    ## Run 162 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323455  max resid 0.03328598 
    ## Run 163 stress 0.08340291 
    ## Run 164 stress 0.06942778 
    ## ... Procrustes: rmse 1.761383e-05  max resid 3.401794e-05 
    ## ... Similar to previous best
    ## Run 165 stress 0.07970526 
    ## Run 166 stress 0.08340289 
    ## Run 167 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326591  max resid 0.03335365 
    ## Run 168 stress 0.07250816 
    ## Run 169 stress 0.08448442 
    ## Run 170 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325067  max resid 0.03332047 
    ## Run 171 stress 0.06978198 
    ## ... Procrustes: rmse 0.0131708  max resid 0.03315245 
    ## Run 172 stress 0.07970526 
    ## Run 173 stress 0.07428314 
    ## Run 174 stress 0.07250812 
    ## Run 175 stress 0.08448438 
    ## Run 176 stress 0.07428312 
    ## Run 177 stress 0.06942776 
    ## ... Procrustes: rmse 1.743952e-05  max resid 4.500145e-05 
    ## ... Similar to previous best
    ## Run 178 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325372  max resid 0.03332916 
    ## Run 179 stress 0.08340299 
    ## Run 180 stress 0.07970526 
    ## Run 181 stress 0.07428315 
    ## Run 182 stress 0.07844926 
    ## Run 183 stress 0.07250812 
    ## Run 184 stress 0.07428314 
    ## Run 185 stress 0.07970525 
    ## Run 186 stress 0.07428313 
    ## Run 187 stress 0.06978196 
    ## ... Procrustes: rmse 0.01317467  max resid 0.03315715 
    ## Run 188 stress 0.07970526 
    ## Run 189 stress 0.08340289 
    ## Run 190 stress 0.06942776 
    ## ... Procrustes: rmse 1.230671e-05  max resid 3.170555e-05 
    ## ... Similar to previous best
    ## Run 191 stress 0.07428313 
    ## Run 192 stress 0.06942776 
    ## ... Procrustes: rmse 1.623918e-05  max resid 3.984474e-05 
    ## ... Similar to previous best
    ## Run 193 stress 0.07250812 
    ## Run 194 stress 0.06942777 
    ## ... Procrustes: rmse 3.986469e-05  max resid 0.0001029264 
    ## ... Similar to previous best
    ## Run 195 stress 0.07970525 
    ## Run 196 stress 0.07970525 
    ## Run 197 stress 0.06942777 
    ## ... Procrustes: rmse 4.659103e-05  max resid 0.0001197577 
    ## ... Similar to previous best
    ## Run 198 stress 0.07970525 
    ## Run 199 stress 0.07250812 
    ## Run 200 stress 0.08448435 
    ## Run 201 stress 0.08340298 
    ## Run 202 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327438  max resid 0.0333715 
    ## Run 203 stress 0.07970526 
    ## Run 204 stress 0.07428313 
    ## Run 205 stress 0.08340286 
    ## Run 206 stress 0.07844948 
    ## Run 207 stress 0.07428314 
    ## Run 208 stress 0.07428314 
    ## Run 209 stress 0.07428317 
    ## Run 210 stress 0.07428314 
    ## Run 211 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323286  max resid 0.03328473 
    ## Run 212 stress 0.07250812 
    ## Run 213 stress 0.0834029 
    ## Run 214 stress 0.07970526 
    ## Run 215 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322788  max resid 0.03327181 
    ## Run 216 stress 0.08340294 
    ## Run 217 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326577  max resid 0.03335098 
    ## Run 218 stress 0.07250815 
    ## Run 219 stress 0.07250813 
    ## Run 220 stress 0.07428317 
    ## Run 221 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322941  max resid 0.03327765 
    ## Run 222 stress 0.07428314 
    ## Run 223 stress 0.08340287 
    ## Run 224 stress 0.08340293 
    ## Run 225 stress 0.06978191 
    ## ... Procrustes: rmse 0.0132109  max resid 0.03323321 
    ## Run 226 stress 0.08448444 
    ## Run 227 stress 0.083403 
    ## Run 228 stress 0.07428315 
    ## Run 229 stress 0.07428312 
    ## Run 230 stress 0.07970525 
    ## Run 231 stress 0.08340289 
    ## Run 232 stress 0.07428314 
    ## Run 233 stress 0.07970525 
    ## Run 234 stress 0.07250813 
    ## Run 235 stress 0.07970525 
    ## Run 236 stress 0.06942777 
    ## ... Procrustes: rmse 4.651459e-05  max resid 0.0001196733 
    ## ... Similar to previous best
    ## Run 237 stress 0.07970526 
    ## Run 238 stress 0.07428321 
    ## Run 239 stress 0.07970525 
    ## Run 240 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325305  max resid 0.03332691 
    ## Run 241 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323873  max resid 0.03329478 
    ## Run 242 stress 0.06978201 
    ## ... Procrustes: rmse 0.0132848  max resid 0.03339036 
    ## Run 243 stress 0.08448446 
    ## Run 244 stress 0.07428318 
    ## Run 245 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327024  max resid 0.03336081 
    ## Run 246 stress 0.06942776 
    ## ... Procrustes: rmse 1.234082e-05  max resid 3.175224e-05 
    ## ... Similar to previous best
    ## Run 247 stress 0.07970525 
    ## Run 248 stress 0.08448445 
    ## Run 249 stress 0.07428313 
    ## Run 250 stress 0.07428315 
    ## Run 251 stress 0.07970526 
    ## Run 252 stress 0.06978204 
    ## ... Procrustes: rmse 0.01328345  max resid 0.03338383 
    ## Run 253 stress 0.06978194 
    ## ... Procrustes: rmse 0.01323498  max resid 0.03328383 
    ## Run 254 stress 0.07250813 
    ## Run 255 stress 0.07428318 
    ## Run 256 stress 0.07250812 
    ## Run 257 stress 0.07428313 
    ## Run 258 stress 0.07970526 
    ## Run 259 stress 0.07250812 
    ## Run 260 stress 0.07970525 
    ## Run 261 stress 0.07428313 
    ## Run 262 stress 0.06942776 
    ## ... Procrustes: rmse 9.470404e-06  max resid 2.468648e-05 
    ## ... Similar to previous best
    ## Run 263 stress 0.07970525 
    ## Run 264 stress 0.08448445 
    ## Run 265 stress 0.07428318 
    ## Run 266 stress 0.06942776 
    ## ... Procrustes: rmse 1.110876e-05  max resid 2.384522e-05 
    ## ... Similar to previous best
    ## Run 267 stress 0.07970526 
    ## Run 268 stress 0.07970525 
    ## Run 269 stress 0.0742832 
    ## Run 270 stress 0.07250813 
    ## Run 271 stress 0.07428313 
    ## Run 272 stress 0.06942776 
    ## ... Procrustes: rmse 3.140554e-05  max resid 8.092545e-05 
    ## ... Similar to previous best
    ## Run 273 stress 0.07428317 
    ## Run 274 stress 0.07250813 
    ## Run 275 stress 0.07970525 
    ## Run 276 stress 0.07428314 
    ## Run 277 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324235  max resid 0.03330283 
    ## Run 278 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324894  max resid 0.03331629 
    ## Run 279 stress 0.08340288 
    ## Run 280 stress 0.07970525 
    ## Run 281 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326525  max resid 0.0333517 
    ## Run 282 stress 0.07428313 
    ## Run 283 stress 0.07250812 
    ## Run 284 stress 0.07250812 
    ## Run 285 stress 0.07428316 
    ## Run 286 stress 0.08340291 
    ## Run 287 stress 0.06978198 
    ## ... Procrustes: rmse 0.0131694  max resid 0.03314895 
    ## Run 288 stress 0.07428313 
    ## Run 289 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325369  max resid 0.03332686 
    ## Run 290 stress 0.06942776 
    ## ... Procrustes: rmse 2.229292e-05  max resid 5.782549e-05 
    ## ... Similar to previous best
    ## Run 291 stress 0.06942776 
    ## ... Procrustes: rmse 1.568389e-05  max resid 4.04328e-05 
    ## ... Similar to previous best
    ## Run 292 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325443  max resid 0.03332802 
    ## Run 293 stress 0.07250812 
    ## Run 294 stress 0.08340286 
    ## Run 295 stress 0.06942776 
    ## ... Procrustes: rmse 1.544291e-05  max resid 3.981741e-05 
    ## ... Similar to previous best
    ## Run 296 stress 0.07250814 
    ## Run 297 stress 0.06942777 
    ## ... Procrustes: rmse 1.397479e-05  max resid 4.385131e-05 
    ## ... Similar to previous best
    ## Run 298 stress 0.06942777 
    ## ... Procrustes: rmse 4.104928e-05  max resid 0.0001056589 
    ## ... Similar to previous best
    ## Run 299 stress 0.06978194 
    ## ... Procrustes: rmse 0.01320141  max resid 0.03322074 
    ## Run 300 stress 0.07250812 
    ## Run 301 stress 0.07428314 
    ## Run 302 stress 0.07250814 
    ## Run 303 stress 0.07970526 
    ## Run 304 stress 0.08448441 
    ## Run 305 stress 0.07428315 
    ## Run 306 stress 0.07428313 
    ## Run 307 stress 0.06942776 
    ## ... Procrustes: rmse 1.881312e-05  max resid 4.849362e-05 
    ## ... Similar to previous best
    ## Run 308 stress 0.07428313 
    ## Run 309 stress 0.07250813 
    ## Run 310 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325273  max resid 0.03332492 
    ## Run 311 stress 0.06978191 
    ## ... Procrustes: rmse 0.01322585  max resid 0.03326787 
    ## Run 312 stress 0.07428313 
    ## Run 313 stress 0.08340298 
    ## Run 314 stress 0.07250812 
    ## Run 315 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324541  max resid 0.03331105 
    ## Run 316 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323309  max resid 0.03328374 
    ## Run 317 stress 0.08340288 
    ## Run 318 stress 0.08340286 
    ## Run 319 stress 0.07250812 
    ## Run 320 stress 0.07250812 
    ## Run 321 stress 0.08448454 
    ## Run 322 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322852  max resid 0.03327821 
    ## Run 323 stress 0.08340287 
    ## Run 324 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327628  max resid 0.03337346 
    ## Run 325 stress 0.07844931 
    ## Run 326 stress 0.07250812 
    ## Run 327 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324453  max resid 0.03330597 
    ## Run 328 stress 0.07428313 
    ## Run 329 stress 0.08340287 
    ## Run 330 stress 0.08340291 
    ## Run 331 stress 0.06942777 
    ## ... Procrustes: rmse 4.27912e-05  max resid 0.0001102688 
    ## ... Similar to previous best
    ## Run 332 stress 0.07428313 
    ## Run 333 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322405  max resid 0.03326506 
    ## Run 334 stress 0.07250812 
    ## Run 335 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323509  max resid 0.03328706 
    ## Run 336 stress 0.069782 
    ## ... Procrustes: rmse 0.01328229  max resid 0.03338413 
    ## Run 337 stress 0.08448439 
    ## Run 338 stress 0.07428316 
    ## Run 339 stress 0.07428316 
    ## Run 340 stress 0.07970525 
    ## Run 341 stress 0.06978194 
    ## ... Procrustes: rmse 0.01318806  max resid 0.03318493 
    ## Run 342 stress 0.06978202 
    ## ... Procrustes: rmse 0.01328679  max resid 0.03339607 
    ## Run 343 stress 0.07428313 
    ## Run 344 stress 0.06942776 
    ## ... Procrustes: rmse 2.813781e-06  max resid 7.230612e-06 
    ## ... Similar to previous best
    ## Run 345 stress 0.07428313 
    ## Run 346 stress 0.07250812 
    ## Run 347 stress 0.08340291 
    ## Run 348 stress 0.08340287 
    ## Run 349 stress 0.07844952 
    ## Run 350 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326176  max resid 0.03334188 
    ## Run 351 stress 0.06942776 
    ## ... Procrustes: rmse 3.661703e-06  max resid 1.041387e-05 
    ## ... Similar to previous best
    ## Run 352 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 7.99961e-06  max resid 2.064138e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.07428313 
    ## Run 354 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326625  max resid 0.03335162 
    ## Run 355 stress 0.06942776 
    ## ... Procrustes: rmse 1.009031e-05  max resid 2.645616e-05 
    ## ... Similar to previous best
    ## Run 356 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326229  max resid 0.03334523 
    ## Run 357 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 7.488588e-07  max resid 2.08318e-06 
    ## ... Similar to previous best
    ## Run 358 stress 0.07428314 
    ## Run 359 stress 0.08340291 
    ## Run 360 stress 0.08340293 
    ## Run 361 stress 0.06942777 
    ## ... Procrustes: rmse 4.330491e-05  max resid 0.0001117947 
    ## ... Similar to previous best
    ## Run 362 stress 0.0834029 
    ## Run 363 stress 0.06942779 
    ## ... Procrustes: rmse 7.767124e-05  max resid 0.0002002115 
    ## ... Similar to previous best
    ## Run 364 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324312  max resid 0.03330387 
    ## Run 365 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319327  max resid 0.03320011 
    ## Run 366 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320293  max resid 0.03322004 
    ## Run 367 stress 0.06978193 
    ## ... Procrustes: rmse 0.01323767  max resid 0.03329108 
    ## Run 368 stress 0.08340293 
    ## Run 369 stress 0.07250815 
    ## Run 370 stress 0.07250812 
    ## Run 371 stress 0.07970526 
    ## Run 372 stress 0.07970525 
    ## Run 373 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326667  max resid 0.03335265 
    ## Run 374 stress 0.08340291 
    ## Run 375 stress 0.07428315 
    ## Run 376 stress 0.07970526 
    ## Run 377 stress 0.07250812 
    ## Run 378 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 4.195903e-06  max resid 1.027609e-05 
    ## ... Similar to previous best
    ## Run 379 stress 0.07250812 
    ## Run 380 stress 0.06942776 
    ## ... Procrustes: rmse 2.950024e-06  max resid 6.649896e-06 
    ## ... Similar to previous best
    ## Run 381 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326551  max resid 0.0333482 
    ## Run 382 stress 0.06978197 
    ## ... Procrustes: rmse 0.01318195  max resid 0.03318047 
    ## Run 383 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320634  max resid 0.03322521 
    ## Run 384 stress 0.07428313 
    ## Run 385 stress 0.07970525 
    ## Run 386 stress 0.06942776 
    ## ... Procrustes: rmse 1.313195e-05  max resid 3.651663e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.07428313 
    ## Run 388 stress 0.07250812 
    ## Run 389 stress 0.07844948 
    ## Run 390 stress 0.07428324 
    ## Run 391 stress 0.07428314 
    ## Run 392 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323518  max resid 0.03328515 
    ## Run 393 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324899  max resid 0.03331882 
    ## Run 394 stress 0.07428315 
    ## Run 395 stress 0.06942776 
    ## ... Procrustes: rmse 8.778275e-06  max resid 2.831027e-05 
    ## ... Similar to previous best
    ## Run 396 stress 0.06942776 
    ## ... Procrustes: rmse 4.686668e-06  max resid 1.11265e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.08340287 
    ## Run 398 stress 0.07428314 
    ## Run 399 stress 0.07428313 
    ## Run 400 stress 0.07250812 
    ## Run 401 stress 0.07250812 
    ## Run 402 stress 0.07970525 
    ## Run 403 stress 0.06942778 
    ## ... Procrustes: rmse 3.464503e-05  max resid 9.220544e-05 
    ## ... Similar to previous best
    ## Run 404 stress 0.06978199 
    ## ... Procrustes: rmse 0.01327842  max resid 0.03337389 
    ## Run 405 stress 0.07250814 
    ## Run 406 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326869  max resid 0.03335198 
    ## Run 407 stress 0.07970526 
    ## Run 408 stress 0.06942778 
    ## ... Procrustes: rmse 5.666213e-05  max resid 0.0001449776 
    ## ... Similar to previous best
    ## Run 409 stress 0.08340291 
    ## Run 410 stress 0.08340286 
    ## Run 411 stress 0.07970528 
    ## Run 412 stress 0.06942778 
    ## ... Procrustes: rmse 5.580057e-05  max resid 0.0001438077 
    ## ... Similar to previous best
    ## Run 413 stress 0.08340289 
    ## Run 414 stress 0.07428312 
    ## Run 415 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327308  max resid 0.03336546 
    ## Run 416 stress 0.08340286 
    ## Run 417 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327153  max resid 0.03336361 
    ## Run 418 stress 0.07428315 
    ## Run 419 stress 0.06978198 
    ## ... Procrustes: rmse 0.0131699  max resid 0.03314482 
    ## Run 420 stress 0.07970525 
    ## Run 421 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325724  max resid 0.03333481 
    ## Run 422 stress 0.08340289 
    ## Run 423 stress 0.07970525 
    ## Run 424 stress 0.06978199 
    ## ... Procrustes: rmse 0.0132722  max resid 0.03336791 
    ## Run 425 stress 0.07250812 
    ## Run 426 stress 0.06942777 
    ## ... Procrustes: rmse 1.682436e-05  max resid 4.552759e-05 
    ## ... Similar to previous best
    ## Run 427 stress 0.07970525 
    ## Run 428 stress 0.07970526 
    ## Run 429 stress 0.07428313 
    ## Run 430 stress 0.07970525 
    ## Run 431 stress 0.08448438 
    ## Run 432 stress 0.07970526 
    ## Run 433 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322961  max resid 0.03327503 
    ## Run 434 stress 0.08340287 
    ## Run 435 stress 0.07970525 
    ## Run 436 stress 0.07428313 
    ## Run 437 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323662  max resid 0.03329051 
    ## Run 438 stress 0.07970525 
    ## Run 439 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324814  max resid 0.03331331 
    ## Run 440 stress 0.08340287 
    ## Run 441 stress 0.07428315 
    ## Run 442 stress 0.07970526 
    ## Run 443 stress 0.08340287 
    ## Run 444 stress 0.06942776 
    ## ... Procrustes: rmse 1.507565e-05  max resid 4.025654e-05 
    ## ... Similar to previous best
    ## Run 445 stress 0.07428313 
    ## Run 446 stress 0.07250813 
    ## Run 447 stress 0.06978192 
    ## ... Procrustes: rmse 0.01323321  max resid 0.03328661 
    ## Run 448 stress 0.06978199 
    ## ... Procrustes: rmse 0.01326903  max resid 0.03335399 
    ## Run 449 stress 0.07970526 
    ## Run 450 stress 0.06942776 
    ## ... Procrustes: rmse 1.437596e-05  max resid 3.824458e-05 
    ## ... Similar to previous best
    ## Run 451 stress 0.07428313 
    ## Run 452 stress 0.07428312 
    ## Run 453 stress 0.07970525 
    ## Run 454 stress 0.07844946 
    ## Run 455 stress 0.07428313 
    ## Run 456 stress 0.0844844 
    ## Run 457 stress 0.07970526 
    ## Run 458 stress 0.07970525 
    ## Run 459 stress 0.08340291 
    ## Run 460 stress 0.07428319 
    ## Run 461 stress 0.07970525 
    ## Run 462 stress 0.08448444 
    ## Run 463 stress 0.06978199 
    ## ... Procrustes: rmse 0.01316384  max resid 0.03313439 
    ## Run 464 stress 0.07428315 
    ## Run 465 stress 0.06942776 
    ## ... Procrustes: rmse 1.660861e-05  max resid 4.252724e-05 
    ## ... Similar to previous best
    ## Run 466 stress 0.08340286 
    ## Run 467 stress 0.06978203 
    ## ... Procrustes: rmse 0.01315642  max resid 0.03311868 
    ## Run 468 stress 0.07250812 
    ## Run 469 stress 0.07970525 
    ## Run 470 stress 0.08448458 
    ## Run 471 stress 0.069782 
    ## ... Procrustes: rmse 0.01316347  max resid 0.03313357 
    ## Run 472 stress 0.06942776 
    ## ... Procrustes: rmse 6.584953e-06  max resid 1.960466e-05 
    ## ... Similar to previous best
    ## Run 473 stress 0.07250813 
    ## Run 474 stress 0.06942776 
    ## ... Procrustes: rmse 1.27441e-05  max resid 3.3054e-05 
    ## ... Similar to previous best
    ## Run 475 stress 0.07844956 
    ## Run 476 stress 0.07428314 
    ## Run 477 stress 0.07250813 
    ## Run 478 stress 0.08448443 
    ## Run 479 stress 0.07970526 
    ## Run 480 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324135  max resid 0.03329872 
    ## Run 481 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323504  max resid 0.03328623 
    ## Run 482 stress 0.08340289 
    ## Run 483 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322501  max resid 0.03326695 
    ## Run 484 stress 0.07250812 
    ## Run 485 stress 0.07970525 
    ## Run 486 stress 0.08340286 
    ## Run 487 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326196  max resid 0.03334338 
    ## Run 488 stress 0.06978193 
    ## ... Procrustes: rmse 0.01324914  max resid 0.03331396 
    ## Run 489 stress 0.07970525 
    ## Run 490 stress 0.07428314 
    ## Run 491 stress 0.06942776 
    ## ... Procrustes: rmse 8.928776e-06  max resid 2.324675e-05 
    ## ... Similar to previous best
    ## Run 492 stress 0.07428321 
    ## Run 493 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321239  max resid 0.03323834 
    ## Run 494 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326749  max resid 0.03335461 
    ## Run 495 stress 0.08340287 
    ## Run 496 stress 0.07428321 
    ## Run 497 stress 0.07970525 
    ## Run 498 stress 0.0742832 
    ## Run 499 stress 0.07428313 
    ## Run 500 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322903  max resid 0.0332741 
    ## *** Best solution repeated 15 times

``` r
# Mixed lakes
SD_beta_geo_M_NMDS <- metaMDS(SD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.001279543 
    ## Run 1 stress 0.1491957 
    ## Run 2 stress 0.1491957 
    ## Run 3 stress 9.53108e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04310352  max resid 0.05938471 
    ## Run 4 stress 0.1491957 
    ## Run 5 stress 0.0003642748 
    ## ... Procrustes: rmse 0.0139761  max resid 0.01900139 
    ## Run 6 stress 9.840878e-05 
    ## ... Procrustes: rmse 0.000211579  max resid 0.0004254032 
    ## ... Similar to previous best
    ## Run 7 stress 0.1990774 
    ## Run 8 stress 0.1491957 
    ## Run 9 stress 8.820077e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001323144  max resid 0.0002316592 
    ## ... Similar to previous best
    ## Run 10 stress 9.723146e-05 
    ## ... Procrustes: rmse 0.0002056977  max resid 0.0004147001 
    ## ... Similar to previous best
    ## Run 11 stress 0.2570109 
    ## Run 12 stress 0.001173459 
    ## Run 13 stress 0.1491957 
    ## Run 14 stress 0.001214688 
    ## Run 15 stress 0.0004574326 
    ## ... Procrustes: rmse 0.01561245  max resid 0.02137042 
    ## Run 16 stress 0.0004870663 
    ## ... Procrustes: rmse 0.01609855  max resid 0.02204033 
    ## Run 17 stress 9.245134e-05 
    ## ... Procrustes: rmse 2.237011e-05  max resid 3.973051e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.0004792971 
    ## ... Procrustes: rmse 0.01595305  max resid 0.02184059 
    ## Run 19 stress 9.234381e-05 
    ## ... Procrustes: rmse 0.0001819672  max resid 0.0003566155 
    ## ... Similar to previous best
    ## Run 20 stress 0.1491957 
    ## Run 21 stress 0.1990774 
    ## Run 22 stress 9.884964e-05 
    ## ... Procrustes: rmse 0.0001730292  max resid 0.0003848319 
    ## ... Similar to previous best
    ## Run 23 stress 9.538217e-05 
    ## ... Procrustes: rmse 0.0001326676  max resid 0.000221409 
    ## ... Similar to previous best
    ## Run 24 stress 0.0004854673 
    ## ... Procrustes: rmse 0.02657591  max resid 0.03654198 
    ## Run 25 stress 8.640689e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001220787  max resid 0.0002154569 
    ## ... Similar to previous best
    ## Run 26 stress 9.19687e-05 
    ## ... Procrustes: rmse 1.2428e-05  max resid 2.099469e-05 
    ## ... Similar to previous best
    ## Run 27 stress 0.0001155055 
    ## ... Procrustes: rmse 0.007854421  max resid 0.01056277 
    ## Run 28 stress 0.001313622 
    ## Run 29 stress 9.444895e-05 
    ## ... Procrustes: rmse 0.0001721722  max resid 0.0002688388 
    ## ... Similar to previous best
    ## Run 30 stress 8.801874e-05 
    ## ... Procrustes: rmse 0.0002192989  max resid 0.0004081751 
    ## ... Similar to previous best
    ## Run 31 stress 0.001216334 
    ## Run 32 stress 9.798545e-05 
    ## ... Procrustes: rmse 0.006561293  max resid 0.008793825 
    ## Run 33 stress 9.376043e-05 
    ## ... Procrustes: rmse 0.0002410645  max resid 0.0004743237 
    ## ... Similar to previous best
    ## Run 34 stress 0.0002682826 
    ## ... Procrustes: rmse 0.01200731  max resid 0.01629423 
    ## Run 35 stress 0.1491957 
    ## Run 36 stress 0.001426281 
    ## Run 37 stress 0.1491957 
    ## Run 38 stress 0.2696744 
    ## Run 39 stress 0.1990774 
    ## Run 40 stress 0.0004750397 
    ## ... Procrustes: rmse 0.01596555  max resid 0.02175519 
    ## Run 41 stress 0.1491957 
    ## Run 42 stress 0.0004898902 
    ## ... Procrustes: rmse 0.02656568  max resid 0.03661427 
    ## Run 43 stress 9.702653e-05 
    ## ... Procrustes: rmse 0.0001921347  max resid 0.0003163885 
    ## ... Similar to previous best
    ## Run 44 stress 0.1990774 
    ## Run 45 stress 0.0005351838 
    ## ... Procrustes: rmse 0.02782396  max resid 0.03835918 
    ## Run 46 stress 9.729036e-05 
    ## ... Procrustes: rmse 7.663038e-05  max resid 0.0001461713 
    ## ... Similar to previous best
    ## Run 47 stress 9.906375e-05 
    ## ... Procrustes: rmse 0.0002111744  max resid 0.0004652847 
    ## ... Similar to previous best
    ## Run 48 stress 0.001429239 
    ## Run 49 stress 0.2848519 
    ## Run 50 stress 0.001287104 
    ## Run 51 stress 9.001198e-05 
    ## ... Procrustes: rmse 0.0001176704  max resid 0.0002038861 
    ## ... Similar to previous best
    ## Run 52 stress 0.0001702881 
    ## ... Procrustes: rmse 0.01564378  max resid 0.02150104 
    ## Run 53 stress 0.2832379 
    ## Run 54 stress 8.224797e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001780881  max resid 0.0002614928 
    ## ... Similar to previous best
    ## Run 55 stress 0.1990774 
    ## Run 56 stress 9.612724e-05 
    ## ... Procrustes: rmse 0.006277782  max resid 0.00867115 
    ## Run 57 stress 0.1990774 
    ## Run 58 stress 9.860423e-05 
    ## ... Procrustes: rmse 0.000201431  max resid 0.000296746 
    ## ... Similar to previous best
    ## Run 59 stress 0.001301194 
    ## Run 60 stress 0.00113917 
    ## Run 61 stress 0.0004920427 
    ## ... Procrustes: rmse 0.01619995  max resid 0.02234653 
    ## Run 62 stress 9.221589e-05 
    ## ... Procrustes: rmse 0.000117444  max resid 0.0001906516 
    ## ... Similar to previous best
    ## Run 63 stress 0.2842805 
    ## Run 64 stress 9.241881e-05 
    ## ... Procrustes: rmse 0.0001818296  max resid 0.0003356979 
    ## ... Similar to previous best
    ## Run 65 stress 0.001307369 
    ## Run 66 stress 0.1491957 
    ## Run 67 stress 0.1776012 
    ## Run 68 stress 0.0005970882 
    ## Run 69 stress 0.001195124 
    ## Run 70 stress 0.001318126 
    ## Run 71 stress 0.1491957 
    ## Run 72 stress 0.001351274 
    ## Run 73 stress 0.001084321 
    ## Run 74 stress 9.478679e-05 
    ## ... Procrustes: rmse 0.0001094466  max resid 0.0001825096 
    ## ... Similar to previous best
    ## Run 75 stress 0.001288921 
    ## Run 76 stress 0.0004379742 
    ## ... Procrustes: rmse 0.01527699  max resid 0.02107097 
    ## Run 77 stress 0.001179846 
    ## Run 78 stress 0.000741062 
    ## Run 79 stress 0.00102263 
    ## Run 80 stress 0.001496165 
    ## Run 81 stress 8.422682e-05 
    ## ... Procrustes: rmse 0.0001091774  max resid 0.0001874806 
    ## ... Similar to previous best
    ## Run 82 stress 0.001333055 
    ## Run 83 stress 0.001201387 
    ## Run 84 stress 9.524589e-05 
    ## ... Procrustes: rmse 0.0001765914  max resid 0.0002571895 
    ## ... Similar to previous best
    ## Run 85 stress 0.2848513 
    ## Run 86 stress 0.1990774 
    ## Run 87 stress 8.335915e-05 
    ## ... Procrustes: rmse 0.0001556945  max resid 0.0002540357 
    ## ... Similar to previous best
    ## Run 88 stress 0.001280946 
    ## Run 89 stress 0.1491957 
    ## Run 90 stress 9.685376e-05 
    ## ... Procrustes: rmse 0.006766541  max resid 0.009339909 
    ## Run 91 stress 0.0003755995 
    ## ... Procrustes: rmse 0.01371825  max resid 0.01892328 
    ## Run 92 stress 0.1990774 
    ## Run 93 stress 0.001167894 
    ## Run 94 stress 0.001295086 
    ## Run 95 stress 9.078642e-05 
    ## ... Procrustes: rmse 0.0001120792  max resid 0.000225936 
    ## ... Similar to previous best
    ## Run 96 stress 0.001404523 
    ## Run 97 stress 0.1776018 
    ## Run 98 stress 0.001434831 
    ## Run 99 stress 0.001246167 
    ## Run 100 stress 9.763472e-05 
    ## ... Procrustes: rmse 0.0001838043  max resid 0.000242491 
    ## ... Similar to previous best
    ## Run 101 stress 0.001323602 
    ## Run 102 stress 9.47286e-05 
    ## ... Procrustes: rmse 0.005228691  max resid 0.007228123 
    ## Run 103 stress 0.001285392 
    ## Run 104 stress 8.822044e-05 
    ## ... Procrustes: rmse 0.0001799779  max resid 0.0002695646 
    ## ... Similar to previous best
    ## Run 105 stress 0.001251599 
    ## Run 106 stress 9.454193e-05 
    ## ... Procrustes: rmse 0.0005032218  max resid 0.0008888427 
    ## ... Similar to previous best
    ## Run 107 stress 0.001262287 
    ## Run 108 stress 0.0004895937 
    ## ... Procrustes: rmse 0.01607143  max resid 0.02218775 
    ## Run 109 stress 0.001422513 
    ## Run 110 stress 0.1990774 
    ## Run 111 stress 9.95149e-05 
    ## ... Procrustes: rmse 0.0001802808  max resid 0.0002387128 
    ## ... Similar to previous best
    ## Run 112 stress 0.3083098 
    ## Run 113 stress 0.1776014 
    ## Run 114 stress 0.1776019 
    ## Run 115 stress 0.00130846 
    ## Run 116 stress 8.80237e-05 
    ## ... Procrustes: rmse 0.0001751522  max resid 0.0002593455 
    ## ... Similar to previous best
    ## Run 117 stress 0.001268236 
    ## Run 118 stress 0.001392935 
    ## Run 119 stress 9.844173e-05 
    ## ... Procrustes: rmse 0.01074631  max resid 0.01500588 
    ## Run 120 stress 0.001129338 
    ## Run 121 stress 8.672977e-05 
    ## ... Procrustes: rmse 0.0001099236  max resid 0.0001727857 
    ## ... Similar to previous best
    ## Run 122 stress 0.0004736601 
    ## ... Procrustes: rmse 0.01587125  max resid 0.02189142 
    ## Run 123 stress 9.506439e-05 
    ## ... Procrustes: rmse 0.0001437024  max resid 0.0002298794 
    ## ... Similar to previous best
    ## Run 124 stress 0.0004516523 
    ## ... Procrustes: rmse 0.01551769  max resid 0.02140319 
    ## Run 125 stress 0.3083098 
    ## Run 126 stress 0.0004225034 
    ## ... Procrustes: rmse 0.01500473  max resid 0.02069556 
    ## Run 127 stress 0.001353675 
    ## Run 128 stress 0.1990774 
    ## Run 129 stress 8.797733e-05 
    ## ... Procrustes: rmse 0.005601264  max resid 0.007738743 
    ## Run 130 stress 9.674669e-05 
    ## ... Procrustes: rmse 0.0001237075  max resid 0.0001931713 
    ## ... Similar to previous best
    ## Run 131 stress 0.0004930326 
    ## ... Procrustes: rmse 0.01620788  max resid 0.02235436 
    ## Run 132 stress 9.530576e-05 
    ## ... Procrustes: rmse 0.0001399223  max resid 0.0002484124 
    ## ... Similar to previous best
    ## Run 133 stress 0.001336326 
    ## Run 134 stress 0.1776016 
    ## Run 135 stress 0.0004384579 
    ## ... Procrustes: rmse 0.01525166  max resid 0.02103616 
    ## Run 136 stress 0.001294809 
    ## Run 137 stress 0.1990774 
    ## Run 138 stress 0.0007778841 
    ## Run 139 stress 9.928553e-05 
    ## ... Procrustes: rmse 0.0001528737  max resid 0.0002425635 
    ## ... Similar to previous best
    ## Run 140 stress 0.1491957 
    ## Run 141 stress 0.2852153 
    ## Run 142 stress 9.064906e-05 
    ## ... Procrustes: rmse 0.0001467574  max resid 0.0002570453 
    ## ... Similar to previous best
    ## Run 143 stress 9.418258e-05 
    ## ... Procrustes: rmse 0.0001848641  max resid 0.0003390455 
    ## ... Similar to previous best
    ## Run 144 stress 0.001283314 
    ## Run 145 stress 0.1491957 
    ## Run 146 stress 0.001323202 
    ## Run 147 stress 9.195532e-05 
    ## ... Procrustes: rmse 0.0001442668  max resid 0.0002285337 
    ## ... Similar to previous best
    ## Run 148 stress 0.0005074327 
    ## ... Procrustes: rmse 0.01645211  max resid 0.02269217 
    ## Run 149 stress 0.1491957 
    ## Run 150 stress 0.001388783 
    ## Run 151 stress 0.177601 
    ## Run 152 stress 0.1491957 
    ## Run 153 stress 0.2494706 
    ## Run 154 stress 9.135826e-05 
    ## ... Procrustes: rmse 0.0001499638  max resid 0.0002566944 
    ## ... Similar to previous best
    ## Run 155 stress 0.001465498 
    ## Run 156 stress 0.3083098 
    ## Run 157 stress 0.00104854 
    ## Run 158 stress 0.0009974647 
    ## Run 159 stress 9.434032e-05 
    ## ... Procrustes: rmse 0.0001843604  max resid 0.0003706863 
    ## ... Similar to previous best
    ## Run 160 stress 0.001220979 
    ## Run 161 stress 0.001365967 
    ## Run 162 stress 0.0005256166 
    ## ... Procrustes: rmse 0.02770105  max resid 0.03843001 
    ## Run 163 stress 9.418898e-05 
    ## ... Procrustes: rmse 0.0001863617  max resid 0.0003418612 
    ## ... Similar to previous best
    ## Run 164 stress 0.2520602 
    ## Run 165 stress 0.001313073 
    ## Run 166 stress 0.0003444551 
    ## ... Procrustes: rmse 0.0134363  max resid 0.01854007 
    ## Run 167 stress 0.000482332 
    ## ... Procrustes: rmse 0.0160372  max resid 0.02212178 
    ## Run 168 stress 0.0004687971 
    ## ... Procrustes: rmse 0.0261557  max resid 0.03628864 
    ## Run 169 stress 0.001439863 
    ## Run 170 stress 0.001191124 
    ## Run 171 stress 0.001409541 
    ## Run 172 stress 7.979567e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0006520196  max resid 0.0009097473 
    ## ... Similar to previous best
    ## Run 173 stress 8.974416e-05 
    ## ... Procrustes: rmse 0.001914717  max resid 0.002640862 
    ## ... Similar to previous best
    ## Run 174 stress 0.001346792 
    ## Run 175 stress 0.0004853575 
    ## ... Procrustes: rmse 0.01543231  max resid 0.02126603 
    ## Run 176 stress 0.00128072 
    ## Run 177 stress 9.912107e-05 
    ## ... Procrustes: rmse 0.0006777025  max resid 0.00102443 
    ## ... Similar to previous best
    ## Run 178 stress 0.001447571 
    ## Run 179 stress 0.2797318 
    ## Run 180 stress 0.0003207418 
    ## ... Procrustes: rmse 0.01208488  max resid 0.01664428 
    ## Run 181 stress 9.688904e-05 
    ## ... Procrustes: rmse 0.0007328267  max resid 0.001186953 
    ## ... Similar to previous best
    ## Run 182 stress 0.001314496 
    ## Run 183 stress 0.0007514229 
    ## Run 184 stress 9.231046e-05 
    ## ... Procrustes: rmse 0.0006645928  max resid 0.001137188 
    ## ... Similar to previous best
    ## Run 185 stress 8.693926e-05 
    ## ... Procrustes: rmse 0.0007320052  max resid 0.00109212 
    ## ... Similar to previous best
    ## Run 186 stress 0.00127944 
    ## Run 187 stress 9.059002e-05 
    ## ... Procrustes: rmse 0.0006738421  max resid 0.001153564 
    ## ... Similar to previous best
    ## Run 188 stress 0.0001423619 
    ## ... Procrustes: rmse 0.008034435  max resid 0.0110578 
    ## Run 189 stress 0.00136944 
    ## Run 190 stress 0.00134236 
    ## Run 191 stress 0.0004596532 
    ## ... Procrustes: rmse 0.01501116  max resid 0.02068504 
    ## Run 192 stress 9.984815e-05 
    ## ... Procrustes: rmse 0.0009094924  max resid 0.001257933 
    ## ... Similar to previous best
    ## Run 193 stress 8.431555e-05 
    ## ... Procrustes: rmse 0.0007098103  max resid 0.001153008 
    ## ... Similar to previous best
    ## Run 194 stress 0.3082926 
    ## Run 195 stress 0.001291225 
    ## Run 196 stress 0.001308483 
    ## Run 197 stress 0.0003871756 
    ## ... Procrustes: rmse 0.01370448  max resid 0.01888218 
    ## Run 198 stress 9.134157e-05 
    ## ... Procrustes: rmse 0.0006774953  max resid 0.001153296 
    ## ... Similar to previous best
    ## Run 199 stress 9.634016e-05 
    ## ... Procrustes: rmse 0.0007374861  max resid 0.001110385 
    ## ... Similar to previous best
    ## Run 200 stress 0.001321799 
    ## Run 201 stress 0.2842805 
    ## Run 202 stress 0.3082927 
    ## Run 203 stress 0.001309395 
    ## Run 204 stress 0.00137557 
    ## Run 205 stress 0.0004044711 
    ## ... Procrustes: rmse 0.01402784  max resid 0.01932877 
    ## Run 206 stress 0.1776025 
    ## Run 207 stress 9.092133e-05 
    ## ... Procrustes: rmse 0.0006958014  max resid 0.001003599 
    ## ... Similar to previous best
    ## Run 208 stress 0.001338917 
    ## Run 209 stress 8.188347e-05 
    ## ... Procrustes: rmse 0.0009528274  max resid 0.00216763 
    ## ... Similar to previous best
    ## Run 210 stress 8.828025e-05 
    ## ... Procrustes: rmse 0.0007200747  max resid 0.0009910216 
    ## ... Similar to previous best
    ## Run 211 stress 9.137231e-05 
    ## ... Procrustes: rmse 0.00112859  max resid 0.002535414 
    ## ... Similar to previous best
    ## Run 212 stress 8.694417e-05 
    ## ... Procrustes: rmse 0.0007047797  max resid 0.001022308 
    ## ... Similar to previous best
    ## Run 213 stress 9.59869e-05 
    ## ... Procrustes: rmse 0.0007328464  max resid 0.001181461 
    ## ... Similar to previous best
    ## Run 214 stress 7.619911e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0006819805  max resid 0.001115357 
    ## ... Similar to previous best
    ## Run 215 stress 0.1776019 
    ## Run 216 stress 0.1491957 
    ## Run 217 stress 9.213349e-05 
    ## ... Procrustes: rmse 0.0001815518  max resid 0.0003182076 
    ## ... Similar to previous best
    ## Run 218 stress 0.00141318 
    ## Run 219 stress 0.001034503 
    ## Run 220 stress 9.375328e-05 
    ## ... Procrustes: rmse 0.0001914685  max resid 0.0003781146 
    ## ... Similar to previous best
    ## Run 221 stress 9.952148e-05 
    ## ... Procrustes: rmse 4.648282e-05  max resid 9.155485e-05 
    ## ... Similar to previous best
    ## Run 222 stress 0.2696744 
    ## Run 223 stress 9.319355e-05 
    ## ... Procrustes: rmse 0.0001080631  max resid 0.0002327791 
    ## ... Similar to previous best
    ## Run 224 stress 9.897863e-05 
    ## ... Procrustes: rmse 0.0002100702  max resid 0.0004505147 
    ## ... Similar to previous best
    ## Run 225 stress 0.1491957 
    ## Run 226 stress 0.001382546 
    ## Run 227 stress 9.997726e-05 
    ## ... Procrustes: rmse 0.007277131  max resid 0.009871878 
    ## Run 228 stress 0.001176866 
    ## Run 229 stress 0.1491957 
    ## Run 230 stress 0.0005004854 
    ## ... Procrustes: rmse 0.01636083  max resid 0.02240688 
    ## Run 231 stress 0.001382457 
    ## Run 232 stress 0.2509018 
    ## Run 233 stress 0.001477539 
    ## Run 234 stress 0.001369426 
    ## Run 235 stress 9.516491e-05 
    ## ... Procrustes: rmse 9.417432e-05  max resid 0.0002132143 
    ## ... Similar to previous best
    ## Run 236 stress 9.910812e-05 
    ## ... Procrustes: rmse 4.943099e-05  max resid 8.906915e-05 
    ## ... Similar to previous best
    ## Run 237 stress 0.1990774 
    ## Run 238 stress 9.031695e-05 
    ## ... Procrustes: rmse 2.864682e-05  max resid 4.741612e-05 
    ## ... Similar to previous best
    ## Run 239 stress 9.213019e-05 
    ## ... Procrustes: rmse 4.635448e-05  max resid 8.176758e-05 
    ## ... Similar to previous best
    ## Run 240 stress 9.405278e-05 
    ## ... Procrustes: rmse 0.004705743  max resid 0.00628255 
    ## ... Similar to previous best
    ## Run 241 stress 9.58824e-05 
    ## ... Procrustes: rmse 0.0001250331  max resid 0.0002091723 
    ## ... Similar to previous best
    ## Run 242 stress 9.421942e-05 
    ## ... Procrustes: rmse 4.072234e-05  max resid 6.47985e-05 
    ## ... Similar to previous best
    ## Run 243 stress 0.1491957 
    ## Run 244 stress 0.2528294 
    ## Run 245 stress 9.783078e-05 
    ## ... Procrustes: rmse 0.0001849282  max resid 0.0003493953 
    ## ... Similar to previous best
    ## Run 246 stress 0.0004797253 
    ## ... Procrustes: rmse 0.0160062  max resid 0.02191205 
    ## Run 247 stress 0.2842805 
    ## Run 248 stress 9.580188e-05 
    ## ... Procrustes: rmse 0.0001842922  max resid 0.0003560483 
    ## ... Similar to previous best
    ## Run 249 stress 9.241856e-05 
    ## ... Procrustes: rmse 0.002452068  max resid 0.003381659 
    ## ... Similar to previous best
    ## Run 250 stress 8.42586e-05 
    ## ... Procrustes: rmse 0.005327688  max resid 0.007293481 
    ## Run 251 stress 0.001213122 
    ## Run 252 stress 0.0003981328 
    ## ... Procrustes: rmse 0.01458442  max resid 0.0199559 
    ## Run 253 stress 9.651619e-05 
    ## ... Procrustes: rmse 2.57282e-05  max resid 4.515206e-05 
    ## ... Similar to previous best
    ## Run 254 stress 0.1990774 
    ## Run 255 stress 0.001257128 
    ## Run 256 stress 0.00116018 
    ## Run 257 stress 0.0004688754 
    ## ... Procrustes: rmse 0.01583385  max resid 0.02167992 
    ## Run 258 stress 0.001302813 
    ## Run 259 stress 9.837262e-05 
    ## ... Procrustes: rmse 0.0001143463  max resid 0.0001675942 
    ## ... Similar to previous best
    ## Run 260 stress 0.001224563 
    ## Run 261 stress 9.32965e-05 
    ## ... Procrustes: rmse 0.0001231661  max resid 0.0002041029 
    ## ... Similar to previous best
    ## Run 262 stress 0.0002964478 
    ## ... Procrustes: rmse 0.02075056  max resid 0.02851023 
    ## Run 263 stress 0.0004647721 
    ## ... Procrustes: rmse 0.01574299  max resid 0.02155465 
    ## Run 264 stress 0.0004391748 
    ## ... Procrustes: rmse 0.01530404  max resid 0.02094886 
    ## Run 265 stress 9.646635e-05 
    ## ... Procrustes: rmse 9.704024e-05  max resid 0.0002186666 
    ## ... Similar to previous best
    ## Run 266 stress 0.0009479189 
    ## Run 267 stress 0.2586411 
    ## Run 268 stress 0.001354184 
    ## Run 269 stress 0.001514104 
    ## Run 270 stress 0.001402151 
    ## Run 271 stress 0.3050283 
    ## Run 272 stress 8.804697e-05 
    ## ... Procrustes: rmse 0.0001152785  max resid 0.0001898729 
    ## ... Similar to previous best
    ## Run 273 stress 9.542538e-05 
    ## ... Procrustes: rmse 0.0001859032  max resid 0.0003724436 
    ## ... Similar to previous best
    ## Run 274 stress 0.00132986 
    ## Run 275 stress 0.3120129 
    ## Run 276 stress 0.1491957 
    ## Run 277 stress 9.138916e-05 
    ## ... Procrustes: rmse 0.0001820141  max resid 0.0003762229 
    ## ... Similar to previous best
    ## Run 278 stress 0.0005067143 
    ## ... Procrustes: rmse 0.01645947  max resid 0.02254333 
    ## Run 279 stress 0.001315871 
    ## Run 280 stress 8.478815e-05 
    ## ... Procrustes: rmse 0.0001367922  max resid 0.0002306245 
    ## ... Similar to previous best
    ## Run 281 stress 0.001363233 
    ## Run 282 stress 0.2829397 
    ## Run 283 stress 9.488501e-05 
    ## ... Procrustes: rmse 0.0001953908  max resid 0.000399543 
    ## ... Similar to previous best
    ## Run 284 stress 8.872702e-05 
    ## ... Procrustes: rmse 2.780779e-05  max resid 4.560437e-05 
    ## ... Similar to previous best
    ## Run 285 stress 0.1491957 
    ## Run 286 stress 0.0002485436 
    ## ... Procrustes: rmse 0.01151381  max resid 0.01571881 
    ## Run 287 stress 0.001308398 
    ## Run 288 stress 0.001114139 
    ## Run 289 stress 0.0004326178 
    ## ... Procrustes: rmse 0.02507989  max resid 0.03450332 
    ## Run 290 stress 0.1491957 
    ## Run 291 stress 9.307998e-05 
    ## ... Procrustes: rmse 0.000187254  max resid 0.0003717978 
    ## ... Similar to previous best
    ## Run 292 stress 0.0004942853 
    ## ... Procrustes: rmse 0.0268133  max resid 0.0369049 
    ## Run 293 stress 0.001311145 
    ## Run 294 stress 0.001214523 
    ## Run 295 stress 9.189996e-05 
    ## ... Procrustes: rmse 0.0001204926  max resid 0.0001990526 
    ## ... Similar to previous best
    ## Run 296 stress 0.1776012 
    ## Run 297 stress 0.001175935 
    ## Run 298 stress 5.856329e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001500125  max resid 0.0003048229 
    ## ... Similar to previous best
    ## Run 299 stress 0.001256609 
    ## Run 300 stress 0.1491957 
    ## Run 301 stress 9.306511e-05 
    ## ... Procrustes: rmse 0.0001782116  max resid 0.0003469964 
    ## ... Similar to previous best
    ## Run 302 stress 0.001239344 
    ## Run 303 stress 8.704138e-05 
    ## ... Procrustes: rmse 0.0001686759  max resid 0.0003243864 
    ## ... Similar to previous best
    ## Run 304 stress 7.346394e-05 
    ## ... Procrustes: rmse 0.0007734247  max resid 0.001061142 
    ## ... Similar to previous best
    ## Run 305 stress 9.520105e-05 
    ## ... Procrustes: rmse 0.0001879038  max resid 0.00034334 
    ## ... Similar to previous best
    ## Run 306 stress 0.1990774 
    ## Run 307 stress 0.0004930905 
    ## ... Procrustes: rmse 0.01622726  max resid 0.02233451 
    ## Run 308 stress 0.0002629958 
    ## ... Procrustes: rmse 0.01183101  max resid 0.01626825 
    ## Run 309 stress 0.0004042182 
    ## ... Procrustes: rmse 0.01468381  max resid 0.02020499 
    ## Run 310 stress 0.0002968784 
    ## ... Procrustes: rmse 0.02082385  max resid 0.02862547 
    ## Run 311 stress 0.001197036 
    ## Run 312 stress 9.759771e-05 
    ## ... Procrustes: rmse 0.0001298207  max resid 0.0001663395 
    ## ... Similar to previous best
    ## Run 313 stress 9.44018e-05 
    ## ... Procrustes: rmse 0.0001770932  max resid 0.0003257295 
    ## ... Similar to previous best
    ## Run 314 stress 0.001235086 
    ## Run 315 stress 0.001243417 
    ## Run 316 stress 0.001222078 
    ## Run 317 stress 8.273769e-05 
    ## ... Procrustes: rmse 0.0006519903  max resid 0.0009736598 
    ## ... Similar to previous best
    ## Run 318 stress 0.1990774 
    ## Run 319 stress 0.2854395 
    ## Run 320 stress 9.725584e-05 
    ## ... Procrustes: rmse 0.0001234188  max resid 0.0001613404 
    ## ... Similar to previous best
    ## Run 321 stress 9.359473e-05 
    ## ... Procrustes: rmse 0.0001791322  max resid 0.0003484407 
    ## ... Similar to previous best
    ## Run 322 stress 0.1491957 
    ## Run 323 stress 0.001322202 
    ## Run 324 stress 0.001188165 
    ## Run 325 stress 9.332177e-05 
    ## ... Procrustes: rmse 0.0001238113  max resid 0.0002063362 
    ## ... Similar to previous best
    ## Run 326 stress 0.0001097719 
    ## ... Procrustes: rmse 0.007623143  max resid 0.01046138 
    ## Run 327 stress 0.1491957 
    ## Run 328 stress 0.001389922 
    ## Run 329 stress 9.333779e-05 
    ## ... Procrustes: rmse 0.0002028106  max resid 0.0003728013 
    ## ... Similar to previous best
    ## Run 330 stress 0.2612448 
    ## Run 331 stress 9.952199e-05 
    ## ... Procrustes: rmse 0.0001324032  max resid 0.0001724107 
    ## ... Similar to previous best
    ## Run 332 stress 0.001319208 
    ## Run 333 stress 0.0003522363 
    ## ... Procrustes: rmse 0.01367645  max resid 0.01881528 
    ## Run 334 stress 0.0005052968 
    ## ... Procrustes: rmse 0.01642783  max resid 0.02261162 
    ## Run 335 stress 0.1491957 
    ## Run 336 stress 9.177255e-05 
    ## ... Procrustes: rmse 0.000115333  max resid 0.0001868724 
    ## ... Similar to previous best
    ## Run 337 stress 0.0004924918 
    ## ... Procrustes: rmse 0.01620889  max resid 0.02230876 
    ## Run 338 stress 9.636329e-05 
    ## ... Procrustes: rmse 0.01001424  max resid 0.01371565 
    ## Run 339 stress 8.52068e-05 
    ## ... Procrustes: rmse 0.0001466771  max resid 0.0002073485 
    ## ... Similar to previous best
    ## Run 340 stress 0.3083098 
    ## Run 341 stress 9.993922e-05 
    ## ... Procrustes: rmse 0.007270651  max resid 0.009974957 
    ## Run 342 stress 0.001129306 
    ## Run 343 stress 0.001362026 
    ## Run 344 stress 0.001272731 
    ## Run 345 stress 7.12823e-05 
    ## ... Procrustes: rmse 0.0001312037  max resid 0.0002630967 
    ## ... Similar to previous best
    ## Run 346 stress 0.1491957 
    ## Run 347 stress 7.623825e-05 
    ## ... Procrustes: rmse 0.003555091  max resid 0.004779922 
    ## ... Similar to previous best
    ## Run 348 stress 0.3083098 
    ## Run 349 stress 9.439311e-05 
    ## ... Procrustes: rmse 0.0001768941  max resid 0.0003266691 
    ## ... Similar to previous best
    ## Run 350 stress 9.508158e-05 
    ## ... Procrustes: rmse 0.0002097558  max resid 0.0003606692 
    ## ... Similar to previous best
    ## Run 351 stress 9.693154e-05 
    ## ... Procrustes: rmse 0.0002123699  max resid 0.0003646416 
    ## ... Similar to previous best
    ## Run 352 stress 0.1776015 
    ## Run 353 stress 0.1776024 
    ## Run 354 stress 0.1491957 
    ## Run 355 stress 9.632199e-05 
    ## ... Procrustes: rmse 0.0002092214  max resid 0.0003697739 
    ## ... Similar to previous best
    ## Run 356 stress 9.782872e-05 
    ## ... Procrustes: rmse 0.0002112096  max resid 0.000378776 
    ## ... Similar to previous best
    ## Run 357 stress 0.0004558674 
    ## ... Procrustes: rmse 0.0154764  max resid 0.02129822 
    ## Run 358 stress 0.001275796 
    ## Run 359 stress 0.001243855 
    ## Run 360 stress 0.001284366 
    ## Run 361 stress 8.3216e-05 
    ## ... Procrustes: rmse 0.000107224  max resid 0.0001868963 
    ## ... Similar to previous best
    ## Run 362 stress 0.001193016 
    ## Run 363 stress 0.001321975 
    ## Run 364 stress 0.001403294 
    ## Run 365 stress 9.756867e-05 
    ## ... Procrustes: rmse 0.0001801624  max resid 0.0003589845 
    ## ... Similar to previous best
    ## Run 366 stress 0.001247175 
    ## Run 367 stress 8.593004e-05 
    ## ... Procrustes: rmse 0.0001546092  max resid 0.0002165388 
    ## ... Similar to previous best
    ## Run 368 stress 0.001219104 
    ## Run 369 stress 0.001439453 
    ## Run 370 stress 9.458846e-05 
    ## ... Procrustes: rmse 0.0002074514  max resid 0.0003659274 
    ## ... Similar to previous best
    ## Run 371 stress 0.0004622406 
    ## ... Procrustes: rmse 0.01568794  max resid 0.02159255 
    ## Run 372 stress 9.851048e-05 
    ## ... Procrustes: rmse 0.0002206642  max resid 0.0003781883 
    ## ... Similar to previous best
    ## Run 373 stress 0.000639291 
    ## Run 374 stress 8.944636e-05 
    ## ... Procrustes: rmse 0.0001064748  max resid 0.0001665901 
    ## ... Similar to previous best
    ## Run 375 stress 0.1491957 
    ## Run 376 stress 9.353215e-05 
    ## ... Procrustes: rmse 0.0001752226  max resid 0.0003577953 
    ## ... Similar to previous best
    ## Run 377 stress 0.001248578 
    ## Run 378 stress 0.0004899953 
    ## ... Procrustes: rmse 0.01616405  max resid 0.0222498 
    ## Run 379 stress 9.49953e-05 
    ## ... Procrustes: rmse 0.0002075186  max resid 0.00035851 
    ## ... Similar to previous best
    ## Run 380 stress 0.3082762 
    ## Run 381 stress 9.384504e-05 
    ## ... Procrustes: rmse 0.0001166918  max resid 0.0001494155 
    ## ... Similar to previous best
    ## Run 382 stress 9.70644e-05 
    ## ... Procrustes: rmse 0.0001832631  max resid 0.0003516657 
    ## ... Similar to previous best
    ## Run 383 stress 9.956042e-05 
    ## ... Procrustes: rmse 0.0001349361  max resid 0.000203727 
    ## ... Similar to previous best
    ## Run 384 stress 0.1776011 
    ## Run 385 stress 9.304153e-05 
    ## ... Procrustes: rmse 0.0001237576  max resid 0.0002195076 
    ## ... Similar to previous best
    ## Run 386 stress 0.1776022 
    ## Run 387 stress 9.477605e-05 
    ## ... Procrustes: rmse 0.0001192576  max resid 0.0002001096 
    ## ... Similar to previous best
    ## Run 388 stress 0.2842805 
    ## Run 389 stress 0.2848514 
    ## Run 390 stress 0.001195486 
    ## Run 391 stress 0.001150185 
    ## Run 392 stress 0.2842805 
    ## Run 393 stress 0.001319851 
    ## Run 394 stress 0.0007703631 
    ## Run 395 stress 0.0004817361 
    ## ... Procrustes: rmse 0.01603974  max resid 0.02207632 
    ## Run 396 stress 9.539574e-05 
    ## ... Procrustes: rmse 0.0001767432  max resid 0.0003358184 
    ## ... Similar to previous best
    ## Run 397 stress 9.57553e-05 
    ## ... Procrustes: rmse 0.0001203735  max resid 0.0001625666 
    ## ... Similar to previous best
    ## Run 398 stress 9.584995e-05 
    ## ... Procrustes: rmse 0.000213992  max resid 0.000368283 
    ## ... Similar to previous best
    ## Run 399 stress 9.36532e-05 
    ## ... Procrustes: rmse 0.0001651682  max resid 0.0002440175 
    ## ... Similar to previous best
    ## Run 400 stress 0.001314418 
    ## Run 401 stress 0.001227006 
    ## Run 402 stress 0.0004192002 
    ## ... Procrustes: rmse 0.01495886  max resid 0.0205841 
    ## Run 403 stress 9.646685e-05 
    ## ... Procrustes: rmse 0.0001619538  max resid 0.0002496595 
    ## ... Similar to previous best
    ## Run 404 stress 9.353216e-05 
    ## ... Procrustes: rmse 0.0001271106  max resid 0.0002056726 
    ## ... Similar to previous best
    ## Run 405 stress 9.537872e-05 
    ## ... Procrustes: rmse 0.000163955  max resid 0.0002555666 
    ## ... Similar to previous best
    ## Run 406 stress 7.123051e-05 
    ## ... Procrustes: rmse 0.0001661261  max resid 0.0003132278 
    ## ... Similar to previous best
    ## Run 407 stress 9.726422e-05 
    ## ... Procrustes: rmse 0.0001577325  max resid 0.0002192786 
    ## ... Similar to previous best
    ## Run 408 stress 0.1491957 
    ## Run 409 stress 9.3907e-05 
    ## ... Procrustes: rmse 0.0001053305  max resid 0.0001591285 
    ## ... Similar to previous best
    ## Run 410 stress 0.0004677953 
    ## ... Procrustes: rmse 0.0153397  max resid 0.02111089 
    ## Run 411 stress 0.3083099 
    ## Run 412 stress 0.1990774 
    ## Run 413 stress 0.1491957 
    ## Run 414 stress 0.00109926 
    ## Run 415 stress 0.001399536 
    ## Run 416 stress 9.604642e-05 
    ## ... Procrustes: rmse 0.0001251045  max resid 0.0002269503 
    ## ... Similar to previous best
    ## Run 417 stress 0.001297975 
    ## Run 418 stress 7.510085e-05 
    ## ... Procrustes: rmse 0.0003016029  max resid 0.0004550286 
    ## ... Similar to previous best
    ## Run 419 stress 9.999867e-05 
    ## ... Procrustes: rmse 0.007272688  max resid 0.009977628 
    ## Run 420 stress 0.1990774 
    ## Run 421 stress 0.001186656 
    ## Run 422 stress 0.0004964293 
    ## ... Procrustes: rmse 0.01628421  max resid 0.02241288 
    ## Run 423 stress 9.63704e-05 
    ## ... Procrustes: rmse 0.0001683726  max resid 0.0003432343 
    ## ... Similar to previous best
    ## Run 424 stress 0.1491957 
    ## Run 425 stress 8.558302e-05 
    ## ... Procrustes: rmse 0.0001064562  max resid 0.0001563985 
    ## ... Similar to previous best
    ## Run 426 stress 0.0004764137 
    ## ... Procrustes: rmse 0.02638065  max resid 0.03631846 
    ## Run 427 stress 0.2852131 
    ## Run 428 stress 0.000405009 
    ## ... Procrustes: rmse 0.01469862  max resid 0.02022541 
    ## Run 429 stress 0.001423024 
    ## Run 430 stress 4.0983e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002782659  max resid 0.0004372145 
    ## ... Similar to previous best
    ## Run 431 stress 0.001252478 
    ## Run 432 stress 8.211268e-05 
    ## ... Procrustes: rmse 0.0001742138  max resid 0.0003284794 
    ## ... Similar to previous best
    ## Run 433 stress 0.0005015505 
    ## ... Procrustes: rmse 0.01637264  max resid 0.02309258 
    ## Run 434 stress 0.00121625 
    ## Run 435 stress 9.727861e-05 
    ## ... Procrustes: rmse 0.0002523119  max resid 0.0005610801 
    ## ... Similar to previous best
    ## Run 436 stress 0.1776022 
    ## Run 437 stress 9.938692e-05 
    ## ... Procrustes: rmse 0.0002465987  max resid 0.0005250412 
    ## ... Similar to previous best
    ## Run 438 stress 0.001381438 
    ## Run 439 stress 0.1990774 
    ## Run 440 stress 9.935568e-05 
    ## ... Procrustes: rmse 0.0003006066  max resid 0.0005637518 
    ## ... Similar to previous best
    ## Run 441 stress 0.001334357 
    ## Run 442 stress 9.465182e-05 
    ## ... Procrustes: rmse 0.0002148322  max resid 0.0003524767 
    ## ... Similar to previous best
    ## Run 443 stress 9.278456e-05 
    ## ... Procrustes: rmse 0.0001918716  max resid 0.0003182188 
    ## ... Similar to previous best
    ## Run 444 stress 0.1491957 
    ## Run 445 stress 0.0004853081 
    ## ... Procrustes: rmse 0.01605576  max resid 0.02265599 
    ## Run 446 stress 0.2848523 
    ## Run 447 stress 9.909314e-05 
    ## ... Procrustes: rmse 0.0002061007  max resid 0.0004158302 
    ## ... Similar to previous best
    ## Run 448 stress 0.001219403 
    ## Run 449 stress 9.0242e-05 
    ## ... Procrustes: rmse 0.0002919722  max resid 0.0005433136 
    ## ... Similar to previous best
    ## Run 450 stress 0.001208487 
    ## Run 451 stress 0.2373687 
    ## Run 452 stress 0.001321721 
    ## Run 453 stress 9.110829e-05 
    ## ... Procrustes: rmse 0.000239646  max resid 0.0005168559 
    ## ... Similar to previous best
    ## Run 454 stress 9.8499e-05 
    ## ... Procrustes: rmse 0.0002552066  max resid 0.000564622 
    ## ... Similar to previous best
    ## Run 455 stress 0.001239004 
    ## Run 456 stress 8.279455e-05 
    ## ... Procrustes: rmse 0.0002172901  max resid 0.0003824671 
    ## ... Similar to previous best
    ## Run 457 stress 0.1990774 
    ## Run 458 stress 9.635017e-05 
    ## ... Procrustes: rmse 0.0002522145  max resid 0.0005589217 
    ## ... Similar to previous best
    ## Run 459 stress 0.001207271 
    ## Run 460 stress 0.001181715 
    ## Run 461 stress 0.0002338257 
    ## ... Procrustes: rmse 0.01116857  max resid 0.01591114 
    ## Run 462 stress 0.001357404 
    ## Run 463 stress 9.724656e-05 
    ## ... Procrustes: rmse 0.003310706  max resid 0.005078068 
    ## ... Similar to previous best
    ## Run 464 stress 0.000474433 
    ## ... Procrustes: rmse 0.01592979  max resid 0.02248151 
    ## Run 465 stress 9.229989e-05 
    ## ... Procrustes: rmse 0.0002102746  max resid 0.0004139185 
    ## ... Similar to previous best
    ## Run 466 stress 0.0004692837 
    ## ... Procrustes: rmse 0.01584046  max resid 0.02235867 
    ## Run 467 stress 0.1491957 
    ## Run 468 stress 9.373114e-05 
    ## ... Procrustes: rmse 0.0001715548  max resid 0.0003344395 
    ## ... Similar to previous best
    ## Run 469 stress 0.001373052 
    ## Run 470 stress 0.001276687 
    ## Run 471 stress 0.1491957 
    ## Run 472 stress 0.1990774 
    ## Run 473 stress 0.001312219 
    ## Run 474 stress 9.295212e-05 
    ## ... Procrustes: rmse 0.0002101025  max resid 0.0004280618 
    ## ... Similar to previous best
    ## Run 475 stress 9.42855e-05 
    ## ... Procrustes: rmse 0.0001701723  max resid 0.0003281466 
    ## ... Similar to previous best
    ## Run 476 stress 8.986719e-05 
    ## ... Procrustes: rmse 0.0003089571  max resid 0.0005276593 
    ## ... Similar to previous best
    ## Run 477 stress 0.001295974 
    ## Run 478 stress 0.001197521 
    ## Run 479 stress 9.437453e-05 
    ## ... Procrustes: rmse 0.0002699742  max resid 0.0005108436 
    ## ... Similar to previous best
    ## Run 480 stress 0.001223597 
    ## Run 481 stress 0.1776019 
    ## Run 482 stress 8.808349e-05 
    ## ... Procrustes: rmse 0.0002110884  max resid 0.0003308886 
    ## ... Similar to previous best
    ## Run 483 stress 0.3083098 
    ## Run 484 stress 0.0004763399 
    ## ... Procrustes: rmse 0.01582303  max resid 0.02233511 
    ## Run 485 stress 0.001419761 
    ## Run 486 stress 0.1491957 
    ## Run 487 stress 0.001039397 
    ## Run 488 stress 0.001273323 
    ## Run 489 stress 0.1491957 
    ## Run 490 stress 0.1990774 
    ## Run 491 stress 9.598285e-05 
    ## ... Procrustes: rmse 0.0002729598  max resid 0.0005213462 
    ## ... Similar to previous best
    ## Run 492 stress 0.001111219 
    ## Run 493 stress 0.000477352 
    ## ... Procrustes: rmse 0.01596875  max resid 0.02254903 
    ## Run 494 stress 0.001355931 
    ## Run 495 stress 0.001309869 
    ## Run 496 stress 0.0010713 
    ## Run 497 stress 8.850303e-05 
    ## ... Procrustes: rmse 0.0002706094  max resid 0.0005086188 
    ## ... Similar to previous best
    ## Run 498 stress 0.0001002548 
    ## ... Procrustes: rmse 0.01183075  max resid 0.01629487 
    ## Run 499 stress 0.1990774 
    ## Run 500 stress 8.581403e-05 
    ## ... Procrustes: rmse 0.0001847303  max resid 0.0003223908 
    ## ... Similar to previous best
    ## *** Best solution repeated 24 times

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
    ## env[surveyed_sites_env, c(34)]  0.997840 -0.065695 0.7388   0.02 *
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
    ## env[surveyed_sites_env, c(34)]  0.997840 -0.065695 0.7388   0.02 *
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
    ## temperature_median -0.48765 -0.87304 0.0586  0.825  
    ## salinity_median     0.99784 -0.06569 0.7388  0.034 *
    ## oxygen_median       0.79309  0.60911 0.5361  0.816  
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
    ## temperature_median -0.48765 -0.87304 0.0586  1.000
    ## salinity_median     0.99784 -0.06569 0.7388  0.102
    ## oxygen_median       0.79309  0.60911 0.5361  1.000
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
    ## temperature_median -0.56551 -0.82474 0.0778  0.795  
    ## salinity_median     0.99174  0.12829 0.7458  0.036 *
    ## oxygen_median       0.94280  0.33337 0.3552  0.886  
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
    ## temperature_median -0.56551 -0.82474 0.0778  1.000
    ## salinity_median     0.99174  0.12829 0.7458  0.108
    ## oxygen_median       0.94280  0.33337 0.3552  1.000
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
    ## temperature_median -0.10614  0.99435 0.0705  0.704   
    ## salinity_median    -0.78732 -0.61654 0.7842  0.010 **
    ## oxygen_median      -0.96711  0.25435 0.7188  0.062 . 
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
    ## salinity_median    -0.78732 -0.61654 0.7842  0.030 *
    ## oxygen_median      -0.96711  0.25435 0.7188  0.186  
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
    ## temperature_median  0.0003827 -1.0000000 0.0340  0.844
    ## salinity_median    -0.0027304 -1.0000000 0.4836  0.888
    ## oxygen_median      -0.0108511 -0.9999400 0.7114  0.973
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
    ## temperature_median  0.0003827 -1.0000000 0.0340      1
    ## salinity_median    -0.0027304 -1.0000000 0.4836      1
    ## oxygen_median      -0.0108511 -0.9999400 0.7114      1
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
    ## temperature_median -0.0041007  0.9999900 0.0301  0.900  
    ## salinity_median     0.0036814 -0.9999900 0.6142  0.070 .
    ## oxygen_median       0.0130802 -0.9999100 0.6096  0.105  
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
    ## temperature_median -0.0041007  0.9999900 0.0301  1.000
    ## salinity_median     0.0036814 -0.9999900 0.6142  0.210
    ## oxygen_median       0.0130802 -0.9999100 0.6096  0.315
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
    ## temperature_median       -0.35881 -0.93341 0.1964  0.592
    ## salinity_median          -0.70036  0.71379 0.5991  0.124
    ## oxygen_median             0.98543 -0.17006 0.5185  0.157
    ## distance_to_ocean_mean_m  0.01326 -0.99991 0.1726  0.596
    ## max_depth                -0.49998 -0.86604 0.0383  0.892
    ## logArea                  -0.26513 -0.96421 0.2144  0.555
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
    ## salinity_median          -0.70036  0.71379 0.5991  0.744
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
    ## distance_to_ocean_min_m -0.99980  0.02020 0.6139  0.125
    ## max_depth               -0.17719 -0.98418 0.1186  0.777
    ## logArea                  0.26555 -0.96410 0.2082  0.116
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
    ## distance_to_ocean_min_m -0.99980  0.02020 0.6139  0.375
    ## max_depth               -0.17719 -0.98418 0.1186  1.000
    ## logArea                  0.26555 -0.96410 0.2082  0.348
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
    ## distance_to_ocean_min_m -0.99980  0.02020 0.6139  0.142
    ## max_depth               -0.17719 -0.98418 0.1186  0.766
    ## logArea                  0.26555 -0.96410 0.2082  0.117
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
    ## distance_to_ocean_min_m -0.99980  0.02020 0.6139  0.426
    ## max_depth               -0.17719 -0.98418 0.1186  1.000
    ## logArea                  0.26555 -0.96410 0.2082  0.351
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
    ## distance_to_ocean_min_m -0.94032  0.34028 0.5170  0.152
    ## max_depth               -0.21482 -0.97665 0.1851  0.650
    ## logArea                 -0.06538  0.99786 0.0118  0.958
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
    ## distance_to_ocean_min_m -0.94032  0.34028 0.5170  0.456
    ## max_depth               -0.21482 -0.97665 0.1851  1.000
    ## logArea                 -0.06538  0.99786 0.0118  1.000
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
    ## distance_to_ocean_min_m  0.86176 -0.50731 0.3198  0.077 .
    ## max_depth               -0.68924 -0.72453 0.5689  0.018 *
    ## logArea                 -0.84908  0.52827 0.2993  0.240  
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
    ## distance_to_ocean_min_m  0.86176 -0.50731 0.3198  0.231  
    ## max_depth               -0.68924 -0.72453 0.5689  0.054 .
    ## logArea                 -0.84908  0.52827 0.2993  0.720  
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
    ## distance_to_ocean_min_m  0.86962  0.49371 0.6793  0.192
    ## max_depth                0.68534  0.72822 0.0510  0.954
    ## logArea                 -0.52043  0.85391 0.2187  0.545
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
    ## distance_to_ocean_min_m  0.86962  0.49371 0.6793  0.576
    ## max_depth                0.68534  0.72822 0.0510  1.000
    ## logArea                 -0.52043  0.85391 0.2187  1.000
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
    ## distance_to_ocean_min_m -0.00074386  1.00000000 0.6632  0.054 . 
    ## max_depth                0.00185967 -1.00000000 0.8305  0.018 * 
    ## logArea                  0.00153151  1.00000000 0.8246  0.003 **
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
    ## distance_to_ocean_min_m -0.00074386  1.00000000 0.6632  0.162   
    ## max_depth                0.00185967 -1.00000000 0.8305  0.054 . 
    ## logArea                  0.00153151  1.00000000 0.8246  0.009 **
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
    ##       Significance: 0.255 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.293 0.325 0.348 0.375 
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
    ## 0.578 0.612 0.626 0.655 
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
    ##       Significance: 0.51 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.586 0.631 0.672 0.708 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_mant_pv <- rbind(SD_beta_env_mant_t$signif, SD_beta_env_mant_s$signif, SD_beta_env_mant_o$signif)
SD_beta_env_mant_pv <- SD_beta_env_mant_pv[,1]
(SD_beta_env_mant_pv <- p.adjust(SD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.765 0.018 1.000

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
    ##       Significance: 0.352 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.340 0.369 0.398 0.426 
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
    ##       Significance: 0.015 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.548 0.588 0.637 0.675 
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
    ##       Significance: 0.606 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.493 0.540 0.565 0.593 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_MS_mant_pv <- rbind(SD_beta_env_MS_mant_t$signif, SD_beta_env_MS_mant_s$signif, SD_beta_env_MS_mant_o$signif)
SD_beta_env_MS_mant_pv <- SD_beta_env_MS_mant_pv[,1]
(SD_beta_env_MS_mant_pv <- p.adjust(SD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.045 1.000

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
    ##       Significance: 0.75 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.123 0.161 0.209 0.405 
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
    ##       Significance: 0.046 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.349 0.413 0.510 0.560 
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
    ##       Significance: 0.024 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.315 0.430 0.523 0.615 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_OM_mant_pv <- rbind(SD_beta_env_OM_mant_t$signif, SD_beta_env_OM_mant_s$signif, SD_beta_env_OM_mant_o$signif)
SD_beta_env_OM_mant_pv <- SD_beta_env_OM_mant_pv[,1]
(SD_beta_env_OM_mant_pv <- p.adjust(SD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.138 0.072

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
    ##       Significance: 0.228 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0270 0.0489 0.0610 0.0802 
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
    ##       Significance: 0.014 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.459 0.504 0.535 0.555 
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
    ##       Significance: 0.369 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.710 0.734 0.747 0.757 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_SO_mant_pv <- rbind(SD_beta_env_SO_mant_t$signif, SD_beta_env_SO_mant_s$signif, SD_beta_env_SO_mant_o$signif)
SD_beta_env_SO_mant_pv <- SD_beta_env_SO_mant_pv[,1]
(SD_beta_env_SO_mant_pv <- p.adjust(SD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 0.684 0.042 1.000

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
    ##       Significance: 0.757 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.260 0.354 0.459 0.697 
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
    ##       Significance: 0.162 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.281 0.345 0.380 0.443 
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
    ## 0.238 0.389 0.476 0.559 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_M_mant_pv <- rbind(SD_beta_env_M_mant_t$signif, SD_beta_env_M_mant_s$signif, SD_beta_env_M_mant_o$signif)
SD_beta_env_M_mant_pv <- SD_beta_env_M_mant_pv[,1]
(SD_beta_env_M_mant_pv <- p.adjust(SD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.486 0.252

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
    ##       Significance: 0.253 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.220 0.288 0.333 0.371 
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
    ##       Significance: 0.458 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.538 0.564 0.582 0.615 
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
    ##       Significance: 0.635 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.201 0.228 0.256 0.280 
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
    ##       Significance: 0.39 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0801 0.0960 0.1077 0.1200 
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
    ##       Significance: 0.378 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.447 0.491 0.538 0.576 
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
    ##       Significance: 0.737 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.223 0.263 0.293 0.320 
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
    ##       Significance: 0.752 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0765 0.0999 0.1146 0.1290 
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
    ##       Significance: 0.588 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.229 0.267 0.292 0.328 
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
    ##       Significance: 0.054 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.177 0.238 0.294 0.347 
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
    ##       Significance: 0.103 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.203 0.237 0.282 0.322 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_OM_mant_pv <- rbind( SD_beta_geo_OM_mant_dmean$signif, SD_beta_geo_OM_mant_md$signif, SD_beta_geo_OM_mant_la$signif)
SD_beta_geo_OM_mant_pv <- SD_beta_geo_OM_mant_pv[,1]
(SD_beta_geo_OM_mant_pv <- p.adjust(SD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.162 0.309

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
    ##       Significance: 0.458 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.674 0.705 0.725 0.746 
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
    ##       Significance: 0.614 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.139 0.169 0.195 0.238 
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
    ##       Significance: 0.208 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.275 0.298 0.314 0.323 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_SO_mant_pv <- rbind( SD_beta_geo_SO_mant_dmean$signif, SD_beta_geo_SO_mant_md$signif, SD_beta_geo_SO_mant_la$signif)
SD_beta_geo_SO_mant_pv <- SD_beta_geo_SO_mant_pv[,1]
(SD_beta_geo_SO_mant_pv <- p.adjust(SD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 1.000 0.624

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
    ##       Significance: 0.452 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.203 0.376 0.471 0.765 
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
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.274 0.420 0.513 0.582 
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
    ##       Significance: 0.029 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.246 0.322 0.438 0.532 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_M_mant_pv <- rbind( SD_beta_geo_M_mant_dmean$signif, SD_beta_geo_M_mant_md$signif, SD_beta_geo_M_mant_la$signif)
SD_beta_geo_M_mant_pv <- SD_beta_geo_M_mant_pv[,1]
(SD_beta_geo_M_mant_pv <- p.adjust(SD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.003 0.087

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
