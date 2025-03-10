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
    ## 2     Stratified vs Ocean  1 1.3093089 3.913330 0.2459152   0.002      0.012
    ## 3 Stratified vs Reference  1 0.6261283 1.909357 0.2143091   0.100      0.600
    ## 4          Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.041      0.246
    ## 5      Mixed vs Reference  1 0.5979203 1.966332 0.2193017   0.128      0.768
    ## 6      Ocean vs Reference  1 0.5599217 1.628213 0.2456489   0.131      0.786
    ##   sig
    ## 1   *
    ## 2   .
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
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.046      0.138

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
    ## 1 Stratified vs Mixed  1 1.3139866 4.158150 0.2289964   0.002      0.006   *
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
    ## 1 Stratified vs Mixed  1 1.4219030 4.731564 0.2827927   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.3260066 4.147587 0.2931657   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.050      0.150

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
    ## Run 1 stress 0.1050338 
    ## ... Procrustes: rmse 0.008873622  max resid 0.02482444 
    ## Run 2 stress 0.1058713 
    ## Run 3 stress 0.1058707 
    ## Run 4 stress 0.1083501 
    ## Run 5 stress 0.1051923 
    ## ... Procrustes: rmse 0.02214318  max resid 0.07895199 
    ## Run 6 stress 0.1138984 
    ## Run 7 stress 0.1082722 
    ## Run 8 stress 0.1063128 
    ## Run 9 stress 0.1084965 
    ## Run 10 stress 0.1050144 
    ## ... Procrustes: rmse 0.006661032  max resid 0.02376114 
    ## Run 11 stress 0.1165305 
    ## Run 12 stress 0.1049791 
    ## ... New best solution
    ## ... Procrustes: rmse 0.006693877  max resid 0.02559062 
    ## Run 13 stress 0.1309381 
    ## Run 14 stress 0.113594 
    ## Run 15 stress 0.1084964 
    ## Run 16 stress 0.1079382 
    ## Run 17 stress 0.1050338 
    ## ... Procrustes: rmse 0.006255406  max resid 0.02186054 
    ## Run 18 stress 0.1063128 
    ## Run 19 stress 0.1051924 
    ## ... Procrustes: rmse 0.02283569  max resid 0.07726194 
    ## Run 20 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283886  max resid 0.07730991 
    ## Run 21 stress 0.1088547 
    ## Run 22 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283825  max resid 0.07729377 
    ## Run 23 stress 0.113594 
    ## Run 24 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180148  max resid 0.07685829 
    ## Run 25 stress 0.1097988 
    ## Run 26 stress 0.1082722 
    ## Run 27 stress 0.1138391 
    ## Run 28 stress 0.1160543 
    ## Run 29 stress 0.1049903 
    ## ... Procrustes: rmse 0.006720314  max resid 0.02577257 
    ## Run 30 stress 0.1049791 
    ## ... Procrustes: rmse 1.716136e-05  max resid 3.808914e-05 
    ## ... Similar to previous best
    ## Run 31 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283648  max resid 0.07727732 
    ## Run 32 stress 0.1050144 
    ## ... Procrustes: rmse 0.00965929  max resid 0.02607231 
    ## Run 33 stress 0.1138396 
    ## Run 34 stress 0.1082722 
    ## Run 35 stress 0.1084965 
    ## Run 36 stress 0.1079381 
    ## Run 37 stress 0.1063128 
    ## Run 38 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180129  max resid 0.07684606 
    ## Run 39 stress 0.1084964 
    ## Run 40 stress 0.113594 
    ## Run 41 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181221  max resid 0.07689101 
    ## Run 42 stress 0.1082722 
    ## Run 43 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180308  max resid 0.07686227 
    ## Run 44 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179964  max resid 0.07682431 
    ## Run 45 stress 0.1088547 
    ## Run 46 stress 0.1097987 
    ## Run 47 stress 0.1084965 
    ## Run 48 stress 0.1050338 
    ## ... Procrustes: rmse 0.006239115  max resid 0.02179011 
    ## Run 49 stress 0.113594 
    ## Run 50 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283888  max resid 0.07730211 
    ## Run 51 stress 0.1050144 
    ## ... Procrustes: rmse 0.009677723  max resid 0.0261698 
    ## Run 52 stress 0.1063128 
    ## Run 53 stress 0.1050338 
    ## ... Procrustes: rmse 0.006250481  max resid 0.02183625 
    ## Run 54 stress 0.1049791 
    ## ... Procrustes: rmse 2.753299e-05  max resid 8.235069e-05 
    ## ... Similar to previous best
    ## Run 55 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180023  max resid 0.07685182 
    ## Run 56 stress 0.1084964 
    ## Run 57 stress 0.1088547 
    ## Run 58 stress 0.113594 
    ## Run 59 stress 0.1079384 
    ## Run 60 stress 0.1063128 
    ## Run 61 stress 0.1049791 
    ## ... Procrustes: rmse 1.223992e-05  max resid 2.406936e-05 
    ## ... Similar to previous best
    ## Run 62 stress 0.1049903 
    ## ... Procrustes: rmse 0.00672011  max resid 0.02576706 
    ## Run 63 stress 0.1135597 
    ## Run 64 stress 0.1083993 
    ## Run 65 stress 0.107938 
    ## Run 66 stress 0.107938 
    ## Run 67 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180092  max resid 0.0768541 
    ## Run 68 stress 0.1135596 
    ## Run 69 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283777  max resid 0.0772979 
    ## Run 70 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180234  max resid 0.07685929 
    ## Run 71 stress 0.1138393 
    ## Run 72 stress 0.1084965 
    ## Run 73 stress 0.1050144 
    ## ... Procrustes: rmse 0.009657003  max resid 0.02605141 
    ## Run 74 stress 0.1083994 
    ## Run 75 stress 0.1136647 
    ## Run 76 stress 0.1102722 
    ## Run 77 stress 0.1138391 
    ## Run 78 stress 0.1050144 
    ## ... Procrustes: rmse 0.009661776  max resid 0.02609282 
    ## Run 79 stress 0.1063128 
    ## Run 80 stress 0.1088548 
    ## Run 81 stress 0.1135596 
    ## Run 82 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283725  max resid 0.07732124 
    ## Run 83 stress 0.1160546 
    ## Run 84 stress 0.1083503 
    ## Run 85 stress 0.1160552 
    ## Run 86 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283564  max resid 0.0772976 
    ## Run 87 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284193  max resid 0.07733195 
    ## Run 88 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180607  max resid 0.07687104 
    ## Run 89 stress 0.1049791 
    ## ... Procrustes: rmse 1.970944e-06  max resid 6.075642e-06 
    ## ... Similar to previous best
    ## Run 90 stress 0.113594 
    ## Run 91 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283864  max resid 0.07732142 
    ## Run 92 stress 0.1136646 
    ## Run 93 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283794  max resid 0.07728944 
    ## Run 94 stress 0.1051923 
    ## ... Procrustes: rmse 0.022838  max resid 0.07731542 
    ## Run 95 stress 0.1082723 
    ## Run 96 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180123  max resid 0.07685588 
    ## Run 97 stress 0.1136646 
    ## Run 98 stress 0.1063128 
    ## Run 99 stress 0.1050144 
    ## ... Procrustes: rmse 0.009659577  max resid 0.02607251 
    ## Run 100 stress 0.1050144 
    ## ... Procrustes: rmse 0.009660831  max resid 0.02608859 
    ## Run 101 stress 0.1051924 
    ## ... Procrustes: rmse 0.02283823  max resid 0.07728494 
    ## Run 102 stress 0.1082722 
    ## Run 103 stress 0.1160551 
    ## Run 104 stress 0.1136648 
    ## Run 105 stress 0.1135597 
    ## Run 106 stress 0.1138399 
    ## Run 107 stress 0.1063128 
    ## Run 108 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283937  max resid 0.07730409 
    ## Run 109 stress 0.1049792 
    ## ... Procrustes: rmse 3.81878e-05  max resid 0.0001044023 
    ## ... Similar to previous best
    ## Run 110 stress 0.1136646 
    ## Run 111 stress 0.1135597 
    ## Run 112 stress 0.1079384 
    ## Run 113 stress 0.1082722 
    ## Run 114 stress 0.1050144 
    ## ... Procrustes: rmse 0.009658582  max resid 0.02605114 
    ## Run 115 stress 0.1063128 
    ## Run 116 stress 0.1063128 
    ## Run 117 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218004  max resid 0.07682466 
    ## Run 118 stress 0.1135597 
    ## Run 119 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180192  max resid 0.07686074 
    ## Run 120 stress 0.1050144 
    ## ... Procrustes: rmse 0.009663785  max resid 0.02610338 
    ## Run 121 stress 0.1063128 
    ## Run 122 stress 0.1049791 
    ## ... Procrustes: rmse 1.020014e-05  max resid 2.86219e-05 
    ## ... Similar to previous best
    ## Run 123 stress 0.10835 
    ## Run 124 stress 0.1136646 
    ## Run 125 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284406  max resid 0.07725599 
    ## Run 126 stress 0.1097989 
    ## Run 127 stress 0.1063128 
    ## Run 128 stress 0.1063128 
    ## Run 129 stress 0.113594 
    ## Run 130 stress 0.1063128 
    ## Run 131 stress 0.1160539 
    ## Run 132 stress 0.1102724 
    ## Run 133 stress 0.1049903 
    ## ... Procrustes: rmse 0.006769198  max resid 0.02599243 
    ## Run 134 stress 0.1136646 
    ## Run 135 stress 0.1063128 
    ## Run 136 stress 0.1063128 
    ## Run 137 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180022  max resid 0.07683403 
    ## Run 138 stress 0.1135597 
    ## Run 139 stress 0.1084965 
    ## Run 140 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180537  max resid 0.07686771 
    ## Run 141 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179816  max resid 0.0768459 
    ## Run 142 stress 0.1063128 
    ## Run 143 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180109  max resid 0.07685262 
    ## Run 144 stress 0.1083996 
    ## Run 145 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283837  max resid 0.07729659 
    ## Run 146 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179962  max resid 0.07683432 
    ## Run 147 stress 0.1138398 
    ## Run 148 stress 0.1135941 
    ## Run 149 stress 0.1135596 
    ## Run 150 stress 0.1102718 
    ## Run 151 stress 0.1079384 
    ## Run 152 stress 0.1063128 
    ## Run 153 stress 0.1063128 
    ## Run 154 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283803  max resid 0.07729539 
    ## Run 155 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180259  max resid 0.07686008 
    ## Run 156 stress 0.1097989 
    ## Run 157 stress 0.1135597 
    ## Run 158 stress 0.1050338 
    ## ... Procrustes: rmse 0.006258377  max resid 0.02187224 
    ## Run 159 stress 0.1160542 
    ## Run 160 stress 0.1135596 
    ## Run 161 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180151  max resid 0.07686095 
    ## Run 162 stress 0.1136647 
    ## Run 163 stress 0.1063128 
    ## Run 164 stress 0.1083996 
    ## Run 165 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180627  max resid 0.07687399 
    ## Run 166 stress 0.1135597 
    ## Run 167 stress 0.1049903 
    ## ... Procrustes: rmse 0.006721239  max resid 0.02577219 
    ## Run 168 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228375  max resid 0.07730249 
    ## Run 169 stress 0.1050338 
    ## ... Procrustes: rmse 0.006261618  max resid 0.0218864 
    ## Run 170 stress 0.1063128 
    ## Run 171 stress 0.113594 
    ## Run 172 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181294  max resid 0.0768936 
    ## Run 173 stress 0.1063128 
    ## Run 174 stress 0.1063128 
    ## Run 175 stress 0.1049791 
    ## ... Procrustes: rmse 1.229941e-05  max resid 2.678125e-05 
    ## ... Similar to previous best
    ## Run 176 stress 0.1051403 
    ## ... Procrustes: rmse 0.02179948  max resid 0.0768507 
    ## Run 177 stress 0.1063128 
    ## Run 178 stress 0.1063128 
    ## Run 179 stress 0.1050144 
    ## ... Procrustes: rmse 0.009654072  max resid 0.02602798 
    ## Run 180 stress 0.113594 
    ## Run 181 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284401  max resid 0.07728326 
    ## Run 182 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180045  max resid 0.07685345 
    ## Run 183 stress 0.1063128 
    ## Run 184 stress 0.1079385 
    ## Run 185 stress 0.1138395 
    ## Run 186 stress 0.1083996 
    ## Run 187 stress 0.1051402 
    ## ... Procrustes: rmse 0.02178904  max resid 0.07681101 
    ## Run 188 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180189  max resid 0.07685964 
    ## Run 189 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283664  max resid 0.07726544 
    ## Run 190 stress 0.1135596 
    ## Run 191 stress 0.1079381 
    ## Run 192 stress 0.1063128 
    ## Run 193 stress 0.1138393 
    ## Run 194 stress 0.1088547 
    ## Run 195 stress 0.1058423 
    ## Run 196 stress 0.1084965 
    ## Run 197 stress 0.1135596 
    ## Run 198 stress 0.1050338 
    ## ... Procrustes: rmse 0.006252627  max resid 0.02184757 
    ## Run 199 stress 0.1082722 
    ## Run 200 stress 0.1097987 
    ## Run 201 stress 0.1079382 
    ## Run 202 stress 0.1063128 
    ## Run 203 stress 0.1138397 
    ## Run 204 stress 0.1063128 
    ## Run 205 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283913  max resid 0.07730937 
    ## Run 206 stress 0.1138396 
    ## Run 207 stress 0.1165303 
    ## Run 208 stress 0.1049903 
    ## ... Procrustes: rmse 0.006791712  max resid 0.02610052 
    ## Run 209 stress 0.1063128 
    ## Run 210 stress 0.1049903 
    ## ... Procrustes: rmse 0.006776202  max resid 0.02602679 
    ## Run 211 stress 0.1063128 
    ## Run 212 stress 0.1138984 
    ## Run 213 stress 0.113594 
    ## Run 214 stress 0.1083499 
    ## Run 215 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180599  max resid 0.0768715 
    ## Run 216 stress 0.1083502 
    ## Run 217 stress 0.1079381 
    ## Run 218 stress 0.1082722 
    ## Run 219 stress 0.1084966 
    ## Run 220 stress 0.11356 
    ## Run 221 stress 0.1088548 
    ## Run 222 stress 0.1058423 
    ## Run 223 stress 0.1084964 
    ## Run 224 stress 0.1084964 
    ## Run 225 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228386  max resid 0.07730541 
    ## Run 226 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228362  max resid 0.07727568 
    ## Run 227 stress 0.1083502 
    ## Run 228 stress 0.1084964 
    ## Run 229 stress 0.1097988 
    ## Run 230 stress 0.1063128 
    ## Run 231 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179381  max resid 0.07682785 
    ## Run 232 stress 0.1079383 
    ## Run 233 stress 0.1097987 
    ## Run 234 stress 0.1309383 
    ## Run 235 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180266  max resid 0.07686009 
    ## Run 236 stress 0.1063128 
    ## Run 237 stress 0.1138394 
    ## Run 238 stress 0.1063128 
    ## Run 239 stress 0.1084964 
    ## Run 240 stress 0.1097989 
    ## Run 241 stress 0.1063128 
    ## Run 242 stress 0.1084965 
    ## Run 243 stress 0.1088547 
    ## Run 244 stress 0.1049903 
    ## ... Procrustes: rmse 0.006750356  max resid 0.02590473 
    ## Run 245 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218115  max resid 0.07689339 
    ## Run 246 stress 0.1063128 
    ## Run 247 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180675  max resid 0.07687287 
    ## Run 248 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180754  max resid 0.07687538 
    ## Run 249 stress 0.1063128 
    ## Run 250 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180253  max resid 0.07686158 
    ## Run 251 stress 0.1058423 
    ## Run 252 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180175  max resid 0.07685419 
    ## Run 253 stress 0.1063128 
    ## Run 254 stress 0.1136646 
    ## Run 255 stress 0.1136646 
    ## Run 256 stress 0.1058708 
    ## Run 257 stress 0.1084964 
    ## Run 258 stress 0.1135597 
    ## Run 259 stress 0.1079381 
    ## Run 260 stress 0.1088547 
    ## Run 261 stress 0.1084964 
    ## Run 262 stress 0.1088548 
    ## Run 263 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283786  max resid 0.07729949 
    ## Run 264 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180698  max resid 0.07687802 
    ## Run 265 stress 0.1135599 
    ## Run 266 stress 0.1097987 
    ## Run 267 stress 0.1063128 
    ## Run 268 stress 0.1058423 
    ## Run 269 stress 0.1160543 
    ## Run 270 stress 0.1138392 
    ## Run 271 stress 0.113594 
    ## Run 272 stress 0.1088547 
    ## Run 273 stress 0.1051402 
    ## ... Procrustes: rmse 0.021807  max resid 0.07687458 
    ## Run 274 stress 0.113594 
    ## Run 275 stress 0.1082722 
    ## Run 276 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179992  max resid 0.07685121 
    ## Run 277 stress 0.1083501 
    ## Run 278 stress 0.1050144 
    ## ... Procrustes: rmse 0.009668351  max resid 0.02612217 
    ## Run 279 stress 0.1063128 
    ## Run 280 stress 0.1050144 
    ## ... Procrustes: rmse 0.009657153  max resid 0.02606689 
    ## Run 281 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283972  max resid 0.07731553 
    ## Run 282 stress 0.1136646 
    ## Run 283 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180452  max resid 0.07687219 
    ## Run 284 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180356  max resid 0.0768631 
    ## Run 285 stress 0.1063128 
    ## Run 286 stress 0.1050144 
    ## ... Procrustes: rmse 0.009653922  max resid 0.0260262 
    ## Run 287 stress 0.1135598 
    ## Run 288 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283703  max resid 0.07728088 
    ## Run 289 stress 0.1051923 
    ## ... Procrustes: rmse 0.02284  max resid 0.07732321 
    ## Run 290 stress 0.1097989 
    ## Run 291 stress 0.1136647 
    ## Run 292 stress 0.1135598 
    ## Run 293 stress 0.1079384 
    ## Run 294 stress 0.1136646 
    ## Run 295 stress 0.10835 
    ## Run 296 stress 0.1135597 
    ## Run 297 stress 0.1084964 
    ## Run 298 stress 0.1050338 
    ## ... Procrustes: rmse 0.006253089  max resid 0.02184954 
    ## Run 299 stress 0.1063128 
    ## Run 300 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228381  max resid 0.07729773 
    ## Run 301 stress 0.1063128 
    ## Run 302 stress 0.1135596 
    ## Run 303 stress 0.113594 
    ## Run 304 stress 0.1079381 
    ## Run 305 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180189  max resid 0.07685672 
    ## Run 306 stress 0.1088547 
    ## Run 307 stress 0.1165298 
    ## Run 308 stress 0.1084965 
    ## Run 309 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180118  max resid 0.07684548 
    ## Run 310 stress 0.113594 
    ## Run 311 stress 0.1082722 
    ## Run 312 stress 0.1097987 
    ## Run 313 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180183  max resid 0.07685819 
    ## Run 314 stress 0.1138393 
    ## Run 315 stress 0.1050144 
    ## ... Procrustes: rmse 0.009661174  max resid 0.02608469 
    ## Run 316 stress 0.113594 
    ## Run 317 stress 0.1136646 
    ## Run 318 stress 0.1136646 
    ## Run 319 stress 0.1097987 
    ## Run 320 stress 0.1082722 
    ## Run 321 stress 0.1138397 
    ## Run 322 stress 0.1088548 
    ## Run 323 stress 0.1309385 
    ## Run 324 stress 0.1050144 
    ## ... Procrustes: rmse 0.009686441  max resid 0.02623866 
    ## Run 325 stress 0.1138395 
    ## Run 326 stress 0.113594 
    ## Run 327 stress 0.1102721 
    ## Run 328 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180007  max resid 0.07683431 
    ## Run 329 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180197  max resid 0.07686048 
    ## Run 330 stress 0.1083499 
    ## Run 331 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283768  max resid 0.07728913 
    ## Run 332 stress 0.1050338 
    ## ... Procrustes: rmse 0.006259475  max resid 0.0218778 
    ## Run 333 stress 0.1058423 
    ## Run 334 stress 0.1082722 
    ## Run 335 stress 0.1050338 
    ## ... Procrustes: rmse 0.006253057  max resid 0.02185239 
    ## Run 336 stress 0.1050338 
    ## ... Procrustes: rmse 0.006266625  max resid 0.02190388 
    ## Run 337 stress 0.1102724 
    ## Run 338 stress 0.1083998 
    ## Run 339 stress 0.1063128 
    ## Run 340 stress 0.1049791 
    ## ... Procrustes: rmse 1.0022e-05  max resid 2.870208e-05 
    ## ... Similar to previous best
    ## Run 341 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180338  max resid 0.07686215 
    ## Run 342 stress 0.1063128 
    ## Run 343 stress 0.1049791 
    ## ... Procrustes: rmse 4.166192e-05  max resid 0.0001175641 
    ## ... Similar to previous best
    ## Run 344 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180542  max resid 0.07686759 
    ## Run 345 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180861  max resid 0.07688026 
    ## Run 346 stress 0.1135597 
    ## Run 347 stress 0.1050144 
    ## ... Procrustes: rmse 0.009634441  max resid 0.02593181 
    ## Run 348 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283956  max resid 0.07730884 
    ## Run 349 stress 0.113899 
    ## Run 350 stress 0.1309402 
    ## Run 351 stress 0.1058423 
    ## Run 352 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180437  max resid 0.07686521 
    ## Run 353 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284322  max resid 0.07727844 
    ## Run 354 stress 0.1051402 
    ## ... Procrustes: rmse 0.0218037  max resid 0.07685958 
    ## Run 355 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180703  max resid 0.07687493 
    ## Run 356 stress 0.1050338 
    ## ... Procrustes: rmse 0.006274082  max resid 0.02193164 
    ## Run 357 stress 0.1049791 
    ## ... Procrustes: rmse 3.872742e-06  max resid 1.245145e-05 
    ## ... Similar to previous best
    ## Run 358 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180186  max resid 0.07685791 
    ## Run 359 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180216  max resid 0.07685903 
    ## Run 360 stress 0.1079381 
    ## Run 361 stress 0.113594 
    ## Run 362 stress 0.1135597 
    ## Run 363 stress 0.1136646 
    ## Run 364 stress 0.1136646 
    ## Run 365 stress 0.1049791 
    ## ... Procrustes: rmse 1.653559e-05  max resid 4.841141e-05 
    ## ... Similar to previous best
    ## Run 366 stress 0.1136646 
    ## Run 367 stress 0.1050338 
    ## ... Procrustes: rmse 0.006258631  max resid 0.0218724 
    ## Run 368 stress 0.1084966 
    ## Run 369 stress 0.1084965 
    ## Run 370 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179942  max resid 0.07682333 
    ## Run 371 stress 0.1097987 
    ## Run 372 stress 0.1136646 
    ## Run 373 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180212  max resid 0.07686689 
    ## Run 374 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283814  max resid 0.0773062 
    ## Run 375 stress 0.1051924 
    ## ... Procrustes: rmse 0.02284536  max resid 0.0772858 
    ## Run 376 stress 0.1135941 
    ## Run 377 stress 0.1050144 
    ## ... Procrustes: rmse 0.009663594  max resid 0.02611739 
    ## Run 378 stress 0.1063128 
    ## Run 379 stress 0.1083499 
    ## Run 380 stress 0.1136646 
    ## Run 381 stress 0.1135597 
    ## Run 382 stress 0.1097989 
    ## Run 383 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180307  max resid 0.07686235 
    ## Run 384 stress 0.1138392 
    ## Run 385 stress 0.1160544 
    ## Run 386 stress 0.1136647 
    ## Run 387 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180204  max resid 0.07686027 
    ## Run 388 stress 0.1083501 
    ## Run 389 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283821  max resid 0.07729348 
    ## Run 390 stress 0.1079381 
    ## Run 391 stress 0.1160555 
    ## Run 392 stress 0.1083502 
    ## Run 393 stress 0.1102719 
    ## Run 394 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180264  max resid 0.07686197 
    ## Run 395 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180259  max resid 0.07685971 
    ## Run 396 stress 0.1082723 
    ## Run 397 stress 0.1136646 
    ## Run 398 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283829  max resid 0.07730729 
    ## Run 399 stress 0.1083501 
    ## Run 400 stress 0.1135942 
    ## Run 401 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180502  max resid 0.07686364 
    ## Run 402 stress 0.113594 
    ## Run 403 stress 0.1097988 
    ## Run 404 stress 0.116055 
    ## Run 405 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179931  max resid 0.07682085 
    ## Run 406 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180287  max resid 0.07686187 
    ## Run 407 stress 0.107938 
    ## Run 408 stress 0.1049791 
    ## ... Procrustes: rmse 1.725429e-05  max resid 3.929654e-05 
    ## ... Similar to previous best
    ## Run 409 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180918  max resid 0.07688197 
    ## Run 410 stress 0.1088547 
    ## Run 411 stress 0.1102722 
    ## Run 412 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283823  max resid 0.07726149 
    ## Run 413 stress 0.1049903 
    ## ... Procrustes: rmse 0.006727513  max resid 0.02579814 
    ## Run 414 stress 0.1135597 
    ## Run 415 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180383  max resid 0.07688777 
    ## Run 416 stress 0.1079382 
    ## Run 417 stress 0.1136646 
    ## Run 418 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283788  max resid 0.07729515 
    ## Run 419 stress 0.1135596 
    ## Run 420 stress 0.1084965 
    ## Run 421 stress 0.1088547 
    ## Run 422 stress 0.1136646 
    ## Run 423 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180148  max resid 0.07684797 
    ## Run 424 stress 0.1049903 
    ## ... Procrustes: rmse 0.006774452  max resid 0.02601866 
    ## Run 425 stress 0.1160545 
    ## Run 426 stress 0.1135596 
    ## Run 427 stress 0.113594 
    ## Run 428 stress 0.1138393 
    ## Run 429 stress 0.1063128 
    ## Run 430 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180391  max resid 0.07686377 
    ## Run 431 stress 0.1160538 
    ## Run 432 stress 0.1097987 
    ## Run 433 stress 0.1083503 
    ## Run 434 stress 0.1097989 
    ## Run 435 stress 0.1063128 
    ## Run 436 stress 0.113899 
    ## Run 437 stress 0.1084964 
    ## Run 438 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228358  max resid 0.07726914 
    ## Run 439 stress 0.1079381 
    ## Run 440 stress 0.1160549 
    ## Run 441 stress 0.1084965 
    ## Run 442 stress 0.1160549 
    ## Run 443 stress 0.1136646 
    ## Run 444 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283786  max resid 0.07729142 
    ## Run 445 stress 0.1088547 
    ## Run 446 stress 0.10835 
    ## Run 447 stress 0.113594 
    ## Run 448 stress 0.1050144 
    ## ... Procrustes: rmse 0.009656419  max resid 0.02604881 
    ## Run 449 stress 0.1079384 
    ## Run 450 stress 0.1136646 
    ## Run 451 stress 0.1136646 
    ## Run 452 stress 0.113594 
    ## Run 453 stress 0.1082722 
    ## Run 454 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180383  max resid 0.07686612 
    ## Run 455 stress 0.1135598 
    ## Run 456 stress 0.1063128 
    ## Run 457 stress 0.1138397 
    ## Run 458 stress 0.1082722 
    ## Run 459 stress 0.113594 
    ## Run 460 stress 0.1050338 
    ## ... Procrustes: rmse 0.006253436  max resid 0.02184722 
    ## Run 461 stress 0.1136646 
    ## Run 462 stress 0.1063128 
    ## Run 463 stress 0.1049903 
    ## ... Procrustes: rmse 0.006762411  max resid 0.02596129 
    ## Run 464 stress 0.1049903 
    ## ... Procrustes: rmse 0.006723738  max resid 0.0257831 
    ## Run 465 stress 0.1135596 
    ## Run 466 stress 0.1138394 
    ## Run 467 stress 0.1136646 
    ## Run 468 stress 0.1083501 
    ## Run 469 stress 0.1051923 
    ## ... Procrustes: rmse 0.0228378  max resid 0.07729692 
    ## Run 470 stress 0.1051923 
    ## ... Procrustes: rmse 0.02283808  max resid 0.07730348 
    ## Run 471 stress 0.1049903 
    ## ... Procrustes: rmse 0.006743213  max resid 0.02586989 
    ## Run 472 stress 0.1136646 
    ## Run 473 stress 0.1138394 
    ## Run 474 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179671  max resid 0.07684017 
    ## Run 475 stress 0.1058423 
    ## Run 476 stress 0.113594 
    ## Run 477 stress 0.1083501 
    ## Run 478 stress 0.10835 
    ## Run 479 stress 0.1136646 
    ## Run 480 stress 0.1051402 
    ## ... Procrustes: rmse 0.02179671  max resid 0.07683948 
    ## Run 481 stress 0.1136647 
    ## Run 482 stress 0.1097988 
    ## Run 483 stress 0.1083995 
    ## Run 484 stress 0.1088547 
    ## Run 485 stress 0.1102722 
    ## Run 486 stress 0.1063128 
    ## Run 487 stress 0.107938 
    ## Run 488 stress 0.1084965 
    ## Run 489 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180211  max resid 0.07685986 
    ## Run 490 stress 0.1058423 
    ## Run 491 stress 0.1063128 
    ## Run 492 stress 0.1097987 
    ## Run 493 stress 0.1051402 
    ## ... Procrustes: rmse 0.02181007  max resid 0.07688487 
    ## Run 494 stress 0.1084965 
    ## Run 495 stress 0.1160543 
    ## Run 496 stress 0.1063128 
    ## Run 497 stress 0.1063128 
    ## Run 498 stress 0.1051402 
    ## ... Procrustes: rmse 0.02180083  max resid 0.07685437 
    ## Run 499 stress 0.1084964 
    ## Run 500 stress 0.1063128 
    ## *** Best solution repeated 12 times

``` r
# Surveyed sites 
SD_beta_NMDS <- metaMDS(SD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09021173 
    ## Run 1 stress 0.1064187 
    ## Run 2 stress 0.09109111 
    ## Run 3 stress 0.1064191 
    ## Run 4 stress 0.1092561 
    ## Run 5 stress 0.1062559 
    ## Run 6 stress 0.1071322 
    ## Run 7 stress 0.08946329 
    ## ... New best solution
    ## ... Procrustes: rmse 0.02371736  max resid 0.07465998 
    ## Run 8 stress 0.09178292 
    ## Run 9 stress 0.0892607 
    ## ... New best solution
    ## ... Procrustes: rmse 0.03751439  max resid 0.1170224 
    ## Run 10 stress 0.0908734 
    ## Run 11 stress 0.1111256 
    ## Run 12 stress 0.08938964 
    ## ... Procrustes: rmse 0.03598786  max resid 0.1183542 
    ## Run 13 stress 0.1091234 
    ## Run 14 stress 0.09018707 
    ## Run 15 stress 0.0904461 
    ## Run 16 stress 0.1096139 
    ## Run 17 stress 0.1075421 
    ## Run 18 stress 0.08938542 
    ## ... Procrustes: rmse 0.01174898  max resid 0.03959615 
    ## Run 19 stress 0.08951722 
    ## ... Procrustes: rmse 0.03510494  max resid 0.1157435 
    ## Run 20 stress 0.08938545 
    ## ... Procrustes: rmse 0.01182601  max resid 0.03960693 
    ## Run 21 stress 0.08946652 
    ## ... Procrustes: rmse 0.03325184  max resid 0.1164133 
    ## Run 22 stress 0.1071322 
    ## Run 23 stress 0.09039136 
    ## Run 24 stress 0.1092127 
    ## Run 25 stress 0.1079019 
    ## Run 26 stress 0.08938967 
    ## ... Procrustes: rmse 0.03599186  max resid 0.118367 
    ## Run 27 stress 0.1087866 
    ## Run 28 stress 0.09503427 
    ## Run 29 stress 0.1056904 
    ## Run 30 stress 0.09044609 
    ## Run 31 stress 0.0894633 
    ## ... Procrustes: rmse 0.03754412  max resid 0.1178315 
    ## Run 32 stress 0.09109108 
    ## Run 33 stress 0.1067489 
    ## Run 34 stress 0.109188 
    ## Run 35 stress 0.1060431 
    ## Run 36 stress 0.08938982 
    ## ... Procrustes: rmse 0.03595188  max resid 0.118302 
    ## Run 37 stress 0.08938962 
    ## ... Procrustes: rmse 0.03596  max resid 0.1183078 
    ## Run 38 stress 0.1062599 
    ## Run 39 stress 0.1052652 
    ## Run 40 stress 0.08946332 
    ## ... Procrustes: rmse 0.03750855  max resid 0.1177876 
    ## Run 41 stress 0.1109762 
    ## Run 42 stress 0.1075421 
    ## Run 43 stress 0.09039133 
    ## Run 44 stress 0.08926076 
    ## ... Procrustes: rmse 6.949834e-05  max resid 0.0001864216 
    ## ... Similar to previous best
    ## Run 45 stress 0.090886 
    ## Run 46 stress 0.08946671 
    ## ... Procrustes: rmse 0.03319159  max resid 0.1163209 
    ## Run 47 stress 0.09109116 
    ## Run 48 stress 0.1061298 
    ## Run 49 stress 0.1104204 
    ## Run 50 stress 0.09228506 
    ## Run 51 stress 0.0893855 
    ## ... Procrustes: rmse 0.01172965  max resid 0.03967659 
    ## Run 52 stress 0.08926065 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002505589  max resid 0.000791078 
    ## ... Similar to previous best
    ## Run 53 stress 0.1052653 
    ## Run 54 stress 0.109256 
    ## Run 55 stress 0.10613 
    ## Run 56 stress 0.1052649 
    ## Run 57 stress 0.08946663 
    ## ... Procrustes: rmse 0.03315627  max resid 0.1162512 
    ## Run 58 stress 0.1060428 
    ## Run 59 stress 0.08938967 
    ## ... Procrustes: rmse 0.03590309  max resid 0.1182185 
    ## Run 60 stress 0.1056897 
    ## Run 61 stress 0.08938546 
    ## ... Procrustes: rmse 0.01186326  max resid 0.0395416 
    ## Run 62 stress 0.1052651 
    ## Run 63 stress 0.0904461 
    ## Run 64 stress 0.0910911 
    ## Run 65 stress 0.1071321 
    ## Run 66 stress 0.1076305 
    ## Run 67 stress 0.08951729 
    ## ... Procrustes: rmse 0.03498916  max resid 0.1155668 
    ## Run 68 stress 0.08926068 
    ## ... Procrustes: rmse 0.0002419168  max resid 0.0006040736 
    ## ... Similar to previous best
    ## Run 69 stress 0.0893897 
    ## ... Procrustes: rmse 0.03589822  max resid 0.1182108 
    ## Run 70 stress 0.08926089 
    ## ... Procrustes: rmse 0.0004195952  max resid 0.001291297 
    ## ... Similar to previous best
    ## Run 71 stress 0.1064427 
    ## Run 72 stress 0.1052649 
    ## Run 73 stress 0.1075422 
    ## Run 74 stress 0.08938965 
    ## ... Procrustes: rmse 0.03592758  max resid 0.1182586 
    ## Run 75 stress 0.08938962 
    ## ... Procrustes: rmse 0.03591884  max resid 0.1182433 
    ## Run 76 stress 0.0910911 
    ## Run 77 stress 0.1076302 
    ## Run 78 stress 0.08926076 
    ## ... Procrustes: rmse 0.0003357345  max resid 0.0009828273 
    ## ... Similar to previous best
    ## Run 79 stress 0.1065372 
    ## Run 80 stress 0.1074434 
    ## Run 81 stress 0.09018707 
    ## Run 82 stress 0.09109108 
    ## Run 83 stress 0.1071322 
    ## Run 84 stress 0.09178286 
    ## Run 85 stress 0.09503424 
    ## Run 86 stress 0.1074434 
    ## Run 87 stress 0.08926082 
    ## ... Procrustes: rmse 0.0003774175  max resid 0.00114511 
    ## ... Similar to previous best
    ## Run 88 stress 0.09228488 
    ## Run 89 stress 0.08938547 
    ## ... Procrustes: rmse 0.01175093  max resid 0.03952814 
    ## Run 90 stress 0.08946651 
    ## ... Procrustes: rmse 0.03319994  max resid 0.1163267 
    ## Run 91 stress 0.08938967 
    ## ... Procrustes: rmse 0.03593655  max resid 0.1182679 
    ## Run 92 stress 0.106751 
    ## Run 93 stress 0.09039132 
    ## Run 94 stress 0.1056895 
    ## Run 95 stress 0.1067514 
    ## Run 96 stress 0.08938961 
    ## ... Procrustes: rmse 0.03592959  max resid 0.1182603 
    ## Run 97 stress 0.1096131 
    ## Run 98 stress 0.09130089 
    ## Run 99 stress 0.1085126 
    ## Run 100 stress 0.1074432 
    ## Run 101 stress 0.1074434 
    ## Run 102 stress 0.0894666 
    ## ... Procrustes: rmse 0.03316249  max resid 0.116267 
    ## Run 103 stress 0.1060093 
    ## Run 104 stress 0.08938974 
    ## ... Procrustes: rmse 0.03589244  max resid 0.1182014 
    ## Run 105 stress 0.08946331 
    ## ... Procrustes: rmse 0.03747006  max resid 0.117714 
    ## Run 106 stress 0.1064199 
    ## Run 107 stress 0.1086561 
    ## Run 108 stress 0.09018707 
    ## Run 109 stress 0.08938537 
    ## ... Procrustes: rmse 0.0118199  max resid 0.0395918 
    ## Run 110 stress 0.08926077 
    ## ... Procrustes: rmse 0.0003227947  max resid 0.001010131 
    ## ... Similar to previous best
    ## Run 111 stress 0.09109112 
    ## Run 112 stress 0.1076303 
    ## Run 113 stress 0.10613 
    ## Run 114 stress 0.09262352 
    ## Run 115 stress 0.09087339 
    ## Run 116 stress 0.1108231 
    ## Run 117 stress 0.1087565 
    ## Run 118 stress 0.09039129 
    ## Run 119 stress 0.1092126 
    ## Run 120 stress 0.08926088 
    ## ... Procrustes: rmse 0.0003983239  max resid 0.001218505 
    ## ... Similar to previous best
    ## Run 121 stress 0.09099531 
    ## Run 122 stress 0.0892608 
    ## ... Procrustes: rmse 0.000363965  max resid 0.001071077 
    ## ... Similar to previous best
    ## Run 123 stress 0.1056897 
    ## Run 124 stress 0.1091449 
    ## Run 125 stress 0.1071322 
    ## Run 126 stress 0.1067378 
    ## Run 127 stress 0.1096133 
    ## Run 128 stress 0.09109116 
    ## Run 129 stress 0.1067378 
    ## Run 130 stress 0.08938539 
    ## ... Procrustes: rmse 0.01182203  max resid 0.03961208 
    ## Run 131 stress 0.08938962 
    ## ... Procrustes: rmse 0.03592782  max resid 0.1182561 
    ## Run 132 stress 0.1087563 
    ## Run 133 stress 0.09021165 
    ## Run 134 stress 0.08938542 
    ## ... Procrustes: rmse 0.01184276  max resid 0.03963253 
    ## Run 135 stress 0.109213 
    ## Run 136 stress 0.0892607 
    ## ... Procrustes: rmse 5.763398e-05  max resid 0.0001668036 
    ## ... Similar to previous best
    ## Run 137 stress 0.09612818 
    ## Run 138 stress 0.09087365 
    ## Run 139 stress 0.09130086 
    ## Run 140 stress 0.1084535 
    ## Run 141 stress 0.1092127 
    ## Run 142 stress 0.08926066 
    ## ... Procrustes: rmse 1.918472e-05  max resid 5.196776e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.1071322 
    ## Run 144 stress 0.1118951 
    ## Run 145 stress 0.106419 
    ## Run 146 stress 0.09094226 
    ## Run 147 stress 0.1060089 
    ## Run 148 stress 0.106009 
    ## Run 149 stress 0.1067367 
    ## Run 150 stress 0.08938961 
    ## ... Procrustes: rmse 0.03593318  max resid 0.1182676 
    ## Run 151 stress 0.08938542 
    ## ... Procrustes: rmse 0.01182668  max resid 0.03961714 
    ## Run 152 stress 0.1108633 
    ## Run 153 stress 0.09088615 
    ## Run 154 stress 0.1092096 
    ## Run 155 stress 0.1071321 
    ## Run 156 stress 0.1074432 
    ## Run 157 stress 0.09180895 
    ## Run 158 stress 0.1067896 
    ## Run 159 stress 0.08926109 
    ## ... Procrustes: rmse 0.0005241754  max resid 0.001653022 
    ## ... Similar to previous best
    ## Run 160 stress 0.09039139 
    ## Run 161 stress 0.1111381 
    ## Run 162 stress 0.1092093 
    ## Run 163 stress 0.08938538 
    ## ... Procrustes: rmse 0.01181632  max resid 0.03955822 
    ## Run 164 stress 0.1076304 
    ## Run 165 stress 0.1065348 
    ## Run 166 stress 0.1092127 
    ## Run 167 stress 0.1092364 
    ## Run 168 stress 0.08946328 
    ## ... Procrustes: rmse 0.03748911  max resid 0.1177415 
    ## Run 169 stress 0.08938969 
    ## ... Procrustes: rmse 0.03590152  max resid 0.1182115 
    ## Run 170 stress 0.08946651 
    ## ... Procrustes: rmse 0.0331999  max resid 0.1163328 
    ## Run 171 stress 0.110864 
    ## Run 172 stress 0.1105335 
    ## Run 173 stress 0.08938965 
    ## ... Procrustes: rmse 0.0359081  max resid 0.1182263 
    ## Run 174 stress 0.09088604 
    ## Run 175 stress 0.09228483 
    ## Run 176 stress 0.09087339 
    ## Run 177 stress 0.09021166 
    ## Run 178 stress 0.09503437 
    ## Run 179 stress 0.09018707 
    ## Run 180 stress 0.09018709 
    ## Run 181 stress 0.09099566 
    ## Run 182 stress 0.1061303 
    ## Run 183 stress 0.1063192 
    ## Run 184 stress 0.1060431 
    ## Run 185 stress 0.09503435 
    ## Run 186 stress 0.1062562 
    ## Run 187 stress 0.1087573 
    ## Run 188 stress 0.0908734 
    ## Run 189 stress 0.1061305 
    ## Run 190 stress 0.09503422 
    ## Run 191 stress 0.1068784 
    ## Run 192 stress 0.09503431 
    ## Run 193 stress 0.1052649 
    ## Run 194 stress 0.1071321 
    ## Run 195 stress 0.08946666 
    ## ... Procrustes: rmse 0.033157  max resid 0.1162598 
    ## Run 196 stress 0.1092127 
    ## Run 197 stress 0.09503431 
    ## Run 198 stress 0.1052651 
    ## Run 199 stress 0.09592132 
    ## Run 200 stress 0.106783 
    ## Run 201 stress 0.09775539 
    ## Run 202 stress 0.08946653 
    ## ... Procrustes: rmse 0.03320611  max resid 0.1163537 
    ## Run 203 stress 0.1092126 
    ## Run 204 stress 0.1088367 
    ## Run 205 stress 0.1101441 
    ## Run 206 stress 0.09039136 
    ## Run 207 stress 0.08946328 
    ## ... Procrustes: rmse 0.03749438  max resid 0.1177624 
    ## Run 208 stress 0.09130087 
    ## Run 209 stress 0.1067374 
    ## Run 210 stress 0.0892608 
    ## ... Procrustes: rmse 0.0003668947  max resid 0.001080854 
    ## ... Similar to previous best
    ## Run 211 stress 0.0910911 
    ## Run 212 stress 0.1087572 
    ## Run 213 stress 0.1075422 
    ## Run 214 stress 0.0894633 
    ## ... Procrustes: rmse 0.03747288  max resid 0.1177189 
    ## Run 215 stress 0.08938561 
    ## ... Procrustes: rmse 0.01171936  max resid 0.03947213 
    ## Run 216 stress 0.09087362 
    ## Run 217 stress 0.1091824 
    ## Run 218 stress 0.1091238 
    ## Run 219 stress 0.1056902 
    ## Run 220 stress 0.1104473 
    ## Run 221 stress 0.08946653 
    ## ... Procrustes: rmse 0.03317697  max resid 0.1162928 
    ## Run 222 stress 0.09087345 
    ## Run 223 stress 0.09109108 
    ## Run 224 stress 0.0893854 
    ## ... Procrustes: rmse 0.01183535  max resid 0.03958683 
    ## Run 225 stress 0.1074434 
    ## Run 226 stress 0.09178305 
    ## Run 227 stress 0.1091449 
    ## Run 228 stress 0.08946335 
    ## ... Procrustes: rmse 0.03751693  max resid 0.1178012 
    ## Run 229 stress 0.1075421 
    ## Run 230 stress 0.08946333 
    ## ... Procrustes: rmse 0.03746473  max resid 0.1177083 
    ## Run 231 stress 0.09503463 
    ## Run 232 stress 0.08946652 
    ## ... Procrustes: rmse 0.03318281  max resid 0.1163034 
    ## Run 233 stress 0.1079018 
    ## Run 234 stress 0.08946653 
    ## ... Procrustes: rmse 0.03321217  max resid 0.1163627 
    ## Run 235 stress 0.0908735 
    ## Run 236 stress 0.1071321 
    ## Run 237 stress 0.1092562 
    ## Run 238 stress 0.08926104 
    ## ... Procrustes: rmse 0.0005116088  max resid 0.001566366 
    ## ... Similar to previous best
    ## Run 239 stress 0.08926099 
    ## ... Procrustes: rmse 0.0004853805  max resid 0.001519917 
    ## ... Similar to previous best
    ## Run 240 stress 0.1074433 
    ## Run 241 stress 0.1074434 
    ## Run 242 stress 0.09018712 
    ## Run 243 stress 0.1092948 
    ## Run 244 stress 0.1089907 
    ## Run 245 stress 0.1061298 
    ## Run 246 stress 0.0892607 
    ## ... Procrustes: rmse 0.0002733583  max resid 0.0007840606 
    ## ... Similar to previous best
    ## Run 247 stress 0.08946652 
    ## ... Procrustes: rmse 0.03318208  max resid 0.1163008 
    ## Run 248 stress 0.1060428 
    ## Run 249 stress 0.1071322 
    ## Run 250 stress 0.1067423 
    ## Run 251 stress 0.1087568 
    ## Run 252 stress 0.1066194 
    ## Run 253 stress 0.0893855 
    ## ... Procrustes: rmse 0.01176509  max resid 0.03949182 
    ## Run 254 stress 0.105265 
    ## Run 255 stress 0.09503421 
    ## Run 256 stress 0.1092209 
    ## Run 257 stress 0.08938538 
    ## ... Procrustes: rmse 0.01182042  max resid 0.03959618 
    ## Run 258 stress 0.1087568 
    ## Run 259 stress 0.09503445 
    ## Run 260 stress 0.107902 
    ## Run 261 stress 0.09039135 
    ## Run 262 stress 0.08926071 
    ## ... Procrustes: rmse 0.0002847728  max resid 0.0008303935 
    ## ... Similar to previous best
    ## Run 263 stress 0.08938965 
    ## ... Procrustes: rmse 0.03590855  max resid 0.1182264 
    ## Run 264 stress 0.1056905 
    ## Run 265 stress 0.09021171 
    ## Run 266 stress 0.1056904 
    ## Run 267 stress 0.09503414 
    ## Run 268 stress 0.08938542 
    ## ... Procrustes: rmse 0.01184611  max resid 0.03955084 
    ## Run 269 stress 0.111138 
    ## Run 270 stress 0.09087341 
    ## Run 271 stress 0.08938541 
    ## ... Procrustes: rmse 0.01184472  max resid 0.03957108 
    ## Run 272 stress 0.09130085 
    ## Run 273 stress 0.1074434 
    ## Run 274 stress 0.1062565 
    ## Run 275 stress 0.0893854 
    ## ... Procrustes: rmse 0.01181248  max resid 0.03959877 
    ## Run 276 stress 0.09099586 
    ## Run 277 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001441751  max resid 0.0003248746 
    ## ... Similar to previous best
    ## Run 278 stress 0.09039137 
    ## Run 279 stress 0.1056892 
    ## Run 280 stress 0.1056906 
    ## Run 281 stress 0.08926095 
    ## ... Procrustes: rmse 0.0004084812  max resid 0.0011011 
    ## ... Similar to previous best
    ## Run 282 stress 0.0902117 
    ## Run 283 stress 0.1071322 
    ## Run 284 stress 0.08938961 
    ## ... Procrustes: rmse 0.03596527  max resid 0.1183055 
    ## Run 285 stress 0.08926099 
    ## ... Procrustes: rmse 0.0004321444  max resid 0.00125082 
    ## ... Similar to previous best
    ## Run 286 stress 0.09503415 
    ## Run 287 stress 0.08946658 
    ## ... Procrustes: rmse 0.0332019  max resid 0.1163074 
    ## Run 288 stress 0.08946661 
    ## ... Procrustes: rmse 0.03319872  max resid 0.1163075 
    ## Run 289 stress 0.08926101 
    ## ... Procrustes: rmse 0.0003591618  max resid 0.00113547 
    ## ... Similar to previous best
    ## Run 290 stress 0.08926074 
    ## ... Procrustes: rmse 0.0001885216  max resid 0.0006168842 
    ## ... Similar to previous best
    ## Run 291 stress 0.1061302 
    ## Run 292 stress 0.08938545 
    ## ... Procrustes: rmse 0.0118538  max resid 0.03965981 
    ## Run 293 stress 0.1056903 
    ## Run 294 stress 0.09088623 
    ## Run 295 stress 0.0893854 
    ## ... Procrustes: rmse 0.01182003  max resid 0.03974288 
    ## Run 296 stress 0.1099387 
    ## Run 297 stress 0.1052648 
    ## Run 298 stress 0.1074432 
    ## Run 299 stress 0.1087582 
    ## Run 300 stress 0.08938543 
    ## ... Procrustes: rmse 0.01188493  max resid 0.03970331 
    ## Run 301 stress 0.1052647 
    ## Run 302 stress 0.09503453 
    ## Run 303 stress 0.09099542 
    ## Run 304 stress 0.1101446 
    ## Run 305 stress 0.1075421 
    ## Run 306 stress 0.1092128 
    ## Run 307 stress 0.1088391 
    ## Run 308 stress 0.09088612 
    ## Run 309 stress 0.1068389 
    ## Run 310 stress 0.08926105 
    ## ... Procrustes: rmse 0.0003820106  max resid 0.001276081 
    ## ... Similar to previous best
    ## Run 311 stress 0.1075422 
    ## Run 312 stress 0.1089198 
    ## Run 313 stress 0.1056905 
    ## Run 314 stress 0.0893856 
    ## ... Procrustes: rmse 0.01191744  max resid 0.03966652 
    ## Run 315 stress 0.1062596 
    ## Run 316 stress 0.1071322 
    ## Run 317 stress 0.09088615 
    ## Run 318 stress 0.1092127 
    ## Run 319 stress 0.1111382 
    ## Run 320 stress 0.08938542 
    ## ... Procrustes: rmse 0.01183464  max resid 0.03975291 
    ## Run 321 stress 0.09775588 
    ## Run 322 stress 0.09130086 
    ## Run 323 stress 0.08938538 
    ## ... Procrustes: rmse 0.01185996  max resid 0.03975162 
    ## Run 324 stress 0.1056892 
    ## Run 325 stress 0.1074431 
    ## Run 326 stress 0.1079015 
    ## Run 327 stress 0.09021166 
    ## Run 328 stress 0.1052652 
    ## Run 329 stress 0.08938963 
    ## ... Procrustes: rmse 0.03594494  max resid 0.1182769 
    ## Run 330 stress 0.1075422 
    ## Run 331 stress 0.0904462 
    ## Run 332 stress 0.08946335 
    ## ... Procrustes: rmse 0.03749935  max resid 0.1177342 
    ## Run 333 stress 0.09503417 
    ## Run 334 stress 0.1065352 
    ## Run 335 stress 0.1065366 
    ## Run 336 stress 0.1052654 
    ## Run 337 stress 0.0903913 
    ## Run 338 stress 0.1091506 
    ## Run 339 stress 0.08946653 
    ## ... Procrustes: rmse 0.03321246  max resid 0.1163251 
    ## Run 340 stress 0.1071321 
    ## Run 341 stress 0.08938969 
    ## ... Procrustes: rmse 0.03593362  max resid 0.1182555 
    ## Run 342 stress 0.1065356 
    ## Run 343 stress 0.1071321 
    ## Run 344 stress 0.1062595 
    ## Run 345 stress 0.1086562 
    ## Run 346 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183133  max resid 0.03972079 
    ## Run 347 stress 0.08926082 
    ## ... Procrustes: rmse 0.0002558175  max resid 0.0008293805 
    ## ... Similar to previous best
    ## Run 348 stress 0.1071323 
    ## Run 349 stress 0.09087349 
    ## Run 350 stress 0.08926089 
    ## ... Procrustes: rmse 0.0003016079  max resid 0.0009997254 
    ## ... Similar to previous best
    ## Run 351 stress 0.1075422 
    ## Run 352 stress 0.09130085 
    ## Run 353 stress 0.09503433 
    ## Run 354 stress 0.08926065 
    ## ... Procrustes: rmse 0.0001424252  max resid 0.0003237188 
    ## ... Similar to previous best
    ## Run 355 stress 0.1108316 
    ## Run 356 stress 0.08946331 
    ## ... Procrustes: rmse 0.0375093  max resid 0.117754 
    ## Run 357 stress 0.1056903 
    ## Run 358 stress 0.0893898 
    ## ... Procrustes: rmse 0.03591624  max resid 0.1182298 
    ## Run 359 stress 0.09021166 
    ## Run 360 stress 0.08926086 
    ## ... Procrustes: rmse 0.0002799202  max resid 0.0009292525 
    ## ... Similar to previous best
    ## Run 361 stress 0.08946653 
    ## ... Procrustes: rmse 0.03321382  max resid 0.1163296 
    ## Run 362 stress 0.110146 
    ## Run 363 stress 0.1092091 
    ## Run 364 stress 0.1056895 
    ## Run 365 stress 0.1074432 
    ## Run 366 stress 0.08938539 
    ## ... Procrustes: rmse 0.01181938  max resid 0.03972226 
    ## Run 367 stress 0.09503437 
    ## Run 368 stress 0.08926075 
    ## ... Procrustes: rmse 0.0001950582  max resid 0.0006619505 
    ## ... Similar to previous best
    ## Run 369 stress 0.1061302 
    ## Run 370 stress 0.1092095 
    ## Run 371 stress 0.08938544 
    ## ... Procrustes: rmse 0.01187623  max resid 0.03978793 
    ## Run 372 stress 0.1056899 
    ## Run 373 stress 0.08938542 
    ## ... Procrustes: rmse 0.011867  max resid 0.03977451 
    ## Run 374 stress 0.1076302 
    ## Run 375 stress 0.08938961 
    ## ... Procrustes: rmse 0.03596232  max resid 0.1182943 
    ## Run 376 stress 0.08938549 
    ## ... Procrustes: rmse 0.01181611  max resid 0.03973656 
    ## Run 377 stress 0.1062597 
    ## Run 378 stress 0.08938538 
    ## ... Procrustes: rmse 0.01186597  max resid 0.0397285 
    ## Run 379 stress 0.08926085 
    ## ... Procrustes: rmse 0.0003397393  max resid 0.0009689325 
    ## ... Similar to previous best
    ## Run 380 stress 0.09130097 
    ## Run 381 stress 0.08938553 
    ## ... Procrustes: rmse 0.0117733  max resid 0.03965472 
    ## Run 382 stress 0.0950342 
    ## Run 383 stress 0.1101447 
    ## Run 384 stress 0.09130086 
    ## Run 385 stress 0.1080638 
    ## Run 386 stress 0.08938964 
    ## ... Procrustes: rmse 0.03594298  max resid 0.1182678 
    ## Run 387 stress 0.08938974 
    ## ... Procrustes: rmse 0.03592394  max resid 0.1182369 
    ## Run 388 stress 0.08926105 
    ## ... Procrustes: rmse 0.0004603097  max resid 0.001331594 
    ## ... Similar to previous best
    ## Run 389 stress 0.09503462 
    ## Run 390 stress 0.1067378 
    ## Run 391 stress 0.09087339 
    ## Run 392 stress 0.1052649 
    ## Run 393 stress 0.09018716 
    ## Run 394 stress 0.1108454 
    ## Run 395 stress 0.1052647 
    ## Run 396 stress 0.08938554 
    ## ... Procrustes: rmse 0.01179951  max resid 0.03961346 
    ## Run 397 stress 0.09088625 
    ## Run 398 stress 0.09180883 
    ## Run 399 stress 0.08946336 
    ## ... Procrustes: rmse 0.03749523  max resid 0.1177301 
    ## Run 400 stress 0.0894665 
    ## ... Procrustes: rmse 0.03323169  max resid 0.1163619 
    ## Run 401 stress 0.08938539 
    ## ... Procrustes: rmse 0.0118264  max resid 0.03973063 
    ## Run 402 stress 0.1061302 
    ## Run 403 stress 0.09592167 
    ## Run 404 stress 0.1056901 
    ## Run 405 stress 0.09612828 
    ## Run 406 stress 0.1091447 
    ## Run 407 stress 0.0892607 
    ## ... Procrustes: rmse 0.0001505588  max resid 0.0005084119 
    ## ... Similar to previous best
    ## Run 408 stress 0.1056894 
    ## Run 409 stress 0.1074432 
    ## Run 410 stress 0.09039147 
    ## Run 411 stress 0.08938963 
    ## ... Procrustes: rmse 0.03597674  max resid 0.1183244 
    ## Run 412 stress 0.09503456 
    ## Run 413 stress 0.09099594 
    ## Run 414 stress 0.08938539 
    ## ... Procrustes: rmse 0.01185626  max resid 0.03976058 
    ## Run 415 stress 0.08926073 
    ## ... Procrustes: rmse 0.0001743289  max resid 0.0005796411 
    ## ... Similar to previous best
    ## Run 416 stress 0.09130092 
    ## Run 417 stress 0.09592133 
    ## Run 418 stress 0.08926072 
    ## ... Procrustes: rmse 0.0001632576  max resid 0.0005361514 
    ## ... Similar to previous best
    ## Run 419 stress 0.09109109 
    ## Run 420 stress 0.1088395 
    ## Run 421 stress 0.09503438 
    ## Run 422 stress 0.0903913 
    ## Run 423 stress 0.08946337 
    ## ... Procrustes: rmse 0.03749341  max resid 0.1177225 
    ## Run 424 stress 0.09018707 
    ## Run 425 stress 0.08946676 
    ## ... Procrustes: rmse 0.03320981  max resid 0.1163276 
    ## Run 426 stress 0.08926068 
    ## ... Procrustes: rmse 0.000121459  max resid 0.0003982415 
    ## ... Similar to previous best
    ## Run 427 stress 0.08926069 
    ## ... Procrustes: rmse 0.0001273798  max resid 0.0004235965 
    ## ... Similar to previous best
    ## Run 428 stress 0.09044611 
    ## Run 429 stress 0.08938966 
    ## ... Procrustes: rmse 0.03597968  max resid 0.1183287 
    ## Run 430 stress 0.09130092 
    ## Run 431 stress 0.1064191 
    ## Run 432 stress 0.09087346 
    ## Run 433 stress 0.08938538 
    ## ... Procrustes: rmse 0.01185666  max resid 0.03975257 
    ## Run 434 stress 0.1052651 
    ## Run 435 stress 0.09044622 
    ## Run 436 stress 0.08946652 
    ## ... Procrustes: rmse 0.03322238  max resid 0.1163422 
    ## Run 437 stress 0.08946328 
    ## ... Procrustes: rmse 0.03753159  max resid 0.117789 
    ## Run 438 stress 0.1092129 
    ## Run 439 stress 0.1074433 
    ## Run 440 stress 0.08926087 
    ## ... Procrustes: rmse 0.000283962  max resid 0.000905371 
    ## ... Similar to previous best
    ## Run 441 stress 0.08938962 
    ## ... Procrustes: rmse 0.03595057  max resid 0.1182789 
    ## Run 442 stress 0.1129739 
    ## Run 443 stress 0.1085124 
    ## Run 444 stress 0.1086568 
    ## Run 445 stress 0.09021165 
    ## Run 446 stress 0.1062568 
    ## Run 447 stress 0.09088604 
    ## Run 448 stress 0.08946662 
    ## ... Procrustes: rmse 0.0331952  max resid 0.1163007 
    ## Run 449 stress 0.1066197 
    ## Run 450 stress 0.1061304 
    ## Run 451 stress 0.09594368 
    ## Run 452 stress 0.111895 
    ## Run 453 stress 0.08926074 
    ## ... Procrustes: rmse 0.0002063285  max resid 0.0005735585 
    ## ... Similar to previous best
    ## Run 454 stress 0.09039132 
    ## Run 455 stress 0.08946652 
    ## ... Procrustes: rmse 0.03324515  max resid 0.1163942 
    ## Run 456 stress 0.09018707 
    ## Run 457 stress 0.09503437 
    ## Run 458 stress 0.08938968 
    ## ... Procrustes: rmse 0.03593432  max resid 0.1182551 
    ## Run 459 stress 0.08946334 
    ## ... Procrustes: rmse 0.03749859  max resid 0.1177361 
    ## Run 460 stress 0.08938968 
    ## ... Procrustes: rmse 0.03593344  max resid 0.1182553 
    ## Run 461 stress 0.08938542 
    ## ... Procrustes: rmse 0.01182249  max resid 0.03968089 
    ## Run 462 stress 0.1075421 
    ## Run 463 stress 0.08938962 
    ## ... Procrustes: rmse 0.03595024  max resid 0.1182801 
    ## Run 464 stress 0.09130095 
    ## Run 465 stress 0.09087344 
    ## Run 466 stress 0.09039145 
    ## Run 467 stress 0.1074434 
    ## Run 468 stress 0.1074432 
    ## Run 469 stress 0.1079014 
    ## Run 470 stress 0.08938546 
    ## ... Procrustes: rmse 0.01191052  max resid 0.03973984 
    ## Run 471 stress 0.1074432 
    ## Run 472 stress 0.08938988 
    ## ... Procrustes: rmse 0.03591198  max resid 0.11822 
    ## Run 473 stress 0.1087571 
    ## Run 474 stress 0.08926086 
    ## ... Procrustes: rmse 0.0002715052  max resid 0.0009033109 
    ## ... Similar to previous best
    ## Run 475 stress 0.08946651 
    ## ... Procrustes: rmse 0.03322458  max resid 0.1163407 
    ## Run 476 stress 0.08938538 
    ## ... Procrustes: rmse 0.0118518  max resid 0.03972576 
    ## Run 477 stress 0.08926075 
    ## ... Procrustes: rmse 0.0001850538  max resid 0.0006239621 
    ## ... Similar to previous best
    ## Run 478 stress 0.09178296 
    ## Run 479 stress 0.08938961 
    ## ... Procrustes: rmse 0.0359654  max resid 0.1182961 
    ## Run 480 stress 0.0892609 
    ## ... Procrustes: rmse 0.0003087848  max resid 0.0009947839 
    ## ... Similar to previous best
    ## Run 481 stress 0.1065384 
    ## Run 482 stress 0.09503446 
    ## Run 483 stress 0.09262352 
    ## Run 484 stress 0.1061304 
    ## Run 485 stress 0.1056898 
    ## Run 486 stress 0.08938554 
    ## ... Procrustes: rmse 0.01177412  max resid 0.03962605 
    ## Run 487 stress 0.08946656 
    ## ... Procrustes: rmse 0.0332049  max resid 0.1163148 
    ## Run 488 stress 0.09503418 
    ## Run 489 stress 0.08926069 
    ## ... Procrustes: rmse 0.0002039842  max resid 0.0005752898 
    ## ... Similar to previous best
    ## Run 490 stress 0.1096137 
    ## Run 491 stress 0.1064426 
    ## Run 492 stress 0.0894634 
    ## ... Procrustes: rmse 0.03749034  max resid 0.1177282 
    ## Run 493 stress 0.1067592 
    ## Run 494 stress 0.1056904 
    ## Run 495 stress 0.1086563 
    ## Run 496 stress 0.1091452 
    ## Run 497 stress 0.1087866 
    ## Run 498 stress 0.111137 
    ## Run 499 stress 0.08938961 
    ## ... Procrustes: rmse 0.03596776  max resid 0.1183064 
    ## Run 500 stress 0.08926104 
    ## ... Procrustes: rmse 0.0003594599  max resid 0.001153156 
    ## ... Similar to previous best
    ## *** Best solution repeated 25 times

``` r
round(SD_beta_NMDS$stress, digits = 2)
```

    ## [1] 0.09

``` r
SD_beta_rep_NMDS <- metaMDS(SD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.3323292 
    ## Run 1 stress 0.3374587 
    ## Run 2 stress 0.3304626 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1868569  max resid 0.2825389 
    ## Run 3 stress 0.3293322 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1772297  max resid 0.377294 
    ## Run 4 stress 0.3471919 
    ## Run 5 stress 0.3495262 
    ## Run 6 stress 0.3365812 
    ## Run 7 stress 0.3514166 
    ## Run 8 stress 0.3282563 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1978292  max resid 0.3270921 
    ## Run 9 stress 0.3433329 
    ## Run 10 stress 0.3377184 
    ## Run 11 stress 0.3357619 
    ## Run 12 stress 0.3507661 
    ## Run 13 stress 0.3462721 
    ## Run 14 stress 0.3302371 
    ## Run 15 stress 0.3418379 
    ## Run 16 stress 0.3496452 
    ## Run 17 stress 0.3390302 
    ## Run 18 stress 0.3283939 
    ## ... Procrustes: rmse 0.1820434  max resid 0.3741733 
    ## Run 19 stress 0.3306595 
    ## Run 20 stress 0.3311174 
    ## Run 21 stress 0.328576 
    ## ... Procrustes: rmse 0.1888224  max resid 0.3546306 
    ## Run 22 stress 0.3318327 
    ## Run 23 stress 0.3335063 
    ## Run 24 stress 0.3347502 
    ## Run 25 stress 0.3255379 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1899395  max resid 0.3376353 
    ## Run 26 stress 0.3346995 
    ## Run 27 stress 0.3446364 
    ## Run 28 stress 0.3394767 
    ## Run 29 stress 0.3453831 
    ## Run 30 stress 0.3359016 
    ## Run 31 stress 0.3399736 
    ## Run 32 stress 0.3292382 
    ## Run 33 stress 0.3363418 
    ## Run 34 stress 0.3341911 
    ## Run 35 stress 0.3401912 
    ## Run 36 stress 0.3320451 
    ## Run 37 stress 0.3406347 
    ## Run 38 stress 0.3379009 
    ## Run 39 stress 0.3467312 
    ## Run 40 stress 0.335347 
    ## Run 41 stress 0.3278788 
    ## Run 42 stress 0.330939 
    ## Run 43 stress 0.345339 
    ## Run 44 stress 0.3324765 
    ## Run 45 stress 0.3371328 
    ## Run 46 stress 0.3373262 
    ## Run 47 stress 0.3469454 
    ## Run 48 stress 0.3525745 
    ## Run 49 stress 0.3313391 
    ## Run 50 stress 0.3329789 
    ## Run 51 stress 0.3315559 
    ## Run 52 stress 0.3362819 
    ## Run 53 stress 0.3349173 
    ## Run 54 stress 0.3367645 
    ## Run 55 stress 0.3420115 
    ## Run 56 stress 0.3332133 
    ## Run 57 stress 0.339724 
    ## Run 58 stress 0.3443216 
    ## Run 59 stress 0.3426353 
    ## Run 60 stress 0.3414892 
    ## Run 61 stress 0.3428015 
    ## Run 62 stress 0.3517217 
    ## Run 63 stress 0.3419504 
    ## Run 64 stress 0.3431947 
    ## Run 65 stress 0.3325768 
    ## Run 66 stress 0.3473107 
    ## Run 67 stress 0.3314542 
    ## Run 68 stress 0.3297643 
    ## Run 69 stress 0.335273 
    ## Run 70 stress 0.3308598 
    ## Run 71 stress 0.3324887 
    ## Run 72 stress 0.3405252 
    ## Run 73 stress 0.3285186 
    ## Run 74 stress 0.3314551 
    ## Run 75 stress 0.3265131 
    ## Run 76 stress 0.3412994 
    ## Run 77 stress 0.3473815 
    ## Run 78 stress 0.3373812 
    ## Run 79 stress 0.3446912 
    ## Run 80 stress 0.3335991 
    ## Run 81 stress 0.3319928 
    ## Run 82 stress 0.3310225 
    ## Run 83 stress 0.3342829 
    ## Run 84 stress 0.3424738 
    ## Run 85 stress 0.3431198 
    ## Run 86 stress 0.3359715 
    ## Run 87 stress 0.3408717 
    ## Run 88 stress 0.3303487 
    ## Run 89 stress 0.3292254 
    ## Run 90 stress 0.3502864 
    ## Run 91 stress 0.3386647 
    ## Run 92 stress 0.3347975 
    ## Run 93 stress 0.3345973 
    ## Run 94 stress 0.3286777 
    ## Run 95 stress 0.3507049 
    ## Run 96 stress 0.3268172 
    ## Run 97 stress 0.3374519 
    ## Run 98 stress 0.3393357 
    ## Run 99 stress 0.3463873 
    ## Run 100 stress 0.3357862 
    ## Run 101 stress 0.3393089 
    ## Run 102 stress 0.3394224 
    ## Run 103 stress 0.3297138 
    ## Run 104 stress 0.3342721 
    ## Run 105 stress 0.3509975 
    ## Run 106 stress 0.339178 
    ## Run 107 stress 0.3410785 
    ## Run 108 stress 0.3343129 
    ## Run 109 stress 0.3336061 
    ## Run 110 stress 0.340256 
    ## Run 111 stress 0.3333631 
    ## Run 112 stress 0.3312211 
    ## Run 113 stress 0.3325757 
    ## Run 114 stress 0.3358991 
    ## Run 115 stress 0.3345698 
    ## Run 116 stress 0.3340912 
    ## Run 117 stress 0.3348723 
    ## Run 118 stress 0.3464642 
    ## Run 119 stress 0.3394255 
    ## Run 120 stress 0.3297213 
    ## Run 121 stress 0.3486019 
    ## Run 122 stress 0.330901 
    ## Run 123 stress 0.3431353 
    ## Run 124 stress 0.3363719 
    ## Run 125 stress 0.3351771 
    ## Run 126 stress 0.3395945 
    ## Run 127 stress 0.3391387 
    ## Run 128 stress 0.338897 
    ## Run 129 stress 0.3467168 
    ## Run 130 stress 0.3472192 
    ## Run 131 stress 0.3449669 
    ## Run 132 stress 0.3460263 
    ## Run 133 stress 0.3424119 
    ## Run 134 stress 0.3311428 
    ## Run 135 stress 0.3329756 
    ## Run 136 stress 0.3361361 
    ## Run 137 stress 0.3446763 
    ## Run 138 stress 0.3361701 
    ## Run 139 stress 0.3310015 
    ## Run 140 stress 0.343922 
    ## Run 141 stress 0.3289562 
    ## Run 142 stress 0.3319723 
    ## Run 143 stress 0.3289073 
    ## Run 144 stress 0.3365578 
    ## Run 145 stress 0.3461835 
    ## Run 146 stress 0.3501881 
    ## Run 147 stress 0.3463 
    ## Run 148 stress 0.3381474 
    ## Run 149 stress 0.330178 
    ## Run 150 stress 0.3307928 
    ## Run 151 stress 0.3461947 
    ## Run 152 stress 0.3358273 
    ## Run 153 stress 0.3431965 
    ## Run 154 stress 0.329298 
    ## Run 155 stress 0.3440004 
    ## Run 156 stress 0.3345873 
    ## Run 157 stress 0.3410813 
    ## Run 158 stress 0.3308137 
    ## Run 159 stress 0.3423264 
    ## Run 160 stress 0.3449852 
    ## Run 161 stress 0.338846 
    ## Run 162 stress 0.3339096 
    ## Run 163 stress 0.3325641 
    ## Run 164 stress 0.3364045 
    ## Run 165 stress 0.3302825 
    ## Run 166 stress 0.3363466 
    ## Run 167 stress 0.3336599 
    ## Run 168 stress 0.3410409 
    ## Run 169 stress 0.332538 
    ## Run 170 stress 0.3269533 
    ## Run 171 stress 0.3271227 
    ## Run 172 stress 0.3306337 
    ## Run 173 stress 0.3311542 
    ## Run 174 stress 0.3300904 
    ## Run 175 stress 0.3360437 
    ## Run 176 stress 0.3330155 
    ## Run 177 stress 0.3352303 
    ## Run 178 stress 0.3341488 
    ## Run 179 stress 0.3395215 
    ## Run 180 stress 0.3560964 
    ## Run 181 stress 0.3266085 
    ## Run 182 stress 0.3320518 
    ## Run 183 stress 0.3406797 
    ## Run 184 stress 0.3385804 
    ## Run 185 stress 0.3373552 
    ## Run 186 stress 0.3383297 
    ## Run 187 stress 0.3320551 
    ## Run 188 stress 0.3304925 
    ## Run 189 stress 0.3434719 
    ## Run 190 stress 0.3260641 
    ## Run 191 stress 0.3331841 
    ## Run 192 stress 0.3352353 
    ## Run 193 stress 0.3288735 
    ## Run 194 stress 0.335085 
    ## Run 195 stress 0.3327203 
    ## Run 196 stress 0.3309972 
    ## Run 197 stress 0.3448301 
    ## Run 198 stress 0.3372636 
    ## Run 199 stress 0.3279755 
    ## Run 200 stress 0.3334712 
    ## Run 201 stress 0.3336889 
    ## Run 202 stress 0.3516778 
    ## Run 203 stress 0.34253 
    ## Run 204 stress 0.3418001 
    ## Run 205 stress 0.330153 
    ## Run 206 stress 0.3367964 
    ## Run 207 stress 0.3496774 
    ## Run 208 stress 0.3412266 
    ## Run 209 stress 0.3406059 
    ## Run 210 stress 0.350247 
    ## Run 211 stress 0.3328689 
    ## Run 212 stress 0.33295 
    ## Run 213 stress 0.3402065 
    ## Run 214 stress 0.342235 
    ## Run 215 stress 0.3377962 
    ## Run 216 stress 0.3496966 
    ## Run 217 stress 0.3379217 
    ## Run 218 stress 0.3313395 
    ## Run 219 stress 0.3412316 
    ## Run 220 stress 0.3513894 
    ## Run 221 stress 0.327128 
    ## Run 222 stress 0.339779 
    ## Run 223 stress 0.3432489 
    ## Run 224 stress 0.3384842 
    ## Run 225 stress 0.3284672 
    ## Run 226 stress 0.3415005 
    ## Run 227 stress 0.3347693 
    ## Run 228 stress 0.3369383 
    ## Run 229 stress 0.3424779 
    ## Run 230 stress 0.3340221 
    ## Run 231 stress 0.3470759 
    ## Run 232 stress 0.3322799 
    ## Run 233 stress 0.3347234 
    ## Run 234 stress 0.3426058 
    ## Run 235 stress 0.3445125 
    ## Run 236 stress 0.3325376 
    ## Run 237 stress 0.343451 
    ## Run 238 stress 0.33064 
    ## Run 239 stress 0.3427548 
    ## Run 240 stress 0.3482146 
    ## Run 241 stress 0.3328357 
    ## Run 242 stress 0.3303179 
    ## Run 243 stress 0.3313032 
    ## Run 244 stress 0.3592894 
    ## Run 245 stress 0.3452798 
    ## Run 246 stress 0.3439624 
    ## Run 247 stress 0.3324203 
    ## Run 248 stress 0.3335273 
    ## Run 249 stress 0.3415459 
    ## Run 250 stress 0.3684501 
    ## Run 251 stress 0.3320626 
    ## Run 252 stress 0.3295891 
    ## Run 253 stress 0.3455099 
    ## Run 254 stress 0.3360222 
    ## Run 255 stress 0.3296237 
    ## Run 256 stress 0.3500379 
    ## Run 257 stress 0.3480235 
    ## Run 258 stress 0.3356225 
    ## Run 259 stress 0.3320725 
    ## Run 260 stress 0.3343834 
    ## Run 261 stress 0.3486721 
    ## Run 262 stress 0.3302675 
    ## Run 263 stress 0.3423278 
    ## Run 264 stress 0.3322627 
    ## Run 265 stress 0.3340726 
    ## Run 266 stress 0.3326724 
    ## Run 267 stress 0.3361804 
    ## Run 268 stress 0.3319649 
    ## Run 269 stress 0.3807933 
    ## Run 270 stress 0.3392032 
    ## Run 271 stress 0.3328769 
    ## Run 272 stress 0.3335199 
    ## Run 273 stress 0.3454595 
    ## Run 274 stress 0.3365733 
    ## Run 275 stress 0.3499939 
    ## Run 276 stress 0.3353411 
    ## Run 277 stress 0.333617 
    ## Run 278 stress 0.3466297 
    ## Run 279 stress 0.3462835 
    ## Run 280 stress 0.3443423 
    ## Run 281 stress 0.3325467 
    ## Run 282 stress 0.331568 
    ## Run 283 stress 0.3415195 
    ## Run 284 stress 0.3430396 
    ## Run 285 stress 0.3388116 
    ## Run 286 stress 0.3308492 
    ## Run 287 stress 0.3367184 
    ## Run 288 stress 0.3366657 
    ## Run 289 stress 0.3370952 
    ## Run 290 stress 0.3468821 
    ## Run 291 stress 0.3469965 
    ## Run 292 stress 0.328631 
    ## Run 293 stress 0.3303248 
    ## Run 294 stress 0.3411562 
    ## Run 295 stress 0.330505 
    ## Run 296 stress 0.3278904 
    ## Run 297 stress 0.341819 
    ## Run 298 stress 0.3286736 
    ## Run 299 stress 0.338591 
    ## Run 300 stress 0.3371715 
    ## Run 301 stress 0.3327832 
    ## Run 302 stress 0.3337372 
    ## Run 303 stress 0.3335982 
    ## Run 304 stress 0.335544 
    ## Run 305 stress 0.3452327 
    ## Run 306 stress 0.3505386 
    ## Run 307 stress 0.3302616 
    ## Run 308 stress 0.3372089 
    ## Run 309 stress 0.3446425 
    ## Run 310 stress 0.3332939 
    ## Run 311 stress 0.3335448 
    ## Run 312 stress 0.3401395 
    ## Run 313 stress 0.3464608 
    ## Run 314 stress 0.3491709 
    ## Run 315 stress 0.3294534 
    ## Run 316 stress 0.3287702 
    ## Run 317 stress 0.3413863 
    ## Run 318 stress 0.3361423 
    ## Run 319 stress 0.334012 
    ## Run 320 stress 0.3417003 
    ## Run 321 stress 0.3373231 
    ## Run 322 stress 0.3324813 
    ## Run 323 stress 0.3338798 
    ## Run 324 stress 0.3320837 
    ## Run 325 stress 0.3454029 
    ## Run 326 stress 0.3314261 
    ## Run 327 stress 0.3420607 
    ## Run 328 stress 0.3315601 
    ## Run 329 stress 0.3352447 
    ## Run 330 stress 0.3316166 
    ## Run 331 stress 0.332046 
    ## Run 332 stress 0.3326492 
    ## Run 333 stress 0.343118 
    ## Run 334 stress 0.33815 
    ## Run 335 stress 0.348657 
    ## Run 336 stress 0.3506133 
    ## Run 337 stress 0.3363378 
    ## Run 338 stress 0.3302756 
    ## Run 339 stress 0.329542 
    ## Run 340 stress 0.3340257 
    ## Run 341 stress 0.3421112 
    ## Run 342 stress 0.3313745 
    ## Run 343 stress 0.3359862 
    ## Run 344 stress 0.3365247 
    ## Run 345 stress 0.3375106 
    ## Run 346 stress 0.3310217 
    ## Run 347 stress 0.334621 
    ## Run 348 stress 0.3402416 
    ## Run 349 stress 0.3364105 
    ## Run 350 stress 0.3398497 
    ## Run 351 stress 0.3300458 
    ## Run 352 stress 0.3297406 
    ## Run 353 stress 0.3393518 
    ## Run 354 stress 0.3479432 
    ## Run 355 stress 0.3402516 
    ## Run 356 stress 0.3417563 
    ## Run 357 stress 0.33614 
    ## Run 358 stress 0.3314338 
    ## Run 359 stress 0.3342448 
    ## Run 360 stress 0.3270597 
    ## Run 361 stress 0.3289269 
    ## Run 362 stress 0.3432875 
    ## Run 363 stress 0.3315823 
    ## Run 364 stress 0.3291059 
    ## Run 365 stress 0.3295111 
    ## Run 366 stress 0.3268804 
    ## Run 367 stress 0.327956 
    ## Run 368 stress 0.3344971 
    ## Run 369 stress 0.3317512 
    ## Run 370 stress 0.3275414 
    ## Run 371 stress 0.3318535 
    ## Run 372 stress 0.3278498 
    ## Run 373 stress 0.3459625 
    ## Run 374 stress 0.3378761 
    ## Run 375 stress 0.3336836 
    ## Run 376 stress 0.3370201 
    ## Run 377 stress 0.3320231 
    ## Run 378 stress 0.3381702 
    ## Run 379 stress 0.3369518 
    ## Run 380 stress 0.3328147 
    ## Run 381 stress 0.3340939 
    ## Run 382 stress 0.3469133 
    ## Run 383 stress 0.328618 
    ## Run 384 stress 0.3363268 
    ## Run 385 stress 0.3345581 
    ## Run 386 stress 0.3434398 
    ## Run 387 stress 0.336933 
    ## Run 388 stress 0.337234 
    ## Run 389 stress 0.3388396 
    ## Run 390 stress 0.3390268 
    ## Run 391 stress 0.3416631 
    ## Run 392 stress 0.3408362 
    ## Run 393 stress 0.3455151 
    ## Run 394 stress 0.3307112 
    ## Run 395 stress 0.3359904 
    ## Run 396 stress 0.3395555 
    ## Run 397 stress 0.3392149 
    ## Run 398 stress 0.3295898 
    ## Run 399 stress 0.3284923 
    ## Run 400 stress 0.3389374 
    ## Run 401 stress 0.3328127 
    ## Run 402 stress 0.3324572 
    ## Run 403 stress 0.3332701 
    ## Run 404 stress 0.3391643 
    ## Run 405 stress 0.3325717 
    ## Run 406 stress 0.3294136 
    ## Run 407 stress 0.3397984 
    ## Run 408 stress 0.3389974 
    ## Run 409 stress 0.3395284 
    ## Run 410 stress 0.3390418 
    ## Run 411 stress 0.3388197 
    ## Run 412 stress 0.3349329 
    ## Run 413 stress 0.3434593 
    ## Run 414 stress 0.3271025 
    ## Run 415 stress 0.3340753 
    ## Run 416 stress 0.3356581 
    ## Run 417 stress 0.3375356 
    ## Run 418 stress 0.3491537 
    ## Run 419 stress 0.3345265 
    ## Run 420 stress 0.3364743 
    ## Run 421 stress 0.336236 
    ## Run 422 stress 0.3294119 
    ## Run 423 stress 0.3291973 
    ## Run 424 stress 0.3353341 
    ## Run 425 stress 0.3439553 
    ## Run 426 stress 0.3307457 
    ## Run 427 stress 0.3481436 
    ## Run 428 stress 0.3319824 
    ## Run 429 stress 0.3362632 
    ## Run 430 stress 0.3366428 
    ## Run 431 stress 0.3448369 
    ## Run 432 stress 0.3340627 
    ## Run 433 stress 0.3388219 
    ## Run 434 stress 0.3362473 
    ## Run 435 stress 0.3368746 
    ## Run 436 stress 0.3345569 
    ## Run 437 stress 0.3376495 
    ## Run 438 stress 0.3449376 
    ## Run 439 stress 0.3352273 
    ## Run 440 stress 0.3298814 
    ## Run 441 stress 0.3347414 
    ## Run 442 stress 0.3374325 
    ## Run 443 stress 0.3308058 
    ## Run 444 stress 0.3312351 
    ## Run 445 stress 0.3411965 
    ## Run 446 stress 0.3376636 
    ## Run 447 stress 0.3450358 
    ## Run 448 stress 0.3604492 
    ## Run 449 stress 0.3358272 
    ## Run 450 stress 0.3289625 
    ## Run 451 stress 0.3306582 
    ## Run 452 stress 0.3356172 
    ## Run 453 stress 0.3288559 
    ## Run 454 stress 0.3354663 
    ## Run 455 stress 0.3253155 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1720512  max resid 0.3571475 
    ## Run 456 stress 0.3401791 
    ## Run 457 stress 0.3479178 
    ## Run 458 stress 0.3339159 
    ## Run 459 stress 0.3383316 
    ## Run 460 stress 0.3278212 
    ## Run 461 stress 0.3333486 
    ## Run 462 stress 0.3312661 
    ## Run 463 stress 0.3301728 
    ## Run 464 stress 0.3320754 
    ## Run 465 stress 0.3448441 
    ## Run 466 stress 0.3337868 
    ## Run 467 stress 0.3304085 
    ## Run 468 stress 0.3344642 
    ## Run 469 stress 0.3323514 
    ## Run 470 stress 0.3273547 
    ## Run 471 stress 0.347632 
    ## Run 472 stress 0.3306916 
    ## Run 473 stress 0.3373246 
    ## Run 474 stress 0.3314271 
    ## Run 475 stress 0.3381743 
    ## Run 476 stress 0.3378909 
    ## Run 477 stress 0.3288642 
    ## Run 478 stress 0.3287773 
    ## Run 479 stress 0.3441666 
    ## Run 480 stress 0.333302 
    ## Run 481 stress 0.3375701 
    ## Run 482 stress 0.3402537 
    ## Run 483 stress 0.332011 
    ## Run 484 stress 0.3406407 
    ## Run 485 stress 0.3404648 
    ## Run 486 stress 0.3462942 
    ## Run 487 stress 0.3309664 
    ## Run 488 stress 0.3284106 
    ## Run 489 stress 0.3383546 
    ## Run 490 stress 0.34065 
    ## Run 491 stress 0.3489587 
    ## Run 492 stress 0.3326469 
    ## Run 493 stress 0.3355108 
    ## Run 494 stress 0.3295959 
    ## Run 495 stress 0.3354833 
    ## Run 496 stress 0.3375275 
    ## Run 497 stress 0.3465619 
    ## Run 498 stress 0.3359059 
    ## Run 499 stress 0.3406537 
    ## Run 500 stress 0.3319088 
    ## *** Best solution was not repeated -- monoMDS stopping criteria:
    ##    500: stress ratio > sratmax

``` r
SD_beta_ric_NMDS <- metaMDS(SD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.01879552 
    ## Run 1 stress 0.02492201 
    ## Run 2 stress 0.01879709 
    ## ... Procrustes: rmse 0.001644742  max resid 0.003385414 
    ## ... Similar to previous best
    ## Run 3 stress 0.01879559 
    ## ... Procrustes: rmse 0.00116931  max resid 0.002408309 
    ## ... Similar to previous best
    ## Run 4 stress 0.01879583 
    ## ... Procrustes: rmse 0.00128036  max resid 0.002636832 
    ## ... Similar to previous best
    ## Run 5 stress 0.01894987 
    ## ... Procrustes: rmse 0.009940781  max resid 0.02023923 
    ## Run 6 stress 0.0250947 
    ## Run 7 stress 0.01879582 
    ## ... Procrustes: rmse 0.0001386014  max resid 0.0002845011 
    ## ... Similar to previous best
    ## Run 8 stress 0.0249218 
    ## Run 9 stress 0.01879608 
    ## ... Procrustes: rmse 0.0002299086  max resid 0.0004718691 
    ## ... Similar to previous best
    ## Run 10 stress 0.0187974 
    ## ... Procrustes: rmse 0.001767768  max resid 0.00364036 
    ## ... Similar to previous best
    ## Run 11 stress 0.018832 
    ## ... Procrustes: rmse 0.003260717  max resid 0.00630706 
    ## ... Similar to previous best
    ## Run 12 stress 0.02009799 
    ## Run 13 stress 0.0248612 
    ## Run 14 stress 0.01894266 
    ## ... Procrustes: rmse 0.009922339  max resid 0.02036356 
    ## Run 15 stress 0.02486128 
    ## Run 16 stress 0.02492197 
    ## Run 17 stress 0.01879572 
    ## ... Procrustes: rmse 0.001219018  max resid 0.002511196 
    ## ... Similar to previous best
    ## Run 18 stress 0.357918 
    ## Run 19 stress 0.01889014 
    ## ... Procrustes: rmse 0.008205323  max resid 0.01685147 
    ## Run 20 stress 0.01879505 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003704121  max resid 0.0007646839 
    ## ... Similar to previous best
    ## Run 21 stress 0.02509467 
    ## Run 22 stress 0.02509459 
    ## Run 23 stress 0.02520897 
    ## Run 24 stress 0.0188487 
    ## ... Procrustes: rmse 0.004590973  max resid 0.009118893 
    ## ... Similar to previous best
    ## Run 25 stress 0.01879593 
    ## ... Procrustes: rmse 0.0005465223  max resid 0.001126677 
    ## ... Similar to previous best
    ## Run 26 stress 0.02486145 
    ## Run 27 stress 0.02486104 
    ## Run 28 stress 0.01889181 
    ## ... Procrustes: rmse 0.007899776  max resid 0.0162255 
    ## Run 29 stress 0.02520904 
    ## Run 30 stress 0.01883687 
    ## ... Procrustes: rmse 0.003518545  max resid 0.006841661 
    ## ... Similar to previous best
    ## Run 31 stress 0.0187959 
    ## ... Procrustes: rmse 0.0009337712  max resid 0.001920965 
    ## ... Similar to previous best
    ## Run 32 stress 0.02509442 
    ## Run 33 stress 0.02509446 
    ## Run 34 stress 0.02504771 
    ## Run 35 stress 0.01879601 
    ## ... Procrustes: rmse 0.0009600457  max resid 0.001975106 
    ## ... Similar to previous best
    ## Run 36 stress 0.0187958 
    ## ... Procrustes: rmse 0.0008939968  max resid 0.001838999 
    ## ... Similar to previous best
    ## Run 37 stress 0.01879579 
    ## ... Procrustes: rmse 0.0008905404  max resid 0.001831867 
    ## ... Similar to previous best
    ## Run 38 stress 0.0187999 
    ## ... Procrustes: rmse 0.00188291  max resid 0.0038896 
    ## ... Similar to previous best
    ## Run 39 stress 0.02520898 
    ## Run 40 stress 0.02492172 
    ## Run 41 stress 0.01879576 
    ## ... Procrustes: rmse 0.0008766955  max resid 0.001803409 
    ## ... Similar to previous best
    ## Run 42 stress 0.02509474 
    ## Run 43 stress 0.02509492 
    ## Run 44 stress 0.02492178 
    ## Run 45 stress 0.01883071 
    ## ... Procrustes: rmse 0.002703507  max resid 0.005105477 
    ## ... Similar to previous best
    ## Run 46 stress 0.02492163 
    ## Run 47 stress 0.02509474 
    ## Run 48 stress 0.01879596 
    ## ... Procrustes: rmse 0.0009514018  max resid 0.001957299 
    ## ... Similar to previous best
    ## Run 49 stress 0.02520909 
    ## Run 50 stress 0.01879579 
    ## ... Procrustes: rmse 0.0004914519  max resid 0.001013715 
    ## ... Similar to previous best
    ## Run 51 stress 0.02509458 
    ## Run 52 stress 0.01881894 
    ## ... Procrustes: rmse 0.003927501  max resid 0.008078376 
    ## ... Similar to previous best
    ## Run 53 stress 0.3149075 
    ## Run 54 stress 0.02520884 
    ## Run 55 stress 0.02509439 
    ## Run 56 stress 0.02509455 
    ## Run 57 stress 0.01880827 
    ## ... Procrustes: rmse 0.003028127  max resid 0.006234068 
    ## ... Similar to previous best
    ## Run 58 stress 0.01879999 
    ## ... Procrustes: rmse 0.001929625  max resid 0.003972383 
    ## ... Similar to previous best
    ## Run 59 stress 0.0249218 
    ## Run 60 stress 0.01906756 
    ## ... Procrustes: rmse 0.0131656  max resid 0.02695988 
    ## Run 61 stress 0.3654815 
    ## Run 62 stress 0.01884592 
    ## ... Procrustes: rmse 0.00431498  max resid 0.008531685 
    ## ... Similar to previous best
    ## Run 63 stress 0.01879874 
    ## ... Procrustes: rmse 0.001694638  max resid 0.003491877 
    ## ... Similar to previous best
    ## Run 64 stress 0.02509459 
    ## Run 65 stress 0.01879573 
    ## ... Procrustes: rmse 0.0008625328  max resid 0.001774197 
    ## ... Similar to previous best
    ## Run 66 stress 0.02492196 
    ## Run 67 stress 0.02492176 
    ## Run 68 stress 0.01885285 
    ## ... Procrustes: rmse 0.004933892  max resid 0.009838728 
    ## ... Similar to previous best
    ## Run 69 stress 0.02509476 
    ## Run 70 stress 0.02486107 
    ## Run 71 stress 0.01880013 
    ## ... Procrustes: rmse 0.001955097  max resid 0.004024823 
    ## ... Similar to previous best
    ## Run 72 stress 0.02509429 
    ## Run 73 stress 0.01879569 
    ## ... Procrustes: rmse 0.0008386262  max resid 0.001725567 
    ## ... Similar to previous best
    ## Run 74 stress 0.01880059 
    ## ... Procrustes: rmse 0.001960924  max resid 0.00404057 
    ## ... Similar to previous best
    ## Run 75 stress 0.01879594 
    ## ... Procrustes: rmse 0.0009450318  max resid 0.001943937 
    ## ... Similar to previous best
    ## Run 76 stress 0.02492214 
    ## Run 77 stress 0.01879557 
    ## ... Procrustes: rmse 0.0003950765  max resid 0.0008156293 
    ## ... Similar to previous best
    ## Run 78 stress 0.01883168 
    ## ... Procrustes: rmse 0.002864704  max resid 0.005445528 
    ## ... Similar to previous best
    ## Run 79 stress 0.01916654 
    ## ... Procrustes: rmse 0.01540682  max resid 0.03160692 
    ## Run 80 stress 0.01879562 
    ## ... Procrustes: rmse 0.0004221312  max resid 0.0008712225 
    ## ... Similar to previous best
    ## Run 81 stress 0.0250944 
    ## Run 82 stress 0.01883853 
    ## ... Procrustes: rmse 0.005345976  max resid 0.0109922 
    ## Run 83 stress 0.02509476 
    ## Run 84 stress 0.01879681 
    ## ... Procrustes: rmse 0.001241719  max resid 0.002554962 
    ## ... Similar to previous best
    ## Run 85 stress 0.02509457 
    ## Run 86 stress 0.01879554 
    ## ... Procrustes: rmse 0.0007568552  max resid 0.001557187 
    ## ... Similar to previous best
    ## Run 87 stress 0.02486133 
    ## Run 88 stress 0.01993098 
    ## Run 89 stress 0.01879597 
    ## ... Procrustes: rmse 0.0009652739  max resid 0.001985763 
    ## ... Similar to previous best
    ## Run 90 stress 0.02509477 
    ## Run 91 stress 0.01879576 
    ## ... Procrustes: rmse 0.0008714417  max resid 0.001792445 
    ## ... Similar to previous best
    ## Run 92 stress 0.02492183 
    ## Run 93 stress 0.01879571 
    ## ... Procrustes: rmse 0.0008577371  max resid 0.001764343 
    ## ... Similar to previous best
    ## Run 94 stress 0.01886652 
    ## ... Procrustes: rmse 0.005842321  max resid 0.01173801 
    ## Run 95 stress 0.01879591 
    ## ... Procrustes: rmse 0.0009298591  max resid 0.001912665 
    ## ... Similar to previous best
    ## Run 96 stress 0.02509444 
    ## Run 97 stress 0.0245686 
    ## Run 98 stress 0.0252092 
    ## Run 99 stress 0.01887751 
    ## ... Procrustes: rmse 0.005238519  max resid 0.01046085 
    ## Run 100 stress 0.0189742 
    ## ... Procrustes: rmse 0.01036656  max resid 0.02111619 
    ## Run 101 stress 0.02492164 
    ## Run 102 stress 0.0250946 
    ## Run 103 stress 0.01883764 
    ## ... Procrustes: rmse 0.003602491  max resid 0.007020766 
    ## ... Similar to previous best
    ## Run 104 stress 0.01879575 
    ## ... Procrustes: rmse 0.0004771383  max resid 0.0009841052 
    ## ... Similar to previous best
    ## Run 105 stress 0.0187958 
    ## ... Procrustes: rmse 0.0008910199  max resid 0.001832958 
    ## ... Similar to previous best
    ## Run 106 stress 0.02509442 
    ## Run 107 stress 0.01879578 
    ## ... Procrustes: rmse 0.0004928959  max resid 0.001016637 
    ## ... Similar to previous best
    ## Run 108 stress 0.02509478 
    ## Run 109 stress 0.0187958 
    ## ... Procrustes: rmse 0.0008973804  max resid 0.001845964 
    ## ... Similar to previous best
    ## Run 110 stress 0.01879537 
    ## ... Procrustes: rmse 0.0006756106  max resid 0.001389273 
    ## ... Similar to previous best
    ## Run 111 stress 0.02520885 
    ## Run 112 stress 0.01880298 
    ## ... Procrustes: rmse 0.002384858  max resid 0.004909709 
    ## ... Similar to previous best
    ## Run 113 stress 0.02509442 
    ## Run 114 stress 0.01891395 
    ## ... Procrustes: rmse 0.006663767  max resid 0.0137028 
    ## Run 115 stress 0.01884594 
    ## ... Procrustes: rmse 0.005611853  max resid 0.01153149 
    ## Run 116 stress 0.0187956 
    ## ... Procrustes: rmse 0.0008019498  max resid 0.00164951 
    ## ... Similar to previous best
    ## Run 117 stress 0.02493543 
    ## Run 118 stress 0.02509469 
    ## Run 119 stress 0.01879606 
    ## ... Procrustes: rmse 0.0009829028  max resid 0.002021836 
    ## ... Similar to previous best
    ## Run 120 stress 0.02509459 
    ## Run 121 stress 0.02509432 
    ## Run 122 stress 0.02509464 
    ## Run 123 stress 0.01879568 
    ## ... Procrustes: rmse 0.0004509924  max resid 0.0009305268 
    ## ... Similar to previous best
    ## Run 124 stress 0.02509478 
    ## Run 125 stress 0.01879632 
    ## ... Procrustes: rmse 0.001088837  max resid 0.002240191 
    ## ... Similar to previous best
    ## Run 126 stress 0.01883659 
    ## ... Procrustes: rmse 0.003486459  max resid 0.006773767 
    ## ... Similar to previous best
    ## Run 127 stress 0.0187955 
    ## ... Procrustes: rmse 0.0003602612  max resid 0.00074411 
    ## ... Similar to previous best
    ## Run 128 stress 0.02509455 
    ## Run 129 stress 0.0250947 
    ## Run 130 stress 0.01879565 
    ## ... Procrustes: rmse 0.0008294877  max resid 0.001706159 
    ## ... Similar to previous best
    ## Run 131 stress 0.01879609 
    ## ... Procrustes: rmse 0.0009990625  max resid 0.002055121 
    ## ... Similar to previous best
    ## Run 132 stress 0.01879544 
    ## ... Procrustes: rmse 0.0007075026  max resid 0.00145481 
    ## ... Similar to previous best
    ## Run 133 stress 0.02486128 
    ## Run 134 stress 0.01879617 
    ## ... Procrustes: rmse 0.001020932  max resid 0.002100073 
    ## ... Similar to previous best
    ## Run 135 stress 0.01883283 
    ## ... Procrustes: rmse 0.003009907  max resid 0.005755036 
    ## ... Similar to previous best
    ## Run 136 stress 0.01879547 
    ## ... Procrustes: rmse 0.0007385336  max resid 0.001518855 
    ## ... Similar to previous best
    ## Run 137 stress 0.02486127 
    ## Run 138 stress 0.01895437 
    ## ... Procrustes: rmse 0.009734984  max resid 0.01981644 
    ## Run 139 stress 0.01879574 
    ## ... Procrustes: rmse 0.0008691184  max resid 0.001787815 
    ## ... Similar to previous best
    ## Run 140 stress 0.02520884 
    ## Run 141 stress 0.02520889 
    ## Run 142 stress 0.01879696 
    ## ... Procrustes: rmse 0.001282261  max resid 0.002638512 
    ## ... Similar to previous best
    ## Run 143 stress 0.02509471 
    ## Run 144 stress 0.01879579 
    ## ... Procrustes: rmse 0.0005006168  max resid 0.001032509 
    ## ... Similar to previous best
    ## Run 145 stress 0.01879837 
    ## ... Procrustes: rmse 0.001621175  max resid 0.003336537 
    ## ... Similar to previous best
    ## Run 146 stress 0.02509447 
    ## Run 147 stress 0.02509467 
    ## Run 148 stress 0.02492194 
    ## Run 149 stress 0.01883461 
    ## ... Procrustes: rmse 0.003254946  max resid 0.006279026 
    ## ... Similar to previous best
    ## Run 150 stress 0.01879585 
    ## ... Procrustes: rmse 0.0005155066  max resid 0.001063129 
    ## ... Similar to previous best
    ## Run 151 stress 0.01879543 
    ## ... Procrustes: rmse 0.0007129602  max resid 0.001466255 
    ## ... Similar to previous best
    ## Run 152 stress 0.02509464 
    ## Run 153 stress 0.02492171 
    ## Run 154 stress 0.0248614 
    ## Run 155 stress 0.02486107 
    ## Run 156 stress 0.01886893 
    ## ... Procrustes: rmse 0.006919904  max resid 0.01421765 
    ## Run 157 stress 0.01879589 
    ## ... Procrustes: rmse 0.0005298318  max resid 0.001092655 
    ## ... Similar to previous best
    ## Run 158 stress 0.02509461 
    ## Run 159 stress 0.01879567 
    ## ... Procrustes: rmse 0.0004434969  max resid 0.0009151404 
    ## ... Similar to previous best
    ## Run 160 stress 0.01879573 
    ## ... Procrustes: rmse 0.0008619387  max resid 0.001773046 
    ## ... Similar to previous best
    ## Run 161 stress 0.01879556 
    ## ... Procrustes: rmse 0.0007827293  max resid 0.001609934 
    ## ... Similar to previous best
    ## Run 162 stress 0.01879578 
    ## ... Procrustes: rmse 0.0004916364  max resid 0.001014 
    ## ... Similar to previous best
    ## Run 163 stress 0.02509449 
    ## Run 164 stress 0.01879547 
    ## ... Procrustes: rmse 0.0007322075  max resid 0.001505878 
    ## ... Similar to previous best
    ## Run 165 stress 0.02016686 
    ## Run 166 stress 0.01880106 
    ## ... Procrustes: rmse 0.002022418  max resid 0.004164115 
    ## ... Similar to previous best
    ## Run 167 stress 0.02492198 
    ## Run 168 stress 0.01879578 
    ## ... Procrustes: rmse 0.0008855463  max resid 0.001821594 
    ## ... Similar to previous best
    ## Run 169 stress 0.02492185 
    ## Run 170 stress 0.02509453 
    ## Run 171 stress 0.01884449 
    ## ... Procrustes: rmse 0.005689678  max resid 0.01169678 
    ## Run 172 stress 0.02509469 
    ## Run 173 stress 0.02492171 
    ## Run 174 stress 0.01887317 
    ## ... Procrustes: rmse 0.007090255  max resid 0.01457 
    ## Run 175 stress 0.02486125 
    ## Run 176 stress 0.02492621 
    ## Run 177 stress 0.0250946 
    ## Run 178 stress 0.0187953 
    ## ... Procrustes: rmse 0.0006276527  max resid 0.001290511 
    ## ... Similar to previous best
    ## Run 179 stress 0.01879563 
    ## ... Procrustes: rmse 0.0008169227  max resid 0.001680272 
    ## ... Similar to previous best
    ## Run 180 stress 0.01879581 
    ## ... Procrustes: rmse 0.0009017502  max resid 0.001854992 
    ## ... Similar to previous best
    ## Run 181 stress 0.02509429 
    ## Run 182 stress 0.01891098 
    ## ... Procrustes: rmse 0.007889632  max resid 0.01600482 
    ## Run 183 stress 0.0187956 
    ## ... Procrustes: rmse 0.000803518  max resid 0.001652701 
    ## ... Similar to previous best
    ## Run 184 stress 0.0187958 
    ## ... Procrustes: rmse 0.0008948345  max resid 0.001840728 
    ## ... Similar to previous best
    ## Run 185 stress 0.02486128 
    ## Run 186 stress 0.01879595 
    ## ... Procrustes: rmse 0.0009455911  max resid 0.001945038 
    ## ... Similar to previous best
    ## Run 187 stress 0.02509456 
    ## Run 188 stress 0.0189467 
    ## ... Procrustes: rmse 0.009873228  max resid 0.02025843 
    ## Run 189 stress 0.01879579 
    ## ... Procrustes: rmse 0.0008906724  max resid 0.001832208 
    ## ... Similar to previous best
    ## Run 190 stress 0.01881772 
    ## ... Procrustes: rmse 0.00354099  max resid 0.007290258 
    ## ... Similar to previous best
    ## Run 191 stress 0.01880732 
    ## ... Procrustes: rmse 0.002918737  max resid 0.00600922 
    ## ... Similar to previous best
    ## Run 192 stress 0.02492185 
    ## Run 193 stress 0.02520904 
    ## Run 194 stress 0.01879563 
    ## ... Procrustes: rmse 0.0008197948  max resid 0.001686219 
    ## ... Similar to previous best
    ## Run 195 stress 0.02486118 
    ## Run 196 stress 0.0252091 
    ## Run 197 stress 0.02509441 
    ## Run 198 stress 0.02509442 
    ## Run 199 stress 0.01879558 
    ## ... Procrustes: rmse 0.0007934449  max resid 0.001631984 
    ## ... Similar to previous best
    ## Run 200 stress 0.02486116 
    ## Run 201 stress 0.02520875 
    ## Run 202 stress 0.01879562 
    ## ... Procrustes: rmse 0.0008142878  max resid 0.001674875 
    ## ... Similar to previous best
    ## Run 203 stress 0.01879895 
    ## ... Procrustes: rmse 0.001734984  max resid 0.003570998 
    ## ... Similar to previous best
    ## Run 204 stress 0.01879527 
    ## ... Procrustes: rmse 0.0006064718  max resid 0.001246805 
    ## ... Similar to previous best
    ## Run 205 stress 0.01879591 
    ## ... Procrustes: rmse 0.0009423834  max resid 0.001938656 
    ## ... Similar to previous best
    ## Run 206 stress 0.01879568 
    ## ... Procrustes: rmse 0.0008429118  max resid 0.00173379 
    ## ... Similar to previous best
    ## Run 207 stress 0.02492201 
    ## Run 208 stress 0.02509442 
    ## Run 209 stress 0.01884005 
    ## ... Procrustes: rmse 0.004631868  max resid 0.009514212 
    ## ... Similar to previous best
    ## Run 210 stress 0.02492189 
    ## Run 211 stress 0.01882307 
    ## ... Procrustes: rmse 0.004305603  max resid 0.008855514 
    ## ... Similar to previous best
    ## Run 212 stress 0.01879567 
    ## ... Procrustes: rmse 0.0004430408  max resid 0.000914213 
    ## ... Similar to previous best
    ## Run 213 stress 0.02520909 
    ## Run 214 stress 0.02492193 
    ## Run 215 stress 0.01882053 
    ## ... Procrustes: rmse 0.004075926  max resid 0.008386081 
    ## ... Similar to previous best
    ## Run 216 stress 0.02509456 
    ## Run 217 stress 0.02509456 
    ## Run 218 stress 0.01879566 
    ## ... Procrustes: rmse 0.0008262545  max resid 0.001699435 
    ## ... Similar to previous best
    ## Run 219 stress 0.02520878 
    ## Run 220 stress 0.02492184 
    ## Run 221 stress 0.01880094 
    ## ... Procrustes: rmse 0.0020888  max resid 0.00430027 
    ## ... Similar to previous best
    ## Run 222 stress 0.01880058 
    ## ... Procrustes: rmse 0.001888987  max resid 0.003888673 
    ## ... Similar to previous best
    ## Run 223 stress 0.02486102 
    ## Run 224 stress 0.01923542 
    ## ... Procrustes: rmse 0.01607603  max resid 0.0330584 
    ## Run 225 stress 0.01879594 
    ## ... Procrustes: rmse 0.0005513825  max resid 0.001136659 
    ## ... Similar to previous best
    ## Run 226 stress 0.02493211 
    ## Run 227 stress 0.01879895 
    ## ... Procrustes: rmse 0.001617124  max resid 0.003330469 
    ## ... Similar to previous best
    ## Run 228 stress 0.01884712 
    ## ... Procrustes: rmse 0.005473753  max resid 0.01125373 
    ## Run 229 stress 0.02509449 
    ## Run 230 stress 0.01890171 
    ## ... Procrustes: rmse 0.008293026  max resid 0.01702989 
    ## Run 231 stress 0.0187956 
    ## ... Procrustes: rmse 0.000802516  max resid 0.001650632 
    ## ... Similar to previous best
    ## Run 232 stress 0.01879764 
    ## ... Procrustes: rmse 0.001436431  max resid 0.002955075 
    ## ... Similar to previous best
    ## Run 233 stress 0.02509453 
    ## Run 234 stress 0.01879591 
    ## ... Procrustes: rmse 0.0009417948  max resid 0.001937433 
    ## ... Similar to previous best
    ## Run 235 stress 0.01879928 
    ## ... Procrustes: rmse 0.001735041  max resid 0.003577693 
    ## ... Similar to previous best
    ## Run 236 stress 0.01879602 
    ## ... Procrustes: rmse 0.0009635997  max resid 0.001982431 
    ## ... Similar to previous best
    ## Run 237 stress 0.02509451 
    ## Run 238 stress 0.01879854 
    ## ... Procrustes: rmse 0.001634503  max resid 0.00336606 
    ## ... Similar to previous best
    ## Run 239 stress 0.02492196 
    ## Run 240 stress 0.01879581 
    ## ... Procrustes: rmse 0.000871073  max resid 0.001791491 
    ## ... Similar to previous best
    ## Run 241 stress 0.01879598 
    ## ... Procrustes: rmse 0.0009578639  max resid 0.001970337 
    ## ... Similar to previous best
    ## Run 242 stress 0.01880055 
    ## ... Procrustes: rmse 0.002024057  max resid 0.00416696 
    ## ... Similar to previous best
    ## Run 243 stress 0.01879541 
    ## ... Procrustes: rmse 0.0007034748  max resid 0.001446745 
    ## ... Similar to previous best
    ## Run 244 stress 0.02486121 
    ## Run 245 stress 0.02486133 
    ## Run 246 stress 0.01879598 
    ## ... Procrustes: rmse 0.0009689969  max resid 0.00199345 
    ## ... Similar to previous best
    ## Run 247 stress 0.02492192 
    ## Run 248 stress 0.02509457 
    ## Run 249 stress 0.01879588 
    ## ... Procrustes: rmse 0.0009292726  max resid 0.001911678 
    ## ... Similar to previous best
    ## Run 250 stress 0.01879565 
    ## ... Procrustes: rmse 0.0008260896  max resid 0.001699135 
    ## ... Similar to previous best
    ## Run 251 stress 0.01879588 
    ## ... Procrustes: rmse 0.0005250287  max resid 0.001082657 
    ## ... Similar to previous best
    ## Run 252 stress 0.02486131 
    ## Run 253 stress 0.01879602 
    ## ... Procrustes: rmse 0.000971828  max resid 0.001999354 
    ## ... Similar to previous best
    ## Run 254 stress 0.01883894 
    ## ... Procrustes: rmse 0.005343041  max resid 0.01098767 
    ## Run 255 stress 0.01879726 
    ## ... Procrustes: rmse 0.00136458  max resid 0.002807947 
    ## ... Similar to previous best
    ## Run 256 stress 0.02492184 
    ## Run 257 stress 0.01891929 
    ## ... Procrustes: rmse 0.008387354  max resid 0.01702772 
    ## Run 258 stress 0.0249216 
    ## Run 259 stress 0.0250944 
    ## Run 260 stress 0.01884848 
    ## ... Procrustes: rmse 0.004161808  max resid 0.008200679 
    ## ... Similar to previous best
    ## Run 261 stress 0.02509451 
    ## Run 262 stress 0.01925356 
    ## ... Procrustes: rmse 0.01579622  max resid 0.03249437 
    ## Run 263 stress 0.01879582 
    ## ... Procrustes: rmse 0.0009052964  max resid 0.001862279 
    ## ... Similar to previous best
    ## Run 264 stress 0.01879602 
    ## ... Procrustes: rmse 0.0009711896  max resid 0.001997736 
    ## ... Similar to previous best
    ## Run 265 stress 0.02509479 
    ## Run 266 stress 0.02509469 
    ## Run 267 stress 0.01879602 
    ## ... Procrustes: rmse 0.0009610143  max resid 0.001977103 
    ## ... Similar to previous best
    ## Run 268 stress 0.02486117 
    ## Run 269 stress 0.02509456 
    ## Run 270 stress 0.02509454 
    ## Run 271 stress 0.02492182 
    ## Run 272 stress 0.02509457 
    ## Run 273 stress 0.0249217 
    ## Run 274 stress 0.02520871 
    ## Run 275 stress 0.02509449 
    ## Run 276 stress 0.02509435 
    ## Run 277 stress 0.01879601 
    ## ... Procrustes: rmse 0.0009748313  max resid 0.002005318 
    ## ... Similar to previous best
    ## Run 278 stress 0.0187955 
    ## ... Procrustes: rmse 0.0007545634  max resid 0.001551874 
    ## ... Similar to previous best
    ## Run 279 stress 0.01879551 
    ## ... Procrustes: rmse 0.000367728  max resid 0.0007594315 
    ## ... Similar to previous best
    ## Run 280 stress 0.01879801 
    ## ... Procrustes: rmse 0.001538238  max resid 0.003161051 
    ## ... Similar to previous best
    ## Run 281 stress 0.01879599 
    ## ... Procrustes: rmse 0.0009604436  max resid 0.001975917 
    ## ... Similar to previous best
    ## Run 282 stress 0.2241413 
    ## Run 283 stress 0.01879558 
    ## ... Procrustes: rmse 0.000396124  max resid 0.000817486 
    ## ... Similar to previous best
    ## Run 284 stress 0.01879585 
    ## ... Procrustes: rmse 0.0009174928  max resid 0.001887395 
    ## ... Similar to previous best
    ## Run 285 stress 0.01879574 
    ## ... Procrustes: rmse 0.0008672278  max resid 0.001783905 
    ## ... Similar to previous best
    ## Run 286 stress 0.02492158 
    ## Run 287 stress 0.01879555 
    ## ... Procrustes: rmse 0.0007761837  max resid 0.001596412 
    ## ... Similar to previous best
    ## Run 288 stress 0.01894114 
    ## ... Procrustes: rmse 0.009261348  max resid 0.01883866 
    ## Run 289 stress 0.0248612 
    ## Run 290 stress 0.01879584 
    ## ... Procrustes: rmse 0.0005115889  max resid 0.001055047 
    ## ... Similar to previous best
    ## Run 291 stress 0.0250947 
    ## Run 292 stress 0.02486106 
    ## Run 293 stress 0.0194599 
    ## Run 294 stress 0.02509476 
    ## Run 295 stress 0.0188839 
    ## ... Procrustes: rmse 0.007580113  max resid 0.01557048 
    ## Run 296 stress 0.02486126 
    ## Run 297 stress 0.01879563 
    ## ... Procrustes: rmse 0.0008185519  max resid 0.00168368 
    ## ... Similar to previous best
    ## Run 298 stress 0.01879604 
    ## ... Procrustes: rmse 0.0009845472  max resid 0.002025527 
    ## ... Similar to previous best
    ## Run 299 stress 0.01879551 
    ## ... Procrustes: rmse 0.0007585204  max resid 0.00156004 
    ## ... Similar to previous best
    ## Run 300 stress 0.01879563 
    ## ... Procrustes: rmse 0.0008194581  max resid 0.001685549 
    ## ... Similar to previous best
    ## Run 301 stress 0.02520895 
    ## Run 302 stress 0.01881551 
    ## ... Procrustes: rmse 0.003253375  max resid 0.006690327 
    ## ... Similar to previous best
    ## Run 303 stress 0.01879559 
    ## ... Procrustes: rmse 0.0004028476  max resid 0.0008315058 
    ## ... Similar to previous best
    ## Run 304 stress 0.02486114 
    ## Run 305 stress 0.02520881 
    ## Run 306 stress 0.01891682 
    ## ... Procrustes: rmse 0.008857195  max resid 0.01818193 
    ## Run 307 stress 0.01879959 
    ## ... Procrustes: rmse 0.001831395  max resid 0.00378002 
    ## ... Similar to previous best
    ## Run 308 stress 0.01879592 
    ## ... Procrustes: rmse 0.0009442922  max resid 0.001942593 
    ## ... Similar to previous best
    ## Run 309 stress 0.02520891 
    ## Run 310 stress 0.01879592 
    ## ... Procrustes: rmse 0.0005513993  max resid 0.001136774 
    ## ... Similar to previous best
    ## Run 311 stress 0.01879576 
    ## ... Procrustes: rmse 0.0004779713  max resid 0.0009863629 
    ## ... Similar to previous best
    ## Run 312 stress 0.02509445 
    ## Run 313 stress 0.01883127 
    ## ... Procrustes: rmse 0.00278401  max resid 0.005275774 
    ## ... Similar to previous best
    ## Run 314 stress 0.02486123 
    ## Run 315 stress 0.01884169 
    ## ... Procrustes: rmse 0.005528423  max resid 0.01136711 
    ## Run 316 stress 0.02492184 
    ## Run 317 stress 0.01889611 
    ## ... Procrustes: rmse 0.007325285  max resid 0.01482856 
    ## Run 318 stress 0.02492347 
    ## Run 319 stress 0.01879551 
    ## ... Procrustes: rmse 0.0007573259  max resid 0.00155756 
    ## ... Similar to previous best
    ## Run 320 stress 0.02486132 
    ## Run 321 stress 0.01879564 
    ## ... Procrustes: rmse 0.000822473  max resid 0.001691756 
    ## ... Similar to previous best
    ## Run 322 stress 0.01881748 
    ## ... Procrustes: rmse 0.003889041  max resid 0.008002347 
    ## ... Similar to previous best
    ## Run 323 stress 0.02509438 
    ## Run 324 stress 0.02520887 
    ## Run 325 stress 0.01879562 
    ## ... Procrustes: rmse 0.0008095861  max resid 0.001665133 
    ## ... Similar to previous best
    ## Run 326 stress 0.0248612 
    ## Run 327 stress 0.01879588 
    ## ... Procrustes: rmse 0.0009302211  max resid 0.001913593 
    ## ... Similar to previous best
    ## Run 328 stress 0.0187954 
    ## ... Procrustes: rmse 0.0006949648  max resid 0.001429142 
    ## ... Similar to previous best
    ## Run 329 stress 0.01879607 
    ## ... Procrustes: rmse 0.0009914647  max resid 0.002039788 
    ## ... Similar to previous best
    ## Run 330 stress 0.0189599 
    ## ... Procrustes: rmse 0.009925102  max resid 0.02020711 
    ## Run 331 stress 0.0187983 
    ## ... Procrustes: rmse 0.001445302  max resid 0.002972254 
    ## ... Similar to previous best
    ## Run 332 stress 0.018829 
    ## ... Procrustes: rmse 0.002459015  max resid 0.004639894 
    ## ... Similar to previous best
    ## Run 333 stress 0.01883196 
    ## ... Procrustes: rmse 0.002900629  max resid 0.005522926 
    ## ... Similar to previous best
    ## Run 334 stress 0.01879571 
    ## ... Procrustes: rmse 0.000850639  max resid 0.001749644 
    ## ... Similar to previous best
    ## Run 335 stress 0.01883125 
    ## ... Procrustes: rmse 0.002796739  max resid 0.005301742 
    ## ... Similar to previous best
    ## Run 336 stress 0.02492805 
    ## Run 337 stress 0.0187957 
    ## ... Procrustes: rmse 0.0008469066  max resid 0.00174212 
    ## ... Similar to previous best
    ## Run 338 stress 0.02486136 
    ## Run 339 stress 0.01879624 
    ## ... Procrustes: rmse 0.001062409  max resid 0.002185691 
    ## ... Similar to previous best
    ## Run 340 stress 0.01879567 
    ## ... Procrustes: rmse 0.0008391583  max resid 0.001726066 
    ## ... Similar to previous best
    ## Run 341 stress 0.01884501 
    ## ... Procrustes: rmse 0.003570898  max resid 0.006950926 
    ## ... Similar to previous best
    ## Run 342 stress 0.02492191 
    ## Run 343 stress 0.02509449 
    ## Run 344 stress 0.01879601 
    ## ... Procrustes: rmse 0.0005713806  max resid 0.001177866 
    ## ... Similar to previous best
    ## Run 345 stress 0.01879601 
    ## ... Procrustes: rmse 0.0009786385  max resid 0.002013293 
    ## ... Similar to previous best
    ## Run 346 stress 0.02509437 
    ## Run 347 stress 0.01879581 
    ## ... Procrustes: rmse 0.0004979208  max resid 0.001026817 
    ## ... Similar to previous best
    ## Run 348 stress 0.02492163 
    ## Run 349 stress 0.0248612 
    ## Run 350 stress 0.01890844 
    ## ... Procrustes: rmse 0.007954431  max resid 0.01613193 
    ## Run 351 stress 0.02509478 
    ## Run 352 stress 0.01879595 
    ## ... Procrustes: rmse 0.0009550671  max resid 0.00196479 
    ## ... Similar to previous best
    ## Run 353 stress 0.01879591 
    ## ... Procrustes: rmse 0.0009413678  max resid 0.001936544 
    ## ... Similar to previous best
    ## Run 354 stress 0.02520876 
    ## Run 355 stress 0.0248613 
    ## Run 356 stress 0.02509444 
    ## Run 357 stress 0.01883245 
    ## ... Procrustes: rmse 0.002965512  max resid 0.005660459 
    ## ... Similar to previous best
    ## Run 358 stress 0.02509461 
    ## Run 359 stress 0.01879597 
    ## ... Procrustes: rmse 0.0009640709  max resid 0.001983348 
    ## ... Similar to previous best
    ## Run 360 stress 0.0187957 
    ## ... Procrustes: rmse 0.0004553837  max resid 0.0009395224 
    ## ... Similar to previous best
    ## Run 361 stress 0.02520888 
    ## Run 362 stress 0.02509479 
    ## Run 363 stress 0.0187954 
    ## ... Procrustes: rmse 0.000695349  max resid 0.001429901 
    ## ... Similar to previous best
    ## Run 364 stress 0.01879573 
    ## ... Procrustes: rmse 0.0004727132  max resid 0.0009751789 
    ## ... Similar to previous best
    ## Run 365 stress 0.01884655 
    ## ... Procrustes: rmse 0.004439486  max resid 0.008796969 
    ## ... Similar to previous best
    ## Run 366 stress 0.02509468 
    ## Run 367 stress 0.01883446 
    ## ... Procrustes: rmse 0.00323237  max resid 0.006231196 
    ## ... Similar to previous best
    ## Run 368 stress 0.02509473 
    ## Run 369 stress 0.01891097 
    ## ... Procrustes: rmse 0.00864354  max resid 0.01774491 
    ## Run 370 stress 0.01887624 
    ## ... Procrustes: rmse 0.00638824  max resid 0.01287707 
    ## Run 371 stress 0.01887952 
    ## ... Procrustes: rmse 0.006577369  max resid 0.01327069 
    ## Run 372 stress 0.02486108 
    ## Run 373 stress 0.01887916 
    ## ... Procrustes: rmse 0.007370311  max resid 0.01514117 
    ## Run 374 stress 0.02492204 
    ## Run 375 stress 0.02492185 
    ## Run 376 stress 0.02509464 
    ## Run 377 stress 0.01880697 
    ## ... Procrustes: rmse 0.002885982  max resid 0.005941027 
    ## ... Similar to previous best
    ## Run 378 stress 0.02486148 
    ## Run 379 stress 0.02492179 
    ## Run 380 stress 0.02509467 
    ## Run 381 stress 0.01883715 
    ## ... Procrustes: rmse 0.003545935  max resid 0.006899574 
    ## ... Similar to previous best
    ## Run 382 stress 0.01879588 
    ## ... Procrustes: rmse 0.0005257487  max resid 0.001083978 
    ## ... Similar to previous best
    ## Run 383 stress 0.01883209 
    ## ... Procrustes: rmse 0.002912088  max resid 0.005544345 
    ## ... Similar to previous best
    ## Run 384 stress 0.01950112 
    ## Run 385 stress 0.02492194 
    ## Run 386 stress 0.01898691 
    ## ... Procrustes: rmse 0.01015291  max resid 0.02084503 
    ## Run 387 stress 0.02509428 
    ## Run 388 stress 0.01882846 
    ## ... Procrustes: rmse 0.002409114  max resid 0.004555687 
    ## ... Similar to previous best
    ## Run 389 stress 0.02509452 
    ## Run 390 stress 0.02492177 
    ## Run 391 stress 0.02509452 
    ## Run 392 stress 0.02509468 
    ## Run 393 stress 0.01879563 
    ## ... Procrustes: rmse 0.0004175997  max resid 0.0008617414 
    ## ... Similar to previous best
    ## Run 394 stress 0.01879561 
    ## ... Procrustes: rmse 0.0008023354  max resid 0.001650139 
    ## ... Similar to previous best
    ## Run 395 stress 0.02492202 
    ## Run 396 stress 0.02492183 
    ## Run 397 stress 0.02486117 
    ## Run 398 stress 0.01879556 
    ## ... Procrustes: rmse 0.0007834914  max resid 0.001611466 
    ## ... Similar to previous best
    ## Run 399 stress 0.01879568 
    ## ... Procrustes: rmse 0.0008388331  max resid 0.001725368 
    ## ... Similar to previous best
    ## Run 400 stress 0.01879577 
    ## ... Procrustes: rmse 0.0008788194  max resid 0.001807824 
    ## ... Similar to previous best
    ## Run 401 stress 0.2261475 
    ## Run 402 stress 0.02509464 
    ## Run 403 stress 0.01879563 
    ## ... Procrustes: rmse 0.0008190145  max resid 0.00168443 
    ## ... Similar to previous best
    ## Run 404 stress 0.01884356 
    ## ... Procrustes: rmse 0.005627312  max resid 0.01157045 
    ## Run 405 stress 0.01879997 
    ## ... Procrustes: rmse 0.001926382  max resid 0.003965661 
    ## ... Similar to previous best
    ## Run 406 stress 0.02492176 
    ## Run 407 stress 0.01879557 
    ## ... Procrustes: rmse 0.0007901727  max resid 0.001625237 
    ## ... Similar to previous best
    ## Run 408 stress 0.01884379 
    ## ... Procrustes: rmse 0.004203735  max resid 0.008298556 
    ## ... Similar to previous best
    ## Run 409 stress 0.02492182 
    ## Run 410 stress 0.01880895 
    ## ... Procrustes: rmse 0.003100717  max resid 0.00638326 
    ## ... Similar to previous best
    ## Run 411 stress 0.02492208 
    ## Run 412 stress 0.02492168 
    ## Run 413 stress 0.01882997 
    ## ... Procrustes: rmse 0.002602394  max resid 0.004930681 
    ## ... Similar to previous best
    ## Run 414 stress 0.02520891 
    ## Run 415 stress 0.01879624 
    ## ... Procrustes: rmse 0.001064111  max resid 0.002189275 
    ## ... Similar to previous best
    ## Run 416 stress 0.01879524 
    ## ... Procrustes: rmse 0.000588123  max resid 0.001209016 
    ## ... Similar to previous best
    ## Run 417 stress 0.02492173 
    ## Run 418 stress 0.02492175 
    ## Run 419 stress 0.01879565 
    ## ... Procrustes: rmse 0.0004221986  max resid 0.0008713639 
    ## ... Similar to previous best
    ## Run 420 stress 0.02509483 
    ## Run 421 stress 0.02520893 
    ## Run 422 stress 0.0187955 
    ## ... Procrustes: rmse 0.0003599392  max resid 0.0007434532 
    ## ... Similar to previous best
    ## Run 423 stress 0.01879616 
    ## ... Procrustes: rmse 0.001018251  max resid 0.002094558 
    ## ... Similar to previous best
    ## Run 424 stress 0.01879586 
    ## ... Procrustes: rmse 0.0009130199  max resid 0.001878239 
    ## ... Similar to previous best
    ## Run 425 stress 0.02509483 
    ## Run 426 stress 0.02509479 
    ## Run 427 stress 0.02520875 
    ## Run 428 stress 0.02509454 
    ## Run 429 stress 0.0250947 
    ## Run 430 stress 0.01880115 
    ## ... Procrustes: rmse 0.002119459  max resid 0.004363597 
    ## ... Similar to previous best
    ## Run 431 stress 0.01879567 
    ## ... Procrustes: rmse 0.0008358718  max resid 0.001719362 
    ## ... Similar to previous best
    ## Run 432 stress 0.01879595 
    ## ... Procrustes: rmse 0.0009453715  max resid 0.001944885 
    ## ... Similar to previous best
    ## Run 433 stress 0.01884817 
    ## ... Procrustes: rmse 0.004572248  max resid 0.009077814 
    ## ... Similar to previous best
    ## Run 434 stress 0.01879568 
    ## ... Procrustes: rmse 0.0008207427  max resid 0.001687925 
    ## ... Similar to previous best
    ## Run 435 stress 0.0187958 
    ## ... Procrustes: rmse 0.0008865441  max resid 0.00182375 
    ## ... Similar to previous best
    ## Run 436 stress 0.01879591 
    ## ... Procrustes: rmse 0.0009190979  max resid 0.001890424 
    ## ... Similar to previous best
    ## Run 437 stress 0.01880469 
    ## ... Procrustes: rmse 0.002613857  max resid 0.005381613 
    ## ... Similar to previous best
    ## Run 438 stress 0.02509441 
    ## Run 439 stress 0.01892618 
    ## ... Procrustes: rmse 0.008633906  max resid 0.01754052 
    ## Run 440 stress 0.01895221 
    ## ... Procrustes: rmse 0.009662013  max resid 0.0196644 
    ## Run 441 stress 0.01879597 
    ## ... Procrustes: rmse 0.0009637178  max resid 0.001982572 
    ## ... Similar to previous best
    ## Run 442 stress 0.02509458 
    ## Run 443 stress 0.01879595 
    ## ... Procrustes: rmse 0.0005554866  max resid 0.001145091 
    ## ... Similar to previous best
    ## Run 444 stress 0.02486128 
    ## Run 445 stress 0.02509454 
    ## Run 446 stress 0.01889161 
    ## ... Procrustes: rmse 0.007191776  max resid 0.01455028 
    ## Run 447 stress 0.02509471 
    ## Run 448 stress 0.02509458 
    ## Run 449 stress 0.02509464 
    ## Run 450 stress 0.02520879 
    ## Run 451 stress 0.02486134 
    ## Run 452 stress 0.02520886 
    ## Run 453 stress 0.01885157 
    ## ... Procrustes: rmse 0.00483849  max resid 0.009638483 
    ## ... Similar to previous best
    ## Run 454 stress 0.02486118 
    ## Run 455 stress 0.02492182 
    ## Run 456 stress 0.01880085 
    ## ... Procrustes: rmse 0.002042994  max resid 0.004206536 
    ## ... Similar to previous best
    ## Run 457 stress 0.01925176 
    ## ... Procrustes: rmse 0.01636293  max resid 0.03360185 
    ## Run 458 stress 0.02492205 
    ## Run 459 stress 0.01879565 
    ## ... Procrustes: rmse 0.0008280166  max resid 0.001703123 
    ## ... Similar to previous best
    ## Run 460 stress 0.02492154 
    ## Run 461 stress 0.02486144 
    ## Run 462 stress 0.0187956 
    ## ... Procrustes: rmse 0.0007946467  max resid 0.001634482 
    ## ... Similar to previous best
    ## Run 463 stress 0.02509443 
    ## Run 464 stress 0.01879573 
    ## ... Procrustes: rmse 0.0008649793  max resid 0.001779231 
    ## ... Similar to previous best
    ## Run 465 stress 0.01882748 
    ## ... Procrustes: rmse 0.004636104  max resid 0.009535706 
    ## ... Similar to previous best
    ## Run 466 stress 0.02486126 
    ## Run 467 stress 0.02492203 
    ## Run 468 stress 0.01879585 
    ## ... Procrustes: rmse 0.000914167  max resid 0.001880582 
    ## ... Similar to previous best
    ## Run 469 stress 0.02509455 
    ## Run 470 stress 0.01879579 
    ## ... Procrustes: rmse 0.0008789595  max resid 0.001808164 
    ## ... Similar to previous best
    ## Run 471 stress 0.01879573 
    ## ... Procrustes: rmse 0.0008628444  max resid 0.001773933 
    ## ... Similar to previous best
    ## Run 472 stress 0.01882043 
    ## ... Procrustes: rmse 0.0041253  max resid 0.008486759 
    ## ... Similar to previous best
    ## Run 473 stress 0.01879594 
    ## ... Procrustes: rmse 0.0009438573  max resid 0.001941524 
    ## ... Similar to previous best
    ## Run 474 stress 0.01879888 
    ## ... Procrustes: rmse 0.00165615  max resid 0.003409081 
    ## ... Similar to previous best
    ## Run 475 stress 0.01879561 
    ## ... Procrustes: rmse 0.000799487  max resid 0.001644467 
    ## ... Similar to previous best
    ## Run 476 stress 0.02492764 
    ## Run 477 stress 0.01883801 
    ## ... Procrustes: rmse 0.003574781  max resid 0.006956387 
    ## ... Similar to previous best
    ## Run 478 stress 0.02509458 
    ## Run 479 stress 0.01879585 
    ## ... Procrustes: rmse 0.0009044745  max resid 0.001860491 
    ## ... Similar to previous best
    ## Run 480 stress 0.01887404 
    ## ... Procrustes: rmse 0.00708596  max resid 0.01455872 
    ## Run 481 stress 0.01886563 
    ## ... Procrustes: rmse 0.005739575  max resid 0.01152126 
    ## Run 482 stress 0.02492162 
    ## Run 483 stress 0.01879758 
    ## ... Procrustes: rmse 0.00136355  max resid 0.002806016 
    ## ... Similar to previous best
    ## Run 484 stress 0.01884132 
    ## ... Procrustes: rmse 0.003974563  max resid 0.007812676 
    ## ... Similar to previous best
    ## Run 485 stress 0.0249219 
    ## Run 486 stress 0.01879593 
    ## ... Procrustes: rmse 0.000949751  max resid 0.001953821 
    ## ... Similar to previous best
    ## Run 487 stress 0.01879545 
    ## ... Procrustes: rmse 0.0007213602  max resid 0.001483752 
    ## ... Similar to previous best
    ## Run 488 stress 0.3068721 
    ## Run 489 stress 0.01879598 
    ## ... Procrustes: rmse 0.0009669614  max resid 0.00198924 
    ## ... Similar to previous best
    ## Run 490 stress 0.024922 
    ## Run 491 stress 0.0250947 
    ## Run 492 stress 0.0187958 
    ## ... Procrustes: rmse 0.0008946442  max resid 0.001840325 
    ## ... Similar to previous best
    ## Run 493 stress 0.02492166 
    ## Run 494 stress 0.01879581 
    ## ... Procrustes: rmse 0.0008995805  max resid 0.001850546 
    ## ... Similar to previous best
    ## Run 495 stress 0.01880349 
    ## ... Procrustes: rmse 0.002297206  max resid 0.004727101 
    ## ... Similar to previous best
    ## Run 496 stress 0.01879576 
    ## ... Procrustes: rmse 0.0008791252  max resid 0.001808415 
    ## ... Similar to previous best
    ## Run 497 stress 0.02509473 
    ## Run 498 stress 0.02520911 
    ## Run 499 stress 0.01879576 
    ## ... Procrustes: rmse 0.0004847233  max resid 0.0009998624 
    ## ... Similar to previous best
    ## Run 500 stress 0.02520908 
    ## *** Best solution repeated 218 times

``` r
# Mixed and stratified lakes
SD_beta_MS_NMDS <- metaMDS(SD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09030412 
    ## Run 2 stress 0.0850347 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001544593  max resid 0.003996185 
    ## ... Similar to previous best
    ## Run 3 stress 0.09539194 
    ## Run 4 stress 0.09145326 
    ## Run 5 stress 0.09908184 
    ## Run 6 stress 0.09416135 
    ## Run 7 stress 0.09908187 
    ## Run 8 stress 0.08440278 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01420246  max resid 0.04313529 
    ## Run 9 stress 0.09464397 
    ## Run 10 stress 0.09286083 
    ## Run 11 stress 0.09030401 
    ## Run 12 stress 0.08440266 
    ## ... New best solution
    ## ... Procrustes: rmse 7.878242e-05  max resid 0.0001397726 
    ## ... Similar to previous best
    ## Run 13 stress 0.08440255 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001192009  max resid 0.0002117599 
    ## ... Similar to previous best
    ## Run 14 stress 0.09721197 
    ## Run 15 stress 0.09268384 
    ## Run 16 stress 0.08773469 
    ## Run 17 stress 0.08503502 
    ## Run 18 stress 0.0844028 
    ## ... Procrustes: rmse 0.0001712239  max resid 0.0003619811 
    ## ... Similar to previous best
    ## Run 19 stress 0.08773482 
    ## Run 20 stress 0.08440256 
    ## ... Procrustes: rmse 4.049664e-05  max resid 7.553224e-05 
    ## ... Similar to previous best
    ## Run 21 stress 0.09030398 
    ## Run 22 stress 0.090304 
    ## Run 23 stress 0.09321784 
    ## Run 24 stress 0.1038039 
    ## Run 25 stress 0.09286083 
    ## Run 26 stress 0.09374312 
    ## Run 27 stress 0.09374261 
    ## Run 28 stress 0.08440258 
    ## ... Procrustes: rmse 3.776095e-05  max resid 6.876018e-05 
    ## ... Similar to previous best
    ## Run 29 stress 0.0903041 
    ## Run 30 stress 0.100461 
    ## Run 31 stress 0.09337249 
    ## Run 32 stress 0.0916895 
    ## Run 33 stress 0.08773473 
    ## Run 34 stress 0.08973874 
    ## Run 35 stress 0.09380586 
    ## Run 36 stress 0.1038041 
    ## Run 37 stress 0.09337267 
    ## Run 38 stress 0.09168932 
    ## Run 39 stress 0.09381846 
    ## Run 40 stress 0.0926839 
    ## Run 41 stress 0.09539192 
    ## Run 42 stress 0.08440268 
    ## ... Procrustes: rmse 0.0001334499  max resid 0.0002386055 
    ## ... Similar to previous best
    ## Run 43 stress 0.09374271 
    ## Run 44 stress 0.09159091 
    ## Run 45 stress 0.09159086 
    ## Run 46 stress 0.09308955 
    ## Run 47 stress 0.09030399 
    ## Run 48 stress 0.08440256 
    ## ... Procrustes: rmse 4.83975e-05  max resid 0.0001087549 
    ## ... Similar to previous best
    ## Run 49 stress 0.08973873 
    ## Run 50 stress 0.09374215 
    ## Run 51 stress 0.08773479 
    ## Run 52 stress 0.08440261 
    ## ... Procrustes: rmse 9.341141e-05  max resid 0.0001723534 
    ## ... Similar to previous best
    ## Run 53 stress 0.08973863 
    ## Run 54 stress 0.08503502 
    ## Run 55 stress 0.09407973 
    ## Run 56 stress 0.08503631 
    ## Run 57 stress 0.09159083 
    ## Run 58 stress 0.09407986 
    ## Run 59 stress 0.09159084 
    ## Run 60 stress 0.09030403 
    ## Run 61 stress 0.1026215 
    ## Run 62 stress 0.09408008 
    ## Run 63 stress 0.09535489 
    ## Run 64 stress 0.08503617 
    ## Run 65 stress 0.08503466 
    ## Run 66 stress 0.09159086 
    ## Run 67 stress 0.09407985 
    ## Run 68 stress 0.105285 
    ## Run 69 stress 0.09374259 
    ## Run 70 stress 0.08440257 
    ## ... Procrustes: rmse 1.058461e-05  max resid 2.18943e-05 
    ## ... Similar to previous best
    ## Run 71 stress 0.08503519 
    ## Run 72 stress 0.08503775 
    ## Run 73 stress 0.1026217 
    ## Run 74 stress 0.09308955 
    ## Run 75 stress 0.08973881 
    ## Run 76 stress 0.09374291 
    ## Run 77 stress 0.09308973 
    ## Run 78 stress 0.09380578 
    ## Run 79 stress 0.09337236 
    ## Run 80 stress 0.09407983 
    ## Run 81 stress 0.09465905 
    ## Run 82 stress 0.09465907 
    ## Run 83 stress 0.08773482 
    ## Run 84 stress 0.08973889 
    ## Run 85 stress 0.09030395 
    ## Run 86 stress 0.09465908 
    ## Run 87 stress 0.09030407 
    ## Run 88 stress 0.09337265 
    ## Run 89 stress 0.08440271 
    ## ... Procrustes: rmse 0.0001483765  max resid 0.0002660096 
    ## ... Similar to previous best
    ## Run 90 stress 0.1038045 
    ## Run 91 stress 0.09969964 
    ## Run 92 stress 0.0926839 
    ## Run 93 stress 0.09381844 
    ## Run 94 stress 0.09407992 
    ## Run 95 stress 0.09145321 
    ## Run 96 stress 0.09669863 
    ## Run 97 stress 0.0940799 
    ## Run 98 stress 0.09407964 
    ## Run 99 stress 0.1052852 
    ## Run 100 stress 0.08773483 
    ## Run 101 stress 0.09308977 
    ## Run 102 stress 0.2868406 
    ## Run 103 stress 0.1038041 
    ## Run 104 stress 0.1038053 
    ## Run 105 stress 0.09969966 
    ## Run 106 stress 0.08773476 
    ## Run 107 stress 0.09374294 
    ## Run 108 stress 0.09760824 
    ## Run 109 stress 0.09286084 
    ## Run 110 stress 0.09408002 
    ## Run 111 stress 0.08503527 
    ## Run 112 stress 0.0937429 
    ## Run 113 stress 0.09268347 
    ## Run 114 stress 0.09760825 
    ## Run 115 stress 0.09760814 
    ## Run 116 stress 0.1038052 
    ## Run 117 stress 0.09590479 
    ## Run 118 stress 0.0930896 
    ## Run 119 stress 0.08973862 
    ## Run 120 stress 0.0850352 
    ## Run 121 stress 0.08503644 
    ## Run 122 stress 0.09268413 
    ## Run 123 stress 0.08773471 
    ## Run 124 stress 0.09408005 
    ## Run 125 stress 0.09145331 
    ## Run 126 stress 0.08503613 
    ## Run 127 stress 0.09168939 
    ## Run 128 stress 0.09969966 
    ## Run 129 stress 0.09268352 
    ## Run 130 stress 0.09308946 
    ## Run 131 stress 0.09407989 
    ## Run 132 stress 0.2528099 
    ## Run 133 stress 0.09407985 
    ## Run 134 stress 0.09168936 
    ## Run 135 stress 0.09145321 
    ## Run 136 stress 0.08773486 
    ## Run 137 stress 0.090304 
    ## Run 138 stress 0.09337229 
    ## Run 139 stress 0.09159099 
    ## Run 140 stress 0.08773483 
    ## Run 141 stress 0.08440257 
    ## ... Procrustes: rmse 2.189376e-05  max resid 6.847898e-05 
    ## ... Similar to previous best
    ## Run 142 stress 0.09168928 
    ## Run 143 stress 0.1038044 
    ## Run 144 stress 0.09969969 
    ## Run 145 stress 0.0926841 
    ## Run 146 stress 0.09374248 
    ## Run 147 stress 0.09407112 
    ## Run 148 stress 0.09286098 
    ## Run 149 stress 0.09030405 
    ## Run 150 stress 0.09969995 
    ## Run 151 stress 0.08440256 
    ## ... Procrustes: rmse 2.68638e-05  max resid 5.624008e-05 
    ## ... Similar to previous best
    ## Run 152 stress 0.0850351 
    ## Run 153 stress 0.09465904 
    ## Run 154 stress 0.09286086 
    ## Run 155 stress 0.09030407 
    ## Run 156 stress 0.09268345 
    ## Run 157 stress 0.08503492 
    ## Run 158 stress 0.08440256 
    ## ... Procrustes: rmse 3.357929e-05  max resid 6.792563e-05 
    ## ... Similar to previous best
    ## Run 159 stress 0.09268398 
    ## Run 160 stress 0.09539197 
    ## Run 161 stress 0.08440261 
    ## ... Procrustes: rmse 9.505675e-05  max resid 0.0001790074 
    ## ... Similar to previous best
    ## Run 162 stress 0.09159084 
    ## Run 163 stress 0.0946591 
    ## Run 164 stress 0.09145325 
    ## Run 165 stress 0.09030413 
    ## Run 166 stress 0.09145328 
    ## Run 167 stress 0.09030418 
    ## Run 168 stress 0.09030414 
    ## Run 169 stress 0.09030405 
    ## Run 170 stress 0.08440262 
    ## ... Procrustes: rmse 9.246557e-05  max resid 0.0001661606 
    ## ... Similar to previous best
    ## Run 171 stress 0.09416131 
    ## Run 172 stress 0.09268354 
    ## Run 173 stress 0.08773466 
    ## Run 174 stress 0.08503502 
    ## Run 175 stress 0.09400545 
    ## Run 176 stress 0.09403414 
    ## Run 177 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001688667  max resid 0.0003189829 
    ## ... Similar to previous best
    ## Run 178 stress 0.09407995 
    ## Run 179 stress 0.09145331 
    ## Run 180 stress 0.09030394 
    ## Run 181 stress 0.09374304 
    ## Run 182 stress 0.09969988 
    ## Run 183 stress 0.09268358 
    ## Run 184 stress 0.08440253 
    ## ... Procrustes: rmse 0.000127807  max resid 0.0002645085 
    ## ... Similar to previous best
    ## Run 185 stress 0.09030398 
    ## Run 186 stress 0.09337266 
    ## Run 187 stress 0.0938058 
    ## Run 188 stress 0.09464429 
    ## Run 189 stress 0.08503584 
    ## Run 190 stress 0.09268357 
    ## Run 191 stress 0.08973881 
    ## Run 192 stress 0.09030411 
    ## Run 193 stress 0.0930895 
    ## Run 194 stress 0.09407975 
    ## Run 195 stress 0.09268337 
    ## Run 196 stress 0.09712929 
    ## Run 197 stress 0.09286083 
    ## Run 198 stress 0.09145321 
    ## Run 199 stress 0.2475294 
    ## Run 200 stress 0.09145323 
    ## Run 201 stress 0.08973889 
    ## Run 202 stress 0.08440261 
    ## ... Procrustes: rmse 0.0002471461  max resid 0.0004888426 
    ## ... Similar to previous best
    ## Run 203 stress 0.1038038 
    ## Run 204 stress 0.08503521 
    ## Run 205 stress 0.09416132 
    ## Run 206 stress 0.08503524 
    ## Run 207 stress 0.09465905 
    ## Run 208 stress 0.08973872 
    ## Run 209 stress 0.08773485 
    ## Run 210 stress 0.08503492 
    ## Run 211 stress 0.09465905 
    ## Run 212 stress 0.09030395 
    ## Run 213 stress 0.09268318 
    ## Run 214 stress 0.08773478 
    ## Run 215 stress 0.09159083 
    ## Run 216 stress 0.1095795 
    ## Run 217 stress 0.09464499 
    ## Run 218 stress 0.09337299 
    ## Run 219 stress 0.09535487 
    ## Run 220 stress 0.08973864 
    ## Run 221 stress 0.09268375 
    ## Run 222 stress 0.09030412 
    ## Run 223 stress 0.08503472 
    ## Run 224 stress 0.09145321 
    ## Run 225 stress 0.0940797 
    ## Run 226 stress 0.09145329 
    ## Run 227 stress 0.09030393 
    ## Run 228 stress 0.1038033 
    ## Run 229 stress 0.08973866 
    ## Run 230 stress 0.1026216 
    ## Run 231 stress 0.09416132 
    ## Run 232 stress 0.0928609 
    ## Run 233 stress 0.09286084 
    ## Run 234 stress 0.09159087 
    ## Run 235 stress 0.09400517 
    ## Run 236 stress 0.08503463 
    ## Run 237 stress 0.08503598 
    ## Run 238 stress 0.09464403 
    ## Run 239 stress 0.09159083 
    ## Run 240 stress 0.09286084 
    ## Run 241 stress 0.09159087 
    ## Run 242 stress 0.09465908 
    ## Run 243 stress 0.09590879 
    ## Run 244 stress 0.09286083 
    ## Run 245 stress 0.09380588 
    ## Run 246 stress 0.09030405 
    ## Run 247 stress 0.09030397 
    ## Run 248 stress 0.08973892 
    ## Run 249 stress 0.1038051 
    ## Run 250 stress 0.09030411 
    ## Run 251 stress 0.09268346 
    ## Run 252 stress 0.09465904 
    ## Run 253 stress 0.09590875 
    ## Run 254 stress 0.09030395 
    ## Run 255 stress 0.08973867 
    ## Run 256 stress 0.09908184 
    ## Run 257 stress 0.09168953 
    ## Run 258 stress 0.09159084 
    ## Run 259 stress 0.09408 
    ## Run 260 stress 0.09030411 
    ## Run 261 stress 0.090304 
    ## Run 262 stress 0.09407977 
    ## Run 263 stress 0.09286086 
    ## Run 264 stress 0.09268391 
    ## Run 265 stress 0.09465904 
    ## Run 266 stress 0.08973882 
    ## Run 267 stress 0.09268332 
    ## Run 268 stress 0.09407985 
    ## Run 269 stress 0.09286084 
    ## Run 270 stress 0.09030404 
    ## Run 271 stress 0.09268339 
    ## Run 272 stress 0.09321655 
    ## Run 273 stress 0.08973862 
    ## Run 274 stress 0.08973894 
    ## Run 275 stress 0.09337238 
    ## Run 276 stress 0.08503637 
    ## Run 277 stress 0.09381844 
    ## Run 278 stress 0.09464449 
    ## Run 279 stress 0.09308995 
    ## Run 280 stress 0.09168938 
    ## Run 281 stress 0.09408002 
    ## Run 282 stress 0.08440261 
    ## ... Procrustes: rmse 0.0002363684  max resid 0.0004449046 
    ## ... Similar to previous best
    ## Run 283 stress 0.09030409 
    ## Run 284 stress 0.09445467 
    ## Run 285 stress 0.09308949 
    ## Run 286 stress 0.08503657 
    ## Run 287 stress 0.09465908 
    ## Run 288 stress 0.08773486 
    ## Run 289 stress 0.09268369 
    ## Run 290 stress 0.08440259 
    ## ... Procrustes: rmse 0.0002090308  max resid 0.0004424583 
    ## ... Similar to previous best
    ## Run 291 stress 0.09030404 
    ## Run 292 stress 0.09408009 
    ## Run 293 stress 0.09145339 
    ## Run 294 stress 0.09159083 
    ## Run 295 stress 0.09408006 
    ## Run 296 stress 0.08503497 
    ## Run 297 stress 0.094034 
    ## Run 298 stress 0.09969969 
    ## Run 299 stress 0.2911124 
    ## Run 300 stress 0.09159083 
    ## Run 301 stress 0.09374336 
    ## Run 302 stress 0.09535486 
    ## Run 303 stress 0.09908188 
    ## Run 304 stress 0.09381844 
    ## Run 305 stress 0.08503459 
    ## Run 306 stress 0.08503657 
    ## Run 307 stress 0.09712977 
    ## Run 308 stress 0.09308976 
    ## Run 309 stress 0.08440253 
    ## ... Procrustes: rmse 0.0001543007  max resid 0.0002954023 
    ## ... Similar to previous best
    ## Run 310 stress 0.09286086 
    ## Run 311 stress 0.09030394 
    ## Run 312 stress 0.09465904 
    ## Run 313 stress 0.08973891 
    ## Run 314 stress 0.09535447 
    ## Run 315 stress 0.08503493 
    ## Run 316 stress 0.08773468 
    ## Run 317 stress 0.09159084 
    ## Run 318 stress 0.09159083 
    ## Run 319 stress 0.1053423 
    ## Run 320 stress 0.09337259 
    ## Run 321 stress 0.09337229 
    ## Run 322 stress 0.08503458 
    ## Run 323 stress 0.08773474 
    ## Run 324 stress 0.09380587 
    ## Run 325 stress 0.09721198 
    ## Run 326 stress 0.08973868 
    ## Run 327 stress 0.08773482 
    ## Run 328 stress 0.09030394 
    ## Run 329 stress 0.1038047 
    ## Run 330 stress 0.09286087 
    ## Run 331 stress 0.08973885 
    ## Run 332 stress 0.09760821 
    ## Run 333 stress 0.09286083 
    ## Run 334 stress 0.09030402 
    ## Run 335 stress 0.08440263 
    ## ... Procrustes: rmse 0.0002645176  max resid 0.0005083197 
    ## ... Similar to previous best
    ## Run 336 stress 0.09030417 
    ## Run 337 stress 0.09465908 
    ## Run 338 stress 0.09030396 
    ## Run 339 stress 0.08773467 
    ## Run 340 stress 0.1006365 
    ## Run 341 stress 0.09145322 
    ## Run 342 stress 0.08503488 
    ## Run 343 stress 0.1038032 
    ## Run 344 stress 0.09539203 
    ## Run 345 stress 0.09590757 
    ## Run 346 stress 0.08973871 
    ## Run 347 stress 0.09760818 
    ## Run 348 stress 0.09374327 
    ## Run 349 stress 0.09465914 
    ## Run 350 stress 0.09159085 
    ## Run 351 stress 0.09669919 
    ## Run 352 stress 0.09286096 
    ## Run 353 stress 0.09268327 
    ## Run 354 stress 0.2488812 
    ## Run 355 stress 0.1112172 
    ## Run 356 stress 0.08503682 
    ## Run 357 stress 0.09145336 
    ## Run 358 stress 0.09407995 
    ## Run 359 stress 0.0930895 
    ## Run 360 stress 0.09535423 
    ## Run 361 stress 0.09539198 
    ## Run 362 stress 0.09407983 
    ## Run 363 stress 0.09308966 
    ## Run 364 stress 0.09168955 
    ## Run 365 stress 0.09030407 
    ## Run 366 stress 0.09374211 
    ## Run 367 stress 0.08773469 
    ## Run 368 stress 0.09404394 
    ## Run 369 stress 0.08440261 
    ## ... Procrustes: rmse 0.000136658  max resid 0.0002795738 
    ## ... Similar to previous best
    ## Run 370 stress 0.09321806 
    ## Run 371 stress 0.09337275 
    ## Run 372 stress 0.09445465 
    ## Run 373 stress 0.1038048 
    ## Run 374 stress 0.09030396 
    ## Run 375 stress 0.09712977 
    ## Run 376 stress 0.0915909 
    ## Run 377 stress 0.09408012 
    ## Run 378 stress 0.09030396 
    ## Run 379 stress 0.09268323 
    ## Run 380 stress 0.09539192 
    ## Run 381 stress 0.09908181 
    ## Run 382 stress 0.09416132 
    ## Run 383 stress 0.2914682 
    ## Run 384 stress 0.0903041 
    ## Run 385 stress 0.0953545 
    ## Run 386 stress 0.09030397 
    ## Run 387 stress 0.08440261 
    ## ... Procrustes: rmse 0.0002392426  max resid 0.0004393958 
    ## ... Similar to previous best
    ## Run 388 stress 0.09030419 
    ## Run 389 stress 0.09030417 
    ## Run 390 stress 0.09286085 
    ## Run 391 stress 0.09465904 
    ## Run 392 stress 0.09145321 
    ## Run 393 stress 0.09407964 
    ## Run 394 stress 0.09159088 
    ## Run 395 stress 0.09381852 
    ## Run 396 stress 0.08503469 
    ## Run 397 stress 0.08973878 
    ## Run 398 stress 0.09159083 
    ## Run 399 stress 0.09464514 
    ## Run 400 stress 0.3304983 
    ## Run 401 stress 0.09030398 
    ## Run 402 stress 0.08773494 
    ## Run 403 stress 0.09535465 
    ## Run 404 stress 0.09030393 
    ## Run 405 stress 0.09465907 
    ## Run 406 stress 0.09030395 
    ## Run 407 stress 0.09030413 
    ## Run 408 stress 0.08440273 
    ## ... Procrustes: rmse 0.0002994685  max resid 0.0005517862 
    ## ... Similar to previous best
    ## Run 409 stress 0.09407997 
    ## Run 410 stress 0.09407992 
    ## Run 411 stress 0.09374298 
    ## Run 412 stress 0.09268365 
    ## Run 413 stress 0.0941615 
    ## Run 414 stress 0.09030412 
    ## Run 415 stress 0.09337323 
    ## Run 416 stress 0.09337267 
    ## Run 417 stress 0.09268333 
    ## Run 418 stress 0.08440257 
    ## ... Procrustes: rmse 0.0002044536  max resid 0.0003871115 
    ## ... Similar to previous best
    ## Run 419 stress 0.08440272 
    ## ... Procrustes: rmse 0.0003296299  max resid 0.0006093944 
    ## ... Similar to previous best
    ## Run 420 stress 0.09030394 
    ## Run 421 stress 0.09760823 
    ## Run 422 stress 0.0940797 
    ## Run 423 stress 0.09030403 
    ## Run 424 stress 0.09145328 
    ## Run 425 stress 0.08503516 
    ## Run 426 stress 0.09159101 
    ## Run 427 stress 0.09168957 
    ## Run 428 stress 0.1038084 
    ## Run 429 stress 0.09168942 
    ## Run 430 stress 0.09286083 
    ## Run 431 stress 0.09268369 
    ## Run 432 stress 0.08503482 
    ## Run 433 stress 0.09407985 
    ## Run 434 stress 0.09535495 
    ## Run 435 stress 0.0903041 
    ## Run 436 stress 0.09374258 
    ## Run 437 stress 0.08773468 
    ## Run 438 stress 0.09590891 
    ## Run 439 stress 0.08773493 
    ## Run 440 stress 0.08973867 
    ## Run 441 stress 0.09590632 
    ## Run 442 stress 0.09159084 
    ## Run 443 stress 0.09286088 
    ## Run 444 stress 0.08440258 
    ## ... Procrustes: rmse 0.000213996  max resid 0.0004059202 
    ## ... Similar to previous best
    ## Run 445 stress 0.09030396 
    ## Run 446 stress 0.09416133 
    ## Run 447 stress 0.09168934 
    ## Run 448 stress 0.09969977 
    ## Run 449 stress 0.09286087 
    ## Run 450 stress 0.09407997 
    ## Run 451 stress 0.08440254 
    ## ... Procrustes: rmse 5.928765e-05  max resid 0.0001312361 
    ## ... Similar to previous best
    ## Run 452 stress 0.09464495 
    ## Run 453 stress 0.09721203 
    ## Run 454 stress 0.09407988 
    ## Run 455 stress 0.09321705 
    ## Run 456 stress 0.0850346 
    ## Run 457 stress 0.09381846 
    ## Run 458 stress 0.3189848 
    ## Run 459 stress 0.2465515 
    ## Run 460 stress 0.09535487 
    ## Run 461 stress 0.08773471 
    ## Run 462 stress 0.08973885 
    ## Run 463 stress 0.09969971 
    ## Run 464 stress 0.09407982 
    ## Run 465 stress 0.09308951 
    ## Run 466 stress 0.09168941 
    ## Run 467 stress 0.09407984 
    ## Run 468 stress 0.08773471 
    ## Run 469 stress 0.09145322 
    ## Run 470 stress 0.08973869 
    ## Run 471 stress 0.08973892 
    ## Run 472 stress 0.1054542 
    ## Run 473 stress 0.09159093 
    ## Run 474 stress 0.08503526 
    ## Run 475 stress 0.08440257 
    ## ... Procrustes: rmse 0.0001748164  max resid 0.0003606606 
    ## ... Similar to previous best
    ## Run 476 stress 0.09535502 
    ## Run 477 stress 0.1038039 
    ## Run 478 stress 0.09381846 
    ## Run 479 stress 0.09321902 
    ## Run 480 stress 0.1053435 
    ## Run 481 stress 0.09030399 
    ## Run 482 stress 0.09030396 
    ## Run 483 stress 0.08503582 
    ## Run 484 stress 0.09030402 
    ## Run 485 stress 0.0877347 
    ## Run 486 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 7.385426e-05  max resid 0.0001388692 
    ## ... Similar to previous best
    ## Run 487 stress 0.09159084 
    ## Run 488 stress 0.09447325 
    ## Run 489 stress 0.09380579 
    ## Run 490 stress 0.09416147 
    ## Run 491 stress 0.08503459 
    ## Run 492 stress 0.09400536 
    ## Run 493 stress 0.09407977 
    ## Run 494 stress 0.09465908 
    ## Run 495 stress 0.08503484 
    ## Run 496 stress 0.09030395 
    ## Run 497 stress 0.100461 
    ## Run 498 stress 0.09969973 
    ## Run 499 stress 0.0897387 
    ## Run 500 stress 0.09159085 
    ## *** Best solution repeated 1 times

``` r
# Ocean sites and mixed lakes
SD_beta_OM_NMDS <- metaMDS(SD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07365783 
    ## Run 1 stress 0.08000467 
    ## Run 2 stress 0.08288044 
    ## Run 3 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001091344  max resid 0.0002601108 
    ## ... Similar to previous best
    ## Run 4 stress 0.08288058 
    ## Run 5 stress 0.07629243 
    ## Run 6 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 1.370736e-05  max resid 3.187845e-05 
    ## ... Similar to previous best
    ## Run 7 stress 0.0800047 
    ## Run 8 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744728  max resid 0.05102681 
    ## Run 9 stress 0.07629246 
    ## Run 10 stress 0.3276486 
    ## Run 11 stress 0.0762925 
    ## Run 12 stress 0.08000468 
    ## Run 13 stress 0.07629235 
    ## Run 14 stress 0.08000468 
    ## Run 15 stress 0.08000468 
    ## Run 16 stress 0.08000467 
    ## Run 17 stress 0.07365783 
    ## ... Procrustes: rmse 3.453399e-05  max resid 8.18353e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.08288038 
    ## Run 19 stress 0.0828806 
    ## Run 20 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174466  max resid 0.05103085 
    ## Run 21 stress 0.08233773 
    ## Run 22 stress 0.07629235 
    ## Run 23 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 5.618748e-06  max resid 1.210869e-05 
    ## ... Similar to previous best
    ## Run 24 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745083  max resid 0.05105097 
    ## Run 25 stress 0.08000468 
    ## Run 26 stress 0.08233766 
    ## Run 27 stress 0.07365785 
    ## ... Procrustes: rmse 8.573822e-05  max resid 0.0001999116 
    ## ... Similar to previous best
    ## Run 28 stress 0.07629236 
    ## Run 29 stress 0.07365783 
    ## ... Procrustes: rmse 1.345985e-05  max resid 2.862139e-05 
    ## ... Similar to previous best
    ## Run 30 stress 0.07629243 
    ## Run 31 stress 0.07629234 
    ## Run 32 stress 0.07629234 
    ## Run 33 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744682  max resid 0.05104917 
    ## Run 34 stress 0.08288038 
    ## Run 35 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001717584  max resid 0.0004047696 
    ## ... Similar to previous best
    ## Run 36 stress 0.07732925 
    ## Run 37 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001105003  max resid 0.0002572499 
    ## ... Similar to previous best
    ## Run 38 stress 0.07365785 
    ## ... Procrustes: rmse 5.309745e-05  max resid 0.0001176011 
    ## ... Similar to previous best
    ## Run 39 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744843  max resid 0.05105982 
    ## Run 40 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744493  max resid 0.05102877 
    ## Run 41 stress 0.07365786 
    ## ... Procrustes: rmse 9.323741e-05  max resid 0.0002189651 
    ## ... Similar to previous best
    ## Run 42 stress 0.08000469 
    ## Run 43 stress 0.08000467 
    ## Run 44 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744513  max resid 0.05104058 
    ## Run 45 stress 0.07629233 
    ## Run 46 stress 0.08288034 
    ## Run 47 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744629  max resid 0.05104909 
    ## Run 48 stress 0.08000469 
    ## Run 49 stress 0.07732936 
    ## Run 50 stress 0.07365784 
    ## ... Procrustes: rmse 5.757358e-05  max resid 0.0001318881 
    ## ... Similar to previous best
    ## Run 51 stress 0.07732927 
    ## Run 52 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743068  max resid 0.05098627 
    ## Run 53 stress 0.07629246 
    ## Run 54 stress 0.0800047 
    ## Run 55 stress 0.08000469 
    ## Run 56 stress 0.07629234 
    ## Run 57 stress 0.07629239 
    ## Run 58 stress 0.07365783 
    ## ... Procrustes: rmse 1.538577e-05  max resid 3.612461e-05 
    ## ... Similar to previous best
    ## Run 59 stress 0.07732928 
    ## Run 60 stress 0.08000471 
    ## Run 61 stress 0.07365785 
    ## ... Procrustes: rmse 5.356784e-05  max resid 0.0001258532 
    ## ... Similar to previous best
    ## Run 62 stress 0.08000467 
    ## Run 63 stress 0.07629233 
    ## Run 64 stress 0.07629236 
    ## Run 65 stress 0.07365784 
    ## ... Procrustes: rmse 2.551245e-05  max resid 4.662727e-05 
    ## ... Similar to previous best
    ## Run 66 stress 0.07365783 
    ## ... Procrustes: rmse 1.317678e-05  max resid 3.098574e-05 
    ## ... Similar to previous best
    ## Run 67 stress 0.07629235 
    ## Run 68 stress 0.07629235 
    ## Run 69 stress 0.08000467 
    ## Run 70 stress 0.08288037 
    ## Run 71 stress 0.07365783 
    ## ... Procrustes: rmse 1.589459e-05  max resid 3.304047e-05 
    ## ... Similar to previous best
    ## Run 72 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001021807  max resid 0.0002368045 
    ## ... Similar to previous best
    ## Run 73 stress 0.07629243 
    ## Run 74 stress 0.08000467 
    ## Run 75 stress 0.349555 
    ## Run 76 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744415  max resid 0.05103626 
    ## Run 77 stress 0.07732932 
    ## Run 78 stress 0.07629245 
    ## Run 79 stress 0.07629243 
    ## Run 80 stress 0.07365783 
    ## ... Procrustes: rmse 3.327297e-05  max resid 7.656093e-05 
    ## ... Similar to previous best
    ## Run 81 stress 0.07732927 
    ## Run 82 stress 0.07732927 
    ## Run 83 stress 0.08000468 
    ## Run 84 stress 0.07365783 
    ## ... Procrustes: rmse 3.172185e-05  max resid 7.279963e-05 
    ## ... Similar to previous best
    ## Run 85 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001582592  max resid 0.0003739537 
    ## ... Similar to previous best
    ## Run 86 stress 0.0737823 
    ## ... Procrustes: rmse 0.01747339  max resid 0.05118562 
    ## Run 87 stress 0.07365793 
    ## ... Procrustes: rmse 0.0001462142  max resid 0.0003310081 
    ## ... Similar to previous best
    ## Run 88 stress 0.07629245 
    ## Run 89 stress 0.07629243 
    ## Run 90 stress 0.07629239 
    ## Run 91 stress 0.08000469 
    ## Run 92 stress 0.0737823 
    ## ... Procrustes: rmse 0.01745117  max resid 0.05109497 
    ## Run 93 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 1.022551e-05  max resid 2.137601e-05 
    ## ... Similar to previous best
    ## Run 94 stress 0.08233762 
    ## Run 95 stress 0.07629241 
    ## Run 96 stress 0.07629233 
    ## Run 97 stress 0.07378232 
    ## ... Procrustes: rmse 0.01743334  max resid 0.05102968 
    ## Run 98 stress 0.07378234 
    ## ... Procrustes: rmse 0.01742189  max resid 0.05094996 
    ## Run 99 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744313  max resid 0.05102028 
    ## Run 100 stress 0.08233756 
    ## Run 101 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744066  max resid 0.05102181 
    ## Run 102 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744484  max resid 0.0510677 
    ## Run 103 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744562  max resid 0.05105331 
    ## Run 104 stress 0.08000474 
    ## Run 105 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174444  max resid 0.05105356 
    ## Run 106 stress 0.07365784 
    ## ... Procrustes: rmse 4.008544e-05  max resid 9.538933e-05 
    ## ... Similar to previous best
    ## Run 107 stress 0.08233765 
    ## Run 108 stress 0.07365787 
    ## ... Procrustes: rmse 9.791232e-05  max resid 0.0002302468 
    ## ... Similar to previous best
    ## Run 109 stress 0.08233769 
    ## Run 110 stress 0.07732926 
    ## Run 111 stress 0.07365785 
    ## ... Procrustes: rmse 6.529155e-05  max resid 0.0001545043 
    ## ... Similar to previous best
    ## Run 112 stress 0.07365783 
    ## ... Procrustes: rmse 6.773623e-06  max resid 1.517415e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.08000468 
    ## Run 114 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174446  max resid 0.0510315 
    ## Run 115 stress 0.07732926 
    ## Run 116 stress 0.07365784 
    ## ... Procrustes: rmse 6.693231e-05  max resid 0.000157862 
    ## ... Similar to previous best
    ## Run 117 stress 0.07365783 
    ## ... Procrustes: rmse 1.470484e-05  max resid 3.670931e-05 
    ## ... Similar to previous best
    ## Run 118 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174446  max resid 0.05105253 
    ## Run 119 stress 0.07629237 
    ## Run 120 stress 0.07365784 
    ## ... Procrustes: rmse 2.065556e-05  max resid 4.685167e-05 
    ## ... Similar to previous best
    ## Run 121 stress 0.07365784 
    ## ... Procrustes: rmse 9.607618e-06  max resid 1.992075e-05 
    ## ... Similar to previous best
    ## Run 122 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744495  max resid 0.05105035 
    ## Run 123 stress 0.07378233 
    ## ... Procrustes: rmse 0.01744855  max resid 0.05102036 
    ## Run 124 stress 0.08000468 
    ## Run 125 stress 0.07365786 
    ## ... Procrustes: rmse 8.42819e-05  max resid 0.0001997761 
    ## ... Similar to previous best
    ## Run 126 stress 0.07365783 
    ## ... Procrustes: rmse 1.891556e-05  max resid 4.419093e-05 
    ## ... Similar to previous best
    ## Run 127 stress 0.08000472 
    ## Run 128 stress 0.07732925 
    ## Run 129 stress 0.07629243 
    ## Run 130 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001324946  max resid 0.0003145326 
    ## ... Similar to previous best
    ## Run 131 stress 0.08000472 
    ## Run 132 stress 0.07629235 
    ## Run 133 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744123  max resid 0.05101317 
    ## Run 134 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744287  max resid 0.05105241 
    ## Run 135 stress 0.07629235 
    ## Run 136 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001430388  max resid 0.0003368462 
    ## ... Similar to previous best
    ## Run 137 stress 0.2646835 
    ## Run 138 stress 0.08288056 
    ## Run 139 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744925  max resid 0.051054 
    ## Run 140 stress 0.07629255 
    ## Run 141 stress 0.07365787 
    ## ... Procrustes: rmse 9.598805e-05  max resid 0.0002254587 
    ## ... Similar to previous best
    ## Run 142 stress 0.07629238 
    ## Run 143 stress 0.07629247 
    ## Run 144 stress 0.08000469 
    ## Run 145 stress 0.2560005 
    ## Run 146 stress 0.07365785 
    ## ... Procrustes: rmse 6.589844e-05  max resid 0.0001557171 
    ## ... Similar to previous best
    ## Run 147 stress 0.07629237 
    ## Run 148 stress 0.08288039 
    ## Run 149 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001311573  max resid 0.0003107562 
    ## ... Similar to previous best
    ## Run 150 stress 0.0762924 
    ## Run 151 stress 0.07365783 
    ## ... Procrustes: rmse 3.245758e-05  max resid 7.44437e-05 
    ## ... Similar to previous best
    ## Run 152 stress 0.08000473 
    ## Run 153 stress 0.07365784 
    ## ... Procrustes: rmse 6.575441e-05  max resid 0.0001564328 
    ## ... Similar to previous best
    ## Run 154 stress 0.07365783 
    ## ... Procrustes: rmse 2.373239e-05  max resid 5.443369e-05 
    ## ... Similar to previous best
    ## Run 155 stress 0.08000468 
    ## Run 156 stress 0.07365785 
    ## ... Procrustes: rmse 5.733328e-05  max resid 0.0001340467 
    ## ... Similar to previous best
    ## Run 157 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744292  max resid 0.051053 
    ## Run 158 stress 0.07365783 
    ## ... Procrustes: rmse 1.339931e-05  max resid 3.157573e-05 
    ## ... Similar to previous best
    ## Run 159 stress 0.07629234 
    ## Run 160 stress 0.08000469 
    ## Run 161 stress 0.07629238 
    ## Run 162 stress 0.07629234 
    ## Run 163 stress 0.08000467 
    ## Run 164 stress 0.07629251 
    ## Run 165 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745539  max resid 0.05110626 
    ## Run 166 stress 0.0823377 
    ## Run 167 stress 0.07378233 
    ## ... Procrustes: rmse 0.01743825  max resid 0.05097164 
    ## Run 168 stress 0.08288053 
    ## Run 169 stress 0.08000467 
    ## Run 170 stress 0.07365784 
    ## ... Procrustes: rmse 1.759093e-05  max resid 3.191277e-05 
    ## ... Similar to previous best
    ## Run 171 stress 0.07629233 
    ## Run 172 stress 0.07365784 
    ## ... Procrustes: rmse 4.649635e-05  max resid 0.0001043407 
    ## ... Similar to previous best
    ## Run 173 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743718  max resid 0.05100813 
    ## Run 174 stress 0.07732927 
    ## Run 175 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744335  max resid 0.05105406 
    ## Run 176 stress 0.07732934 
    ## Run 177 stress 0.07629243 
    ## Run 178 stress 0.07629234 
    ## Run 179 stress 0.07629235 
    ## Run 180 stress 0.08233755 
    ## Run 181 stress 0.0736579 
    ## ... Procrustes: rmse 3.244334e-05  max resid 6.609031e-05 
    ## ... Similar to previous best
    ## Run 182 stress 0.07629236 
    ## Run 183 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743755  max resid 0.05101323 
    ## Run 184 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001495537  max resid 0.0003530252 
    ## ... Similar to previous best
    ## Run 185 stress 0.08233757 
    ## Run 186 stress 0.07365785 
    ## ... Procrustes: rmse 6.099743e-05  max resid 0.0001446747 
    ## ... Similar to previous best
    ## Run 187 stress 0.08233759 
    ## Run 188 stress 0.07365783 
    ## ... Procrustes: rmse 2.445309e-05  max resid 5.681375e-05 
    ## ... Similar to previous best
    ## Run 189 stress 0.07629235 
    ## Run 190 stress 0.07629234 
    ## Run 191 stress 0.2437516 
    ## Run 192 stress 0.07365784 
    ## ... Procrustes: rmse 5.18832e-05  max resid 0.0001221437 
    ## ... Similar to previous best
    ## Run 193 stress 0.08000476 
    ## Run 194 stress 0.07365785 
    ## ... Procrustes: rmse 5.582458e-05  max resid 0.0001302387 
    ## ... Similar to previous best
    ## Run 195 stress 0.07365784 
    ## ... Procrustes: rmse 3.241757e-05  max resid 7.685273e-05 
    ## ... Similar to previous best
    ## Run 196 stress 0.07378231 
    ## ... Procrustes: rmse 0.01745741  max resid 0.05105979 
    ## Run 197 stress 0.08000471 
    ## Run 198 stress 0.07365783 
    ## ... Procrustes: rmse 2.520151e-05  max resid 5.929424e-05 
    ## ... Similar to previous best
    ## Run 199 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 1.88732e-06  max resid 3.185249e-06 
    ## ... Similar to previous best
    ## Run 200 stress 0.07365795 
    ## ... Procrustes: rmse 0.0001742625  max resid 0.0004148773 
    ## ... Similar to previous best
    ## Run 201 stress 0.07365783 
    ## ... Procrustes: rmse 1.507217e-05  max resid 3.653573e-05 
    ## ... Similar to previous best
    ## Run 202 stress 0.07365785 
    ## ... Procrustes: rmse 7.046546e-05  max resid 0.0001602527 
    ## ... Similar to previous best
    ## Run 203 stress 0.07629239 
    ## Run 204 stress 0.08233757 
    ## Run 205 stress 0.07732926 
    ## Run 206 stress 0.07629239 
    ## Run 207 stress 0.08233765 
    ## Run 208 stress 0.0823375 
    ## Run 209 stress 0.07365784 
    ## ... Procrustes: rmse 1.385056e-05  max resid 2.569625e-05 
    ## ... Similar to previous best
    ## Run 210 stress 0.07629246 
    ## Run 211 stress 0.3578485 
    ## Run 212 stress 0.07629238 
    ## Run 213 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174451  max resid 0.05102517 
    ## Run 214 stress 0.08000474 
    ## Run 215 stress 0.08233762 
    ## Run 216 stress 0.07629237 
    ## Run 217 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174511  max resid 0.05107665 
    ## Run 218 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001135122  max resid 0.0002702802 
    ## ... Similar to previous best
    ## Run 219 stress 0.08000468 
    ## Run 220 stress 0.08288049 
    ## Run 221 stress 0.3423054 
    ## Run 222 stress 0.08000469 
    ## Run 223 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744037  max resid 0.05102339 
    ## Run 224 stress 0.07365783 
    ## ... Procrustes: rmse 5.381907e-06  max resid 8.048042e-06 
    ## ... Similar to previous best
    ## Run 225 stress 0.0773294 
    ## Run 226 stress 0.07629236 
    ## Run 227 stress 0.08000469 
    ## Run 228 stress 0.08288035 
    ## Run 229 stress 0.07629241 
    ## Run 230 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001095632  max resid 0.0002588063 
    ## ... Similar to previous best
    ## Run 231 stress 0.07629234 
    ## Run 232 stress 0.0828805 
    ## Run 233 stress 0.07629247 
    ## Run 234 stress 0.07732935 
    ## Run 235 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001186014  max resid 0.0002783567 
    ## ... Similar to previous best
    ## Run 236 stress 0.07629236 
    ## Run 237 stress 0.08000471 
    ## Run 238 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174432  max resid 0.05102977 
    ## Run 239 stress 0.08233765 
    ## Run 240 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744495  max resid 0.05102525 
    ## Run 241 stress 0.07629243 
    ## Run 242 stress 0.08233767 
    ## Run 243 stress 0.2409739 
    ## Run 244 stress 0.07629233 
    ## Run 245 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001387482  max resid 0.0003264454 
    ## ... Similar to previous best
    ## Run 246 stress 0.07629233 
    ## Run 247 stress 0.2676493 
    ## Run 248 stress 0.07629234 
    ## Run 249 stress 0.07365784 
    ## ... Procrustes: rmse 2.493136e-05  max resid 4.672805e-05 
    ## ... Similar to previous best
    ## Run 250 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001603624  max resid 0.0003732034 
    ## ... Similar to previous best
    ## Run 251 stress 0.08233766 
    ## Run 252 stress 0.07629238 
    ## Run 253 stress 0.07629233 
    ## Run 254 stress 0.07629235 
    ## Run 255 stress 0.07629242 
    ## Run 256 stress 0.07629234 
    ## Run 257 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001525059  max resid 0.0003599601 
    ## ... Similar to previous best
    ## Run 258 stress 0.07365793 
    ## ... Procrustes: rmse 0.000161333  max resid 0.0003840382 
    ## ... Similar to previous best
    ## Run 259 stress 0.08000468 
    ## Run 260 stress 0.07629236 
    ## Run 261 stress 0.07365784 
    ## ... Procrustes: rmse 3.252471e-05  max resid 7.54208e-05 
    ## ... Similar to previous best
    ## Run 262 stress 0.0800047 
    ## Run 263 stress 0.08000471 
    ## Run 264 stress 0.08233765 
    ## Run 265 stress 0.08000469 
    ## Run 266 stress 0.07629247 
    ## Run 267 stress 0.07365784 
    ## ... Procrustes: rmse 5.318993e-05  max resid 0.0001234066 
    ## ... Similar to previous best
    ## Run 268 stress 0.07732932 
    ## Run 269 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744856  max resid 0.05105133 
    ## Run 270 stress 0.07365786 
    ## ... Procrustes: rmse 8.270862e-05  max resid 0.0001951839 
    ## ... Similar to previous best
    ## Run 271 stress 0.0762924 
    ## Run 272 stress 0.07365783 
    ## ... Procrustes: rmse 3.642103e-05  max resid 8.609158e-05 
    ## ... Similar to previous best
    ## Run 273 stress 0.07629243 
    ## Run 274 stress 0.07732925 
    ## Run 275 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744481  max resid 0.0510354 
    ## Run 276 stress 0.07365783 
    ## ... Procrustes: rmse 1.171828e-05  max resid 2.567315e-05 
    ## ... Similar to previous best
    ## Run 277 stress 0.07629237 
    ## Run 278 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745176  max resid 0.05107798 
    ## Run 279 stress 0.07365784 
    ## ... Procrustes: rmse 4.472433e-05  max resid 0.0001039768 
    ## ... Similar to previous best
    ## Run 280 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001371231  max resid 0.000324355 
    ## ... Similar to previous best
    ## Run 281 stress 0.08288046 
    ## Run 282 stress 0.08288039 
    ## Run 283 stress 0.07629233 
    ## Run 284 stress 0.07629236 
    ## Run 285 stress 0.07365784 
    ## ... Procrustes: rmse 5.080041e-05  max resid 0.0001214531 
    ## ... Similar to previous best
    ## Run 286 stress 0.07629236 
    ## Run 287 stress 0.08000473 
    ## Run 288 stress 0.07365783 
    ## ... Procrustes: rmse 3.007028e-05  max resid 7.146449e-05 
    ## ... Similar to previous best
    ## Run 289 stress 0.07365783 
    ## ... Procrustes: rmse 1.538023e-05  max resid 3.510653e-05 
    ## ... Similar to previous best
    ## Run 290 stress 0.07629235 
    ## Run 291 stress 0.07629241 
    ## Run 292 stress 0.08288054 
    ## Run 293 stress 0.07629236 
    ## Run 294 stress 0.07629238 
    ## Run 295 stress 0.08000468 
    ## Run 296 stress 0.08000468 
    ## Run 297 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745641  max resid 0.05110275 
    ## Run 298 stress 0.07629243 
    ## Run 299 stress 0.07365785 
    ## ... Procrustes: rmse 6.636519e-05  max resid 0.0001559057 
    ## ... Similar to previous best
    ## Run 300 stress 0.08288061 
    ## Run 301 stress 0.08000468 
    ## Run 302 stress 0.3223601 
    ## Run 303 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744716  max resid 0.05105272 
    ## Run 304 stress 0.08000467 
    ## Run 305 stress 0.08000474 
    ## Run 306 stress 0.07378227 
    ## ... Procrustes: rmse 0.01742846  max resid 0.05096839 
    ## Run 307 stress 0.07629238 
    ## Run 308 stress 0.08000472 
    ## Run 309 stress 0.07629234 
    ## Run 310 stress 0.08000467 
    ## Run 311 stress 0.08000468 
    ## Run 312 stress 0.07629236 
    ## Run 313 stress 0.07629235 
    ## Run 314 stress 0.0736579 
    ## ... Procrustes: rmse 0.000146441  max resid 0.0003441551 
    ## ... Similar to previous best
    ## Run 315 stress 0.3407909 
    ## Run 316 stress 0.07365784 
    ## ... Procrustes: rmse 4.482132e-05  max resid 0.0001055348 
    ## ... Similar to previous best
    ## Run 317 stress 0.08000476 
    ## Run 318 stress 0.08000468 
    ## Run 319 stress 0.07629237 
    ## Run 320 stress 0.08288037 
    ## Run 321 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744315  max resid 0.05105344 
    ## Run 322 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744891  max resid 0.05108161 
    ## Run 323 stress 0.07629241 
    ## Run 324 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 5.569201e-06  max resid 1.238291e-05 
    ## ... Similar to previous best
    ## Run 325 stress 0.07365784 
    ## ... Procrustes: rmse 3.619612e-05  max resid 8.649452e-05 
    ## ... Similar to previous best
    ## Run 326 stress 0.07629235 
    ## Run 327 stress 0.08000472 
    ## Run 328 stress 0.08288049 
    ## Run 329 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744683  max resid 0.05104823 
    ## Run 330 stress 0.0800047 
    ## Run 331 stress 0.07365787 
    ## ... Procrustes: rmse 9.780664e-05  max resid 0.0002335349 
    ## ... Similar to previous best
    ## Run 332 stress 0.08000476 
    ## Run 333 stress 0.07365784 
    ## ... Procrustes: rmse 4.031167e-05  max resid 9.57882e-05 
    ## ... Similar to previous best
    ## Run 334 stress 0.07365783 
    ## ... Procrustes: rmse 3.211741e-05  max resid 7.173338e-05 
    ## ... Similar to previous best
    ## Run 335 stress 0.07365784 
    ## ... Procrustes: rmse 4.598846e-05  max resid 0.0001076102 
    ## ... Similar to previous best
    ## Run 336 stress 0.07629246 
    ## Run 337 stress 0.08000474 
    ## Run 338 stress 0.07629234 
    ## Run 339 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001074802  max resid 0.0002524672 
    ## ... Similar to previous best
    ## Run 340 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001087159  max resid 0.0002531012 
    ## ... Similar to previous best
    ## Run 341 stress 0.07732932 
    ## Run 342 stress 0.07629237 
    ## Run 343 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744763  max resid 0.0510568 
    ## Run 344 stress 0.07629236 
    ## Run 345 stress 0.08288054 
    ## Run 346 stress 0.08000473 
    ## Run 347 stress 0.07629234 
    ## Run 348 stress 0.07629243 
    ## Run 349 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744586  max resid 0.05104084 
    ## Run 350 stress 0.07365788 
    ## ... Procrustes: rmse 4.764154e-05  max resid 0.0001060923 
    ## ... Similar to previous best
    ## Run 351 stress 0.07365787 
    ## ... Procrustes: rmse 8.562806e-05  max resid 0.0001977153 
    ## ... Similar to previous best
    ## Run 352 stress 0.07365786 
    ## ... Procrustes: rmse 8.316334e-05  max resid 0.0001966609 
    ## ... Similar to previous best
    ## Run 353 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744567  max resid 0.05105839 
    ## Run 354 stress 0.07629235 
    ## Run 355 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744422  max resid 0.05102328 
    ## Run 356 stress 0.08000473 
    ## Run 357 stress 0.07365784 
    ## ... Procrustes: rmse 5.831449e-05  max resid 0.0001372503 
    ## ... Similar to previous best
    ## Run 358 stress 0.07365783 
    ## ... Procrustes: rmse 4.945539e-06  max resid 9.220576e-06 
    ## ... Similar to previous best
    ## Run 359 stress 0.08000471 
    ## Run 360 stress 0.07365786 
    ## ... Procrustes: rmse 9.014841e-05  max resid 0.0002089426 
    ## ... Similar to previous best
    ## Run 361 stress 0.08000472 
    ## Run 362 stress 0.07629237 
    ## Run 363 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744692  max resid 0.05102841 
    ## Run 364 stress 0.07365783 
    ## ... Procrustes: rmse 1.879281e-05  max resid 4.419915e-05 
    ## ... Similar to previous best
    ## Run 365 stress 0.07365799 
    ## ... Procrustes: rmse 0.0001880316  max resid 0.0004506763 
    ## ... Similar to previous best
    ## Run 366 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743234  max resid 0.05098796 
    ## Run 367 stress 0.07378233 
    ## ... Procrustes: rmse 0.01744897  max resid 0.05108776 
    ## Run 368 stress 0.07365786 
    ## ... Procrustes: rmse 9.176821e-05  max resid 0.0002147208 
    ## ... Similar to previous best
    ## Run 369 stress 0.07629236 
    ## Run 370 stress 0.08000472 
    ## Run 371 stress 0.07365787 
    ## ... Procrustes: rmse 6.204402e-05  max resid 0.0001449368 
    ## ... Similar to previous best
    ## Run 372 stress 0.07365783 
    ## ... Procrustes: rmse 5.201903e-06  max resid 1.012807e-05 
    ## ... Similar to previous best
    ## Run 373 stress 0.08000468 
    ## Run 374 stress 0.07732926 
    ## Run 375 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001255606  max resid 0.000297781 
    ## ... Similar to previous best
    ## Run 376 stress 0.07629234 
    ## Run 377 stress 0.08000467 
    ## Run 378 stress 0.07365784 
    ## ... Procrustes: rmse 2.347793e-05  max resid 5.459104e-05 
    ## ... Similar to previous best
    ## Run 379 stress 0.08000467 
    ## Run 380 stress 0.07629253 
    ## Run 381 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745928  max resid 0.0511114 
    ## Run 382 stress 0.07732935 
    ## Run 383 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744418  max resid 0.05102102 
    ## Run 384 stress 0.07629252 
    ## Run 385 stress 0.07365783 
    ## ... Procrustes: rmse 2.256607e-05  max resid 5.249225e-05 
    ## ... Similar to previous best
    ## Run 386 stress 0.07629241 
    ## Run 387 stress 0.07378228 
    ## ... Procrustes: rmse 0.01747831  max resid 0.05116507 
    ## Run 388 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745174  max resid 0.05107488 
    ## Run 389 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745056  max resid 0.05106109 
    ## Run 390 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744957  max resid 0.0510673 
    ## Run 391 stress 0.08000468 
    ## Run 392 stress 0.08000471 
    ## Run 393 stress 0.08000467 
    ## Run 394 stress 0.07365786 
    ## ... Procrustes: rmse 9.223942e-05  max resid 0.0002195555 
    ## ... Similar to previous best
    ## Run 395 stress 0.07629241 
    ## Run 396 stress 0.08000471 
    ## Run 397 stress 0.07732932 
    ## Run 398 stress 0.07629248 
    ## Run 399 stress 0.08233762 
    ## Run 400 stress 0.08233757 
    ## Run 401 stress 0.07365785 
    ## ... Procrustes: rmse 7.087686e-05  max resid 0.0001625209 
    ## ... Similar to previous best
    ## Run 402 stress 0.08233761 
    ## Run 403 stress 0.07378228 
    ## ... Procrustes: rmse 0.017445  max resid 0.05105698 
    ## Run 404 stress 0.08000468 
    ## Run 405 stress 0.07629234 
    ## Run 406 stress 0.07629237 
    ## Run 407 stress 0.08000482 
    ## Run 408 stress 0.07732938 
    ## Run 409 stress 0.0737823 
    ## ... Procrustes: rmse 0.0174422  max resid 0.05105616 
    ## Run 410 stress 0.07732948 
    ## Run 411 stress 0.07629233 
    ## Run 412 stress 0.07365785 
    ## ... Procrustes: rmse 5.137873e-05  max resid 0.0001219952 
    ## ... Similar to previous best
    ## Run 413 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744889  max resid 0.05106035 
    ## Run 414 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744407  max resid 0.05103233 
    ## Run 415 stress 0.07629249 
    ## Run 416 stress 0.07629242 
    ## Run 417 stress 0.07732932 
    ## Run 418 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744878  max resid 0.05105728 
    ## Run 419 stress 0.07629245 
    ## Run 420 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745174  max resid 0.05105417 
    ## Run 421 stress 0.08000469 
    ## Run 422 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174416  max resid 0.05101336 
    ## Run 423 stress 0.08000471 
    ## Run 424 stress 0.07629237 
    ## Run 425 stress 0.07732933 
    ## Run 426 stress 0.07629234 
    ## Run 427 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744524  max resid 0.05105812 
    ## Run 428 stress 0.08288059 
    ## Run 429 stress 0.07365786 
    ## ... Procrustes: rmse 8.25939e-05  max resid 0.0001915044 
    ## ... Similar to previous best
    ## Run 430 stress 0.07732929 
    ## Run 431 stress 0.07365783 
    ## ... Procrustes: rmse 7.383003e-06  max resid 1.437888e-05 
    ## ... Similar to previous best
    ## Run 432 stress 0.07365783 
    ## ... Procrustes: rmse 1.499758e-05  max resid 3.404022e-05 
    ## ... Similar to previous best
    ## Run 433 stress 0.08233762 
    ## Run 434 stress 0.07629248 
    ## Run 435 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001076473  max resid 0.0002544495 
    ## ... Similar to previous best
    ## Run 436 stress 0.0773294 
    ## Run 437 stress 0.0762924 
    ## Run 438 stress 0.0762924 
    ## Run 439 stress 0.08233763 
    ## Run 440 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744775  max resid 0.05105795 
    ## Run 441 stress 0.0762924 
    ## Run 442 stress 0.08000469 
    ## Run 443 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001381652  max resid 0.0003261486 
    ## ... Similar to previous best
    ## Run 444 stress 0.07732942 
    ## Run 445 stress 0.07732933 
    ## Run 446 stress 0.08233751 
    ## Run 447 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001161267  max resid 0.0002642695 
    ## ... Similar to previous best
    ## Run 448 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744498  max resid 0.05106138 
    ## Run 449 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744566  max resid 0.05103475 
    ## Run 450 stress 0.0800047 
    ## Run 451 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744475  max resid 0.05105805 
    ## Run 452 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744289  max resid 0.05103021 
    ## Run 453 stress 0.07629239 
    ## Run 454 stress 0.08000469 
    ## Run 455 stress 0.08000467 
    ## Run 456 stress 0.07629236 
    ## Run 457 stress 0.08288042 
    ## Run 458 stress 0.08000471 
    ## Run 459 stress 0.07629236 
    ## Run 460 stress 0.0737823 
    ## ... Procrustes: rmse 0.01746622  max resid 0.05115475 
    ## Run 461 stress 0.07365783 
    ## ... Procrustes: rmse 1.356434e-05  max resid 2.983305e-05 
    ## ... Similar to previous best
    ## Run 462 stress 0.07629238 
    ## Run 463 stress 0.07365785 
    ## ... Procrustes: rmse 7.384671e-05  max resid 0.0001757239 
    ## ... Similar to previous best
    ## Run 464 stress 0.08233763 
    ## Run 465 stress 0.08000467 
    ## Run 466 stress 0.07365784 
    ## ... Procrustes: rmse 5.185599e-05  max resid 0.0001239084 
    ## ... Similar to previous best
    ## Run 467 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745917  max resid 0.05111221 
    ## Run 468 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174358  max resid 0.05100936 
    ## Run 469 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744416  max resid 0.05103767 
    ## Run 470 stress 0.07378232 
    ## ... Procrustes: rmse 0.01744006  max resid 0.05105497 
    ## Run 471 stress 0.08000474 
    ## Run 472 stress 0.08233759 
    ## Run 473 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001491421  max resid 0.000352965 
    ## ... Similar to previous best
    ## Run 474 stress 0.08000469 
    ## Run 475 stress 0.07365785 
    ## ... Procrustes: rmse 6.981224e-05  max resid 0.0001655685 
    ## ... Similar to previous best
    ## Run 476 stress 0.08000481 
    ## Run 477 stress 0.07629235 
    ## Run 478 stress 0.07365786 
    ## ... Procrustes: rmse 8.029701e-05  max resid 0.0001818412 
    ## ... Similar to previous best
    ## Run 479 stress 0.07365785 
    ## ... Procrustes: rmse 7.009473e-05  max resid 0.000164031 
    ## ... Similar to previous best
    ## Run 480 stress 0.07629239 
    ## Run 481 stress 0.07732925 
    ## Run 482 stress 0.07629238 
    ## Run 483 stress 0.07732937 
    ## Run 484 stress 0.07629235 
    ## Run 485 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001374613  max resid 0.0003204308 
    ## ... Similar to previous best
    ## Run 486 stress 0.08233763 
    ## Run 487 stress 0.07365785 
    ## ... Procrustes: rmse 5.838867e-05  max resid 0.0001390478 
    ## ... Similar to previous best
    ## Run 488 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744697  max resid 0.05104871 
    ## Run 489 stress 0.07629241 
    ## Run 490 stress 0.07365784 
    ## ... Procrustes: rmse 5.91017e-05  max resid 0.0001411815 
    ## ... Similar to previous best
    ## Run 491 stress 0.07732936 
    ## Run 492 stress 0.08000468 
    ## Run 493 stress 0.08233763 
    ## Run 494 stress 0.07629239 
    ## Run 495 stress 0.07629233 
    ## Run 496 stress 0.07629245 
    ## Run 497 stress 0.07365783 
    ## ... Procrustes: rmse 1.290688e-05  max resid 2.786295e-05 
    ## ... Similar to previous best
    ## Run 498 stress 0.0800047 
    ## Run 499 stress 0.07629236 
    ## Run 500 stress 0.08000468 
    ## *** Best solution repeated 42 times

``` r
# Stratified lakes and ocean sites
SD_beta_SO_NMDS <- metaMDS(SD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07970526 
    ## Run 1 stress 0.07250812 
    ## ... New best solution
    ## ... Procrustes: rmse 0.08994047  max resid 0.2576279 
    ## Run 2 stress 0.07428313 
    ## Run 3 stress 0.0844846 
    ## Run 4 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04030603  max resid 0.1229625 
    ## Run 5 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321251  max resid 0.03324174 
    ## Run 6 stress 0.06978194 
    ## ... Procrustes: rmse 0.013253  max resid 0.03332816 
    ## Run 7 stress 0.07250815 
    ## Run 8 stress 0.08340287 
    ## Run 9 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 3.090505e-05  max resid 7.965766e-05 
    ## ... Similar to previous best
    ## Run 10 stress 0.07428319 
    ## Run 11 stress 0.06942776 
    ## ... Procrustes: rmse 7.391559e-06  max resid 1.79487e-05 
    ## ... Similar to previous best
    ## Run 12 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324824  max resid 0.03331486 
    ## Run 13 stress 0.07428314 
    ## Run 14 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 4.704297e-06  max resid 1.174329e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.08340287 
    ## Run 16 stress 0.07428321 
    ## Run 17 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325014  max resid 0.03332002 
    ## Run 18 stress 0.08340287 
    ## Run 19 stress 0.08340287 
    ## Run 20 stress 0.07250814 
    ## Run 21 stress 0.06942776 
    ## ... Procrustes: rmse 2.444566e-06  max resid 6.46141e-06 
    ## ... Similar to previous best
    ## Run 22 stress 0.07970525 
    ## Run 23 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324031  max resid 0.03329773 
    ## Run 24 stress 0.07250812 
    ## Run 25 stress 0.07428315 
    ## Run 26 stress 0.0844844 
    ## Run 27 stress 0.06942778 
    ## ... Procrustes: rmse 5.54523e-05  max resid 0.0001426913 
    ## ... Similar to previous best
    ## Run 28 stress 0.07250813 
    ## Run 29 stress 0.07970525 
    ## Run 30 stress 0.06942776 
    ## ... Procrustes: rmse 7.853244e-06  max resid 2.026018e-05 
    ## ... Similar to previous best
    ## Run 31 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326826  max resid 0.03335822 
    ## Run 32 stress 0.08340288 
    ## Run 33 stress 0.08340293 
    ## Run 34 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322877  max resid 0.0332742 
    ## Run 35 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326756  max resid 0.03335337 
    ## Run 36 stress 0.07428313 
    ## Run 37 stress 0.07428314 
    ## Run 38 stress 0.08340294 
    ## Run 39 stress 0.07250812 
    ## Run 40 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321738  max resid 0.03325104 
    ## Run 41 stress 0.08448451 
    ## Run 42 stress 0.07250814 
    ## Run 43 stress 0.07250813 
    ## Run 44 stress 0.08340286 
    ## Run 45 stress 0.07844928 
    ## Run 46 stress 0.08340288 
    ## Run 47 stress 0.07428313 
    ## Run 48 stress 0.08340293 
    ## Run 49 stress 0.07428313 
    ## Run 50 stress 0.07428317 
    ## Run 51 stress 0.069782 
    ## ... Procrustes: rmse 0.01328201  max resid 0.03338607 
    ## Run 52 stress 0.07428318 
    ## Run 53 stress 0.07970527 
    ## Run 54 stress 0.08340286 
    ## Run 55 stress 0.08340288 
    ## Run 56 stress 0.06942776 
    ## ... Procrustes: rmse 2.366709e-06  max resid 6.275173e-06 
    ## ... Similar to previous best
    ## Run 57 stress 0.07428316 
    ## Run 58 stress 0.07428314 
    ## Run 59 stress 0.07970525 
    ## Run 60 stress 0.07428313 
    ## Run 61 stress 0.07250812 
    ## Run 62 stress 0.08340297 
    ## Run 63 stress 0.07250812 
    ## Run 64 stress 0.08448438 
    ## Run 65 stress 0.0844844 
    ## Run 66 stress 0.06978202 
    ## ... Procrustes: rmse 0.01315719  max resid 0.03311902 
    ## Run 67 stress 0.07428321 
    ## Run 68 stress 0.07428313 
    ## Run 69 stress 0.08340286 
    ## Run 70 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324285  max resid 0.03330446 
    ## Run 71 stress 0.07250813 
    ## Run 72 stress 0.07428316 
    ## Run 73 stress 0.06942776 
    ## ... Procrustes: rmse 3.512552e-06  max resid 8.313663e-06 
    ## ... Similar to previous best
    ## Run 74 stress 0.06942776 
    ## ... Procrustes: rmse 3.032536e-05  max resid 7.844959e-05 
    ## ... Similar to previous best
    ## Run 75 stress 0.07970526 
    ## Run 76 stress 0.06942776 
    ## ... Procrustes: rmse 1.933243e-05  max resid 4.985019e-05 
    ## ... Similar to previous best
    ## Run 77 stress 0.07250815 
    ## Run 78 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327329  max resid 0.03336626 
    ## Run 79 stress 0.06978191 
    ## ... Procrustes: rmse 0.013242  max resid 0.03330347 
    ## Run 80 stress 0.06978198 
    ## ... Procrustes: rmse 0.01326772  max resid 0.03335243 
    ## Run 81 stress 0.07250812 
    ## Run 82 stress 0.08448436 
    ## Run 83 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323134  max resid 0.03327993 
    ## Run 84 stress 0.07970526 
    ## Run 85 stress 0.06942776 
    ## ... Procrustes: rmse 2.972864e-05  max resid 7.687812e-05 
    ## ... Similar to previous best
    ## Run 86 stress 0.06978202 
    ## ... Procrustes: rmse 0.01315794  max resid 0.03312029 
    ## Run 87 stress 0.07250813 
    ## Run 88 stress 0.07250815 
    ## Run 89 stress 0.08448438 
    ## Run 90 stress 0.06978192 
    ## ... Procrustes: rmse 0.01323863  max resid 0.03329884 
    ## Run 91 stress 0.08340288 
    ## Run 92 stress 0.07428314 
    ## Run 93 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325017  max resid 0.03331845 
    ## Run 94 stress 0.07428315 
    ## Run 95 stress 0.06942776 
    ## ... Procrustes: rmse 3.539999e-06  max resid 1.139857e-05 
    ## ... Similar to previous best
    ## Run 96 stress 0.07844913 
    ## Run 97 stress 0.07970526 
    ## Run 98 stress 0.08448444 
    ## Run 99 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323357  max resid 0.03328707 
    ## Run 100 stress 0.08340286 
    ## Run 101 stress 0.07428313 
    ## Run 102 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323795  max resid 0.03329553 
    ## Run 103 stress 0.08340298 
    ## Run 104 stress 0.0784495 
    ## Run 105 stress 0.06942777 
    ## ... Procrustes: rmse 3.877704e-05  max resid 0.0001007023 
    ## ... Similar to previous best
    ## Run 106 stress 0.07428313 
    ## Run 107 stress 0.07428315 
    ## Run 108 stress 0.07428313 
    ## Run 109 stress 0.06942777 
    ## ... Procrustes: rmse 3.635735e-05  max resid 9.384094e-05 
    ## ... Similar to previous best
    ## Run 110 stress 0.07250812 
    ## Run 111 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322428  max resid 0.03326579 
    ## Run 112 stress 0.06942776 
    ## ... Procrustes: rmse 1.462715e-05  max resid 3.677122e-05 
    ## ... Similar to previous best
    ## Run 113 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322843  max resid 0.0332733 
    ## Run 114 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320396  max resid 0.03322161 
    ## Run 115 stress 0.06978195 
    ## ... Procrustes: rmse 0.01317678  max resid 0.03316355 
    ## Run 116 stress 0.07428318 
    ## Run 117 stress 0.06942777 
    ## ... Procrustes: rmse 4.727669e-05  max resid 0.0001219058 
    ## ... Similar to previous best
    ## Run 118 stress 0.07970525 
    ## Run 119 stress 0.07428315 
    ## Run 120 stress 0.07250812 
    ## Run 121 stress 0.07428313 
    ## Run 122 stress 0.06942776 
    ## ... Procrustes: rmse 8.100062e-07  max resid 1.372422e-06 
    ## ... Similar to previous best
    ## Run 123 stress 0.08340289 
    ## Run 124 stress 0.07250812 
    ## Run 125 stress 0.07250813 
    ## Run 126 stress 0.07428313 
    ## Run 127 stress 0.08448452 
    ## Run 128 stress 0.06942776 
    ## ... Procrustes: rmse 5.859354e-06  max resid 1.541427e-05 
    ## ... Similar to previous best
    ## Run 129 stress 0.08340287 
    ## Run 130 stress 0.08340288 
    ## Run 131 stress 0.07428321 
    ## Run 132 stress 0.07428313 
    ## Run 133 stress 0.06942776 
    ## ... Procrustes: rmse 2.265011e-05  max resid 5.928105e-05 
    ## ... Similar to previous best
    ## Run 134 stress 0.0697819 
    ## ... Procrustes: rmse 0.01324019  max resid 0.0332981 
    ## Run 135 stress 0.07428313 
    ## Run 136 stress 0.07428317 
    ## Run 137 stress 0.07844948 
    ## Run 138 stress 0.06942776 
    ## ... Procrustes: rmse 2.651667e-05  max resid 6.846731e-05 
    ## ... Similar to previous best
    ## Run 139 stress 0.06942776 
    ## ... Procrustes: rmse 3.613008e-06  max resid 1.015503e-05 
    ## ... Similar to previous best
    ## Run 140 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323456  max resid 0.03328781 
    ## Run 141 stress 0.06942776 
    ## ... Procrustes: rmse 1.080601e-05  max resid 2.791086e-05 
    ## ... Similar to previous best
    ## Run 142 stress 0.06942776 
    ## ... Procrustes: rmse 2.571262e-05  max resid 6.621079e-05 
    ## ... Similar to previous best
    ## Run 143 stress 0.08340296 
    ## Run 144 stress 0.08340289 
    ## Run 145 stress 0.06942776 
    ## ... Procrustes: rmse 1.680345e-05  max resid 4.380247e-05 
    ## ... Similar to previous best
    ## Run 146 stress 0.07428316 
    ## Run 147 stress 0.08340286 
    ## Run 148 stress 0.0742832 
    ## Run 149 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324699  max resid 0.03331223 
    ## Run 150 stress 0.07428315 
    ## Run 151 stress 0.08340287 
    ## Run 152 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326645  max resid 0.03335268 
    ## Run 153 stress 0.07970525 
    ## Run 154 stress 0.08448438 
    ## Run 155 stress 0.07970525 
    ## Run 156 stress 0.07250812 
    ## Run 157 stress 0.06942777 
    ## ... Procrustes: rmse 2.542835e-05  max resid 6.715144e-05 
    ## ... Similar to previous best
    ## Run 158 stress 0.07250812 
    ## Run 159 stress 0.07250812 
    ## Run 160 stress 0.07250813 
    ## Run 161 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320643  max resid 0.03322848 
    ## Run 162 stress 0.07250812 
    ## Run 163 stress 0.07250813 
    ## Run 164 stress 0.07250813 
    ## Run 165 stress 0.07428313 
    ## Run 166 stress 0.07970525 
    ## Run 167 stress 0.06978197 
    ## ... Procrustes: rmse 0.0131805  max resid 0.03316734 
    ## Run 168 stress 0.08340287 
    ## Run 169 stress 0.06942776 
    ## ... Procrustes: rmse 1.825223e-05  max resid 4.731916e-05 
    ## ... Similar to previous best
    ## Run 170 stress 0.07250813 
    ## Run 171 stress 0.07428315 
    ## Run 172 stress 0.07428318 
    ## Run 173 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325112  max resid 0.03332363 
    ## Run 174 stress 0.07250813 
    ## Run 175 stress 0.07428313 
    ## Run 176 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322269  max resid 0.03325966 
    ## Run 177 stress 0.07844958 
    ## Run 178 stress 0.07428313 
    ## Run 179 stress 0.07428315 
    ## Run 180 stress 0.0834029 
    ## Run 181 stress 0.06978193 
    ## ... Procrustes: rmse 0.01324363  max resid 0.03330298 
    ## Run 182 stress 0.07428313 
    ## Run 183 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 3.202851e-06  max resid 7.740828e-06 
    ## ... Similar to previous best
    ## Run 184 stress 0.07970526 
    ## Run 185 stress 0.07250812 
    ## Run 186 stress 0.07428318 
    ## Run 187 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132298  max resid 0.03327478 
    ## Run 188 stress 0.07428316 
    ## Run 189 stress 0.07970525 
    ## Run 190 stress 0.07970528 
    ## Run 191 stress 0.06942776 
    ## ... Procrustes: rmse 1.661021e-05  max resid 4.229427e-05 
    ## ... Similar to previous best
    ## Run 192 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321495  max resid 0.03324546 
    ## Run 193 stress 0.08340287 
    ## Run 194 stress 0.08448446 
    ## Run 195 stress 0.07428314 
    ## Run 196 stress 0.07428317 
    ## Run 197 stress 0.08340288 
    ## Run 198 stress 0.0697819 
    ## ... Procrustes: rmse 0.01320478  max resid 0.03322185 
    ## Run 199 stress 0.06942776 
    ## ... Procrustes: rmse 2.130534e-05  max resid 5.515531e-05 
    ## ... Similar to previous best
    ## Run 200 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325518  max resid 0.03332809 
    ## Run 201 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324795  max resid 0.03331452 
    ## Run 202 stress 0.07428315 
    ## Run 203 stress 0.07250814 
    ## Run 204 stress 0.06942776 
    ## ... Procrustes: rmse 7.29602e-06  max resid 1.859164e-05 
    ## ... Similar to previous best
    ## Run 205 stress 0.07250812 
    ## Run 206 stress 0.07970525 
    ## Run 207 stress 0.07428315 
    ## Run 208 stress 0.07250812 
    ## Run 209 stress 0.07428316 
    ## Run 210 stress 0.07970525 
    ## Run 211 stress 0.08340287 
    ## Run 212 stress 0.07428313 
    ## Run 213 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325008  max resid 0.03331828 
    ## Run 214 stress 0.07970525 
    ## Run 215 stress 0.07428318 
    ## Run 216 stress 0.08340287 
    ## Run 217 stress 0.06942777 
    ## ... Procrustes: rmse 3.560393e-05  max resid 9.158175e-05 
    ## ... Similar to previous best
    ## Run 218 stress 0.06942776 
    ## ... Procrustes: rmse 5.626486e-06  max resid 1.384979e-05 
    ## ... Similar to previous best
    ## Run 219 stress 0.08448438 
    ## Run 220 stress 0.08448437 
    ## Run 221 stress 0.08340294 
    ## Run 222 stress 0.07250812 
    ## Run 223 stress 0.06942776 
    ## ... Procrustes: rmse 1.818325e-06  max resid 4.970952e-06 
    ## ... Similar to previous best
    ## Run 224 stress 0.08340298 
    ## Run 225 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326113  max resid 0.03334537 
    ## Run 226 stress 0.08340294 
    ## Run 227 stress 0.07250812 
    ## Run 228 stress 0.08340287 
    ## Run 229 stress 0.07844948 
    ## Run 230 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326909  max resid 0.03335347 
    ## Run 231 stress 0.07428315 
    ## Run 232 stress 0.07428314 
    ## Run 233 stress 0.06978194 
    ## ... Procrustes: rmse 0.01321471  max resid 0.03324873 
    ## Run 234 stress 0.07428313 
    ## Run 235 stress 0.07970525 
    ## Run 236 stress 0.07428313 
    ## Run 237 stress 0.06978196 
    ## ... Procrustes: rmse 0.01325907  max resid 0.03333381 
    ## Run 238 stress 0.06978195 
    ## ... Procrustes: rmse 0.01317895  max resid 0.03316651 
    ## Run 239 stress 0.07250812 
    ## Run 240 stress 0.06942776 
    ## ... Procrustes: rmse 3.087858e-05  max resid 7.914334e-05 
    ## ... Similar to previous best
    ## Run 241 stress 0.06942776 
    ## ... Procrustes: rmse 7.201468e-06  max resid 2.005744e-05 
    ## ... Similar to previous best
    ## Run 242 stress 0.06978197 
    ## ... Procrustes: rmse 0.01326558  max resid 0.03335519 
    ## Run 243 stress 0.07970526 
    ## Run 244 stress 0.07428315 
    ## Run 245 stress 0.06942776 
    ## ... Procrustes: rmse 1.725435e-05  max resid 4.423411e-05 
    ## ... Similar to previous best
    ## Run 246 stress 0.07970525 
    ## Run 247 stress 0.08340287 
    ## Run 248 stress 0.07970526 
    ## Run 249 stress 0.06942776 
    ## ... Procrustes: rmse 6.236138e-06  max resid 1.137994e-05 
    ## ... Similar to previous best
    ## Run 250 stress 0.07970526 
    ## Run 251 stress 0.07970525 
    ## Run 252 stress 0.07250812 
    ## Run 253 stress 0.07428313 
    ## Run 254 stress 0.07428313 
    ## Run 255 stress 0.07250813 
    ## Run 256 stress 0.07428319 
    ## Run 257 stress 0.06942776 
    ## ... Procrustes: rmse 8.651002e-06  max resid 2.619482e-05 
    ## ... Similar to previous best
    ## Run 258 stress 0.08448447 
    ## Run 259 stress 0.06942776 
    ## ... Procrustes: rmse 1.549185e-05  max resid 3.946245e-05 
    ## ... Similar to previous best
    ## Run 260 stress 0.07970525 
    ## Run 261 stress 0.07428316 
    ## Run 262 stress 0.06942777 
    ## ... Procrustes: rmse 5.440853e-05  max resid 0.000139396 
    ## ... Similar to previous best
    ## Run 263 stress 0.07428313 
    ## Run 264 stress 0.06942776 
    ## ... Procrustes: rmse 2.853643e-05  max resid 7.421005e-05 
    ## ... Similar to previous best
    ## Run 265 stress 0.0742832 
    ## Run 266 stress 0.06942777 
    ## ... Procrustes: rmse 2.854077e-05  max resid 7.360697e-05 
    ## ... Similar to previous best
    ## Run 267 stress 0.07970525 
    ## Run 268 stress 0.07250814 
    ## Run 269 stress 0.08448437 
    ## Run 270 stress 0.0834029 
    ## Run 271 stress 0.08340295 
    ## Run 272 stress 0.07428314 
    ## Run 273 stress 0.07970526 
    ## Run 274 stress 0.07428313 
    ## Run 275 stress 0.07428313 
    ## Run 276 stress 0.07428316 
    ## Run 277 stress 0.06942777 
    ## ... Procrustes: rmse 3.289775e-05  max resid 8.543425e-05 
    ## ... Similar to previous best
    ## Run 278 stress 0.06942776 
    ## ... Procrustes: rmse 9.89414e-06  max resid 2.675983e-05 
    ## ... Similar to previous best
    ## Run 279 stress 0.07428313 
    ## Run 280 stress 0.06978194 
    ## ... Procrustes: rmse 0.01318391  max resid 0.03317804 
    ## Run 281 stress 0.06942781 
    ## ... Procrustes: rmse 3.736868e-05  max resid 0.0001068185 
    ## ... Similar to previous best
    ## Run 282 stress 0.08340289 
    ## Run 283 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326841  max resid 0.03335789 
    ## Run 284 stress 0.07428313 
    ## Run 285 stress 0.08340287 
    ## Run 286 stress 0.08340292 
    ## Run 287 stress 0.07970525 
    ## Run 288 stress 0.06942776 
    ## ... Procrustes: rmse 2.098035e-06  max resid 6.67915e-06 
    ## ... Similar to previous best
    ## Run 289 stress 0.07250812 
    ## Run 290 stress 0.07250817 
    ## Run 291 stress 0.07970525 
    ## Run 292 stress 0.07428315 
    ## Run 293 stress 0.07970525 
    ## Run 294 stress 0.07428315 
    ## Run 295 stress 0.08340291 
    ## Run 296 stress 0.07250814 
    ## Run 297 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325858  max resid 0.03333691 
    ## Run 298 stress 0.07250812 
    ## Run 299 stress 0.07250812 
    ## Run 300 stress 0.07250812 
    ## Run 301 stress 0.08448437 
    ## Run 302 stress 0.07250812 
    ## Run 303 stress 0.08340289 
    ## Run 304 stress 0.08448442 
    ## Run 305 stress 0.06942776 
    ## ... Procrustes: rmse 8.598981e-06  max resid 2.166145e-05 
    ## ... Similar to previous best
    ## Run 306 stress 0.08340287 
    ## Run 307 stress 0.08448442 
    ## Run 308 stress 0.08340288 
    ## Run 309 stress 0.07428314 
    ## Run 310 stress 0.07250812 
    ## Run 311 stress 0.07428313 
    ## Run 312 stress 0.07428314 
    ## Run 313 stress 0.08340289 
    ## Run 314 stress 0.0834029 
    ## Run 315 stress 0.08340288 
    ## Run 316 stress 0.07428315 
    ## Run 317 stress 0.08340296 
    ## Run 318 stress 0.07428317 
    ## Run 319 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325126  max resid 0.03332213 
    ## Run 320 stress 0.07250814 
    ## Run 321 stress 0.07250813 
    ## Run 322 stress 0.08340308 
    ## Run 323 stress 0.0834029 
    ## Run 324 stress 0.07970525 
    ## Run 325 stress 0.0742832 
    ## Run 326 stress 0.07970525 
    ## Run 327 stress 0.08340288 
    ## Run 328 stress 0.06978192 
    ## ... Procrustes: rmse 0.01321409  max resid 0.03324809 
    ## Run 329 stress 0.07970525 
    ## Run 330 stress 0.06942776 
    ## ... Procrustes: rmse 2.76175e-06  max resid 7.383583e-06 
    ## ... Similar to previous best
    ## Run 331 stress 0.07428313 
    ## Run 332 stress 0.07428313 
    ## Run 333 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132398  max resid 0.03329808 
    ## Run 334 stress 0.0834029 
    ## Run 335 stress 0.07428333 
    ## Run 336 stress 0.07970526 
    ## Run 337 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323633  max resid 0.03328916 
    ## Run 338 stress 0.07250813 
    ## Run 339 stress 0.07250812 
    ## Run 340 stress 0.08448449 
    ## Run 341 stress 0.07428313 
    ## Run 342 stress 0.07250812 
    ## Run 343 stress 0.07250813 
    ## Run 344 stress 0.08448439 
    ## Run 345 stress 0.07428313 
    ## Run 346 stress 0.07844936 
    ## Run 347 stress 0.07428315 
    ## Run 348 stress 0.06942776 
    ## ... Procrustes: rmse 1.973606e-05  max resid 5.091303e-05 
    ## ... Similar to previous best
    ## Run 349 stress 0.07428313 
    ## Run 350 stress 0.08340288 
    ## Run 351 stress 0.08340291 
    ## Run 352 stress 0.07428314 
    ## Run 353 stress 0.06942776 
    ## ... Procrustes: rmse 1.764669e-05  max resid 4.498789e-05 
    ## ... Similar to previous best
    ## Run 354 stress 0.06942777 
    ## ... Procrustes: rmse 2.610008e-05  max resid 7.254208e-05 
    ## ... Similar to previous best
    ## Run 355 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326365  max resid 0.03334786 
    ## Run 356 stress 0.07428316 
    ## Run 357 stress 0.06942777 
    ## ... Procrustes: rmse 5.550198e-05  max resid 0.0001425227 
    ## ... Similar to previous best
    ## Run 358 stress 0.07970525 
    ## Run 359 stress 0.07250813 
    ## Run 360 stress 0.08448454 
    ## Run 361 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323373  max resid 0.03328511 
    ## Run 362 stress 0.07428315 
    ## Run 363 stress 0.07428319 
    ## Run 364 stress 0.07970525 
    ## Run 365 stress 0.08340287 
    ## Run 366 stress 0.07428315 
    ## Run 367 stress 0.07250812 
    ## Run 368 stress 0.06942777 
    ## ... Procrustes: rmse 2.215684e-05  max resid 6.18604e-05 
    ## ... Similar to previous best
    ## Run 369 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322157  max resid 0.03326085 
    ## Run 370 stress 0.06942777 
    ## ... Procrustes: rmse 2.007932e-05  max resid 5.297235e-05 
    ## ... Similar to previous best
    ## Run 371 stress 0.08340287 
    ## Run 372 stress 0.07250812 
    ## Run 373 stress 0.08448449 
    ## Run 374 stress 0.07970525 
    ## Run 375 stress 0.08340291 
    ## Run 376 stress 0.07844915 
    ## Run 377 stress 0.07250815 
    ## Run 378 stress 0.07428316 
    ## Run 379 stress 0.08340288 
    ## Run 380 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322096  max resid 0.03326347 
    ## Run 381 stress 0.06978192 
    ## ... Procrustes: rmse 0.01323809  max resid 0.03329871 
    ## Run 382 stress 0.08340287 
    ## Run 383 stress 0.07428317 
    ## Run 384 stress 0.06978204 
    ## ... Procrustes: rmse 0.01329433  max resid 0.03341 
    ## Run 385 stress 0.07428325 
    ## Run 386 stress 0.2584132 
    ## Run 387 stress 0.06942777 
    ## ... Procrustes: rmse 3.669162e-05  max resid 9.416724e-05 
    ## ... Similar to previous best
    ## Run 388 stress 0.06942777 
    ## ... Procrustes: rmse 4.029083e-05  max resid 0.0001039444 
    ## ... Similar to previous best
    ## Run 389 stress 0.07428322 
    ## Run 390 stress 0.08340287 
    ## Run 391 stress 0.07970525 
    ## Run 392 stress 0.06942776 
    ## ... Procrustes: rmse 6.840498e-06  max resid 1.723658e-05 
    ## ... Similar to previous best
    ## Run 393 stress 0.07970526 
    ## Run 394 stress 0.07970526 
    ## Run 395 stress 0.07970526 
    ## Run 396 stress 0.08340288 
    ## Run 397 stress 0.07250812 
    ## Run 398 stress 0.06942776 
    ## ... Procrustes: rmse 3.999514e-06  max resid 1.048961e-05 
    ## ... Similar to previous best
    ## Run 399 stress 0.07428313 
    ## Run 400 stress 0.07428317 
    ## Run 401 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322163  max resid 0.03326785 
    ## Run 402 stress 0.08340294 
    ## Run 403 stress 0.07428319 
    ## Run 404 stress 0.07250812 
    ## Run 405 stress 0.07970525 
    ## Run 406 stress 0.08340291 
    ## Run 407 stress 0.0697819 
    ## ... Procrustes: rmse 0.01324199  max resid 0.03330278 
    ## Run 408 stress 0.07250812 
    ## Run 409 stress 0.06942776 
    ## ... Procrustes: rmse 2.564559e-06  max resid 5.203299e-06 
    ## ... Similar to previous best
    ## Run 410 stress 0.08340292 
    ## Run 411 stress 0.08340297 
    ## Run 412 stress 0.07250812 
    ## Run 413 stress 0.07250813 
    ## Run 414 stress 0.06942776 
    ## ... Procrustes: rmse 1.470095e-05  max resid 3.368668e-05 
    ## ... Similar to previous best
    ## Run 415 stress 0.06942777 
    ## ... Procrustes: rmse 3.689276e-05  max resid 9.472955e-05 
    ## ... Similar to previous best
    ## Run 416 stress 0.08448438 
    ## Run 417 stress 0.08340288 
    ## Run 418 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322672  max resid 0.03327048 
    ## Run 419 stress 0.08340295 
    ## Run 420 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319689  max resid 0.03320709 
    ## Run 421 stress 0.07250813 
    ## Run 422 stress 0.06978192 
    ## ... Procrustes: rmse 0.01320118  max resid 0.03321828 
    ## Run 423 stress 0.07844934 
    ## Run 424 stress 0.07250812 
    ## Run 425 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325185  max resid 0.0333211 
    ## Run 426 stress 0.08340287 
    ## Run 427 stress 0.07970525 
    ## Run 428 stress 0.069782 
    ## ... Procrustes: rmse 0.01316283  max resid 0.03313285 
    ## Run 429 stress 0.07428314 
    ## Run 430 stress 0.07250812 
    ## Run 431 stress 0.07428313 
    ## Run 432 stress 0.07250812 
    ## Run 433 stress 0.08448436 
    ## Run 434 stress 0.07970525 
    ## Run 435 stress 0.07970525 
    ## Run 436 stress 0.06942776 
    ## ... Procrustes: rmse 1.070443e-05  max resid 2.682118e-05 
    ## ... Similar to previous best
    ## Run 437 stress 0.08340287 
    ## Run 438 stress 0.06942776 
    ## ... Procrustes: rmse 9.009876e-06  max resid 2.503405e-05 
    ## ... Similar to previous best
    ## Run 439 stress 0.07250813 
    ## Run 440 stress 0.07250815 
    ## Run 441 stress 0.06942776 
    ## ... Procrustes: rmse 6.038508e-06  max resid 1.518454e-05 
    ## ... Similar to previous best
    ## Run 442 stress 0.08448442 
    ## Run 443 stress 0.07428313 
    ## Run 444 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322068  max resid 0.03325752 
    ## Run 445 stress 0.07428313 
    ## Run 446 stress 0.07970525 
    ## Run 447 stress 0.06942776 
    ## ... Procrustes: rmse 1.287841e-05  max resid 3.294615e-05 
    ## ... Similar to previous best
    ## Run 448 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326557  max resid 0.0333514 
    ## Run 449 stress 0.07428313 
    ## Run 450 stress 0.06942777 
    ## ... Procrustes: rmse 4.44753e-05  max resid 0.0001156722 
    ## ... Similar to previous best
    ## Run 451 stress 0.07428321 
    ## Run 452 stress 0.07250812 
    ## Run 453 stress 0.06978191 
    ## ... Procrustes: rmse 0.0132461  max resid 0.03331118 
    ## Run 454 stress 0.08340288 
    ## Run 455 stress 0.083403 
    ## Run 456 stress 0.08340295 
    ## Run 457 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324244  max resid 0.03330405 
    ## Run 458 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323372  max resid 0.03328536 
    ## Run 459 stress 0.07428317 
    ## Run 460 stress 0.06942778 
    ## ... Procrustes: rmse 1.983506e-05  max resid 3.623515e-05 
    ## ... Similar to previous best
    ## Run 461 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324216  max resid 0.03330473 
    ## Run 462 stress 0.07428316 
    ## Run 463 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326825  max resid 0.03335856 
    ## Run 464 stress 0.06978201 
    ## ... Procrustes: rmse 0.01328489  max resid 0.03338791 
    ## Run 465 stress 0.06978192 
    ## ... Procrustes: rmse 0.0132529  max resid 0.03332459 
    ## Run 466 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323325  max resid 0.0332859 
    ## Run 467 stress 0.08340288 
    ## Run 468 stress 0.08340286 
    ## Run 469 stress 0.06942776 
    ## ... Procrustes: rmse 1.690955e-05  max resid 4.326677e-05 
    ## ... Similar to previous best
    ## Run 470 stress 0.08340287 
    ## Run 471 stress 0.06942776 
    ## ... Procrustes: rmse 2.460061e-05  max resid 6.291579e-05 
    ## ... Similar to previous best
    ## Run 472 stress 0.07250812 
    ## Run 473 stress 0.07428318 
    ## Run 474 stress 0.07428314 
    ## Run 475 stress 0.07250813 
    ## Run 476 stress 0.08340289 
    ## Run 477 stress 0.07844954 
    ## Run 478 stress 0.06942777 
    ## ... Procrustes: rmse 4.930336e-05  max resid 0.0001264323 
    ## ... Similar to previous best
    ## Run 479 stress 0.07428313 
    ## Run 480 stress 0.06978193 
    ## ... Procrustes: rmse 0.01320302  max resid 0.03322347 
    ## Run 481 stress 0.08340294 
    ## Run 482 stress 0.07428313 
    ## Run 483 stress 0.06942778 
    ## ... Procrustes: rmse 5.975864e-05  max resid 0.0001534153 
    ## ... Similar to previous best
    ## Run 484 stress 0.08340292 
    ## Run 485 stress 0.07970525 
    ## Run 486 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325042  max resid 0.03332184 
    ## Run 487 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325043  max resid 0.0333207 
    ## Run 488 stress 0.08340289 
    ## Run 489 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323208  max resid 0.03328227 
    ## Run 490 stress 0.08340287 
    ## Run 491 stress 0.08340299 
    ## Run 492 stress 0.08340287 
    ## Run 493 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325562  max resid 0.0333313 
    ## Run 494 stress 0.07428312 
    ## Run 495 stress 0.07970525 
    ## Run 496 stress 0.07844943 
    ## Run 497 stress 0.07970526 
    ## Run 498 stress 0.08340296 
    ## Run 499 stress 0.08340288 
    ## Run 500 stress 0.06978195 
    ## ... Procrustes: rmse 0.01317888  max resid 0.03316731 
    ## *** Best solution repeated 45 times

``` r
# Stratified lakes
SD_beta_S_NMDS <- metaMDS(SD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.1410449 
    ## Run 1 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 0.2673134  max resid 0.5425959 
    ## Run 2 stress 0.2691259 
    ## Run 3 stress 0.1434197 
    ## Run 4 stress 0.2005461 
    ## Run 5 stress 0.1572668 
    ## Run 6 stress 0.1583016 
    ## Run 7 stress 0.1572668 
    ## Run 8 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 7.791469e-07  max resid 1.506537e-06 
    ## ... Similar to previous best
    ## Run 9 stress 0.1320258 
    ## ... Procrustes: rmse 7.500148e-07  max resid 1.381181e-06 
    ## ... Similar to previous best
    ## Run 10 stress 0.1693717 
    ## Run 11 stress 0.1771969 
    ## Run 12 stress 0.1383682 
    ## Run 13 stress 0.2008162 
    ## Run 14 stress 0.3083098 
    ## Run 15 stress 0.2823279 
    ## Run 16 stress 0.1320258 
    ## ... Procrustes: rmse 2.12401e-06  max resid 4.454918e-06 
    ## ... Similar to previous best
    ## Run 17 stress 0.1410451 
    ## Run 18 stress 0.1410451 
    ## Run 19 stress 0.1970303 
    ## Run 20 stress 0.1320258 
    ## ... Procrustes: rmse 1.837546e-06  max resid 3.855089e-06 
    ## ... Similar to previous best
    ## Run 21 stress 0.1383681 
    ## Run 22 stress 0.1581642 
    ## Run 23 stress 0.1383681 
    ## Run 24 stress 0.1383682 
    ## Run 25 stress 0.2373121 
    ## Run 26 stress 0.1771969 
    ## Run 27 stress 0.1383681 
    ## Run 28 stress 0.1771969 
    ## Run 29 stress 0.1383681 
    ## Run 30 stress 0.1415299 
    ## Run 31 stress 0.1693717 
    ## Run 32 stress 0.1320258 
    ## ... Procrustes: rmse 1.452815e-06  max resid 3.001194e-06 
    ## ... Similar to previous best
    ## Run 33 stress 0.1383681 
    ## Run 34 stress 0.1628606 
    ## Run 35 stress 0.2550044 
    ## Run 36 stress 0.1383681 
    ## Run 37 stress 0.1583016 
    ## Run 38 stress 0.2422742 
    ## Run 39 stress 0.1551863 
    ## Run 40 stress 0.1383681 
    ## Run 41 stress 0.1407298 
    ## Run 42 stress 0.1434197 
    ## Run 43 stress 0.1383681 
    ## Run 44 stress 0.1383681 
    ## Run 45 stress 0.1410448 
    ## Run 46 stress 0.1830332 
    ## Run 47 stress 0.1771969 
    ## Run 48 stress 0.1628606 
    ## Run 49 stress 0.1320258 
    ## ... Procrustes: rmse 8.321668e-07  max resid 1.281928e-06 
    ## ... Similar to previous best
    ## Run 50 stress 0.2385029 
    ## Run 51 stress 0.1410449 
    ## Run 52 stress 0.1320258 
    ## ... Procrustes: rmse 1.110886e-06  max resid 1.856216e-06 
    ## ... Similar to previous best
    ## Run 53 stress 0.1383681 
    ## Run 54 stress 0.1410449 
    ## Run 55 stress 0.1572668 
    ## Run 56 stress 0.1320258 
    ## ... Procrustes: rmse 9.860274e-07  max resid 1.96714e-06 
    ## ... Similar to previous best
    ## Run 57 stress 0.1830332 
    ## Run 58 stress 0.1434197 
    ## Run 59 stress 0.1572668 
    ## Run 60 stress 0.1383681 
    ## Run 61 stress 0.140721 
    ## Run 62 stress 0.1320258 
    ## ... Procrustes: rmse 4.520974e-07  max resid 9.135821e-07 
    ## ... Similar to previous best
    ## Run 63 stress 0.1320258 
    ## ... Procrustes: rmse 9.165274e-07  max resid 1.908712e-06 
    ## ... Similar to previous best
    ## Run 64 stress 0.1320258 
    ## ... Procrustes: rmse 2.041636e-06  max resid 4.205846e-06 
    ## ... Similar to previous best
    ## Run 65 stress 0.1583016 
    ## Run 66 stress 0.1796865 
    ## Run 67 stress 0.1551863 
    ## Run 68 stress 0.1572668 
    ## Run 69 stress 0.1407298 
    ## Run 70 stress 0.1434197 
    ## Run 71 stress 0.2422742 
    ## Run 72 stress 0.1407298 
    ## Run 73 stress 0.2224299 
    ## Run 74 stress 0.1802752 
    ## Run 75 stress 0.1320258 
    ## ... Procrustes: rmse 3.457904e-07  max resid 6.334868e-07 
    ## ... Similar to previous best
    ## Run 76 stress 0.1572668 
    ## Run 77 stress 0.1551863 
    ## Run 78 stress 0.1551863 
    ## Run 79 stress 0.1407298 
    ## Run 80 stress 0.1693717 
    ## Run 81 stress 0.1320258 
    ## ... Procrustes: rmse 2.088948e-06  max resid 4.214561e-06 
    ## ... Similar to previous best
    ## Run 82 stress 0.1320258 
    ## ... Procrustes: rmse 1.673107e-06  max resid 3.41419e-06 
    ## ... Similar to previous best
    ## Run 83 stress 0.1693717 
    ## Run 84 stress 0.2384387 
    ## Run 85 stress 0.1407298 
    ## Run 86 stress 0.1572668 
    ## Run 87 stress 0.1320258 
    ## ... Procrustes: rmse 2.163384e-06  max resid 4.463036e-06 
    ## ... Similar to previous best
    ## Run 88 stress 0.1572668 
    ## Run 89 stress 0.1320258 
    ## ... Procrustes: rmse 1.114541e-06  max resid 2.277257e-06 
    ## ... Similar to previous best
    ## Run 90 stress 0.1410449 
    ## Run 91 stress 0.141045 
    ## Run 92 stress 0.1962476 
    ## Run 93 stress 0.141045 
    ## Run 94 stress 0.1320258 
    ## ... Procrustes: rmse 2.787289e-06  max resid 5.738293e-06 
    ## ... Similar to previous best
    ## Run 95 stress 0.1407298 
    ## Run 96 stress 0.2227526 
    ## Run 97 stress 0.1415299 
    ## Run 98 stress 0.1583016 
    ## Run 99 stress 0.2510079 
    ## Run 100 stress 0.2416997 
    ## Run 101 stress 0.1771966 
    ## Run 102 stress 0.1771969 
    ## Run 103 stress 0.1410448 
    ## Run 104 stress 0.1581642 
    ## Run 105 stress 0.1410456 
    ## Run 106 stress 0.1383681 
    ## Run 107 stress 0.1572668 
    ## Run 108 stress 0.1407208 
    ## Run 109 stress 0.1320258 
    ## ... Procrustes: rmse 1.473954e-06  max resid 2.252175e-06 
    ## ... Similar to previous best
    ## Run 110 stress 0.1551863 
    ## Run 111 stress 0.2008158 
    ## Run 112 stress 0.2384387 
    ## Run 113 stress 0.1415299 
    ## Run 114 stress 0.141045 
    ## Run 115 stress 0.1551863 
    ## Run 116 stress 0.1598154 
    ## Run 117 stress 0.1383682 
    ## Run 118 stress 0.2422667 
    ## Run 119 stress 0.1771969 
    ## Run 120 stress 0.1776793 
    ## Run 121 stress 0.1572668 
    ## Run 122 stress 0.1383681 
    ## Run 123 stress 0.1320258 
    ## ... Procrustes: rmse 1.721608e-06  max resid 2.983876e-06 
    ## ... Similar to previous best
    ## Run 124 stress 0.1407207 
    ## Run 125 stress 0.1434197 
    ## Run 126 stress 0.1320258 
    ## ... Procrustes: rmse 2.580978e-06  max resid 5.326364e-06 
    ## ... Similar to previous best
    ## Run 127 stress 0.1383681 
    ## Run 128 stress 0.1415299 
    ## Run 129 stress 0.2008163 
    ## Run 130 stress 0.2842805 
    ## Run 131 stress 0.1410449 
    ## Run 132 stress 0.200816 
    ## Run 133 stress 0.1407298 
    ## Run 134 stress 0.1320258 
    ## ... Procrustes: rmse 1.107095e-06  max resid 2.307389e-06 
    ## ... Similar to previous best
    ## Run 135 stress 0.1320258 
    ## ... Procrustes: rmse 4.755257e-07  max resid 7.505665e-07 
    ## ... Similar to previous best
    ## Run 136 stress 0.141045 
    ## Run 137 stress 0.1929519 
    ## Run 138 stress 0.1830332 
    ## Run 139 stress 0.1407298 
    ## Run 140 stress 0.140721 
    ## Run 141 stress 0.1572668 
    ## Run 142 stress 0.1383681 
    ## Run 143 stress 0.1407298 
    ## Run 144 stress 0.2385029 
    ## Run 145 stress 0.1410455 
    ## Run 146 stress 0.1407298 
    ## Run 147 stress 0.1693717 
    ## Run 148 stress 0.1383682 
    ## Run 149 stress 0.1383681 
    ## Run 150 stress 0.2373121 
    ## Run 151 stress 0.1415299 
    ## Run 152 stress 0.1693717 
    ## Run 153 stress 0.1383681 
    ## Run 154 stress 0.1693717 
    ## Run 155 stress 0.1383681 
    ## Run 156 stress 0.2385029 
    ## Run 157 stress 0.1410449 
    ## Run 158 stress 0.2578648 
    ## Run 159 stress 0.2074935 
    ## Run 160 stress 0.2578649 
    ## Run 161 stress 0.2578639 
    ## Run 162 stress 0.3012631 
    ## Run 163 stress 0.2550045 
    ## Run 164 stress 0.1415299 
    ## Run 165 stress 0.1383682 
    ## Run 166 stress 0.1572668 
    ## Run 167 stress 0.1407207 
    ## Run 168 stress 0.1693717 
    ## Run 169 stress 0.1693717 
    ## Run 170 stress 0.2422667 
    ## Run 171 stress 0.1415299 
    ## Run 172 stress 0.1407298 
    ## Run 173 stress 0.1415299 
    ## Run 174 stress 0.1383681 
    ## Run 175 stress 0.1771969 
    ## Run 176 stress 0.2629507 
    ## Run 177 stress 0.2422667 
    ## Run 178 stress 0.2496387 
    ## Run 179 stress 0.1551863 
    ## Run 180 stress 0.1693717 
    ## Run 181 stress 0.1407207 
    ## Run 182 stress 0.1407207 
    ## Run 183 stress 0.1383681 
    ## Run 184 stress 0.1407298 
    ## Run 185 stress 0.1628606 
    ## Run 186 stress 0.1970303 
    ## Run 187 stress 0.1693717 
    ## Run 188 stress 0.1771969 
    ## Run 189 stress 0.1796865 
    ## Run 190 stress 0.1628606 
    ## Run 191 stress 0.1551863 
    ## Run 192 stress 0.1410448 
    ## Run 193 stress 0.1796865 
    ## Run 194 stress 0.1410453 
    ## Run 195 stress 0.141045 
    ## Run 196 stress 0.1320258 
    ## ... Procrustes: rmse 1.180805e-06  max resid 2.389015e-06 
    ## ... Similar to previous best
    ## Run 197 stress 0.1410449 
    ## Run 198 stress 0.1320258 
    ## ... Procrustes: rmse 1.707622e-06  max resid 3.490156e-06 
    ## ... Similar to previous best
    ## Run 199 stress 0.1383681 
    ## Run 200 stress 0.1320258 
    ## ... Procrustes: rmse 3.126138e-07  max resid 4.828986e-07 
    ## ... Similar to previous best
    ## Run 201 stress 0.1445811 
    ## Run 202 stress 0.1407206 
    ## Run 203 stress 0.1445811 
    ## Run 204 stress 0.1320258 
    ## ... Procrustes: rmse 7.440516e-07  max resid 1.448049e-06 
    ## ... Similar to previous best
    ## Run 205 stress 0.1929519 
    ## Run 206 stress 0.1415299 
    ## Run 207 stress 0.1407298 
    ## Run 208 stress 0.1410448 
    ## Run 209 stress 0.1410449 
    ## Run 210 stress 0.2842805 
    ## Run 211 stress 0.1407298 
    ## Run 212 stress 0.2416997 
    ## Run 213 stress 0.1410449 
    ## Run 214 stress 0.1693717 
    ## Run 215 stress 0.1320258 
    ## ... Procrustes: rmse 3.075696e-06  max resid 5.444813e-06 
    ## ... Similar to previous best
    ## Run 216 stress 0.1434197 
    ## Run 217 stress 0.1830332 
    ## Run 218 stress 0.1415299 
    ## Run 219 stress 0.1383681 
    ## Run 220 stress 0.2422742 
    ## Run 221 stress 0.1407207 
    ## Run 222 stress 0.1410448 
    ## Run 223 stress 0.1415299 
    ## Run 224 stress 0.1776793 
    ## Run 225 stress 0.1572668 
    ## Run 226 stress 0.1771969 
    ## Run 227 stress 0.1320258 
    ## ... Procrustes: rmse 1.039429e-06  max resid 2.073156e-06 
    ## ... Similar to previous best
    ## Run 228 stress 0.1693717 
    ## Run 229 stress 0.1551863 
    ## Run 230 stress 0.1410452 
    ## Run 231 stress 0.1320258 
    ## ... Procrustes: rmse 2.396807e-06  max resid 4.909756e-06 
    ## ... Similar to previous best
    ## Run 232 stress 0.2385029 
    ## Run 233 stress 0.1320258 
    ## ... Procrustes: rmse 3.613428e-06  max resid 7.2536e-06 
    ## ... Similar to previous best
    ## Run 234 stress 0.1410448 
    ## Run 235 stress 0.2005461 
    ## Run 236 stress 0.1383681 
    ## Run 237 stress 0.1410448 
    ## Run 238 stress 0.1320258 
    ## ... Procrustes: rmse 1.195975e-06  max resid 2.4895e-06 
    ## ... Similar to previous best
    ## Run 239 stress 0.2224299 
    ## Run 240 stress 0.1410449 
    ## Run 241 stress 0.1581642 
    ## Run 242 stress 0.1693717 
    ## Run 243 stress 0.1410449 
    ## Run 244 stress 0.1415299 
    ## Run 245 stress 0.2227526 
    ## Run 246 stress 0.1415299 
    ## Run 247 stress 0.1407298 
    ## Run 248 stress 0.1383682 
    ## Run 249 stress 0.1410449 
    ## Run 250 stress 0.1320258 
    ## ... Procrustes: rmse 2.423971e-06  max resid 4.785728e-06 
    ## ... Similar to previous best
    ## Run 251 stress 0.2224298 
    ## Run 252 stress 0.1383682 
    ## Run 253 stress 0.1320258 
    ## ... Procrustes: rmse 1.564258e-06  max resid 3.16183e-06 
    ## ... Similar to previous best
    ## Run 254 stress 0.1410454 
    ## Run 255 stress 0.1407208 
    ## Run 256 stress 0.1410451 
    ## Run 257 stress 0.1383682 
    ## Run 258 stress 0.1598154 
    ## Run 259 stress 0.1383682 
    ## Run 260 stress 0.1320258 
    ## ... Procrustes: rmse 1.91677e-06  max resid 3.931051e-06 
    ## ... Similar to previous best
    ## Run 261 stress 0.1572668 
    ## Run 262 stress 0.1693717 
    ## Run 263 stress 0.1598154 
    ## Run 264 stress 0.1776793 
    ## Run 265 stress 0.1434197 
    ## Run 266 stress 0.1445811 
    ## Run 267 stress 0.1383681 
    ## Run 268 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 5.475797e-07  max resid 1.043935e-06 
    ## ... Similar to previous best
    ## Run 269 stress 0.1434197 
    ## Run 270 stress 0.1771969 
    ## Run 271 stress 0.1320258 
    ## ... Procrustes: rmse 5.662272e-06  max resid 1.160165e-05 
    ## ... Similar to previous best
    ## Run 272 stress 0.1407298 
    ## Run 273 stress 0.1407298 
    ## Run 274 stress 0.1410449 
    ## Run 275 stress 0.1320258 
    ## ... Procrustes: rmse 1.692586e-06  max resid 2.813925e-06 
    ## ... Similar to previous best
    ## Run 276 stress 0.2538018 
    ## Run 277 stress 0.1407298 
    ## Run 278 stress 0.1551863 
    ## Run 279 stress 0.1551863 
    ## Run 280 stress 0.1598154 
    ## Run 281 stress 0.1572668 
    ## Run 282 stress 0.1383682 
    ## Run 283 stress 0.1415299 
    ## Run 284 stress 0.1410449 
    ## Run 285 stress 0.1320258 
    ## ... Procrustes: rmse 1.817407e-06  max resid 3.706219e-06 
    ## ... Similar to previous best
    ## Run 286 stress 0.1434197 
    ## Run 287 stress 0.1320258 
    ## ... Procrustes: rmse 1.020634e-06  max resid 2.013595e-06 
    ## ... Similar to previous best
    ## Run 288 stress 0.1415299 
    ## Run 289 stress 0.1581642 
    ## Run 290 stress 0.1320258 
    ## ... Procrustes: rmse 8.444996e-07  max resid 1.653848e-06 
    ## ... Similar to previous best
    ## Run 291 stress 0.2501068 
    ## Run 292 stress 0.1434197 
    ## Run 293 stress 0.1415299 
    ## Run 294 stress 0.1771969 
    ## Run 295 stress 0.1445811 
    ## Run 296 stress 0.1663544 
    ## Run 297 stress 0.1771969 
    ## Run 298 stress 0.1587623 
    ## Run 299 stress 0.1407298 
    ## Run 300 stress 0.1320258 
    ## ... Procrustes: rmse 4.068249e-07  max resid 7.404225e-07 
    ## ... Similar to previous best
    ## Run 301 stress 0.1415299 
    ## Run 302 stress 0.1771966 
    ## Run 303 stress 0.2848524 
    ## Run 304 stress 0.1771969 
    ## Run 305 stress 0.1410448 
    ## Run 306 stress 0.1320258 
    ## ... Procrustes: rmse 1.100481e-06  max resid 2.258318e-06 
    ## ... Similar to previous best
    ## Run 307 stress 0.1802752 
    ## Run 308 stress 0.1320258 
    ## ... Procrustes: rmse 1.864261e-06  max resid 4.23512e-06 
    ## ... Similar to previous best
    ## Run 309 stress 0.1320258 
    ## ... Procrustes: rmse 6.455832e-07  max resid 1.287028e-06 
    ## ... Similar to previous best
    ## Run 310 stress 0.2385029 
    ## Run 311 stress 0.1320258 
    ## ... Procrustes: rmse 6.370791e-07  max resid 1.246754e-06 
    ## ... Similar to previous best
    ## Run 312 stress 0.1587623 
    ## Run 313 stress 0.1776793 
    ## Run 314 stress 0.1383682 
    ## Run 315 stress 0.1802752 
    ## Run 316 stress 0.1830332 
    ## Run 317 stress 0.1383682 
    ## Run 318 stress 0.1320258 
    ## ... Procrustes: rmse 1.326517e-06  max resid 2.708936e-06 
    ## ... Similar to previous best
    ## Run 319 stress 0.1415299 
    ## Run 320 stress 0.1320258 
    ## ... Procrustes: rmse 1.86225e-06  max resid 3.861941e-06 
    ## ... Similar to previous best
    ## Run 321 stress 0.2074935 
    ## Run 322 stress 0.1572668 
    ## Run 323 stress 0.2842805 
    ## Run 324 stress 0.1962476 
    ## Run 325 stress 0.1583016 
    ## Run 326 stress 0.1693717 
    ## Run 327 stress 0.1445811 
    ## Run 328 stress 0.1320258 
    ## ... Procrustes: rmse 2.258012e-06  max resid 4.647818e-06 
    ## ... Similar to previous best
    ## Run 329 stress 0.1407298 
    ## Run 330 stress 0.1572668 
    ## Run 331 stress 0.1830332 
    ## Run 332 stress 0.1410453 
    ## Run 333 stress 0.1383681 
    ## Run 334 stress 0.1415299 
    ## Run 335 stress 0.1572668 
    ## Run 336 stress 0.1415299 
    ## Run 337 stress 0.2385029 
    ## Run 338 stress 0.1320258 
    ## ... Procrustes: rmse 2.11822e-06  max resid 4.316135e-06 
    ## ... Similar to previous best
    ## Run 339 stress 0.1693717 
    ## Run 340 stress 0.1410448 
    ## Run 341 stress 0.1320258 
    ## ... Procrustes: rmse 2.857441e-06  max resid 5.811486e-06 
    ## ... Similar to previous best
    ## Run 342 stress 0.1320258 
    ## ... New best solution
    ## ... Procrustes: rmse 1.676795e-07  max resid 2.712936e-07 
    ## ... Similar to previous best
    ## Run 343 stress 0.1415299 
    ## Run 344 stress 0.2005461 
    ## Run 345 stress 0.1434197 
    ## Run 346 stress 0.1415299 
    ## Run 347 stress 0.1383681 
    ## Run 348 stress 0.1383682 
    ## Run 349 stress 0.1383681 
    ## Run 350 stress 0.1410451 
    ## Run 351 stress 0.1434197 
    ## Run 352 stress 0.1415299 
    ## Run 353 stress 0.1383681 
    ## Run 354 stress 0.1383682 
    ## Run 355 stress 0.1383682 
    ## Run 356 stress 0.2015531 
    ## Run 357 stress 0.1320258 
    ## ... Procrustes: rmse 9.968476e-07  max resid 1.822323e-06 
    ## ... Similar to previous best
    ## Run 358 stress 0.1415299 
    ## Run 359 stress 0.1415299 
    ## Run 360 stress 0.1415299 
    ## Run 361 stress 0.2385029 
    ## Run 362 stress 0.1320258 
    ## ... Procrustes: rmse 1.069775e-06  max resid 1.836967e-06 
    ## ... Similar to previous best
    ## Run 363 stress 0.1320258 
    ## ... Procrustes: rmse 7.379649e-07  max resid 1.375764e-06 
    ## ... Similar to previous best
    ## Run 364 stress 0.1320258 
    ## ... Procrustes: rmse 1.941445e-07  max resid 3.317063e-07 
    ## ... Similar to previous best
    ## Run 365 stress 0.1572668 
    ## Run 366 stress 0.1415299 
    ## Run 367 stress 0.1383681 
    ## Run 368 stress 0.1383682 
    ## Run 369 stress 0.2227526 
    ## Run 370 stress 0.1445811 
    ## Run 371 stress 0.1693717 
    ## Run 372 stress 0.2373121 
    ## Run 373 stress 0.1383681 
    ## Run 374 stress 0.1796865 
    ## Run 375 stress 0.1320258 
    ## ... Procrustes: rmse 4.249444e-07  max resid 7.81099e-07 
    ## ... Similar to previous best
    ## Run 376 stress 0.1407207 
    ## Run 377 stress 0.1551863 
    ## Run 378 stress 0.1410449 
    ## Run 379 stress 0.1410452 
    ## Run 380 stress 0.1320258 
    ## ... Procrustes: rmse 1.740233e-06  max resid 3.607669e-06 
    ## ... Similar to previous best
    ## Run 381 stress 0.1410449 
    ## Run 382 stress 0.1383682 
    ## Run 383 stress 0.1410448 
    ## Run 384 stress 0.1407298 
    ## Run 385 stress 0.1598154 
    ## Run 386 stress 0.1771966 
    ## Run 387 stress 0.1383681 
    ## Run 388 stress 0.2711712 
    ## Run 389 stress 0.1407298 
    ## Run 390 stress 0.1962476 
    ## Run 391 stress 0.2015531 
    ## Run 392 stress 0.141045 
    ## Run 393 stress 0.1407298 
    ## Run 394 stress 0.1663544 
    ## Run 395 stress 0.1434197 
    ## Run 396 stress 0.2008163 
    ## Run 397 stress 0.1551863 
    ## Run 398 stress 0.1320258 
    ## ... Procrustes: rmse 9.855065e-07  max resid 1.97491e-06 
    ## ... Similar to previous best
    ## Run 399 stress 0.1383682 
    ## Run 400 stress 0.1410449 
    ## Run 401 stress 0.1407298 
    ## Run 402 stress 0.1970303 
    ## Run 403 stress 0.1802587 
    ## Run 404 stress 0.1415299 
    ## Run 405 stress 0.1407298 
    ## Run 406 stress 0.2775523 
    ## Run 407 stress 0.2520602 
    ## Run 408 stress 0.1410449 
    ## Run 409 stress 0.2085297 
    ## Run 410 stress 0.1572668 
    ## Run 411 stress 0.1320258 
    ## ... Procrustes: rmse 7.619289e-07  max resid 1.556717e-06 
    ## ... Similar to previous best
    ## Run 412 stress 0.1407207 
    ## Run 413 stress 0.1998955 
    ## Run 414 stress 0.1830332 
    ## Run 415 stress 0.1407298 
    ## Run 416 stress 0.1970303 
    ## Run 417 stress 0.1445811 
    ## Run 418 stress 0.1693717 
    ## Run 419 stress 0.1415299 
    ## Run 420 stress 0.1410453 
    ## Run 421 stress 0.1320258 
    ## ... Procrustes: rmse 9.756279e-07  max resid 1.457619e-06 
    ## ... Similar to previous best
    ## Run 422 stress 0.1445811 
    ## Run 423 stress 0.1663544 
    ## Run 424 stress 0.1445811 
    ## Run 425 stress 0.1383681 
    ## Run 426 stress 0.1407298 
    ## Run 427 stress 0.1434197 
    ## Run 428 stress 0.1551863 
    ## Run 429 stress 0.1830332 
    ## Run 430 stress 0.1572668 
    ## Run 431 stress 0.1407207 
    ## Run 432 stress 0.283815 
    ## Run 433 stress 0.1320258 
    ## ... Procrustes: rmse 8.201225e-07  max resid 1.704266e-06 
    ## ... Similar to previous best
    ## Run 434 stress 0.1410448 
    ## Run 435 stress 0.1410451 
    ## Run 436 stress 0.2538018 
    ## Run 437 stress 0.1407298 
    ## Run 438 stress 0.2520602 
    ## Run 439 stress 0.1383682 
    ## Run 440 stress 0.1383682 
    ## Run 441 stress 0.1771966 
    ## Run 442 stress 0.1551863 
    ## Run 443 stress 0.2538018 
    ## Run 444 stress 0.1320258 
    ## ... Procrustes: rmse 4.376859e-07  max resid 7.146995e-07 
    ## ... Similar to previous best
    ## Run 445 stress 0.2749582 
    ## Run 446 stress 0.1383681 
    ## Run 447 stress 0.1383681 
    ## Run 448 stress 0.2008164 
    ## Run 449 stress 0.1415299 
    ## Run 450 stress 0.1572668 
    ## Run 451 stress 0.1407298 
    ## Run 452 stress 0.1410449 
    ## Run 453 stress 0.2008159 
    ## Run 454 stress 0.1771969 
    ## Run 455 stress 0.1830332 
    ## Run 456 stress 0.1407207 
    ## Run 457 stress 0.2842805 
    ## Run 458 stress 0.2775534 
    ## Run 459 stress 0.1410448 
    ## Run 460 stress 0.1693717 
    ## Run 461 stress 0.1445811 
    ## Run 462 stress 0.1383681 
    ## Run 463 stress 0.1572668 
    ## Run 464 stress 0.2600627 
    ## Run 465 stress 0.1383682 
    ## Run 466 stress 0.1771969 
    ## Run 467 stress 0.1383681 
    ## Run 468 stress 0.1320258 
    ## ... Procrustes: rmse 1.306828e-06  max resid 2.64606e-06 
    ## ... Similar to previous best
    ## Run 469 stress 0.1383682 
    ## Run 470 stress 0.1320258 
    ## ... Procrustes: rmse 1.170192e-06  max resid 2.626537e-06 
    ## ... Similar to previous best
    ## Run 471 stress 0.1771966 
    ## Run 472 stress 0.1407207 
    ## Run 473 stress 0.2514078 
    ## Run 474 stress 0.1572668 
    ## Run 475 stress 0.1410448 
    ## Run 476 stress 0.1320258 
    ## ... Procrustes: rmse 2.081027e-06  max resid 4.322231e-06 
    ## ... Similar to previous best
    ## Run 477 stress 0.141045 
    ## Run 478 stress 0.257865 
    ## Run 479 stress 0.1320258 
    ## ... Procrustes: rmse 2.019753e-06  max resid 4.082944e-06 
    ## ... Similar to previous best
    ## Run 480 stress 0.1320258 
    ## ... Procrustes: rmse 2.820416e-07  max resid 6.193136e-07 
    ## ... Similar to previous best
    ## Run 481 stress 0.1445811 
    ## Run 482 stress 0.1407298 
    ## Run 483 stress 0.1415299 
    ## Run 484 stress 0.2008162 
    ## Run 485 stress 0.1572668 
    ## Run 486 stress 0.1320258 
    ## ... Procrustes: rmse 6.341315e-07  max resid 9.768864e-07 
    ## ... Similar to previous best
    ## Run 487 stress 0.1410449 
    ## Run 488 stress 0.1415299 
    ## Run 489 stress 0.1572668 
    ## Run 490 stress 0.1693717 
    ## Run 491 stress 0.1970303 
    ## Run 492 stress 0.1771969 
    ## Run 493 stress 0.1551863 
    ## Run 494 stress 0.1320258 
    ## ... Procrustes: rmse 7.572884e-07  max resid 1.363516e-06 
    ## ... Similar to previous best
    ## Run 495 stress 0.1320258 
    ## ... Procrustes: rmse 1.859439e-06  max resid 3.811113e-06 
    ## ... Similar to previous best
    ## Run 496 stress 0.1663544 
    ## Run 497 stress 0.1830332 
    ## Run 498 stress 0.1572668 
    ## Run 499 stress 0.2384387 
    ## Run 500 stress 0.1383682 
    ## *** Best solution repeated 20 times

``` r
### Environmental
# Surveyed sites 
SD_beta_env_NMDS <- metaMDS(SD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07669178 
    ## Run 1 stress 0.08968589 
    ## Run 2 stress 0.08514304 
    ## Run 3 stress 0.08337152 
    ## Run 4 stress 0.08265485 
    ## Run 5 stress 0.07871737 
    ## Run 6 stress 0.08233217 
    ## Run 7 stress 0.08337201 
    ## Run 8 stress 0.07960692 
    ## Run 9 stress 0.08459584 
    ## Run 10 stress 0.08144897 
    ## Run 11 stress 0.08250125 
    ## Run 12 stress 0.08097747 
    ## Run 13 stress 0.07868344 
    ## Run 14 stress 0.0810493 
    ## Run 15 stress 0.0818212 
    ## Run 16 stress 0.08182131 
    ## Run 17 stress 0.08218525 
    ## Run 18 stress 0.08097744 
    ## Run 19 stress 0.07943342 
    ## Run 20 stress 0.08329748 
    ## Run 21 stress 0.08265437 
    ## Run 22 stress 0.08265425 
    ## Run 23 stress 0.07669156 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004498418  max resid 0.000780648 
    ## ... Similar to previous best
    ## Run 24 stress 0.0766917 
    ## ... Procrustes: rmse 0.0001425953  max resid 0.0002480094 
    ## ... Similar to previous best
    ## Run 25 stress 0.08175692 
    ## Run 26 stress 0.08485154 
    ## Run 27 stress 0.08335426 
    ## Run 28 stress 0.08144909 
    ## Run 29 stress 0.08175693 
    ## Run 30 stress 0.08189377 
    ## Run 31 stress 0.08012945 
    ## Run 32 stress 0.0796066 
    ## Run 33 stress 0.08250183 
    ## Run 34 stress 0.08568909 
    ## Run 35 stress 0.08104893 
    ## Run 36 stress 0.08003192 
    ## Run 37 stress 0.08265474 
    ## Run 38 stress 0.07960684 
    ## Run 39 stress 0.08898789 
    ## Run 40 stress 0.08305515 
    ## Run 41 stress 0.08097745 
    ## Run 42 stress 0.08360266 
    ## Run 43 stress 0.08421026 
    ## Run 44 stress 0.08898776 
    ## Run 45 stress 0.08314396 
    ## Run 46 stress 0.08175703 
    ## Run 47 stress 0.08527017 
    ## Run 48 stress 0.0817571 
    ## Run 49 stress 0.08189401 
    ## Run 50 stress 0.08364077 
    ## Run 51 stress 0.08175694 
    ## Run 52 stress 0.08097757 
    ## Run 53 stress 0.08362753 
    ## Run 54 stress 0.08252124 
    ## Run 55 stress 0.08485135 
    ## Run 56 stress 0.08898771 
    ## Run 57 stress 0.0831441 
    ## Run 58 stress 0.08243366 
    ## Run 59 stress 0.08175689 
    ## Run 60 stress 0.08553082 
    ## Run 61 stress 0.08682653 
    ## Run 62 stress 0.08241503 
    ## Run 63 stress 0.08175687 
    ## Run 64 stress 0.08527033 
    ## Run 65 stress 0.081757 
    ## Run 66 stress 0.08243397 
    ## Run 67 stress 0.07868399 
    ## Run 68 stress 0.08252792 
    ## Run 69 stress 0.07868345 
    ## Run 70 stress 0.08929064 
    ## Run 71 stress 0.08265394 
    ## Run 72 stress 0.08218525 
    ## Run 73 stress 0.08189386 
    ## Run 74 stress 0.08097738 
    ## Run 75 stress 0.07669162 
    ## ... Procrustes: rmse 0.0003424918  max resid 0.0005882023 
    ## ... Similar to previous best
    ## Run 76 stress 0.08189393 
    ## Run 77 stress 0.07868339 
    ## Run 78 stress 0.08175702 
    ## Run 79 stress 0.08682734 
    ## Run 80 stress 0.08175696 
    ## Run 81 stress 0.08218513 
    ## Run 82 stress 0.08337179 
    ## Run 83 stress 0.08252122 
    ## Run 84 stress 0.08104883 
    ## Run 85 stress 0.08283579 
    ## Run 86 stress 0.08013485 
    ## Run 87 stress 0.08012971 
    ## Run 88 stress 0.08305459 
    ## Run 89 stress 0.0824339 
    ## Run 90 stress 0.08175687 
    ## Run 91 stress 0.07731533 
    ## Run 92 stress 0.08104901 
    ## Run 93 stress 0.08182126 
    ## Run 94 stress 0.08317333 
    ## Run 95 stress 0.08217348 
    ## Run 96 stress 0.07960657 
    ## Run 97 stress 0.08477451 
    ## Run 98 stress 0.08104883 
    ## Run 99 stress 0.07948983 
    ## Run 100 stress 0.08012975 
    ## Run 101 stress 0.08175697 
    ## Run 102 stress 0.08527016 
    ## Run 103 stress 0.08252135 
    ## Run 104 stress 0.08217352 
    ## Run 105 stress 0.08182124 
    ## Run 106 stress 0.08097755 
    ## Run 107 stress 0.0787177 
    ## Run 108 stress 0.08175697 
    ## Run 109 stress 0.0809776 
    ## Run 110 stress 0.0830546 
    ## Run 111 stress 0.08175687 
    ## Run 112 stress 0.08684113 
    ## Run 113 stress 0.0826536 
    ## Run 114 stress 0.08003198 
    ## Run 115 stress 0.08104896 
    ## Run 116 stress 0.08486961 
    ## Run 117 stress 0.08305313 
    ## Run 118 stress 0.08317308 
    ## Run 119 stress 0.07871727 
    ## Run 120 stress 0.08375706 
    ## Run 121 stress 0.08379861 
    ## Run 122 stress 0.08957565 
    ## Run 123 stress 0.08527047 
    ## Run 124 stress 0.08433694 
    ## Run 125 stress 0.08233267 
    ## Run 126 stress 0.08421039 
    ## Run 127 stress 0.0766916 
    ## ... Procrustes: rmse 0.0002138673  max resid 0.0003937657 
    ## ... Similar to previous best
    ## Run 128 stress 0.08433686 
    ## Run 129 stress 0.08175711 
    ## Run 130 stress 0.08477427 
    ## Run 131 stress 0.08241508 
    ## Run 132 stress 0.07871741 
    ## Run 133 stress 0.08003207 
    ## Run 134 stress 0.08003219 
    ## Run 135 stress 0.08305462 
    ## Run 136 stress 0.08433698 
    ## Run 137 stress 0.08252118 
    ## Run 138 stress 0.08314402 
    ## Run 139 stress 0.08317369 
    ## Run 140 stress 0.07951473 
    ## Run 141 stress 0.08337154 
    ## Run 142 stress 0.0766916 
    ## ... Procrustes: rmse 0.0003281795  max resid 0.000590764 
    ## ... Similar to previous best
    ## Run 143 stress 0.08175705 
    ## Run 144 stress 0.08314374 
    ## Run 145 stress 0.07868342 
    ## Run 146 stress 0.08104888 
    ## Run 147 stress 0.08182125 
    ## Run 148 stress 0.07948968 
    ## Run 149 stress 0.08252119 
    ## Run 150 stress 0.08104911 
    ## Run 151 stress 0.08375721 
    ## Run 152 stress 0.08360268 
    ## Run 153 stress 0.08317343 
    ## Run 154 stress 0.08252764 
    ## Run 155 stress 0.08375748 
    ## Run 156 stress 0.08265396 
    ## Run 157 stress 0.08012956 
    ## Run 158 stress 0.08355807 
    ## Run 159 stress 0.08957505 
    ## Run 160 stress 0.08104889 
    ## Run 161 stress 0.08360257 
    ## Run 162 stress 0.084337 
    ## Run 163 stress 0.08317359 
    ## Run 164 stress 0.08104919 
    ## Run 165 stress 0.07871747 
    ## Run 166 stress 0.08421015 
    ## Run 167 stress 0.08218529 
    ## Run 168 stress 0.08104906 
    ## Run 169 stress 0.08104908 
    ## Run 170 stress 0.0810491 
    ## Run 171 stress 0.08957496 
    ## Run 172 stress 0.08500374 
    ## Run 173 stress 0.08175703 
    ## Run 174 stress 0.08317313 
    ## Run 175 stress 0.081049 
    ## Run 176 stress 0.08548311 
    ## Run 177 stress 0.08097739 
    ## Run 178 stress 0.07960658 
    ## Run 179 stress 0.07871745 
    ## Run 180 stress 0.07871727 
    ## Run 181 stress 0.08493076 
    ## Run 182 stress 0.07960687 
    ## Run 183 stress 0.08012976 
    ## Run 184 stress 0.08175713 
    ## Run 185 stress 0.08485134 
    ## Run 186 stress 0.08175689 
    ## Run 187 stress 0.07943349 
    ## Run 188 stress 0.08243392 
    ## Run 189 stress 0.08684164 
    ## Run 190 stress 0.08305449 
    ## Run 191 stress 0.08175691 
    ## Run 192 stress 0.07871748 
    ## Run 193 stress 0.08335398 
    ## Run 194 stress 0.08241501 
    ## Run 195 stress 0.08314418 
    ## Run 196 stress 0.08684348 
    ## Run 197 stress 0.07669152 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001034226  max resid 0.0001760882 
    ## ... Similar to previous best
    ## Run 198 stress 0.08470964 
    ## Run 199 stress 0.08265348 
    ## Run 200 stress 0.07960675 
    ## Run 201 stress 0.08337174 
    ## Run 202 stress 0.0809776 
    ## Run 203 stress 0.0801297 
    ## Run 204 stress 0.08274907 
    ## Run 205 stress 0.08527039 
    ## Run 206 stress 0.08097744 
    ## Run 207 stress 0.08553085 
    ## Run 208 stress 0.07868374 
    ## Run 209 stress 0.08927512 
    ## Run 210 stress 0.08252116 
    ## Run 211 stress 0.08175692 
    ## Run 212 stress 0.07871751 
    ## Run 213 stress 0.07951495 
    ## Run 214 stress 0.08530162 
    ## Run 215 stress 0.07731529 
    ## Run 216 stress 0.07868324 
    ## Run 217 stress 0.08242197 
    ## Run 218 stress 0.08243369 
    ## Run 219 stress 0.08252125 
    ## Run 220 stress 0.08144901 
    ## Run 221 stress 0.08175691 
    ## Run 222 stress 0.08477436 
    ## Run 223 stress 0.08265349 
    ## Run 224 stress 0.08012938 
    ## Run 225 stress 0.08243402 
    ## Run 226 stress 0.08684158 
    ## Run 227 stress 0.08013555 
    ## Run 228 stress 0.08527036 
    ## Run 229 stress 0.0826548 
    ## Run 230 stress 0.0801296 
    ## Run 231 stress 0.08003212 
    ## Run 232 stress 0.08175702 
    ## Run 233 stress 0.08175697 
    ## Run 234 stress 0.08370287 
    ## Run 235 stress 0.08175721 
    ## Run 236 stress 0.07669154 
    ## ... Procrustes: rmse 0.0001770932  max resid 0.0003687926 
    ## ... Similar to previous best
    ## Run 237 stress 0.08527021 
    ## Run 238 stress 0.08527031 
    ## Run 239 stress 0.08097775 
    ## Run 240 stress 0.08003213 
    ## Run 241 stress 0.08428709 
    ## Run 242 stress 0.07949025 
    ## Run 243 stress 0.08243375 
    ## Run 244 stress 0.08104884 
    ## Run 245 stress 0.08012963 
    ## Run 246 stress 0.07948985 
    ## Run 247 stress 0.08250152 
    ## Run 248 stress 0.0766918 
    ## ... Procrustes: rmse 0.0003737292  max resid 0.0007094807 
    ## ... Similar to previous best
    ## Run 249 stress 0.07669153 
    ## ... Procrustes: rmse 4.540265e-05  max resid 9.336804e-05 
    ## ... Similar to previous best
    ## Run 250 stress 0.0817569 
    ## Run 251 stress 0.08175711 
    ## Run 252 stress 0.08530177 
    ## Run 253 stress 0.08499324 
    ## Run 254 stress 0.08233206 
    ## Run 255 stress 0.08364069 
    ## Run 256 stress 0.07949013 
    ## Run 257 stress 0.07871736 
    ## Run 258 stress 0.08568765 
    ## Run 259 stress 0.08470918 
    ## Run 260 stress 0.08175706 
    ## Run 261 stress 0.08265376 
    ## Run 262 stress 0.08305283 
    ## Run 263 stress 0.07669153 
    ## ... Procrustes: rmse 0.0001517515  max resid 0.0003075657 
    ## ... Similar to previous best
    ## Run 264 stress 0.08241483 
    ## Run 265 stress 0.08876527 
    ## Run 266 stress 0.08252122 
    ## Run 267 stress 0.0853015 
    ## Run 268 stress 0.08453233 
    ## Run 269 stress 0.08433684 
    ## Run 270 stress 0.08249633 
    ## Run 271 stress 0.07669176 
    ## ... Procrustes: rmse 0.0002520828  max resid 0.0004282769 
    ## ... Similar to previous best
    ## Run 272 stress 0.08274946 
    ## Run 273 stress 0.07948985 
    ## Run 274 stress 0.08175712 
    ## Run 275 stress 0.08218517 
    ## Run 276 stress 0.08443288 
    ## Run 277 stress 0.08175709 
    ## Run 278 stress 0.08175688 
    ## Run 279 stress 0.07669164 
    ## ... Procrustes: rmse 0.0002767732  max resid 0.0005594514 
    ## ... Similar to previous best
    ## Run 280 stress 0.0773154 
    ## Run 281 stress 0.07871769 
    ## Run 282 stress 0.08097748 
    ## Run 283 stress 0.08104927 
    ## Run 284 stress 0.08175699 
    ## Run 285 stress 0.08433688 
    ## Run 286 stress 0.08957515 
    ## Run 287 stress 0.08360266 
    ## Run 288 stress 0.0818211 
    ## Run 289 stress 0.08265387 
    ## Run 290 stress 0.08003209 
    ## Run 291 stress 0.07871733 
    ## Run 292 stress 0.08453188 
    ## Run 293 stress 0.08012952 
    ## Run 294 stress 0.08233216 
    ## Run 295 stress 0.08553071 
    ## Run 296 stress 0.08104906 
    ## Run 297 stress 0.07871742 
    ## Run 298 stress 0.08667542 
    ## Run 299 stress 0.07731513 
    ## Run 300 stress 0.08421043 
    ## Run 301 stress 0.07871726 
    ## Run 302 stress 0.08097758 
    ## Run 303 stress 0.08097774 
    ## Run 304 stress 0.08144956 
    ## Run 305 stress 0.08530164 
    ## Run 306 stress 0.08012978 
    ## Run 307 stress 0.07943371 
    ## Run 308 stress 0.08175697 
    ## Run 309 stress 0.07868353 
    ## Run 310 stress 0.08553061 
    ## Run 311 stress 0.0931125 
    ## Run 312 stress 0.08405601 
    ## Run 313 stress 0.07868364 
    ## Run 314 stress 0.08335412 
    ## Run 315 stress 0.07669171 
    ## ... Procrustes: rmse 0.0003190304  max resid 0.0006121426 
    ## ... Similar to previous best
    ## Run 316 stress 0.08443358 
    ## Run 317 stress 0.08003188 
    ## Run 318 stress 0.08527048 
    ## Run 319 stress 0.08950139 
    ## Run 320 stress 0.08387139 
    ## Run 321 stress 0.08265462 
    ## Run 322 stress 0.08314382 
    ## Run 323 stress 0.08189375 
    ## Run 324 stress 0.08364081 
    ## Run 325 stress 0.08104912 
    ## Run 326 stress 0.08252127 
    ## Run 327 stress 0.08252129 
    ## Run 328 stress 0.08362768 
    ## Run 329 stress 0.08249597 
    ## Run 330 stress 0.08370303 
    ## Run 331 stress 0.08418573 
    ## Run 332 stress 0.08443251 
    ## Run 333 stress 0.0817571 
    ## Run 334 stress 0.08684111 
    ## Run 335 stress 0.08317353 
    ## Run 336 stress 0.08459594 
    ## Run 337 stress 0.07871727 
    ## Run 338 stress 0.08097745 
    ## Run 339 stress 0.08252132 
    ## Run 340 stress 0.07669155 
    ## ... Procrustes: rmse 0.0001448998  max resid 0.0002676739 
    ## ... Similar to previous best
    ## Run 341 stress 0.07669174 
    ## ... Procrustes: rmse 0.0003302653  max resid 0.0006489663 
    ## ... Similar to previous best
    ## Run 342 stress 0.07951512 
    ## Run 343 stress 0.07948974 
    ## Run 344 stress 0.08097755 
    ## Run 345 stress 0.08097761 
    ## Run 346 stress 0.08250138 
    ## Run 347 stress 0.0828355 
    ## Run 348 stress 0.08217354 
    ## Run 349 stress 0.08175701 
    ## Run 350 stress 0.0855307 
    ## Run 351 stress 0.08274935 
    ## Run 352 stress 0.08265375 
    ## Run 353 stress 0.08175696 
    ## Run 354 stress 0.07669173 
    ## ... Procrustes: rmse 0.0002154775  max resid 0.0003929843 
    ## ... Similar to previous best
    ## Run 355 stress 0.08835165 
    ## Run 356 stress 0.08314408 
    ## Run 357 stress 0.08487596 
    ## Run 358 stress 0.08422689 
    ## Run 359 stress 0.08104911 
    ## Run 360 stress 0.08433687 
    ## Run 361 stress 0.09044883 
    ## Run 362 stress 0.08097755 
    ## Run 363 stress 0.07948976 
    ## Run 364 stress 0.08835136 
    ## Run 365 stress 0.08433685 
    ## Run 366 stress 0.08175693 
    ## Run 367 stress 0.0883514 
    ## Run 368 stress 0.08104883 
    ## Run 369 stress 0.08241484 
    ## Run 370 stress 0.08314385 
    ## Run 371 stress 0.09063312 
    ## Run 372 stress 0.08217349 
    ## Run 373 stress 0.08329758 
    ## Run 374 stress 0.09044898 
    ## Run 375 stress 0.08097766 
    ## Run 376 stress 0.08175688 
    ## Run 377 stress 0.08249622 
    ## Run 378 stress 0.08527026 
    ## Run 379 stress 0.07871755 
    ## Run 380 stress 0.07949008 
    ## Run 381 stress 0.0796066 
    ## Run 382 stress 0.08097746 
    ## Run 383 stress 0.08459573 
    ## Run 384 stress 0.0809777 
    ## Run 385 stress 0.07951506 
    ## Run 386 stress 0.08500384 
    ## Run 387 stress 0.07731519 
    ## Run 388 stress 0.08375732 
    ## Run 389 stress 0.08265423 
    ## Run 390 stress 0.08584214 
    ## Run 391 stress 0.07943324 
    ## Run 392 stress 0.08530151 
    ## Run 393 stress 0.0817569 
    ## Run 394 stress 0.08433685 
    ## Run 395 stress 0.0827492 
    ## Run 396 stress 0.07948962 
    ## Run 397 stress 0.08317346 
    ## Run 398 stress 0.08379946 
    ## Run 399 stress 0.08182115 
    ## Run 400 stress 0.08486939 
    ## Run 401 stress 0.08514241 
    ## Run 402 stress 0.08684136 
    ## Run 403 stress 0.08842086 
    ## Run 404 stress 0.08003201 
    ## Run 405 stress 0.08252145 
    ## Run 406 stress 0.08252122 
    ## Run 407 stress 0.08252149 
    ## Run 408 stress 0.07871754 
    ## Run 409 stress 0.08182114 
    ## Run 410 stress 0.07731517 
    ## Run 411 stress 0.07669158 
    ## ... Procrustes: rmse 0.0001546118  max resid 0.000270998 
    ## ... Similar to previous best
    ## Run 412 stress 0.07868345 
    ## Run 413 stress 0.08615932 
    ## Run 414 stress 0.08471581 
    ## Run 415 stress 0.07943328 
    ## Run 416 stress 0.07871747 
    ## Run 417 stress 0.08217348 
    ## Run 418 stress 0.08684083 
    ## Run 419 stress 0.08433693 
    ## Run 420 stress 0.08335398 
    ## Run 421 stress 0.08249612 
    ## Run 422 stress 0.08314377 
    ## Run 423 stress 0.08217342 
    ## Run 424 stress 0.07868369 
    ## Run 425 stress 0.08189407 
    ## Run 426 stress 0.08968621 
    ## Run 427 stress 0.08097761 
    ## Run 428 stress 0.08314381 
    ## Run 429 stress 0.08252117 
    ## Run 430 stress 0.08097746 
    ## Run 431 stress 0.08835166 
    ## Run 432 stress 0.08682716 
    ## Run 433 stress 0.08252133 
    ## Run 434 stress 0.08104887 
    ## Run 435 stress 0.08003214 
    ## Run 436 stress 0.08314376 
    ## Run 437 stress 0.0818938 
    ## Run 438 stress 0.08527043 
    ## Run 439 stress 0.07868395 
    ## Run 440 stress 0.08684131 
    ## Run 441 stress 0.08175706 
    ## Run 442 stress 0.08104886 
    ## Run 443 stress 0.08175695 
    ## Run 444 stress 0.08012948 
    ## Run 445 stress 0.07731513 
    ## Run 446 stress 0.08553072 
    ## Run 447 stress 0.08487617 
    ## Run 448 stress 0.07948996 
    ## Run 449 stress 0.08375738 
    ## Run 450 stress 0.0766916 
    ## ... Procrustes: rmse 0.0002077588  max resid 0.0003669403 
    ## ... Similar to previous best
    ## Run 451 stress 0.07669152 
    ## ... New best solution
    ## ... Procrustes: rmse 2.792338e-05  max resid 4.811125e-05 
    ## ... Similar to previous best
    ## Run 452 stress 0.08249639 
    ## Run 453 stress 0.08283558 
    ## Run 454 stress 0.08957609 
    ## Run 455 stress 0.08337168 
    ## Run 456 stress 0.08433714 
    ## Run 457 stress 0.08305329 
    ## Run 458 stress 0.08835143 
    ## Run 459 stress 0.0830547 
    ## Run 460 stress 0.08527023 
    ## Run 461 stress 0.08252141 
    ## Run 462 stress 0.08241482 
    ## Run 463 stress 0.08250126 
    ## Run 464 stress 0.0817569 
    ## Run 465 stress 0.07868337 
    ## Run 466 stress 0.08175697 
    ## Run 467 stress 0.08252151 
    ## Run 468 stress 0.08104928 
    ## Run 469 stress 0.07871734 
    ## Run 470 stress 0.08305361 
    ## Run 471 stress 0.0801297 
    ## Run 472 stress 0.07669163 
    ## ... Procrustes: rmse 0.0002106409  max resid 0.0003586055 
    ## ... Similar to previous best
    ## Run 473 stress 0.08265374 
    ## Run 474 stress 0.08471584 
    ## Run 475 stress 0.08003199 
    ## Run 476 stress 0.08265446 
    ## Run 477 stress 0.08360282 
    ## Run 478 stress 0.08485144 
    ## Run 479 stress 0.08217355 
    ## Run 480 stress 0.08360267 
    ## Run 481 stress 0.08405579 
    ## Run 482 stress 0.08379864 
    ## Run 483 stress 0.08175699 
    ## Run 484 stress 0.08265462 
    ## Run 485 stress 0.08314433 
    ## Run 486 stress 0.08012988 
    ## Run 487 stress 0.08243378 
    ## Run 488 stress 0.0817569 
    ## Run 489 stress 0.08337196 
    ## Run 490 stress 0.08929125 
    ## Run 491 stress 0.08443211 
    ## Run 492 stress 0.0826546 
    ## Run 493 stress 0.0806395 
    ## Run 494 stress 0.07669181 
    ## ... Procrustes: rmse 0.0003415229  max resid 0.0005769111 
    ## ... Similar to previous best
    ## Run 495 stress 0.07871742 
    ## Run 496 stress 0.08927451 
    ## Run 497 stress 0.07669156 
    ## ... Procrustes: rmse 0.0001530509  max resid 0.0002737731 
    ## ... Similar to previous best
    ## Run 498 stress 0.08274893 
    ## Run 499 stress 0.08684137 
    ## Run 500 stress 0.08337144 
    ## *** Best solution repeated 4 times

``` r
# Mixed and stratified lakes
SD_beta_env_MS_NMDS <- metaMDS(SD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.09407979 
    ## Run 2 stress 0.2913191 
    ## Run 3 stress 0.09381847 
    ## Run 4 stress 0.09168955 
    ## Run 5 stress 0.1038038 
    ## Run 6 stress 0.08503604 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001143094  max resid 0.0002573003 
    ## ... Similar to previous best
    ## Run 7 stress 0.09159083 
    ## Run 8 stress 0.09030409 
    ## Run 9 stress 0.09268323 
    ## Run 10 stress 0.09168934 
    ## Run 11 stress 0.09268344 
    ## Run 12 stress 0.1072864 
    ## Run 13 stress 0.09308954 
    ## Run 14 stress 0.09308974 
    ## Run 15 stress 0.09159085 
    ## Run 16 stress 0.08773473 
    ## Run 17 stress 0.09407989 
    ## Run 18 stress 0.08773466 
    ## Run 19 stress 0.08503463 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001400013  max resid 0.003556733 
    ## ... Similar to previous best
    ## Run 20 stress 0.09168927 
    ## Run 21 stress 0.09268354 
    ## Run 22 stress 0.09416131 
    ## Run 23 stress 0.08503549 
    ## ... Procrustes: rmse 0.001148311  max resid 0.002838248 
    ## ... Similar to previous best
    ## Run 24 stress 0.09030421 
    ## Run 25 stress 0.0930896 
    ## Run 26 stress 0.09539195 
    ## Run 27 stress 0.09969986 
    ## Run 28 stress 0.09159095 
    ## Run 29 stress 0.09407978 
    ## Run 30 stress 0.08503563 
    ## ... Procrustes: rmse 0.001215816  max resid 0.00306428 
    ## ... Similar to previous best
    ## Run 31 stress 0.09407966 
    ## Run 32 stress 0.09159083 
    ## Run 33 stress 0.09539205 
    ## Run 34 stress 0.09416135 
    ## Run 35 stress 0.3193146 
    ## Run 36 stress 0.09030398 
    ## Run 37 stress 0.09381864 
    ## Run 38 stress 0.09407965 
    ## Run 39 stress 0.08973879 
    ## Run 40 stress 0.08973864 
    ## Run 41 stress 0.09030402 
    ## Run 42 stress 0.09168936 
    ## Run 43 stress 0.08503509 
    ## ... Procrustes: rmse 0.0008684028  max resid 0.002125387 
    ## ... Similar to previous best
    ## Run 44 stress 0.09465905 
    ## Run 45 stress 0.09416132 
    ## Run 46 stress 0.09321658 
    ## Run 47 stress 0.09268396 
    ## Run 48 stress 0.09407971 
    ## Run 49 stress 0.09590968 
    ## Run 50 stress 0.09030394 
    ## Run 51 stress 0.09168927 
    ## Run 52 stress 0.09590943 
    ## Run 53 stress 0.0928609 
    ## Run 54 stress 0.09381847 
    ## Run 55 stress 0.09535488 
    ## Run 56 stress 0.08773482 
    ## Run 57 stress 0.08440269 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0141998  max resid 0.04302515 
    ## Run 58 stress 0.08973878 
    ## Run 59 stress 0.09159087 
    ## Run 60 stress 0.09268368 
    ## Run 61 stress 0.0966984 
    ## Run 62 stress 0.08440252 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003174049  max resid 0.0005854016 
    ## ... Similar to previous best
    ## Run 63 stress 0.08773471 
    ## Run 64 stress 0.09407989 
    ## Run 65 stress 0.097212 
    ## Run 66 stress 0.09416131 
    ## Run 67 stress 0.09465905 
    ## Run 68 stress 0.09464389 
    ## Run 69 stress 0.08773487 
    ## Run 70 stress 0.0844026 
    ## ... Procrustes: rmse 0.0002391341  max resid 0.0004479353 
    ## ... Similar to previous best
    ## Run 71 stress 0.09760814 
    ## Run 72 stress 0.09145321 
    ## Run 73 stress 0.0850348 
    ## Run 74 stress 0.09268324 
    ## Run 75 stress 0.08440274 
    ## ... Procrustes: rmse 0.0003484752  max resid 0.0006326589 
    ## ... Similar to previous best
    ## Run 76 stress 0.09308982 
    ## Run 77 stress 0.09407967 
    ## Run 78 stress 0.3067549 
    ## Run 79 stress 0.09268339 
    ## Run 80 stress 0.09969964 
    ## Run 81 stress 0.08440252 
    ## ... Procrustes: rmse 0.0001331927  max resid 0.0002543817 
    ## ... Similar to previous best
    ## Run 82 stress 0.09464469 
    ## Run 83 stress 0.09168935 
    ## Run 84 stress 0.09407994 
    ## Run 85 stress 0.08773465 
    ## Run 86 stress 0.08773481 
    ## Run 87 stress 0.1052851 
    ## Run 88 stress 0.09465904 
    ## Run 89 stress 0.0916894 
    ## Run 90 stress 0.09465912 
    ## Run 91 stress 0.09465913 
    ## Run 92 stress 0.08773472 
    ## Run 93 stress 0.09400529 
    ## Run 94 stress 0.09337269 
    ## Run 95 stress 0.08503494 
    ## Run 96 stress 0.09416151 
    ## Run 97 stress 0.09308947 
    ## Run 98 stress 0.09445513 
    ## Run 99 stress 0.0877347 
    ## Run 100 stress 0.08503648 
    ## Run 101 stress 0.1112169 
    ## Run 102 stress 0.09374258 
    ## Run 103 stress 0.08440263 
    ## ... Procrustes: rmse 0.0001496423  max resid 0.0003131173 
    ## ... Similar to previous best
    ## Run 104 stress 0.08503523 
    ## Run 105 stress 0.09337242 
    ## Run 106 stress 0.09030404 
    ## Run 107 stress 0.08440263 
    ## ... Procrustes: rmse 0.0002545859  max resid 0.0005384412 
    ## ... Similar to previous best
    ## Run 108 stress 0.09535511 
    ## Run 109 stress 0.08973872 
    ## Run 110 stress 0.08503495 
    ## Run 111 stress 0.09969984 
    ## Run 112 stress 0.09286083 
    ## Run 113 stress 0.09337265 
    ## Run 114 stress 0.08973889 
    ## Run 115 stress 0.09465905 
    ## Run 116 stress 0.09268366 
    ## Run 117 stress 0.08440259 
    ## ... Procrustes: rmse 0.00023439  max resid 0.0004423732 
    ## ... Similar to previous best
    ## Run 118 stress 0.09407998 
    ## Run 119 stress 0.09465905 
    ## Run 120 stress 0.09464493 
    ## Run 121 stress 0.09374344 
    ## Run 122 stress 0.09030396 
    ## Run 123 stress 0.09286085 
    ## Run 124 stress 0.09407992 
    ## Run 125 stress 0.09464474 
    ## Run 126 stress 0.09030403 
    ## Run 127 stress 0.09159086 
    ## Run 128 stress 0.09417823 
    ## Run 129 stress 0.09407967 
    ## Run 130 stress 0.09969967 
    ## Run 131 stress 0.1052851 
    ## Run 132 stress 0.08773465 
    ## Run 133 stress 0.09400522 
    ## Run 134 stress 0.09268355 
    ## Run 135 stress 0.08503479 
    ## Run 136 stress 0.09407981 
    ## Run 137 stress 0.09308953 
    ## Run 138 stress 0.09159091 
    ## Run 139 stress 0.08440255 
    ## ... Procrustes: rmse 4.570672e-05  max resid 0.0001019171 
    ## ... Similar to previous best
    ## Run 140 stress 0.09030398 
    ## Run 141 stress 0.08503496 
    ## Run 142 stress 0.0946591 
    ## Run 143 stress 0.09168945 
    ## Run 144 stress 0.09030408 
    ## Run 145 stress 0.09407967 
    ## Run 146 stress 0.08503489 
    ## Run 147 stress 0.09159085 
    ## Run 148 stress 0.08503565 
    ## Run 149 stress 0.09407981 
    ## Run 150 stress 0.09416133 
    ## Run 151 stress 0.09445496 
    ## Run 152 stress 0.09145325 
    ## Run 153 stress 0.0932191 
    ## Run 154 stress 0.09408015 
    ## Run 155 stress 0.0916894 
    ## Run 156 stress 0.0940711 
    ## Run 157 stress 0.09465905 
    ## Run 158 stress 0.09308985 
    ## Run 159 stress 0.09721199 
    ## Run 160 stress 0.09030398 
    ## Run 161 stress 0.09030396 
    ## Run 162 stress 0.0877347 
    ## Run 163 stress 0.09465909 
    ## Run 164 stress 0.08503477 
    ## Run 165 stress 0.09337229 
    ## Run 166 stress 0.09539195 
    ## Run 167 stress 0.0940052 
    ## Run 168 stress 0.3052587 
    ## Run 169 stress 0.090304 
    ## Run 170 stress 0.09465906 
    ## Run 171 stress 0.09400527 
    ## Run 172 stress 0.1038037 
    ## Run 173 stress 0.08503575 
    ## Run 174 stress 0.08503587 
    ## Run 175 stress 0.09407968 
    ## Run 176 stress 0.09407135 
    ## Run 177 stress 0.09159087 
    ## Run 178 stress 0.09159092 
    ## Run 179 stress 0.0937421 
    ## Run 180 stress 0.09416136 
    ## Run 181 stress 0.08773476 
    ## Run 182 stress 0.1038044 
    ## Run 183 stress 0.08503497 
    ## Run 184 stress 0.0897387 
    ## Run 185 stress 0.09407974 
    ## Run 186 stress 0.0937433 
    ## Run 187 stress 0.09465905 
    ## Run 188 stress 0.09030396 
    ## Run 189 stress 0.09407976 
    ## Run 190 stress 0.09321663 
    ## Run 191 stress 0.09465904 
    ## Run 192 stress 0.08503511 
    ## Run 193 stress 0.09030393 
    ## Run 194 stress 0.09168957 
    ## Run 195 stress 0.09407969 
    ## Run 196 stress 0.09268374 
    ## Run 197 stress 0.09337295 
    ## Run 198 stress 0.09445498 
    ## Run 199 stress 0.2450179 
    ## Run 200 stress 0.09374284 
    ## Run 201 stress 0.08503456 
    ## Run 202 stress 0.093743 
    ## Run 203 stress 0.09159084 
    ## Run 204 stress 0.09374228 
    ## Run 205 stress 0.09030394 
    ## Run 206 stress 0.09447324 
    ## Run 207 stress 0.09417788 
    ## Run 208 stress 0.09447297 
    ## Run 209 stress 0.08973873 
    ## Run 210 stress 0.09465905 
    ## Run 211 stress 0.09969969 
    ## Run 212 stress 0.09286085 
    ## Run 213 stress 0.09145324 
    ## Run 214 stress 0.09407984 
    ## Run 215 stress 0.09407973 
    ## Run 216 stress 0.0997001 
    ## Run 217 stress 0.08440253 
    ## ... Procrustes: rmse 0.0001463081  max resid 0.0002884153 
    ## ... Similar to previous best
    ## Run 218 stress 0.09337287 
    ## Run 219 stress 0.09030402 
    ## Run 220 stress 0.1053421 
    ## Run 221 stress 0.1004611 
    ## Run 222 stress 0.08503644 
    ## Run 223 stress 0.08773466 
    ## Run 224 stress 0.08973894 
    ## Run 225 stress 0.09407987 
    ## Run 226 stress 0.08503455 
    ## Run 227 stress 0.09407439 
    ## Run 228 stress 0.08440261 
    ## ... Procrustes: rmse 0.0002291911  max resid 0.0004926299 
    ## ... Similar to previous best
    ## Run 229 stress 0.09407137 
    ## Run 230 stress 0.09407966 
    ## Run 231 stress 0.09337301 
    ## Run 232 stress 0.09337281 
    ## Run 233 stress 0.1026216 
    ## Run 234 stress 0.09337276 
    ## Run 235 stress 0.09159089 
    ## Run 236 stress 0.09168937 
    ## Run 237 stress 0.09337273 
    ## Run 238 stress 0.09159098 
    ## Run 239 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001890033  max resid 0.000376708 
    ## ... Similar to previous best
    ## Run 240 stress 0.08440262 
    ## ... Procrustes: rmse 0.0001423865  max resid 0.000282464 
    ## ... Similar to previous best
    ## Run 241 stress 0.08440267 
    ## ... Procrustes: rmse 0.0003041394  max resid 0.0005702918 
    ## ... Similar to previous best
    ## Run 242 stress 0.09268335 
    ## Run 243 stress 0.2466508 
    ## Run 244 stress 0.09464472 
    ## Run 245 stress 0.09380575 
    ## Run 246 stress 0.09465904 
    ## Run 247 stress 0.09445465 
    ## Run 248 stress 0.1072865 
    ## Run 249 stress 0.09268351 
    ## Run 250 stress 0.09030403 
    ## Run 251 stress 0.09268405 
    ## Run 252 stress 0.08973891 
    ## Run 253 stress 0.09407965 
    ## Run 254 stress 0.09374184 
    ## Run 255 stress 0.08503508 
    ## Run 256 stress 0.08973868 
    ## Run 257 stress 0.09308955 
    ## Run 258 stress 0.09760833 
    ## Run 259 stress 0.2465515 
    ## Run 260 stress 0.09286094 
    ## Run 261 stress 0.0933728 
    ## Run 262 stress 0.09337241 
    ## Run 263 stress 0.09145325 
    ## Run 264 stress 0.09380587 
    ## Run 265 stress 0.1072862 
    ## Run 266 stress 0.09030393 
    ## Run 267 stress 0.08773466 
    ## Run 268 stress 0.09286083 
    ## Run 269 stress 0.1026215 
    ## Run 270 stress 0.09445498 
    ## Run 271 stress 0.08973874 
    ## Run 272 stress 0.09030404 
    ## Run 273 stress 0.1038042 
    ## Run 274 stress 0.09030396 
    ## Run 275 stress 0.09145329 
    ## Run 276 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001270251  max resid 0.000241816 
    ## ... Similar to previous best
    ## Run 277 stress 0.09407512 
    ## Run 278 stress 0.0877347 
    ## Run 279 stress 0.09286084 
    ## Run 280 stress 0.08973872 
    ## Run 281 stress 0.0850359 
    ## Run 282 stress 0.08440256 
    ## ... Procrustes: rmse 7.53328e-05  max resid 0.0001591818 
    ## ... Similar to previous best
    ## Run 283 stress 0.09760815 
    ## Run 284 stress 0.09030407 
    ## Run 285 stress 0.09416132 
    ## Run 286 stress 0.09145332 
    ## Run 287 stress 0.09145322 
    ## Run 288 stress 0.08503463 
    ## Run 289 stress 0.09286084 
    ## Run 290 stress 0.09030393 
    ## Run 291 stress 0.09030402 
    ## Run 292 stress 0.09145326 
    ## Run 293 stress 0.09539195 
    ## Run 294 stress 0.09465907 
    ## Run 295 stress 0.08773466 
    ## Run 296 stress 0.0938186 
    ## Run 297 stress 0.08503463 
    ## Run 298 stress 0.0996998 
    ## Run 299 stress 0.09145326 
    ## Run 300 stress 0.1026217 
    ## Run 301 stress 0.09168936 
    ## Run 302 stress 0.09590894 
    ## Run 303 stress 0.09465904 
    ## Run 304 stress 0.09407967 
    ## Run 305 stress 0.09977533 
    ## Run 306 stress 0.09669839 
    ## Run 307 stress 0.09268324 
    ## Run 308 stress 0.08503551 
    ## Run 309 stress 0.09321864 
    ## Run 310 stress 0.09416143 
    ## Run 311 stress 0.08503605 
    ## Run 312 stress 0.08973869 
    ## Run 313 stress 0.0850347 
    ## Run 314 stress 0.09159084 
    ## Run 315 stress 0.08973867 
    ## Run 316 stress 0.09145334 
    ## Run 317 stress 0.09168955 
    ## Run 318 stress 0.09407993 
    ## Run 319 stress 0.09407997 
    ## Run 320 stress 0.09407979 
    ## Run 321 stress 0.09337281 
    ## Run 322 stress 0.08503483 
    ## Run 323 stress 0.09159084 
    ## Run 324 stress 0.08973865 
    ## Run 325 stress 0.0938185 
    ## Run 326 stress 0.08503605 
    ## Run 327 stress 0.08773466 
    ## Run 328 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 4.694117e-06  max resid 8.262041e-06 
    ## ... Similar to previous best
    ## Run 329 stress 0.09407123 
    ## Run 330 stress 0.09465906 
    ## Run 331 stress 0.2471333 
    ## Run 332 stress 0.08773466 
    ## Run 333 stress 0.09286093 
    ## Run 334 stress 0.09030397 
    ## Run 335 stress 0.09030408 
    ## Run 336 stress 0.09465906 
    ## Run 337 stress 0.09465905 
    ## Run 338 stress 0.09380576 
    ## Run 339 stress 0.09308979 
    ## Run 340 stress 0.09380576 
    ## Run 341 stress 0.08973893 
    ## Run 342 stress 0.09286084 
    ## Run 343 stress 0.103805 
    ## Run 344 stress 0.08503504 
    ## Run 345 stress 0.09168935 
    ## Run 346 stress 0.08440259 
    ## ... Procrustes: rmse 0.0002170952  max resid 0.0004550848 
    ## ... Similar to previous best
    ## Run 347 stress 0.08503536 
    ## Run 348 stress 0.09286084 
    ## Run 349 stress 0.08973862 
    ## Run 350 stress 0.09416136 
    ## Run 351 stress 0.1038041 
    ## Run 352 stress 0.08503486 
    ## Run 353 stress 0.09969976 
    ## Run 354 stress 0.09030399 
    ## Run 355 stress 0.09400515 
    ## Run 356 stress 0.08973862 
    ## Run 357 stress 0.08973864 
    ## Run 358 stress 0.09374344 
    ## Run 359 stress 0.2510662 
    ## Run 360 stress 0.09337257 
    ## Run 361 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001903864  max resid 0.0003558251 
    ## ... Similar to previous best
    ## Run 362 stress 0.09268438 
    ## Run 363 stress 0.09168941 
    ## Run 364 stress 0.09416136 
    ## Run 365 stress 0.09337274 
    ## Run 366 stress 0.08773481 
    ## Run 367 stress 0.09030393 
    ## Run 368 stress 0.08503561 
    ## Run 369 stress 0.08773483 
    ## Run 370 stress 0.09030397 
    ## Run 371 stress 0.09535423 
    ## Run 372 stress 0.08503476 
    ## Run 373 stress 0.08973883 
    ## Run 374 stress 0.09381844 
    ## Run 375 stress 0.09465907 
    ## Run 376 stress 0.08503669 
    ## Run 377 stress 0.09308982 
    ## Run 378 stress 0.08773466 
    ## Run 379 stress 0.09030395 
    ## Run 380 stress 0.08440261 
    ## ... Procrustes: rmse 0.0001183975  max resid 0.000336496 
    ## ... Similar to previous best
    ## Run 381 stress 0.09145324 
    ## Run 382 stress 0.09407973 
    ## Run 383 stress 0.09030395 
    ## Run 384 stress 0.09374317 
    ## Run 385 stress 0.09374298 
    ## Run 386 stress 0.08773478 
    ## Run 387 stress 0.09465904 
    ## Run 388 stress 0.09969968 
    ## Run 389 stress 0.09308946 
    ## Run 390 stress 0.08973862 
    ## Run 391 stress 0.08973877 
    ## Run 392 stress 0.09535517 
    ## Run 393 stress 0.09590506 
    ## Run 394 stress 0.09381844 
    ## Run 395 stress 0.09145322 
    ## Run 396 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 1.747327e-05  max resid 3.797089e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.0930898 
    ## Run 398 stress 0.08440268 
    ## ... Procrustes: rmse 0.0002038008  max resid 0.0004034167 
    ## ... Similar to previous best
    ## Run 399 stress 0.09030398 
    ## Run 400 stress 0.08503499 
    ## Run 401 stress 0.2772439 
    ## Run 402 stress 0.09400544 
    ## Run 403 stress 0.09969991 
    ## Run 404 stress 0.09268324 
    ## Run 405 stress 0.09380585 
    ## Run 406 stress 0.09030393 
    ## Run 407 stress 0.09286084 
    ## Run 408 stress 0.09268354 
    ## Run 409 stress 0.09407975 
    ## Run 410 stress 0.08503475 
    ## Run 411 stress 0.09030404 
    ## Run 412 stress 0.09381844 
    ## Run 413 stress 0.09416146 
    ## Run 414 stress 0.09159083 
    ## Run 415 stress 0.09407996 
    ## Run 416 stress 0.09610753 
    ## Run 417 stress 0.09159087 
    ## Run 418 stress 0.08503496 
    ## Run 419 stress 0.09539197 
    ## Run 420 stress 0.09407987 
    ## Run 421 stress 0.08503575 
    ## Run 422 stress 0.08773466 
    ## Run 423 stress 0.09159085 
    ## Run 424 stress 0.08503497 
    ## Run 425 stress 0.09374289 
    ## Run 426 stress 0.09168934 
    ## Run 427 stress 0.08973872 
    ## Run 428 stress 0.09159087 
    ## Run 429 stress 0.09168949 
    ## Run 430 stress 0.09286089 
    ## Run 431 stress 0.08503462 
    ## Run 432 stress 0.09268343 
    ## Run 433 stress 0.08440262 
    ## ... Procrustes: rmse 0.000246174  max resid 0.0004645553 
    ## ... Similar to previous best
    ## Run 434 stress 0.08440259 
    ## ... Procrustes: rmse 0.0002152261  max resid 0.0003955803 
    ## ... Similar to previous best
    ## Run 435 stress 0.1038034 
    ## Run 436 stress 0.09145325 
    ## Run 437 stress 0.08973866 
    ## Run 438 stress 0.09337279 
    ## Run 439 stress 0.09969981 
    ## Run 440 stress 0.09374319 
    ## Run 441 stress 0.0877347 
    ## Run 442 stress 0.08503504 
    ## Run 443 stress 0.08503508 
    ## Run 444 stress 0.09380578 
    ## Run 445 stress 0.08773484 
    ## Run 446 stress 0.09380577 
    ## Run 447 stress 0.09407969 
    ## Run 448 stress 0.09407991 
    ## Run 449 stress 0.0877349 
    ## Run 450 stress 0.09030393 
    ## Run 451 stress 0.09337294 
    ## Run 452 stress 0.09380588 
    ## Run 453 stress 0.09380591 
    ## Run 454 stress 0.0940799 
    ## Run 455 stress 0.09159088 
    ## Run 456 stress 0.09465916 
    ## Run 457 stress 0.09337322 
    ## Run 458 stress 0.08440279 
    ## ... Procrustes: rmse 0.0003590937  max resid 0.000647125 
    ## ... Similar to previous best
    ## Run 459 stress 0.08440267 
    ## ... Procrustes: rmse 0.0002753163  max resid 0.0005595978 
    ## ... Similar to previous best
    ## Run 460 stress 0.09969967 
    ## Run 461 stress 0.09159083 
    ## Run 462 stress 0.1052847 
    ## Run 463 stress 0.1052848 
    ## Run 464 stress 0.08773488 
    ## Run 465 stress 0.09445363 
    ## Run 466 stress 0.09969989 
    ## Run 467 stress 0.08503503 
    ## Run 468 stress 0.09145321 
    ## Run 469 stress 0.1038052 
    ## Run 470 stress 0.0850352 
    ## Run 471 stress 0.2857908 
    ## Run 472 stress 0.09464451 
    ## Run 473 stress 0.08973885 
    ## Run 474 stress 0.09590835 
    ## Run 475 stress 0.09407975 
    ## Run 476 stress 0.09159092 
    ## Run 477 stress 0.1026215 
    ## Run 478 stress 0.1038052 
    ## Run 479 stress 0.09445486 
    ## Run 480 stress 0.09407981 
    ## Run 481 stress 0.1038039 
    ## Run 482 stress 0.09145335 
    ## Run 483 stress 0.09145325 
    ## Run 484 stress 0.1052849 
    ## Run 485 stress 0.09400521 
    ## Run 486 stress 0.0897387 
    ## Run 487 stress 0.08773471 
    ## Run 488 stress 0.09374285 
    ## Run 489 stress 0.0940736 
    ## Run 490 stress 0.09417773 
    ## Run 491 stress 0.09407974 
    ## Run 492 stress 0.0937429 
    ## Run 493 stress 0.0953541 
    ## Run 494 stress 0.1052849 
    ## Run 495 stress 0.09268382 
    ## Run 496 stress 0.09465906 
    ## Run 497 stress 0.2743247 
    ## Run 498 stress 0.09321558 
    ## Run 499 stress 0.08973878 
    ## Run 500 stress 0.09268352 
    ## *** Best solution repeated 6 times

``` r
# Ocean sites and mixed lakes
SD_beta_env_OM_NMDS <- metaMDS(SD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.06623458 
    ## Run 1 stress 0.06477855 
    ## ... New best solution
    ## ... Procrustes: rmse 0.1519853  max resid 0.2644943 
    ## Run 2 stress 0.06623457 
    ## Run 3 stress 0.08226167 
    ## Run 4 stress 0.07411945 
    ## Run 5 stress 0.064779 
    ## ... Procrustes: rmse 0.0003303203  max resid 0.0005524336 
    ## ... Similar to previous best
    ## Run 6 stress 0.07411944 
    ## Run 7 stress 0.06623458 
    ## Run 8 stress 0.08226163 
    ## Run 9 stress 0.07411941 
    ## Run 10 stress 0.08226176 
    ## Run 11 stress 0.06623459 
    ## Run 12 stress 0.06477983 
    ## ... Procrustes: rmse 0.002415179  max resid 0.004058825 
    ## ... Similar to previous best
    ## Run 13 stress 0.06477879 
    ## ... Procrustes: rmse 0.001685672  max resid 0.002838168 
    ## ... Similar to previous best
    ## Run 14 stress 0.06124512 
    ## ... New best solution
    ## ... Procrustes: rmse 0.09173862  max resid 0.2199517 
    ## Run 15 stress 0.07166974 
    ## Run 16 stress 0.06124527 
    ## ... Procrustes: rmse 0.0009076285  max resid 0.002415224 
    ## ... Similar to previous best
    ## Run 17 stress 0.06124511 
    ## ... New best solution
    ## ... Procrustes: rmse 5.776458e-06  max resid 1.304584e-05 
    ## ... Similar to previous best
    ## Run 18 stress 0.06623458 
    ## Run 19 stress 0.07166974 
    ## Run 20 stress 0.3252872 
    ## Run 21 stress 0.07166969 
    ## Run 22 stress 0.06124503 
    ## ... New best solution
    ## ... Procrustes: rmse 7.211241e-05  max resid 0.000190978 
    ## ... Similar to previous best
    ## Run 23 stress 0.06623458 
    ## Run 24 stress 0.07166985 
    ## Run 25 stress 0.06477938 
    ## Run 26 stress 0.07411942 
    ## Run 27 stress 0.06477889 
    ## Run 28 stress 0.06477963 
    ## Run 29 stress 0.0741194 
    ## Run 30 stress 0.06623462 
    ## Run 31 stress 0.06124511 
    ## ... Procrustes: rmse 7.054569e-05  max resid 0.0001866235 
    ## ... Similar to previous best
    ## Run 32 stress 0.07166988 
    ## Run 33 stress 0.07411941 
    ## Run 34 stress 0.07166969 
    ## Run 35 stress 0.0612452 
    ## ... Procrustes: rmse 0.0001334007  max resid 0.0003529631 
    ## ... Similar to previous best
    ## Run 36 stress 0.06477868 
    ## Run 37 stress 0.06623459 
    ## Run 38 stress 0.06124533 
    ## ... Procrustes: rmse 0.0001732963  max resid 0.0004465936 
    ## ... Similar to previous best
    ## Run 39 stress 0.07166977 
    ## Run 40 stress 0.07166981 
    ## Run 41 stress 0.06124506 
    ## ... Procrustes: rmse 2.621458e-05  max resid 6.366353e-05 
    ## ... Similar to previous best
    ## Run 42 stress 0.06123361 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01509874  max resid 0.04199709 
    ## Run 43 stress 0.07411942 
    ## Run 44 stress 0.06623465 
    ## Run 45 stress 0.06124506 
    ## ... Procrustes: rmse 0.0151201  max resid 0.04186978 
    ## Run 46 stress 0.06123361 
    ## ... New best solution
    ## ... Procrustes: rmse 5.775158e-05  max resid 7.210896e-05 
    ## ... Similar to previous best
    ## Run 47 stress 0.06477887 
    ## Run 48 stress 0.0647798 
    ## Run 49 stress 0.07411946 
    ## Run 50 stress 0.07411941 
    ## Run 51 stress 0.0647795 
    ## Run 52 stress 0.061245 
    ## ... Procrustes: rmse 0.01507281  max resid 0.04174412 
    ## Run 53 stress 0.07411939 
    ## Run 54 stress 0.0647783 
    ## Run 55 stress 0.06477841 
    ## Run 56 stress 0.06477863 
    ## Run 57 stress 0.07411949 
    ## Run 58 stress 0.07411941 
    ## Run 59 stress 0.06477934 
    ## Run 60 stress 0.06623458 
    ## Run 61 stress 0.06477879 
    ## Run 62 stress 0.06124506 
    ## ... Procrustes: rmse 0.01512874  max resid 0.04189404 
    ## Run 63 stress 0.06477927 
    ## Run 64 stress 0.0822616 
    ## Run 65 stress 0.07411938 
    ## Run 66 stress 0.08226156 
    ## Run 67 stress 0.06623457 
    ## Run 68 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 2.977613e-05  max resid 4.603933e-05 
    ## ... Similar to previous best
    ## Run 69 stress 0.06477967 
    ## Run 70 stress 0.06477917 
    ## Run 71 stress 0.07411948 
    ## Run 72 stress 0.07166986 
    ## Run 73 stress 0.07411951 
    ## Run 74 stress 0.07166974 
    ## Run 75 stress 0.0612336 
    ## ... New best solution
    ## ... Procrustes: rmse 3.21966e-06  max resid 6.281765e-06 
    ## ... Similar to previous best
    ## Run 76 stress 0.07411938 
    ## Run 77 stress 0.08226159 
    ## Run 78 stress 0.07411945 
    ## Run 79 stress 0.3157929 
    ## Run 80 stress 0.06477832 
    ## Run 81 stress 0.06623458 
    ## Run 82 stress 0.07166984 
    ## Run 83 stress 0.06123361 
    ## ... Procrustes: rmse 4.081732e-05  max resid 5.151863e-05 
    ## ... Similar to previous best
    ## Run 84 stress 0.06623457 
    ## Run 85 stress 0.0647794 
    ## Run 86 stress 0.06123362 
    ## ... Procrustes: rmse 1.636952e-05  max resid 3.547823e-05 
    ## ... Similar to previous best
    ## Run 87 stress 0.06623458 
    ## Run 88 stress 0.06477841 
    ## Run 89 stress 0.06623462 
    ## Run 90 stress 0.06623458 
    ## Run 91 stress 0.06623463 
    ## Run 92 stress 0.0612336 
    ## ... Procrustes: rmse 2.759163e-06  max resid 3.81696e-06 
    ## ... Similar to previous best
    ## Run 93 stress 0.06623464 
    ## Run 94 stress 0.07411938 
    ## Run 95 stress 0.06477843 
    ## Run 96 stress 0.07411942 
    ## Run 97 stress 0.06477984 
    ## Run 98 stress 0.0716697 
    ## Run 99 stress 0.06124526 
    ## ... Procrustes: rmse 0.01523742  max resid 0.04218595 
    ## Run 100 stress 0.0662346 
    ## Run 101 stress 0.07411944 
    ## Run 102 stress 0.07166973 
    ## Run 103 stress 0.06124501 
    ## ... Procrustes: rmse 0.01508419  max resid 0.04177243 
    ## Run 104 stress 0.06477827 
    ## Run 105 stress 0.07411941 
    ## Run 106 stress 0.06477846 
    ## Run 107 stress 0.08226159 
    ## Run 108 stress 0.08226169 
    ## Run 109 stress 0.2345379 
    ## Run 110 stress 0.07411938 
    ## Run 111 stress 0.06477872 
    ## Run 112 stress 0.06623457 
    ## Run 113 stress 0.06623458 
    ## Run 114 stress 0.08226155 
    ## Run 115 stress 0.07166973 
    ## Run 116 stress 0.07411952 
    ## Run 117 stress 0.08226156 
    ## Run 118 stress 0.06623457 
    ## Run 119 stress 0.06623458 
    ## Run 120 stress 0.07411945 
    ## Run 121 stress 0.0612336 
    ## ... Procrustes: rmse 8.450269e-06  max resid 1.079999e-05 
    ## ... Similar to previous best
    ## Run 122 stress 0.06477871 
    ## Run 123 stress 0.06123363 
    ## ... Procrustes: rmse 8.669457e-05  max resid 0.0001204347 
    ## ... Similar to previous best
    ## Run 124 stress 0.06124504 
    ## ... Procrustes: rmse 0.01511483  max resid 0.04185424 
    ## Run 125 stress 0.06124524 
    ## ... Procrustes: rmse 0.01524605  max resid 0.04220821 
    ## Run 126 stress 0.06623458 
    ## Run 127 stress 0.06123361 
    ## ... Procrustes: rmse 3.634634e-05  max resid 6.52621e-05 
    ## ... Similar to previous best
    ## Run 128 stress 0.07411951 
    ## Run 129 stress 0.2592302 
    ## Run 130 stress 0.07411949 
    ## Run 131 stress 0.08226174 
    ## Run 132 stress 0.06623461 
    ## Run 133 stress 0.06477906 
    ## Run 134 stress 0.06623458 
    ## Run 135 stress 0.06124522 
    ## ... Procrustes: rmse 0.01525596  max resid 0.04223334 
    ## Run 136 stress 0.2394588 
    ## Run 137 stress 0.06124494 
    ## ... Procrustes: rmse 0.01454585  max resid 0.04033091 
    ## Run 138 stress 0.07166974 
    ## Run 139 stress 0.0662346 
    ## Run 140 stress 0.06477839 
    ## Run 141 stress 0.06623458 
    ## Run 142 stress 0.06477896 
    ## Run 143 stress 0.07166965 
    ## Run 144 stress 0.07411948 
    ## Run 145 stress 0.08226163 
    ## Run 146 stress 0.07411949 
    ## Run 147 stress 0.07411939 
    ## Run 148 stress 0.0612336 
    ## ... Procrustes: rmse 1.300371e-05  max resid 3.067548e-05 
    ## ... Similar to previous best
    ## Run 149 stress 0.06124503 
    ## ... Procrustes: rmse 0.01510273  max resid 0.04182167 
    ## Run 150 stress 0.06477864 
    ## Run 151 stress 0.07411939 
    ## Run 152 stress 0.07411948 
    ## Run 153 stress 0.06124526 
    ## ... Procrustes: rmse 0.01527081  max resid 0.04227068 
    ## Run 154 stress 0.06623466 
    ## Run 155 stress 0.0662346 
    ## Run 156 stress 0.06124507 
    ## ... Procrustes: rmse 0.01514691  max resid 0.04194135 
    ## Run 157 stress 0.06477845 
    ## Run 158 stress 0.07411944 
    ## Run 159 stress 0.06123361 
    ## ... Procrustes: rmse 1.781519e-05  max resid 3.928901e-05 
    ## ... Similar to previous best
    ## Run 160 stress 0.06123362 
    ## ... Procrustes: rmse 5.306393e-05  max resid 0.0001142109 
    ## ... Similar to previous best
    ## Run 161 stress 0.07411942 
    ## Run 162 stress 0.07166987 
    ## Run 163 stress 0.06123362 
    ## ... Procrustes: rmse 5.20652e-05  max resid 6.518658e-05 
    ## ... Similar to previous best
    ## Run 164 stress 0.06124501 
    ## ... Procrustes: rmse 0.01509162  max resid 0.04179332 
    ## Run 165 stress 0.07411938 
    ## Run 166 stress 0.0716697 
    ## Run 167 stress 0.06477914 
    ## Run 168 stress 0.06623457 
    ## Run 169 stress 0.07411941 
    ## Run 170 stress 0.07166982 
    ## Run 171 stress 0.06477843 
    ## Run 172 stress 0.08226164 
    ## Run 173 stress 0.06623457 
    ## Run 174 stress 0.06124495 
    ## ... Procrustes: rmse 0.01502109  max resid 0.0416042 
    ## Run 175 stress 0.08226156 
    ## Run 176 stress 0.06623458 
    ## Run 177 stress 0.07166967 
    ## Run 178 stress 0.07411942 
    ## Run 179 stress 0.06124522 
    ## ... Procrustes: rmse 0.01428658  max resid 0.03963386 
    ## Run 180 stress 0.06477994 
    ## Run 181 stress 0.06623457 
    ## Run 182 stress 0.0612336 
    ## ... Procrustes: rmse 1.305178e-05  max resid 1.86867e-05 
    ## ... Similar to previous best
    ## Run 183 stress 0.08226155 
    ## Run 184 stress 0.07166984 
    ## Run 185 stress 0.3407343 
    ## Run 186 stress 0.06623467 
    ## Run 187 stress 0.07411939 
    ## Run 188 stress 0.06477862 
    ## Run 189 stress 0.06623458 
    ## Run 190 stress 0.06477965 
    ## Run 191 stress 0.07166979 
    ## Run 192 stress 0.06623458 
    ## Run 193 stress 0.0612336 
    ## ... Procrustes: rmse 1.959118e-05  max resid 2.468421e-05 
    ## ... Similar to previous best
    ## Run 194 stress 0.06623463 
    ## Run 195 stress 0.08226163 
    ## Run 196 stress 0.07166971 
    ## Run 197 stress 0.07411942 
    ## Run 198 stress 0.07411947 
    ## Run 199 stress 0.07166966 
    ## Run 200 stress 0.06623464 
    ## Run 201 stress 0.06124494 
    ## ... Procrustes: rmse 0.01500878  max resid 0.04157194 
    ## Run 202 stress 0.07411938 
    ## Run 203 stress 0.0647794 
    ## Run 204 stress 0.06124512 
    ## ... Procrustes: rmse 0.01514775  max resid 0.04194071 
    ## Run 205 stress 0.06477996 
    ## Run 206 stress 0.06477854 
    ## Run 207 stress 0.06124494 
    ## ... Procrustes: rmse 0.01499713  max resid 0.04154068 
    ## Run 208 stress 0.07411946 
    ## Run 209 stress 0.06124508 
    ## ... Procrustes: rmse 0.01515204  max resid 0.04195507 
    ## Run 210 stress 0.07411944 
    ## Run 211 stress 0.06477945 
    ## Run 212 stress 0.06477972 
    ## Run 213 stress 0.0822616 
    ## Run 214 stress 0.0647785 
    ## Run 215 stress 0.07411939 
    ## Run 216 stress 0.0612336 
    ## ... Procrustes: rmse 9.88487e-06  max resid 1.328782e-05 
    ## ... Similar to previous best
    ## Run 217 stress 0.07166973 
    ## Run 218 stress 0.06123361 
    ## ... Procrustes: rmse 1.664931e-05  max resid 2.240063e-05 
    ## ... Similar to previous best
    ## Run 219 stress 0.07411943 
    ## Run 220 stress 0.06623457 
    ## Run 221 stress 0.07166969 
    ## Run 222 stress 0.06477831 
    ## Run 223 stress 0.07411943 
    ## Run 224 stress 0.06477841 
    ## Run 225 stress 0.07411955 
    ## Run 226 stress 0.06123364 
    ## ... Procrustes: rmse 2.74868e-05  max resid 6.322773e-05 
    ## ... Similar to previous best
    ## Run 227 stress 0.09300802 
    ## Run 228 stress 0.07411948 
    ## Run 229 stress 0.07166964 
    ## Run 230 stress 0.07411941 
    ## Run 231 stress 0.06124516 
    ## ... Procrustes: rmse 0.01520793  max resid 0.04210536 
    ## Run 232 stress 0.06477975 
    ## Run 233 stress 0.07411939 
    ## Run 234 stress 0.06477926 
    ## Run 235 stress 0.07411948 
    ## Run 236 stress 0.06124522 
    ## ... Procrustes: rmse 0.01428459  max resid 0.03962748 
    ## Run 237 stress 0.07166966 
    ## Run 238 stress 0.07411953 
    ## Run 239 stress 0.07166966 
    ## Run 240 stress 0.06477855 
    ## Run 241 stress 0.07166964 
    ## Run 242 stress 0.06477939 
    ## Run 243 stress 0.07411938 
    ## Run 244 stress 0.0741194 
    ## Run 245 stress 0.07411945 
    ## Run 246 stress 0.0612336 
    ## ... Procrustes: rmse 1.47447e-05  max resid 2.153314e-05 
    ## ... Similar to previous best
    ## Run 247 stress 0.06623461 
    ## Run 248 stress 0.06124488 
    ## ... Procrustes: rmse 0.01487061  max resid 0.04120185 
    ## Run 249 stress 0.06623462 
    ## Run 250 stress 0.07411947 
    ## Run 251 stress 0.06124498 
    ## ... Procrustes: rmse 0.01505303  max resid 0.04168947 
    ## Run 252 stress 0.06477904 
    ## Run 253 stress 0.06124502 
    ## ... Procrustes: rmse 0.01507484  max resid 0.04174969 
    ## Run 254 stress 0.0741194 
    ## Run 255 stress 0.06124505 
    ## ... Procrustes: rmse 0.0144246  max resid 0.04000539 
    ## Run 256 stress 0.06124532 
    ## ... Procrustes: rmse 0.0152571  max resid 0.04223148 
    ## Run 257 stress 0.0716699 
    ## Run 258 stress 0.0741194 
    ## Run 259 stress 0.07411948 
    ## Run 260 stress 0.06623457 
    ## Run 261 stress 0.06124495 
    ## ... Procrustes: rmse 0.01502683  max resid 0.04161965 
    ## Run 262 stress 0.07411946 
    ## Run 263 stress 0.07411947 
    ## Run 264 stress 0.06477938 
    ## Run 265 stress 0.07411948 
    ## Run 266 stress 0.08226169 
    ## Run 267 stress 0.07166966 
    ## Run 268 stress 0.0647796 
    ## Run 269 stress 0.06623464 
    ## Run 270 stress 0.06623457 
    ## Run 271 stress 0.06623462 
    ## Run 272 stress 0.06477919 
    ## Run 273 stress 0.06477953 
    ## Run 274 stress 0.06623458 
    ## Run 275 stress 0.07411958 
    ## Run 276 stress 0.06623462 
    ## Run 277 stress 0.06623462 
    ## Run 278 stress 0.06477899 
    ## Run 279 stress 0.0662346 
    ## Run 280 stress 0.07166971 
    ## Run 281 stress 0.2970709 
    ## Run 282 stress 0.06477873 
    ## Run 283 stress 0.06623458 
    ## Run 284 stress 0.06477896 
    ## Run 285 stress 0.07166975 
    ## Run 286 stress 0.07411938 
    ## Run 287 stress 0.2361811 
    ## Run 288 stress 0.07166966 
    ## Run 289 stress 0.06623457 
    ## Run 290 stress 0.0612336 
    ## ... Procrustes: rmse 2.077514e-06  max resid 2.934154e-06 
    ## ... Similar to previous best
    ## Run 291 stress 0.07411944 
    ## Run 292 stress 0.0662346 
    ## Run 293 stress 0.07411941 
    ## Run 294 stress 0.0612336 
    ## ... Procrustes: rmse 1.473956e-05  max resid 3.207805e-05 
    ## ... Similar to previous best
    ## Run 295 stress 0.06124498 
    ## ... Procrustes: rmse 0.01505128  max resid 0.04168414 
    ## Run 296 stress 0.07411943 
    ## Run 297 stress 0.06623459 
    ## Run 298 stress 0.06123362 
    ## ... Procrustes: rmse 4.860256e-05  max resid 8.579583e-05 
    ## ... Similar to previous best
    ## Run 299 stress 0.07411947 
    ## Run 300 stress 0.06623461 
    ## Run 301 stress 0.06477853 
    ## Run 302 stress 0.0741194 
    ## Run 303 stress 0.07411943 
    ## Run 304 stress 0.07411941 
    ## Run 305 stress 0.06623457 
    ## Run 306 stress 0.06623457 
    ## Run 307 stress 0.06477956 
    ## Run 308 stress 0.07166973 
    ## Run 309 stress 0.06623466 
    ## Run 310 stress 0.06623462 
    ## Run 311 stress 0.07411942 
    ## Run 312 stress 0.06124502 
    ## ... Procrustes: rmse 0.01509042  max resid 0.04179106 
    ## Run 313 stress 0.06477859 
    ## Run 314 stress 0.0741194 
    ## Run 315 stress 0.08226164 
    ## Run 316 stress 0.0612336 
    ## ... Procrustes: rmse 1.599095e-05  max resid 3.253697e-05 
    ## ... Similar to previous best
    ## Run 317 stress 0.07166971 
    ## Run 318 stress 0.08226159 
    ## Run 319 stress 0.06623459 
    ## Run 320 stress 0.06124495 
    ## ... Procrustes: rmse 0.01500973  max resid 0.04157353 
    ## Run 321 stress 0.0612452 
    ## ... Procrustes: rmse 0.01522084  max resid 0.0421363 
    ## Run 322 stress 0.06124507 
    ## ... Procrustes: rmse 0.01514461  max resid 0.04193485 
    ## Run 323 stress 0.07411939 
    ## Run 324 stress 0.07166986 
    ## Run 325 stress 0.06623458 
    ## Run 326 stress 0.07411939 
    ## Run 327 stress 0.08226159 
    ## Run 328 stress 0.06477862 
    ## Run 329 stress 0.3266727 
    ## Run 330 stress 0.06623457 
    ## Run 331 stress 0.06123361 
    ## ... Procrustes: rmse 1.513221e-05  max resid 1.973814e-05 
    ## ... Similar to previous best
    ## Run 332 stress 0.06477891 
    ## Run 333 stress 0.07166965 
    ## Run 334 stress 0.0716698 
    ## Run 335 stress 0.06623458 
    ## Run 336 stress 0.06623458 
    ## Run 337 stress 0.06477942 
    ## Run 338 stress 0.08226161 
    ## Run 339 stress 0.06477957 
    ## Run 340 stress 0.06623459 
    ## Run 341 stress 0.0741194 
    ## Run 342 stress 0.08226164 
    ## Run 343 stress 0.06123361 
    ## ... Procrustes: rmse 4.044063e-05  max resid 5.120787e-05 
    ## ... Similar to previous best
    ## Run 344 stress 0.06623462 
    ## Run 345 stress 0.327167 
    ## Run 346 stress 0.06477991 
    ## Run 347 stress 0.06623465 
    ## Run 348 stress 0.06623458 
    ## Run 349 stress 0.0612336 
    ## ... Procrustes: rmse 6.491974e-06  max resid 1.153444e-05 
    ## ... Similar to previous best
    ## Run 350 stress 0.06124525 
    ## ... Procrustes: rmse 0.01524629  max resid 0.04220918 
    ## Run 351 stress 0.06623465 
    ## Run 352 stress 0.0741195 
    ## Run 353 stress 0.07166965 
    ## Run 354 stress 0.06124513 
    ## ... Procrustes: rmse 0.01519715  max resid 0.04207492 
    ## Run 355 stress 0.07166965 
    ## Run 356 stress 0.06124511 
    ## ... Procrustes: rmse 0.01517658  max resid 0.04201978 
    ## Run 357 stress 0.07411952 
    ## Run 358 stress 0.06124528 
    ## ... Procrustes: rmse 0.01529255  max resid 0.04233071 
    ## Run 359 stress 0.2585391 
    ## Run 360 stress 0.0716699 
    ## Run 361 stress 0.07166964 
    ## Run 362 stress 0.06124509 
    ## ... Procrustes: rmse 0.01513836  max resid 0.04191565 
    ## Run 363 stress 0.06623459 
    ## Run 364 stress 0.0822616 
    ## Run 365 stress 0.2570957 
    ## Run 366 stress 0.07166964 
    ## Run 367 stress 0.07411939 
    ## Run 368 stress 0.06623462 
    ## Run 369 stress 0.07166964 
    ## Run 370 stress 0.06124517 
    ## ... Procrustes: rmse 0.01522132  max resid 0.04213945 
    ## Run 371 stress 0.06477849 
    ## Run 372 stress 0.06124507 
    ## ... Procrustes: rmse 0.01514173  max resid 0.04192617 
    ## Run 373 stress 0.0612451 
    ## ... Procrustes: rmse 0.01516078  max resid 0.04197941 
    ## Run 374 stress 0.07411938 
    ## Run 375 stress 0.07411946 
    ## Run 376 stress 0.07411948 
    ## Run 377 stress 0.06623459 
    ## Run 378 stress 0.06477928 
    ## Run 379 stress 0.06477881 
    ## Run 380 stress 0.07411939 
    ## Run 381 stress 0.08226175 
    ## Run 382 stress 0.0662346 
    ## Run 383 stress 0.07411938 
    ## Run 384 stress 0.06623463 
    ## Run 385 stress 0.07411949 
    ## Run 386 stress 0.0612449 
    ## ... Procrustes: rmse 0.0149371  max resid 0.04137917 
    ## Run 387 stress 0.06123361 
    ## ... Procrustes: rmse 4.316411e-05  max resid 5.287921e-05 
    ## ... Similar to previous best
    ## Run 388 stress 0.06623458 
    ## Run 389 stress 0.06124489 
    ## ... Procrustes: rmse 0.01493297  max resid 0.04136865 
    ## Run 390 stress 0.07166981 
    ## Run 391 stress 0.06623458 
    ## Run 392 stress 0.06477888 
    ## Run 393 stress 0.08226162 
    ## Run 394 stress 0.07166969 
    ## Run 395 stress 0.06477934 
    ## Run 396 stress 0.0612336 
    ## ... Procrustes: rmse 4.539708e-06  max resid 1.019625e-05 
    ## ... Similar to previous best
    ## Run 397 stress 0.07411941 
    ## Run 398 stress 0.07411954 
    ## Run 399 stress 0.06477973 
    ## Run 400 stress 0.0946895 
    ## Run 401 stress 0.07166985 
    ## Run 402 stress 0.07411942 
    ## Run 403 stress 0.06623463 
    ## Run 404 stress 0.06477844 
    ## Run 405 stress 0.06623463 
    ## Run 406 stress 0.06124497 
    ## ... Procrustes: rmse 0.01503978  max resid 0.04165418 
    ## Run 407 stress 0.06623458 
    ## Run 408 stress 0.08226157 
    ## Run 409 stress 0.07166966 
    ## Run 410 stress 0.0612336 
    ## ... Procrustes: rmse 7.937374e-06  max resid 1.337337e-05 
    ## ... Similar to previous best
    ## Run 411 stress 0.08226159 
    ## Run 412 stress 0.06623462 
    ## Run 413 stress 0.07411941 
    ## Run 414 stress 0.07166968 
    ## Run 415 stress 0.06477865 
    ## Run 416 stress 0.07411942 
    ## Run 417 stress 0.06124505 
    ## ... Procrustes: rmse 0.01512991  max resid 0.04189531 
    ## Run 418 stress 0.07411938 
    ## Run 419 stress 0.06477856 
    ## Run 420 stress 0.06477888 
    ## Run 421 stress 0.0612336 
    ## ... Procrustes: rmse 1.893549e-05  max resid 2.42851e-05 
    ## ... Similar to previous best
    ## Run 422 stress 0.2518595 
    ## Run 423 stress 0.0612336 
    ## ... Procrustes: rmse 6.155427e-06  max resid 9.376327e-06 
    ## ... Similar to previous best
    ## Run 424 stress 0.07411953 
    ## Run 425 stress 0.06623461 
    ## Run 426 stress 0.06477933 
    ## Run 427 stress 0.06124488 
    ## ... Procrustes: rmse 0.01490335  max resid 0.04128855 
    ## Run 428 stress 0.07166964 
    ## Run 429 stress 0.06477929 
    ## Run 430 stress 0.06477904 
    ## Run 431 stress 0.07166973 
    ## Run 432 stress 0.07166987 
    ## Run 433 stress 0.07166982 
    ## Run 434 stress 0.06124505 
    ## ... Procrustes: rmse 0.01511756  max resid 0.04186335 
    ## Run 435 stress 0.2324715 
    ## Run 436 stress 0.07166972 
    ## Run 437 stress 0.07166983 
    ## Run 438 stress 0.07411939 
    ## Run 439 stress 0.07166982 
    ## Run 440 stress 0.0612336 
    ## ... Procrustes: rmse 1.344305e-05  max resid 1.631575e-05 
    ## ... Similar to previous best
    ## Run 441 stress 0.07166971 
    ## Run 442 stress 0.06477838 
    ## Run 443 stress 0.06477889 
    ## Run 444 stress 0.07411947 
    ## Run 445 stress 0.06623457 
    ## Run 446 stress 0.06477912 
    ## Run 447 stress 0.06477838 
    ## Run 448 stress 0.07411938 
    ## Run 449 stress 0.0612336 
    ## ... Procrustes: rmse 1.849545e-05  max resid 2.408676e-05 
    ## ... Similar to previous best
    ## Run 450 stress 0.07411938 
    ## Run 451 stress 0.07166988 
    ## Run 452 stress 0.06477877 
    ## Run 453 stress 0.06623459 
    ## Run 454 stress 0.06623458 
    ## Run 455 stress 0.06124512 
    ## ... Procrustes: rmse 0.01518334  max resid 0.04203721 
    ## Run 456 stress 0.07411938 
    ## Run 457 stress 0.06123364 
    ## ... Procrustes: rmse 9.57605e-05  max resid 0.0001341217 
    ## ... Similar to previous best
    ## Run 458 stress 0.07166967 
    ## Run 459 stress 0.06623461 
    ## Run 460 stress 0.07166982 
    ## Run 461 stress 0.07411943 
    ## Run 462 stress 0.06623462 
    ## Run 463 stress 0.06124486 
    ## ... Procrustes: rmse 0.01482005  max resid 0.04106518 
    ## Run 464 stress 0.3119309 
    ## Run 465 stress 0.06477861 
    ## Run 466 stress 0.08226171 
    ## Run 467 stress 0.06477844 
    ## Run 468 stress 0.06477892 
    ## Run 469 stress 0.3157923 
    ## Run 470 stress 0.06124522 
    ## ... Procrustes: rmse 0.0152557  max resid 0.04223247 
    ## Run 471 stress 0.07411946 
    ## Run 472 stress 0.07411945 
    ## Run 473 stress 0.0741194 
    ## Run 474 stress 0.2381016 
    ## Run 475 stress 0.06477947 
    ## Run 476 stress 0.07166978 
    ## Run 477 stress 0.06623458 
    ## Run 478 stress 0.3407301 
    ## Run 479 stress 0.08226159 
    ## Run 480 stress 0.0647794 
    ## Run 481 stress 0.06477867 
    ## Run 482 stress 0.06623458 
    ## Run 483 stress 0.0612336 
    ## ... Procrustes: rmse 7.378157e-06  max resid 9.27337e-06 
    ## ... Similar to previous best
    ## Run 484 stress 0.08226167 
    ## Run 485 stress 0.07411944 
    ## Run 486 stress 0.06477886 
    ## Run 487 stress 0.06123362 
    ## ... Procrustes: rmse 3.726033e-05  max resid 8.685544e-05 
    ## ... Similar to previous best
    ## Run 488 stress 0.06623457 
    ## Run 489 stress 0.2394588 
    ## Run 490 stress 0.07166981 
    ## Run 491 stress 0.2361811 
    ## Run 492 stress 0.06477859 
    ## Run 493 stress 0.08226163 
    ## Run 494 stress 0.07166965 
    ## Run 495 stress 0.06477837 
    ## Run 496 stress 0.07166971 
    ## Run 497 stress 0.07166976 
    ## Run 498 stress 0.07166981 
    ## Run 499 stress 0.06124517 
    ## ... Procrustes: rmse 0.0152197  max resid 0.04213477 
    ## Run 500 stress 0.2361811 
    ## *** Best solution repeated 34 times

``` r
# Stratified lakes and ocean sites
SD_beta_env_SO_NMDS <- metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 9.59483e-05 
    ## Run 1 stress 9.455847e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001443333  max resid 0.0002935966 
    ## ... Similar to previous best
    ## Run 2 stress 9.411163e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002043952  max resid 0.0003646894 
    ## ... Similar to previous best
    ## Run 3 stress 9.307744e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001767588  max resid 0.0004454481 
    ## ... Similar to previous best
    ## Run 4 stress 9.669822e-05 
    ## ... Procrustes: rmse 0.0002084912  max resid 0.0004650377 
    ## ... Similar to previous best
    ## Run 5 stress 8.891625e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000171131  max resid 0.0002699598 
    ## ... Similar to previous best
    ## Run 6 stress 9.266009e-05 
    ## ... Procrustes: rmse 0.0001029227  max resid 0.0001806406 
    ## ... Similar to previous best
    ## Run 7 stress 9.560392e-05 
    ## ... Procrustes: rmse 0.0001929886  max resid 0.0002998487 
    ## ... Similar to previous best
    ## Run 8 stress 9.92592e-05 
    ## ... Procrustes: rmse 0.0002214553  max resid 0.0003229699 
    ## ... Similar to previous best
    ## Run 9 stress 9.601467e-05 
    ## ... Procrustes: rmse 0.0001563044  max resid 0.0002516388 
    ## ... Similar to previous best
    ## Run 10 stress 9.306288e-05 
    ## ... Procrustes: rmse 0.0001508612  max resid 0.0002283014 
    ## ... Similar to previous best
    ## Run 11 stress 9.860299e-05 
    ## ... Procrustes: rmse 0.0001724979  max resid 0.0003835142 
    ## ... Similar to previous best
    ## Run 12 stress 9.748299e-05 
    ## ... Procrustes: rmse 0.0001415091  max resid 0.0002523283 
    ## ... Similar to previous best
    ## Run 13 stress 9.80262e-05 
    ## ... Procrustes: rmse 0.0001574431  max resid 0.0003635146 
    ## ... Similar to previous best
    ## Run 14 stress 9.62671e-05 
    ## ... Procrustes: rmse 0.000191129  max resid 0.000391342 
    ## ... Similar to previous best
    ## Run 15 stress 8.724366e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001480807  max resid 0.0003402862 
    ## ... Similar to previous best
    ## Run 16 stress 9.600817e-05 
    ## ... Procrustes: rmse 0.0001808388  max resid 0.0002676301 
    ## ... Similar to previous best
    ## Run 17 stress 9.877251e-05 
    ## ... Procrustes: rmse 0.0001828315  max resid 0.0003742204 
    ## ... Similar to previous best
    ## Run 18 stress 9.420667e-05 
    ## ... Procrustes: rmse 0.000170883  max resid 0.000368362 
    ## ... Similar to previous best
    ## Run 19 stress 8.460642e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001058656  max resid 0.0002272886 
    ## ... Similar to previous best
    ## Run 20 stress 9.625627e-05 
    ## ... Procrustes: rmse 0.0002093575  max resid 0.0003927585 
    ## ... Similar to previous best
    ## Run 21 stress 9.711649e-05 
    ## ... Procrustes: rmse 0.0001401306  max resid 0.0002867668 
    ## ... Similar to previous best
    ## Run 22 stress 9.905063e-05 
    ## ... Procrustes: rmse 0.0002179253  max resid 0.0003958363 
    ## ... Similar to previous best
    ## Run 23 stress 9.639847e-05 
    ## ... Procrustes: rmse 0.000219322  max resid 0.0003947442 
    ## ... Similar to previous best
    ## Run 24 stress 9.043128e-05 
    ## ... Procrustes: rmse 0.0001585433  max resid 0.0003627034 
    ## ... Similar to previous best
    ## Run 25 stress 9.879227e-05 
    ## ... Procrustes: rmse 0.0001920482  max resid 0.0003694536 
    ## ... Similar to previous best
    ## Run 26 stress 9.562586e-05 
    ## ... Procrustes: rmse 0.0001698123  max resid 0.0003017524 
    ## ... Similar to previous best
    ## Run 27 stress 8.402542e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001543614  max resid 0.0002428775 
    ## ... Similar to previous best
    ## Run 28 stress 9.288763e-05 
    ## ... Procrustes: rmse 0.0001651537  max resid 0.0003121072 
    ## ... Similar to previous best
    ## Run 29 stress 9.131959e-05 
    ## ... Procrustes: rmse 0.0001702599  max resid 0.0003069595 
    ## ... Similar to previous best
    ## Run 30 stress 9.281136e-05 
    ## ... Procrustes: rmse 0.0001221327  max resid 0.0002389288 
    ## ... Similar to previous best
    ## Run 31 stress 9.000392e-05 
    ## ... Procrustes: rmse 0.0001307518  max resid 0.0003225383 
    ## ... Similar to previous best
    ## Run 32 stress 9.785087e-05 
    ## ... Procrustes: rmse 0.0001725174  max resid 0.000273779 
    ## ... Similar to previous best
    ## Run 33 stress 9.701078e-05 
    ## ... Procrustes: rmse 0.0001254419  max resid 0.0002580844 
    ## ... Similar to previous best
    ## Run 34 stress 9.310346e-05 
    ## ... Procrustes: rmse 0.0001683245  max resid 0.0003097519 
    ## ... Similar to previous best
    ## Run 35 stress 9.475704e-05 
    ## ... Procrustes: rmse 0.0001682172  max resid 0.0004218254 
    ## ... Similar to previous best
    ## Run 36 stress 9.731708e-05 
    ## ... Procrustes: rmse 3.264527e-05  max resid 5.343998e-05 
    ## ... Similar to previous best
    ## Run 37 stress 9.67596e-05 
    ## ... Procrustes: rmse 4.452161e-05  max resid 0.0001049383 
    ## ... Similar to previous best
    ## Run 38 stress 9.330914e-05 
    ## ... Procrustes: rmse 0.0001931843  max resid 0.0004000606 
    ## ... Similar to previous best
    ## Run 39 stress 9.421592e-05 
    ## ... Procrustes: rmse 0.0001714127  max resid 0.0004240388 
    ## ... Similar to previous best
    ## Run 40 stress 9.742081e-05 
    ## ... Procrustes: rmse 0.0001384379  max resid 0.0003365722 
    ## ... Similar to previous best
    ## Run 41 stress 9.656447e-05 
    ## ... Procrustes: rmse 0.0001673147  max resid 0.000267071 
    ## ... Similar to previous best
    ## Run 42 stress 9.368289e-05 
    ## ... Procrustes: rmse 0.0001963248  max resid 0.0004286616 
    ## ... Similar to previous best
    ## Run 43 stress 9.919176e-05 
    ## ... Procrustes: rmse 0.0001718849  max resid 0.0004032921 
    ## ... Similar to previous best
    ## Run 44 stress 9.26123e-05 
    ## ... Procrustes: rmse 0.0001185468  max resid 0.0002524491 
    ## ... Similar to previous best
    ## Run 45 stress 9.557455e-05 
    ## ... Procrustes: rmse 0.000162093  max resid 0.0003298136 
    ## ... Similar to previous best
    ## Run 46 stress 9.551216e-05 
    ## ... Procrustes: rmse 0.0001980091  max resid 0.0003608028 
    ## ... Similar to previous best
    ## Run 47 stress 9.035287e-05 
    ## ... Procrustes: rmse 0.0001257833  max resid 0.0002643044 
    ## ... Similar to previous best
    ## Run 48 stress 9.98322e-05 
    ## ... Procrustes: rmse 0.0002060193  max resid 0.0003532284 
    ## ... Similar to previous best
    ## Run 49 stress 9.499464e-05 
    ## ... Procrustes: rmse 2.398565e-05  max resid 4.522667e-05 
    ## ... Similar to previous best
    ## Run 50 stress 9.375509e-05 
    ## ... Procrustes: rmse 0.000134253  max resid 0.0003286451 
    ## ... Similar to previous best
    ## Run 51 stress 8.931708e-05 
    ## ... Procrustes: rmse 0.0001187339  max resid 0.0002564043 
    ## ... Similar to previous best
    ## Run 52 stress 9.46867e-05 
    ## ... Procrustes: rmse 0.0001941618  max resid 0.0004255914 
    ## ... Similar to previous best
    ## Run 53 stress 9.616492e-05 
    ## ... Procrustes: rmse 3.004443e-05  max resid 4.633363e-05 
    ## ... Similar to previous best
    ## Run 54 stress 9.912896e-05 
    ## ... Procrustes: rmse 0.0001418082  max resid 0.0003393672 
    ## ... Similar to previous best
    ## Run 55 stress 9.45852e-05 
    ## ... Procrustes: rmse 0.0002032908  max resid 0.0004310406 
    ## ... Similar to previous best
    ## Run 56 stress 9.111103e-05 
    ## ... Procrustes: rmse 0.0001669283  max resid 0.0004153249 
    ## ... Similar to previous best
    ## Run 57 stress 9.754871e-05 
    ## ... Procrustes: rmse 0.0001666058  max resid 0.0003379202 
    ## ... Similar to previous best
    ## Run 58 stress 9.512474e-05 
    ## ... Procrustes: rmse 0.0001346975  max resid 0.0003300561 
    ## ... Similar to previous best
    ## Run 59 stress 9.502472e-05 
    ## ... Procrustes: rmse 0.0001952774  max resid 0.0003860255 
    ## ... Similar to previous best
    ## Run 60 stress 9.69423e-05 
    ## ... Procrustes: rmse 6.607953e-05  max resid 0.0001375191 
    ## ... Similar to previous best
    ## Run 61 stress 9.470401e-05 
    ## ... Procrustes: rmse 0.0001724199  max resid 0.0004240569 
    ## ... Similar to previous best
    ## Run 62 stress 9.760873e-05 
    ## ... Procrustes: rmse 0.0001393702  max resid 0.0003168855 
    ## ... Similar to previous best
    ## Run 63 stress 9.970798e-05 
    ## ... Procrustes: rmse 0.0001389697  max resid 0.0003359363 
    ## ... Similar to previous best
    ## Run 64 stress 9.984759e-05 
    ## ... Procrustes: rmse 0.0001427994  max resid 0.0003416057 
    ## ... Similar to previous best
    ## Run 65 stress 9.708497e-05 
    ## ... Procrustes: rmse 0.000203766  max resid 0.0004001227 
    ## ... Similar to previous best
    ## Run 66 stress 9.955665e-05 
    ## ... Procrustes: rmse 0.0001424476  max resid 0.0003421156 
    ## ... Similar to previous best
    ## Run 67 stress 9.978696e-05 
    ## ... Procrustes: rmse 0.0002108318  max resid 0.0003568588 
    ## ... Similar to previous best
    ## Run 68 stress 9.228539e-05 
    ## ... Procrustes: rmse 0.0001724966  max resid 0.0003459117 
    ## ... Similar to previous best
    ## Run 69 stress 9.559428e-05 
    ## ... Procrustes: rmse 0.0001993662  max resid 0.0003929546 
    ## ... Similar to previous best
    ## Run 70 stress 9.551463e-05 
    ## ... Procrustes: rmse 0.0001277619  max resid 0.0002224988 
    ## ... Similar to previous best
    ## Run 71 stress 9.163447e-05 
    ## ... Procrustes: rmse 8.73858e-05  max resid 0.000153363 
    ## ... Similar to previous best
    ## Run 72 stress 9.948369e-05 
    ## ... Procrustes: rmse 0.0002078837  max resid 0.0003474614 
    ## ... Similar to previous best
    ## Run 73 stress 8.631643e-05 
    ## ... Procrustes: rmse 0.0001283361  max resid 0.0003186986 
    ## ... Similar to previous best
    ## Run 74 stress 9.73531e-05 
    ## ... Procrustes: rmse 0.000130076  max resid 0.0002660222 
    ## ... Similar to previous best
    ## Run 75 stress 9.537926e-05 
    ## ... Procrustes: rmse 0.0001912296  max resid 0.0003479111 
    ## ... Similar to previous best
    ## Run 76 stress 9.955952e-05 
    ## ... Procrustes: rmse 0.0002054251  max resid 0.0003620061 
    ## ... Similar to previous best
    ## Run 77 stress 9.069071e-05 
    ## ... Procrustes: rmse 0.0001337344  max resid 0.0003220301 
    ## ... Similar to previous best
    ## Run 78 stress 9.625243e-05 
    ## ... Procrustes: rmse 0.0002078939  max resid 0.0003770002 
    ## ... Similar to previous best
    ## Run 79 stress 9.43522e-05 
    ## ... Procrustes: rmse 0.000130907  max resid 0.0003188632 
    ## ... Similar to previous best
    ## Run 80 stress 9.764172e-05 
    ## ... Procrustes: rmse 2.446296e-05  max resid 4.450332e-05 
    ## ... Similar to previous best
    ## Run 81 stress 9.732611e-05 
    ## ... Procrustes: rmse 0.0002044609  max resid 0.0003183692 
    ## ... Similar to previous best
    ## Run 82 stress 8.596233e-05 
    ## ... Procrustes: rmse 0.0002099514  max resid 0.0003411903 
    ## ... Similar to previous best
    ## Run 83 stress 9.647804e-05 
    ## ... Procrustes: rmse 3.472645e-05  max resid 6.487574e-05 
    ## ... Similar to previous best
    ## Run 84 stress 8.901817e-05 
    ## ... Procrustes: rmse 8.570181e-05  max resid 0.0002077554 
    ## ... Similar to previous best
    ## Run 85 stress 9.657595e-05 
    ## ... Procrustes: rmse 0.0001313881  max resid 0.000266041 
    ## ... Similar to previous best
    ## Run 86 stress 9.258794e-05 
    ## ... Procrustes: rmse 0.0001162144  max resid 0.0002301379 
    ## ... Similar to previous best
    ## Run 87 stress 9.636139e-05 
    ## ... Procrustes: rmse 0.0001982249  max resid 0.0004093983 
    ## ... Similar to previous best
    ## Run 88 stress 9.423101e-05 
    ## ... Procrustes: rmse 0.0001331579  max resid 0.0003257354 
    ## ... Similar to previous best
    ## Run 89 stress 9.479653e-05 
    ## ... Procrustes: rmse 0.000155177  max resid 0.0002522108 
    ## ... Similar to previous best
    ## Run 90 stress 9.793773e-05 
    ## ... Procrustes: rmse 3.167593e-05  max resid 5.032691e-05 
    ## ... Similar to previous best
    ## Run 91 stress 9.726563e-05 
    ## ... Procrustes: rmse 0.0001615147  max resid 0.0002596691 
    ## ... Similar to previous best
    ## Run 92 stress 9.571764e-05 
    ## ... Procrustes: rmse 0.000167028  max resid 0.0003012138 
    ## ... Similar to previous best
    ## Run 93 stress 9.332454e-05 
    ## ... Procrustes: rmse 0.0001758325  max resid 0.0004037932 
    ## ... Similar to previous best
    ## Run 94 stress 9.690245e-05 
    ## ... Procrustes: rmse 0.0001368667  max resid 0.0003298052 
    ## ... Similar to previous best
    ## Run 95 stress 9.294546e-05 
    ## ... Procrustes: rmse 9.324895e-05  max resid 0.0001954456 
    ## ... Similar to previous best
    ## Run 96 stress 9.318519e-05 
    ## ... Procrustes: rmse 0.0001349022  max resid 0.0003340633 
    ## ... Similar to previous best
    ## Run 97 stress 9.059695e-05 
    ## ... Procrustes: rmse 0.0001960379  max resid 0.0003367408 
    ## ... Similar to previous best
    ## Run 98 stress 9.583616e-05 
    ## ... Procrustes: rmse 2.896557e-05  max resid 4.272422e-05 
    ## ... Similar to previous best
    ## Run 99 stress 9.393535e-05 
    ## ... Procrustes: rmse 0.0001971287  max resid 0.0004302324 
    ## ... Similar to previous best
    ## Run 100 stress 9.720536e-05 
    ## ... Procrustes: rmse 0.0001661817  max resid 0.0003376544 
    ## ... Similar to previous best
    ## Run 101 stress 9.166726e-05 
    ## ... Procrustes: rmse 0.0001293305  max resid 0.0003199365 
    ## ... Similar to previous best
    ## Run 102 stress 9.705555e-05 
    ## ... Procrustes: rmse 0.0002010285  max resid 0.0003572257 
    ## ... Similar to previous best
    ## Run 103 stress 9.440921e-05 
    ## ... Procrustes: rmse 9.500733e-05  max resid 0.0002022245 
    ## ... Similar to previous best
    ## Run 104 stress 9.372342e-05 
    ## ... Procrustes: rmse 0.0001660347  max resid 0.0003080207 
    ## ... Similar to previous best
    ## Run 105 stress 7.608425e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001762485  max resid 0.0003123027 
    ## ... Similar to previous best
    ## Run 106 stress 9.654072e-05 
    ## ... Procrustes: rmse 6.227813e-05  max resid 0.0001038337 
    ## ... Similar to previous best
    ## Run 107 stress 9.741213e-05 
    ## ... Procrustes: rmse 8.942363e-05  max resid 0.0001525443 
    ## ... Similar to previous best
    ## Run 108 stress 9.254293e-05 
    ## ... Procrustes: rmse 7.01002e-05  max resid 0.0001360264 
    ## ... Similar to previous best
    ## Run 109 stress 9.30041e-05 
    ## ... Procrustes: rmse 6.146969e-05  max resid 0.0001161762 
    ## ... Similar to previous best
    ## Run 110 stress 9.964452e-05 
    ## ... Procrustes: rmse 0.0001817791  max resid 0.0003499273 
    ## ... Similar to previous best
    ## Run 111 stress 9.778108e-05 
    ## ... Procrustes: rmse 0.0001159718  max resid 0.000198664 
    ## ... Similar to previous best
    ## Run 112 stress 9.821393e-05 
    ## ... Procrustes: rmse 0.0001586393  max resid 0.000335031 
    ## ... Similar to previous best
    ## Run 113 stress 9.637443e-05 
    ## ... Procrustes: rmse 0.0001567166  max resid 0.0003054072 
    ## ... Similar to previous best
    ## Run 114 stress 9.478392e-05 
    ## ... Procrustes: rmse 0.0002197838  max resid 0.000362603 
    ## ... Similar to previous best
    ## Run 115 stress 9.654135e-05 
    ## ... Procrustes: rmse 0.0001528774  max resid 0.0003059549 
    ## ... Similar to previous best
    ## Run 116 stress 9.33775e-05 
    ## ... Procrustes: rmse 0.0001278736  max resid 0.0002575595 
    ## ... Similar to previous best
    ## Run 117 stress 9.075055e-05 
    ## ... Procrustes: rmse 0.0001821998  max resid 0.0003215393 
    ## ... Similar to previous best
    ## Run 118 stress 9.418814e-05 
    ## ... Procrustes: rmse 0.0001811216  max resid 0.0003239594 
    ## ... Similar to previous best
    ## Run 119 stress 8.690762e-05 
    ## ... Procrustes: rmse 0.0001143913  max resid 0.0001971633 
    ## ... Similar to previous best
    ## Run 120 stress 9.575336e-05 
    ## ... Procrustes: rmse 0.0001280738  max resid 0.0002598122 
    ## ... Similar to previous best
    ## Run 121 stress 6.84554e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001617037  max resid 0.0002766423 
    ## ... Similar to previous best
    ## Run 122 stress 9.032897e-05 
    ## ... Procrustes: rmse 0.0001383149  max resid 0.0002612843 
    ## ... Similar to previous best
    ## Run 123 stress 9.584009e-05 
    ## ... Procrustes: rmse 0.0001677448  max resid 0.0003733603 
    ## ... Similar to previous best
    ## Run 124 stress 9.801103e-05 
    ## ... Procrustes: rmse 0.0002110756  max resid 0.000414268 
    ## ... Similar to previous best
    ## Run 125 stress 9.755786e-05 
    ## ... Procrustes: rmse 0.000127644  max resid 0.0002782959 
    ## ... Similar to previous best
    ## Run 126 stress 9.921329e-05 
    ## ... Procrustes: rmse 0.0001626808  max resid 0.0002842446 
    ## ... Similar to previous best
    ## Run 127 stress 8.247948e-05 
    ## ... Procrustes: rmse 0.0001110778  max resid 0.0002606743 
    ## ... Similar to previous best
    ## Run 128 stress 9.887605e-05 
    ## ... Procrustes: rmse 0.0002115896  max resid 0.0004141266 
    ## ... Similar to previous best
    ## Run 129 stress 9.252997e-05 
    ## ... Procrustes: rmse 0.0001348198  max resid 0.0002925974 
    ## ... Similar to previous best
    ## Run 130 stress 8.829778e-05 
    ## ... Procrustes: rmse 0.0001571038  max resid 0.0002367451 
    ## ... Similar to previous best
    ## Run 131 stress 9.548417e-05 
    ## ... Procrustes: rmse 0.0001932863  max resid 0.0003689214 
    ## ... Similar to previous best
    ## Run 132 stress 9.494977e-05 
    ## ... Procrustes: rmse 0.0001638896  max resid 0.0003672721 
    ## ... Similar to previous best
    ## Run 133 stress 9.270122e-05 
    ## ... Procrustes: rmse 0.0001523995  max resid 0.0002353133 
    ## ... Similar to previous best
    ## Run 134 stress 9.894877e-05 
    ## ... Procrustes: rmse 0.0001587905  max resid 0.0003484404 
    ## ... Similar to previous best
    ## Run 135 stress 9.163621e-05 
    ## ... Procrustes: rmse 0.0001676294  max resid 0.0003175693 
    ## ... Similar to previous best
    ## Run 136 stress 9.340011e-05 
    ## ... Procrustes: rmse 7.107627e-05  max resid 9.518962e-05 
    ## ... Similar to previous best
    ## Run 137 stress 9.997751e-05 
    ## ... Procrustes: rmse 0.000152118  max resid 0.0002645327 
    ## ... Similar to previous best
    ## Run 138 stress 9.564183e-05 
    ## ... Procrustes: rmse 0.0001671839  max resid 0.0003757028 
    ## ... Similar to previous best
    ## Run 139 stress 9.303409e-05 
    ## ... Procrustes: rmse 0.0001530789  max resid 0.0002337864 
    ## ... Similar to previous best
    ## Run 140 stress 9.909035e-05 
    ## ... Procrustes: rmse 0.0001625377  max resid 0.0003041786 
    ## ... Similar to previous best
    ## Run 141 stress 9.56393e-05 
    ## ... Procrustes: rmse 0.0001950642  max resid 0.0003183652 
    ## ... Similar to previous best
    ## Run 142 stress 9.940568e-05 
    ## ... Procrustes: rmse 0.0001978777  max resid 0.0003222023 
    ## ... Similar to previous best
    ## Run 143 stress 9.920076e-05 
    ## ... Procrustes: rmse 0.0001281073  max resid 0.0001877916 
    ## ... Similar to previous best
    ## Run 144 stress 9.402419e-05 
    ## ... Procrustes: rmse 0.0001461608  max resid 0.000234053 
    ## ... Similar to previous best
    ## Run 145 stress 9.692129e-05 
    ## ... Procrustes: rmse 0.0002043301  max resid 0.0003183006 
    ## ... Similar to previous best
    ## Run 146 stress 9.715112e-05 
    ## ... Procrustes: rmse 0.0001502372  max resid 0.0002855654 
    ## ... Similar to previous best
    ## Run 147 stress 9.755305e-05 
    ## ... Procrustes: rmse 0.0001698568  max resid 0.0003842335 
    ## ... Similar to previous best
    ## Run 148 stress 8.35932e-05 
    ## ... Procrustes: rmse 0.000133639  max resid 0.0002588665 
    ## ... Similar to previous best
    ## Run 149 stress 8.313866e-05 
    ## ... Procrustes: rmse 9.753766e-05  max resid 0.0002406926 
    ## ... Similar to previous best
    ## Run 150 stress 9.940611e-05 
    ## ... Procrustes: rmse 7.837208e-05  max resid 0.0001125444 
    ## ... Similar to previous best
    ## Run 151 stress 9.335694e-05 
    ## ... Procrustes: rmse 0.0001618583  max resid 0.0003651571 
    ## ... Similar to previous best
    ## Run 152 stress 9.867131e-05 
    ## ... Procrustes: rmse 0.0001545569  max resid 0.0002412323 
    ## ... Similar to previous best
    ## Run 153 stress 9.407657e-05 
    ## ... Procrustes: rmse 0.0001879482  max resid 0.0003789425 
    ## ... Similar to previous best
    ## Run 154 stress 8.157217e-05 
    ## ... Procrustes: rmse 0.0001268484  max resid 0.0002255537 
    ## ... Similar to previous best
    ## Run 155 stress 9.547053e-05 
    ## ... Procrustes: rmse 0.0001734017  max resid 0.0003494899 
    ## ... Similar to previous best
    ## Run 156 stress 9.663957e-05 
    ## ... Procrustes: rmse 0.0001378763  max resid 0.0002779311 
    ## ... Similar to previous best
    ## Run 157 stress 9.283142e-05 
    ## ... Procrustes: rmse 0.0001389682  max resid 0.0002428478 
    ## ... Similar to previous best
    ## Run 158 stress 9.78567e-05 
    ## ... Procrustes: rmse 0.0001802901  max resid 0.000325408 
    ## ... Similar to previous best
    ## Run 159 stress 9.97529e-05 
    ## ... Procrustes: rmse 0.0001547692  max resid 0.0003302491 
    ## ... Similar to previous best
    ## Run 160 stress 8.358304e-05 
    ## ... Procrustes: rmse 0.0001470125  max resid 0.0003438059 
    ## ... Similar to previous best
    ## Run 161 stress 9.367769e-05 
    ## ... Procrustes: rmse 0.0001331676  max resid 0.0002889598 
    ## ... Similar to previous best
    ## Run 162 stress 9.890481e-05 
    ## ... Procrustes: rmse 0.0001449294  max resid 0.0002867028 
    ## ... Similar to previous best
    ## Run 163 stress 9.559847e-05 
    ## ... Procrustes: rmse 0.0001757317  max resid 0.0003230444 
    ## ... Similar to previous best
    ## Run 164 stress 9.715262e-05 
    ## ... Procrustes: rmse 0.0001389435  max resid 0.000282641 
    ## ... Similar to previous best
    ## Run 165 stress 9.191306e-05 
    ## ... Procrustes: rmse 0.0001175121  max resid 0.0002149525 
    ## ... Similar to previous best
    ## Run 166 stress 9.330033e-05 
    ## ... Procrustes: rmse 0.0001855179  max resid 0.0003093834 
    ## ... Similar to previous best
    ## Run 167 stress 9.512198e-05 
    ## ... Procrustes: rmse 0.0001446942  max resid 0.0002614607 
    ## ... Similar to previous best
    ## Run 168 stress 9.805073e-05 
    ## ... Procrustes: rmse 0.0001629892  max resid 0.0002450164 
    ## ... Similar to previous best
    ## Run 169 stress 9.205906e-05 
    ## ... Procrustes: rmse 0.0001625383  max resid 0.0002660476 
    ## ... Similar to previous best
    ## Run 170 stress 9.576749e-05 
    ## ... Procrustes: rmse 0.0001594917  max resid 0.0002448383 
    ## ... Similar to previous best
    ## Run 171 stress 9.368866e-05 
    ## ... Procrustes: rmse 7.305418e-05  max resid 0.0001024591 
    ## ... Similar to previous best
    ## Run 172 stress 9.345304e-05 
    ## ... Procrustes: rmse 0.0001599922  max resid 0.0002725414 
    ## ... Similar to previous best
    ## Run 173 stress 8.909401e-05 
    ## ... Procrustes: rmse 0.0001445358  max resid 0.0002669593 
    ## ... Similar to previous best
    ## Run 174 stress 9.041002e-05 
    ## ... Procrustes: rmse 0.0001317338  max resid 0.0002908391 
    ## ... Similar to previous best
    ## Run 175 stress 8.972582e-05 
    ## ... Procrustes: rmse 0.0001299242  max resid 0.0002212676 
    ## ... Similar to previous best
    ## Run 176 stress 9.805069e-05 
    ## ... Procrustes: rmse 0.00019371  max resid 0.0003176995 
    ## ... Similar to previous best
    ## Run 177 stress 9.410816e-05 
    ## ... Procrustes: rmse 0.0001556805  max resid 0.0003610116 
    ## ... Similar to previous best
    ## Run 178 stress 9.191387e-05 
    ## ... Procrustes: rmse 9.672023e-05  max resid 0.0001968928 
    ## ... Similar to previous best
    ## Run 179 stress 9.619749e-05 
    ## ... Procrustes: rmse 0.0001504918  max resid 0.0002883698 
    ## ... Similar to previous best
    ## Run 180 stress 9.331043e-05 
    ## ... Procrustes: rmse 0.0001517729  max resid 0.0002656255 
    ## ... Similar to previous best
    ## Run 181 stress 0.2280258 
    ## Run 182 stress 9.984261e-05 
    ## ... Procrustes: rmse 0.0002145222  max resid 0.0004174381 
    ## ... Similar to previous best
    ## Run 183 stress 9.792507e-05 
    ## ... Procrustes: rmse 0.0001132647  max resid 0.0002778867 
    ## ... Similar to previous best
    ## Run 184 stress 9.819366e-05 
    ## ... Procrustes: rmse 0.0001445215  max resid 0.0003030584 
    ## ... Similar to previous best
    ## Run 185 stress 9.878389e-05 
    ## ... Procrustes: rmse 0.0001542727  max resid 0.0003238357 
    ## ... Similar to previous best
    ## Run 186 stress 9.722039e-05 
    ## ... Procrustes: rmse 0.0001752882  max resid 0.0002626506 
    ## ... Similar to previous best
    ## Run 187 stress 9.927454e-05 
    ## ... Procrustes: rmse 8.018079e-05  max resid 0.0001142643 
    ## ... Similar to previous best
    ## Run 188 stress 9.929745e-05 
    ## ... Procrustes: rmse 0.0002025193  max resid 0.0003256308 
    ## ... Similar to previous best
    ## Run 189 stress 9.62028e-05 
    ## ... Procrustes: rmse 0.0001938408  max resid 0.0003690344 
    ## ... Similar to previous best
    ## Run 190 stress 9.817795e-05 
    ## ... Procrustes: rmse 8.091768e-05  max resid 0.0001214788 
    ## ... Similar to previous best
    ## Run 191 stress 9.442356e-05 
    ## ... Procrustes: rmse 0.0001540602  max resid 0.000237706 
    ## ... Similar to previous best
    ## Run 192 stress 9.885819e-05 
    ## ... Procrustes: rmse 0.0002094325  max resid 0.0003225447 
    ## ... Similar to previous best
    ## Run 193 stress 9.920088e-05 
    ## ... Procrustes: rmse 8.098429e-05  max resid 0.0001309949 
    ## ... Similar to previous best
    ## Run 194 stress 9.637978e-05 
    ## ... Procrustes: rmse 0.0001914725  max resid 0.0003814048 
    ## ... Similar to previous best
    ## Run 195 stress 9.752963e-05 
    ## ... Procrustes: rmse 0.0001354981  max resid 0.0002912698 
    ## ... Similar to previous best
    ## Run 196 stress 9.686972e-05 
    ## ... Procrustes: rmse 0.0001781652  max resid 0.0003580439 
    ## ... Similar to previous best
    ## Run 197 stress 9.55846e-05 
    ## ... Procrustes: rmse 6.862316e-05  max resid 9.247714e-05 
    ## ... Similar to previous best
    ## Run 198 stress 9.293136e-05 
    ## ... Procrustes: rmse 0.0001488396  max resid 0.0002892699 
    ## ... Similar to previous best
    ## Run 199 stress 9.856946e-05 
    ## ... Procrustes: rmse 0.0001551372  max resid 0.0002836456 
    ## ... Similar to previous best
    ## Run 200 stress 9.445679e-05 
    ## ... Procrustes: rmse 0.0001573158  max resid 0.0003589247 
    ## ... Similar to previous best
    ## Run 201 stress 9.861572e-05 
    ## ... Procrustes: rmse 0.0001478971  max resid 0.0002384996 
    ## ... Similar to previous best
    ## Run 202 stress 9.39506e-05 
    ## ... Procrustes: rmse 7.096655e-05  max resid 9.267447e-05 
    ## ... Similar to previous best
    ## Run 203 stress 8.196295e-05 
    ## ... Procrustes: rmse 0.00011523  max resid 0.0002730058 
    ## ... Similar to previous best
    ## Run 204 stress 9.339592e-05 
    ## ... Procrustes: rmse 0.0001755781  max resid 0.0002753535 
    ## ... Similar to previous best
    ## Run 205 stress 9.777293e-05 
    ## ... Procrustes: rmse 8.132177e-05  max resid 0.0001119587 
    ## ... Similar to previous best
    ## Run 206 stress 8.934633e-05 
    ## ... Procrustes: rmse 0.0001937505  max resid 0.0003905471 
    ## ... Similar to previous best
    ## Run 207 stress 9.870565e-05 
    ## ... Procrustes: rmse 0.000171238  max resid 0.0002567967 
    ## ... Similar to previous best
    ## Run 208 stress 9.626364e-05 
    ## ... Procrustes: rmse 0.0001685773  max resid 0.0003749868 
    ## ... Similar to previous best
    ## Run 209 stress 9.871289e-05 
    ## ... Procrustes: rmse 0.000144316  max resid 0.0003025778 
    ## ... Similar to previous best
    ## Run 210 stress 9.391986e-05 
    ## ... Procrustes: rmse 0.0001507073  max resid 0.0002917136 
    ## ... Similar to previous best
    ## Run 211 stress 9.955726e-05 
    ## ... Procrustes: rmse 0.0001801323  max resid 0.0002703218 
    ## ... Similar to previous best
    ## Run 212 stress 9.603954e-05 
    ## ... Procrustes: rmse 0.0001632057  max resid 0.0003556289 
    ## ... Similar to previous best
    ## Run 213 stress 9.15923e-05 
    ## ... Procrustes: rmse 0.0001418896  max resid 0.0002655169 
    ## ... Similar to previous best
    ## Run 214 stress 9.632727e-05 
    ## ... Procrustes: rmse 0.0001510561  max resid 0.0002732355 
    ## ... Similar to previous best
    ## Run 215 stress 8.814164e-05 
    ## ... Procrustes: rmse 0.000103759  max resid 0.0002867299 
    ## ... Similar to previous best
    ## Run 216 stress 9.479947e-05 
    ## ... Procrustes: rmse 0.000192941  max resid 0.0003684367 
    ## ... Similar to previous best
    ## Run 217 stress 9.655435e-05 
    ## ... Procrustes: rmse 7.064602e-05  max resid 0.0001126886 
    ## ... Similar to previous best
    ## Run 218 stress 9.477567e-05 
    ## ... Procrustes: rmse 0.0001386159  max resid 0.0003006861 
    ## ... Similar to previous best
    ## Run 219 stress 9.838639e-05 
    ## ... Procrustes: rmse 0.000185019  max resid 0.0003841841 
    ## ... Similar to previous best
    ## Run 220 stress 9.516696e-05 
    ## ... Procrustes: rmse 0.0001849102  max resid 0.0003090399 
    ## ... Similar to previous best
    ## Run 221 stress 9.590315e-05 
    ## ... Procrustes: rmse 0.0001977657  max resid 0.0003934132 
    ## ... Similar to previous best
    ## Run 222 stress 9.904078e-05 
    ## ... Procrustes: rmse 0.0001875165  max resid 0.0003810243 
    ## ... Similar to previous best
    ## Run 223 stress 8.941165e-05 
    ## ... Procrustes: rmse 0.0001756003  max resid 0.0002959516 
    ## ... Similar to previous best
    ## Run 224 stress 9.842223e-05 
    ## ... Procrustes: rmse 0.0001550212  max resid 0.0002939385 
    ## ... Similar to previous best
    ## Run 225 stress 8.270233e-05 
    ## ... Procrustes: rmse 0.0001008933  max resid 0.0002483193 
    ## ... Similar to previous best
    ## Run 226 stress 9.732332e-05 
    ## ... Procrustes: rmse 0.0001031747  max resid 0.0002156665 
    ## ... Similar to previous best
    ## Run 227 stress 9.345201e-05 
    ## ... Procrustes: rmse 9.891767e-05  max resid 0.0002102748 
    ## ... Similar to previous best
    ## Run 228 stress 9.697471e-05 
    ## ... Procrustes: rmse 0.0001365689  max resid 0.0002257832 
    ## ... Similar to previous best
    ## Run 229 stress 9.685344e-05 
    ## ... Procrustes: rmse 0.0001412412  max resid 0.0002998977 
    ## ... Similar to previous best
    ## Run 230 stress 9.829788e-05 
    ## ... Procrustes: rmse 0.0001458722  max resid 0.0002561616 
    ## ... Similar to previous best
    ## Run 231 stress 9.276003e-05 
    ## ... Procrustes: rmse 0.0001536017  max resid 0.0002433422 
    ## ... Similar to previous best
    ## Run 232 stress 9.765774e-05 
    ## ... Procrustes: rmse 0.0001424628  max resid 0.0002988919 
    ## ... Similar to previous best
    ## Run 233 stress 8.935554e-05 
    ## ... Procrustes: rmse 0.0001574233  max resid 0.0003000051 
    ## ... Similar to previous best
    ## Run 234 stress 9.83464e-05 
    ## ... Procrustes: rmse 8.080043e-05  max resid 0.0001053845 
    ## ... Similar to previous best
    ## Run 235 stress 9.095259e-05 
    ## ... Procrustes: rmse 6.312205e-05  max resid 8.97612e-05 
    ## ... Similar to previous best
    ## Run 236 stress 8.945619e-05 
    ## ... Procrustes: rmse 8.732345e-05  max resid 0.0001998285 
    ## ... Similar to previous best
    ## Run 237 stress 9.105949e-05 
    ## ... Procrustes: rmse 0.0001331515  max resid 0.0002527516 
    ## ... Similar to previous best
    ## Run 238 stress 9.778851e-05 
    ## ... Procrustes: rmse 0.000117237  max resid 0.0002952527 
    ## ... Similar to previous best
    ## Run 239 stress 9.191966e-05 
    ## ... Procrustes: rmse 0.0001323199  max resid 0.0002910007 
    ## ... Similar to previous best
    ## Run 240 stress 9.998333e-05 
    ## ... Procrustes: rmse 0.0001542075  max resid 0.0003169982 
    ## ... Similar to previous best
    ## Run 241 stress 9.465878e-05 
    ## ... Procrustes: rmse 0.0001836905  max resid 0.0003739233 
    ## ... Similar to previous best
    ## Run 242 stress 9.50112e-05 
    ## ... Procrustes: rmse 0.0001226295  max resid 0.0002812522 
    ## ... Similar to previous best
    ## Run 243 stress 9.138855e-05 
    ## ... Procrustes: rmse 0.0001400052  max resid 0.0002699949 
    ## ... Similar to previous best
    ## Run 244 stress 9.940016e-05 
    ## ... Procrustes: rmse 0.0002087218  max resid 0.0004109008 
    ## ... Similar to previous best
    ## Run 245 stress 8.888819e-05 
    ## ... Procrustes: rmse 0.0001359739  max resid 0.0002617123 
    ## ... Similar to previous best
    ## Run 246 stress 9.428846e-05 
    ## ... Procrustes: rmse 0.000136163  max resid 0.0002947401 
    ## ... Similar to previous best
    ## Run 247 stress 9.309838e-05 
    ## ... Procrustes: rmse 6.60415e-05  max resid 8.860294e-05 
    ## ... Similar to previous best
    ## Run 248 stress 9.477251e-05 
    ## ... Procrustes: rmse 7.423822e-05  max resid 0.000128177 
    ## ... Similar to previous best
    ## Run 249 stress 9.48133e-05 
    ## ... Procrustes: rmse 6.848122e-05  max resid 9.83431e-05 
    ## ... Similar to previous best
    ## Run 250 stress 9.251948e-05 
    ## ... Procrustes: rmse 0.0001590535  max resid 0.0003169028 
    ## ... Similar to previous best
    ## Run 251 stress 8.946313e-05 
    ## ... Procrustes: rmse 0.0001077939  max resid 0.0001908094 
    ## ... Similar to previous best
    ## Run 252 stress 9.577215e-05 
    ## ... Procrustes: rmse 0.0001355393  max resid 0.0002954322 
    ## ... Similar to previous best
    ## Run 253 stress 8.627419e-05 
    ## ... Procrustes: rmse 0.0001022543  max resid 0.0002330377 
    ## ... Similar to previous best
    ## Run 254 stress 9.513634e-05 
    ## ... Procrustes: rmse 0.0001123432  max resid 0.0002616899 
    ## ... Similar to previous best
    ## Run 255 stress 9.79505e-05 
    ## ... Procrustes: rmse 0.0001401476  max resid 0.0002993819 
    ## ... Similar to previous best
    ## Run 256 stress 9.96967e-05 
    ## ... Procrustes: rmse 0.0002000865  max resid 0.0003495948 
    ## ... Similar to previous best
    ## Run 257 stress 9.304499e-05 
    ## ... Procrustes: rmse 0.0001356708  max resid 0.0002957442 
    ## ... Similar to previous best
    ## Run 258 stress 9.417357e-05 
    ## ... Procrustes: rmse 0.000136233  max resid 0.0002959367 
    ## ... Similar to previous best
    ## Run 259 stress 8.620938e-05 
    ## ... Procrustes: rmse 0.0001877843  max resid 0.0003822104 
    ## ... Similar to previous best
    ## Run 260 stress 9.898721e-05 
    ## ... Procrustes: rmse 0.0001406313  max resid 0.0002972627 
    ## ... Similar to previous best
    ## Run 261 stress 8.454588e-05 
    ## ... Procrustes: rmse 5.116136e-05  max resid 7.181788e-05 
    ## ... Similar to previous best
    ## Run 262 stress 9.580847e-05 
    ## ... Procrustes: rmse 7.460026e-05  max resid 0.000100392 
    ## ... Similar to previous best
    ## Run 263 stress 9.137729e-05 
    ## ... Procrustes: rmse 0.0001399947  max resid 0.0002638508 
    ## ... Similar to previous best
    ## Run 264 stress 9.609024e-05 
    ## ... Procrustes: rmse 0.0001138445  max resid 0.0001900133 
    ## ... Similar to previous best
    ## Run 265 stress 9.703469e-05 
    ## ... Procrustes: rmse 0.0002060772  max resid 0.0003194935 
    ## ... Similar to previous best
    ## Run 266 stress 9.745833e-05 
    ## ... Procrustes: rmse 0.0001755991  max resid 0.0002640075 
    ## ... Similar to previous best
    ## Run 267 stress 8.962359e-05 
    ## ... Procrustes: rmse 5.65328e-05  max resid 8.138884e-05 
    ## ... Similar to previous best
    ## Run 268 stress 9.640553e-05 
    ## ... Procrustes: rmse 0.000157662  max resid 0.000342726 
    ## ... Similar to previous best
    ## Run 269 stress 8.661347e-05 
    ## ... Procrustes: rmse 8.341739e-05  max resid 0.0001845162 
    ## ... Similar to previous best
    ## Run 270 stress 9.7866e-05 
    ## ... Procrustes: rmse 0.0001441146  max resid 0.0003023087 
    ## ... Similar to previous best
    ## Run 271 stress 9.579202e-05 
    ## ... Procrustes: rmse 0.0001021231  max resid 0.0002140646 
    ## ... Similar to previous best
    ## Run 272 stress 9.892295e-05 
    ## ... Procrustes: rmse 0.0001452555  max resid 0.0003027372 
    ## ... Similar to previous best
    ## Run 273 stress 9.982143e-05 
    ## ... Procrustes: rmse 0.0001579119  max resid 0.000304652 
    ## ... Similar to previous best
    ## Run 274 stress 9.379528e-05 
    ## ... Procrustes: rmse 0.000202093  max resid 0.0004008064 
    ## ... Similar to previous best
    ## Run 275 stress 9.965194e-05 
    ## ... Procrustes: rmse 0.0001515979  max resid 0.0002388209 
    ## ... Similar to previous best
    ## Run 276 stress 9.400055e-05 
    ## ... Procrustes: rmse 0.0001356765  max resid 0.0002927326 
    ## ... Similar to previous best
    ## Run 277 stress 9.825092e-05 
    ## ... Procrustes: rmse 0.0001768131  max resid 0.0002646433 
    ## ... Similar to previous best
    ## Run 278 stress 9.504475e-05 
    ## ... Procrustes: rmse 0.0001569486  max resid 0.0003523487 
    ## ... Similar to previous best
    ## Run 279 stress 9.716584e-05 
    ## ... Procrustes: rmse 0.0001452356  max resid 0.0002355032 
    ## ... Similar to previous best
    ## Run 280 stress 8.763995e-05 
    ## ... Procrustes: rmse 8.580857e-05  max resid 0.0001895345 
    ## ... Similar to previous best
    ## Run 281 stress 9.548568e-05 
    ## ... Procrustes: rmse 0.0001709737  max resid 0.000356374 
    ## ... Similar to previous best
    ## Run 282 stress 9.525212e-05 
    ## ... Procrustes: rmse 0.0001272796  max resid 0.0002850402 
    ## ... Similar to previous best
    ## Run 283 stress 9.752187e-05 
    ## ... Procrustes: rmse 0.0001998167  max resid 0.0003228312 
    ## ... Similar to previous best
    ## Run 284 stress 9.738438e-05 
    ## ... Procrustes: rmse 0.0001544898  max resid 0.0002972757 
    ## ... Similar to previous best
    ## Run 285 stress 9.522077e-05 
    ## ... Procrustes: rmse 7.367335e-05  max resid 0.0001015294 
    ## ... Similar to previous best
    ## Run 286 stress 9.402087e-05 
    ## ... Procrustes: rmse 9.772472e-05  max resid 0.0002074596 
    ## ... Similar to previous best
    ## Run 287 stress 9.456392e-05 
    ## ... Procrustes: rmse 0.0001881188  max resid 0.0003784054 
    ## ... Similar to previous best
    ## Run 288 stress 9.903929e-05 
    ## ... Procrustes: rmse 0.0002005062  max resid 0.0003768124 
    ## ... Similar to previous best
    ## Run 289 stress 9.96968e-05 
    ## ... Procrustes: rmse 0.00020234  max resid 0.0003874782 
    ## ... Similar to previous best
    ## Run 290 stress 8.818229e-05 
    ## ... Procrustes: rmse 0.0001436688  max resid 0.0002801864 
    ## ... Similar to previous best
    ## Run 291 stress 9.637996e-05 
    ## ... Procrustes: rmse 0.0001410619  max resid 0.0003013631 
    ## ... Similar to previous best
    ## Run 292 stress 9.218383e-05 
    ## ... Procrustes: rmse 0.0001458992  max resid 0.0002610364 
    ## ... Similar to previous best
    ## Run 293 stress 9.033917e-05 
    ## ... Procrustes: rmse 0.0001236126  max resid 0.0002763882 
    ## ... Similar to previous best
    ## Run 294 stress 9.566045e-05 
    ## ... Procrustes: rmse 0.0001709692  max resid 0.0003039848 
    ## ... Similar to previous best
    ## Run 295 stress 9.873756e-05 
    ## ... Procrustes: rmse 0.0001926096  max resid 0.0003056096 
    ## ... Similar to previous best
    ## Run 296 stress 8.182561e-05 
    ## ... Procrustes: rmse 0.0001168857  max resid 0.0002650777 
    ## ... Similar to previous best
    ## Run 297 stress 9.861514e-05 
    ## ... Procrustes: rmse 0.0001692646  max resid 0.0003548221 
    ## ... Similar to previous best
    ## Run 298 stress 9.85864e-05 
    ## ... Procrustes: rmse 0.0001701938  max resid 0.0003510795 
    ## ... Similar to previous best
    ## Run 299 stress 9.025307e-05 
    ## ... Procrustes: rmse 0.00012726  max resid 0.0002852421 
    ## ... Similar to previous best
    ## Run 300 stress 9.404845e-05 
    ## ... Procrustes: rmse 0.0001661786  max resid 0.0002485019 
    ## ... Similar to previous best
    ## Run 301 stress 9.63471e-05 
    ## ... Procrustes: rmse 0.0001729077  max resid 0.0002601143 
    ## ... Similar to previous best
    ## Run 302 stress 9.287069e-05 
    ## ... Procrustes: rmse 0.000163179  max resid 0.0003035361 
    ## ... Similar to previous best
    ## Run 303 stress 9.715676e-05 
    ## ... Procrustes: rmse 0.0001575409  max resid 0.0002373725 
    ## ... Similar to previous best
    ## Run 304 stress 9.672818e-05 
    ## ... Procrustes: rmse 0.0002078831  max resid 0.0004085029 
    ## ... Similar to previous best
    ## Run 305 stress 9.352698e-05 
    ## ... Procrustes: rmse 0.0001776475  max resid 0.0003531746 
    ## ... Similar to previous best
    ## Run 306 stress 9.52076e-05 
    ## ... Procrustes: rmse 0.0001075484  max resid 0.0002756148 
    ## ... Similar to previous best
    ## Run 307 stress 9.094355e-05 
    ## ... Procrustes: rmse 0.0001734307  max resid 0.000265261 
    ## ... Similar to previous best
    ## Run 308 stress 9.012125e-05 
    ## ... Procrustes: rmse 0.0001932565  max resid 0.0003890288 
    ## ... Similar to previous best
    ## Run 309 stress 0.2282325 
    ## Run 310 stress 9.801914e-05 
    ## ... Procrustes: rmse 0.0002108354  max resid 0.0004132595 
    ## ... Similar to previous best
    ## Run 311 stress 8.90232e-05 
    ## ... Procrustes: rmse 0.0001275676  max resid 0.0002831202 
    ## ... Similar to previous best
    ## Run 312 stress 9.441764e-05 
    ## ... Procrustes: rmse 0.000137365  max resid 0.0002949141 
    ## ... Similar to previous best
    ## Run 313 stress 8.564259e-05 
    ## ... Procrustes: rmse 0.0001448032  max resid 0.0002929993 
    ## ... Similar to previous best
    ## Run 314 stress 9.441106e-05 
    ## ... Procrustes: rmse 0.0001920088  max resid 0.0003388825 
    ## ... Similar to previous best
    ## Run 315 stress 9.679027e-05 
    ## ... Procrustes: rmse 0.0001579226  max resid 0.0003184827 
    ## ... Similar to previous best
    ## Run 316 stress 9.262837e-05 
    ## ... Procrustes: rmse 6.510956e-05  max resid 9.95656e-05 
    ## ... Similar to previous best
    ## Run 317 stress 9.871372e-05 
    ## ... Procrustes: rmse 0.0001455311  max resid 0.000304547 
    ## ... Similar to previous best
    ## Run 318 stress 8.296584e-05 
    ## ... Procrustes: rmse 0.0001436328  max resid 0.0003029868 
    ## ... Similar to previous best
    ## Run 319 stress 8.216551e-05 
    ## ... Procrustes: rmse 4.567762e-05  max resid 6.496826e-05 
    ## ... Similar to previous best
    ## Run 320 stress 9.012084e-05 
    ## ... Procrustes: rmse 0.0001879474  max resid 0.0003006693 
    ## ... Similar to previous best
    ## Run 321 stress 9.30499e-05 
    ## ... Procrustes: rmse 0.0001885856  max resid 0.000365489 
    ## ... Similar to previous best
    ## Run 322 stress 9.483467e-05 
    ## ... Procrustes: rmse 0.0001252291  max resid 0.0001901015 
    ## ... Similar to previous best
    ## Run 323 stress 9.929505e-05 
    ## ... Procrustes: rmse 7.464933e-05  max resid 0.0001087527 
    ## ... Similar to previous best
    ## Run 324 stress 9.837719e-05 
    ## ... Procrustes: rmse 0.0001720915  max resid 0.0003967359 
    ## ... Similar to previous best
    ## Run 325 stress 9.672789e-05 
    ## ... Procrustes: rmse 0.0001412459  max resid 0.0003006412 
    ## ... Similar to previous best
    ## Run 326 stress 9.592411e-05 
    ## ... Procrustes: rmse 0.0001649829  max resid 0.0003691169 
    ## ... Similar to previous best
    ## Run 327 stress 9.894734e-05 
    ## ... Procrustes: rmse 7.739282e-05  max resid 0.0001054471 
    ## ... Similar to previous best
    ## Run 328 stress 9.937513e-05 
    ## ... Procrustes: rmse 0.0002083608  max resid 0.0003221985 
    ## ... Similar to previous best
    ## Run 329 stress 9.562255e-05 
    ## ... Procrustes: rmse 0.0001376836  max resid 0.0002957725 
    ## ... Similar to previous best
    ## Run 330 stress 9.404193e-05 
    ## ... Procrustes: rmse 0.000191752  max resid 0.0003391477 
    ## ... Similar to previous best
    ## Run 331 stress 9.156304e-05 
    ## ... Procrustes: rmse 5.712601e-05  max resid 0.0001050933 
    ## ... Similar to previous best
    ## Run 332 stress 9.622069e-05 
    ## ... Procrustes: rmse 0.0001141745  max resid 0.0002608107 
    ## ... Similar to previous best
    ## Run 333 stress 9.560696e-05 
    ## ... Procrustes: rmse 0.0001570571  max resid 0.0003598425 
    ## ... Similar to previous best
    ## Run 334 stress 8.760445e-05 
    ## ... Procrustes: rmse 0.0001281791  max resid 0.0002720193 
    ## ... Similar to previous best
    ## Run 335 stress 9.882872e-05 
    ## ... Procrustes: rmse 8.374307e-05  max resid 0.0001183574 
    ## ... Similar to previous best
    ## Run 336 stress 9.318148e-05 
    ## ... Procrustes: rmse 0.0001925563  max resid 0.0003418602 
    ## ... Similar to previous best
    ## Run 337 stress 9.350053e-05 
    ## ... Procrustes: rmse 0.0001054546  max resid 0.0001658059 
    ## ... Similar to previous best
    ## Run 338 stress 9.049952e-05 
    ## ... Procrustes: rmse 0.0001313243  max resid 0.0002888353 
    ## ... Similar to previous best
    ## Run 339 stress 9.928154e-05 
    ## ... Procrustes: rmse 0.0001458577  max resid 0.0003054771 
    ## ... Similar to previous best
    ## Run 340 stress 9.148732e-05 
    ## ... Procrustes: rmse 6.653694e-05  max resid 9.238686e-05 
    ## ... Similar to previous best
    ## Run 341 stress 9.699956e-05 
    ## ... Procrustes: rmse 0.0002046997  max resid 0.0003177989 
    ## ... Similar to previous best
    ## Run 342 stress 8.929282e-05 
    ## ... Procrustes: rmse 0.0001245387  max resid 0.0003088862 
    ## ... Similar to previous best
    ## Run 343 stress 8.916846e-05 
    ## ... Procrustes: rmse 0.0001297357  max resid 0.0002883139 
    ## ... Similar to previous best
    ## Run 344 stress 9.535098e-05 
    ## ... Procrustes: rmse 0.0001990981  max resid 0.0003128714 
    ## ... Similar to previous best
    ## Run 345 stress 9.260622e-05 
    ## ... Procrustes: rmse 0.0001492453  max resid 0.0002739497 
    ## ... Similar to previous best
    ## Run 346 stress 9.188566e-05 
    ## ... Procrustes: rmse 0.0001394191  max resid 0.0002387536 
    ## ... Similar to previous best
    ## Run 347 stress 9.499092e-05 
    ## ... Procrustes: rmse 0.0001881609  max resid 0.0003582733 
    ## ... Similar to previous best
    ## Run 348 stress 9.302268e-05 
    ## ... Procrustes: rmse 0.0001387638  max resid 0.0002628655 
    ## ... Similar to previous best
    ## Run 349 stress 9.910634e-05 
    ## ... Procrustes: rmse 0.0002028615  max resid 0.0003259562 
    ## ... Similar to previous best
    ## Run 350 stress 7.787346e-05 
    ## ... Procrustes: rmse 0.0001596275  max resid 0.0002824627 
    ## ... Similar to previous best
    ## Run 351 stress 7.83697e-05 
    ## ... Procrustes: rmse 0.0001620532  max resid 0.0002833154 
    ## ... Similar to previous best
    ## Run 352 stress 9.929643e-05 
    ## ... Procrustes: rmse 8.214277e-05  max resid 0.0001146469 
    ## ... Similar to previous best
    ## Run 353 stress 9.996398e-05 
    ## ... Procrustes: rmse 8.620462e-05  max resid 0.0001189492 
    ## ... Similar to previous best
    ## Run 354 stress 9.683649e-05 
    ## ... Procrustes: rmse 0.0001673921  max resid 0.0003762411 
    ## ... Similar to previous best
    ## Run 355 stress 8.955136e-05 
    ## ... Procrustes: rmse 0.0001794087  max resid 0.000337818 
    ## ... Similar to previous best
    ## Run 356 stress 9.47753e-05 
    ## ... Procrustes: rmse 0.0001308937  max resid 0.0002410696 
    ## ... Similar to previous best
    ## Run 357 stress 8.574048e-05 
    ## ... Procrustes: rmse 0.0001504579  max resid 0.0003234439 
    ## ... Similar to previous best
    ## Run 358 stress 8.775626e-05 
    ## ... Procrustes: rmse 0.0001303726  max resid 0.0002549399 
    ## ... Similar to previous best
    ## Run 359 stress 7.220569e-05 
    ## ... Procrustes: rmse 5.6111e-05  max resid 8.456206e-05 
    ## ... Similar to previous best
    ## Run 360 stress 9.371596e-05 
    ## ... Procrustes: rmse 0.000120688  max resid 0.000230811 
    ## ... Similar to previous best
    ## Run 361 stress 9.506567e-05 
    ## ... Procrustes: rmse 0.0001677254  max resid 0.0002512463 
    ## ... Similar to previous best
    ## Run 362 stress 9.453657e-05 
    ## ... Procrustes: rmse 0.0001386014  max resid 0.0002969476 
    ## ... Similar to previous best
    ## Run 363 stress 9.314527e-05 
    ## ... Procrustes: rmse 0.0001982278  max resid 0.0003125099 
    ## ... Similar to previous best
    ## Run 364 stress 9.398647e-05 
    ## ... Procrustes: rmse 0.0001789612  max resid 0.0002690573 
    ## ... Similar to previous best
    ## Run 365 stress 9.9458e-05 
    ## ... Procrustes: rmse 0.0001478485  max resid 0.0002341998 
    ## ... Similar to previous best
    ## Run 366 stress 9.647911e-05 
    ## ... Procrustes: rmse 8.077934e-05  max resid 0.0001716717 
    ## ... Similar to previous best
    ## Run 367 stress 9.938315e-05 
    ## ... Procrustes: rmse 0.0002026872  max resid 0.0003254978 
    ## ... Similar to previous best
    ## Run 368 stress 8.221016e-05 
    ## ... Procrustes: rmse 0.000176996  max resid 0.0003653965 
    ## ... Similar to previous best
    ## Run 369 stress 9.243418e-05 
    ## ... Procrustes: rmse 0.0001705455  max resid 0.0003153554 
    ## ... Similar to previous best
    ## Run 370 stress 9.924814e-05 
    ## ... Procrustes: rmse 0.0001739987  max resid 0.0003717052 
    ## ... Similar to previous best
    ## Run 371 stress 9.622968e-05 
    ## ... Procrustes: rmse 0.0001816077  max resid 0.0003583226 
    ## ... Similar to previous best
    ## Run 372 stress 6.561038e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 6.352203e-05  max resid 0.0001391777 
    ## ... Similar to previous best
    ## Run 373 stress 9.471356e-05 
    ## ... Procrustes: rmse 7.449357e-05  max resid 0.0001087583 
    ## ... Similar to previous best
    ## Run 374 stress 9.560509e-05 
    ## ... Procrustes: rmse 0.0001972753  max resid 0.00039968 
    ## ... Similar to previous best
    ## Run 375 stress 9.722389e-05 
    ## ... Procrustes: rmse 0.0001180952  max resid 0.0001759816 
    ## ... Similar to previous best
    ## Run 376 stress 9.595462e-05 
    ## ... Procrustes: rmse 0.0002111554  max resid 0.0003967439 
    ## ... Similar to previous best
    ## Run 377 stress 9.917869e-05 
    ## ... Procrustes: rmse 0.0001363534  max resid 0.0002877881 
    ## ... Similar to previous best
    ## Run 378 stress 9.905557e-05 
    ## ... Procrustes: rmse 0.0002259172  max resid 0.0003413876 
    ## ... Similar to previous best
    ## Run 379 stress 9.078445e-05 
    ## ... Procrustes: rmse 0.0001372725  max resid 0.0002838152 
    ## ... Similar to previous best
    ## Run 380 stress 9.921788e-05 
    ## ... Procrustes: rmse 0.0001137522  max resid 0.0002121437 
    ## ... Similar to previous best
    ## Run 381 stress 9.355579e-05 
    ## ... Procrustes: rmse 0.0001602667  max resid 0.0002887971 
    ## ... Similar to previous best
    ## Run 382 stress 7.986623e-05 
    ## ... Procrustes: rmse 4.105089e-05  max resid 7.930645e-05 
    ## ... Similar to previous best
    ## Run 383 stress 9.80003e-05 
    ## ... Procrustes: rmse 0.0001155335  max resid 0.0002427091 
    ## ... Similar to previous best
    ## Run 384 stress 9.313406e-05 
    ## ... Procrustes: rmse 0.0001747222  max resid 0.0003877697 
    ## ... Similar to previous best
    ## Run 385 stress 9.43507e-05 
    ## ... Procrustes: rmse 0.0001580689  max resid 0.0002871679 
    ## ... Similar to previous best
    ## Run 386 stress 9.593968e-05 
    ## ... Procrustes: rmse 0.0002046333  max resid 0.0003854424 
    ## ... Similar to previous best
    ## Run 387 stress 9.635123e-05 
    ## ... Procrustes: rmse 0.0001205378  max resid 0.0001928718 
    ## ... Similar to previous best
    ## Run 388 stress 9.690144e-05 
    ## ... Procrustes: rmse 0.0002110027  max resid 0.0003993633 
    ## ... Similar to previous best
    ## Run 389 stress 9.288633e-05 
    ## ... Procrustes: rmse 0.000221265  max resid 0.0004147836 
    ## ... Similar to previous best
    ## Run 390 stress 9.172886e-05 
    ## ... Procrustes: rmse 0.0001084546  max resid 0.0001848905 
    ## ... Similar to previous best
    ## Run 391 stress 9.499041e-05 
    ## ... Procrustes: rmse 0.0002056438  max resid 0.0003883894 
    ## ... Similar to previous best
    ## Run 392 stress 9.343054e-05 
    ## ... Procrustes: rmse 0.0001392225  max resid 0.0002782127 
    ## ... Similar to previous best
    ## Run 393 stress 9.838275e-05 
    ## ... Procrustes: rmse 0.0001862034  max resid 0.000341422 
    ## ... Similar to previous best
    ## Run 394 stress 9.737474e-05 
    ## ... Procrustes: rmse 0.000199785  max resid 0.0003713679 
    ## ... Similar to previous best
    ## Run 395 stress 9.938989e-05 
    ## ... Procrustes: rmse 0.0001856044  max resid 0.0003891748 
    ## ... Similar to previous best
    ## Run 396 stress 8.520192e-05 
    ## ... Procrustes: rmse 6.116517e-05  max resid 0.0001277547 
    ## ... Similar to previous best
    ## Run 397 stress 9.812245e-05 
    ## ... Procrustes: rmse 0.0001858229  max resid 0.000399558 
    ## ... Similar to previous best
    ## Run 398 stress 7.725428e-05 
    ## ... Procrustes: rmse 0.0001804819  max resid 0.0002701155 
    ## ... Similar to previous best
    ## Run 399 stress 9.032809e-05 
    ## ... Procrustes: rmse 0.0001923637  max resid 0.0003719791 
    ## ... Similar to previous best
    ## Run 400 stress 9.758707e-05 
    ## ... Procrustes: rmse 0.0001814138  max resid 0.0003402034 
    ## ... Similar to previous best
    ## Run 401 stress 9.949708e-05 
    ## ... Procrustes: rmse 6.756807e-05  max resid 0.0001095551 
    ## ... Similar to previous best
    ## Run 402 stress 9.828918e-05 
    ## ... Procrustes: rmse 0.0001844559  max resid 0.0002580063 
    ## ... Similar to previous best
    ## Run 403 stress 6.889061e-05 
    ## ... Procrustes: rmse 0.0001398356  max resid 0.0002762832 
    ## ... Similar to previous best
    ## Run 404 stress 9.516469e-05 
    ## ... Procrustes: rmse 0.0001934515  max resid 0.0003657813 
    ## ... Similar to previous best
    ## Run 405 stress 8.299394e-05 
    ## ... Procrustes: rmse 0.0001343898  max resid 0.0002657073 
    ## ... Similar to previous best
    ## Run 406 stress 9.611461e-05 
    ## ... Procrustes: rmse 0.0002115048  max resid 0.0003975147 
    ## ... Similar to previous best
    ## Run 407 stress 9.921843e-05 
    ## ... Procrustes: rmse 0.0002053522  max resid 0.0003210564 
    ## ... Similar to previous best
    ## Run 408 stress 9.82516e-05 
    ## ... Procrustes: rmse 0.0001801106  max resid 0.000327265 
    ## ... Similar to previous best
    ## Run 409 stress 9.874438e-05 
    ## ... Procrustes: rmse 0.0001657954  max resid 0.0003561579 
    ## ... Similar to previous best
    ## Run 410 stress 9.869618e-05 
    ## ... Procrustes: rmse 0.000225095  max resid 0.0003397933 
    ## ... Similar to previous best
    ## Run 411 stress 9.836297e-05 
    ## ... Procrustes: rmse 0.0002030735  max resid 0.0003078024 
    ## ... Similar to previous best
    ## Run 412 stress 9.386242e-05 
    ## ... Procrustes: rmse 0.000110378  max resid 0.0001782637 
    ## ... Similar to previous best
    ## Run 413 stress 9.50984e-05 
    ## ... Procrustes: rmse 0.0001895696  max resid 0.0003591739 
    ## ... Similar to previous best
    ## Run 414 stress 9.676563e-05 
    ## ... Procrustes: rmse 0.0001948799  max resid 0.0003685379 
    ## ... Similar to previous best
    ## Run 415 stress 9.436037e-05 
    ## ... Procrustes: rmse 0.0001854675  max resid 0.0003186735 
    ## ... Similar to previous best
    ## Run 416 stress 9.917869e-05 
    ## ... Procrustes: rmse 0.000181718  max resid 0.0003837475 
    ## ... Similar to previous best
    ## Run 417 stress 9.442462e-05 
    ## ... Procrustes: rmse 0.0001628813  max resid 0.0002920704 
    ## ... Similar to previous best
    ## Run 418 stress 9.747711e-05 
    ## ... Procrustes: rmse 0.0001470887  max resid 0.0002805006 
    ## ... Similar to previous best
    ## Run 419 stress 9.613739e-05 
    ## ... Procrustes: rmse 0.0002027867  max resid 0.0003897512 
    ## ... Similar to previous best
    ## Run 420 stress 8.750831e-05 
    ## ... Procrustes: rmse 0.0001939043  max resid 0.0002861818 
    ## ... Similar to previous best
    ## Run 421 stress 9.94647e-05 
    ## ... Procrustes: rmse 0.0001895605  max resid 0.0003266835 
    ## ... Similar to previous best
    ## Run 422 stress 9.880082e-05 
    ## ... Procrustes: rmse 0.000170571  max resid 0.0003000335 
    ## ... Similar to previous best
    ## Run 423 stress 9.754518e-05 
    ## ... Procrustes: rmse 0.000139111  max resid 0.0002748588 
    ## ... Similar to previous best
    ## Run 424 stress 8.75688e-05 
    ## ... Procrustes: rmse 0.0001007697  max resid 0.0001774477 
    ## ... Similar to previous best
    ## Run 425 stress 9.308166e-05 
    ## ... Procrustes: rmse 0.000109562  max resid 0.0001797627 
    ## ... Similar to previous best
    ## Run 426 stress 9.947689e-05 
    ## ... Procrustes: rmse 0.0001517623  max resid 0.0002805314 
    ## ... Similar to previous best
    ## Run 427 stress 8.431598e-05 
    ## ... Procrustes: rmse 0.0001930104  max resid 0.0003050268 
    ## ... Similar to previous best
    ## Run 428 stress 9.656529e-05 
    ## ... Procrustes: rmse 0.0002225183  max resid 0.0003214787 
    ## ... Similar to previous best
    ## Run 429 stress 7.092357e-05 
    ## ... Procrustes: rmse 0.0001019864  max resid 0.0001794509 
    ## ... Similar to previous best
    ## Run 430 stress 9.868185e-05 
    ## ... Procrustes: rmse 0.0001701847  max resid 0.0002994088 
    ## ... Similar to previous best
    ## Run 431 stress 9.934527e-05 
    ## ... Procrustes: rmse 0.0001954509  max resid 0.0003351381 
    ## ... Similar to previous best
    ## Run 432 stress 9.443957e-05 
    ## ... Procrustes: rmse 0.00011433  max resid 0.000184643 
    ## ... Similar to previous best
    ## Run 433 stress 9.611173e-05 
    ## ... Procrustes: rmse 0.0001201445  max resid 0.0002029081 
    ## ... Similar to previous best
    ## Run 434 stress 9.439926e-05 
    ## ... Procrustes: rmse 0.0001627494  max resid 0.0002907286 
    ## ... Similar to previous best
    ## Run 435 stress 9.031407e-05 
    ## ... Procrustes: rmse 0.00014657  max resid 0.0002634756 
    ## ... Similar to previous best
    ## Run 436 stress 9.698908e-05 
    ## ... Procrustes: rmse 0.0001873185  max resid 0.0003646257 
    ## ... Similar to previous best
    ## Run 437 stress 9.911852e-05 
    ## ... Procrustes: rmse 0.0001490841  max resid 0.0002899485 
    ## ... Similar to previous best
    ## Run 438 stress 9.242832e-05 
    ## ... Procrustes: rmse 0.0001970953  max resid 0.0003790837 
    ## ... Similar to previous best
    ## Run 439 stress 8.989568e-05 
    ## ... Procrustes: rmse 0.0002076592  max resid 0.000298905 
    ## ... Similar to previous best
    ## Run 440 stress 9.926622e-05 
    ## ... Procrustes: rmse 6.600796e-05  max resid 0.0001105804 
    ## ... Similar to previous best
    ## Run 441 stress 9.336128e-05 
    ## ... Procrustes: rmse 0.0001695754  max resid 0.0003611496 
    ## ... Similar to previous best
    ## Run 442 stress 9.424062e-05 
    ## ... Procrustes: rmse 0.0001971451  max resid 0.0003817871 
    ## ... Similar to previous best
    ## Run 443 stress 8.88764e-05 
    ## ... Procrustes: rmse 9.188265e-05  max resid 0.0001497373 
    ## ... Similar to previous best
    ## Run 444 stress 9.757906e-05 
    ## ... Procrustes: rmse 0.0001750706  max resid 0.000369368 
    ## ... Similar to previous best
    ## Run 445 stress 9.225652e-05 
    ## ... Procrustes: rmse 0.0002013986  max resid 0.0003831182 
    ## ... Similar to previous best
    ## Run 446 stress 9.613314e-05 
    ## ... Procrustes: rmse 0.0001242502  max resid 0.0001786357 
    ## ... Similar to previous best
    ## Run 447 stress 9.748927e-05 
    ## ... Procrustes: rmse 0.0001846835  max resid 0.0003188758 
    ## ... Similar to previous best
    ## Run 448 stress 8.951046e-05 
    ## ... Procrustes: rmse 0.0001027985  max resid 0.0002705143 
    ## ... Similar to previous best
    ## Run 449 stress 9.441609e-05 
    ## ... Procrustes: rmse 0.0001337713  max resid 0.0002707538 
    ## ... Similar to previous best
    ## Run 450 stress 8.137674e-05 
    ## ... Procrustes: rmse 0.000129545  max resid 0.0002435367 
    ## ... Similar to previous best
    ## Run 451 stress 8.992915e-05 
    ## ... Procrustes: rmse 9.787777e-05  max resid 0.0002373367 
    ## ... Similar to previous best
    ## Run 452 stress 8.588007e-05 
    ## ... Procrustes: rmse 9.654009e-05  max resid 0.000170338 
    ## ... Similar to previous best
    ## Run 453 stress 9.685781e-05 
    ## ... Procrustes: rmse 9.392507e-05  max resid 0.0002096179 
    ## ... Similar to previous best
    ## Run 454 stress 9.712463e-05 
    ## ... Procrustes: rmse 0.00020411  max resid 0.0003913964 
    ## ... Similar to previous best
    ## Run 455 stress 9.174485e-05 
    ## ... Procrustes: rmse 0.0001565784  max resid 0.0002850973 
    ## ... Similar to previous best
    ## Run 456 stress 9.314811e-05 
    ## ... Procrustes: rmse 0.0001685358  max resid 0.0003412245 
    ## ... Similar to previous best
    ## Run 457 stress 9.00462e-05 
    ## ... Procrustes: rmse 0.0002029546  max resid 0.0002952228 
    ## ... Similar to previous best
    ## Run 458 stress 9.194765e-05 
    ## ... Procrustes: rmse 0.0001058708  max resid 0.0001793804 
    ## ... Similar to previous best
    ## Run 459 stress 9.627058e-05 
    ## ... Procrustes: rmse 0.000167438  max resid 0.0002911377 
    ## ... Similar to previous best
    ## Run 460 stress 9.948748e-05 
    ## ... Procrustes: rmse 0.0001231351  max resid 0.000196761 
    ## ... Similar to previous best
    ## Run 461 stress 8.921014e-05 
    ## ... Procrustes: rmse 0.0001276525  max resid 0.0002339314 
    ## ... Similar to previous best
    ## Run 462 stress 9.474523e-05 
    ## ... Procrustes: rmse 0.0001763718  max resid 0.0003065866 
    ## ... Similar to previous best
    ## Run 463 stress 9.991494e-05 
    ## ... Procrustes: rmse 0.0002255615  max resid 0.0003407702 
    ## ... Similar to previous best
    ## Run 464 stress 9.738056e-05 
    ## ... Procrustes: rmse 0.0002085876  max resid 0.0003952474 
    ## ... Similar to previous best
    ## Run 465 stress 8.803905e-05 
    ## ... Procrustes: rmse 0.0001490833  max resid 0.000276275 
    ## ... Similar to previous best
    ## Run 466 stress 9.959516e-05 
    ## ... Procrustes: rmse 0.0001781699  max resid 0.0002804791 
    ## ... Similar to previous best
    ## Run 467 stress 7.476346e-05 
    ## ... Procrustes: rmse 0.0001097169  max resid 0.0002050042 
    ## ... Similar to previous best
    ## Run 468 stress 9.907443e-05 
    ## ... Procrustes: rmse 0.0002170264  max resid 0.0004044848 
    ## ... Similar to previous best
    ## Run 469 stress 9.363489e-05 
    ## ... Procrustes: rmse 7.550204e-05  max resid 0.0001301429 
    ## ... Similar to previous best
    ## Run 470 stress 8.879926e-05 
    ## ... Procrustes: rmse 0.0001980321  max resid 0.0003922194 
    ## ... Similar to previous best
    ## Run 471 stress 9.953193e-05 
    ## ... Procrustes: rmse 0.000194507  max resid 0.0003699738 
    ## ... Similar to previous best
    ## Run 472 stress 9.605281e-05 
    ## ... Procrustes: rmse 0.0001818628  max resid 0.0003185819 
    ## ... Similar to previous best
    ## Run 473 stress 9.683724e-05 
    ## ... Procrustes: rmse 0.000145138  max resid 0.0002697414 
    ## ... Similar to previous best
    ## Run 474 stress 9.690193e-05 
    ## ... Procrustes: rmse 0.0002122108  max resid 0.0004006006 
    ## ... Similar to previous best
    ## Run 475 stress 9.904533e-05 
    ## ... Procrustes: rmse 0.0002198588  max resid 0.0004214014 
    ## ... Similar to previous best
    ## Run 476 stress 9.937698e-05 
    ## ... Procrustes: rmse 0.0002124451  max resid 0.0004039955 
    ## ... Similar to previous best
    ## Run 477 stress 9.958501e-05 
    ## ... Procrustes: rmse 0.0001890963  max resid 0.000325053 
    ## ... Similar to previous best
    ## Run 478 stress 9.748765e-05 
    ## ... Procrustes: rmse 0.0002135504  max resid 0.0004119483 
    ## ... Similar to previous best
    ## Run 479 stress 9.673728e-05 
    ## ... Procrustes: rmse 0.0002082382  max resid 0.0003958853 
    ## ... Similar to previous best
    ## Run 480 stress 9.790713e-05 
    ## ... Procrustes: rmse 0.0001197237  max resid 0.0001882409 
    ## ... Similar to previous best
    ## Run 481 stress 9.636345e-05 
    ## ... Procrustes: rmse 0.0002111467  max resid 0.0004106396 
    ## ... Similar to previous best
    ## Run 482 stress 9.410962e-05 
    ## ... Procrustes: rmse 0.00012595  max resid 0.0002717538 
    ## ... Similar to previous best
    ## Run 483 stress 9.990486e-05 
    ## ... Procrustes: rmse 0.000172746  max resid 0.0002819683 
    ## ... Similar to previous best
    ## Run 484 stress 9.903776e-05 
    ## ... Procrustes: rmse 0.0001709168  max resid 0.0002992216 
    ## ... Similar to previous best
    ## Run 485 stress 9.749677e-05 
    ## ... Procrustes: rmse 0.0002066682  max resid 0.0004047102 
    ## ... Similar to previous best
    ## Run 486 stress 9.831573e-05 
    ## ... Procrustes: rmse 0.0001858441  max resid 0.0003408688 
    ## ... Similar to previous best
    ## Run 487 stress 9.963467e-05 
    ## ... Procrustes: rmse 0.0002016019  max resid 0.0003771288 
    ## ... Similar to previous best
    ## Run 488 stress 9.956652e-05 
    ## ... Procrustes: rmse 0.0002109814  max resid 0.000412003 
    ## ... Similar to previous best
    ## Run 489 stress 9.037462e-05 
    ## ... Procrustes: rmse 0.000163975  max resid 0.0003171146 
    ## ... Similar to previous best
    ## Run 490 stress 9.098704e-05 
    ## ... Procrustes: rmse 0.0001401724  max resid 0.0002872912 
    ## ... Similar to previous best
    ## Run 491 stress 9.913318e-05 
    ## ... Procrustes: rmse 0.0002084561  max resid 0.0003242696 
    ## ... Similar to previous best
    ## Run 492 stress 9.77772e-05 
    ## ... Procrustes: rmse 0.0002142124  max resid 0.0004016329 
    ## ... Similar to previous best
    ## Run 493 stress 9.401196e-05 
    ## ... Procrustes: rmse 0.0001364585  max resid 0.0002865086 
    ## ... Similar to previous best
    ## Run 494 stress 9.531611e-05 
    ## ... Procrustes: rmse 0.0001605353  max resid 0.0002898944 
    ## ... Similar to previous best
    ## Run 495 stress 8.872394e-05 
    ## ... Procrustes: rmse 0.0001771625  max resid 0.0003342324 
    ## ... Similar to previous best
    ## Run 496 stress 9.609065e-05 
    ## ... Procrustes: rmse 0.0001042516  max resid 0.0001696599 
    ## ... Similar to previous best
    ## Run 497 stress 9.984666e-05 
    ## ... Procrustes: rmse 0.0002104797  max resid 0.000404655 
    ## ... Similar to previous best
    ## Run 498 stress 9.847417e-05 
    ## ... Procrustes: rmse 0.0001495969  max resid 0.0002832329 
    ## ... Similar to previous best
    ## Run 499 stress 9.975475e-05 
    ## ... Procrustes: rmse 0.0002182695  max resid 0.0004064621 
    ## ... Similar to previous best
    ## Run 500 stress 8.614658e-05 
    ## ... Procrustes: rmse 0.000181301  max resid 0.0003283087 
    ## ... Similar to previous best
    ## *** Best solution repeated 129 times

    ## Warning in metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Mixed lakes
SD_beta_env_M_NMDS <- metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.001279543 
    ## Run 1 stress 0.1491957 
    ## Run 2 stress 9.762649e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04320589  max resid 0.05937883 
    ## Run 3 stress 0.1491957 
    ## Run 4 stress 0.000356625 
    ## ... Procrustes: rmse 0.01374414  max resid 0.01882348 
    ## Run 5 stress 0.0004295945 
    ## ... Procrustes: rmse 0.01509518  max resid 0.02068768 
    ## Run 6 stress 4.804319e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.00426535  max resid 0.005650903 
    ## ... Similar to previous best
    ## Run 7 stress 0.1990774 
    ## Run 8 stress 0.0004937168 
    ## ... Procrustes: rmse 0.01557862  max resid 0.02836298 
    ## Run 9 stress 9.516078e-05 
    ## ... Procrustes: rmse 0.004291777  max resid 0.005636141 
    ## ... Similar to previous best
    ## Run 10 stress 0.001384264 
    ## Run 11 stress 0.1491957 
    ## Run 12 stress 0.001237562 
    ## Run 13 stress 0.2835383 
    ## Run 14 stress 0.001401045 
    ## Run 15 stress 5.62416e-05 
    ## ... Procrustes: rmse 0.00428589  max resid 0.005837404 
    ## ... Similar to previous best
    ## Run 16 stress 9.778188e-05 
    ## ... Procrustes: rmse 0.004227215  max resid 0.005628186 
    ## ... Similar to previous best
    ## Run 17 stress 0.001208643 
    ## Run 18 stress 9.326802e-05 
    ## ... Procrustes: rmse 0.004263052  max resid 0.005660176 
    ## ... Similar to previous best
    ## Run 19 stress 0.0008373243 
    ## Run 20 stress 8.350957e-05 
    ## ... Procrustes: rmse 0.004298382  max resid 0.0059269 
    ## ... Similar to previous best
    ## Run 21 stress 0.0002846142 
    ## ... Procrustes: rmse 0.01188079  max resid 0.02301773 
    ## Run 22 stress 9.431482e-05 
    ## ... Procrustes: rmse 0.003976016  max resid 0.005495077 
    ## ... Similar to previous best
    ## Run 23 stress 9.085513e-05 
    ## ... Procrustes: rmse 0.004216624  max resid 0.005713156 
    ## ... Similar to previous best
    ## Run 24 stress 0.001301258 
    ## Run 25 stress 0.1491957 
    ## Run 26 stress 0.0004670354 
    ## ... Procrustes: rmse 0.0152096  max resid 0.02783262 
    ## Run 27 stress 0.0003022224 
    ## ... Procrustes: rmse 0.01225363  max resid 0.02355677 
    ## Run 28 stress 0.1776014 
    ## Run 29 stress 0.001332924 
    ## Run 30 stress 9.235058e-05 
    ## ... Procrustes: rmse 0.004236583  max resid 0.005654263 
    ## ... Similar to previous best
    ## Run 31 stress 9.348312e-05 
    ## ... Procrustes: rmse 0.004280716  max resid 0.005898168 
    ## ... Similar to previous best
    ## Run 32 stress 0.0006298121 
    ## Run 33 stress 0.001304269 
    ## Run 34 stress 0.0001028583 
    ## ... Procrustes: rmse 0.007466462  max resid 0.0162252 
    ## Run 35 stress 9.292771e-05 
    ## ... Procrustes: rmse 0.004194562  max resid 0.005669187 
    ## ... Similar to previous best
    ## Run 36 stress 0.1990774 
    ## Run 37 stress 0.001345982 
    ## Run 38 stress 0.0004727584 
    ## ... Procrustes: rmse 0.01528812  max resid 0.02794481 
    ## Run 39 stress 0.001188554 
    ## Run 40 stress 0.0003575309 
    ## ... Procrustes: rmse 0.01330739  max resid 0.02509205 
    ## Run 41 stress 0.001147208 
    ## Run 42 stress 0.2848446 
    ## Run 43 stress 0.0004657721 
    ## ... Procrustes: rmse 0.0151878  max resid 0.02780422 
    ## Run 44 stress 0.001310095 
    ## Run 45 stress 0.001284587 
    ## Run 46 stress 0.001454226 
    ## Run 47 stress 0.3083099 
    ## Run 48 stress 0.001348134 
    ## Run 49 stress 8.542267e-05 
    ## ... Procrustes: rmse 0.004258454  max resid 0.005731958 
    ## ... Similar to previous best
    ## Run 50 stress 0.001122843 
    ## Run 51 stress 0.177601 
    ## Run 52 stress 0.001139778 
    ## Run 53 stress 0.0003253132 
    ## ... Procrustes: rmse 0.01751536  max resid 0.02417947 
    ## Run 54 stress 9.300536e-05 
    ## ... Procrustes: rmse 0.004207622  max resid 0.005661339 
    ## ... Similar to previous best
    ## Run 55 stress 0.001365328 
    ## Run 56 stress 9.294189e-05 
    ## ... Procrustes: rmse 0.004159219  max resid 0.005649464 
    ## ... Similar to previous best
    ## Run 57 stress 0.0004422944 
    ## ... Procrustes: rmse 0.01480018  max resid 0.02724643 
    ## Run 58 stress 9.698808e-05 
    ## ... Procrustes: rmse 0.004188252  max resid 0.005654281 
    ## ... Similar to previous best
    ## Run 59 stress 0.2568281 
    ## Run 60 stress 0.001307219 
    ## Run 61 stress 0.001373172 
    ## Run 62 stress 0.2528294 
    ## Run 63 stress 0.0005063563 
    ## ... Procrustes: rmse 0.01583965  max resid 0.02873403 
    ## Run 64 stress 9.190961e-05 
    ## ... Procrustes: rmse 0.004167898  max resid 0.005670404 
    ## ... Similar to previous best
    ## Run 65 stress 0.001302507 
    ## Run 66 stress 0.001056076 
    ## Run 67 stress 0.0003620909 
    ## ... Procrustes: rmse 0.01339428  max resid 0.02521773 
    ## Run 68 stress 9.821059e-05 
    ## ... Procrustes: rmse 0.004218013  max resid 0.005627463 
    ## ... Similar to previous best
    ## Run 69 stress 9.527855e-05 
    ## ... Procrustes: rmse 0.004270097  max resid 0.005649631 
    ## ... Similar to previous best
    ## Run 70 stress 0.00123041 
    ## Run 71 stress 0.001335925 
    ## Run 72 stress 0.1491957 
    ## Run 73 stress 7.352163e-05 
    ## ... Procrustes: rmse 0.003246514  max resid 0.004437743 
    ## ... Similar to previous best
    ## Run 74 stress 0.0004732307 
    ## ... Procrustes: rmse 0.01530076  max resid 0.02797057 
    ## Run 75 stress 0.1491957 
    ## Run 76 stress 0.1776019 
    ## Run 77 stress 0.0003459327 
    ## ... Procrustes: rmse 0.01819764  max resid 0.02512338 
    ## Run 78 stress 0.0004800318 
    ## ... Procrustes: rmse 0.01541793  max resid 0.02813251 
    ## Run 79 stress 0.1491957 
    ## Run 80 stress 0.0004516404 
    ## ... Procrustes: rmse 0.01495574  max resid 0.02747022 
    ## Run 81 stress 0.001358299 
    ## Run 82 stress 0.1776019 
    ## Run 83 stress 8.588112e-05 
    ## ... Procrustes: rmse 0.004166104  max resid 0.005664836 
    ## ... Similar to previous best
    ## Run 84 stress 8.813869e-05 
    ## ... Procrustes: rmse 0.004241335  max resid 0.005691419 
    ## ... Similar to previous best
    ## Run 85 stress 7.749049e-05 
    ## ... Procrustes: rmse 0.004207879  max resid 0.005711599 
    ## ... Similar to previous best
    ## Run 86 stress 0.001337179 
    ## Run 87 stress 9.799616e-05 
    ## ... Procrustes: rmse 0.004195885  max resid 0.005648555 
    ## ... Similar to previous best
    ## Run 88 stress 9.046353e-05 
    ## ... Procrustes: rmse 0.002203647  max resid 0.003030944 
    ## ... Similar to previous best
    ## Run 89 stress 9.24748e-05 
    ## ... Procrustes: rmse 0.004197194  max resid 0.005666561 
    ## ... Similar to previous best
    ## Run 90 stress 0.001163378 
    ## Run 91 stress 0.1491957 
    ## Run 92 stress 7.219471e-05 
    ## ... Procrustes: rmse 0.004193798  max resid 0.00572597 
    ## ... Similar to previous best
    ## Run 93 stress 0.1491957 
    ## Run 94 stress 0.001271236 
    ## Run 95 stress 9.289643e-05 
    ## ... Procrustes: rmse 0.004217763  max resid 0.005662015 
    ## ... Similar to previous best
    ## Run 96 stress 0.001394181 
    ## Run 97 stress 7.916229e-05 
    ## ... Procrustes: rmse 0.004190907  max resid 0.005714228 
    ## ... Similar to previous best
    ## Run 98 stress 0.001287337 
    ## Run 99 stress 0.0008163161 
    ## Run 100 stress 8.143367e-05 
    ## ... Procrustes: rmse 0.004161234  max resid 0.007097907 
    ## ... Similar to previous best
    ## Run 101 stress 9.323684e-05 
    ## ... Procrustes: rmse 0.004240124  max resid 0.00567048 
    ## ... Similar to previous best
    ## Run 102 stress 9.605299e-05 
    ## ... Procrustes: rmse 0.004275215  max resid 0.00567189 
    ## ... Similar to previous best
    ## Run 103 stress 0.001299112 
    ## Run 104 stress 9.638351e-05 
    ## ... Procrustes: rmse 0.004183013  max resid 0.00564722 
    ## ... Similar to previous best
    ## Run 105 stress 8.651245e-05 
    ## ... Procrustes: rmse 0.004253731  max resid 0.005698109 
    ## ... Similar to previous best
    ## Run 106 stress 9.274631e-05 
    ## ... Procrustes: rmse 0.0041999  max resid 0.005682251 
    ## ... Similar to previous best
    ## Run 107 stress 0.001443484 
    ## Run 108 stress 9.931764e-05 
    ## ... Procrustes: rmse 0.00416588  max resid 0.005665091 
    ## ... Similar to previous best
    ## Run 109 stress 0.0004647503 
    ## ... Procrustes: rmse 0.01517183  max resid 0.02777889 
    ## Run 110 stress 8.519618e-05 
    ## ... Procrustes: rmse 0.004185795  max resid 0.005691848 
    ## ... Similar to previous best
    ## Run 111 stress 0.001228269 
    ## Run 112 stress 0.2842805 
    ## Run 113 stress 0.001205343 
    ## Run 114 stress 0.2528294 
    ## Run 115 stress 0.1990774 
    ## Run 116 stress 9.942907e-05 
    ## ... Procrustes: rmse 0.004214597  max resid 0.005619369 
    ## ... Similar to previous best
    ## Run 117 stress 0.1491957 
    ## Run 118 stress 9.145132e-05 
    ## ... Procrustes: rmse 0.00428874  max resid 0.005657054 
    ## ... Similar to previous best
    ## Run 119 stress 0.0004145284 
    ## ... Procrustes: rmse 0.01405197  max resid 0.02616961 
    ## Run 120 stress 0.0002325387 
    ## ... Procrustes: rmse 0.01069921  max resid 0.02126079 
    ## Run 121 stress 0.001165792 
    ## Run 122 stress 0.0004552682 
    ## ... Procrustes: rmse 0.01500176  max resid 0.02753512 
    ## Run 123 stress 9.216127e-05 
    ## ... Procrustes: rmse 0.004223945  max resid 0.005645461 
    ## ... Similar to previous best
    ## Run 124 stress 0.0005129369 
    ## ... Procrustes: rmse 0.0154659  max resid 0.02820028 
    ## Run 125 stress 9.533581e-05 
    ## ... Procrustes: rmse 0.004197562  max resid 0.005656546 
    ## ... Similar to previous best
    ## Run 126 stress 0.2440272 
    ## Run 127 stress 0.2842805 
    ## Run 128 stress 0.001183192 
    ## Run 129 stress 9.027152e-05 
    ## ... Procrustes: rmse 0.004175686  max resid 0.005686214 
    ## ... Similar to previous best
    ## Run 130 stress 9.338598e-05 
    ## ... Procrustes: rmse 0.004256988  max resid 0.005660038 
    ## ... Similar to previous best
    ## Run 131 stress 0.001343457 
    ## Run 132 stress 9.675601e-05 
    ## ... Procrustes: rmse 0.004228406  max resid 0.005631983 
    ## ... Similar to previous best
    ## Run 133 stress 0.001370426 
    ## Run 134 stress 9.96975e-05 
    ## ... Procrustes: rmse 0.006921008  max resid 0.009562922 
    ## Run 135 stress 9.435755e-05 
    ## ... Procrustes: rmse 0.004219251  max resid 0.005644543 
    ## ... Similar to previous best
    ## Run 136 stress 9.632929e-05 
    ## ... Procrustes: rmse 0.004218927  max resid 0.005652606 
    ## ... Similar to previous best
    ## Run 137 stress 0.001254679 
    ## Run 138 stress 0.001249216 
    ## Run 139 stress 0.3119769 
    ## Run 140 stress 0.0005001353 
    ## ... Procrustes: rmse 0.01574163  max resid 0.02859393 
    ## Run 141 stress 9.511041e-05 
    ## ... Procrustes: rmse 0.004162658  max resid 0.005657644 
    ## ... Similar to previous best
    ## Run 142 stress 0.0001223643 
    ## ... Procrustes: rmse 0.008037987  max resid 0.01715684 
    ## Run 143 stress 0.3071914 
    ## Run 144 stress 0.001174809 
    ## Run 145 stress 0.001149814 
    ## Run 146 stress 9.752622e-05 
    ## ... Procrustes: rmse 0.004223538  max resid 0.005635062 
    ## ... Similar to previous best
    ## Run 147 stress 0.00083599 
    ## Run 148 stress 0.3083098 
    ## Run 149 stress 0.001346018 
    ## Run 150 stress 8.855702e-05 
    ## ... Procrustes: rmse 0.004199797  max resid 0.005667214 
    ## ... Similar to previous best
    ## Run 151 stress 0.0005038981 
    ## ... Procrustes: rmse 0.01580121  max resid 0.02867905 
    ## Run 152 stress 0.001296193 
    ## Run 153 stress 0.001301512 
    ## Run 154 stress 0.001320177 
    ## Run 155 stress 0.001286364 
    ## Run 156 stress 8.601798e-05 
    ## ... Procrustes: rmse 0.004229361  max resid 0.0056725 
    ## ... Similar to previous best
    ## Run 157 stress 0.257864 
    ## Run 158 stress 9.143555e-05 
    ## ... Procrustes: rmse 0.004179723  max resid 0.005672928 
    ## ... Similar to previous best
    ## Run 159 stress 0.001359077 
    ## Run 160 stress 0.00115971 
    ## Run 161 stress 0.001387657 
    ## Run 162 stress 0.001336683 
    ## Run 163 stress 8.698676e-05 
    ## ... Procrustes: rmse 0.00419441  max resid 0.005674132 
    ## ... Similar to previous best
    ## Run 164 stress 0.001362447 
    ## Run 165 stress 8.19126e-05 
    ## ... Procrustes: rmse 0.004207974  max resid 0.005690218 
    ## ... Similar to previous best
    ## Run 166 stress 0.1990774 
    ## Run 167 stress 8.028737e-05 
    ## ... Procrustes: rmse 0.004040222  max resid 0.00551666 
    ## ... Similar to previous best
    ## Run 168 stress 0.0009881192 
    ## Run 169 stress 9.74946e-05 
    ## ... Procrustes: rmse 0.004152321  max resid 0.005635924 
    ## ... Similar to previous best
    ## Run 170 stress 8.526878e-05 
    ## ... Procrustes: rmse 0.004269298  max resid 0.00568361 
    ## ... Similar to previous best
    ## Run 171 stress 0.1990774 
    ## Run 172 stress 0.00045381 
    ## ... Procrustes: rmse 0.01495119  max resid 0.02744181 
    ## Run 173 stress 9.393168e-05 
    ## ... Procrustes: rmse 0.005269119  max resid 0.01216449 
    ## Run 174 stress 0.001230057 
    ## Run 175 stress 0.001380046 
    ## Run 176 stress 0.001388376 
    ## Run 177 stress 0.0004803589 
    ## ... Procrustes: rmse 0.01542332  max resid 0.02814119 
    ## Run 178 stress 9.62518e-05 
    ## ... Procrustes: rmse 0.004155223  max resid 0.005636952 
    ## ... Similar to previous best
    ## Run 179 stress 0.0003074311 
    ## ... Procrustes: rmse 0.01235351  max resid 0.02370309 
    ## Run 180 stress 8.939364e-05 
    ## ... Procrustes: rmse 0.00416172  max resid 0.005654467 
    ## ... Similar to previous best
    ## Run 181 stress 9.904525e-05 
    ## ... Procrustes: rmse 0.004191268  max resid 0.005709189 
    ## ... Similar to previous best
    ## Run 182 stress 8.221292e-05 
    ## ... Procrustes: rmse 0.004228092  max resid 0.006540885 
    ## ... Similar to previous best
    ## Run 183 stress 9.854883e-05 
    ## ... Procrustes: rmse 0.004214194  max resid 0.005624822 
    ## ... Similar to previous best
    ## Run 184 stress 0.001354909 
    ## Run 185 stress 0.2797318 
    ## Run 186 stress 0.001346089 
    ## Run 187 stress 8.675306e-05 
    ## ... Procrustes: rmse 0.00425849  max resid 0.006228267 
    ## ... Similar to previous best
    ## Run 188 stress 0.1491957 
    ## Run 189 stress 0.2373687 
    ## Run 190 stress 9.784084e-05 
    ## ... Procrustes: rmse 0.004148246  max resid 0.005624622 
    ## ... Similar to previous best
    ## Run 191 stress 9.07979e-05 
    ## ... Procrustes: rmse 0.004210999  max resid 0.005667486 
    ## ... Similar to previous best
    ## Run 192 stress 0.1491957 
    ## Run 193 stress 0.001332156 
    ## Run 194 stress 0.001247024 
    ## Run 195 stress 0.0006663784 
    ## Run 196 stress 9.544321e-05 
    ## ... Procrustes: rmse 0.004230705  max resid 0.005639075 
    ## ... Similar to previous best
    ## Run 197 stress 9.909164e-05 
    ## ... Procrustes: rmse 0.00422002  max resid 0.005694542 
    ## ... Similar to previous best
    ## Run 198 stress 0.001142417 
    ## Run 199 stress 0.0005239165 
    ## ... Procrustes: rmse 0.0155813  max resid 0.02836549 
    ## Run 200 stress 0.001351893 
    ## Run 201 stress 0.001280113 
    ## Run 202 stress 0.001012953 
    ## Run 203 stress 0.001363667 
    ## Run 204 stress 0.2852163 
    ## Run 205 stress 0.0002426586 
    ## ... Procrustes: rmse 0.01452774  max resid 0.02004726 
    ## Run 206 stress 9.02853e-05 
    ## ... Procrustes: rmse 0.004268613  max resid 0.005661971 
    ## ... Similar to previous best
    ## Run 207 stress 9.154285e-05 
    ## ... Procrustes: rmse 0.004221998  max resid 0.005659602 
    ## ... Similar to previous best
    ## Run 208 stress 9.834107e-05 
    ## ... Procrustes: rmse 0.004152638  max resid 0.005635404 
    ## ... Similar to previous best
    ## Run 209 stress 0.1491957 
    ## Run 210 stress 0.2848522 
    ## Run 211 stress 0.001223379 
    ## Run 212 stress 8.506153e-05 
    ## ... Procrustes: rmse 0.004125682  max resid 0.005619744 
    ## ... Similar to previous best
    ## Run 213 stress 0.001086027 
    ## Run 214 stress 9.988884e-05 
    ## ... Procrustes: rmse 0.00422503  max resid 0.00562132 
    ## ... Similar to previous best
    ## Run 215 stress 0.0004865057 
    ## ... Procrustes: rmse 0.01552415  max resid 0.02828391 
    ## Run 216 stress 9.600819e-05 
    ## ... Procrustes: rmse 0.00427114  max resid 0.005652548 
    ## ... Similar to previous best
    ## Run 217 stress 8.99076e-05 
    ## ... Procrustes: rmse 0.004216758  max resid 0.005776893 
    ## ... Similar to previous best
    ## Run 218 stress 0.1491957 
    ## Run 219 stress 0.001221745 
    ## Run 220 stress 7.947772e-05 
    ## ... Procrustes: rmse 0.004243352  max resid 0.005689691 
    ## ... Similar to previous best
    ## Run 221 stress 0.0004783168 
    ## ... Procrustes: rmse 0.01539267  max resid 0.02809483 
    ## Run 222 stress 0.001387909 
    ## Run 223 stress 0.1491957 
    ## Run 224 stress 0.1990774 
    ## Run 225 stress 0.0004729311 
    ## ... Procrustes: rmse 0.0153044  max resid 0.02797284 
    ## Run 226 stress 0.1990774 
    ## Run 227 stress 0.2797324 
    ## Run 228 stress 9.524617e-05 
    ## ... Procrustes: rmse 0.004157119  max resid 0.005644816 
    ## ... Similar to previous best
    ## Run 229 stress 0.1491957 
    ## Run 230 stress 0.001190186 
    ## Run 231 stress 0.2842805 
    ## Run 232 stress 0.001381879 
    ## Run 233 stress 9.001414e-05 
    ## ... Procrustes: rmse 0.004262383  max resid 0.005713975 
    ## ... Similar to previous best
    ## Run 234 stress 0.001352285 
    ## Run 235 stress 0.001272314 
    ## Run 236 stress 9.599056e-05 
    ## ... Procrustes: rmse 0.004292452  max resid 0.005638877 
    ## ... Similar to previous best
    ## Run 237 stress 0.0004011434 
    ## ... Procrustes: rmse 0.01967077  max resid 0.02716334 
    ## Run 238 stress 0.0005055353 
    ## ... Procrustes: rmse 0.01578979  max resid 0.02866339 
    ## Run 239 stress 9.263841e-05 
    ## ... Procrustes: rmse 0.004292302  max resid 0.005644092 
    ## ... Similar to previous best
    ## Run 240 stress 9.060163e-05 
    ## ... Procrustes: rmse 0.004222881  max resid 0.00565172 
    ## ... Similar to previous best
    ## Run 241 stress 0.2361461 
    ## Run 242 stress 0.001316696 
    ## Run 243 stress 0.1776019 
    ## Run 244 stress 9.362917e-05 
    ## ... Procrustes: rmse 0.004228113  max resid 0.00564058 
    ## ... Similar to previous best
    ## Run 245 stress 0.0009108137 
    ## Run 246 stress 0.1491957 
    ## Run 247 stress 0.3083098 
    ## Run 248 stress 9.469496e-05 
    ## ... Procrustes: rmse 0.004229471  max resid 0.005650677 
    ## ... Similar to previous best
    ## Run 249 stress 0.1491957 
    ## Run 250 stress 0.001266043 
    ## Run 251 stress 0.0004860931 
    ## ... Procrustes: rmse 0.01551672  max resid 0.02827194 
    ## Run 252 stress 0.2851281 
    ## Run 253 stress 9.919878e-05 
    ## ... Procrustes: rmse 0.004275839  max resid 0.005975443 
    ## ... Similar to previous best
    ## Run 254 stress 0.001296038 
    ## Run 255 stress 0.1491957 
    ## Run 256 stress 0.0004243711 
    ## ... Procrustes: rmse 0.0144904  max resid 0.02680054 
    ## Run 257 stress 0.1491957 
    ## Run 258 stress 0.1990774 
    ## Run 259 stress 0.3120129 
    ## Run 260 stress 7.025392e-05 
    ## ... Procrustes: rmse 0.004242358  max resid 0.00576465 
    ## ... Similar to previous best
    ## Run 261 stress 0.0003468659 
    ## ... Procrustes: rmse 0.01310278  max resid 0.02479546 
    ## Run 262 stress 9.215166e-05 
    ## ... Procrustes: rmse 0.004229353  max resid 0.005656054 
    ## ... Similar to previous best
    ## Run 263 stress 8.940359e-05 
    ## ... Procrustes: rmse 0.004197463  max resid 0.005839426 
    ## ... Similar to previous best
    ## Run 264 stress 9.323314e-05 
    ## ... Procrustes: rmse 0.004188483  max resid 0.005657772 
    ## ... Similar to previous best
    ## Run 265 stress 8.949942e-05 
    ## ... Procrustes: rmse 0.00418537  max resid 0.00568068 
    ## ... Similar to previous best
    ## Run 266 stress 9.257475e-05 
    ## ... Procrustes: rmse 0.00415976  max resid 0.005649224 
    ## ... Similar to previous best
    ## Run 267 stress 0.1491957 
    ## Run 268 stress 9.336499e-05 
    ## ... Procrustes: rmse 0.004234463  max resid 0.005647939 
    ## ... Similar to previous best
    ## Run 269 stress 0.2520602 
    ## Run 270 stress 0.0004328042 
    ## ... Procrustes: rmse 0.01463833  max resid 0.02701554 
    ## Run 271 stress 0.001161096 
    ## Run 272 stress 9.200176e-05 
    ## ... Procrustes: rmse 0.004297442  max resid 0.006318225 
    ## ... Similar to previous best
    ## Run 273 stress 0.0002442042 
    ## ... Procrustes: rmse 0.01460155  max resid 0.02014946 
    ## Run 274 stress 0.001287076 
    ## Run 275 stress 0.001309442 
    ## Run 276 stress 0.001404128 
    ## Run 277 stress 0.0004746065 
    ## ... Procrustes: rmse 0.02205269  max resid 0.0304604 
    ## Run 278 stress 0.001337305 
    ## Run 279 stress 8.557037e-05 
    ## ... Procrustes: rmse 0.004237576  max resid 0.005671555 
    ## ... Similar to previous best
    ## Run 280 stress 0.1491957 
    ## Run 281 stress 0.0006271444 
    ## Run 282 stress 9.297351e-05 
    ## ... Procrustes: rmse 0.004255239  max resid 0.005657112 
    ## ... Similar to previous best
    ## Run 283 stress 0.001336545 
    ## Run 284 stress 9.171603e-05 
    ## ... Procrustes: rmse 0.00423472  max resid 0.005655835 
    ## ... Similar to previous best
    ## Run 285 stress 8.764071e-05 
    ## ... Procrustes: rmse 0.004190193  max resid 0.005679518 
    ## ... Similar to previous best
    ## Run 286 stress 0.001222809 
    ## Run 287 stress 0.001169971 
    ## Run 288 stress 0.000432952 
    ## ... Procrustes: rmse 0.01464312  max resid 0.02702045 
    ## Run 289 stress 9.022461e-05 
    ## ... Procrustes: rmse 0.003862866  max resid 0.005256946 
    ## ... Similar to previous best
    ## Run 290 stress 0.2822095 
    ## Run 291 stress 7.486424e-05 
    ## ... Procrustes: rmse 0.004249515  max resid 0.005726017 
    ## ... Similar to previous best
    ## Run 292 stress 8.800712e-05 
    ## ... Procrustes: rmse 0.004269054  max resid 0.005674117 
    ## ... Similar to previous best
    ## Run 293 stress 0.0001051053 
    ## ... Procrustes: rmse 0.007533762  max resid 0.01633643 
    ## Run 294 stress 8.930871e-05 
    ## ... Procrustes: rmse 0.004259056  max resid 0.005664638 
    ## ... Similar to previous best
    ## Run 295 stress 0.001216415 
    ## Run 296 stress 8.980826e-05 
    ## ... Procrustes: rmse 0.00380683  max resid 0.005177603 
    ## ... Similar to previous best
    ## Run 297 stress 0.2520602 
    ## Run 298 stress 0.001389847 
    ## Run 299 stress 0.2570107 
    ## Run 300 stress 0.001481734 
    ## Run 301 stress 0.1491957 
    ## Run 302 stress 0.0004865938 
    ## ... Procrustes: rmse 0.01552531  max resid 0.02828501 
    ## Run 303 stress 0.1491957 
    ## Run 304 stress 9.050487e-05 
    ## ... Procrustes: rmse 0.004263348  max resid 0.00566429 
    ## ... Similar to previous best
    ## Run 305 stress 8.649336e-05 
    ## ... Procrustes: rmse 0.004518682  max resid 0.01016258 
    ## Run 306 stress 9.861328e-05 
    ## ... Procrustes: rmse 0.004149179  max resid 0.005626332 
    ## ... Similar to previous best
    ## Run 307 stress 0.1491957 
    ## Run 308 stress 0.0003233 
    ## ... Procrustes: rmse 0.012666  max resid 0.02415926 
    ## Run 309 stress 0.1491957 
    ## Run 310 stress 9.497577e-05 
    ## ... Procrustes: rmse 0.004215026  max resid 0.005652702 
    ## ... Similar to previous best
    ## Run 311 stress 9.309234e-05 
    ## ... Procrustes: rmse 0.004226777  max resid 0.005645784 
    ## ... Similar to previous best
    ## Run 312 stress 7.786228e-05 
    ## ... Procrustes: rmse 0.003380945  max resid 0.004613401 
    ## ... Similar to previous best
    ## Run 313 stress 0.0004469817 
    ## ... Procrustes: rmse 0.01487816  max resid 0.02735757 
    ## Run 314 stress 0.001341818 
    ## Run 315 stress 9.462142e-05 
    ## ... Procrustes: rmse 0.004221819  max resid 0.00564195 
    ## ... Similar to previous best
    ## Run 316 stress 0.001299058 
    ## Run 317 stress 9.961569e-05 
    ## ... Procrustes: rmse 0.004291763  max resid 0.005621737 
    ## ... Similar to previous best
    ## Run 318 stress 7.911359e-05 
    ## ... Procrustes: rmse 0.004180889  max resid 0.005693709 
    ## ... Similar to previous best
    ## Run 319 stress 8.080442e-05 
    ## ... Procrustes: rmse 0.004208717  max resid 0.005709414 
    ## ... Similar to previous best
    ## Run 320 stress 8.832751e-05 
    ## ... Procrustes: rmse 0.004220382  max resid 0.005665882 
    ## ... Similar to previous best
    ## Run 321 stress 0.001410468 
    ## Run 322 stress 0.2520602 
    ## Run 323 stress 0.2373687 
    ## Run 324 stress 0.0004321095 
    ## ... Procrustes: rmse 0.01462845  max resid 0.02700012 
    ## Run 325 stress 0.001207352 
    ## Run 326 stress 9.663605e-05 
    ## ... Procrustes: rmse 0.004151612  max resid 0.005631605 
    ## ... Similar to previous best
    ## Run 327 stress 0.001370226 
    ## Run 328 stress 7.787799e-05 
    ## ... Procrustes: rmse 0.004250707  max resid 0.005728344 
    ## ... Similar to previous best
    ## Run 329 stress 0.001413401 
    ## Run 330 stress 0.1491957 
    ## Run 331 stress 0.001098377 
    ## Run 332 stress 9.599212e-05 
    ## ... Procrustes: rmse 0.004279978  max resid 0.005958958 
    ## ... Similar to previous best
    ## Run 333 stress 0.2848512 
    ## Run 334 stress 0.1990774 
    ## Run 335 stress 0.00105949 
    ## Run 336 stress 0.1491957 
    ## Run 337 stress 9.57755e-05 
    ## ... Procrustes: rmse 0.004217109  max resid 0.005638569 
    ## ... Similar to previous best
    ## Run 338 stress 0.1491957 
    ## Run 339 stress 9.77485e-05 
    ## ... Procrustes: rmse 0.003014176  max resid 0.004136581 
    ## ... Similar to previous best
    ## Run 340 stress 0.0004939477 
    ## ... Procrustes: rmse 0.01563885  max resid 0.02844724 
    ## Run 341 stress 0.001462056 
    ## Run 342 stress 0.00139222 
    ## Run 343 stress 0.2509018 
    ## Run 344 stress 0.001021317 
    ## Run 345 stress 9.397505e-05 
    ## ... Procrustes: rmse 0.004231942  max resid 0.005644779 
    ## ... Similar to previous best
    ## Run 346 stress 0.3083098 
    ## Run 347 stress 0.0001298622 
    ## ... Procrustes: rmse 0.00824824  max resid 0.01749324 
    ## Run 348 stress 0.1776017 
    ## Run 349 stress 0.001184044 
    ## Run 350 stress 0.0004400787 
    ## ... Procrustes: rmse 0.01475009  max resid 0.02717374 
    ## Run 351 stress 0.0004963723 
    ## ... Procrustes: rmse 0.01567593  max resid 0.02850091 
    ## Run 352 stress 9.735614e-05 
    ## ... Procrustes: rmse 0.004181095  max resid 0.005641161 
    ## ... Similar to previous best
    ## Run 353 stress 8.677154e-05 
    ## ... Procrustes: rmse 0.004288619  max resid 0.005668165 
    ## ... Similar to previous best
    ## Run 354 stress 0.1491957 
    ## Run 355 stress 0.2520602 
    ## Run 356 stress 9.665588e-05 
    ## ... Procrustes: rmse 0.004214767  max resid 0.005633999 
    ## ... Similar to previous best
    ## Run 357 stress 9.191955e-05 
    ## ... Procrustes: rmse 0.004212815  max resid 0.005686045 
    ## ... Similar to previous best
    ## Run 358 stress 9.304164e-05 
    ## ... Procrustes: rmse 0.004231454  max resid 0.005653573 
    ## ... Similar to previous best
    ## Run 359 stress 0.001378994 
    ## Run 360 stress 6.640782e-05 
    ## ... Procrustes: rmse 0.003156085  max resid 0.004314617 
    ## ... Similar to previous best
    ## Run 361 stress 0.0005264906 
    ## ... Procrustes: rmse 0.01598392  max resid 0.02894636 
    ## Run 362 stress 0.001453249 
    ## Run 363 stress 0.001177707 
    ## Run 364 stress 9.679752e-05 
    ## ... Procrustes: rmse 0.004221307  max resid 0.005631176 
    ## ... Similar to previous best
    ## Run 365 stress 9.230486e-05 
    ## ... Procrustes: rmse 0.004291788  max resid 0.005644899 
    ## ... Similar to previous best
    ## Run 366 stress 0.001422543 
    ## Run 367 stress 9.779622e-05 
    ## ... Procrustes: rmse 0.00419206  max resid 0.005642734 
    ## ... Similar to previous best
    ## Run 368 stress 0.001157605 
    ## Run 369 stress 9.410914e-05 
    ## ... Procrustes: rmse 0.004168555  max resid 0.005666914 
    ## ... Similar to previous best
    ## Run 370 stress 9.221016e-05 
    ## ... Procrustes: rmse 0.004218602  max resid 0.005648372 
    ## ... Similar to previous best
    ## Run 371 stress 0.0003408822 
    ## ... Procrustes: rmse 0.01801927  max resid 0.02487682 
    ## Run 372 stress 9.150775e-05 
    ## ... Procrustes: rmse 0.004259032  max resid 0.005663415 
    ## ... Similar to previous best
    ## Run 373 stress 0.001335982 
    ## Run 374 stress 0.001283384 
    ## Run 375 stress 0.001040528 
    ## Run 376 stress 0.001126151 
    ## Run 377 stress 9.686108e-05 
    ## ... Procrustes: rmse 0.001704051  max resid 0.002331213 
    ## ... Similar to previous best
    ## Run 378 stress 0.1776021 
    ## Run 379 stress 0.1491957 
    ## Run 380 stress 0.001287109 
    ## Run 381 stress 0.0003981388 
    ## ... Procrustes: rmse 0.01983655  max resid 0.02739166 
    ## Run 382 stress 8.963341e-05 
    ## ... Procrustes: rmse 0.004259348  max resid 0.005700339 
    ## ... Similar to previous best
    ## Run 383 stress 8.060349e-05 
    ## ... Procrustes: rmse 0.004564767  max resid 0.01032752 
    ## Run 384 stress 0.0004797434 
    ## ... Procrustes: rmse 0.0154139  max resid 0.02812534 
    ## Run 385 stress 8.32222e-05 
    ## ... Procrustes: rmse 0.004178683  max resid 0.005687493 
    ## ... Similar to previous best
    ## Run 386 stress 0.1491957 
    ## Run 387 stress 0.1491957 
    ## Run 388 stress 0.0004217885 
    ## ... Procrustes: rmse 0.01445057  max resid 0.02674329 
    ## Run 389 stress 8.187806e-05 
    ## ... Procrustes: rmse 0.004273024  max resid 0.005707447 
    ## ... Similar to previous best
    ## Run 390 stress 0.0001167283 
    ## ... Procrustes: rmse 0.007876238  max resid 0.0168959 
    ## Run 391 stress 0.1990774 
    ## Run 392 stress 0.001297405 
    ## Run 393 stress 9.617654e-05 
    ## ... Procrustes: rmse 0.004284275  max resid 0.005938882 
    ## ... Similar to previous best
    ## Run 394 stress 9.762926e-05 
    ## ... Procrustes: rmse 0.004221368  max resid 0.00562605 
    ## ... Similar to previous best
    ## Run 395 stress 9.125509e-05 
    ## ... Procrustes: rmse 0.00424861  max resid 0.005667446 
    ## ... Similar to previous best
    ## Run 396 stress 9.068348e-05 
    ## ... Procrustes: rmse 0.004259727  max resid 0.005668971 
    ## ... Similar to previous best
    ## Run 397 stress 0.2568281 
    ## Run 398 stress 0.001297229 
    ## Run 399 stress 0.001315044 
    ## Run 400 stress 8.605461e-05 
    ## ... Procrustes: rmse 0.004170746  max resid 0.005672541 
    ## ... Similar to previous best
    ## Run 401 stress 9.935029e-05 
    ## ... Procrustes: rmse 0.004201589  max resid 0.005636117 
    ## ... Similar to previous best
    ## Run 402 stress 0.001238196 
    ## Run 403 stress 0.001064927 
    ## Run 404 stress 9.817967e-05 
    ## ... Procrustes: rmse 0.004255453  max resid 0.005655386 
    ## ... Similar to previous best
    ## Run 405 stress 8.967831e-05 
    ## ... Procrustes: rmse 0.004233001  max resid 0.005675307 
    ## ... Similar to previous best
    ## Run 406 stress 0.1491957 
    ## Run 407 stress 0.00110512 
    ## Run 408 stress 8.956935e-05 
    ## ... Procrustes: rmse 0.004193176  max resid 0.00567195 
    ## ... Similar to previous best
    ## Run 409 stress 9.906552e-05 
    ## ... Procrustes: rmse 0.007272702  max resid 0.01590531 
    ## Run 410 stress 0.0003007369 
    ## ... Procrustes: rmse 0.01222421  max resid 0.02351351 
    ## Run 411 stress 7.964213e-05 
    ## ... Procrustes: rmse 0.004193227  max resid 0.005717888 
    ## ... Similar to previous best
    ## Run 412 stress 8.618508e-05 
    ## ... Procrustes: rmse 0.00146217  max resid 0.001999076 
    ## ... Similar to previous best
    ## Run 413 stress 8.585683e-05 
    ## ... Procrustes: rmse 0.004144461  max resid 0.007881182 
    ## ... Similar to previous best
    ## Run 414 stress 0.0007537202 
    ## Run 415 stress 0.0005259583 
    ## ... Procrustes: rmse 0.01614579  max resid 0.02917106 
    ## Run 416 stress 0.1491957 
    ## Run 417 stress 9.315712e-05 
    ## ... Procrustes: rmse 0.004229524  max resid 0.005684962 
    ## ... Similar to previous best
    ## Run 418 stress 0.1990774 
    ## Run 419 stress 9.547938e-05 
    ## ... Procrustes: rmse 0.004184196  max resid 0.005650773 
    ## ... Similar to previous best
    ## Run 420 stress 9.493446e-05 
    ## ... Procrustes: rmse 0.004264079  max resid 0.005663765 
    ## ... Similar to previous best
    ## Run 421 stress 0.00110623 
    ## Run 422 stress 0.001197682 
    ## Run 423 stress 0.001371717 
    ## Run 424 stress 0.1990774 
    ## Run 425 stress 8.49153e-05 
    ## ... Procrustes: rmse 0.004231731  max resid 0.005747122 
    ## ... Similar to previous best
    ## Run 426 stress 0.0004438025 
    ## ... Procrustes: rmse 0.01477701  max resid 0.02721652 
    ## Run 427 stress 8.823656e-05 
    ## ... Procrustes: rmse 0.004163411  max resid 0.00565766 
    ## ... Similar to previous best
    ## Run 428 stress 0.1776016 
    ## Run 429 stress 8.706572e-05 
    ## ... Procrustes: rmse 0.004300082  max resid 0.005977439 
    ## ... Similar to previous best
    ## Run 430 stress 0.2528294 
    ## Run 431 stress 0.0004826695 
    ## ... Procrustes: rmse 0.02227439  max resid 0.03076725 
    ## Run 432 stress 0.0004858674 
    ## ... Procrustes: rmse 0.01551306  max resid 0.02826724 
    ## Run 433 stress 0.2848519 
    ## Run 434 stress 0.2852162 
    ## Run 435 stress 9.768802e-05 
    ## ... Procrustes: rmse 0.004148656  max resid 0.005626072 
    ## ... Similar to previous best
    ## Run 436 stress 0.0004617159 
    ## ... Procrustes: rmse 0.01512134  max resid 0.02770647 
    ## Run 437 stress 9.816764e-05 
    ## ... Procrustes: rmse 0.0042872  max resid 0.006022025 
    ## ... Similar to previous best
    ## Run 438 stress 9.47485e-05 
    ## ... Procrustes: rmse 0.004156163  max resid 0.005642736 
    ## ... Similar to previous best
    ## Run 439 stress 0.0001301158 
    ## ... Procrustes: rmse 0.009411605  max resid 0.01297731 
    ## Run 440 stress 0.1776013 
    ## Run 441 stress 0.0001321369 
    ## ... Procrustes: rmse 0.008310645  max resid 0.01759257 
    ## Run 442 stress 0.0003660517 
    ## ... Procrustes: rmse 0.01346744  max resid 0.02532437 
    ## Run 443 stress 0.00103224 
    ## Run 444 stress 0.001374638 
    ## Run 445 stress 0.0004335913 
    ## ... Procrustes: rmse 0.01464481  max resid 0.02702316 
    ## Run 446 stress 0.0001305654 
    ## ... Procrustes: rmse 0.008268953  max resid 0.01752629 
    ## Run 447 stress 0.001207067 
    ## Run 448 stress 0.001270452 
    ## Run 449 stress 0.001352581 
    ## Run 450 stress 0.00125862 
    ## Run 451 stress 9.838678e-05 
    ## ... Procrustes: rmse 0.006366393  max resid 0.01433504 
    ## Run 452 stress 0.001220827 
    ## Run 453 stress 0.001398292 
    ## Run 454 stress 0.1776014 
    ## Run 455 stress 0.001403837 
    ## Run 456 stress 0.001329095 
    ## Run 457 stress 0.1990774 
    ## Run 458 stress 0.0008444508 
    ## Run 459 stress 0.1491957 
    ## Run 460 stress 8.435511e-05 
    ## ... Procrustes: rmse 0.004284452  max resid 0.005938615 
    ## ... Similar to previous best
    ## Run 461 stress 0.001258047 
    ## Run 462 stress 0.3083098 
    ## Run 463 stress 0.001330899 
    ## Run 464 stress 9.464846e-05 
    ## ... Procrustes: rmse 0.004219461  max resid 0.005652443 
    ## ... Similar to previous best
    ## Run 465 stress 0.1990774 
    ## Run 466 stress 0.001397263 
    ## Run 467 stress 0.001177381 
    ## Run 468 stress 9.902534e-05 
    ## ... Procrustes: rmse 0.004167558  max resid 0.005650796 
    ## ... Similar to previous best
    ## Run 469 stress 0.001002468 
    ## Run 470 stress 0.1491957 
    ## Run 471 stress 9.320777e-05 
    ## ... Procrustes: rmse 0.004226471  max resid 0.005644317 
    ## ... Similar to previous best
    ## Run 472 stress 8.551201e-05 
    ## ... Procrustes: rmse 0.005394789  max resid 0.01243538 
    ## Run 473 stress 0.1990774 
    ## Run 474 stress 0.001275544 
    ## Run 475 stress 0.001425237 
    ## Run 476 stress 0.1990774 
    ## Run 477 stress 0.0009049052 
    ## Run 478 stress 0.001143568 
    ## Run 479 stress 0.001379025 
    ## Run 480 stress 0.2528294 
    ## Run 481 stress 0.1990774 
    ## Run 482 stress 0.001328855 
    ## Run 483 stress 0.1776019 
    ## Run 484 stress 9.179438e-05 
    ## ... Procrustes: rmse 0.00422482  max resid 0.005649648 
    ## ... Similar to previous best
    ## Run 485 stress 9.190852e-05 
    ## ... Procrustes: rmse 0.004163453  max resid 0.005660033 
    ## ... Similar to previous best
    ## Run 486 stress 0.001440676 
    ## Run 487 stress 0.001214749 
    ## Run 488 stress 9.534572e-05 
    ## ... Procrustes: rmse 0.004230946  max resid 0.005649661 
    ## ... Similar to previous best
    ## Run 489 stress 0.001348409 
    ## Run 490 stress 0.00134157 
    ## Run 491 stress 0.2842805 
    ## Run 492 stress 9.466554e-05 
    ## ... Procrustes: rmse 0.004226802  max resid 0.005650557 
    ## ... Similar to previous best
    ## Run 493 stress 9.418802e-05 
    ## ... Procrustes: rmse 0.004221228  max resid 0.005641675 
    ## ... Similar to previous best
    ## Run 494 stress 0.0008538723 
    ## Run 495 stress 0.001306695 
    ## Run 496 stress 0.0004796433 
    ## ... Procrustes: rmse 0.0154139  max resid 0.02812546 
    ## Run 497 stress 0.0002846932 
    ## ... Procrustes: rmse 0.01513853  max resid 0.02088634 
    ## Run 498 stress 0.001291126 
    ## Run 499 stress 0.0004249708 
    ## ... Procrustes: rmse 0.01450737  max resid 0.02682506 
    ## Run 500 stress 9.143504e-05 
    ## ... Procrustes: rmse 0.004183525  max resid 0.005673764 
    ## ... Similar to previous best
    ## *** Best solution repeated 161 times

    ## Warning in metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
SD_beta_geo_NMDS <- metaMDS(SD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.09021173 
    ## Run 1 stress 0.1071322 
    ## Run 2 stress 0.09503447 
    ## Run 3 stress 0.09044608 
    ## ... Procrustes: rmse 0.01009228  max resid 0.03434769 
    ## Run 4 stress 0.1080637 
    ## Run 5 stress 0.09130086 
    ## Run 6 stress 0.08926078 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04363911  max resid 0.1170214 
    ## Run 7 stress 0.091783 
    ## Run 8 stress 0.1071321 
    ## Run 9 stress 0.1110415 
    ## Run 10 stress 0.0893854 
    ## ... Procrustes: rmse 0.01176039  max resid 0.03960219 
    ## Run 11 stress 0.1104469 
    ## Run 12 stress 0.08926074 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0003605091  max resid 0.001139885 
    ## ... Similar to previous best
    ## Run 13 stress 0.08938961 
    ## ... Procrustes: rmse 0.03591123  max resid 0.11824 
    ## Run 14 stress 0.08938986 
    ## ... Procrustes: rmse 0.03585916  max resid 0.1181583 
    ## Run 15 stress 0.09178285 
    ## Run 16 stress 0.09503418 
    ## Run 17 stress 0.08926078 
    ## ... Procrustes: rmse 0.0004138195  max resid 0.001176566 
    ## ... Similar to previous best
    ## Run 18 stress 0.1085124 
    ## Run 19 stress 0.08938552 
    ## ... Procrustes: rmse 0.01185375  max resid 0.03942573 
    ## Run 20 stress 0.08926105 
    ## ... Procrustes: rmse 0.0005803086  max resid 0.001694725 
    ## ... Similar to previous best
    ## Run 21 stress 0.09018708 
    ## Run 22 stress 0.09594302 
    ## Run 23 stress 0.1067508 
    ## Run 24 stress 0.08946329 
    ## ... Procrustes: rmse 0.03746582  max resid 0.1177368 
    ## Run 25 stress 0.08938963 
    ## ... Procrustes: rmse 0.03591918  max resid 0.1182624 
    ## Run 26 stress 0.08938963 
    ## ... Procrustes: rmse 0.03591155  max resid 0.1182548 
    ## Run 27 stress 0.08926078 
    ## ... Procrustes: rmse 0.0004051308  max resid 0.001127995 
    ## ... Similar to previous best
    ## Run 28 stress 0.1110873 
    ## Run 29 stress 0.0893854 
    ## ... Procrustes: rmse 0.01180558  max resid 0.03950583 
    ## Run 30 stress 0.08946329 
    ## ... Procrustes: rmse 0.03744984  max resid 0.1177017 
    ## Run 31 stress 0.0918088 
    ## Run 32 stress 0.09503421 
    ## Run 33 stress 0.08938539 
    ## ... Procrustes: rmse 0.01180179  max resid 0.03950239 
    ## Run 34 stress 0.1052648 
    ## Run 35 stress 0.09018707 
    ## Run 36 stress 0.08946329 
    ## ... Procrustes: rmse 0.03746891  max resid 0.1177213 
    ## Run 37 stress 0.09021167 
    ## Run 38 stress 0.09088601 
    ## Run 39 stress 0.1086566 
    ## Run 40 stress 0.08938548 
    ## ... Procrustes: rmse 0.01174279  max resid 0.03951169 
    ## Run 41 stress 0.1061297 
    ## Run 42 stress 0.1087577 
    ## Run 43 stress 0.105265 
    ## Run 44 stress 0.0893856 
    ## ... Procrustes: rmse 0.01173197  max resid 0.03942497 
    ## Run 45 stress 0.1076303 
    ## Run 46 stress 0.1101448 
    ## Run 47 stress 0.09130085 
    ## Run 48 stress 0.1064188 
    ## Run 49 stress 0.08938966 
    ## ... Procrustes: rmse 0.03592614  max resid 0.1182727 
    ## Run 50 stress 0.09612787 
    ## Run 51 stress 0.09180879 
    ## Run 52 stress 0.08938538 
    ## ... Procrustes: rmse 0.01179763  max resid 0.03948993 
    ## Run 53 stress 0.08938961 
    ## ... Procrustes: rmse 0.03591556  max resid 0.1182504 
    ## Run 54 stress 0.09592141 
    ## Run 55 stress 0.1087865 
    ## Run 56 stress 0.08926096 
    ## ... Procrustes: rmse 0.0005388887  max resid 0.001547078 
    ## ... Similar to previous best
    ## Run 57 stress 0.1074434 
    ## Run 58 stress 0.09044625 
    ## Run 59 stress 0.09088622 
    ## Run 60 stress 0.1067508 
    ## Run 61 stress 0.09503453 
    ## Run 62 stress 0.08946335 
    ## ... Procrustes: rmse 0.03743321  max resid 0.1176818 
    ## Run 63 stress 0.08926081 
    ## ... Procrustes: rmse 0.0002024777  max resid 0.0006444046 
    ## ... Similar to previous best
    ## Run 64 stress 0.1052649 
    ## Run 65 stress 0.0950342 
    ## Run 66 stress 0.08938966 
    ## ... Procrustes: rmse 0.03588355  max resid 0.1182007 
    ## Run 67 stress 0.1052648 
    ## Run 68 stress 0.09775571 
    ## Run 69 stress 0.1060432 
    ## Run 70 stress 0.0902117 
    ## Run 71 stress 0.0893856 
    ## ... Procrustes: rmse 0.01172185  max resid 0.03952983 
    ## Run 72 stress 0.09503443 
    ## Run 73 stress 0.08938541 
    ## ... Procrustes: rmse 0.01179952  max resid 0.03951449 
    ## Run 74 stress 0.09503444 
    ## Run 75 stress 0.1091449 
    ## Run 76 stress 0.1079023 
    ## Run 77 stress 0.1067586 
    ## Run 78 stress 0.08926083 
    ## ... Procrustes: rmse 0.0004541331  max resid 0.001308784 
    ## ... Similar to previous best
    ## Run 79 stress 0.1108323 
    ## Run 80 stress 0.09039134 
    ## Run 81 stress 0.08926115 
    ## ... Procrustes: rmse 0.0006339228  max resid 0.001841791 
    ## ... Similar to previous best
    ## Run 82 stress 0.08926107 
    ## ... Procrustes: rmse 0.0005876209  max resid 0.00176831 
    ## ... Similar to previous best
    ## Run 83 stress 0.1102604 
    ## Run 84 stress 0.1086568 
    ## Run 85 stress 0.09503417 
    ## Run 86 stress 0.09503413 
    ## Run 87 stress 0.08926106 
    ## ... Procrustes: rmse 0.0005691288  max resid 0.001644411 
    ## ... Similar to previous best
    ## Run 88 stress 0.1063186 
    ## Run 89 stress 0.08946661 
    ## ... Procrustes: rmse 0.03313861  max resid 0.1162437 
    ## Run 90 stress 0.1067501 
    ## Run 91 stress 0.09130091 
    ## Run 92 stress 0.1074433 
    ## Run 93 stress 0.09594332 
    ## Run 94 stress 0.08946334 
    ## ... Procrustes: rmse 0.03743235  max resid 0.1176739 
    ## Run 95 stress 0.08926082 
    ## ... Procrustes: rmse 0.0001528685  max resid 0.000495879 
    ## ... Similar to previous best
    ## Run 96 stress 0.1087571 
    ## Run 97 stress 0.1052647 
    ## Run 98 stress 0.08946331 
    ## ... Procrustes: rmse 0.03744248  max resid 0.1176895 
    ## Run 99 stress 0.0908735 
    ## Run 100 stress 0.08938967 
    ## ... Procrustes: rmse 0.0358807  max resid 0.1181942 
    ## Run 101 stress 0.1056907 
    ## Run 102 stress 0.1104186 
    ## Run 103 stress 0.106419 
    ## Run 104 stress 0.1092949 
    ## Run 105 stress 0.1105343 
    ## Run 106 stress 0.1063194 
    ## Run 107 stress 0.1052653 
    ## Run 108 stress 0.08938549 
    ## ... Procrustes: rmse 0.01173201  max resid 0.03937682 
    ## Run 109 stress 0.09109108 
    ## Run 110 stress 0.1079018 
    ## Run 111 stress 0.08926064 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002419526  max resid 0.0005843955 
    ## ... Similar to previous best
    ## Run 112 stress 0.08938538 
    ## ... Procrustes: rmse 0.01183745  max resid 0.03971478 
    ## Run 113 stress 0.1056892 
    ## Run 114 stress 0.0893854 
    ## ... Procrustes: rmse 0.01184534  max resid 0.03974092 
    ## Run 115 stress 0.08946651 
    ## ... Procrustes: rmse 0.03321868  max resid 0.1163497 
    ## Run 116 stress 0.1075421 
    ## Run 117 stress 0.08926094 
    ## ... Procrustes: rmse 0.0002906596  max resid 0.0009057823 
    ## ... Similar to previous best
    ## Run 118 stress 0.08926064 
    ## ... New best solution
    ## ... Procrustes: rmse 8.699467e-05  max resid 0.0002908103 
    ## ... Similar to previous best
    ## Run 119 stress 0.1096133 
    ## Run 120 stress 0.08946653 
    ## ... Procrustes: rmse 0.03320174  max resid 0.1163103 
    ## Run 121 stress 0.09018711 
    ## Run 122 stress 0.0908734 
    ## Run 123 stress 0.08938538 
    ## ... Procrustes: rmse 0.0118393  max resid 0.03971169 
    ## Run 124 stress 0.08926065 
    ## ... Procrustes: rmse 0.0001054817  max resid 0.0003564078 
    ## ... Similar to previous best
    ## Run 125 stress 0.09503442 
    ## Run 126 stress 0.1060428 
    ## Run 127 stress 0.09088615 
    ## Run 128 stress 0.1052649 
    ## Run 129 stress 0.08926063 
    ## ... New best solution
    ## ... Procrustes: rmse 2.728512e-05  max resid 7.983594e-05 
    ## ... Similar to previous best
    ## Run 130 stress 0.08938541 
    ## ... Procrustes: rmse 0.0118729  max resid 0.03965131 
    ## Run 131 stress 0.08926089 
    ## ... Procrustes: rmse 0.0003545075  max resid 0.001123118 
    ## ... Similar to previous best
    ## Run 132 stress 0.09592123 
    ## Run 133 stress 0.1076302 
    ## Run 134 stress 0.1087862 
    ## Run 135 stress 0.08946333 
    ## ... Procrustes: rmse 0.03748545  max resid 0.1177229 
    ## Run 136 stress 0.09088631 
    ## Run 137 stress 0.1065384 
    ## Run 138 stress 0.08951722 
    ## ... Procrustes: rmse 0.03507634  max resid 0.1156461 
    ## Run 139 stress 0.08938967 
    ## ... Procrustes: rmse 0.03596258  max resid 0.1183206 
    ## Run 140 stress 0.1096132 
    ## Run 141 stress 0.09503439 
    ## Run 142 stress 0.09592157 
    ## Run 143 stress 0.08926091 
    ## ... Procrustes: rmse 0.000374007  max resid 0.001186916 
    ## ... Similar to previous best
    ## Run 144 stress 0.09018716 
    ## Run 145 stress 0.08946333 
    ## ... Procrustes: rmse 0.0375363  max resid 0.1178068 
    ## Run 146 stress 0.1067585 
    ## Run 147 stress 0.08946328 
    ## ... Procrustes: rmse 0.03750463  max resid 0.1177583 
    ## Run 148 stress 0.1071322 
    ## Run 149 stress 0.1074432 
    ## Run 150 stress 0.08938963 
    ## ... Procrustes: rmse 0.03595867  max resid 0.1183002 
    ## Run 151 stress 0.1065238 
    ## Run 152 stress 0.1061299 
    ## Run 153 stress 0.1071321 
    ## Run 154 stress 0.1074431 
    ## Run 155 stress 0.0901871 
    ## Run 156 stress 0.1075421 
    ## Run 157 stress 0.09503433 
    ## Run 158 stress 0.08926063 
    ## ... Procrustes: rmse 6.147554e-05  max resid 0.0002159832 
    ## ... Similar to previous best
    ## Run 159 stress 0.1067884 
    ## Run 160 stress 0.09039136 
    ## Run 161 stress 0.08946667 
    ## ... Procrustes: rmse 0.0331707  max resid 0.1162673 
    ## Run 162 stress 0.08926089 
    ## ... Procrustes: rmse 0.0003593645  max resid 0.001135662 
    ## ... Similar to previous best
    ## Run 163 stress 0.0892608 
    ## ... Procrustes: rmse 0.0002685207  max resid 0.0009003435 
    ## ... Similar to previous best
    ## Run 164 stress 0.1063193 
    ## Run 165 stress 0.09087366 
    ## Run 166 stress 0.08946341 
    ## ... Procrustes: rmse 0.03746775  max resid 0.1176945 
    ## Run 167 stress 0.09503448 
    ## Run 168 stress 0.1108633 
    ## Run 169 stress 0.08938546 
    ## ... Procrustes: rmse 0.01189618  max resid 0.0396558 
    ## Run 170 stress 0.08938558 
    ## ... Procrustes: rmse 0.01189969  max resid 0.03961786 
    ## Run 171 stress 0.08938537 
    ## ... Procrustes: rmse 0.01184801  max resid 0.0396789 
    ## Run 172 stress 0.1071321 
    ## Run 173 stress 0.08926075 
    ## ... Procrustes: rmse 0.0002522042  max resid 0.0008116338 
    ## ... Similar to previous best
    ## Run 174 stress 0.09021166 
    ## Run 175 stress 0.1074431 
    ## Run 176 stress 0.1061303 
    ## Run 177 stress 0.1067894 
    ## Run 178 stress 0.08946654 
    ## ... Procrustes: rmse 0.03319123  max resid 0.1163001 
    ## Run 179 stress 0.08946333 
    ## ... Procrustes: rmse 0.03753835  max resid 0.1178084 
    ## Run 180 stress 0.09039132 
    ## Run 181 stress 0.1075423 
    ## Run 182 stress 0.09130087 
    ## Run 183 stress 0.08938544 
    ## ... Procrustes: rmse 0.01188915  max resid 0.03965459 
    ## Run 184 stress 0.08946656 
    ## ... Procrustes: rmse 0.03319279  max resid 0.1163008 
    ## Run 185 stress 0.08938539 
    ## ... Procrustes: rmse 0.0118566  max resid 0.03970887 
    ## Run 186 stress 0.1071322 
    ## Run 187 stress 0.1091883 
    ## Run 188 stress 0.08938553 
    ## ... Procrustes: rmse 0.01175382  max resid 0.03958649 
    ## Run 189 stress 0.08951718 
    ## ... Procrustes: rmse 0.03505879  max resid 0.1156599 
    ## Run 190 stress 0.08946331 
    ## ... Procrustes: rmse 0.03752898  max resid 0.1177931 
    ## Run 191 stress 0.1052649 
    ## Run 192 stress 0.08938538 
    ## ... Procrustes: rmse 0.01184542  max resid 0.03966262 
    ## Run 193 stress 0.1086562 
    ## Run 194 stress 0.08938539 
    ## ... Procrustes: rmse 0.01182529  max resid 0.03968243 
    ## Run 195 stress 0.0909959 
    ## Run 196 stress 0.08926077 
    ## ... Procrustes: rmse 0.0002713227  max resid 0.0008674497 
    ## ... Similar to previous best
    ## Run 197 stress 0.09044608 
    ## Run 198 stress 0.08926093 
    ## ... Procrustes: rmse 0.0003354428  max resid 0.0009833488 
    ## ... Similar to previous best
    ## Run 199 stress 0.09039133 
    ## Run 200 stress 0.09592156 
    ## Run 201 stress 0.08946335 
    ## ... Procrustes: rmse 0.03748124  max resid 0.1177107 
    ## Run 202 stress 0.1064193 
    ## Run 203 stress 0.09018715 
    ## Run 204 stress 0.1065385 
    ## Run 205 stress 0.08938962 
    ## ... Procrustes: rmse 0.03595967  max resid 0.1182985 
    ## Run 206 stress 0.08938963 
    ## ... Procrustes: rmse 0.03595997  max resid 0.1183034 
    ## Run 207 stress 0.1060092 
    ## Run 208 stress 0.1075422 
    ## Run 209 stress 0.1075423 
    ## Run 210 stress 0.106419 
    ## Run 211 stress 0.09503455 
    ## Run 212 stress 0.1056905 
    ## Run 213 stress 0.08938549 
    ## ... Procrustes: rmse 0.01177104  max resid 0.03961019 
    ## Run 214 stress 0.08926085 
    ## ... Procrustes: rmse 0.0003333923  max resid 0.001065958 
    ## ... Similar to previous best
    ## Run 215 stress 0.08946651 
    ## ... Procrustes: rmse 0.03322108  max resid 0.1163552 
    ## Run 216 stress 0.08946334 
    ## ... Procrustes: rmse 0.03748118  max resid 0.1177182 
    ## Run 217 stress 0.0908861 
    ## Run 218 stress 0.09594323 
    ## Run 219 stress 0.09021167 
    ## Run 220 stress 0.08938538 
    ## ... Procrustes: rmse 0.01185662  max resid 0.03969005 
    ## Run 221 stress 0.1056894 
    ## Run 222 stress 0.1092209 
    ## Run 223 stress 0.1091825 
    ## Run 224 stress 0.1103709 
    ## Run 225 stress 0.1092091 
    ## Run 226 stress 0.08938968 
    ## ... Procrustes: rmse 0.03591715  max resid 0.1182334 
    ## Run 227 stress 0.09044626 
    ## Run 228 stress 0.1064193 
    ## Run 229 stress 0.1087863 
    ## Run 230 stress 0.08938552 
    ## ... Procrustes: rmse 0.01178714  max resid 0.0397206 
    ## Run 231 stress 0.1076305 
    ## Run 232 stress 0.0903913 
    ## Run 233 stress 0.08946653 
    ## ... Procrustes: rmse 0.03320159  max resid 0.1163169 
    ## Run 234 stress 0.109256 
    ## Run 235 stress 0.08946329 
    ## ... Procrustes: rmse 0.03751044  max resid 0.1177628 
    ## Run 236 stress 0.1061297 
    ## Run 237 stress 0.1056902 
    ## Run 238 stress 0.08951718 
    ## ... Procrustes: rmse 0.03506613  max resid 0.1156741 
    ## Run 239 stress 0.1066194 
    ## Run 240 stress 0.1092561 
    ## Run 241 stress 0.08926088 
    ## ... Procrustes: rmse 0.0003513283  max resid 0.001128181 
    ## ... Similar to previous best
    ## Run 242 stress 0.09130088 
    ## Run 243 stress 0.08938539 
    ## ... Procrustes: rmse 0.01184423  max resid 0.03968094 
    ## Run 244 stress 0.08926079 
    ## ... Procrustes: rmse 0.0002882125  max resid 0.0009171273 
    ## ... Similar to previous best
    ## Run 245 stress 0.08938539 
    ## ... Procrustes: rmse 0.01187112  max resid 0.03968978 
    ## Run 246 stress 0.08938968 
    ## ... Procrustes: rmse 0.03597427  max resid 0.1183241 
    ## Run 247 stress 0.1075421 
    ## Run 248 stress 0.1052649 
    ## Run 249 stress 0.1096134 
    ## Run 250 stress 0.1088391 
    ## Run 251 stress 0.08938964 
    ## ... Procrustes: rmse 0.03592603  max resid 0.1182452 
    ## Run 252 stress 0.08938539 
    ## ... Procrustes: rmse 0.01185735  max resid 0.03969002 
    ## Run 253 stress 0.08938543 
    ## ... Procrustes: rmse 0.01180646  max resid 0.03970291 
    ## Run 254 stress 0.1103725 
    ## Run 255 stress 0.1092128 
    ## Run 256 stress 0.08938561 
    ## ... Procrustes: rmse 0.01177304  max resid 0.03973495 
    ## Run 257 stress 0.08946334 
    ## ... Procrustes: rmse 0.0374827  max resid 0.1177187 
    ## Run 258 stress 0.09178285 
    ## Run 259 stress 0.08938548 
    ## ... Procrustes: rmse 0.01190183  max resid 0.03965283 
    ## Run 260 stress 0.09180882 
    ## Run 261 stress 0.09087344 
    ## Run 262 stress 0.0908734 
    ## Run 263 stress 0.1067489 
    ## Run 264 stress 0.1052646 
    ## Run 265 stress 0.09018707 
    ## Run 266 stress 0.09039131 
    ## Run 267 stress 0.1099391 
    ## Run 268 stress 0.09503419 
    ## Run 269 stress 0.110534 
    ## Run 270 stress 0.09021169 
    ## Run 271 stress 0.09180901 
    ## Run 272 stress 0.1063188 
    ## Run 273 stress 0.08946653 
    ## ... Procrustes: rmse 0.03322459  max resid 0.1163711 
    ## Run 274 stress 0.0903913 
    ## Run 275 stress 0.09087341 
    ## Run 276 stress 0.1056904 
    ## Run 277 stress 0.1052647 
    ## Run 278 stress 0.08938542 
    ## ... Procrustes: rmse 0.01181452  max resid 0.03960859 
    ## Run 279 stress 0.09039133 
    ## Run 280 stress 0.1071322 
    ## Run 281 stress 0.1056895 
    ## Run 282 stress 0.08938537 
    ## ... Procrustes: rmse 0.01183955  max resid 0.03968251 
    ## Run 283 stress 0.08926108 
    ## ... Procrustes: rmse 0.0004645418  max resid 0.001524305 
    ## ... Similar to previous best
    ## Run 284 stress 0.1068767 
    ## Run 285 stress 0.1071321 
    ## Run 286 stress 0.09039132 
    ## Run 287 stress 0.09503455 
    ## Run 288 stress 0.08938554 
    ## ... Procrustes: rmse 0.01179746  max resid 0.03973822 
    ## Run 289 stress 0.09018715 
    ## Run 290 stress 0.1108641 
    ## Run 291 stress 0.09018707 
    ## Run 292 stress 0.08938965 
    ## ... Procrustes: rmse 0.03592344  max resid 0.1182442 
    ## Run 293 stress 0.1061303 
    ## Run 294 stress 0.08938544 
    ## ... Procrustes: rmse 0.01189049  max resid 0.03967318 
    ## Run 295 stress 0.08938543 
    ## ... Procrustes: rmse 0.01182501  max resid 0.03969529 
    ## Run 296 stress 0.08938542 
    ## ... Procrustes: rmse 0.01186471  max resid 0.03972557 
    ## Run 297 stress 0.106043 
    ## Run 298 stress 0.1052647 
    ## Run 299 stress 0.0950342 
    ## Run 300 stress 0.08926105 
    ## ... Procrustes: rmse 0.0004202039  max resid 0.001245901 
    ## ... Similar to previous best
    ## Run 301 stress 0.08926087 
    ## ... Procrustes: rmse 0.0003085331  max resid 0.0009710984 
    ## ... Similar to previous best
    ## Run 302 stress 0.09021168 
    ## Run 303 stress 0.1052649 
    ## Run 304 stress 0.09130093 
    ## Run 305 stress 0.1076306 
    ## Run 306 stress 0.08926099 
    ## ... Procrustes: rmse 0.0004189796  max resid 0.001301528 
    ## ... Similar to previous best
    ## Run 307 stress 0.0892607 
    ## ... Procrustes: rmse 0.0002058091  max resid 0.0006594383 
    ## ... Similar to previous best
    ## Run 308 stress 0.1062565 
    ## Run 309 stress 0.1071321 
    ## Run 310 stress 0.1087565 
    ## Run 311 stress 0.09039147 
    ## Run 312 stress 0.1052647 
    ## Run 313 stress 0.08938547 
    ## ... Procrustes: rmse 0.01177402  max resid 0.03960436 
    ## Run 314 stress 0.1092096 
    ## Run 315 stress 0.1052651 
    ## Run 316 stress 0.0892611 
    ## ... Procrustes: rmse 0.0004755489  max resid 0.001494252 
    ## ... Similar to previous best
    ## Run 317 stress 0.1056904 
    ## Run 318 stress 0.111137 
    ## Run 319 stress 0.1118955 
    ## Run 320 stress 0.0892608 
    ## ... Procrustes: rmse 0.0002852773  max resid 0.0009599053 
    ## ... Similar to previous best
    ## Run 321 stress 0.09087348 
    ## Run 322 stress 0.1071321 
    ## Run 323 stress 0.08946651 
    ## ... Procrustes: rmse 0.03321421  max resid 0.1163544 
    ## Run 324 stress 0.09228502 
    ## Run 325 stress 0.09178296 
    ## Run 326 stress 0.09178293 
    ## Run 327 stress 0.1052647 
    ## Run 328 stress 0.1067496 
    ## Run 329 stress 0.1067516 
    ## Run 330 stress 0.1091883 
    ## Run 331 stress 0.1061302 
    ## Run 332 stress 0.09228503 
    ## Run 333 stress 0.1096133 
    ## Run 334 stress 0.09039136 
    ## Run 335 stress 0.0901871 
    ## Run 336 stress 0.08938553 
    ## ... Procrustes: rmse 0.01188885  max resid 0.03960815 
    ## Run 337 stress 0.09130088 
    ## Run 338 stress 0.08926094 
    ## ... Procrustes: rmse 0.0003908719  max resid 0.001266645 
    ## ... Similar to previous best
    ## Run 339 stress 0.1063191 
    ## Run 340 stress 0.09503432 
    ## Run 341 stress 0.1056896 
    ## Run 342 stress 0.08938554 
    ## ... Procrustes: rmse 0.0118962  max resid 0.03962995 
    ## Run 343 stress 0.08926069 
    ## ... Procrustes: rmse 0.0001603578  max resid 0.0005078685 
    ## ... Similar to previous best
    ## Run 344 stress 0.1074432 
    ## Run 345 stress 0.09262364 
    ## Run 346 stress 0.1076302 
    ## Run 347 stress 0.1091877 
    ## Run 348 stress 0.08938538 
    ## ... Procrustes: rmse 0.01184109  max resid 0.03967524 
    ## Run 349 stress 0.09088604 
    ## Run 350 stress 0.10569 
    ## Run 351 stress 0.08938549 
    ## ... Procrustes: rmse 0.01178711  max resid 0.03964981 
    ## Run 352 stress 0.09180887 
    ## Run 353 stress 0.1076303 
    ## Run 354 stress 0.08938961 
    ## ... Procrustes: rmse 0.03595013  max resid 0.118285 
    ## Run 355 stress 0.1074433 
    ## Run 356 stress 0.08938543 
    ## ... Procrustes: rmse 0.01186502  max resid 0.03971379 
    ## Run 357 stress 0.1088389 
    ## Run 358 stress 0.09503431 
    ## Run 359 stress 0.08938541 
    ## ... Procrustes: rmse 0.01179929  max resid 0.03962969 
    ## Run 360 stress 0.1067376 
    ## Run 361 stress 0.09109109 
    ## Run 362 stress 0.08946328 
    ## ... Procrustes: rmse 0.03750762  max resid 0.1177584 
    ## Run 363 stress 0.09018707 
    ## Run 364 stress 0.09503448 
    ## Run 365 stress 0.09592167 
    ## Run 366 stress 0.08938537 
    ## ... Procrustes: rmse 0.01183767  max resid 0.03967006 
    ## Run 367 stress 0.09039133 
    ## Run 368 stress 0.108837 
    ## Run 369 stress 0.09018712 
    ## Run 370 stress 0.1067586 
    ## Run 371 stress 0.09021166 
    ## Run 372 stress 0.09262375 
    ## Run 373 stress 0.09018708 
    ## Run 374 stress 0.08938968 
    ## ... Procrustes: rmse 0.03594763  max resid 0.1182849 
    ## Run 375 stress 0.1086565 
    ## Run 376 stress 0.1071322 
    ## Run 377 stress 0.08926064 
    ## ... Procrustes: rmse 0.0001090326  max resid 0.0003563163 
    ## ... Similar to previous best
    ## Run 378 stress 0.1052648 
    ## Run 379 stress 0.09099539 
    ## Run 380 stress 0.0908734 
    ## Run 381 stress 0.10569 
    ## Run 382 stress 0.111137 
    ## Run 383 stress 0.08926105 
    ## ... Procrustes: rmse 0.0004528549  max resid 0.001448149 
    ## ... Similar to previous best
    ## Run 384 stress 0.1085127 
    ## Run 385 stress 0.1087862 
    ## Run 386 stress 0.10569 
    ## Run 387 stress 0.09018708 
    ## Run 388 stress 0.1071321 
    ## Run 389 stress 0.09503433 
    ## Run 390 stress 0.08926118 
    ## ... Procrustes: rmse 0.0005096068  max resid 0.001599616 
    ## ... Similar to previous best
    ## Run 391 stress 0.1074432 
    ## Run 392 stress 0.1092128 
    ## Run 393 stress 0.1065335 
    ## Run 394 stress 0.0893854 
    ## ... Procrustes: rmse 0.01185492  max resid 0.03970997 
    ## Run 395 stress 0.08938963 
    ## ... Procrustes: rmse 0.03593447  max resid 0.1182636 
    ## Run 396 stress 0.08946656 
    ## ... Procrustes: rmse 0.03319009  max resid 0.116299 
    ## Run 397 stress 0.1086561 
    ## Run 398 stress 0.08951722 
    ## ... Procrustes: rmse 0.03503029  max resid 0.1156069 
    ## Run 399 stress 0.08926063 
    ## ... Procrustes: rmse 2.098344e-05  max resid 5.894153e-05 
    ## ... Similar to previous best
    ## Run 400 stress 0.08938546 
    ## ... Procrustes: rmse 0.01187454  max resid 0.03973716 
    ## Run 401 stress 0.1060431 
    ## Run 402 stress 0.09088602 
    ## Run 403 stress 0.09130095 
    ## Run 404 stress 0.1092091 
    ## Run 405 stress 0.1074433 
    ## Run 406 stress 0.1062565 
    ## Run 407 stress 0.09039129 
    ## Run 408 stress 0.08926065 
    ## ... Procrustes: rmse 4.30419e-05  max resid 0.000121661 
    ## ... Similar to previous best
    ## Run 409 stress 0.08926074 
    ## ... Procrustes: rmse 0.0001719079  max resid 0.000476905 
    ## ... Similar to previous best
    ## Run 410 stress 0.08926071 
    ## ... Procrustes: rmse 0.0001348722  max resid 0.0003730048 
    ## ... Similar to previous best
    ## Run 411 stress 0.1103725 
    ## Run 412 stress 0.1075422 
    ## Run 413 stress 0.09503436 
    ## Run 414 stress 0.08938541 
    ## ... Procrustes: rmse 0.01184093  max resid 0.0396264 
    ## Run 415 stress 0.09130093 
    ## Run 416 stress 0.09130102 
    ## Run 417 stress 0.08946657 
    ## ... Procrustes: rmse 0.03318651  max resid 0.1162944 
    ## Run 418 stress 0.1087586 
    ## Run 419 stress 0.10569 
    ## Run 420 stress 0.09178283 
    ## Run 421 stress 0.08946653 
    ## ... Procrustes: rmse 0.03322898  max resid 0.1163822 
    ## Run 422 stress 0.09087356 
    ## Run 423 stress 0.08946333 
    ## ... Procrustes: rmse 0.03748349  max resid 0.1177211 
    ## Run 424 stress 0.1056909 
    ## Run 425 stress 0.1063197 
    ## Run 426 stress 0.08946664 
    ## ... Procrustes: rmse 0.03317616  max resid 0.116284 
    ## Run 427 stress 0.08926072 
    ## ... Procrustes: rmse 0.000215475  max resid 0.0006822573 
    ## ... Similar to previous best
    ## Run 428 stress 0.08938539 
    ## ... Procrustes: rmse 0.0118645  max resid 0.03967596 
    ## Run 429 stress 0.08926065 
    ## ... Procrustes: rmse 0.0001211542  max resid 0.0003863825 
    ## ... Similar to previous best
    ## Run 430 stress 0.08926089 
    ## ... Procrustes: rmse 0.0003635461  max resid 0.001144969 
    ## ... Similar to previous best
    ## Run 431 stress 0.08938562 
    ## ... Procrustes: rmse 0.01192964  max resid 0.03963239 
    ## Run 432 stress 0.09503439 
    ## Run 433 stress 0.09021166 
    ## Run 434 stress 0.1092127 
    ## Run 435 stress 0.09044614 
    ## Run 436 stress 0.1091878 
    ## Run 437 stress 0.08946329 
    ## ... Procrustes: rmse 0.03750723  max resid 0.1177606 
    ## Run 438 stress 0.1076305 
    ## Run 439 stress 0.08938542 
    ## ... Procrustes: rmse 0.01187044  max resid 0.03964671 
    ## Run 440 stress 0.09088604 
    ## Run 441 stress 0.0908734 
    ## Run 442 stress 0.1096136 
    ## Run 443 stress 0.08938965 
    ## ... Procrustes: rmse 0.03596116  max resid 0.1183046 
    ## Run 444 stress 0.09262371 
    ## Run 445 stress 0.08926069 
    ## ... Procrustes: rmse 8.057338e-05  max resid 0.0001775616 
    ## ... Similar to previous best
    ## Run 446 stress 0.09178301 
    ## Run 447 stress 0.08938549 
    ## ... Procrustes: rmse 0.01179226  max resid 0.03970141 
    ## Run 448 stress 0.0892607 
    ## ... Procrustes: rmse 0.0001441876  max resid 0.0003932602 
    ## ... Similar to previous best
    ## Run 449 stress 0.1056892 
    ## Run 450 stress 0.09044617 
    ## Run 451 stress 0.1063195 
    ## Run 452 stress 0.09262377 
    ## Run 453 stress 0.1075423 
    ## Run 454 stress 0.08938561 
    ## ... Procrustes: rmse 0.01174726  max resid 0.03958877 
    ## Run 455 stress 0.1067898 
    ## Run 456 stress 0.09021167 
    ## Run 457 stress 0.1091827 
    ## Run 458 stress 0.1087862 
    ## Run 459 stress 0.0950342 
    ## Run 460 stress 0.1110409 
    ## Run 461 stress 0.1091825 
    ## Run 462 stress 0.08938547 
    ## ... Procrustes: rmse 0.01182619  max resid 0.03967502 
    ## Run 463 stress 0.08938539 
    ## ... Procrustes: rmse 0.01183629  max resid 0.03970146 
    ## Run 464 stress 0.1087863 
    ## Run 465 stress 0.1092126 
    ## Run 466 stress 0.09087344 
    ## Run 467 stress 0.09039135 
    ## Run 468 stress 0.09130088 
    ## Run 469 stress 0.0892609 
    ## ... Procrustes: rmse 0.0003680561  max resid 0.001199219 
    ## ... Similar to previous best
    ## Run 470 stress 0.08938978 
    ## ... Procrustes: rmse 0.03590317  max resid 0.1182094 
    ## Run 471 stress 0.1087864 
    ## Run 472 stress 0.108839 
    ## Run 473 stress 0.1074434 
    ## Run 474 stress 0.0903913 
    ## Run 475 stress 0.09592184 
    ## Run 476 stress 0.08946653 
    ## ... Procrustes: rmse 0.0331971  max resid 0.1163109 
    ## Run 477 stress 0.09039134 
    ## Run 478 stress 0.1076304 
    ## Run 479 stress 0.09099547 
    ## Run 480 stress 0.1092128 
    ## Run 481 stress 0.08926088 
    ## ... Procrustes: rmse 0.0003528685  max resid 0.001134954 
    ## ... Similar to previous best
    ## Run 482 stress 0.0903913 
    ## Run 483 stress 0.09088603 
    ## Run 484 stress 0.08946653 
    ## ... Procrustes: rmse 0.03319493  max resid 0.116303 
    ## Run 485 stress 0.1076303 
    ## Run 486 stress 0.1065367 
    ## Run 487 stress 0.08946655 
    ## ... Procrustes: rmse 0.03322966  max resid 0.1163883 
    ## Run 488 stress 0.09612796 
    ## Run 489 stress 0.09503417 
    ## Run 490 stress 0.09018707 
    ## Run 491 stress 0.08938543 
    ## ... Procrustes: rmse 0.01186641  max resid 0.03972584 
    ## Run 492 stress 0.09021167 
    ## Run 493 stress 0.0894665 
    ## ... Procrustes: rmse 0.03321664  max resid 0.116344 
    ## Run 494 stress 0.1087571 
    ## Run 495 stress 0.08926066 
    ## ... Procrustes: rmse 5.149269e-05  max resid 0.0001587693 
    ## ... Similar to previous best
    ## Run 496 stress 0.08951717 
    ## ... Procrustes: rmse 0.03504314  max resid 0.1156212 
    ## Run 497 stress 0.1075423 
    ## Run 498 stress 0.1104473 
    ## Run 499 stress 0.1065377 
    ## Run 500 stress 0.1066194 
    ## *** Best solution repeated 36 times

``` r
# Mixed and stratified lakes
SD_beta_geo_MS_NMDS <- metaMDS(SD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.08503614 
    ## Run 1 stress 0.08440255 
    ## ... New best solution
    ## ... Procrustes: rmse 0.01442994  max resid 0.04297109 
    ## Run 2 stress 0.08503532 
    ## Run 3 stress 0.08773471 
    ## Run 4 stress 0.1006346 
    ## Run 5 stress 0.085035 
    ## Run 6 stress 0.09539197 
    ## Run 7 stress 0.0996997 
    ## Run 8 stress 0.08503496 
    ## Run 9 stress 0.08973862 
    ## Run 10 stress 0.1038045 
    ## Run 11 stress 0.09669871 
    ## Run 12 stress 0.08440263 
    ## ... Procrustes: rmse 6.02631e-05  max resid 0.000117583 
    ## ... Similar to previous best
    ## Run 13 stress 0.08973892 
    ## Run 14 stress 0.08503465 
    ## Run 15 stress 0.09159086 
    ## Run 16 stress 0.09337292 
    ## Run 17 stress 0.0926836 
    ## Run 18 stress 0.0930895 
    ## Run 19 stress 0.09145321 
    ## Run 20 stress 0.09030399 
    ## Run 21 stress 0.09030395 
    ## Run 22 stress 0.08440256 
    ## ... Procrustes: rmse 1.404458e-05  max resid 3.133671e-05 
    ## ... Similar to previous best
    ## Run 23 stress 0.09030394 
    ## Run 24 stress 0.09465905 
    ## Run 25 stress 0.09407973 
    ## Run 26 stress 0.08773475 
    ## Run 27 stress 0.09400516 
    ## Run 28 stress 0.09669878 
    ## Run 29 stress 0.09308967 
    ## Run 30 stress 0.08503471 
    ## Run 31 stress 0.09465909 
    ## Run 32 stress 0.09159083 
    ## Run 33 stress 0.09380594 
    ## Run 34 stress 0.09408017 
    ## Run 35 stress 0.09030396 
    ## Run 36 stress 0.09145322 
    ## Run 37 stress 0.08440268 
    ## ... Procrustes: rmse 0.0001218379  max resid 0.0002365011 
    ## ... Similar to previous best
    ## Run 38 stress 0.09969962 
    ## Run 39 stress 0.09407987 
    ## Run 40 stress 0.103804 
    ## Run 41 stress 0.0903041 
    ## Run 42 stress 0.09337281 
    ## Run 43 stress 0.09465915 
    ## Run 44 stress 0.09030393 
    ## Run 45 stress 0.08503475 
    ## Run 46 stress 0.09268387 
    ## Run 47 stress 0.09145322 
    ## Run 48 stress 0.103805 
    ## Run 49 stress 0.100461 
    ## Run 50 stress 0.08503509 
    ## Run 51 stress 0.09030397 
    ## Run 52 stress 0.0897387 
    ## Run 53 stress 0.09159086 
    ## Run 54 stress 0.09030403 
    ## Run 55 stress 0.09159084 
    ## Run 56 stress 0.09268321 
    ## Run 57 stress 0.09407965 
    ## Run 58 stress 0.09030399 
    ## Run 59 stress 0.09416134 
    ## Run 60 stress 0.08503477 
    ## Run 61 stress 0.08773493 
    ## Run 62 stress 0.09268328 
    ## Run 63 stress 0.09286084 
    ## Run 64 stress 0.09145327 
    ## Run 65 stress 0.09159085 
    ## Run 66 stress 0.1038037 
    ## Run 67 stress 0.08503563 
    ## Run 68 stress 0.09465916 
    ## Run 69 stress 0.09374352 
    ## Run 70 stress 0.09380585 
    ## Run 71 stress 0.0850353 
    ## Run 72 stress 0.09268409 
    ## Run 73 stress 0.09030399 
    ## Run 74 stress 0.09030395 
    ## Run 75 stress 0.0844026 
    ## ... Procrustes: rmse 0.0003010249  max resid 0.0005840153 
    ## ... Similar to previous best
    ## Run 76 stress 0.09337235 
    ## Run 77 stress 0.09407985 
    ## Run 78 stress 0.09030393 
    ## Run 79 stress 0.09159095 
    ## Run 80 stress 0.09408006 
    ## Run 81 stress 0.09408 
    ## Run 82 stress 0.08973873 
    ## Run 83 stress 0.09465906 
    ## Run 84 stress 0.08440251 
    ## ... New best solution
    ## ... Procrustes: rmse 0.000151086  max resid 0.0002838115 
    ## ... Similar to previous best
    ## Run 85 stress 0.09159084 
    ## Run 86 stress 0.09337292 
    ## Run 87 stress 0.08773482 
    ## Run 88 stress 0.09408 
    ## Run 89 stress 0.09445461 
    ## Run 90 stress 0.08503489 
    ## Run 91 stress 0.09407983 
    ## Run 92 stress 0.08503458 
    ## Run 93 stress 0.09760825 
    ## Run 94 stress 0.08973864 
    ## Run 95 stress 0.09407999 
    ## Run 96 stress 0.09407977 
    ## Run 97 stress 0.09030393 
    ## Run 98 stress 0.0928609 
    ## Run 99 stress 0.08503471 
    ## Run 100 stress 0.09030396 
    ## Run 101 stress 0.08773483 
    ## Run 102 stress 0.09407966 
    ## Run 103 stress 0.09407973 
    ## Run 104 stress 0.09168934 
    ## Run 105 stress 0.09030393 
    ## Run 106 stress 0.09030406 
    ## Run 107 stress 0.1052851 
    ## Run 108 stress 0.09308954 
    ## Run 109 stress 0.09969982 
    ## Run 110 stress 0.09159095 
    ## Run 111 stress 0.09159084 
    ## Run 112 stress 0.08440256 
    ## ... Procrustes: rmse 0.0001504052  max resid 0.00031062 
    ## ... Similar to previous best
    ## Run 113 stress 0.09760831 
    ## Run 114 stress 0.09416132 
    ## Run 115 stress 0.09404379 
    ## Run 116 stress 0.09030406 
    ## Run 117 stress 0.09535464 
    ## Run 118 stress 0.08773478 
    ## Run 119 stress 0.09337307 
    ## Run 120 stress 0.1038049 
    ## Run 121 stress 0.09030399 
    ## Run 122 stress 0.09337261 
    ## Run 123 stress 0.09308949 
    ## Run 124 stress 0.09374362 
    ## Run 125 stress 0.090304 
    ## Run 126 stress 0.09969982 
    ## Run 127 stress 0.09969964 
    ## Run 128 stress 0.09539203 
    ## Run 129 stress 0.08773469 
    ## Run 130 stress 0.09030399 
    ## Run 131 stress 0.0928609 
    ## Run 132 stress 0.09407965 
    ## Run 133 stress 0.08440264 
    ## ... Procrustes: rmse 0.0002334752  max resid 0.0004265662 
    ## ... Similar to previous best
    ## Run 134 stress 0.09447299 
    ## Run 135 stress 0.08503591 
    ## Run 136 stress 0.09268393 
    ## Run 137 stress 0.0933726 
    ## Run 138 stress 0.09337244 
    ## Run 139 stress 0.09416135 
    ## Run 140 stress 0.09381851 
    ## Run 141 stress 0.0897387 
    ## Run 142 stress 0.08440261 
    ## ... Procrustes: rmse 0.0002052094  max resid 0.0003784446 
    ## ... Similar to previous best
    ## Run 143 stress 0.09337263 
    ## Run 144 stress 0.09268327 
    ## Run 145 stress 0.09407973 
    ## Run 146 stress 0.08973877 
    ## Run 147 stress 0.093373 
    ## Run 148 stress 0.09416136 
    ## Run 149 stress 0.09465908 
    ## Run 150 stress 0.0903041 
    ## Run 151 stress 0.09286087 
    ## Run 152 stress 0.09030393 
    ## Run 153 stress 0.09030406 
    ## Run 154 stress 0.09145325 
    ## Run 155 stress 0.09168945 
    ## Run 156 stress 0.09337302 
    ## Run 157 stress 0.08773485 
    ## Run 158 stress 0.0930898 
    ## Run 159 stress 0.09408006 
    ## Run 160 stress 0.0877348 
    ## Run 161 stress 0.08440265 
    ## ... Procrustes: rmse 0.0002425241  max resid 0.0004407606 
    ## ... Similar to previous best
    ## Run 162 stress 0.08503496 
    ## Run 163 stress 0.1026215 
    ## Run 164 stress 0.09268404 
    ## Run 165 stress 0.09416139 
    ## Run 166 stress 0.09030398 
    ## Run 167 stress 0.09464468 
    ## Run 168 stress 0.09308979 
    ## Run 169 stress 0.08503636 
    ## Run 170 stress 0.1038049 
    ## Run 171 stress 0.09404374 
    ## Run 172 stress 0.09030397 
    ## Run 173 stress 0.1038041 
    ## Run 174 stress 0.08973867 
    ## Run 175 stress 0.09465904 
    ## Run 176 stress 0.09407977 
    ## Run 177 stress 0.09286084 
    ## Run 178 stress 0.08503539 
    ## Run 179 stress 0.08503491 
    ## Run 180 stress 0.09445488 
    ## Run 181 stress 0.09465905 
    ## Run 182 stress 0.09465908 
    ## Run 183 stress 0.0933724 
    ## Run 184 stress 0.09407967 
    ## Run 185 stress 0.09535468 
    ## Run 186 stress 0.09030393 
    ## Run 187 stress 0.09374336 
    ## Run 188 stress 0.08503637 
    ## Run 189 stress 0.100461 
    ## Run 190 stress 0.09159085 
    ## Run 191 stress 0.09721199 
    ## Run 192 stress 0.1038038 
    ## Run 193 stress 0.09286088 
    ## Run 194 stress 0.09416137 
    ## Run 195 stress 0.09407473 
    ## Run 196 stress 0.08503495 
    ## Run 197 stress 0.09030402 
    ## Run 198 stress 0.09030397 
    ## Run 199 stress 0.08503677 
    ## Run 200 stress 0.09168926 
    ## Run 201 stress 0.09168926 
    ## Run 202 stress 0.09416152 
    ## Run 203 stress 0.09168933 
    ## Run 204 stress 0.09308965 
    ## Run 205 stress 0.08503471 
    ## Run 206 stress 0.09030401 
    ## Run 207 stress 0.09030399 
    ## Run 208 stress 0.09286086 
    ## Run 209 stress 0.08440254 
    ## ... Procrustes: rmse 0.0001294097  max resid 0.0002376396 
    ## ... Similar to previous best
    ## Run 210 stress 0.09400548 
    ## Run 211 stress 0.09380582 
    ## Run 212 stress 0.09374255 
    ## Run 213 stress 0.09145324 
    ## Run 214 stress 0.0844027 
    ## ... Procrustes: rmse 0.0002847755  max resid 0.0005029691 
    ## ... Similar to previous best
    ## Run 215 stress 0.08773473 
    ## Run 216 stress 0.09416139 
    ## Run 217 stress 0.09721198 
    ## Run 218 stress 0.09030407 
    ## Run 219 stress 0.2567612 
    ## Run 220 stress 0.1038043 
    ## Run 221 stress 0.09590665 
    ## Run 222 stress 0.09337282 
    ## Run 223 stress 0.08440261 
    ## ... Procrustes: rmse 0.0002120252  max resid 0.0003926273 
    ## ... Similar to previous best
    ## Run 224 stress 0.1054521 
    ## Run 225 stress 0.09760814 
    ## Run 226 stress 0.08973884 
    ## Run 227 stress 0.09337246 
    ## Run 228 stress 0.09408006 
    ## Run 229 stress 0.1026216 
    ## Run 230 stress 0.08773469 
    ## Run 231 stress 0.09407996 
    ## Run 232 stress 0.09407454 
    ## Run 233 stress 0.2456114 
    ## Run 234 stress 0.0940054 
    ## Run 235 stress 0.09969965 
    ## Run 236 stress 0.08773466 
    ## Run 237 stress 0.09465905 
    ## Run 238 stress 0.09445464 
    ## Run 239 stress 0.09969987 
    ## Run 240 stress 0.09400517 
    ## Run 241 stress 0.09168928 
    ## Run 242 stress 0.09159084 
    ## Run 243 stress 0.08973869 
    ## Run 244 stress 0.09168945 
    ## Run 245 stress 0.09381845 
    ## Run 246 stress 0.09308981 
    ## Run 247 stress 0.09969981 
    ## Run 248 stress 0.0850348 
    ## Run 249 stress 0.09590659 
    ## Run 250 stress 0.09374243 
    ## Run 251 stress 0.08503624 
    ## Run 252 stress 0.08503514 
    ## Run 253 stress 0.08973864 
    ## Run 254 stress 0.09408 
    ## Run 255 stress 0.09407981 
    ## Run 256 stress 0.08973879 
    ## Run 257 stress 0.09268377 
    ## Run 258 stress 0.09465915 
    ## Run 259 stress 0.09030415 
    ## Run 260 stress 0.09400538 
    ## Run 261 stress 0.09145326 
    ## Run 262 stress 0.08773465 
    ## Run 263 stress 0.09030408 
    ## Run 264 stress 0.08503501 
    ## Run 265 stress 0.09308946 
    ## Run 266 stress 0.09381852 
    ## Run 267 stress 0.09464465 
    ## Run 268 stress 0.09159097 
    ## Run 269 stress 0.08973863 
    ## Run 270 stress 0.09374264 
    ## Run 271 stress 0.3025013 
    ## Run 272 stress 0.1026216 
    ## Run 273 stress 0.09030405 
    ## Run 274 stress 0.09030395 
    ## Run 275 stress 0.09145321 
    ## Run 276 stress 0.09145325 
    ## Run 277 stress 0.09465904 
    ## Run 278 stress 0.09159084 
    ## Run 279 stress 0.09268354 
    ## Run 280 stress 0.0914533 
    ## Run 281 stress 0.09145322 
    ## Run 282 stress 0.105452 
    ## Run 283 stress 0.09030394 
    ## Run 284 stress 0.09159105 
    ## Run 285 stress 0.09268386 
    ## Run 286 stress 0.08973872 
    ## Run 287 stress 0.0877349 
    ## Run 288 stress 0.09969985 
    ## Run 289 stress 0.09145327 
    ## Run 290 stress 0.09268369 
    ## Run 291 stress 0.09374286 
    ## Run 292 stress 0.09407969 
    ## Run 293 stress 0.08973881 
    ## Run 294 stress 0.09407115 
    ## Run 295 stress 0.1097012 
    ## Run 296 stress 0.08503516 
    ## Run 297 stress 0.09721197 
    ## Run 298 stress 0.08440264 
    ## ... Procrustes: rmse 0.0002406615  max resid 0.0004335314 
    ## ... Similar to previous best
    ## Run 299 stress 0.1026216 
    ## Run 300 stress 0.09407965 
    ## Run 301 stress 0.1038042 
    ## Run 302 stress 0.08973887 
    ## Run 303 stress 0.09407993 
    ## Run 304 stress 0.09337291 
    ## Run 305 stress 0.09407996 
    ## Run 306 stress 0.09145329 
    ## Run 307 stress 0.3146835 
    ## Run 308 stress 0.08503656 
    ## Run 309 stress 0.08503592 
    ## Run 310 stress 0.09721202 
    ## Run 311 stress 0.09465913 
    ## Run 312 stress 0.09145323 
    ## Run 313 stress 0.0897389 
    ## Run 314 stress 0.09286091 
    ## Run 315 stress 0.09374256 
    ## Run 316 stress 0.08503645 
    ## Run 317 stress 0.08773466 
    ## Run 318 stress 0.08440272 
    ## ... Procrustes: rmse 0.00029052  max resid 0.0005628956 
    ## ... Similar to previous best
    ## Run 319 stress 0.097212 
    ## Run 320 stress 0.09030403 
    ## Run 321 stress 0.09030398 
    ## Run 322 stress 0.09465911 
    ## Run 323 stress 0.08503487 
    ## Run 324 stress 0.08503477 
    ## Run 325 stress 0.09337273 
    ## Run 326 stress 0.09159085 
    ## Run 327 stress 0.09374253 
    ## Run 328 stress 0.08973891 
    ## Run 329 stress 0.09159084 
    ## Run 330 stress 0.09969971 
    ## Run 331 stress 0.08503513 
    ## Run 332 stress 0.0959064 
    ## Run 333 stress 0.09168934 
    ## Run 334 stress 0.08440255 
    ## ... Procrustes: rmse 0.0001247972  max resid 0.0002411524 
    ## ... Similar to previous best
    ## Run 335 stress 0.09286086 
    ## Run 336 stress 0.09374307 
    ## Run 337 stress 0.09407987 
    ## Run 338 stress 0.08440254 
    ## ... Procrustes: rmse 0.0001289801  max resid 0.0002370591 
    ## ... Similar to previous best
    ## Run 339 stress 0.0937436 
    ## Run 340 stress 0.08973864 
    ## Run 341 stress 0.09030404 
    ## Run 342 stress 0.09407989 
    ## Run 343 stress 0.09159092 
    ## Run 344 stress 0.09030422 
    ## Run 345 stress 0.09268342 
    ## Run 346 stress 0.1072863 
    ## Run 347 stress 0.09286083 
    ## Run 348 stress 0.09168948 
    ## Run 349 stress 0.09416134 
    ## Run 350 stress 0.08503482 
    ## Run 351 stress 0.09407108 
    ## Run 352 stress 0.09416133 
    ## Run 353 stress 0.08503501 
    ## Run 354 stress 0.09337289 
    ## Run 355 stress 0.08773475 
    ## Run 356 stress 0.0976082 
    ## Run 357 stress 0.09669877 
    ## Run 358 stress 0.09969989 
    ## Run 359 stress 0.08503657 
    ## Run 360 stress 0.09145322 
    ## Run 361 stress 0.08773478 
    ## Run 362 stress 0.09381851 
    ## Run 363 stress 0.09465905 
    ## Run 364 stress 0.09407987 
    ## Run 365 stress 0.09590557 
    ## Run 366 stress 0.09374375 
    ## Run 367 stress 0.08973862 
    ## Run 368 stress 0.09268321 
    ## Run 369 stress 0.08440269 
    ## ... Procrustes: rmse 0.0002310427  max resid 0.0004736122 
    ## ... Similar to previous best
    ## Run 370 stress 0.1038086 
    ## Run 371 stress 0.09590941 
    ## Run 372 stress 0.09145321 
    ## Run 373 stress 0.08773466 
    ## Run 374 stress 0.09159088 
    ## Run 375 stress 0.09465904 
    ## Run 376 stress 0.09407971 
    ## Run 377 stress 0.0877347 
    ## Run 378 stress 0.09539197 
    ## Run 379 stress 0.09159085 
    ## Run 380 stress 0.0916893 
    ## Run 381 stress 0.08973862 
    ## Run 382 stress 0.09321739 
    ## Run 383 stress 0.08503648 
    ## Run 384 stress 0.08503488 
    ## Run 385 stress 0.09337259 
    ## Run 386 stress 0.09030401 
    ## Run 387 stress 0.08773482 
    ## Run 388 stress 0.09374331 
    ## Run 389 stress 0.0946591 
    ## Run 390 stress 0.0940799 
    ## Run 391 stress 0.09408006 
    ## Run 392 stress 0.09168942 
    ## Run 393 stress 0.0914533 
    ## Run 394 stress 0.08440253 
    ## ... Procrustes: rmse 0.0001202528  max resid 0.0002272184 
    ## ... Similar to previous best
    ## Run 395 stress 0.09030419 
    ## Run 396 stress 0.09465912 
    ## Run 397 stress 0.09268381 
    ## Run 398 stress 0.09381845 
    ## Run 399 stress 0.090304 
    ## Run 400 stress 0.09159088 
    ## Run 401 stress 0.09464485 
    ## Run 402 stress 0.09268357 
    ## Run 403 stress 0.09159083 
    ## Run 404 stress 0.08440258 
    ## ... Procrustes: rmse 0.0001849193  max resid 0.0003387183 
    ## ... Similar to previous best
    ## Run 405 stress 0.08973866 
    ## Run 406 stress 0.0877347 
    ## Run 407 stress 0.09030397 
    ## Run 408 stress 0.08503519 
    ## Run 409 stress 0.09030393 
    ## Run 410 stress 0.2574703 
    ## Run 411 stress 0.09030405 
    ## Run 412 stress 0.09400537 
    ## Run 413 stress 0.08440262 
    ## ... Procrustes: rmse 0.0002211388  max resid 0.000401983 
    ## ... Similar to previous best
    ## Run 414 stress 0.09168939 
    ## Run 415 stress 0.09374349 
    ## Run 416 stress 0.09145324 
    ## Run 417 stress 0.08440264 
    ## ... Procrustes: rmse 0.0002163359  max resid 0.0003905239 
    ## ... Similar to previous best
    ## Run 418 stress 0.09145322 
    ## Run 419 stress 0.09145337 
    ## Run 420 stress 0.09760823 
    ## Run 421 stress 0.09268373 
    ## Run 422 stress 0.09374361 
    ## Run 423 stress 0.08503499 
    ## Run 424 stress 0.08973892 
    ## Run 425 stress 0.1026215 
    ## Run 426 stress 0.09760834 
    ## Run 427 stress 0.09969967 
    ## Run 428 stress 0.09159087 
    ## Run 429 stress 0.09447307 
    ## Run 430 stress 0.09464422 
    ## Run 431 stress 0.09417787 
    ## Run 432 stress 0.09286084 
    ## Run 433 stress 0.09374303 
    ## Run 434 stress 0.09159084 
    ## Run 435 stress 0.09465908 
    ## Run 436 stress 0.08503488 
    ## Run 437 stress 0.0941615 
    ## Run 438 stress 0.08973871 
    ## Run 439 stress 0.08773467 
    ## Run 440 stress 0.08973872 
    ## Run 441 stress 0.09145325 
    ## Run 442 stress 0.09416152 
    ## Run 443 stress 0.09447312 
    ## Run 444 stress 0.1026216 
    ## Run 445 stress 0.09417789 
    ## Run 446 stress 0.09268368 
    ## Run 447 stress 0.2774276 
    ## Run 448 stress 0.09407971 
    ## Run 449 stress 0.09159093 
    ## Run 450 stress 0.08440259 
    ## ... Procrustes: rmse 0.0001515666  max resid 0.0003176895 
    ## ... Similar to previous best
    ## Run 451 stress 0.09030393 
    ## Run 452 stress 0.09145321 
    ## Run 453 stress 0.09268401 
    ## Run 454 stress 0.08973882 
    ## Run 455 stress 0.09030398 
    ## Run 456 stress 0.09408011 
    ## Run 457 stress 0.08773466 
    ## Run 458 stress 0.09030393 
    ## Run 459 stress 0.09408023 
    ## Run 460 stress 0.08440252 
    ## ... Procrustes: rmse 4.053569e-05  max resid 7.444291e-05 
    ## ... Similar to previous best
    ## Run 461 stress 0.08503464 
    ## Run 462 stress 0.09539199 
    ## Run 463 stress 0.08973862 
    ## Run 464 stress 0.09381844 
    ## Run 465 stress 0.090304 
    ## Run 466 stress 0.09407969 
    ## Run 467 stress 0.09445501 
    ## Run 468 stress 0.08440273 
    ## ... Procrustes: rmse 0.000231072  max resid 0.0004994104 
    ## ... Similar to previous best
    ## Run 469 stress 0.09760831 
    ## Run 470 stress 0.09168961 
    ## Run 471 stress 0.09337317 
    ## Run 472 stress 0.09337274 
    ## Run 473 stress 0.09445487 
    ## Run 474 stress 0.09380576 
    ## Run 475 stress 0.09286089 
    ## Run 476 stress 0.08503472 
    ## Run 477 stress 0.08973863 
    ## Run 478 stress 0.0926834 
    ## Run 479 stress 0.09416139 
    ## Run 480 stress 0.103805 
    ## Run 481 stress 0.09610763 
    ## Run 482 stress 0.08503506 
    ## Run 483 stress 0.08973867 
    ## Run 484 stress 0.08503726 
    ## Run 485 stress 0.08973881 
    ## Run 486 stress 0.09407983 
    ## Run 487 stress 0.09030397 
    ## Run 488 stress 0.09417766 
    ## Run 489 stress 0.09908183 
    ## Run 490 stress 0.09159088 
    ## Run 491 stress 0.09374266 
    ## Run 492 stress 0.09268355 
    ## Run 493 stress 0.09337265 
    ## Run 494 stress 0.08773486 
    ## Run 495 stress 0.09465905 
    ## Run 496 stress 0.09168939 
    ## Run 497 stress 0.08973885 
    ## Run 498 stress 0.09308955 
    ## Run 499 stress 0.08973862 
    ## Run 500 stress 0.09030405 
    ## *** Best solution repeated 20 times

``` r
# Ocean sites and mixed lakes
SD_beta_geo_OM_NMDS <- metaMDS(SD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07365783 
    ## Run 1 stress 0.0800047 
    ## Run 2 stress 0.07365785 
    ## ... Procrustes: rmse 7.199295e-05  max resid 0.0001707482 
    ## ... Similar to previous best
    ## Run 3 stress 0.08000474 
    ## Run 4 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745427  max resid 0.05110149 
    ## Run 5 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001116103  max resid 0.0002661578 
    ## ... Similar to previous best
    ## Run 6 stress 0.08288052 
    ## Run 7 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744087  max resid 0.05101074 
    ## Run 8 stress 0.07629233 
    ## Run 9 stress 0.07629246 
    ## Run 10 stress 0.07365783 
    ## ... Procrustes: rmse 1.521608e-05  max resid 3.395649e-05 
    ## ... Similar to previous best
    ## Run 11 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744064  max resid 0.05101072 
    ## Run 12 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 6.223016e-06  max resid 1.3394e-05 
    ## ... Similar to previous best
    ## Run 13 stress 0.07732931 
    ## Run 14 stress 0.07365783 
    ## ... Procrustes: rmse 1.225993e-05  max resid 2.726226e-05 
    ## ... Similar to previous best
    ## Run 15 stress 0.07365784 
    ## ... Procrustes: rmse 4.124043e-05  max resid 9.560816e-05 
    ## ... Similar to previous best
    ## Run 16 stress 0.08000471 
    ## Run 17 stress 0.08000473 
    ## Run 18 stress 0.3495549 
    ## Run 19 stress 0.07378229 
    ## ... Procrustes: rmse 0.01743629  max resid 0.05098557 
    ## Run 20 stress 0.07365783 
    ## ... Procrustes: rmse 1.257562e-05  max resid 2.71727e-05 
    ## ... Similar to previous best
    ## Run 21 stress 0.0800047 
    ## Run 22 stress 0.07732942 
    ## Run 23 stress 0.07629249 
    ## Run 24 stress 0.08000468 
    ## Run 25 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001293507  max resid 0.0003014151 
    ## ... Similar to previous best
    ## Run 26 stress 0.07365784 
    ## ... Procrustes: rmse 3.300202e-05  max resid 7.768904e-05 
    ## ... Similar to previous best
    ## Run 27 stress 0.08000468 
    ## Run 28 stress 0.07629237 
    ## Run 29 stress 0.07365785 
    ## ... Procrustes: rmse 1.680288e-05  max resid 3.383582e-05 
    ## ... Similar to previous best
    ## Run 30 stress 0.07365783 
    ## ... Procrustes: rmse 1.058289e-05  max resid 2.347609e-05 
    ## ... Similar to previous best
    ## Run 31 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 5.502669e-06  max resid 9.967863e-06 
    ## ... Similar to previous best
    ## Run 32 stress 0.07365786 
    ## ... Procrustes: rmse 9.086676e-05  max resid 0.0002147353 
    ## ... Similar to previous best
    ## Run 33 stress 0.07629242 
    ## Run 34 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001111493  max resid 0.0002578643 
    ## ... Similar to previous best
    ## Run 35 stress 0.07365786 
    ## ... Procrustes: rmse 9.645696e-05  max resid 0.0002290405 
    ## ... Similar to previous best
    ## Run 36 stress 0.3495545 
    ## Run 37 stress 0.07365783 
    ## ... Procrustes: rmse 1.833447e-05  max resid 4.387005e-05 
    ## ... Similar to previous best
    ## Run 38 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001154352  max resid 0.0002622377 
    ## ... Similar to previous best
    ## Run 39 stress 0.07365786 
    ## ... Procrustes: rmse 6.765926e-05  max resid 0.0001532653 
    ## ... Similar to previous best
    ## Run 40 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744494  max resid 0.05102138 
    ## Run 41 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745516  max resid 0.05108582 
    ## Run 42 stress 0.07365787 
    ## ... Procrustes: rmse 9.17632e-05  max resid 0.0002129568 
    ## ... Similar to previous best
    ## Run 43 stress 0.07629242 
    ## Run 44 stress 0.07378229 
    ## ... Procrustes: rmse 0.01747262  max resid 0.05113425 
    ## Run 45 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744407  max resid 0.05103127 
    ## Run 46 stress 0.08233757 
    ## Run 47 stress 0.0800047 
    ## Run 48 stress 0.08000472 
    ## Run 49 stress 0.07365784 
    ## ... Procrustes: rmse 5.505665e-05  max resid 0.0001277411 
    ## ... Similar to previous best
    ## Run 50 stress 0.07365784 
    ## ... Procrustes: rmse 6.304888e-05  max resid 0.0001496693 
    ## ... Similar to previous best
    ## Run 51 stress 0.08000468 
    ## Run 52 stress 0.0800047 
    ## Run 53 stress 0.07629248 
    ## Run 54 stress 0.07629242 
    ## Run 55 stress 0.08233755 
    ## Run 56 stress 0.07629237 
    ## Run 57 stress 0.07365783 
    ## ... Procrustes: rmse 2.702122e-06  max resid 6.468296e-06 
    ## ... Similar to previous best
    ## Run 58 stress 0.08233765 
    ## Run 59 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744121  max resid 0.05101143 
    ## Run 60 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744689  max resid 0.05105529 
    ## Run 61 stress 0.07732926 
    ## Run 62 stress 0.07365784 
    ## ... Procrustes: rmse 3.138927e-05  max resid 7.2587e-05 
    ## ... Similar to previous best
    ## Run 63 stress 0.0800047 
    ## Run 64 stress 0.07629246 
    ## Run 65 stress 0.07732928 
    ## Run 66 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001235313  max resid 0.000288803 
    ## ... Similar to previous best
    ## Run 67 stress 0.07629246 
    ## Run 68 stress 0.08288058 
    ## Run 69 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744746  max resid 0.05103857 
    ## Run 70 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744657  max resid 0.05105497 
    ## Run 71 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744  max resid 0.05099424 
    ## Run 72 stress 0.07732928 
    ## Run 73 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001220808  max resid 0.0002901163 
    ## ... Similar to previous best
    ## Run 74 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743651  max resid 0.05100881 
    ## Run 75 stress 0.07629236 
    ## Run 76 stress 0.07365788 
    ## ... Procrustes: rmse 0.000110096  max resid 0.0002590426 
    ## ... Similar to previous best
    ## Run 77 stress 0.07629238 
    ## Run 78 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744596  max resid 0.05104695 
    ## Run 79 stress 0.07629236 
    ## Run 80 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744595  max resid 0.05105204 
    ## Run 81 stress 0.08233777 
    ## Run 82 stress 0.07629237 
    ## Run 83 stress 0.07365783 
    ## ... Procrustes: rmse 5.082192e-06  max resid 1.05509e-05 
    ## ... Similar to previous best
    ## Run 84 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744204  max resid 0.05101537 
    ## Run 85 stress 0.07629233 
    ## Run 86 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744388  max resid 0.05102682 
    ## Run 87 stress 0.07629239 
    ## Run 88 stress 0.07732933 
    ## Run 89 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001220251  max resid 0.0002901647 
    ## ... Similar to previous best
    ## Run 90 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001354605  max resid 0.0003196983 
    ## ... Similar to previous best
    ## Run 91 stress 0.08233764 
    ## Run 92 stress 0.08233769 
    ## Run 93 stress 0.07365796 
    ## ... Procrustes: rmse 0.0001886328  max resid 0.0004474459 
    ## ... Similar to previous best
    ## Run 94 stress 0.07365784 
    ## ... Procrustes: rmse 6.622601e-05  max resid 0.0001566644 
    ## ... Similar to previous best
    ## Run 95 stress 0.07629237 
    ## Run 96 stress 0.07365787 
    ## ... Procrustes: rmse 8.908961e-05  max resid 0.0002040485 
    ## ... Similar to previous best
    ## Run 97 stress 0.08000474 
    ## Run 98 stress 0.0800047 
    ## Run 99 stress 0.08000469 
    ## Run 100 stress 0.08000467 
    ## Run 101 stress 0.07732928 
    ## Run 102 stress 0.07629235 
    ## Run 103 stress 0.07378231 
    ## ... Procrustes: rmse 0.017438  max resid 0.05104503 
    ## Run 104 stress 0.07629234 
    ## Run 105 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745921  max resid 0.05109297 
    ## Run 106 stress 0.08000468 
    ## Run 107 stress 0.07365786 
    ## ... Procrustes: rmse 8.94429e-05  max resid 0.0002114682 
    ## ... Similar to previous best
    ## Run 108 stress 0.07629234 
    ## Run 109 stress 0.08233771 
    ## Run 110 stress 0.07365783 
    ## ... Procrustes: rmse 1.398941e-05  max resid 2.969151e-05 
    ## ... Similar to previous best
    ## Run 111 stress 0.07629234 
    ## Run 112 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744374  max resid 0.05102837 
    ## Run 113 stress 0.07732938 
    ## Run 114 stress 0.07629234 
    ## Run 115 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744856  max resid 0.05106899 
    ## Run 116 stress 0.07365783 
    ## ... Procrustes: rmse 1.356465e-05  max resid 3.196535e-05 
    ## ... Similar to previous best
    ## Run 117 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744647  max resid 0.05105108 
    ## Run 118 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745243  max resid 0.05108623 
    ## Run 119 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174445  max resid 0.05104264 
    ## Run 120 stress 0.07629248 
    ## Run 121 stress 0.07365788 
    ## ... Procrustes: rmse 3.524947e-05  max resid 6.922103e-05 
    ## ... Similar to previous best
    ## Run 122 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001118578  max resid 0.0002649945 
    ## ... Similar to previous best
    ## Run 123 stress 0.08000468 
    ## Run 124 stress 0.07629236 
    ## Run 125 stress 0.3578479 
    ## Run 126 stress 0.07365785 
    ## ... Procrustes: rmse 7.582461e-05  max resid 0.0001803348 
    ## ... Similar to previous best
    ## Run 127 stress 0.07732928 
    ## Run 128 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744431  max resid 0.05103453 
    ## Run 129 stress 0.08000473 
    ## Run 130 stress 0.07732935 
    ## Run 131 stress 0.07365783 
    ## ... Procrustes: rmse 1.916565e-05  max resid 4.539603e-05 
    ## ... Similar to previous best
    ## Run 132 stress 0.07629234 
    ## Run 133 stress 0.07629239 
    ## Run 134 stress 0.07365786 
    ## ... Procrustes: rmse 8.26144e-05  max resid 0.0001942563 
    ## ... Similar to previous best
    ## Run 135 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001394752  max resid 0.0003297984 
    ## ... Similar to previous best
    ## Run 136 stress 0.08000468 
    ## Run 137 stress 0.07365786 
    ## ... Procrustes: rmse 8.203887e-05  max resid 0.0001933887 
    ## ... Similar to previous best
    ## Run 138 stress 0.07629235 
    ## Run 139 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744596  max resid 0.05105496 
    ## Run 140 stress 0.07365784 
    ## ... Procrustes: rmse 4.818057e-05  max resid 0.0001142466 
    ## ... Similar to previous best
    ## Run 141 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745943  max resid 0.05111397 
    ## Run 142 stress 0.07629238 
    ## Run 143 stress 0.07365784 
    ## ... Procrustes: rmse 6.156227e-05  max resid 0.0001449831 
    ## ... Similar to previous best
    ## Run 144 stress 0.08000469 
    ## Run 145 stress 0.07378229 
    ## ... Procrustes: rmse 0.01746315  max resid 0.05113686 
    ## Run 146 stress 0.0737823 
    ## ... Procrustes: rmse 0.01744022  max resid 0.05104864 
    ## Run 147 stress 0.07732927 
    ## Run 148 stress 0.07365785 
    ## ... Procrustes: rmse 6.922608e-05  max resid 0.0001645271 
    ## ... Similar to previous best
    ## Run 149 stress 0.08000478 
    ## Run 150 stress 0.0762924 
    ## Run 151 stress 0.07629235 
    ## Run 152 stress 0.0773294 
    ## Run 153 stress 0.07365786 
    ## ... Procrustes: rmse 9.403627e-05  max resid 0.0002187682 
    ## ... Similar to previous best
    ## Run 154 stress 0.08000467 
    ## Run 155 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001047349  max resid 0.000244588 
    ## ... Similar to previous best
    ## Run 156 stress 0.07365783 
    ## ... Procrustes: rmse 2.633763e-05  max resid 6.318439e-05 
    ## ... Similar to previous best
    ## Run 157 stress 0.08000469 
    ## Run 158 stress 0.07365783 
    ## ... Procrustes: rmse 1.494422e-05  max resid 3.384349e-05 
    ## ... Similar to previous best
    ## Run 159 stress 0.08000469 
    ## Run 160 stress 0.07629251 
    ## Run 161 stress 0.07365783 
    ## ... Procrustes: rmse 1.669225e-05  max resid 3.647318e-05 
    ## ... Similar to previous best
    ## Run 162 stress 0.07629236 
    ## Run 163 stress 0.07365784 
    ## ... Procrustes: rmse 2.883394e-05  max resid 6.535301e-05 
    ## ... Similar to previous best
    ## Run 164 stress 0.07732938 
    ## Run 165 stress 0.07365785 
    ## ... Procrustes: rmse 6.854828e-05  max resid 0.0001613628 
    ## ... Similar to previous best
    ## Run 166 stress 0.07629239 
    ## Run 167 stress 0.07365784 
    ## ... Procrustes: rmse 3.714541e-05  max resid 8.875656e-05 
    ## ... Similar to previous best
    ## Run 168 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744539  max resid 0.05104331 
    ## Run 169 stress 0.0762924 
    ## Run 170 stress 0.0773293 
    ## Run 171 stress 0.07629236 
    ## Run 172 stress 0.07365787 
    ## ... Procrustes: rmse 4.184198e-05  max resid 7.696855e-05 
    ## ... Similar to previous best
    ## Run 173 stress 0.07732927 
    ## Run 174 stress 0.07629235 
    ## Run 175 stress 0.07365786 
    ## ... Procrustes: rmse 7.817031e-05  max resid 0.0001785093 
    ## ... Similar to previous best
    ## Run 176 stress 0.0737823 
    ## ... Procrustes: rmse 0.01743699  max resid 0.05097617 
    ## Run 177 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001013133  max resid 0.000227782 
    ## ... Similar to previous best
    ## Run 178 stress 0.07365784 
    ## ... Procrustes: rmse 3.687535e-05  max resid 8.783271e-05 
    ## ... Similar to previous best
    ## Run 179 stress 0.0773293 
    ## Run 180 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744316  max resid 0.05102152 
    ## Run 181 stress 0.07365792 
    ## ... Procrustes: rmse 0.000149003  max resid 0.0003485664 
    ## ... Similar to previous best
    ## Run 182 stress 0.07365783 
    ## ... Procrustes: rmse 1.307172e-05  max resid 2.988948e-05 
    ## ... Similar to previous best
    ## Run 183 stress 0.07365785 
    ## ... Procrustes: rmse 3.113442e-05  max resid 6.819901e-05 
    ## ... Similar to previous best
    ## Run 184 stress 0.07378232 
    ## ... Procrustes: rmse 0.01742688  max resid 0.0509278 
    ## Run 185 stress 0.07629242 
    ## Run 186 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744361  max resid 0.05105862 
    ## Run 187 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744861  max resid 0.0510595 
    ## Run 188 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744551  max resid 0.0510552 
    ## Run 189 stress 0.07365784 
    ## ... Procrustes: rmse 4.695895e-05  max resid 0.000110061 
    ## ... Similar to previous best
    ## Run 190 stress 0.08000469 
    ## Run 191 stress 0.07629237 
    ## Run 192 stress 0.07732947 
    ## Run 193 stress 0.08000468 
    ## Run 194 stress 0.07365784 
    ## ... Procrustes: rmse 4.035742e-05  max resid 9.348839e-05 
    ## ... Similar to previous best
    ## Run 195 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743439  max resid 0.05099158 
    ## Run 196 stress 0.07629246 
    ## Run 197 stress 0.0762924 
    ## Run 198 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001207427  max resid 0.0002885019 
    ## ... Similar to previous best
    ## Run 199 stress 0.07629243 
    ## Run 200 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174456  max resid 0.05105061 
    ## Run 201 stress 0.08233771 
    ## Run 202 stress 0.07365785 
    ## ... Procrustes: rmse 8.212015e-05  max resid 0.000194231 
    ## ... Similar to previous best
    ## Run 203 stress 0.07365785 
    ## ... Procrustes: rmse 5.182061e-05  max resid 0.0001182027 
    ## ... Similar to previous best
    ## Run 204 stress 0.08000467 
    ## Run 205 stress 0.0823377 
    ## Run 206 stress 0.07365783 
    ## ... Procrustes: rmse 2.735372e-05  max resid 6.24336e-05 
    ## ... Similar to previous best
    ## Run 207 stress 0.08000467 
    ## Run 208 stress 0.07629235 
    ## Run 209 stress 0.07629247 
    ## Run 210 stress 0.07378227 
    ## ... Procrustes: rmse 0.01746075  max resid 0.05110828 
    ## Run 211 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743765  max resid 0.05101566 
    ## Run 212 stress 0.08000474 
    ## Run 213 stress 0.07365784 
    ## ... Procrustes: rmse 1.196268e-05  max resid 2.498739e-05 
    ## ... Similar to previous best
    ## Run 214 stress 0.296091 
    ## Run 215 stress 0.07629236 
    ## Run 216 stress 0.275406 
    ## Run 217 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744622  max resid 0.05103984 
    ## Run 218 stress 0.07629233 
    ## Run 219 stress 0.07629238 
    ## Run 220 stress 0.07365784 
    ## ... Procrustes: rmse 3.290028e-05  max resid 7.674542e-05 
    ## ... Similar to previous best
    ## Run 221 stress 0.08288036 
    ## Run 222 stress 0.07365786 
    ## ... Procrustes: rmse 7.95657e-05  max resid 0.0001868389 
    ## ... Similar to previous best
    ## Run 223 stress 0.08000469 
    ## Run 224 stress 0.07365783 
    ## ... New best solution
    ## ... Procrustes: rmse 2.698418e-06  max resid 5.59165e-06 
    ## ... Similar to previous best
    ## Run 225 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001081405  max resid 0.000249861 
    ## ... Similar to previous best
    ## Run 226 stress 0.07629234 
    ## Run 227 stress 0.07365783 
    ## ... Procrustes: rmse 2.21094e-05  max resid 5.118527e-05 
    ## ... Similar to previous best
    ## Run 228 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001158563  max resid 0.0002658484 
    ## ... Similar to previous best
    ## Run 229 stress 0.08233758 
    ## Run 230 stress 0.07365791 
    ## ... Procrustes: rmse 0.0001500227  max resid 0.0003556114 
    ## ... Similar to previous best
    ## Run 231 stress 0.08000467 
    ## Run 232 stress 0.07365784 
    ## ... Procrustes: rmse 1.843961e-05  max resid 4.080903e-05 
    ## ... Similar to previous best
    ## Run 233 stress 0.08000472 
    ## Run 234 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744509  max resid 0.05104672 
    ## Run 235 stress 0.08000468 
    ## Run 236 stress 0.07732926 
    ## Run 237 stress 0.07732934 
    ## Run 238 stress 0.07365783 
    ## ... Procrustes: rmse 2.104643e-05  max resid 4.98594e-05 
    ## ... Similar to previous best
    ## Run 239 stress 0.07629239 
    ## Run 240 stress 0.08233751 
    ## Run 241 stress 0.08288055 
    ## Run 242 stress 0.08000478 
    ## Run 243 stress 0.07732937 
    ## Run 244 stress 0.07365784 
    ## ... Procrustes: rmse 4.376898e-05  max resid 0.0001029144 
    ## ... Similar to previous best
    ## Run 245 stress 0.0773293 
    ## Run 246 stress 0.07365783 
    ## ... Procrustes: rmse 4.227591e-06  max resid 6.785856e-06 
    ## ... Similar to previous best
    ## Run 247 stress 0.08000468 
    ## Run 248 stress 0.357848 
    ## Run 249 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744612  max resid 0.05105461 
    ## Run 250 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744565  max resid 0.05105485 
    ## Run 251 stress 0.07629237 
    ## Run 252 stress 0.08000469 
    ## Run 253 stress 0.07378228 
    ## ... Procrustes: rmse 0.0174399  max resid 0.05104433 
    ## Run 254 stress 0.0800047 
    ## Run 255 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744833  max resid 0.05107742 
    ## Run 256 stress 0.07629236 
    ## Run 257 stress 0.07365783 
    ## ... Procrustes: rmse 2.614369e-05  max resid 6.153992e-05 
    ## ... Similar to previous best
    ## Run 258 stress 0.07365783 
    ## ... Procrustes: rmse 2.011357e-05  max resid 4.724998e-05 
    ## ... Similar to previous best
    ## Run 259 stress 0.07365784 
    ## ... Procrustes: rmse 4.009211e-05  max resid 9.515321e-05 
    ## ... Similar to previous best
    ## Run 260 stress 0.0762925 
    ## Run 261 stress 0.08233762 
    ## Run 262 stress 0.07629253 
    ## Run 263 stress 0.07365784 
    ## ... Procrustes: rmse 3.70172e-05  max resid 8.844479e-05 
    ## ... Similar to previous best
    ## Run 264 stress 0.0736579 
    ## ... Procrustes: rmse 0.0001414531  max resid 0.0003317608 
    ## ... Similar to previous best
    ## Run 265 stress 0.07378232 
    ## ... Procrustes: rmse 0.01743836  max resid 0.05104904 
    ## Run 266 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744837  max resid 0.05106575 
    ## Run 267 stress 0.07365783 
    ## ... Procrustes: rmse 4.555776e-06  max resid 9.446523e-06 
    ## ... Similar to previous best
    ## Run 268 stress 0.07629233 
    ## Run 269 stress 0.07378233 
    ## ... Procrustes: rmse 0.0174431  max resid 0.05099622 
    ## Run 270 stress 0.07365784 
    ## ... Procrustes: rmse 3.277226e-05  max resid 7.640036e-05 
    ## ... Similar to previous best
    ## Run 271 stress 0.07378229 
    ## ... Procrustes: rmse 0.0174445  max resid 0.05101828 
    ## Run 272 stress 0.07629237 
    ## Run 273 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001199166  max resid 0.0002800544 
    ## ... Similar to previous best
    ## Run 274 stress 0.07365783 
    ## ... Procrustes: rmse 5.08032e-06  max resid 1.05636e-05 
    ## ... Similar to previous best
    ## Run 275 stress 0.08233765 
    ## Run 276 stress 0.08288058 
    ## Run 277 stress 0.07629246 
    ## Run 278 stress 0.07629255 
    ## Run 279 stress 0.08000469 
    ## Run 280 stress 0.08233758 
    ## Run 281 stress 0.07365786 
    ## ... Procrustes: rmse 3.275407e-05  max resid 5.367837e-05 
    ## ... Similar to previous best
    ## Run 282 stress 0.07629235 
    ## Run 283 stress 0.07378233 
    ## ... Procrustes: rmse 0.01743754  max resid 0.05104559 
    ## Run 284 stress 0.07365784 
    ## ... Procrustes: rmse 2.87211e-05  max resid 6.317871e-05 
    ## ... Similar to previous best
    ## Run 285 stress 0.07629235 
    ## Run 286 stress 0.07365784 
    ## ... Procrustes: rmse 3.50904e-05  max resid 8.148196e-05 
    ## ... Similar to previous best
    ## Run 287 stress 0.08000472 
    ## Run 288 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745435  max resid 0.05110607 
    ## Run 289 stress 0.0762924 
    ## Run 290 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001110284  max resid 0.0002629245 
    ## ... Similar to previous best
    ## Run 291 stress 0.08000469 
    ## Run 292 stress 0.07732925 
    ## Run 293 stress 0.2493376 
    ## Run 294 stress 0.08000468 
    ## Run 295 stress 0.08000467 
    ## Run 296 stress 0.07365783 
    ## ... Procrustes: rmse 3.610923e-05  max resid 8.575119e-05 
    ## ... Similar to previous best
    ## Run 297 stress 0.07365787 
    ## ... Procrustes: rmse 5.61296e-05  max resid 0.0001311357 
    ## ... Similar to previous best
    ## Run 298 stress 0.07365783 
    ## ... Procrustes: rmse 3.069892e-05  max resid 7.333129e-05 
    ## ... Similar to previous best
    ## Run 299 stress 0.07629239 
    ## Run 300 stress 0.07629238 
    ## Run 301 stress 0.07629243 
    ## Run 302 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744781  max resid 0.0510643 
    ## Run 303 stress 0.07629234 
    ## Run 304 stress 0.07629237 
    ## Run 305 stress 0.07629241 
    ## Run 306 stress 0.07378226 
    ## ... Procrustes: rmse 0.017442  max resid 0.05102823 
    ## Run 307 stress 0.07629234 
    ## Run 308 stress 0.07629233 
    ## Run 309 stress 0.07378226 
    ## ... Procrustes: rmse 0.01745005  max resid 0.05106262 
    ## Run 310 stress 0.08288039 
    ## Run 311 stress 0.07365785 
    ## ... Procrustes: rmse 2.551908e-05  max resid 5.328274e-05 
    ## ... Similar to previous best
    ## Run 312 stress 0.07732926 
    ## Run 313 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744595  max resid 0.05106412 
    ## Run 314 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744695  max resid 0.05103605 
    ## Run 315 stress 0.07629233 
    ## Run 316 stress 0.08000467 
    ## Run 317 stress 0.08233771 
    ## Run 318 stress 0.07629238 
    ## Run 319 stress 0.07732943 
    ## Run 320 stress 0.07365783 
    ## ... Procrustes: rmse 2.107209e-05  max resid 4.974629e-05 
    ## ... Similar to previous best
    ## Run 321 stress 0.07629241 
    ## Run 322 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174427  max resid 0.05102496 
    ## Run 323 stress 0.07365783 
    ## ... Procrustes: rmse 1.878674e-05  max resid 4.011581e-05 
    ## ... Similar to previous best
    ## Run 324 stress 0.0823376 
    ## Run 325 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744341  max resid 0.05105396 
    ## Run 326 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744451  max resid 0.05105701 
    ## Run 327 stress 0.07365784 
    ## ... Procrustes: rmse 1.880674e-05  max resid 3.782189e-05 
    ## ... Similar to previous best
    ## Run 328 stress 0.07629245 
    ## Run 329 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001043012  max resid 0.0002468583 
    ## ... Similar to previous best
    ## Run 330 stress 0.08000467 
    ## Run 331 stress 0.07732937 
    ## Run 332 stress 0.07365783 
    ## ... Procrustes: rmse 3.17243e-06  max resid 7.074759e-06 
    ## ... Similar to previous best
    ## Run 333 stress 0.08000469 
    ## Run 334 stress 0.08288047 
    ## Run 335 stress 0.07629233 
    ## Run 336 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744487  max resid 0.05104351 
    ## Run 337 stress 0.07365788 
    ## ... Procrustes: rmse 7.016893e-05  max resid 0.0001493612 
    ## ... Similar to previous best
    ## Run 338 stress 0.0773293 
    ## Run 339 stress 0.07378232 
    ## ... Procrustes: rmse 0.01744478  max resid 0.05101098 
    ## Run 340 stress 0.07365783 
    ## ... Procrustes: rmse 6.400235e-06  max resid 1.378678e-05 
    ## ... Similar to previous best
    ## Run 341 stress 0.3411846 
    ## Run 342 stress 0.07629234 
    ## Run 343 stress 0.08000478 
    ## Run 344 stress 0.07365784 
    ## ... Procrustes: rmse 4.628717e-05  max resid 0.000104278 
    ## ... Similar to previous best
    ## Run 345 stress 0.08000468 
    ## Run 346 stress 0.07629234 
    ## Run 347 stress 0.0800047 
    ## Run 348 stress 0.08000468 
    ## Run 349 stress 0.07629247 
    ## Run 350 stress 0.07732934 
    ## Run 351 stress 0.07378228 
    ## ... Procrustes: rmse 0.01743549  max resid 0.0509829 
    ## Run 352 stress 0.08000473 
    ## Run 353 stress 0.07365789 
    ## ... Procrustes: rmse 0.0001260344  max resid 0.0002933354 
    ## ... Similar to previous best
    ## Run 354 stress 0.07629235 
    ## Run 355 stress 0.08000467 
    ## Run 356 stress 0.07365786 
    ## ... Procrustes: rmse 8.35328e-05  max resid 0.0001993903 
    ## ... Similar to previous best
    ## Run 357 stress 0.0762924 
    ## Run 358 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744529  max resid 0.05104887 
    ## Run 359 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745097  max resid 0.05109149 
    ## Run 360 stress 0.07378233 
    ## ... Procrustes: rmse 0.01743489  max resid 0.05104234 
    ## Run 361 stress 0.08233763 
    ## Run 362 stress 0.0762924 
    ## Run 363 stress 0.07378227 
    ## ... Procrustes: rmse 0.01743328  max resid 0.05098831 
    ## Run 364 stress 0.07365783 
    ## ... Procrustes: rmse 6.13709e-06  max resid 1.355161e-05 
    ## ... Similar to previous best
    ## Run 365 stress 0.07629251 
    ## Run 366 stress 0.07378229 
    ## ... Procrustes: rmse 0.0174538  max resid 0.05105977 
    ## Run 367 stress 0.07365785 
    ## ... Procrustes: rmse 6.244568e-05  max resid 0.0001482136 
    ## ... Similar to previous best
    ## Run 368 stress 0.07378227 
    ## ... Procrustes: rmse 0.0174444  max resid 0.05103037 
    ## Run 369 stress 0.07629233 
    ## Run 370 stress 0.0800047 
    ## Run 371 stress 0.07378226 
    ## ... Procrustes: rmse 0.0174427  max resid 0.05103597 
    ## Run 372 stress 0.08000467 
    ## Run 373 stress 0.08000493 
    ## Run 374 stress 0.08233758 
    ## Run 375 stress 0.08000467 
    ## Run 376 stress 0.0800047 
    ## Run 377 stress 0.07378229 
    ## ... Procrustes: rmse 0.0174299  max resid 0.05096405 
    ## Run 378 stress 0.07629243 
    ## Run 379 stress 0.07629236 
    ## Run 380 stress 0.07365796 
    ## ... Procrustes: rmse 0.0001935303  max resid 0.0004535154 
    ## ... Similar to previous best
    ## Run 381 stress 0.0773295 
    ## Run 382 stress 0.07629238 
    ## Run 383 stress 0.07365785 
    ## ... Procrustes: rmse 6.049493e-05  max resid 0.0001440344 
    ## ... Similar to previous best
    ## Run 384 stress 0.07365784 
    ## ... Procrustes: rmse 4.180233e-05  max resid 9.956624e-05 
    ## ... Similar to previous best
    ## Run 385 stress 0.0800047 
    ## Run 386 stress 0.07629243 
    ## Run 387 stress 0.08288034 
    ## Run 388 stress 0.2493376 
    ## Run 389 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744125  max resid 0.05102749 
    ## Run 390 stress 0.08233756 
    ## Run 391 stress 0.07629234 
    ## Run 392 stress 0.07365792 
    ## ... Procrustes: rmse 0.0001492644  max resid 0.0003473266 
    ## ... Similar to previous best
    ## Run 393 stress 0.08000472 
    ## Run 394 stress 0.07378231 
    ## ... Procrustes: rmse 0.01743814  max resid 0.05104845 
    ## Run 395 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744158  max resid 0.05101296 
    ## Run 396 stress 0.07629242 
    ## Run 397 stress 0.07365784 
    ## ... Procrustes: rmse 3.687572e-05  max resid 8.805413e-05 
    ## ... Similar to previous best
    ## Run 398 stress 0.07378229 
    ## ... Procrustes: rmse 0.01743221  max resid 0.05101542 
    ## Run 399 stress 0.08288042 
    ## Run 400 stress 0.3267051 
    ## Run 401 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744853  max resid 0.05107502 
    ## Run 402 stress 0.07629242 
    ## Run 403 stress 0.07629248 
    ## Run 404 stress 0.07365784 
    ## ... Procrustes: rmse 5.85058e-05  max resid 0.0001378811 
    ## ... Similar to previous best
    ## Run 405 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744134  max resid 0.05102845 
    ## Run 406 stress 0.07732931 
    ## Run 407 stress 0.08000471 
    ## Run 408 stress 0.07378231 
    ## ... Procrustes: rmse 0.017458  max resid 0.05112268 
    ## Run 409 stress 0.07629234 
    ## Run 410 stress 0.08000471 
    ## Run 411 stress 0.07629234 
    ## Run 412 stress 0.07629234 
    ## Run 413 stress 0.0823375 
    ## Run 414 stress 0.08288063 
    ## Run 415 stress 0.08000468 
    ## Run 416 stress 0.07365784 
    ## ... Procrustes: rmse 1.510218e-05  max resid 2.99934e-05 
    ## ... Similar to previous best
    ## Run 417 stress 0.07732938 
    ## Run 418 stress 0.07378228 
    ## ... Procrustes: rmse 0.01744137  max resid 0.05100923 
    ## Run 419 stress 0.07378231 
    ## ... Procrustes: rmse 0.01746961  max resid 0.0511724 
    ## Run 420 stress 0.08000468 
    ## Run 421 stress 0.07365783 
    ## ... Procrustes: rmse 2.541342e-05  max resid 6.011785e-05 
    ## ... Similar to previous best
    ## Run 422 stress 0.08233763 
    ## Run 423 stress 0.07365786 
    ## ... Procrustes: rmse 6.869912e-05  max resid 0.0001530737 
    ## ... Similar to previous best
    ## Run 424 stress 0.07629234 
    ## Run 425 stress 0.07365783 
    ## ... Procrustes: rmse 1.83523e-05  max resid 3.817812e-05 
    ## ... Similar to previous best
    ## Run 426 stress 0.07365786 
    ## ... Procrustes: rmse 8.908378e-05  max resid 0.0002095852 
    ## ... Similar to previous best
    ## Run 427 stress 0.07378227 
    ## ... Procrustes: rmse 0.01744222  max resid 0.05103707 
    ## Run 428 stress 0.0800047 
    ## Run 429 stress 0.07365784 
    ## ... Procrustes: rmse 4.964959e-05  max resid 0.0001165172 
    ## ... Similar to previous best
    ## Run 430 stress 0.07365784 
    ## ... Procrustes: rmse 4.346775e-05  max resid 0.00010137 
    ## ... Similar to previous best
    ## Run 431 stress 0.07629244 
    ## Run 432 stress 0.07629236 
    ## Run 433 stress 0.07629238 
    ## Run 434 stress 0.07365787 
    ## ... Procrustes: rmse 8.632645e-05  max resid 0.000205757 
    ## ... Similar to previous best
    ## Run 435 stress 0.07365783 
    ## ... Procrustes: rmse 1.63417e-05  max resid 4.181584e-05 
    ## ... Similar to previous best
    ## Run 436 stress 0.08288062 
    ## Run 437 stress 0.07629246 
    ## Run 438 stress 0.07732928 
    ## Run 439 stress 0.07365783 
    ## ... Procrustes: rmse 2.131512e-05  max resid 4.986074e-05 
    ## ... Similar to previous best
    ## Run 440 stress 0.07732932 
    ## Run 441 stress 0.07378226 
    ## ... Procrustes: rmse 0.01743398  max resid 0.05100193 
    ## Run 442 stress 0.07365785 
    ## ... Procrustes: rmse 6.876707e-05  max resid 0.0001596554 
    ## ... Similar to previous best
    ## Run 443 stress 0.08000468 
    ## Run 444 stress 0.08288045 
    ## Run 445 stress 0.07629251 
    ## Run 446 stress 0.07365787 
    ## ... Procrustes: rmse 9.19998e-05  max resid 0.0002096331 
    ## ... Similar to previous best
    ## Run 447 stress 0.07732925 
    ## Run 448 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744383  max resid 0.05104014 
    ## Run 449 stress 0.07365787 
    ## ... Procrustes: rmse 8.879155e-05  max resid 0.000212464 
    ## ... Similar to previous best
    ## Run 450 stress 0.07378228 
    ## ... Procrustes: rmse 0.01740142  max resid 0.05086971 
    ## Run 451 stress 0.08288042 
    ## Run 452 stress 0.07378231 
    ## ... Procrustes: rmse 0.01748999  max resid 0.05125232 
    ## Run 453 stress 0.07365784 
    ## ... Procrustes: rmse 4.781206e-05  max resid 0.0001118462 
    ## ... Similar to previous best
    ## Run 454 stress 0.07365783 
    ## ... Procrustes: rmse 2.348526e-05  max resid 5.537183e-05 
    ## ... Similar to previous best
    ## Run 455 stress 0.08000468 
    ## Run 456 stress 0.07365787 
    ## ... Procrustes: rmse 0.0001041249  max resid 0.000243137 
    ## ... Similar to previous best
    ## Run 457 stress 0.08000472 
    ## Run 458 stress 0.08000469 
    ## Run 459 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744471  max resid 0.05104128 
    ## Run 460 stress 0.07629234 
    ## Run 461 stress 0.07365785 
    ## ... Procrustes: rmse 6.835131e-05  max resid 0.000161384 
    ## ... Similar to previous best
    ## Run 462 stress 0.07629241 
    ## Run 463 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001089979  max resid 0.0002594104 
    ## ... Similar to previous best
    ## Run 464 stress 0.07365789 
    ## ... Procrustes: rmse 0.000130227  max resid 0.000309176 
    ## ... Similar to previous best
    ## Run 465 stress 0.07629234 
    ## Run 466 stress 0.08000469 
    ## Run 467 stress 0.0762925 
    ## Run 468 stress 0.07365786 
    ## ... Procrustes: rmse 8.497409e-05  max resid 0.0001993979 
    ## ... Similar to previous best
    ## Run 469 stress 0.07365784 
    ## ... Procrustes: rmse 2.631843e-05  max resid 5.88162e-05 
    ## ... Similar to previous best
    ## Run 470 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001155993  max resid 0.0002701903 
    ## ... Similar to previous best
    ## Run 471 stress 0.07365784 
    ## ... Procrustes: rmse 3.470055e-05  max resid 8.211725e-05 
    ## ... Similar to previous best
    ## Run 472 stress 0.07365783 
    ## ... Procrustes: rmse 1.098357e-05  max resid 2.509586e-05 
    ## ... Similar to previous best
    ## Run 473 stress 0.07365786 
    ## ... Procrustes: rmse 8.995356e-05  max resid 0.0002127588 
    ## ... Similar to previous best
    ## Run 474 stress 0.07365786 
    ## ... Procrustes: rmse 8.992381e-05  max resid 0.000213912 
    ## ... Similar to previous best
    ## Run 475 stress 0.07378229 
    ## ... Procrustes: rmse 0.01745325  max resid 0.05105734 
    ## Run 476 stress 0.07629235 
    ## Run 477 stress 0.07365785 
    ## ... Procrustes: rmse 4.474676e-05  max resid 9.596884e-05 
    ## ... Similar to previous best
    ## Run 478 stress 0.07732929 
    ## Run 479 stress 0.0800047 
    ## Run 480 stress 0.07732932 
    ## Run 481 stress 0.07365783 
    ## ... Procrustes: rmse 1.342512e-05  max resid 2.790299e-05 
    ## ... Similar to previous best
    ## Run 482 stress 0.0762924 
    ## Run 483 stress 0.07365788 
    ## ... Procrustes: rmse 0.0001197385  max resid 0.0002820311 
    ## ... Similar to previous best
    ## Run 484 stress 0.07378227 
    ## ... Procrustes: rmse 0.01745191  max resid 0.05107406 
    ## Run 485 stress 0.08233758 
    ## Run 486 stress 0.07365783 
    ## ... Procrustes: rmse 2.831031e-05  max resid 6.698304e-05 
    ## ... Similar to previous best
    ## Run 487 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744687  max resid 0.05104988 
    ## Run 488 stress 0.07365785 
    ## ... Procrustes: rmse 7.840757e-05  max resid 0.0001843567 
    ## ... Similar to previous best
    ## Run 489 stress 0.07732934 
    ## Run 490 stress 0.07378228 
    ## ... Procrustes: rmse 0.01743971  max resid 0.0510029 
    ## Run 491 stress 0.07629235 
    ## Run 492 stress 0.07378226 
    ## ... Procrustes: rmse 0.01744555  max resid 0.0510513 
    ## Run 493 stress 0.07732938 
    ## Run 494 stress 0.08000471 
    ## Run 495 stress 0.07378228 
    ## ... Procrustes: rmse 0.01745445  max resid 0.05107318 
    ## Run 496 stress 0.07365784 
    ## ... Procrustes: rmse 4.598382e-05  max resid 0.0001074404 
    ## ... Similar to previous best
    ## Run 497 stress 0.08000472 
    ## Run 498 stress 0.07378229 
    ## ... Procrustes: rmse 0.01744381  max resid 0.05103925 
    ## Run 499 stress 0.07365788 
    ## ... Procrustes: rmse 0.000106899  max resid 0.0002496311 
    ## ... Similar to previous best
    ## Run 500 stress 0.08288044 
    ## *** Best solution repeated 77 times

``` r
# Stratified lakes and ocean sites
SD_beta_geo_SO_NMDS <- metaMDS(SD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.07970526 
    ## Run 1 stress 0.08340291 
    ## Run 2 stress 0.08448459 
    ## Run 3 stress 0.07250813 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0899346  max resid 0.2576286 
    ## Run 4 stress 0.07844946 
    ## Run 5 stress 0.06942777 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04030482  max resid 0.1229911 
    ## Run 6 stress 0.07428318 
    ## Run 7 stress 0.07250812 
    ## Run 8 stress 0.06978198 
    ## ... Procrustes: rmse 0.0132781  max resid 0.03337955 
    ## Run 9 stress 0.08340292 
    ## Run 10 stress 0.08448443 
    ## Run 11 stress 0.07970527 
    ## Run 12 stress 0.08448443 
    ## Run 13 stress 0.08340286 
    ## Run 14 stress 0.07970525 
    ## Run 15 stress 0.06978197 
    ## ... Procrustes: rmse 0.01317506  max resid 0.03316081 
    ## Run 16 stress 0.07428317 
    ## Run 17 stress 0.07428315 
    ## Run 18 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 3.077415e-05  max resid 8.520557e-05 
    ## ... Similar to previous best
    ## Run 19 stress 0.08340287 
    ## Run 20 stress 0.08340288 
    ## Run 21 stress 0.08340287 
    ## Run 22 stress 0.06942778 
    ## ... Procrustes: rmse 2.015661e-05  max resid 5.759644e-05 
    ## ... Similar to previous best
    ## Run 23 stress 0.07250813 
    ## Run 24 stress 0.07970526 
    ## Run 25 stress 0.07970525 
    ## Run 26 stress 0.07970525 
    ## Run 27 stress 0.06942777 
    ## ... Procrustes: rmse 3.753475e-05  max resid 9.480433e-05 
    ## ... Similar to previous best
    ## Run 28 stress 0.08340297 
    ## Run 29 stress 0.07428313 
    ## Run 30 stress 0.07250812 
    ## Run 31 stress 0.06978192 
    ## ... Procrustes: rmse 0.01324389  max resid 0.03331106 
    ## Run 32 stress 0.07250812 
    ## Run 33 stress 0.07970525 
    ## Run 34 stress 0.07428313 
    ## Run 35 stress 0.06978201 
    ## ... Procrustes: rmse 0.01328407  max resid 0.03339303 
    ## Run 36 stress 0.06942777 
    ## ... Procrustes: rmse 1.549045e-05  max resid 4.359374e-05 
    ## ... Similar to previous best
    ## Run 37 stress 0.07428318 
    ## Run 38 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.685676e-05  max resid 6.971075e-05 
    ## ... Similar to previous best
    ## Run 39 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319719  max resid 0.03320716 
    ## Run 40 stress 0.07844966 
    ## Run 41 stress 0.08340287 
    ## Run 42 stress 0.06978197 
    ## ... Procrustes: rmse 0.0132765  max resid 0.0333723 
    ## Run 43 stress 0.08340287 
    ## Run 44 stress 0.08340286 
    ## Run 45 stress 0.07250813 
    ## Run 46 stress 0.07970526 
    ## Run 47 stress 0.06942777 
    ## ... Procrustes: rmse 6.076357e-05  max resid 0.0001562931 
    ## ... Similar to previous best
    ## Run 48 stress 0.07970525 
    ## Run 49 stress 0.07970526 
    ## Run 50 stress 0.07250813 
    ## Run 51 stress 0.07970526 
    ## Run 52 stress 0.07250815 
    ## Run 53 stress 0.06942777 
    ## ... Procrustes: rmse 7.106093e-05  max resid 0.0001829606 
    ## ... Similar to previous best
    ## Run 54 stress 0.08340288 
    ## Run 55 stress 0.06978193 
    ## ... Procrustes: rmse 0.01319138  max resid 0.03319481 
    ## Run 56 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323337  max resid 0.03328378 
    ## Run 57 stress 0.07970525 
    ## Run 58 stress 0.07970526 
    ## Run 59 stress 0.07428313 
    ## Run 60 stress 0.07970526 
    ## Run 61 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326509  max resid 0.0333466 
    ## Run 62 stress 0.07250812 
    ## Run 63 stress 0.07428315 
    ## Run 64 stress 0.07428318 
    ## Run 65 stress 0.07250812 
    ## Run 66 stress 0.08448439 
    ## Run 67 stress 0.07970525 
    ## Run 68 stress 0.07428314 
    ## Run 69 stress 0.07428319 
    ## Run 70 stress 0.06942776 
    ## ... Procrustes: rmse 4.125122e-05  max resid 0.0001066284 
    ## ... Similar to previous best
    ## Run 71 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326213  max resid 0.03334484 
    ## Run 72 stress 0.07844923 
    ## Run 73 stress 0.08340288 
    ## Run 74 stress 0.07970525 
    ## Run 75 stress 0.07428313 
    ## Run 76 stress 0.08340287 
    ## Run 77 stress 0.0742832 
    ## Run 78 stress 0.07250814 
    ## Run 79 stress 0.07250814 
    ## Run 80 stress 0.07250812 
    ## Run 81 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325762  max resid 0.03333403 
    ## Run 82 stress 0.07970526 
    ## Run 83 stress 0.06942777 
    ## ... Procrustes: rmse 4.952952e-05  max resid 0.0001283956 
    ## ... Similar to previous best
    ## Run 84 stress 0.07970525 
    ## Run 85 stress 0.07428316 
    ## Run 86 stress 0.08340288 
    ## Run 87 stress 0.08340291 
    ## Run 88 stress 0.07250813 
    ## Run 89 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326253  max resid 0.03334607 
    ## Run 90 stress 0.08340289 
    ## Run 91 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321907  max resid 0.03325211 
    ## Run 92 stress 0.08340287 
    ## Run 93 stress 0.07250812 
    ## Run 94 stress 0.07970525 
    ## Run 95 stress 0.08340291 
    ## Run 96 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 9.695845e-06  max resid 2.517033e-05 
    ## ... Similar to previous best
    ## Run 97 stress 0.06942776 
    ## ... Procrustes: rmse 1.312995e-05  max resid 3.366228e-05 
    ## ... Similar to previous best
    ## Run 98 stress 0.08340287 
    ## Run 99 stress 0.07250812 
    ## Run 100 stress 0.07428314 
    ## Run 101 stress 0.07970525 
    ## Run 102 stress 0.07428318 
    ## Run 103 stress 0.06942776 
    ## ... Procrustes: rmse 1.788564e-05  max resid 4.625975e-05 
    ## ... Similar to previous best
    ## Run 104 stress 0.07250812 
    ## Run 105 stress 0.08340291 
    ## Run 106 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 9.100316e-06  max resid 2.340305e-05 
    ## ... Similar to previous best
    ## Run 107 stress 0.06942777 
    ## ... Procrustes: rmse 3.554657e-05  max resid 9.264087e-05 
    ## ... Similar to previous best
    ## Run 108 stress 0.07970525 
    ## Run 109 stress 0.07428319 
    ## Run 110 stress 0.07970526 
    ## Run 111 stress 0.06978191 
    ## ... Procrustes: rmse 0.01321798  max resid 0.03324249 
    ## Run 112 stress 0.07250815 
    ## Run 113 stress 0.07844914 
    ## Run 114 stress 0.06978203 
    ## ... Procrustes: rmse 0.01315496  max resid 0.03311639 
    ## Run 115 stress 0.08340288 
    ## Run 116 stress 0.07250814 
    ## Run 117 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321111  max resid 0.03323198 
    ## Run 118 stress 0.07970525 
    ## Run 119 stress 0.07970526 
    ## Run 120 stress 0.07428313 
    ## Run 121 stress 0.07250813 
    ## Run 122 stress 0.07250814 
    ## Run 123 stress 0.06942776 
    ## ... Procrustes: rmse 9.107563e-06  max resid 2.35363e-05 
    ## ... Similar to previous best
    ## Run 124 stress 0.07970525 
    ## Run 125 stress 0.07428313 
    ## Run 126 stress 0.08340286 
    ## Run 127 stress 0.07250812 
    ## Run 128 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323963  max resid 0.03329818 
    ## Run 129 stress 0.07970526 
    ## Run 130 stress 0.08340287 
    ## Run 131 stress 0.07428313 
    ## Run 132 stress 0.07970525 
    ## Run 133 stress 0.07250813 
    ## Run 134 stress 0.07250813 
    ## Run 135 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325676  max resid 0.03333536 
    ## Run 136 stress 0.0697819 
    ## ... Procrustes: rmse 0.01322811  max resid 0.03327414 
    ## Run 137 stress 0.07970525 
    ## Run 138 stress 0.06978201 
    ## ... Procrustes: rmse 0.0131598  max resid 0.03312599 
    ## Run 139 stress 0.07844912 
    ## Run 140 stress 0.06942777 
    ## ... Procrustes: rmse 2.756412e-05  max resid 7.448398e-05 
    ## ... Similar to previous best
    ## Run 141 stress 0.08340287 
    ## Run 142 stress 0.07970525 
    ## Run 143 stress 0.07250812 
    ## Run 144 stress 0.08340289 
    ## Run 145 stress 0.07428315 
    ## Run 146 stress 0.06978195 
    ## ... Procrustes: rmse 0.0132631  max resid 0.03334736 
    ## Run 147 stress 0.08340287 
    ## Run 148 stress 0.06978194 
    ## ... Procrustes: rmse 0.01318181  max resid 0.03317267 
    ## Run 149 stress 0.07970525 
    ## Run 150 stress 0.07250814 
    ## Run 151 stress 0.07428315 
    ## Run 152 stress 0.08448446 
    ## Run 153 stress 0.08340286 
    ## Run 154 stress 0.07428316 
    ## Run 155 stress 0.07428313 
    ## Run 156 stress 0.06978199 
    ## ... Procrustes: rmse 0.01327712  max resid 0.03337203 
    ## Run 157 stress 0.08448435 
    ## Run 158 stress 0.08340299 
    ## Run 159 stress 0.0844845 
    ## Run 160 stress 0.07970525 
    ## Run 161 stress 0.07250813 
    ## Run 162 stress 0.07250812 
    ## Run 163 stress 0.07250812 
    ## Run 164 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326989  max resid 0.03335902 
    ## Run 165 stress 0.06942776 
    ## ... Procrustes: rmse 2.374047e-05  max resid 6.134507e-05 
    ## ... Similar to previous best
    ## Run 166 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319766  max resid 0.03320836 
    ## Run 167 stress 0.08448435 
    ## Run 168 stress 0.06978193 
    ## ... Procrustes: rmse 0.01318711  max resid 0.03318683 
    ## Run 169 stress 0.0834029 
    ## Run 170 stress 0.07970525 
    ## Run 171 stress 0.06942778 
    ## ... Procrustes: rmse 3.148377e-05  max resid 8.899668e-05 
    ## ... Similar to previous best
    ## Run 172 stress 0.07250812 
    ## Run 173 stress 0.0834029 
    ## Run 174 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322675  max resid 0.03327049 
    ## Run 175 stress 0.07428313 
    ## Run 176 stress 0.08340287 
    ## Run 177 stress 0.07428313 
    ## Run 178 stress 0.07970526 
    ## Run 179 stress 0.06978195 
    ## ... Procrustes: rmse 0.01325666  max resid 0.03333632 
    ## Run 180 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 2.377858e-06  max resid 6.22234e-06 
    ## ... Similar to previous best
    ## Run 181 stress 0.07970525 
    ## Run 182 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321988  max resid 0.03325638 
    ## Run 183 stress 0.08448436 
    ## Run 184 stress 0.07250813 
    ## Run 185 stress 0.06942776 
    ## ... Procrustes: rmse 7.474453e-06  max resid 1.917833e-05 
    ## ... Similar to previous best
    ## Run 186 stress 0.06942776 
    ## ... Procrustes: rmse 7.686809e-06  max resid 2.456622e-05 
    ## ... Similar to previous best
    ## Run 187 stress 0.07250815 
    ## Run 188 stress 0.08340287 
    ## Run 189 stress 0.06978197 
    ## ... Procrustes: rmse 0.01327379  max resid 0.03336797 
    ## Run 190 stress 0.06942776 
    ## ... New best solution
    ## ... Procrustes: rmse 1.416165e-06  max resid 3.401441e-06 
    ## ... Similar to previous best
    ## Run 191 stress 0.07250813 
    ## Run 192 stress 0.07970525 
    ## Run 193 stress 0.08340286 
    ## Run 194 stress 0.07844933 
    ## Run 195 stress 0.07250812 
    ## Run 196 stress 0.08340298 
    ## Run 197 stress 0.07844945 
    ## Run 198 stress 0.07250812 
    ## Run 199 stress 0.07250812 
    ## Run 200 stress 0.07428313 
    ## Run 201 stress 0.06978191 
    ## ... Procrustes: rmse 0.01320027  max resid 0.03321376 
    ## Run 202 stress 0.08340294 
    ## Run 203 stress 0.06942776 
    ## ... Procrustes: rmse 1.29587e-05  max resid 3.560754e-05 
    ## ... Similar to previous best
    ## Run 204 stress 0.07970526 
    ## Run 205 stress 0.07428315 
    ## Run 206 stress 0.07428314 
    ## Run 207 stress 0.06978189 
    ## ... Procrustes: rmse 0.01322022  max resid 0.03325676 
    ## Run 208 stress 0.07970525 
    ## Run 209 stress 0.07428316 
    ## Run 210 stress 0.06978201 
    ## ... Procrustes: rmse 0.01316015  max resid 0.03312751 
    ## Run 211 stress 0.07428312 
    ## Run 212 stress 0.07250814 
    ## Run 213 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323351  max resid 0.03328406 
    ## Run 214 stress 0.07970525 
    ## Run 215 stress 0.07428318 
    ## Run 216 stress 0.08448435 
    ## Run 217 stress 0.07844937 
    ## Run 218 stress 0.07428317 
    ## Run 219 stress 0.07250812 
    ## Run 220 stress 0.07970525 
    ## Run 221 stress 0.07970525 
    ## Run 222 stress 0.06942776 
    ## ... Procrustes: rmse 2.063599e-05  max resid 5.403702e-05 
    ## ... Similar to previous best
    ## Run 223 stress 0.06978194 
    ## ... Procrustes: rmse 0.01326102  max resid 0.03334198 
    ## Run 224 stress 0.07250812 
    ## Run 225 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323099  max resid 0.0332796 
    ## Run 226 stress 0.06942776 
    ## ... Procrustes: rmse 9.93006e-06  max resid 2.478912e-05 
    ## ... Similar to previous best
    ## Run 227 stress 0.07428315 
    ## Run 228 stress 0.07250815 
    ## Run 229 stress 0.07970525 
    ## Run 230 stress 0.06942776 
    ## ... Procrustes: rmse 1.612697e-06  max resid 3.541076e-06 
    ## ... Similar to previous best
    ## Run 231 stress 0.07250812 
    ## Run 232 stress 0.08340291 
    ## Run 233 stress 0.07250812 
    ## Run 234 stress 0.07428313 
    ## Run 235 stress 0.07428313 
    ## Run 236 stress 0.07428313 
    ## Run 237 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323901  max resid 0.03329635 
    ## Run 238 stress 0.07970525 
    ## Run 239 stress 0.07250812 
    ## Run 240 stress 0.07970526 
    ## Run 241 stress 0.07428314 
    ## Run 242 stress 0.08448438 
    ## Run 243 stress 0.06942776 
    ## ... Procrustes: rmse 6.389067e-06  max resid 1.537294e-05 
    ## ... Similar to previous best
    ## Run 244 stress 0.07428313 
    ## Run 245 stress 0.06978196 
    ## ... Procrustes: rmse 0.01317842  max resid 0.03316931 
    ## Run 246 stress 0.07844929 
    ## Run 247 stress 0.06978197 
    ## ... Procrustes: rmse 0.01317263  max resid 0.03315403 
    ## Run 248 stress 0.07844919 
    ## Run 249 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319767  max resid 0.03321046 
    ## Run 250 stress 0.06978191 
    ## ... Procrustes: rmse 0.01320612  max resid 0.03322501 
    ## Run 251 stress 0.06978192 
    ## ... Procrustes: rmse 0.01325168  max resid 0.03332395 
    ## Run 252 stress 0.07250813 
    ## Run 253 stress 0.07428317 
    ## Run 254 stress 0.06942778 
    ## ... Procrustes: rmse 4.363026e-05  max resid 0.0001146529 
    ## ... Similar to previous best
    ## Run 255 stress 0.06942777 
    ## ... Procrustes: rmse 5.274296e-05  max resid 0.000135528 
    ## ... Similar to previous best
    ## Run 256 stress 0.07250813 
    ## Run 257 stress 0.07428312 
    ## Run 258 stress 0.07970525 
    ## Run 259 stress 0.07250816 
    ## Run 260 stress 0.07250813 
    ## Run 261 stress 0.07428313 
    ## Run 262 stress 0.0844844 
    ## Run 263 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325597  max resid 0.0333317 
    ## Run 264 stress 0.08340296 
    ## Run 265 stress 0.08340286 
    ## Run 266 stress 0.07970525 
    ## Run 267 stress 0.08340294 
    ## Run 268 stress 0.07250815 
    ## Run 269 stress 0.06978198 
    ## ... Procrustes: rmse 0.01327111  max resid 0.03335905 
    ## Run 270 stress 0.06942777 
    ## ... Procrustes: rmse 4.660754e-05  max resid 0.0001199709 
    ## ... Similar to previous best
    ## Run 271 stress 0.07970525 
    ## Run 272 stress 0.06942776 
    ## ... Procrustes: rmse 1.033178e-05  max resid 2.982442e-05 
    ## ... Similar to previous best
    ## Run 273 stress 0.07970526 
    ## Run 274 stress 0.08340291 
    ## Run 275 stress 0.07250812 
    ## Run 276 stress 0.07428314 
    ## Run 277 stress 0.08340287 
    ## Run 278 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323011  max resid 0.03327701 
    ## Run 279 stress 0.06942776 
    ## ... Procrustes: rmse 4.541585e-06  max resid 9.152764e-06 
    ## ... Similar to previous best
    ## Run 280 stress 0.07250816 
    ## Run 281 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323656  max resid 0.03329083 
    ## Run 282 stress 0.08340286 
    ## Run 283 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323768  max resid 0.03329381 
    ## Run 284 stress 0.07428313 
    ## Run 285 stress 0.06942776 
    ## ... Procrustes: rmse 2.292922e-05  max resid 5.889559e-05 
    ## ... Similar to previous best
    ## Run 286 stress 0.08340296 
    ## Run 287 stress 0.07970525 
    ## Run 288 stress 0.06942777 
    ## ... Procrustes: rmse 4.702485e-05  max resid 0.0001212797 
    ## ... Similar to previous best
    ## Run 289 stress 0.07970526 
    ## Run 290 stress 0.07428316 
    ## Run 291 stress 0.07428315 
    ## Run 292 stress 0.08340287 
    ## Run 293 stress 0.06942776 
    ## ... Procrustes: rmse 2.664436e-05  max resid 6.8637e-05 
    ## ... Similar to previous best
    ## Run 294 stress 0.07250812 
    ## Run 295 stress 0.07250812 
    ## Run 296 stress 0.07428313 
    ## Run 297 stress 0.08340288 
    ## Run 298 stress 0.07428313 
    ## Run 299 stress 0.08340286 
    ## Run 300 stress 0.06942777 
    ## ... Procrustes: rmse 1.730255e-05  max resid 4.778339e-05 
    ## ... Similar to previous best
    ## Run 301 stress 0.06942776 
    ## ... Procrustes: rmse 1.802521e-06  max resid 4.215915e-06 
    ## ... Similar to previous best
    ## Run 302 stress 0.07970525 
    ## Run 303 stress 0.08448449 
    ## Run 304 stress 0.07970525 
    ## Run 305 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324209  max resid 0.03330218 
    ## Run 306 stress 0.08448466 
    ## Run 307 stress 0.06942776 
    ## ... Procrustes: rmse 3.025407e-05  max resid 7.783287e-05 
    ## ... Similar to previous best
    ## Run 308 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325456  max resid 0.03332761 
    ## Run 309 stress 0.06978193 
    ## ... Procrustes: rmse 0.01325802  max resid 0.03333614 
    ## Run 310 stress 0.07970526 
    ## Run 311 stress 0.08340288 
    ## Run 312 stress 0.07970525 
    ## Run 313 stress 0.07250812 
    ## Run 314 stress 0.06942776 
    ## ... Procrustes: rmse 2.278861e-05  max resid 5.907708e-05 
    ## ... Similar to previous best
    ## Run 315 stress 0.06942776 
    ## ... Procrustes: rmse 4.033412e-06  max resid 8.689443e-06 
    ## ... Similar to previous best
    ## Run 316 stress 0.08448439 
    ## Run 317 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325783  max resid 0.03333504 
    ## Run 318 stress 0.07250813 
    ## Run 319 stress 0.06978189 
    ## ... Procrustes: rmse 0.01321677  max resid 0.03324941 
    ## Run 320 stress 0.08340287 
    ## Run 321 stress 0.07970526 
    ## Run 322 stress 0.0697819 
    ## ... Procrustes: rmse 0.01321631  max resid 0.03324787 
    ## Run 323 stress 0.07428323 
    ## Run 324 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323236  max resid 0.03328809 
    ## Run 325 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323106  max resid 0.03327965 
    ## Run 326 stress 0.06942776 
    ## ... Procrustes: rmse 1.222313e-05  max resid 3.16729e-05 
    ## ... Similar to previous best
    ## Run 327 stress 0.07250812 
    ## Run 328 stress 0.06978193 
    ## ... Procrustes: rmse 0.01323191  max resid 0.03327773 
    ## Run 329 stress 0.06942776 
    ## ... Procrustes: rmse 6.158007e-06  max resid 1.572e-05 
    ## ... Similar to previous best
    ## Run 330 stress 0.07428313 
    ## Run 331 stress 0.07970526 
    ## Run 332 stress 0.08340287 
    ## Run 333 stress 0.07250812 
    ## Run 334 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326882  max resid 0.03335453 
    ## Run 335 stress 0.08340292 
    ## Run 336 stress 0.07250813 
    ## Run 337 stress 0.08448455 
    ## Run 338 stress 0.07428315 
    ## Run 339 stress 0.07428319 
    ## Run 340 stress 0.0697819 
    ## ... Procrustes: rmse 0.0132327  max resid 0.03328432 
    ## Run 341 stress 0.08340287 
    ## Run 342 stress 0.08340287 
    ## Run 343 stress 0.06942776 
    ## ... Procrustes: rmse 1.450119e-05  max resid 3.73013e-05 
    ## ... Similar to previous best
    ## Run 344 stress 0.08340287 
    ## Run 345 stress 0.07250812 
    ## Run 346 stress 0.07970525 
    ## Run 347 stress 0.07428313 
    ## Run 348 stress 0.07250814 
    ## Run 349 stress 0.06942777 
    ## ... Procrustes: rmse 3.782562e-05  max resid 9.728607e-05 
    ## ... Similar to previous best
    ## Run 350 stress 0.07970526 
    ## Run 351 stress 0.07428313 
    ## Run 352 stress 0.06942776 
    ## ... Procrustes: rmse 1.000269e-05  max resid 2.620021e-05 
    ## ... Similar to previous best
    ## Run 353 stress 0.08340289 
    ## Run 354 stress 0.07250813 
    ## Run 355 stress 0.07428314 
    ## Run 356 stress 0.07428313 
    ## Run 357 stress 0.07428315 
    ## Run 358 stress 0.07428314 
    ## Run 359 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326406  max resid 0.03334881 
    ## Run 360 stress 0.08340293 
    ## Run 361 stress 0.08340287 
    ## Run 362 stress 0.08340293 
    ## Run 363 stress 0.06942778 
    ## ... Procrustes: rmse 5.845492e-05  max resid 0.000151138 
    ## ... Similar to previous best
    ## Run 364 stress 0.07428314 
    ## Run 365 stress 0.06978199 
    ## ... Procrustes: rmse 0.01328101  max resid 0.033383 
    ## Run 366 stress 0.08340295 
    ## Run 367 stress 0.07428314 
    ## Run 368 stress 0.08340288 
    ## Run 369 stress 0.07428313 
    ## Run 370 stress 0.07428313 
    ## Run 371 stress 0.08340287 
    ## Run 372 stress 0.07428312 
    ## Run 373 stress 0.07970526 
    ## Run 374 stress 0.08448447 
    ## Run 375 stress 0.06978202 
    ## ... Procrustes: rmse 0.01328673  max resid 0.03339164 
    ## Run 376 stress 0.06942778 
    ## ... Procrustes: rmse 5.554277e-05  max resid 0.0001428196 
    ## ... Similar to previous best
    ## Run 377 stress 0.07970525 
    ## Run 378 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325745  max resid 0.03333308 
    ## Run 379 stress 0.08340288 
    ## Run 380 stress 0.06942776 
    ## ... Procrustes: rmse 5.28461e-06  max resid 1.475819e-05 
    ## ... Similar to previous best
    ## Run 381 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323889  max resid 0.03329619 
    ## Run 382 stress 0.0834029 
    ## Run 383 stress 0.08340297 
    ## Run 384 stress 0.07970525 
    ## Run 385 stress 0.06978189 
    ## ... Procrustes: rmse 0.0132217  max resid 0.03326011 
    ## Run 386 stress 0.06942776 
    ## ... Procrustes: rmse 3.330764e-05  max resid 8.470504e-05 
    ## ... Similar to previous best
    ## Run 387 stress 0.08340289 
    ## Run 388 stress 0.07250814 
    ## Run 389 stress 0.07428319 
    ## Run 390 stress 0.07970525 
    ## Run 391 stress 0.07844957 
    ## Run 392 stress 0.07428317 
    ## Run 393 stress 0.06978196 
    ## ... Procrustes: rmse 0.01317581  max resid 0.0331612 
    ## Run 394 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326825  max resid 0.03335855 
    ## Run 395 stress 0.08340286 
    ## Run 396 stress 0.07428313 
    ## Run 397 stress 0.07970527 
    ## Run 398 stress 0.07428315 
    ## Run 399 stress 0.08340291 
    ## Run 400 stress 0.07428314 
    ## Run 401 stress 0.06942777 
    ## ... Procrustes: rmse 3.996387e-05  max resid 0.000102857 
    ## ... Similar to previous best
    ## Run 402 stress 0.06942776 
    ## ... Procrustes: rmse 3.792754e-06  max resid 8.642949e-06 
    ## ... Similar to previous best
    ## Run 403 stress 0.07428313 
    ## Run 404 stress 0.07428313 
    ## Run 405 stress 0.06942776 
    ## ... Procrustes: rmse 1.812624e-05  max resid 4.779849e-05 
    ## ... Similar to previous best
    ## Run 406 stress 0.07428313 
    ## Run 407 stress 0.07250813 
    ## Run 408 stress 0.07970526 
    ## Run 409 stress 0.06942776 
    ## ... Procrustes: rmse 2.753443e-05  max resid 7.094087e-05 
    ## ... Similar to previous best
    ## Run 410 stress 0.08340287 
    ## Run 411 stress 0.07970526 
    ## Run 412 stress 0.069782 
    ## ... Procrustes: rmse 0.0131644  max resid 0.03313803 
    ## Run 413 stress 0.08340287 
    ## Run 414 stress 0.07428313 
    ## Run 415 stress 0.07844935 
    ## Run 416 stress 0.07428313 
    ## Run 417 stress 0.07428317 
    ## Run 418 stress 0.08340286 
    ## Run 419 stress 0.06942776 
    ## ... Procrustes: rmse 2.83102e-05  max resid 7.301734e-05 
    ## ... Similar to previous best
    ## Run 420 stress 0.08340286 
    ## Run 421 stress 0.06978193 
    ## ... Procrustes: rmse 0.0132549  max resid 0.03332828 
    ## Run 422 stress 0.07428317 
    ## Run 423 stress 0.06942777 
    ## ... Procrustes: rmse 4.158844e-05  max resid 0.0001075187 
    ## ... Similar to previous best
    ## Run 424 stress 0.06978199 
    ## ... Procrustes: rmse 0.01327526  max resid 0.03337375 
    ## Run 425 stress 0.07250813 
    ## Run 426 stress 0.07428315 
    ## Run 427 stress 0.06978191 
    ## ... Procrustes: rmse 0.01323636  max resid 0.03329314 
    ## Run 428 stress 0.06942776 
    ## ... Procrustes: rmse 2.655776e-05  max resid 6.864941e-05 
    ## ... Similar to previous best
    ## Run 429 stress 0.07428322 
    ## Run 430 stress 0.07970525 
    ## Run 431 stress 0.07970527 
    ## Run 432 stress 0.08448437 
    ## Run 433 stress 0.07970526 
    ## Run 434 stress 0.07970526 
    ## Run 435 stress 0.08448437 
    ## Run 436 stress 0.07250812 
    ## Run 437 stress 0.06942776 
    ## ... Procrustes: rmse 3.745109e-06  max resid 9.485468e-06 
    ## ... Similar to previous best
    ## Run 438 stress 0.07250812 
    ## Run 439 stress 0.07428313 
    ## Run 440 stress 0.08340293 
    ## Run 441 stress 0.08340292 
    ## Run 442 stress 0.06942776 
    ## ... Procrustes: rmse 2.885408e-06  max resid 7.531483e-06 
    ## ... Similar to previous best
    ## Run 443 stress 0.06978192 
    ## ... Procrustes: rmse 0.01319509  max resid 0.0332059 
    ## Run 444 stress 0.0834029 
    ## Run 445 stress 0.08340287 
    ## Run 446 stress 0.07428313 
    ## Run 447 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319947  max resid 0.03321382 
    ## Run 448 stress 0.07428317 
    ## Run 449 stress 0.06942777 
    ## ... Procrustes: rmse 5.051902e-05  max resid 0.0001299649 
    ## ... Similar to previous best
    ## Run 450 stress 0.07970525 
    ## Run 451 stress 0.06942776 
    ## ... Procrustes: rmse 2.17787e-05  max resid 5.616861e-05 
    ## ... Similar to previous best
    ## Run 452 stress 0.06942776 
    ## ... Procrustes: rmse 1.665977e-06  max resid 4.237387e-06 
    ## ... Similar to previous best
    ## Run 453 stress 0.08340287 
    ## Run 454 stress 0.07428315 
    ## Run 455 stress 0.08340286 
    ## Run 456 stress 0.07428313 
    ## Run 457 stress 0.06978191 
    ## ... Procrustes: rmse 0.01319507  max resid 0.03320404 
    ## Run 458 stress 0.06978191 
    ## ... Procrustes: rmse 0.0132428  max resid 0.03330339 
    ## Run 459 stress 0.07428313 
    ## Run 460 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323828  max resid 0.03329573 
    ## Run 461 stress 0.06978195 
    ## ... Procrustes: rmse 0.01325792  max resid 0.03333719 
    ## Run 462 stress 0.06978201 
    ## ... Procrustes: rmse 0.01326751  max resid 0.03336095 
    ## Run 463 stress 0.06942777 
    ## ... Procrustes: rmse 5.010497e-05  max resid 0.0001298161 
    ## ... Similar to previous best
    ## Run 464 stress 0.06978196 
    ## ... Procrustes: rmse 0.01326153  max resid 0.03334615 
    ## Run 465 stress 0.06978191 
    ## ... Procrustes: rmse 0.01320028  max resid 0.03320678 
    ## Run 466 stress 0.07428317 
    ## Run 467 stress 0.06942776 
    ## ... Procrustes: rmse 3.127264e-05  max resid 8.085576e-05 
    ## ... Similar to previous best
    ## Run 468 stress 0.08340287 
    ## Run 469 stress 0.06978193 
    ## ... Procrustes: rmse 0.01322777  max resid 0.03326806 
    ## Run 470 stress 0.06978199 
    ## ... Procrustes: rmse 0.01326344  max resid 0.03334298 
    ## Run 471 stress 0.07250812 
    ## Run 472 stress 0.08340286 
    ## Run 473 stress 0.07250812 
    ## Run 474 stress 0.07428312 
    ## Run 475 stress 0.08448441 
    ## Run 476 stress 0.07970525 
    ## Run 477 stress 0.06942776 
    ## ... Procrustes: rmse 7.446369e-06  max resid 2.246962e-05 
    ## ... Similar to previous best
    ## Run 478 stress 0.0697819 
    ## ... Procrustes: rmse 0.01323748  max resid 0.0332921 
    ## Run 479 stress 0.07970525 
    ## Run 480 stress 0.07970526 
    ## Run 481 stress 0.07970526 
    ## Run 482 stress 0.08340287 
    ## Run 483 stress 0.07428316 
    ## Run 484 stress 0.06978194 
    ## ... Procrustes: rmse 0.01325394  max resid 0.03332979 
    ## Run 485 stress 0.08340287 
    ## Run 486 stress 0.08340289 
    ## Run 487 stress 0.07970526 
    ## Run 488 stress 0.06978191 
    ## ... Procrustes: rmse 0.01324342  max resid 0.03330838 
    ## Run 489 stress 0.06978195 
    ## ... Procrustes: rmse 0.01326441  max resid 0.03334606 
    ## Run 490 stress 0.07844917 
    ## Run 491 stress 0.07428313 
    ## Run 492 stress 0.07250813 
    ## Run 493 stress 0.07970525 
    ## Run 494 stress 0.08340292 
    ## Run 495 stress 0.07970525 
    ## Run 496 stress 0.06942776 
    ## ... Procrustes: rmse 1.885517e-05  max resid 4.849736e-05 
    ## ... Similar to previous best
    ## Run 497 stress 0.08448439 
    ## Run 498 stress 0.069782 
    ## ... Procrustes: rmse 0.01328348  max resid 0.03338668 
    ## Run 499 stress 0.08448437 
    ## Run 500 stress 0.06942776 
    ## ... Procrustes: rmse 1.367752e-05  max resid 3.576071e-05 
    ## ... Similar to previous best
    ## *** Best solution repeated 45 times

``` r
# Mixed lakes
SD_beta_geo_M_NMDS <- metaMDS(SD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Run 0 stress 0.001279543 
    ## Run 1 stress 0.1491957 
    ## Run 2 stress 0.001218741 
    ## ... New best solution
    ## ... Procrustes: rmse 0.001053994  max resid 0.001462106 
    ## ... Similar to previous best
    ## Run 3 stress 9.999269e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.04072166  max resid 0.06829503 
    ## Run 4 stress 0.1990774 
    ## Run 5 stress 9.667907e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.006866224  max resid 0.009442798 
    ## Run 6 stress 0.001261278 
    ## Run 7 stress 8.948694e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004376193  max resid 0.0008615734 
    ## ... Similar to previous best
    ## Run 8 stress 8.901942e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001906923  max resid 0.0003924015 
    ## ... Similar to previous best
    ## Run 9 stress 0.001287533 
    ## Run 10 stress 0.2842805 
    ## Run 11 stress 0.001457894 
    ## Run 12 stress 0.0004638115 
    ## ... Procrustes: rmse 0.01547  max resid 0.0212375 
    ## Run 13 stress 9.719129e-05 
    ## ... Procrustes: rmse 0.0002016521  max resid 0.0004157812 
    ## ... Similar to previous best
    ## Run 14 stress 0.001519844 
    ## Run 15 stress 8.780869e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001783891  max resid 0.0002960851 
    ## ... Similar to previous best
    ## Run 16 stress 9.185653e-05 
    ## ... Procrustes: rmse 0.0001756944  max resid 0.0002926951 
    ## ... Similar to previous best
    ## Run 17 stress 0.00132723 
    ## Run 18 stress 0.001358133 
    ## Run 19 stress 0.001076155 
    ## Run 20 stress 0.001211602 
    ## Run 21 stress 0.0003035227 
    ## ... Procrustes: rmse 0.02099097  max resid 0.02905042 
    ## Run 22 stress 8.043144e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001695742  max resid 0.0003450236 
    ## ... Similar to previous best
    ## Run 23 stress 0.3078457 
    ## Run 24 stress 0.001087908 
    ## Run 25 stress 0.2842805 
    ## Run 26 stress 0.001200269 
    ## Run 27 stress 0.1491957 
    ## Run 28 stress 9.385299e-05 
    ## ... Procrustes: rmse 0.0001863117  max resid 0.0004441314 
    ## ... Similar to previous best
    ## Run 29 stress 9.2726e-05 
    ## ... Procrustes: rmse 0.0001608044  max resid 0.000349975 
    ## ... Similar to previous best
    ## Run 30 stress 0.0004808608 
    ## ... Procrustes: rmse 0.01608484  max resid 0.02254599 
    ## Run 31 stress 0.001345244 
    ## Run 32 stress 0.1491957 
    ## Run 33 stress 0.1776022 
    ## Run 34 stress 9.657009e-05 
    ## ... Procrustes: rmse 0.00015385  max resid 0.0003205403 
    ## ... Similar to previous best
    ## Run 35 stress 0.001236825 
    ## Run 36 stress 9.244938e-05 
    ## ... Procrustes: rmse 0.0002192302  max resid 0.0004010798 
    ## ... Similar to previous best
    ## Run 37 stress 0.001223283 
    ## Run 38 stress 0.0004436451 
    ## ... Procrustes: rmse 0.01544585  max resid 0.02166494 
    ## Run 39 stress 0.001139842 
    ## Run 40 stress 9.211743e-05 
    ## ... Procrustes: rmse 0.0001689129  max resid 0.0003773101 
    ## ... Similar to previous best
    ## Run 41 stress 0.0005081309 
    ## ... Procrustes: rmse 0.01653472  max resid 0.02316669 
    ## Run 42 stress 8.649189e-05 
    ## ... Procrustes: rmse 0.0001991487  max resid 0.0003610085 
    ## ... Similar to previous best
    ## Run 43 stress 0.001231825 
    ## Run 44 stress 9.080918e-05 
    ## ... Procrustes: rmse 0.0001723464  max resid 0.000406786 
    ## ... Similar to previous best
    ## Run 45 stress 0.0004777707 
    ## ... Procrustes: rmse 0.01602593  max resid 0.02246465 
    ## Run 46 stress 0.001367154 
    ## Run 47 stress 8.517879e-05 
    ## ... Procrustes: rmse 0.005317296  max resid 0.007699636 
    ## Run 48 stress 0.0008987556 
    ## Run 49 stress 9.765653e-05 
    ## ... Procrustes: rmse 0.0001689418  max resid 0.0003944306 
    ## ... Similar to previous best
    ## Run 50 stress 9.865783e-05 
    ## ... Procrustes: rmse 0.007749005  max resid 0.01069906 
    ## Run 51 stress 8.531039e-05 
    ## ... Procrustes: rmse 0.0001661518  max resid 0.0003657634 
    ## ... Similar to previous best
    ## Run 52 stress 7.469818e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0005790532  max resid 0.0008843069 
    ## ... Similar to previous best
    ## Run 53 stress 9.747719e-05 
    ## ... Procrustes: rmse 0.0005582393  max resid 0.0009813461 
    ## ... Similar to previous best
    ## Run 54 stress 0.2842805 
    ## Run 55 stress 9.701198e-05 
    ## ... Procrustes: rmse 0.01104952  max resid 0.0162259 
    ## Run 56 stress 0.001270319 
    ## Run 57 stress 0.0004227437 
    ## ... Procrustes: rmse 0.01454377  max resid 0.02003923 
    ## Run 58 stress 9.037417e-05 
    ## ... Procrustes: rmse 0.0004761761  max resid 0.0007440503 
    ## ... Similar to previous best
    ## Run 59 stress 0.001354091 
    ## Run 60 stress 8.489741e-05 
    ## ... Procrustes: rmse 0.004575529  max resid 0.007265499 
    ## ... Similar to previous best
    ## Run 61 stress 9.895097e-05 
    ## ... Procrustes: rmse 0.0004818094  max resid 0.0008086357 
    ## ... Similar to previous best
    ## Run 62 stress 0.001347223 
    ## Run 63 stress 8.444915e-05 
    ## ... Procrustes: rmse 0.0002805056  max resid 0.0003843364 
    ## ... Similar to previous best
    ## Run 64 stress 9.116429e-05 
    ## ... Procrustes: rmse 0.003714038  max resid 0.005121707 
    ## ... Similar to previous best
    ## Run 65 stress 7.062994e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0004522729  max resid 0.0006225669 
    ## ... Similar to previous best
    ## Run 66 stress 0.0004018532 
    ## ... Procrustes: rmse 0.01461014  max resid 0.02008142 
    ## Run 67 stress 0.1776011 
    ## Run 68 stress 9.841277e-05 
    ## ... Procrustes: rmse 8.62835e-05  max resid 0.0001079036 
    ## ... Similar to previous best
    ## Run 69 stress 0.001318627 
    ## Run 70 stress 0.0003500502 
    ## ... Procrustes: rmse 0.01359211  max resid 0.01867669 
    ## Run 71 stress 0.1491957 
    ## Run 72 stress 7.930896e-05 
    ## ... Procrustes: rmse 0.001653958  max resid 0.002274623 
    ## ... Similar to previous best
    ## Run 73 stress 0.1990774 
    ## Run 74 stress 0.0004867769 
    ## ... Procrustes: rmse 0.01608901  max resid 0.02212213 
    ## Run 75 stress 0.0004953871 
    ## ... Procrustes: rmse 0.01623131  max resid 0.02231897 
    ## Run 76 stress 0.001298292 
    ## Run 77 stress 0.1491957 
    ## Run 78 stress 0.0001052608 
    ## ... Procrustes: rmse 0.007429144  max resid 0.01017222 
    ## Run 79 stress 0.001361678 
    ## Run 80 stress 0.0007346707 
    ## Run 81 stress 9.677143e-05 
    ## ... Procrustes: rmse 0.0001202394  max resid 0.0002300036 
    ## ... Similar to previous best
    ## Run 82 stress 0.1491957 
    ## Run 83 stress 0.001593937 
    ## Run 84 stress 0.1776017 
    ## Run 85 stress 0.0013903 
    ## Run 86 stress 0.0004558502 
    ## ... Procrustes: rmse 0.01551835  max resid 0.02133287 
    ## Run 87 stress 8.49777e-05 
    ## ... Procrustes: rmse 0.004998445  max resid 0.006902128 
    ## ... Similar to previous best
    ## Run 88 stress 9.968861e-05 
    ## ... Procrustes: rmse 8.238975e-05  max resid 9.849851e-05 
    ## ... Similar to previous best
    ## Run 89 stress 9.443502e-05 
    ## ... Procrustes: rmse 0.0002182766  max resid 0.0004237994 
    ## ... Similar to previous best
    ## Run 90 stress 0.0004625656 
    ## ... Procrustes: rmse 0.01568063  max resid 0.02155903 
    ## Run 91 stress 0.2528294 
    ## Run 92 stress 0.001458243 
    ## Run 93 stress 9.204503e-05 
    ## ... Procrustes: rmse 0.0002198576  max resid 0.0004199945 
    ## ... Similar to previous best
    ## Run 94 stress 0.177601 
    ## Run 95 stress 9.947466e-05 
    ## ... Procrustes: rmse 0.0001266081  max resid 0.0002073288 
    ## ... Similar to previous best
    ## Run 96 stress 0.001340667 
    ## Run 97 stress 9.635189e-05 
    ## ... Procrustes: rmse 0.0002740211  max resid 0.000474886 
    ## ... Similar to previous best
    ## Run 98 stress 0.001309385 
    ## Run 99 stress 9.434951e-05 
    ## ... Procrustes: rmse 0.0001206035  max resid 0.0001947503 
    ## ... Similar to previous best
    ## Run 100 stress 0.1990774 
    ## Run 101 stress 0.1491957 
    ## Run 102 stress 8.801018e-05 
    ## ... Procrustes: rmse 0.0001368054  max resid 0.0002223012 
    ## ... Similar to previous best
    ## Run 103 stress 0.2842805 
    ## Run 104 stress 8.438155e-05 
    ## ... Procrustes: rmse 0.0005245851  max resid 0.0006815316 
    ## ... Similar to previous best
    ## Run 105 stress 0.001334606 
    ## Run 106 stress 7.147873e-05 
    ## ... Procrustes: rmse 8.569711e-05  max resid 0.0001452508 
    ## ... Similar to previous best
    ## Run 107 stress 8.741402e-05 
    ## ... Procrustes: rmse 0.0001123529  max resid 0.0002121932 
    ## ... Similar to previous best
    ## Run 108 stress 0.0004912234 
    ## ... Procrustes: rmse 0.01615871  max resid 0.02222191 
    ## Run 109 stress 0.001233165 
    ## Run 110 stress 0.1491957 
    ## Run 111 stress 9.635422e-05 
    ## ... Procrustes: rmse 0.000221161  max resid 0.0004169338 
    ## ... Similar to previous best
    ## Run 112 stress 0.0006953471 
    ## Run 113 stress 0.0005171845 
    ## ... Procrustes: rmse 0.0165883  max resid 0.02281078 
    ## Run 114 stress 0.2842805 
    ## Run 115 stress 0.00123362 
    ## Run 116 stress 0.1491957 
    ## Run 117 stress 0.001224831 
    ## Run 118 stress 0.1990774 
    ## Run 119 stress 0.0002000767 
    ## ... Procrustes: rmse 0.0170901  max resid 0.02339503 
    ## Run 120 stress 0.001280686 
    ## Run 121 stress 8.943081e-05 
    ## ... Procrustes: rmse 5.655292e-05  max resid 8.43062e-05 
    ## ... Similar to previous best
    ## Run 122 stress 0.1491957 
    ## Run 123 stress 8.870071e-05 
    ## ... Procrustes: rmse 0.0001664538  max resid 0.0002605248 
    ## ... Similar to previous best
    ## Run 124 stress 0.2373687 
    ## Run 125 stress 0.001446055 
    ## Run 126 stress 0.001111652 
    ## Run 127 stress 8.664289e-05 
    ## ... Procrustes: rmse 0.0001726678  max resid 0.0002422942 
    ## ... Similar to previous best
    ## Run 128 stress 0.0004771386 
    ## ... Procrustes: rmse 0.01590416  max resid 0.02186483 
    ## Run 129 stress 0.001333485 
    ## Run 130 stress 0.001352339 
    ## Run 131 stress 9.847006e-05 
    ## ... Procrustes: rmse 0.01122281  max resid 0.01531343 
    ## Run 132 stress 0.1491957 
    ## Run 133 stress 0.0004755115 
    ## ... Procrustes: rmse 0.01578291  max resid 0.02169876 
    ## Run 134 stress 0.0006777516 
    ## Run 135 stress 0.001319441 
    ## Run 136 stress 8.756106e-05 
    ## ... Procrustes: rmse 0.0002117943  max resid 0.0004086474 
    ## ... Similar to previous best
    ## Run 137 stress 0.1776018 
    ## Run 138 stress 9.465124e-05 
    ## ... Procrustes: rmse 0.000220887  max resid 0.0004152475 
    ## ... Similar to previous best
    ## Run 139 stress 0.001333614 
    ## Run 140 stress 0.1776017 
    ## Run 141 stress 0.2520602 
    ## Run 142 stress 9.326097e-05 
    ## ... Procrustes: rmse 0.000213634  max resid 0.0003935037 
    ## ... Similar to previous best
    ## Run 143 stress 0.1990774 
    ## Run 144 stress 0.1491957 
    ## Run 145 stress 0.00121948 
    ## Run 146 stress 0.3075147 
    ## Run 147 stress 9.781337e-05 
    ## ... Procrustes: rmse 0.0002134818  max resid 0.0004165865 
    ## ... Similar to previous best
    ## Run 148 stress 0.1491957 
    ## Run 149 stress 0.0008503953 
    ## Run 150 stress 0.001279523 
    ## Run 151 stress 9.40069e-05 
    ## ... Procrustes: rmse 0.0002230237  max resid 0.0004221593 
    ## ... Similar to previous best
    ## Run 152 stress 9.264414e-05 
    ## ... Procrustes: rmse 0.0001048531  max resid 0.0001846087 
    ## ... Similar to previous best
    ## Run 153 stress 0.0004147262 
    ## ... Procrustes: rmse 0.0148439  max resid 0.02040478 
    ## Run 154 stress 0.0004840566 
    ## ... Procrustes: rmse 0.01604099  max resid 0.02205738 
    ## Run 155 stress 0.001442735 
    ## Run 156 stress 9.122835e-05 
    ## ... Procrustes: rmse 0.0001104125  max resid 0.0002039251 
    ## ... Similar to previous best
    ## Run 157 stress 0.0001330276 
    ## ... Procrustes: rmse 0.01393707  max resid 0.01903728 
    ## Run 158 stress 0.0012002 
    ## Run 159 stress 0.001216806 
    ## Run 160 stress 0.001320095 
    ## Run 161 stress 0.1491957 
    ## Run 162 stress 0.284852 
    ## Run 163 stress 0.001376701 
    ## Run 164 stress 0.0008297168 
    ## Run 165 stress 8.098201e-05 
    ## ... Procrustes: rmse 0.001080168  max resid 0.001511264 
    ## ... Similar to previous best
    ## Run 166 stress 8.209071e-05 
    ## ... Procrustes: rmse 0.0001531474  max resid 0.0002398496 
    ## ... Similar to previous best
    ## Run 167 stress 0.001239397 
    ## Run 168 stress 8.635766e-05 
    ## ... Procrustes: rmse 0.002915433  max resid 0.003965596 
    ## ... Similar to previous best
    ## Run 169 stress 0.001123425 
    ## Run 170 stress 0.001352447 
    ## Run 171 stress 0.1491957 
    ## Run 172 stress 8.82917e-05 
    ## ... Procrustes: rmse 0.0001575443  max resid 0.0002404765 
    ## ... Similar to previous best
    ## Run 173 stress 0.000861039 
    ## Run 174 stress 9.601874e-05 
    ## ... Procrustes: rmse 0.0003476249  max resid 0.0005404755 
    ## ... Similar to previous best
    ## Run 175 stress 0.000972744 
    ## Run 176 stress 9.867693e-05 
    ## ... Procrustes: rmse 0.0001157502  max resid 0.0001979076 
    ## ... Similar to previous best
    ## Run 177 stress 0.001269274 
    ## Run 178 stress 0.1491957 
    ## Run 179 stress 0.001377258 
    ## Run 180 stress 9.360784e-05 
    ## ... Procrustes: rmse 0.0002197629  max resid 0.0004147672 
    ## ... Similar to previous best
    ## Run 181 stress 0.2509018 
    ## Run 182 stress 0.001350625 
    ## Run 183 stress 0.0004538307 
    ## ... Procrustes: rmse 0.02570851  max resid 0.03532083 
    ## Run 184 stress 9.809457e-05 
    ## ... Procrustes: rmse 0.0002112813  max resid 0.0003840254 
    ## ... Similar to previous best
    ## Run 185 stress 0.3083098 
    ## Run 186 stress 0.00130631 
    ## Run 187 stress 0.001331023 
    ## Run 188 stress 9.464845e-05 
    ## ... Procrustes: rmse 0.0001186901  max resid 0.0002296925 
    ## ... Similar to previous best
    ## Run 189 stress 0.0004176641 
    ## ... Procrustes: rmse 0.01489722  max resid 0.02047746 
    ## Run 190 stress 9.105515e-05 
    ## ... Procrustes: rmse 0.0001152717  max resid 0.0001863011 
    ## ... Similar to previous best
    ## Run 191 stress 9.692885e-05 
    ## ... Procrustes: rmse 8.016593e-05  max resid 0.0001351596 
    ## ... Similar to previous best
    ## Run 192 stress 0.001337195 
    ## Run 193 stress 0.001204022 
    ## Run 194 stress 0.001270927 
    ## Run 195 stress 0.0004709176 
    ## ... Procrustes: rmse 0.0158241  max resid 0.02175728 
    ## Run 196 stress 0.000486647 
    ## ... Procrustes: rmse 0.01608764  max resid 0.02212007 
    ## Run 197 stress 0.0013263 
    ## Run 198 stress 8.672259e-05 
    ## ... Procrustes: rmse 0.0002037522  max resid 0.0003899072 
    ## ... Similar to previous best
    ## Run 199 stress 0.0008202471 
    ## Run 200 stress 8.54923e-05 
    ## ... Procrustes: rmse 0.0001157436  max resid 0.0001965767 
    ## ... Similar to previous best
    ## Run 201 stress 0.00110858 
    ## Run 202 stress 7.902954e-05 
    ## ... Procrustes: rmse 0.0001227548  max resid 0.0002160299 
    ## ... Similar to previous best
    ## Run 203 stress 0.00134781 
    ## Run 204 stress 0.1491957 
    ## Run 205 stress 8.99371e-05 
    ## ... Procrustes: rmse 0.0001562722  max resid 0.0002180342 
    ## ... Similar to previous best
    ## Run 206 stress 9.235414e-05 
    ## ... Procrustes: rmse 9.184624e-05  max resid 0.0001607269 
    ## ... Similar to previous best
    ## Run 207 stress 0.3083098 
    ## Run 208 stress 9.595367e-05 
    ## ... Procrustes: rmse 0.000183877  max resid 0.0003968569 
    ## ... Similar to previous best
    ## Run 209 stress 9.870964e-05 
    ## ... Procrustes: rmse 9.470947e-05  max resid 0.0001226529 
    ## ... Similar to previous best
    ## Run 210 stress 9.178937e-05 
    ## ... Procrustes: rmse 0.0002179879  max resid 0.0004214303 
    ## ... Similar to previous best
    ## Run 211 stress 0.001425869 
    ## Run 212 stress 9.123858e-05 
    ## ... Procrustes: rmse 0.0002155901  max resid 0.0004234055 
    ## ... Similar to previous best
    ## Run 213 stress 0.1776011 
    ## Run 214 stress 7.222108e-05 
    ## ... Procrustes: rmse 0.0005021492  max resid 0.0007605577 
    ## ... Similar to previous best
    ## Run 215 stress 0.0004823334 
    ## ... Procrustes: rmse 0.01601637  max resid 0.0220218 
    ## Run 216 stress 0.2852157 
    ## Run 217 stress 0.1491957 
    ## Run 218 stress 0.2568281 
    ## Run 219 stress 0.2842805 
    ## Run 220 stress 6.147011e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0001621372  max resid 0.000312888 
    ## ... Similar to previous best
    ## Run 221 stress 5.041931e-05 
    ## ... New best solution
    ## ... Procrustes: rmse 0.0002094894  max resid 0.0003703157 
    ## ... Similar to previous best
    ## Run 222 stress 6.200413e-05 
    ## ... Procrustes: rmse 0.0002866394  max resid 0.0004197411 
    ## ... Similar to previous best
    ## Run 223 stress 0.001037397 
    ## Run 224 stress 9.755136e-05 
    ## ... Procrustes: rmse 0.0002429539  max resid 0.0005254758 
    ## ... Similar to previous best
    ## Run 225 stress 0.1776011 
    ## Run 226 stress 0.001217543 
    ## Run 227 stress 0.001509814 
    ## Run 228 stress 8.711781e-05 
    ## ... Procrustes: rmse 0.002707962  max resid 0.004212738 
    ## ... Similar to previous best
    ## Run 229 stress 0.001249854 
    ## Run 230 stress 0.0003865706 
    ## ... Procrustes: rmse 0.01437439  max resid 0.02029986 
    ## Run 231 stress 8.115234e-05 
    ## ... Procrustes: rmse 0.000242135  max resid 0.0004618326 
    ## ... Similar to previous best
    ## Run 232 stress 0.001249357 
    ## Run 233 stress 0.00137343 
    ## Run 234 stress 0.0004491463 
    ## ... Procrustes: rmse 0.0154979  max resid 0.02185024 
    ## Run 235 stress 0.00126433 
    ## Run 236 stress 9.12822e-05 
    ## ... Procrustes: rmse 0.0001583331  max resid 0.0002971204 
    ## ... Similar to previous best
    ## Run 237 stress 0.0004866116 
    ## ... Procrustes: rmse 0.01613327  max resid 0.02272589 
    ## Run 238 stress 9.344639e-05 
    ## ... Procrustes: rmse 0.0002755148  max resid 0.000463375 
    ## ... Similar to previous best
    ## Run 239 stress 6.541712e-05 
    ## ... Procrustes: rmse 0.0005796548  max resid 0.0007870097 
    ## ... Similar to previous best
    ## Run 240 stress 8.758082e-05 
    ## ... Procrustes: rmse 0.0002797336  max resid 0.0004030968 
    ## ... Similar to previous best
    ## Run 241 stress 0.001272062 
    ## Run 242 stress 0.001032628 
    ## Run 243 stress 0.001347989 
    ## Run 244 stress 9.932686e-05 
    ## ... Procrustes: rmse 0.0002544526  max resid 0.0003828441 
    ## ... Similar to previous best
    ## Run 245 stress 8.84993e-05 
    ## ... Procrustes: rmse 0.0002340139  max resid 0.0004677089 
    ## ... Similar to previous best
    ## Run 246 stress 0.001215423 
    ## Run 247 stress 9.880442e-05 
    ## ... Procrustes: rmse 0.0002667445  max resid 0.0004872794 
    ## ... Similar to previous best
    ## Run 248 stress 0.001172285 
    ## Run 249 stress 0.001304138 
    ## Run 250 stress 9.968464e-05 
    ## ... Procrustes: rmse 0.0001871731  max resid 0.0002804137 
    ## ... Similar to previous best
    ## Run 251 stress 0.2440272 
    ## Run 252 stress 0.001298628 
    ## Run 253 stress 0.1990774 
    ## Run 254 stress 9.720506e-05 
    ## ... Procrustes: rmse 0.003219582  max resid 0.004442611 
    ## ... Similar to previous best
    ## Run 255 stress 0.0004931938 
    ## ... Procrustes: rmse 0.01624199  max resid 0.02287708 
    ## Run 256 stress 0.001230606 
    ## Run 257 stress 8.708241e-05 
    ## ... Procrustes: rmse 0.0003011134  max resid 0.0005492932 
    ## ... Similar to previous best
    ## Run 258 stress 6.348359e-05 
    ## ... Procrustes: rmse 0.004793631  max resid 0.006647919 
    ## ... Similar to previous best
    ## Run 259 stress 0.0004794273 
    ## ... Procrustes: rmse 0.01601195  max resid 0.02255958 
    ## Run 260 stress 8.741231e-05 
    ## ... Procrustes: rmse 0.000167808  max resid 0.0003189806 
    ## ... Similar to previous best
    ## Run 261 stress 0.1990774 
    ## Run 262 stress 9.264173e-05 
    ## ... Procrustes: rmse 0.0002297491  max resid 0.0005065605 
    ## ... Similar to previous best
    ## Run 263 stress 0.001405824 
    ## Run 264 stress 0.001155484 
    ## Run 265 stress 9.827401e-05 
    ## ... Procrustes: rmse 0.0001576761  max resid 0.0003008276 
    ## ... Similar to previous best
    ## Run 266 stress 0.3083098 
    ## Run 267 stress 0.001282958 
    ## Run 268 stress 0.0004742959 
    ## ... Procrustes: rmse 0.0159272  max resid 0.02244177 
    ## Run 269 stress 0.0004673772 
    ## ... Procrustes: rmse 0.01580995  max resid 0.02228011 
    ## Run 270 stress 0.001190194 
    ## Run 271 stress 0.001375055 
    ## Run 272 stress 0.001125814 
    ## Run 273 stress 0.00138562 
    ## Run 274 stress 0.001109459 
    ## Run 275 stress 8.482036e-05 
    ## ... Procrustes: rmse 0.0002418462  max resid 0.0005096434 
    ## ... Similar to previous best
    ## Run 276 stress 9.061393e-05 
    ## ... Procrustes: rmse 0.0002610561  max resid 0.0004792887 
    ## ... Similar to previous best
    ## Run 277 stress 0.001358225 
    ## Run 278 stress 0.1990774 
    ## Run 279 stress 9.568148e-05 
    ## ... Procrustes: rmse 0.0001964577  max resid 0.0003861338 
    ## ... Similar to previous best
    ## Run 280 stress 0.2848525 
    ## Run 281 stress 0.0004164568 
    ## ... Procrustes: rmse 0.01492167  max resid 0.02105455 
    ## Run 282 stress 9.706128e-05 
    ## ... Procrustes: rmse 0.0002613713  max resid 0.0004859972 
    ## ... Similar to previous best
    ## Run 283 stress 0.001213275 
    ## Run 284 stress 0.0002286402 
    ## ... Procrustes: rmse 0.01802528  max resid 0.02486068 
    ## Run 285 stress 0.001362465 
    ## Run 286 stress 0.001441857 
    ## Run 287 stress 0.001424102 
    ## Run 288 stress 9.823692e-05 
    ## ... Procrustes: rmse 0.0002622969  max resid 0.0004960745 
    ## ... Similar to previous best
    ## Run 289 stress 0.001306009 
    ## Run 290 stress 0.0004840643 
    ## ... Procrustes: rmse 0.01607588  max resid 0.02265157 
    ## Run 291 stress 9.786079e-05 
    ## ... Procrustes: rmse 0.006822189  max resid 0.009886099 
    ## Run 292 stress 0.00121255 
    ## Run 293 stress 0.1491957 
    ## Run 294 stress 0.0006835131 
    ## Run 295 stress 0.1491957 
    ## Run 296 stress 0.2509018 
    ## Run 297 stress 0.0004628722 
    ## ... Procrustes: rmse 0.01573354  max resid 0.02217548 
    ## Run 298 stress 0.001370081 
    ## Run 299 stress 9.180411e-05 
    ## ... Procrustes: rmse 0.0002031259  max resid 0.000342609 
    ## ... Similar to previous best
    ## Run 300 stress 0.00120182 
    ## Run 301 stress 0.001421482 
    ## Run 302 stress 9.720778e-05 
    ## ... Procrustes: rmse 0.0002569807  max resid 0.0004086368 
    ## ... Similar to previous best
    ## Run 303 stress 8.981199e-05 
    ## ... Procrustes: rmse 0.0001616866  max resid 0.0003076665 
    ## ... Similar to previous best
    ## Run 304 stress 9.499283e-05 
    ## ... Procrustes: rmse 0.0001577684  max resid 0.000298431 
    ## ... Similar to previous best
    ## Run 305 stress 9.864192e-05 
    ## ... Procrustes: rmse 0.0002857663  max resid 0.0004885215 
    ## ... Similar to previous best
    ## Run 306 stress 0.001232174 
    ## Run 307 stress 9.536479e-05 
    ## ... Procrustes: rmse 0.0001586225  max resid 0.0003028326 
    ## ... Similar to previous best
    ## Run 308 stress 0.001168881 
    ## Run 309 stress 0.0009894307 
    ## Run 310 stress 0.1491957 
    ## Run 311 stress 0.1990774 
    ## Run 312 stress 0.001228539 
    ## Run 313 stress 0.001364922 
    ## Run 314 stress 0.001205178 
    ## Run 315 stress 0.1491957 
    ## Run 316 stress 0.0004496708 
    ## ... Procrustes: rmse 0.01550583  max resid 0.02186057 
    ## Run 317 stress 0.001321403 
    ## Run 318 stress 0.1491957 
    ## Run 319 stress 0.001221492 
    ## Run 320 stress 0.0004914331 
    ## ... Procrustes: rmse 0.01621361  max resid 0.02283739 
    ## Run 321 stress 8.901991e-05 
    ## ... Procrustes: rmse 0.0002268325  max resid 0.0004856272 
    ## ... Similar to previous best
    ## Run 322 stress 0.1990774 
    ## Run 323 stress 0.1491957 
    ## Run 324 stress 0.1491957 
    ## Run 325 stress 0.2509018 
    ## Run 326 stress 9.576073e-05 
    ## ... Procrustes: rmse 0.0002580089  max resid 0.000515636 
    ## ... Similar to previous best
    ## Run 327 stress 0.0004434122 
    ## ... Procrustes: rmse 0.01539761  max resid 0.0217115 
    ## Run 328 stress 0.001225938 
    ## Run 329 stress 9.054812e-05 
    ## ... Procrustes: rmse 0.0002429706  max resid 0.000501458 
    ## ... Similar to previous best
    ## Run 330 stress 0.1491957 
    ## Run 331 stress 9.889861e-05 
    ## ... Procrustes: rmse 0.0002481649  max resid 0.0005149873 
    ## ... Similar to previous best
    ## Run 332 stress 9.196756e-05 
    ## ... Procrustes: rmse 0.0002551388  max resid 0.0004720053 
    ## ... Similar to previous best
    ## Run 333 stress 0.001263244 
    ## Run 334 stress 0.001253194 
    ## Run 335 stress 0.0003748015 
    ## ... Procrustes: rmse 0.01415335  max resid 0.01999516 
    ## Run 336 stress 0.2842805 
    ## Run 337 stress 0.001393239 
    ## Run 338 stress 0.0004831239 
    ## ... Procrustes: rmse 0.0160586  max resid 0.02262842 
    ## Run 339 stress 8.700678e-05 
    ## ... Procrustes: rmse 0.004959367  max resid 0.006890546 
    ## ... Similar to previous best
    ## Run 340 stress 0.2852157 
    ## Run 341 stress 0.00136327 
    ## Run 342 stress 0.000195985 
    ## ... Procrustes: rmse 0.01022092  max resid 0.01456777 
    ## Run 343 stress 8.477935e-05 
    ## ... Procrustes: rmse 0.0001970645  max resid 0.0003905701 
    ## ... Similar to previous best
    ## Run 344 stress 0.1990774 
    ## Run 345 stress 9.164217e-05 
    ## ... Procrustes: rmse 0.0001572527  max resid 0.0002903858 
    ## ... Similar to previous best
    ## Run 346 stress 0.001065871 
    ## Run 347 stress 0.1990774 
    ## Run 348 stress 9.00432e-05 
    ## ... Procrustes: rmse 0.0001863714  max resid 0.0002737835 
    ## ... Similar to previous best
    ## Run 349 stress 0.001040065 
    ## Run 350 stress 0.0004726531 
    ## ... Procrustes: rmse 0.01589975  max resid 0.02240519 
    ## Run 351 stress 9.953587e-05 
    ## ... Procrustes: rmse 0.0002447895  max resid 0.000523917 
    ## ... Similar to previous best
    ## Run 352 stress 9.618533e-05 
    ## ... Procrustes: rmse 0.001481954  max resid 0.002109023 
    ## ... Similar to previous best
    ## Run 353 stress 0.001509279 
    ## Run 354 stress 0.0012765 
    ## Run 355 stress 0.001262668 
    ## Run 356 stress 0.001155332 
    ## Run 357 stress 0.001200944 
    ## Run 358 stress 0.0004736416 
    ## ... Procrustes: rmse 0.01591646  max resid 0.02242756 
    ## Run 359 stress 0.00120118 
    ## Run 360 stress 9.356681e-05 
    ## ... Procrustes: rmse 0.0002071815  max resid 0.0003954502 
    ## ... Similar to previous best
    ## Run 361 stress 0.1990774 
    ## Run 362 stress 0.001177947 
    ## Run 363 stress 8.523175e-05 
    ## ... Procrustes: rmse 0.001247241  max resid 0.002145217 
    ## ... Similar to previous best
    ## Run 364 stress 0.0006770365 
    ## Run 365 stress 0.1776011 
    ## Run 366 stress 9.689988e-05 
    ## ... Procrustes: rmse 0.0002576735  max resid 0.000474568 
    ## ... Similar to previous best
    ## Run 367 stress 0.001336704 
    ## Run 368 stress 0.001261147 
    ## Run 369 stress 0.001331805 
    ## Run 370 stress 0.0009622561 
    ## Run 371 stress 9.907368e-05 
    ## ... Procrustes: rmse 0.0002680222  max resid 0.0005083033 
    ## ... Similar to previous best
    ## Run 372 stress 9.079256e-05 
    ## ... Procrustes: rmse 0.0002353391  max resid 0.0005099819 
    ## ... Similar to previous best
    ## Run 373 stress 0.1491957 
    ## Run 374 stress 0.2842805 
    ## Run 375 stress 9.040378e-05 
    ## ... Procrustes: rmse 0.0002324128  max resid 0.0005099361 
    ## ... Similar to previous best
    ## Run 376 stress 0.0006936666 
    ## Run 377 stress 0.001327112 
    ## Run 378 stress 9.249184e-05 
    ## ... Procrustes: rmse 0.0002611326  max resid 0.0004922418 
    ## ... Similar to previous best
    ## Run 379 stress 0.001221612 
    ## Run 380 stress 0.2361461 
    ## Run 381 stress 0.001196658 
    ## Run 382 stress 0.00134491 
    ## Run 383 stress 0.2494706 
    ## Run 384 stress 0.1776012 
    ## Run 385 stress 0.001199949 
    ## Run 386 stress 0.001211321 
    ## Run 387 stress 0.2578637 
    ## Run 388 stress 0.1491957 
    ## Run 389 stress 0.001234887 
    ## Run 390 stress 9.903587e-05 
    ## ... Procrustes: rmse 0.01117264  max resid 0.01540663 
    ## Run 391 stress 0.2520602 
    ## Run 392 stress 9.340322e-05 
    ## ... Procrustes: rmse 0.0001423918  max resid 0.0002441436 
    ## ... Similar to previous best
    ## Run 393 stress 0.177602 
    ## Run 394 stress 0.00140058 
    ## Run 395 stress 0.0001474218 
    ## ... Procrustes: rmse 0.008858118  max resid 0.0126866 
    ## Run 396 stress 0.001097254 
    ## Run 397 stress 8.855704e-05 
    ## ... Procrustes: rmse 0.0002835031  max resid 0.0005458324 
    ## ... Similar to previous best
    ## Run 398 stress 0.001329956 
    ## Run 399 stress 0.001304012 
    ## Run 400 stress 9.720561e-05 
    ## ... Procrustes: rmse 0.0001969945  max resid 0.0003612736 
    ## ... Similar to previous best
    ## Run 401 stress 0.001312163 
    ## Run 402 stress 0.0004596436 
    ## ... Procrustes: rmse 0.01567795  max resid 0.02209858 
    ## Run 403 stress 0.001201148 
    ## Run 404 stress 0.2842805 
    ## Run 405 stress 0.001346995 
    ## Run 406 stress 0.2440272 
    ## Run 407 stress 0.1776018 
    ## Run 408 stress 0.001407534 
    ## Run 409 stress 0.2842805 
    ## Run 410 stress 8.386984e-05 
    ## ... Procrustes: rmse 0.0001822965  max resid 0.0003104882 
    ## ... Similar to previous best
    ## Run 411 stress 9.291647e-05 
    ## ... Procrustes: rmse 0.0003030315  max resid 0.000536849 
    ## ... Similar to previous best
    ## Run 412 stress 0.0007093793 
    ## Run 413 stress 0.001339017 
    ## Run 414 stress 0.001188949 
    ## Run 415 stress 0.001323672 
    ## Run 416 stress 0.00107184 
    ## Run 417 stress 0.0004235828 
    ## ... Procrustes: rmse 0.0246302  max resid 0.03399963 
    ## Run 418 stress 0.2696744 
    ## Run 419 stress 0.3083099 
    ## Run 420 stress 9.350209e-05 
    ## ... Procrustes: rmse 0.0002861371  max resid 0.0004843003 
    ## ... Similar to previous best
    ## Run 421 stress 0.2520602 
    ## Run 422 stress 0.1491957 
    ## Run 423 stress 0.1990774 
    ## Run 424 stress 0.0004188924 
    ## ... Procrustes: rmse 0.01496433  max resid 0.02111451 
    ## Run 425 stress 9.983226e-05 
    ## ... Procrustes: rmse 0.0001602842  max resid 0.0002947667 
    ## ... Similar to previous best
    ## Run 426 stress 0.0009498715 
    ## Run 427 stress 0.0007002414 
    ## Run 428 stress 0.001414605 
    ## Run 429 stress 0.001449797 
    ## Run 430 stress 0.001237102 
    ## Run 431 stress 0.001363207 
    ## Run 432 stress 0.001262302 
    ## Run 433 stress 0.1491957 
    ## Run 434 stress 0.2848519 
    ## Run 435 stress 0.1491957 
    ## Run 436 stress 0.0004740762 
    ## ... Procrustes: rmse 0.01592272  max resid 0.02243664 
    ## Run 437 stress 0.1491957 
    ## Run 438 stress 9.884864e-05 
    ## ... Procrustes: rmse 0.003721026  max resid 0.005153467 
    ## ... Similar to previous best
    ## Run 439 stress 8.497055e-05 
    ## ... Procrustes: rmse 0.0002305807  max resid 0.0003956466 
    ## ... Similar to previous best
    ## Run 440 stress 0.1990774 
    ## Run 441 stress 0.0006285275 
    ## Run 442 stress 9.908478e-05 
    ## ... Procrustes: rmse 0.0006426993  max resid 0.001276982 
    ## ... Similar to previous best
    ## Run 443 stress 9.582251e-05 
    ## ... Procrustes: rmse 0.0001559703  max resid 0.0002698089 
    ## ... Similar to previous best
    ## Run 444 stress 0.1776012 
    ## Run 445 stress 0.00117994 
    ## Run 446 stress 0.001385093 
    ## Run 447 stress 0.2834463 
    ## Run 448 stress 0.1491957 
    ## Run 449 stress 0.001178103 
    ## Run 450 stress 8.746469e-05 
    ## ... Procrustes: rmse 0.000275677  max resid 0.0004530623 
    ## ... Similar to previous best
    ## Run 451 stress 0.0004708502 
    ## ... Procrustes: rmse 0.0158695  max resid 0.02236291 
    ## Run 452 stress 0.0005059316 
    ## ... Procrustes: rmse 0.01645085  max resid 0.02316846 
    ## Run 453 stress 7.392607e-05 
    ## ... Procrustes: rmse 0.0005536517  max resid 0.0008026335 
    ## ... Similar to previous best
    ## Run 454 stress 0.0001801606 
    ## ... Procrustes: rmse 0.01597395  max resid 0.02202496 
    ## Run 455 stress 9.061995e-05 
    ## ... Procrustes: rmse 0.0002563923  max resid 0.00047583 
    ## ... Similar to previous best
    ## Run 456 stress 9.872175e-05 
    ## ... Procrustes: rmse 0.0002624402  max resid 0.0004910316 
    ## ... Similar to previous best
    ## Run 457 stress 0.001269224 
    ## Run 458 stress 0.001187738 
    ## Run 459 stress 0.001426827 
    ## Run 460 stress 9.221624e-05 
    ## ... Procrustes: rmse 0.0002045239  max resid 0.0004004959 
    ## ... Similar to previous best
    ## Run 461 stress 0.001242651 
    ## Run 462 stress 9.537418e-05 
    ## ... Procrustes: rmse 0.000259435  max resid 0.0003742061 
    ## ... Similar to previous best
    ## Run 463 stress 0.001418238 
    ## Run 464 stress 0.1491957 
    ## Run 465 stress 9.995235e-05 
    ## ... Procrustes: rmse 0.0002512455  max resid 0.0004400516 
    ## ... Similar to previous best
    ## Run 466 stress 0.000641274 
    ## Run 467 stress 0.1491957 
    ## Run 468 stress 0.001152273 
    ## Run 469 stress 8.211859e-05 
    ## ... Procrustes: rmse 0.0002062119  max resid 0.0003845787 
    ## ... Similar to previous best
    ## Run 470 stress 0.001178002 
    ## Run 471 stress 0.1491957 
    ## Run 472 stress 8.70864e-05 
    ## ... Procrustes: rmse 0.0002422775  max resid 0.0004928318 
    ## ... Similar to previous best
    ## Run 473 stress 0.001228342 
    ## Run 474 stress 0.2440272 
    ## Run 475 stress 0.1491957 
    ## Run 476 stress 0.0004601677 
    ## ... Procrustes: rmse 0.01568768  max resid 0.02211197 
    ## Run 477 stress 0.0004361751 
    ## ... Procrustes: rmse 0.0152659  max resid 0.02153551 
    ## Run 478 stress 0.001213879 
    ## Run 479 stress 0.000215813 
    ## ... Procrustes: rmse 0.01750738  max resid 0.02414459 
    ## Run 480 stress 0.001277405 
    ## Run 481 stress 0.001242181 
    ## Run 482 stress 9.967869e-05 
    ## ... Procrustes: rmse 0.000226539  max resid 0.0005066748 
    ## ... Similar to previous best
    ## Run 483 stress 0.1491957 
    ## Run 484 stress 0.00139282 
    ## Run 485 stress 0.001368088 
    ## Run 486 stress 7.484758e-05 
    ## ... Procrustes: rmse 0.0002460519  max resid 0.0004304222 
    ## ... Similar to previous best
    ## Run 487 stress 0.2440272 
    ## Run 488 stress 0.001268572 
    ## Run 489 stress 0.1491957 
    ## Run 490 stress 0.0007788428 
    ## Run 491 stress 0.0004658937 
    ## ... Procrustes: rmse 0.01577824  max resid 0.02224263 
    ## Run 492 stress 0.0004266015 
    ## ... Procrustes: rmse 0.01509834  max resid 0.02129875 
    ## Run 493 stress 0.001375968 
    ## Run 494 stress 8.758525e-05 
    ## ... Procrustes: rmse 0.0003025357  max resid 0.0005437788 
    ## ... Similar to previous best
    ## Run 495 stress 9.994168e-05 
    ## ... Procrustes: rmse 0.007282132  max resid 0.01051213 
    ## Run 496 stress 0.001241755 
    ## Run 497 stress 9.392671e-05 
    ## ... Procrustes: rmse 0.000233354  max resid 0.0005156851 
    ## ... Similar to previous best
    ## Run 498 stress 9.294501e-05 
    ## ... Procrustes: rmse 0.000286441  max resid 0.0004819631 
    ## ... Similar to previous best
    ## Run 499 stress 0.001209357 
    ## Run 500 stress 0.0002509818 
    ## ... Procrustes: rmse 0.01142694  max resid 0.01622992 
    ## *** Best solution repeated 73 times

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
    ## env[surveyed_sites_env, c(34)]  0.997830 -0.065785 0.7388  0.038 *
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
    ## env[surveyed_sites_env, c(34)]  0.997830 -0.065785 0.7388  0.038 *
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
    ## temperature_median -0.48711 -0.87334 0.0586  0.826  
    ## salinity_median     0.99783 -0.06578 0.7388  0.041 *
    ## oxygen_median       0.79276  0.60953 0.5362  0.822  
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
    ## temperature_median -0.48711 -0.87334 0.0586  1.000
    ## salinity_median     0.99783 -0.06578 0.7388  0.123
    ## oxygen_median       0.79276  0.60953 0.5362  1.000
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
    ## temperature_median -0.56580 -0.82455 0.0778  0.788  
    ## salinity_median     0.99174  0.12828 0.7458  0.035 *
    ## oxygen_median       0.94304  0.33268 0.3551  0.894  
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
    ## temperature_median -0.56580 -0.82455 0.0778  1.000
    ## salinity_median     0.99174  0.12828 0.7458  0.105
    ## oxygen_median       0.94304  0.33268 0.3551  1.000
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
    ## temperature_median -0.10614  0.99435 0.0705  0.716  
    ## salinity_median    -0.78732 -0.61655 0.7842  0.012 *
    ## oxygen_median      -0.96711  0.25434 0.7188  0.051 .
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
    ## salinity_median    -0.78732 -0.61655 0.7842  0.036 *
    ## oxygen_median      -0.96711  0.25434 0.7188  0.153  
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
    ## temperature_median  0.0001064  1.0000000 0.3331  0.088 .
    ## salinity_median    -0.0041863 -0.9999900 0.4855  0.833  
    ## oxygen_median      -0.0027769 -1.0000000 0.7219  0.579  
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
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## temperature_median  0.0001064  1.0000000 0.3331  0.264
    ## salinity_median    -0.0041863 -0.9999900 0.4855  1.000
    ## oxygen_median      -0.0027769 -1.0000000 0.7219  1.000
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
    ##                        NMDS1     NMDS2     r2 Pr(>r)  
    ## temperature_median -0.016725  0.999860 0.0306  0.925  
    ## salinity_median     0.015090 -0.999890 0.6246  0.071 .
    ## oxygen_median       0.054967 -0.998490 0.6093  0.089 .
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
    ##                        NMDS1     NMDS2     r2 Pr(>r)
    ## temperature_median -0.016725  0.999860 0.0306  1.000
    ## salinity_median     0.015090 -0.999890 0.6246  0.213
    ## oxygen_median       0.054967 -0.998490 0.6093  0.267
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
    ## temperature_median       -0.35881 -0.93341 0.1964  0.630
    ## salinity_median          -0.70036  0.71379 0.5991  0.123
    ## oxygen_median             0.98543 -0.17006 0.5185  0.172
    ## distance_to_ocean_mean_m  0.01326 -0.99991 0.1726  0.621
    ## max_depth                -0.49998 -0.86604 0.0383  0.879
    ## logArea                  -0.26513 -0.96421 0.2144  0.557
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
    ## salinity_median          -0.70036  0.71379 0.5991  0.738
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
    ## distance_to_ocean_min_m -0.99980  0.02024 0.6139  0.135
    ## max_depth               -0.17723 -0.98417 0.1185  0.785
    ## logArea                  0.26562 -0.96408 0.2081  0.106
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
    ## distance_to_ocean_min_m -0.99980  0.02024 0.6139  0.405
    ## max_depth               -0.17723 -0.98417 0.1185  1.000
    ## logArea                  0.26562 -0.96408 0.2081  0.318
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
    ## distance_to_ocean_min_m -0.99980  0.02024 0.6139  0.116
    ## max_depth               -0.17723 -0.98417 0.1185  0.783
    ## logArea                  0.26562 -0.96408 0.2081  0.125
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
    ## distance_to_ocean_min_m -0.99980  0.02024 0.6139  0.348
    ## max_depth               -0.17723 -0.98417 0.1185  1.000
    ## logArea                  0.26562 -0.96408 0.2081  0.375
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
    ## distance_to_ocean_min_m -0.94035  0.34022 0.5170  0.164
    ## max_depth               -0.21481 -0.97666 0.1851  0.664
    ## logArea                 -0.06541  0.99786 0.0118  0.959
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
    ## distance_to_ocean_min_m -0.94035  0.34022 0.5170  0.492
    ## max_depth               -0.21481 -0.97666 0.1851  1.000
    ## logArea                 -0.06541  0.99786 0.0118  1.000
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
    ## distance_to_ocean_min_m  0.86175 -0.50733 0.3198  0.088 .
    ## max_depth               -0.68923 -0.72454 0.5689  0.022 *
    ## logArea                 -0.84908  0.52827 0.2993  0.230  
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
    ## distance_to_ocean_min_m  0.86175 -0.50733 0.3198  0.264  
    ## max_depth               -0.68923 -0.72454 0.5689  0.066 .
    ## logArea                 -0.84908  0.52827 0.2993  0.690  
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
    ## distance_to_ocean_min_m  0.86963  0.49371 0.6793  0.195
    ## max_depth                0.68534  0.72822 0.0510  0.962
    ## logArea                 -0.52042  0.85391 0.2187  0.537
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
    ## distance_to_ocean_min_m  0.86963  0.49371 0.6793  0.585
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
    ##                               NMDS1       NMDS2     r2 Pr(>r)   
    ## distance_to_ocean_min_m -0.00070192  1.00000000 0.6553  0.056 . 
    ## max_depth                0.00171321 -1.00000000 0.8337  0.020 * 
    ## logArea                  0.00143052  1.00000000 0.8241  0.005 **
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
    ## distance_to_ocean_min_m -0.00070192  1.00000000 0.6553  0.168  
    ## max_depth                0.00171321 -1.00000000 0.8337  0.060 .
    ## logArea                  0.00143052  1.00000000 0.8241  0.015 *
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
    ##       Significance: 0.244 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.297 0.332 0.359 0.382 
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
    ##       Significance: 0.005 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.564 0.607 0.634 0.654 
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
    ##       Significance: 0.47 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.587 0.620 0.651 0.683 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_mant_pv <- rbind(SD_beta_env_mant_t$signif, SD_beta_env_mant_s$signif, SD_beta_env_mant_o$signif)
SD_beta_env_mant_pv <- SD_beta_env_mant_pv[,1]
(SD_beta_env_mant_pv <- p.adjust(SD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.732 0.015 1.000

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
    ##       Significance: 0.358 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.341 0.374 0.403 0.428 
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
    ##       Significance: 0.008 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.573 0.621 0.643 0.662 
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
    ##       Significance: 0.6 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.475 0.521 0.568 0.598 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_MS_mant_pv <- rbind(SD_beta_env_MS_mant_t$signif, SD_beta_env_MS_mant_s$signif, SD_beta_env_MS_mant_o$signif)
SD_beta_env_MS_mant_pv <- SD_beta_env_MS_mant_pv[,1]
(SD_beta_env_MS_mant_pv <- p.adjust(SD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.024 1.000

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
    ##       Significance: 0.77 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.126 0.177 0.282 0.398 
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
    ##       Significance: 0.036 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.329 0.401 0.464 0.550 
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
    ##       Significance: 0.025 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.330 0.445 0.518 0.608 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_OM_mant_pv <- rbind(SD_beta_env_OM_mant_t$signif, SD_beta_env_OM_mant_s$signif, SD_beta_env_OM_mant_o$signif)
SD_beta_env_OM_mant_pv <- SD_beta_env_OM_mant_pv[,1]
(SD_beta_env_OM_mant_pv <- p.adjust(SD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.108 0.075

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
    ##       Significance: 0.223 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0353 0.0576 0.0740 0.0856 
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
    ##       Significance: 0.018 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.468 0.513 0.540 0.568 
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
    ##       Significance: 0.42 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.713 0.735 0.746 0.765 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_SO_mant_pv <- rbind(SD_beta_env_SO_mant_t$signif, SD_beta_env_SO_mant_s$signif, SD_beta_env_SO_mant_o$signif)
SD_beta_env_SO_mant_pv <- SD_beta_env_SO_mant_pv[,1]
(SD_beta_env_SO_mant_pv <- p.adjust(SD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 0.669 0.054 1.000

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
    ##       Significance: 0.728 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.235 0.345 0.421 0.662 
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
    ##       Significance: 0.143 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.240 0.333 0.388 0.425 
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
    ##       Significance: 0.097 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.257 0.405 0.496 0.628 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_M_mant_pv <- rbind(SD_beta_env_M_mant_t$signif, SD_beta_env_M_mant_s$signif, SD_beta_env_M_mant_o$signif)
SD_beta_env_M_mant_pv <- SD_beta_env_M_mant_pv[,1]
(SD_beta_env_M_mant_pv <- p.adjust(SD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.429 0.291

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
    ## 0.222 0.288 0.332 0.396 
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
    ##       Significance: 0.45 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.554 0.579 0.595 0.626 
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
    ##       Significance: 0.65 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.204 0.230 0.254 0.288 
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
    ##       Significance: 0.427 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0824 0.0975 0.1066 0.1245 
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
    ##       Significance: 0.411 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.450 0.494 0.531 0.559 
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
    ##       Significance: 0.747 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.218 0.265 0.295 0.332 
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
    ##       Significance: 0.75 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0644 0.0859 0.1061 0.1380 
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
    ##       Significance: 0.62 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.235 0.268 0.286 0.353 
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
    ##       Significance: 0.05 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.176 0.234 0.298 0.372 
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
    ##       Significance: 0.118 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.211 0.259 0.289 0.326 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_OM_mant_pv <- rbind( SD_beta_geo_OM_mant_dmean$signif, SD_beta_geo_OM_mant_md$signif, SD_beta_geo_OM_mant_la$signif)
SD_beta_geo_OM_mant_pv <- SD_beta_geo_OM_mant_pv[,1]
(SD_beta_geo_OM_mant_pv <- p.adjust(SD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.150 0.354

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
    ##       Significance: 0.462 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.670 0.700 0.721 0.740 
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
    ##       Significance: 0.597 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.134 0.167 0.188 0.219 
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
    ##       Significance: 0.221 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.280 0.294 0.310 0.323 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_SO_mant_pv <- rbind( SD_beta_geo_SO_mant_dmean$signif, SD_beta_geo_SO_mant_md$signif, SD_beta_geo_SO_mant_la$signif)
SD_beta_geo_SO_mant_pv <- SD_beta_geo_SO_mant_pv[,1]
(SD_beta_geo_SO_mant_pv <- p.adjust(SD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 1.000 0.663

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
    ##       Significance: 0.502 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.250 0.409 0.585 0.787 
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
    ##       Significance: 0.004 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.274 0.423 0.539 0.582 
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
    ##       Significance: 0.033 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.225 0.306 0.458 0.539 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_M_mant_pv <- rbind( SD_beta_geo_M_mant_dmean$signif, SD_beta_geo_M_mant_md$signif, SD_beta_geo_M_mant_la$signif)
SD_beta_geo_M_mant_pv <- SD_beta_geo_M_mant_pv[,1]
(SD_beta_geo_M_mant_pv <- p.adjust(SD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.012 0.099

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
        panel.border = element_rect(fill = NA, size = 14), 
        axis.text = element_text(color = "black"), 
        legend.key = element_blank()) +
  annotate("text", x = -0.7, y = 0.8, size = 5,
           label = paste("Stress: ", round(SD_beta_ref_NMDS$stress, digits = 2))) +
  labs(colour = "Site type:  ", fill = "Site type:  ", tag = "a") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
```

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
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
