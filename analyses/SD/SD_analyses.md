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
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: row_sum
    ##                    Sum Sq Df F value  Pr(>F)  
    ## (Intercept)        1909.3  1  7.2025 0.01781 *
    ## temperature_median  432.9  1  1.6329 0.22208  
    ## salinity_median       2.8  1  0.0104 0.92030  
    ## oxygen_median       340.9  1  1.2859 0.27584  
    ## pH_median          1789.6  1  6.7509 0.02105 *
    ## Residuals          3711.3 14                  
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
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: row_sum
    ##                    Sum Sq Df F value   Pr(>F)   
    ## (Intercept)         289.1  1  0.7882 0.388657   
    ## temperature_median   13.2  1  0.0361 0.851873   
    ## salinity_median    2961.0  1  8.0741 0.012377 * 
    ## oxygen_median      3808.8  1 10.3860 0.005693 **
    ## Residuals          5500.9 15                    
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


mod <- aov(temperature_median ~ Site_type, SR_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2   2.94   1.470   1.092  0.359
    ## Residuals   16  21.54   1.346               
    ## 4 observations deleted due to missingness

``` r
mod <- aov(row_sum ~ oxygen_median + Site_type, SR_env)
car::vif(mod)
```

    ##                   GVIF Df GVIF^(1/(2*Df))
    ## oxygen_median 3.374424  1        1.836961
    ## Site_type     3.374424  2        1.355345

``` r
mod <- aov(salinity_median ~ Site_type, SR_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Site_type    2 147.34   73.67   13.61 0.000353 ***
    ## Residuals   16  86.62    5.41                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 4 observations deleted due to missingness

``` r
mod <- aov(row_sum ~ salinity_median + Site_type, SR_env)
car::vif(mod)
```

    ##                     GVIF Df GVIF^(1/(2*Df))
    ## salinity_median 2.701125  1        1.643510
    ## Site_type       2.701125  2        1.281995

``` r
mod <- aov(oxygen_median ~ Site_type, SR_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Site_type    2 14.440    7.22      19 5.95e-05 ***
    ## Residuals   16  6.081    0.38                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 4 observations deleted due to missingness

``` r
mod <- aov(row_sum ~ oxygen_median + Site_type, SR_env)
car::vif(mod)
```

    ##                   GVIF Df GVIF^(1/(2*Df))
    ## oxygen_median 3.374424  1        1.836961
    ## Site_type     3.374424  2        1.355345

``` r
SR_ss_env <- SR_env[surveyed_sites_env,]
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
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: row_sum
    ##                            Sum Sq Df F value  Pr(>F)  
    ## (Intercept)                 440.2  1  0.9946 0.34467  
    ## volume_m3_w_chemocline      537.0  1  1.2135 0.29922  
    ## volume_m3                   375.3  1  0.8479 0.38115  
    ## surface_area_m2            1629.2  1  3.6815 0.08723 .
    ## distance_to_ocean_min_m       8.8  1  0.0199 0.89103  
    ## distance_to_ocean_mean_m      9.6  1  0.0216 0.88633  
    ## distance_to_ocean_median_m   30.6  1  0.0692 0.79837  
    ## tidal_lag_time_minutes       42.8  1  0.0968 0.76279  
    ## tidal_efficiency              5.7  1  0.0128 0.91243  
    ## perimeter_fromSat           332.3  1  0.7508 0.40874  
    ## max_depth                    29.7  1  0.0670 0.80152  
    ## logArea                    1160.0  1  2.6212 0.13990  
    ## Residuals                  3982.9  9                  
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
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: row_sum
    ##                          Sum Sq Df F value  Pr(>F)  
    ## (Intercept)               115.8  1  0.1720 0.68323  
    ## distance_to_ocean_min_m  4916.0  1  7.3006 0.01459 *
    ## max_depth                  18.2  1  0.0271 0.87116  
    ## logArea                   699.5  1  1.0388 0.32161  
    ## Residuals               12120.7 18                  
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

mod <- aov(distance_to_ocean_min_m ~ Site_type, SR_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Site_type    2  80935   40468   16.49 7.05e-05 ***
    ## Residuals   19  46637    2455                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 1 observation deleted due to missingness

``` r
mod <- aov(row_sum ~ distance_to_ocean_min_m + Site_type, SR_env)
car::vif(mod)
```

    ##                             GVIF Df GVIF^(1/(2*Df))
    ## distance_to_ocean_min_m 2.735421  1        1.653911
    ## Site_type               2.735421  2        1.286045

``` r
mod <- aov(max_depth ~ Site_type, SR_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2    535   267.5   2.235  0.134
    ## Residuals   19   2274   119.7               
    ## 1 observation deleted due to missingness

``` r
mod <- aov(row_sum ~ max_depth + Site_type, SR_env)
car::vif(mod)
```

    ##               GVIF Df GVIF^(1/(2*Df))
    ## max_depth 1.235315  1        1.111447
    ## Site_type 1.235315  2        1.054252

``` r
mod <- aov(logArea ~ Site_type, SR_env)
summary(mod)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## Site_type    2  14.37   7.184   2.341  0.123
    ## Residuals   19  58.32   3.069               
    ## 1 observation deleted due to missingness

``` r
mod <- aov(row_sum ~ logArea + Site_type, SR_env)
car::vif(mod)
```

    ##               GVIF Df GVIF^(1/(2*Df))
    ## logArea   1.246371  1        1.116410
    ## Site_type 1.246371  2        1.056603

``` r
SR_ss_geo <- SR_env[surveyed_sites,]
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
mod <-  lm(row_sum ~ salinity_median + oxygen_median + temperature_median + distance_to_ocean_min_m + max_depth + logArea, SR_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ salinity_median + oxygen_median + temperature_median + 
    ##     distance_to_ocean_min_m + max_depth + logArea, data = SR_env)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -31.899  -6.041   0.311   6.860  27.594 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             -84.35957  153.48785  -0.550   0.5927  
    ## salinity_median           3.86324    2.01538   1.917   0.0794 .
    ## oxygen_median            10.28180    4.31263   2.384   0.0345 *
    ## temperature_median       -4.50254    4.13686  -1.088   0.2978  
    ## distance_to_ocean_min_m  -0.01767    0.09508  -0.186   0.8557  
    ## max_depth                -0.72903    0.55986  -1.302   0.2173  
    ## logArea                  11.20688    3.83364   2.923   0.0128 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 15.9 on 12 degrees of freedom
    ##   (4 observations deleted due to missingness)
    ## Multiple R-squared:  0.8362, Adjusted R-squared:  0.7543 
    ## F-statistic: 10.21 on 6 and 12 DF,  p-value: 0.0003995

``` r
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: row_sum
    ##                          Sum Sq Df F value  Pr(>F)  
    ## (Intercept)               76.37  1  0.3021 0.59266  
    ## salinity_median          928.96  1  3.6744 0.07937 .
    ## oxygen_median           1437.02  1  5.6840 0.03450 *
    ## temperature_median       299.49  1  1.1846 0.29780  
    ## distance_to_ocean_min_m    8.73  1  0.0345 0.85565  
    ## max_depth                428.68  1  1.6956 0.21730  
    ## logArea                 2160.51  1  8.5457 0.01276 *
    ## Residuals               3033.82 12                  
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
mod <-  lm(row_sum ~ oxygen_median + temperature_median + distance_to_ocean_min_m + logArea, SR_env)
summary(mod)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ oxygen_median + temperature_median + distance_to_ocean_min_m + 
    ##     logArea, data = SR_env)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -23.557 -13.492  -1.246   9.123  37.830 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)             152.15529  123.23638   1.235   0.2373  
    ## oxygen_median            12.21960    4.59804   2.658   0.0187 *
    ## temperature_median       -7.57170    4.25822  -1.778   0.0971 .
    ## distance_to_ocean_min_m  -0.16235    0.06355  -2.555   0.0229 *
    ## logArea                   8.10894    2.96536   2.735   0.0161 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 17.55 on 14 degrees of freedom
    ##   (4 observations deleted due to missingness)
    ## Multiple R-squared:  0.7672, Adjusted R-squared:  0.7007 
    ## F-statistic: 11.54 on 4 and 14 DF,  p-value: 0.0002359

``` r
Anova(mod, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: row_sum
    ##                         Sum Sq Df F value  Pr(>F)  
    ## (Intercept)              469.5  1  1.5244 0.23728  
    ## oxygen_median           2175.3  1  7.0627 0.01875 *
    ## temperature_median       973.8  1  3.1618 0.09710 .
    ## distance_to_ocean_min_m 2010.1  1  6.5264 0.02291 *
    ## logArea                 2303.1  1  7.4778 0.01613 *
    ## Residuals               4311.9 14                  
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
par(mfrow=c(2,2)) 
#### Environmental
### SRic ANOVA
## Surveyed sites
# Temperature
# Interaction
SD_SRic_am_temp <- aov(row_sum ~ temperature_median * Site_type, data = SR_env[surveyed_sites_env,])
summary(SD_SRic_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## temperature_median            1   1618    1618   4.099 0.063963 .  
    ## Site_type                     2   9895    4947  12.530 0.000928 ***
    ## temperature_median:Site_type  2   1878     939   2.378 0.131806    
    ## Residuals                    13   5133     395                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_lm_temp <- lm(row_sum ~ temperature_median + Site_type, data = SR_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_lm_temp)
    ## W = 0.9529, p-value = 0.4421

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_lm_temp) ~ SR_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  3.3867 0.05937 .
    ##       16                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.313815          0.0051196     0.097272

``` r
# Plot residuals
plot(SD_SRic_lm_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-1.png)<!-- -->

``` r
# Test relationship
SD_SRic_am_temp <- aov(row_sum ~ temperature_median + Site_type, data = SR_env[surveyed_sites_env,])

# Salinity
# Interaction
SD_SRic_am_sal <- aov(row_sum ~ salinity_median * Site_type, data = SR_env[surveyed_sites_env,])
summary(SD_SRic_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median            1   9173    9173  22.544 0.000381 ***
    ## Site_type                  2   2620    1310   3.220 0.073142 .  
    ## salinity_median:Site_type  2   1442     721   1.772 0.208663    
    ## Residuals                 13   5289     407                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_lm_sal <- lm(row_sum ~ salinity_median + Site_type, data = SR_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_lm_sal)
    ## W = 0.93527, p-value = 0.2164

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_lm_sal) ~ SR_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  4.3382 0.03124 *
    ##       16                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.436669          0.0040095     0.076181

``` r
# Plot residuals
plot(SD_SRic_lm_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-2.png)<!-- -->

``` r
# Test relationship
SD_SRic_am_sal <- aov(row_sum ~ salinity_median + Site_type, data = SR_env[surveyed_sites_env,])

# Oxygen
# Interaction
SD_SRic_am_oxy <- aov(row_sum ~ oxygen_median * Site_type, data = SR_env[surveyed_sites_env,])
summary(SD_SRic_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median            1   9605    9605  29.593 0.000113 ***
    ## Site_type                2   2160    1080   3.328 0.068075 .  
    ## oxygen_median:Site_type  2   2540    1270   3.913 0.046733 *  
    ## Residuals               13   4219     325                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_lm_oxy <- lm(row_sum ~ oxygen_median + Site_type, data = SR_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_lm_oxy)
    ## W = 0.94821, p-value = 0.3682

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_lm_oxy) ~ SR_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  2.2869 0.1338
    ##       16

``` r
# Check for outliers
outlierTest(SD_SRic_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.254152          0.0057649      0.10953

``` r
# Plot residuals
plot(SD_SRic_lm_oxy)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-3.png)<!-- -->

``` r
# Test relationship
SD_SRic_am_oxy <- aov(row_sum ~ oxygen_median + Site_type, data = SR_env[surveyed_sites_env,])

# Anova outputs
summary(SD_SRic_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value  Pr(>F)   
    ## temperature_median  1   1618    1618   3.463 0.08249 . 
    ## Site_type           2   9895    4947  10.585 0.00136 **
    ## Residuals          15   7011     467                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median  1   9173    9173  20.440 0.000406 ***
    ## Site_type        2   2620    1310   2.919 0.084953 .  
    ## Residuals       15   6731     449                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median  1   9605    9605  21.314 0.000336 ***
    ## Site_type      2   2160    1080   2.397 0.124957    
    ## Residuals     15   6759     451                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(SD_SRic_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(SD_SRic_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(SD_SRic_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.494932100 0.008150692          NA 0.002435311 0.509718354          NA
    ## [7] 0.002014076 0.749743578          NA

``` r
## Mixed and stratified lakes
# Temperature
# Interaction
SD_SRic_MS_am_temp <- aov(row_sum ~ temperature_median * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_SRic_MS_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value  Pr(>F)   
    ## temperature_median            1   2185    2185   5.166 0.04221 * 
    ## Site_type                     1   6632    6632  15.684 0.00189 **
    ## temperature_median:Site_type  1   1000    1000   2.366 0.14996   
    ## Residuals                    12   5074     423                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_MS_lm_temp <- lm(row_sum ~ temperature_median + Site_type, data = SR_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_MS_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_MS_lm_temp)
    ## W = 0.92755, p-value = 0.2229

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_MS_lm_temp) ~ SR_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)  
    ## group  1  7.2669 0.0174 *
    ##       14                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_MS_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.305486          0.0062766      0.10043

``` r
# Plot residuals
plot(SD_SRic_MS_lm_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-4.png)<!-- -->

``` r
# Test relationship
SD_SRic_MS_am_temp <- aov(row_sum ~ temperature_median + Site_type, data = SR_env[mixed_stratified_lakes,])

# Salinity
# Interaction
SD_SRic_MS_am_sal <- aov(row_sum ~ salinity_median * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_SRic_MS_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value  Pr(>F)   
    ## salinity_median            1   6790    6790  17.997 0.00114 **
    ## Site_type                  1   2141    2141   5.674 0.03463 * 
    ## salinity_median:Site_type  1   1434    1434   3.800 0.07502 . 
    ## Residuals                 12   4527     377                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_MS_lm_sal <- lm(row_sum ~ salinity_median + Site_type, data = SR_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_MS_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_MS_lm_sal)
    ## W = 0.8911, p-value = 0.05799

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_MS_lm_sal) ~ SR_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  1  10.052 0.006809 **
    ##       14                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_MS_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.582117          0.0037676     0.060281

``` r
# Plot residuals
plot(SD_SRic_MS_lm_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-5.png)<!-- -->

``` r
# Test relationship
SD_SRic_MS_am_sal <- aov(row_sum ~ salinity_median + Site_type, data = SR_env[mixed_stratified_lakes,])

# Oxygen
# Interaction
SD_SRic_MS_am_oxy <- aov(row_sum ~ oxygen_median * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_SRic_MS_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value  Pr(>F)   
    ## oxygen_median            1   6282    6282  17.871 0.00117 **
    ## Site_type                1   2469    2469   7.024 0.02117 * 
    ## oxygen_median:Site_type  1   1923    1923   5.471 0.03745 * 
    ## Residuals               12   4218     351                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_MS_lm_oxy <- lm(row_sum ~ oxygen_median + Site_type, data = SR_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_MS_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_MS_lm_oxy)
    ## W = 0.93808, p-value = 0.3262

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_MS_lm_oxy) ~ SR_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)  
    ## group  1  5.8654 0.0296 *
    ##       14                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_MS_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.365788          0.0056138      0.08982

``` r
# Plot residuals
plot(SD_SRic_MS_lm_oxy)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-6.png)<!-- -->

``` r
# Test relationship
SD_SRic_MS_am_oxy <- aov(row_sum ~ oxygen_median + Site_type, data = SR_env[mixed_stratified_lakes,])

# Anova outputs
summary(SD_SRic_MS_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value  Pr(>F)   
    ## temperature_median  1   2185    2185   4.675 0.04983 * 
    ## Site_type           1   6632    6632  14.193 0.00235 **
    ## Residuals          13   6075     467                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_MS_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value  Pr(>F)   
    ## salinity_median  1   6790    6790  14.808 0.00201 **
    ## Site_type        1   2141    2141   4.669 0.04997 * 
    ## Residuals       13   5961     459                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_MS_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value  Pr(>F)   
    ## oxygen_median  1   6282    6282  13.298 0.00296 **
    ## Site_type      1   2469    2469   5.226 0.03967 * 
    ## Residuals     13   6141     472                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(SD_SRic_MS_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(SD_SRic_MS_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(SD_SRic_MS_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.29900871 0.01409312         NA 0.01208932 0.29979203         NA 0.01773844
    ## [8] 0.23799529         NA

``` r
## Ocean sites and mixed lakes
# Temperature
# Interaction
SD_SRic_OM_am_temp <- aov(row_sum ~ temperature_median * Site_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_SRic_OM_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median            1     60    60.4   0.087  0.777
    ## Site_type                     1    318   318.1   0.458  0.520
    ## temperature_median:Site_type  1   1784  1784.0   2.569  0.153
    ## Residuals                     7   4862   694.5

``` r
# Linear model
SD_SRic_OM_lm_temp <- lm(row_sum ~ temperature_median + Site_type, data = SR_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_OM_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_OM_lm_temp)
    ## W = 0.96592, p-value = 0.8426

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_OM_lm_temp) ~ SR_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.3487 0.5694
    ##        9

``` r
# Check for outliers
outlierTest(SD_SRic_OM_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN  2.63858           0.033492      0.36842

``` r
# Plot residuals
plot(SD_SRic_OM_lm_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-7.png)<!-- -->

``` r
# Test relationship
SD_SRic_OM_am_temp <- aov(row_sum ~ temperature_median + Site_type, data = SR_env[ocean_mixed_sites_env,])

# Salinity
# Interaction
SD_SRic_OM_am_sal <- aov(row_sum ~ salinity_median * Site_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_SRic_OM_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)
    ## salinity_median            1   1694  1693.9   2.283  0.175
    ## Site_type                  1    124   124.5   0.168  0.694
    ## salinity_median:Site_type  1     11    11.3   0.015  0.905
    ## Residuals                  7   5194   742.1

``` r
# Linear model
SD_SRic_OM_lm_sal <- lm(row_sum ~ salinity_median + Site_type, data = SR_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_OM_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_OM_lm_sal)
    ## W = 0.91894, p-value = 0.3099

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_OM_lm_sal) ~ SR_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.2029  0.663
    ##        9

``` r
# Check for outliers
outlierTest(SD_SRic_OM_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.736003           0.029085      0.31994

``` r
# Plot residuals
plot(SD_SRic_OM_lm_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-8.png)<!-- -->

``` r
# Test relationship
SD_SRic_OM_am_sal <- aov(row_sum ~ salinity_median + Site_type, data = SR_env[ocean_mixed_sites_env,])

# Oxygen
# Interaction
SD_SRic_OM_am_oxy <- aov(row_sum ~ oxygen_median * Site_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_SRic_OM_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)  
    ## oxygen_median            1   2353  2353.3   4.022 0.0849 .
    ## Site_type                1    472   471.9   0.806 0.3990  
    ## oxygen_median:Site_type  1    103   102.7   0.176 0.6878  
    ## Residuals                7   4096   585.2                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_OM_lm_oxy <- lm(row_sum ~ oxygen_median + Site_type, data = SR_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_OM_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_OM_lm_oxy)
    ## W = 0.97017, p-value = 0.8883

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_OM_lm_oxy) ~ SR_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  3.0772 0.1133
    ##        9

``` r
# Check for outliers
outlierTest(SD_SRic_OM_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.567212           0.037159      0.40875

``` r
# Plot residuals
plot(SD_SRic_OM_lm_oxy)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-9.png)<!-- -->

``` r
# Test relationship
SD_SRic_OM_am_oxy <- aov(row_sum ~ oxygen_median + Site_type, data = SR_env[ocean_mixed_sites_env,])

# Anova outputs
summary(SD_SRic_OM_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1     60    60.4   0.073  0.794
    ## Site_type           1    318   318.1   0.383  0.553
    ## Residuals           8   6646   830.7

``` r
summary(SD_SRic_OM_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value Pr(>F)
    ## salinity_median  1   1694  1693.9   2.603  0.145
    ## Site_type        1    124   124.5   0.191  0.673
    ## Residuals        8   5206   650.7

``` r
summary(SD_SRic_OM_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)  
    ## oxygen_median  1   2353  2353.3   4.484 0.0671 .
    ## Site_type      1    472   471.9   0.899 0.3708  
    ## Residuals      8   4199   524.9                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(SD_SRic_OM_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(SD_SRic_OM_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(SD_SRic_OM_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 1.0000000        NA 0.8718936 1.0000000        NA 0.4025299
    ## [8] 1.0000000        NA

``` r
## Stratified lakes and ocean sites
# Temperature
# Interaction
SD_SRic_SO_am_temp <- aov(row_sum ~ temperature_median * Site_type, data = SR_env[ocean_stratified_sites_env,])
summary(SD_SRic_SO_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## temperature_median            1     84      84   1.781   0.2238    
    ## Site_type                     1   6992    6992 148.334 5.76e-06 ***
    ## temperature_median:Site_type  1    690     690  14.630   0.0065 ** 
    ## Residuals                     7    330      47                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_SO_lm_temp <- lm(row_sum ~ temperature_median + Site_type, data = SR_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_SO_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_SO_lm_temp)
    ## W = 0.88507, p-value = 0.1206

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_SO_lm_temp) ~ SR_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  1.0308 0.3365
    ##        9

``` r
# Check for outliers
outlierTest(SD_SRic_SO_lm_temp)
```

    ##     rstudent unadjusted p-value Bonferroni p
    ## IBK 4.376517          0.0032482      0.03573

``` r
# Plot residuals
plot(SD_SRic_SO_lm_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-10.png)<!-- -->

``` r
# Test relationship
SD_SRic_SO_am_temp <- aov(row_sum ~ temperature_median + Site_type, data = SR_env[ocean_stratified_sites_env,])

# Salinity
# Interaction
SD_SRic_SO_am_sal <- aov(row_sum ~ salinity_median * Site_type, data = SR_env[ocean_stratified_sites_env,])
summary(SD_SRic_SO_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median            1   4627    4627  37.805 0.000468 ***
    ## Site_type                  1   2603    2603  21.268 0.002450 ** 
    ## salinity_median:Site_type  1      9       9   0.074 0.793928    
    ## Residuals                  7    857     122                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_SO_lm_sal <- lm(row_sum ~ salinity_median + Site_type, data = SR_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_SO_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_SO_lm_sal)
    ## W = 0.89263, p-value = 0.15

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_SO_lm_sal) ~ SR_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  2.1706 0.1748
    ##        9

``` r
# Check for outliers
outlierTest(SD_SRic_SO_lm_sal)
```

    ##     rstudent unadjusted p-value Bonferroni p
    ## IBK 7.420502         0.00014684    0.0016153

``` r
# Plot residuals
plot(SD_SRic_SO_lm_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-11.png)<!-- -->

``` r
# Test relationship
SD_SRic_SO_am_sal <- aov(row_sum ~ salinity_median + Site_type, data = SR_env[ocean_stratified_sites_env,])

# Oxygen
# Interaction
SD_SRic_SO_am_oxy <- aov(row_sum ~ oxygen_median * Site_type, data = SR_env[ocean_stratified_sites_env,])
summary(SD_SRic_SO_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median            1   4589    4589  258.57 8.74e-07 ***
    ## Site_type                1   2491    2491  140.35 6.93e-06 ***
    ## oxygen_median:Site_type  1    892     892   50.26 0.000196 ***
    ## Residuals                7    124      18                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_SO_lm_oxy <- lm(row_sum ~ oxygen_median + Site_type, data = SR_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_SO_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_SO_lm_oxy)
    ## W = 0.90127, p-value = 0.1916

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_SO_lm_oxy) ~ SR_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1   1.841 0.2079
    ##        9

``` r
# Check for outliers
outlierTest(SD_SRic_SO_lm_oxy)
```

    ##     rstudent unadjusted p-value Bonferroni p
    ## IBK  7.03862         0.00020443    0.0022488

``` r
# Plot residuals
plot(SD_SRic_SO_lm_oxy)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-12.png)<!-- -->

``` r
# Test relationship
SD_SRic_SO_am_oxy <- aov(row_sum ~ oxygen_median + Site_type, data = SR_env[ocean_stratified_sites_env,])

# Anova outputs
summary(SD_SRic_SO_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value   Pr(>F)    
    ## temperature_median  1     84      84   0.659     0.44    
    ## Site_type           1   6992    6992  54.863 7.57e-05 ***
    ## Residuals           8   1020     127                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_SO_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median  1   4627    4627   42.76 0.000181 ***
    ## Site_type        1   2603    2603   24.05 0.001187 ** 
    ## Residuals        8    866     108                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_SO_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value  Pr(>F)    
    ## oxygen_median  1   4589    4589   36.13 0.00032 ***
    ## Site_type      1   2491    2491   19.61 0.00220 ** 
    ## Residuals      8   1016     127                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(SD_SRic_SO_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(SD_SRic_SO_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(SD_SRic_SO_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000000 0.0004541492           NA 0.0010834964 0.0071226989
    ## [6]           NA 0.0019174822 0.0132121587           NA

``` r
## Ocean sites
SD_SRic_env_O_lm <- lm(row_sum ~ salinity_median + oxygen_median + temperature_median, SR_env[ocean_sites_env,])
summary(SD_SRic_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ salinity_median + oxygen_median + temperature_median, 
    ##     data = SR_env[ocean_sites_env, ])
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
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 2 and 0 DF,  p-value: NA

``` r
p_values <- summary(SD_SRic_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##     (Intercept) salinity_median   oxygen_median 
    ##             NaN             NaN             NaN

``` r
## Mixed lakes
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
Anova(SD_SRic_env_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: row_sum
    ##                     Sum Sq Df F value Pr(>F)
    ## (Intercept)           0.48  1  0.0006 0.9812
    ## salinity_median     137.16  1  0.1794 0.6937
    ## oxygen_median       975.73  1  1.2762 0.3218
    ## temperature_median  505.41  1  0.6610 0.4618
    ## Residuals          3058.30  4

``` r
p_values <- summary(SD_SRic_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
## Stratified lakes
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
Anova(SD_SRic_env_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: row_sum
    ##                    Sum Sq Df F value Pr(>F)
    ## (Intercept)         1.314  1  0.0622 0.8154
    ## salinity_median    36.742  1  1.7380 0.2578
    ## oxygen_median      10.197  1  0.4824 0.5256
    ## temperature_median  0.008  1  0.0004 0.9852
    ## Residuals          84.560  4

``` r
p_values <- summary(SD_SRic_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
#### Geographical
### SRic ANOVA
## Surveyed sites
# Distance
# Interaction
SD_SRic_am_dist <- aov(row_sum ~ distance_to_ocean_min_m * Site_type, data = SR_env[surveyed_sites,])
summary(SD_SRic_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m            1   6762    6762  11.938 0.00326 **
    ## Site_type                          2   4154    2077   3.667 0.04889 * 
    ## distance_to_ocean_min_m:Site_type  2     32      16   0.028 0.97262   
    ## Residuals                         16   9063     566                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_lm_dist <- lm(row_sum ~ distance_to_ocean_min_m + Site_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_lm_dist)
    ## W = 0.94125, p-value = 0.2101

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_lm_dist) ~ SR_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  4.6736 0.02235 *
    ##       19                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.170152          0.0055954       0.1231

``` r
# Plot residuals
plot(SD_SRic_lm_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-13.png)<!-- -->

``` r
# Test relationship
SD_SRic_am_dist <- aov(row_sum ~ distance_to_ocean_min_m + Site_type, data = SR_env[surveyed_sites,])

# Max depth
# Interaction
SD_SRic_am_mxd <- aov(row_sum ~ max_depth * Site_type, data = SR_env[surveyed_sites,])
summary(SD_SRic_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth            1     33      33   0.093    0.764    
    ## Site_type            2  12617    6308  17.595 9.11e-05 ***
    ## max_depth:Site_type  2   1624     812   2.264    0.136    
    ## Residuals           16   5737     359                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_lm_mxd <- lm(row_sum ~ max_depth + Site_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_lm_mxd)
    ## W = 0.95124, p-value = 0.3342

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_lm_mxd) ~ SR_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  1.0674 0.3636
    ##       19

``` r
# Check for outliers
outlierTest(SD_SRic_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.218891          0.0050385      0.11085

``` r
# Plot residuals
plot(SD_SRic_lm_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-14.png)<!-- -->

``` r
# Test relationship
SD_SRic_am_mxd <- aov(row_sum ~ max_depth + Site_type, data = SR_env[surveyed_sites,])

# Log area
# Interaction
SD_SRic_am_lga <- aov(row_sum ~ logArea * Site_type, data = SR_env[surveyed_sites,])
summary(SD_SRic_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## logArea            1   2166    2166   7.338   0.0155 *  
    ## Site_type          2  10588    5294  17.931 8.21e-05 ***
    ## logArea:Site_type  2   2532    1266   4.287   0.0323 *  
    ## Residuals         16   4724     295                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_lm_lga <- lm(row_sum ~ logArea + Site_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_lm_lga)
    ## W = 0.97073, p-value = 0.7279

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_lm_lga) ~ SR_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  1.0393  0.373
    ##       19

``` r
# Check for outliers
outlierTest(SD_SRic_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OOO -2.994945           0.008142      0.17912

``` r
# Plot residuals
plot(SD_SRic_lm_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-15.png)<!-- -->

``` r
# Test relationship
SD_SRic_am_lga <- aov(row_sum ~ logArea + Site_type, data = SR_env[surveyed_sites,])

# Anova outputs
summary(SD_SRic_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)   
    ## distance_to_ocean_min_m  1   6762    6762  13.383 0.0018 **
    ## Site_type                2   4154    2077   4.111 0.0339 * 
    ## Residuals               18   9095     505                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth    1     33      33   0.082 0.778361    
    ## Site_type    2  12617    6308  15.428 0.000125 ***
    ## Residuals   18   7360     409                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_am_lga)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## logArea      1   2166    2166   5.375 0.032401 *  
    ## Site_type    2  10588    5294  13.134 0.000304 ***
    ## Residuals   18   7256     403                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
dist_p_values <- summary(SD_SRic_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(SD_SRic_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(SD_SRic_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.0107882928 0.2031108319           NA 1.0000000000 0.0007505244
    ## [6]           NA 0.1944063709 0.0018228571           NA

``` r
## Mixed and stratified lakes
# Distance
# Interaction
SD_SRic_MS_am_dist <- aov(row_sum ~ distance_to_ocean_min_m * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_SRic_MS_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)  
    ## distance_to_ocean_min_m            1   4715    4715   9.296 0.0101 *
    ## Site_type                          1   4066    4066   8.017 0.0151 *
    ## distance_to_ocean_min_m:Site_type  1     25      25   0.049 0.8284  
    ## Residuals                         12   6086     507                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_MS_lm_dist <- lm(row_sum ~ distance_to_ocean_min_m + Site_type, data = SR_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_MS_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_MS_lm_dist)
    ## W = 0.89383, p-value = 0.06405

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_MS_lm_dist) ~ SR_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  1  9.1557 0.009074 **
    ##       14                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_MS_lm_dist)
```

    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.858631          0.0022743     0.036388

``` r
# Plot residuals
plot(SD_SRic_MS_lm_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-16.png)<!-- -->

``` r
# Test relationship
SD_SRic_MS_am_dist <- aov(row_sum ~ distance_to_ocean_min_m + Site_type, data = SR_env[mixed_stratified_lakes,])

# Max depth
# Interaction
SD_SRic_MS_am_mxd <- aov(row_sum ~ max_depth * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_SRic_MS_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value  Pr(>F)    
    ## max_depth            1    394     394   1.127 0.30927    
    ## Site_type            1   8958    8958  25.599 0.00028 ***
    ## max_depth:Site_type  1   1340    1340   3.830 0.07401 .  
    ## Residuals           12   4199     350                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_MS_lm_mxd <- lm(row_sum ~ max_depth + Site_type, data = SR_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_MS_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_MS_lm_mxd)
    ## W = 0.93659, p-value = 0.3093

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_MS_lm_mxd) ~ SR_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  3.7709 0.07255 .
    ##       14                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_MS_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.598848          0.0036536     0.058458

``` r
# Plot residuals
plot(SD_SRic_MS_lm_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-17.png)<!-- -->

``` r
# Test relationship
SD_SRic_MS_am_mxd <- aov(row_sum ~ max_depth + Site_type, data = SR_env[mixed_stratified_lakes,])

# Log area
# Interaction
SD_SRic_MS_am_lga <- aov(row_sum ~ logArea * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_SRic_MS_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value  Pr(>F)    
    ## logArea            1    768     768   4.667 0.05167 .  
    ## Site_type          1  10409   10409  63.236   4e-06 ***
    ## logArea:Site_type  1   1739    1739  10.567 0.00695 ** 
    ## Residuals         12   1975     165                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_MS_lm_lga <- lm(row_sum ~ logArea + Site_type, data = SR_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_MS_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_MS_lm_lga)
    ## W = 0.96379, p-value = 0.7307

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_MS_lm_lga) ~ SR_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0996 0.7569
    ##       14

``` r
# Check for outliers
outlierTest(SD_SRic_MS_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.766241            0.01708      0.27329

``` r
# Plot residuals
plot(SD_SRic_MS_lm_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-18.png)<!-- -->

``` r
# Test relationship
SD_SRic_MS_am_lga <- aov(row_sum ~ logArea + Site_type, data = SR_env[mixed_stratified_lakes,])

# Anova outputs
summary(SD_SRic_MS_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m  1   4715    4715   10.03 0.00742 **
    ## Site_type                1   4066    4066    8.65 0.01147 * 
    ## Residuals               13   6111     470                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_MS_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth    1    394     394   0.926 0.353540    
    ## Site_type    1   8958    8958  21.023 0.000511 ***
    ## Residuals   13   5539     426                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_MS_am_lga)
```

    ##             Df Sum Sq Mean Sq F value  Pr(>F)    
    ## logArea      1    768     768   2.689   0.125    
    ## Site_type    1  10409   10409  36.428 4.2e-05 ***
    ## Residuals   13   3715     286                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
dist_p_values <- summary(SD_SRic_MS_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(SD_SRic_MS_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(SD_SRic_MS_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.0445452685 0.0688057806           NA 1.0000000000 0.0030683957
    ## [6]           NA 0.7501227093 0.0002517625           NA

``` r
## Ocean sites and mixed lakes
# Distance
# Interaction
SD_SRic_OM_am_dist <- aov(row_sum ~ distance_to_ocean_min_m * Site_type, data = SR_env[ocean_mixed_sites,])
summary(SD_SRic_OM_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m            1     14    13.9   0.016  0.903
    ## Site_type                          1     73    73.0   0.082  0.780
    ## distance_to_ocean_min_m:Site_type  1     13    13.2   0.015  0.905
    ## Residuals                         10   8895   889.5

``` r
# Linear model
SD_SRic_OM_lm_dist <- lm(row_sum ~ distance_to_ocean_min_m + Site_type, data = SR_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_OM_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_OM_lm_dist)
    ## W = 0.95157, p-value = 0.5853

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_OM_lm_dist) ~ SR_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.1659  0.691
    ##       12

``` r
# Check for outliers
outlierTest(SD_SRic_OM_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.836501           0.017654      0.24716

``` r
# Plot residuals
plot(SD_SRic_OM_lm_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-19.png)<!-- -->

``` r
# Test relationship
SD_SRic_OM_am_dist <- aov(row_sum ~ distance_to_ocean_min_m + Site_type, data = SR_env[ocean_mixed_sites,])

# Max depth
# Interaction
SD_SRic_OM_am_mxd <- aov(row_sum ~ max_depth * Site_type, data = SR_env[ocean_mixed_sites,])
summary(SD_SRic_OM_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value Pr(>F)  
    ## max_depth            1   3439    3439   6.298 0.0309 *
    ## Site_type            1     35      35   0.063 0.8062  
    ## max_depth:Site_type  1     61      61   0.112 0.7446  
    ## Residuals           10   5461     546                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_OM_lm_mxd <- lm(row_sum ~ max_depth + Site_type, data = SR_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_OM_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_OM_lm_mxd)
    ## W = 0.87725, p-value = 0.05311

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_OM_lm_mxd) ~ SR_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.2863 0.6024
    ##       12

``` r
# Check for outliers
outlierTest(SD_SRic_OM_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN  2.91868           0.015335      0.21469

``` r
# Plot residuals
plot(SD_SRic_OM_lm_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-20.png)<!-- -->

``` r
# Test relationship
SD_SRic_OM_am_mxd <- aov(row_sum ~ max_depth + Site_type, data = SR_env[ocean_mixed_sites,])

# Log area
# Interaction
SD_SRic_OM_am_lga <- aov(row_sum ~ logArea * Site_type, data = SR_env[ocean_mixed_sites,])
summary(SD_SRic_OM_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value Pr(>F)  
    ## logArea            1   1826  1826.2   4.105 0.0703 .
    ## Site_type          1    635   634.9   1.427 0.2598  
    ## logArea:Site_type  1   2086  2085.9   4.689 0.0556 .
    ## Residuals         10   4448   444.8                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_OM_lm_lga <- lm(row_sum ~ logArea + Site_type, data = SR_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_OM_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_OM_lm_lga)
    ## W = 0.96682, p-value = 0.8316

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_OM_lm_lga) ~ SR_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0786  0.784
    ##       12

``` r
# Check for outliers
outlierTest(SD_SRic_OM_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## OOO -2.89853           0.015874      0.22223

``` r
# Plot residuals
plot(SD_SRic_OM_lm_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-21.png)<!-- -->

``` r
# Test relationship
SD_SRic_OM_am_lga <- aov(row_sum ~ logArea + Site_type, data = SR_env[ocean_mixed_sites,])

# Anova outputs
summary(SD_SRic_OM_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1     14    13.9   0.017  0.898
    ## Site_type                1     73    73.0   0.090  0.770
    ## Residuals               11   8909   809.9

``` r
summary(SD_SRic_OM_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## max_depth    1   3439    3439   6.851 0.0239 *
    ## Site_type    1     35      35   0.069 0.7976  
    ## Residuals   11   5522     502                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_OM_am_lga)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## logArea      1   1826  1826.2   3.074  0.107
    ## Site_type    1    635   634.9   1.069  0.323
    ## Residuals   11   6534   594.0

``` r
# p-values
dist_p_values <- summary(SD_SRic_OM_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(SD_SRic_OM_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(SD_SRic_OM_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 1.0000000        NA 0.1436420 1.0000000        NA 0.6439722
    ## [8] 1.0000000        NA

``` r
## Stratified lakes and ocean sites
# Distance
# Interaction
SD_SRic_SO_am_dist <- aov(row_sum ~ distance_to_ocean_min_m * Site_type, data = SR_env[ocean_stratified_sites,])
summary(SD_SRic_SO_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m            1   5520    5520  17.550 0.00186 **
    ## Site_type                          1   1631    1631   5.186 0.04600 * 
    ## distance_to_ocean_min_m:Site_type  1      6       6   0.018 0.89653   
    ## Residuals                         10   3145     314                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_SO_lm_dist <- lm(row_sum ~ distance_to_ocean_min_m + Site_type, data = SR_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_SO_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_SO_lm_dist)
    ## W = 0.91486, p-value = 0.1854

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_SO_lm_dist) ~ SR_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  8.6917 0.01219 *
    ##       12                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_SO_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## IBK 2.801545           0.018745      0.26243

``` r
# Plot residuals
plot(SD_SRic_SO_lm_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-22.png)<!-- -->

``` r
# Test relationship
SD_SRic_SO_am_dist <- aov(row_sum ~ distance_to_ocean_min_m + Site_type, data = SR_env[ocean_stratified_sites,])

# Max depth
# Interaction
SD_SRic_SO_am_mxd <- aov(row_sum ~ max_depth * Site_type, data = SR_env[ocean_stratified_sites,])
summary(SD_SRic_SO_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth            1    119     119   0.655   0.4373    
    ## Site_type            1   7497    7497  41.343 7.54e-05 ***
    ## max_depth:Site_type  1    872     872   4.806   0.0531 .  
    ## Residuals           10   1813     181                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_SO_lm_mxd <- lm(row_sum ~ max_depth + Site_type, data = SR_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_SO_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_SO_lm_mxd)
    ## W = 0.96356, p-value = 0.781

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_SO_lm_mxd) ~ SR_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  3.5772 0.08296 .
    ##       12                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_SO_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## IBK  2.14316           0.057727      0.80817

``` r
# Plot residuals
plot(SD_SRic_SO_lm_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-23.png)<!-- -->

``` r
# Test relationship
SD_SRic_SO_am_mxd <- aov(row_sum ~ max_depth + Site_type, data = SR_env[ocean_stratified_sites,])

# Log area
# Interaction
SD_SRic_SO_am_lga <- aov(row_sum ~ logArea * Site_type, data = SR_env[ocean_stratified_sites,])
summary(SD_SRic_SO_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value  Pr(>F)   
    ## logArea            1   1872    1872   6.189 0.03211 * 
    ## Site_type          1   5342    5342  17.665 0.00182 **
    ## logArea:Site_type  1     62      62   0.206 0.65924   
    ## Residuals         10   3024     302                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_SRic_SO_lm_lga <- lm(row_sum ~ logArea + Site_type, data = SR_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_SO_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_SO_lm_lga)
    ## W = 0.95698, p-value = 0.6732

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_SO_lm_lga) ~ SR_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  4.3649 0.05865 .
    ##       12                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_SO_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OOO -3.217237          0.0092182      0.12906

``` r
# Plot residuals
plot(SD_SRic_SO_lm_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-24.png)<!-- -->

``` r
# Test relationship
SD_SRic_SO_am_lga <- aov(row_sum ~ logArea + Site_type, data = SR_env[ocean_stratified_sites,])

# Anova outputs
summary(SD_SRic_SO_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m  1   5520    5520  19.271 0.00108 **
    ## Site_type                1   1631    1631   5.694 0.03610 * 
    ## Residuals               11   3151     286                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_SO_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth    1    119     119   0.486 0.500047    
    ## Site_type    1   7497    7497  30.714 0.000175 ***
    ## Residuals   11   2685     244                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_SRic_SO_am_lga)
```

    ##             Df Sum Sq Mean Sq F value  Pr(>F)   
    ## logArea      1   1872    1872    6.67 0.02547 * 
    ## Site_type    1   5342    5342   19.04 0.00113 **
    ## Residuals   11   3087     281                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
dist_p_values <- summary(SD_SRic_SO_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(SD_SRic_SO_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(SD_SRic_SO_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.006487265 0.216612949          NA 1.000000000 0.001048991          NA
    ## [7] 0.152818514 0.006780374          NA

``` r
## Ocean sites
SD_SRic_geo_O_lm <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[ocean_sites,])
summary(SD_SRic_geo_O_lm)
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
p_values <- summary(SD_SRic_geo_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##                       1                       1                       1 
    ##                 logArea 
    ##                       1

``` r
## Mixed lakes
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
Anova(SD_SRic_geo_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: row_sum
    ##                          Sum Sq Df F value  Pr(>F)  
    ## (Intercept)             1272.73  1  3.1193 0.15212  
    ## distance_to_ocean_min_m   26.23  1  0.0643 0.81235  
    ## max_depth                  1.56  1  0.0038 0.95361  
    ## logArea                 1956.27  1  4.7946 0.09373 .
    ## Residuals               1632.06  4                  
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
## Stratified lakes
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
Anova(SD_SRic_geo_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: row_sum
    ##                          Sum Sq Df F value Pr(>F)
    ## (Intercept)               1.978  1  0.0483 0.8368
    ## distance_to_ocean_min_m 110.787  1  2.7051 0.1754
    ## max_depth                 0.019  1  0.0005 0.9839
    ## logArea                   0.951  1  0.0232 0.8862
    ## Residuals               163.819  4

``` r
p_values <- summary(SD_SRic_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##                 1.00000                 0.70149                 1.00000 
    ##                 logArea 
    ##                 1.00000

``` r
#### Environmental
### logSRic ANOVA
## Surveyed sites
# Temperature
# Interaction
SD_logSRic_am_temp <- aov(log(row_sum) ~ temperature_median * Site_type, data = SR_env[surveyed_sites_env,])
summary(SD_logSRic_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## temperature_median            1  2.599   2.599   5.440 0.036399 *  
    ## Site_type                     2 19.292   9.646  20.188 0.000103 ***
    ## temperature_median:Site_type  2  0.388   0.194   0.406 0.674310    
    ## Residuals                    13  6.211   0.478                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_lm_temp <- lm(log(row_sum) ~ temperature_median + Site_type, data = SR_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_lm_temp)
    ## W = 0.94972, p-value = 0.3909

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_lm_temp) ~ SR_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.7421 0.4918
    ##       16

``` r
# Check for outliers
outlierTest(SD_logSRic_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## HLM 2.038909           0.060805           NA

``` r
# Plot residuals
plot(SD_logSRic_lm_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-25.png)<!-- -->

``` r
# Test relationship
SD_logSRic_am_temp <- aov(log(row_sum) ~ temperature_median + Site_type, data = SR_env[surveyed_sites_env,])

# Salinity
# Interaction
SD_logSRic_am_sal <- aov(log(row_sum) ~ salinity_median * Site_type, data = SR_env[surveyed_sites_env,])
summary(SD_logSRic_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median            1 22.260  22.260  89.390 3.44e-07 ***
    ## Site_type                  2  2.358   1.179   4.735   0.0285 *  
    ## salinity_median:Site_type  2  0.635   0.317   1.275   0.3122    
    ## Residuals                 13  3.237   0.249                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_lm_sal <- lm(log(row_sum) ~ salinity_median + Site_type, data = SR_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_lm_sal)
    ## W = 0.95657, p-value = 0.5069

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_lm_sal) ~ SR_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.8224 0.4571
    ##       16

``` r
# Check for outliers
outlierTest(SD_logSRic_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OTM -1.854307           0.084879           NA

``` r
# Plot residuals
plot(SD_logSRic_lm_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-26.png)<!-- -->

``` r
# Test relationship
SD_logSRic_am_sal <- aov(log(row_sum) ~ salinity_median + Site_type, data = SR_env[surveyed_sites_env,])

# Oxygen
# Interaction
SD_logSRic_am_oxy <- aov(log(row_sum) ~ oxygen_median * Site_type, data = SR_env[surveyed_sites_env,])
summary(SD_logSRic_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median            1 12.415  12.415  43.041 1.82e-05 ***
    ## Site_type                2  9.779   4.890  16.952 0.000239 ***
    ## oxygen_median:Site_type  2  2.547   1.273   4.415 0.034424 *  
    ## Residuals               13  3.750   0.288                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_lm_oxy <- lm(log(row_sum) ~ oxygen_median + Site_type, data = SR_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_lm_oxy)
    ## W = 0.96235, p-value = 0.6194

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_lm_oxy) ~ SR_env[surveyed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.5387 0.5937
    ##       16

``` r
# Check for outliers
outlierTest(SD_logSRic_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## HLM 1.851206           0.085349           NA

``` r
# Plot residuals
plot(SD_logSRic_lm_oxy)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-27.png)<!-- -->

``` r
# Test relationship
SD_logSRic_am_oxy <- aov(log(row_sum) ~ oxygen_median + Site_type, data = SR_env[surveyed_sites_env,])

# Anova outputs
summary(SD_logSRic_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value   Pr(>F)    
    ## temperature_median  1  2.599   2.599   5.908   0.0281 *  
    ## Site_type           2 19.292   9.646  21.924 3.53e-05 ***
    ## Residuals          15  6.600   0.440                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median  1 22.260  22.260  86.230 1.31e-07 ***
    ## Site_type        2  2.358   1.179   4.568   0.0282 *  
    ## Residuals       15  3.872   0.258                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median  1 12.415   12.41   29.58 6.85e-05 ***
    ## Site_type      2  9.779    4.89   11.65 0.000885 ***
    ## Residuals     15  6.297    0.42                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(SD_logSRic_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(SD_logSRic_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(SD_logSRic_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.685435e-01 2.117708e-04           NA 7.873462e-07 1.693822e-01
    ## [6]           NA 4.112906e-04 5.310782e-03           NA

``` r
## Mixed and stratified lakes
# Temperature
# Interaction
SD_logSRic_MS_am_temp <- aov(log(row_sum) ~ temperature_median * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value   Pr(>F)    
    ## temperature_median            1  2.891   2.891   5.598 0.035656 *  
    ## Site_type                     1 14.497  14.497  28.071 0.000189 ***
    ## temperature_median:Site_type  1  0.199   0.199   0.385 0.546604    
    ## Residuals                    12  6.198   0.516                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_MS_lm_temp <- lm(log(row_sum) ~ temperature_median + Site_type, data = SR_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_MS_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_MS_lm_temp)
    ## W = 0.93532, p-value = 0.2956

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_MS_lm_temp) ~ SR_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.1222 0.7319
    ##       14

``` r
# Check for outliers
outlierTest(SD_logSRic_MS_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -2.182316           0.049687        0.795

``` r
# Plot residuals
plot(SD_logSRic_MS_lm_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-28.png)<!-- -->

``` r
# Test relationship
SD_logSRic_MS_am_temp <- aov(log(row_sum) ~ temperature_median + Site_type, data = SR_env[mixed_stratified_lakes,])

# Salinity
# Interaction
SD_logSRic_MS_am_sal <- aov(log(row_sum) ~ salinity_median * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median            1 18.016  18.016  70.319 2.31e-06 ***
    ## Site_type                  1  2.060   2.060   8.039    0.015 *  
    ## salinity_median:Site_type  1  0.635   0.635   2.478    0.141    
    ## Residuals                 12  3.074   0.256                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_MS_lm_sal <- lm(log(row_sum) ~ salinity_median + Site_type, data = SR_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_MS_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_MS_lm_sal)
    ## W = 0.94056, p-value = 0.3559

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_MS_lm_sal) ~ SR_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.5694  0.463
    ##       14

``` r
# Check for outliers
outlierTest(SD_logSRic_MS_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OTM -1.763547            0.10323           NA

``` r
# Plot residuals
plot(SD_logSRic_MS_lm_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-29.png)<!-- -->

``` r
# Test relationship
SD_logSRic_MS_am_sal <- aov(log(row_sum) ~ salinity_median + Site_type, data = SR_env[mixed_stratified_lakes,])

# Oxygen
# Interaction
SD_logSRic_MS_am_oxy <- aov(log(row_sum) ~ oxygen_median * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median            1  7.797   7.797  24.952 0.000312 ***
    ## Site_type                1  9.969   9.969  31.903 0.000108 ***
    ## oxygen_median:Site_type  1  2.269   2.269   7.263 0.019490 *  
    ## Residuals               12  3.750   0.312                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_MS_lm_oxy <- lm(log(row_sum) ~ oxygen_median + Site_type, data = SR_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_MS_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_MS_lm_oxy)
    ## W = 0.9435, p-value = 0.3941

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_MS_lm_oxy) ~ SR_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.1231 0.7309
    ##       14

``` r
# Check for outliers
outlierTest(SD_logSRic_MS_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## HLM 1.714236            0.11217           NA

``` r
# Plot residuals
plot(SD_logSRic_MS_lm_oxy)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-30.png)<!-- -->

``` r
# Test relationship
SD_logSRic_MS_am_oxy <- aov(log(row_sum) ~ oxygen_median + Site_type, data = SR_env[mixed_stratified_lakes,])

# Anova outputs
summary(SD_logSRic_MS_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value   Pr(>F)    
    ## temperature_median  1  2.891   2.891   5.876 0.030668 *  
    ## Site_type           1 14.497  14.497  29.465 0.000115 ***
    ## Residuals          13  6.396   0.492                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_MS_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median  1 18.016  18.016  63.141 2.41e-06 ***
    ## Site_type        1  2.060   2.060   7.218   0.0187 *  
    ## Residuals       13  3.709   0.285                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_MS_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median  1  7.797   7.797   16.84 0.001245 ** 
    ## Site_type      1  9.969   9.969   21.53 0.000463 ***
    ## Residuals     13  6.019   0.463                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(SD_logSRic_MS_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(SD_logSRic_MS_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(SD_logSRic_MS_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.1840095285 0.0006926945           NA 0.0000144461 0.1119775178
    ## [6]           NA 0.0074681089 0.0027764148           NA

``` r
## Ocean sites and mixed lakes
# Temperature
# Interaction
SD_logSRic_OM_am_temp <- aov(log(row_sum) ~ temperature_median * Site_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_OM_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median            1 0.0031  0.0031   0.009  0.927
    ## Site_type                     1 0.2209  0.2209   0.637  0.451
    ## temperature_median:Site_type  1 0.3722  0.3722   1.074  0.335
    ## Residuals                     7 2.4265  0.3466

``` r
# Linear model
SD_logSRic_OM_lm_temp <- lm(log(row_sum) ~ temperature_median + Site_type, data = SR_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_OM_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_OM_lm_temp)
    ## W = 0.91473, p-value = 0.2771

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_OM_lm_temp) ~ SR_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.8649 0.3766
    ##        9

``` r
# Check for outliers
outlierTest(SD_logSRic_OM_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OLO -1.798868            0.11507           NA

``` r
# Plot residuals
plot(SD_logSRic_OM_lm_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-31.png)<!-- -->

``` r
# Test relationship
SD_logSRic_OM_am_temp <- aov(log(row_sum) ~ temperature_median + Site_type, data = SR_env[ocean_mixed_sites_env,])

# Salinity
# Interaction
SD_logSRic_OM_am_sal <- aov(log(row_sum) ~ salinity_median * Site_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_OM_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median            1 1.1092  1.1092   4.203 0.0795 .
    ## Site_type                  1 0.0513  0.0513   0.194 0.6727  
    ## salinity_median:Site_type  1 0.0150  0.0150   0.057 0.8181  
    ## Residuals                  7 1.8472  0.2639                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_OM_lm_sal <- lm(log(row_sum) ~ salinity_median + Site_type, data = SR_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_OM_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_OM_lm_sal)
    ## W = 0.95831, p-value = 0.7505

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_OM_lm_sal) ~ SR_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  1.0283  0.337
    ##        9

``` r
# Check for outliers
outlierTest(SD_logSRic_OM_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 1.867974            0.10399           NA

``` r
# Plot residuals
plot(SD_logSRic_OM_lm_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-32.png)<!-- -->

``` r
# Test relationship
SD_logSRic_OM_am_sal <- aov(log(row_sum) ~ salinity_median + Site_type, data = SR_env[ocean_mixed_sites_env,])

# Oxygen
# Interaction
SD_logSRic_OM_am_oxy <- aov(log(row_sum) ~ oxygen_median * Site_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_OM_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)  
    ## oxygen_median            1 1.1393  1.1393   4.518 0.0711 .
    ## Site_type                1 0.1166  0.1166   0.462 0.5184  
    ## oxygen_median:Site_type  1 0.0014  0.0014   0.005 0.9430  
    ## Residuals                7 1.7653  0.2522                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_OM_lm_oxy <- lm(log(row_sum) ~ oxygen_median + Site_type, data = SR_env[ocean_mixed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_OM_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_OM_lm_oxy)
    ## W = 0.93881, p-value = 0.5066

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_OM_lm_oxy) ~ SR_env[ocean_mixed_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  1  4.1659 0.07164 .
    ##        9                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_logSRic_OM_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OLO -3.035435           0.018969      0.20865

``` r
# Plot residuals
plot(SD_logSRic_OM_lm_oxy)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-33.png)<!-- -->

``` r
# Test relationship
SD_logSRic_OM_am_oxy <- aov(log(row_sum) ~ oxygen_median + Site_type, data = SR_env[ocean_mixed_sites_env,])

# Anova outputs
summary(SD_logSRic_OM_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value Pr(>F)
    ## temperature_median  1 0.0031  0.0031   0.009  0.927
    ## Site_type           1 0.2209  0.2209   0.631  0.450
    ## Residuals           8 2.7986  0.3498

``` r
summary(SD_logSRic_OM_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value Pr(>F)  
    ## salinity_median  1 1.1092  1.1092   4.765 0.0606 .
    ## Site_type        1 0.0513  0.0513   0.220 0.6514  
    ## Residuals        8 1.8622  0.2328                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OM_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)  
    ## oxygen_median  1 1.1393  1.1393   5.159 0.0528 .
    ## Site_type      1 0.1166  0.1166   0.528 0.4881  
    ## Residuals      8 1.7667  0.2208                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(SD_logSRic_OM_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(SD_logSRic_OM_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(SD_logSRic_OM_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000 1.0000000        NA 0.3635379 1.0000000        NA 0.3166768
    ## [8] 1.0000000        NA

``` r
## Stratified lakes and ocean sites
# Temperature
# Interaction
SD_logSRic_SO_am_temp <- aov(log(row_sum) ~ temperature_median * Site_type, data = SR_env[ocean_stratified_sites_env,])
summary(SD_logSRic_SO_am_temp)
```

    ##                              Df Sum Sq Mean Sq F value  Pr(>F)   
    ## temperature_median            1  0.303   0.303   0.558 0.47939   
    ## Site_type                     1 12.076  12.076  22.252 0.00216 **
    ## temperature_median:Site_type  1  0.150   0.150   0.277 0.61486   
    ## Residuals                     7  3.799   0.543                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_SO_lm_temp <- lm(log(row_sum) ~ temperature_median + Site_type, data = SR_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_SO_lm_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_SO_lm_temp)
    ## W = 0.85487, p-value = 0.0494

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_SO_lm_temp) ~ SR_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  1.3348 0.2777
    ##        9

``` r
# Check for outliers
outlierTest(SD_logSRic_SO_lm_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## HLM 2.254231            0.05883      0.64713

``` r
# Plot residuals
plot(SD_logSRic_SO_lm_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-34.png)<!-- -->

``` r
# Test relationship
SD_logSRic_SO_am_temp <- aov(log(row_sum) ~ temperature_median + Site_type, data = SR_env[ocean_stratified_sites_env,])

# Salinity
# Interaction
SD_logSRic_SO_am_sal <- aov(log(row_sum) ~ salinity_median * Site_type, data = SR_env[ocean_stratified_sites_env,])
summary(SD_logSRic_SO_am_sal)
```

    ##                           Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median            1 12.655  12.655  57.041 0.000131 ***
    ## Site_type                  1  2.120   2.120   9.555 0.017539 *  
    ## salinity_median:Site_type  1  0.000   0.000   0.000 0.982859    
    ## Residuals                  7  1.553   0.222                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_SO_lm_sal <- lm(log(row_sum) ~ salinity_median + Site_type, data = SR_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_SO_lm_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_SO_lm_sal)
    ## W = 0.9724, p-value = 0.9097

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_SO_lm_sal) ~ SR_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.5406 0.4809
    ##        9

``` r
# Check for outliers
outlierTest(SD_logSRic_SO_lm_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OTM -2.545489           0.038355      0.42191

``` r
# Plot residuals
plot(SD_logSRic_SO_lm_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-35.png)<!-- -->

``` r
# Test relationship
SD_logSRic_SO_am_sal <- aov(log(row_sum) ~ salinity_median + Site_type, data = SR_env[ocean_stratified_sites_env,])

# Oxygen
# Interaction
SD_logSRic_SO_am_oxy <- aov(log(row_sum) ~ oxygen_median * Site_type, data = SR_env[ocean_stratified_sites_env,])
summary(SD_logSRic_SO_am_oxy)
```

    ##                         Df Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median            1  5.357   5.357  18.894 0.003368 ** 
    ## Site_type                1  8.495   8.495  29.962 0.000932 ***
    ## oxygen_median:Site_type  1  0.492   0.492   1.736 0.229168    
    ## Residuals                7  1.985   0.284                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_SO_lm_oxy <- lm(log(row_sum) ~ oxygen_median + Site_type, data = SR_env[ocean_stratified_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_SO_lm_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_SO_lm_oxy)
    ## W = 0.98379, p-value = 0.9835

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_SO_lm_oxy) ~ SR_env[ocean_stratified_sites_env,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.1765 0.6842
    ##        9

``` r
# Check for outliers
outlierTest(SD_logSRic_SO_lm_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OTM -2.315882           0.053719      0.59091

``` r
# Plot residuals
plot(SD_logSRic_SO_lm_oxy)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-36.png)<!-- -->

``` r
# Test relationship
SD_logSRic_SO_am_oxy <- aov(log(row_sum) ~ oxygen_median + Site_type, data = SR_env[ocean_stratified_sites_env,])

# Anova outputs
summary(SD_logSRic_SO_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value  Pr(>F)   
    ## temperature_median  1  0.303   0.303   0.613 0.45604   
    ## Site_type           1 12.076  12.076  24.463 0.00113 **
    ## Residuals           8  3.949   0.494                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_SO_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median  1 12.655  12.655   65.19 4.09e-05 ***
    ## Site_type        1  2.120   2.120   10.92   0.0108 *  
    ## Residuals        8  1.553   0.194                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_SO_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median  1  5.357   5.357   17.30 0.003167 ** 
    ## Site_type      1  8.495   8.495   27.44 0.000785 ***
    ## Residuals      8  2.477   0.310                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
temp_p_values <- summary(SD_logSRic_SO_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(SD_logSRic_SO_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(SD_logSRic_SO_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.0000000000 0.0067595020           NA 0.0002452472 0.0647207568
    ## [6]           NA 0.0189993093 0.0047099840           NA

``` r
## Ocean sites
SD_logSRic_env_O_lm <- lm(log(row_sum) ~ salinity_median + oxygen_median + temperature_median, SR_env[ocean_sites_env,])
summary(SD_logSRic_env_O_lm)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ salinity_median + oxygen_median + 
    ##     temperature_median, data = SR_env[ocean_sites_env, ])
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
    ## Multiple R-squared:      1,  Adjusted R-squared:    NaN 
    ## F-statistic:   NaN on 2 and 0 DF,  p-value: NA

``` r
p_values <- summary(SD_logSRic_env_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##     (Intercept) salinity_median   oxygen_median 
    ##             NaN             NaN             NaN

``` r
## Mixed lakes
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
Anova(SD_logSRic_env_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: log(row_sum)
    ##                     Sum Sq Df F value Pr(>F)
    ## (Intercept)        0.09631  1  0.2889 0.6194
    ## salinity_median    0.27660  1  0.8298 0.4139
    ## oxygen_median      0.34222  1  1.0267 0.3683
    ## temperature_median 0.02010  1  0.0603 0.8181
    ## Residuals          1.33334  4

``` r
p_values <- summary(SD_logSRic_env_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##                  1                  1                  1                  1

``` r
## Stratified lakes
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
Anova(SD_logSRic_env_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: log(row_sum)
    ##                     Sum Sq Df F value Pr(>F)
    ## (Intercept)        0.08340  1  0.2587 0.6378
    ## salinity_median    0.69323  1  2.1503 0.2164
    ## oxygen_median      0.04881  1  0.1514 0.7170
    ## temperature_median 0.04124  1  0.1279 0.7387
    ## Residuals          1.28953  4

``` r
p_values <- summary(SD_logSRic_env_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##        (Intercept)    salinity_median      oxygen_median temperature_median 
    ##          1.0000000          0.8657002          1.0000000          1.0000000

``` r
#### Geographical
### logSRic ANOVA
## Surveyed sites
# Distance
# Interaction
SD_logSRic_am_dist <- aov(log(row_sum) ~ distance_to_ocean_min_m * Site_type, data = SR_env[surveyed_sites,])
summary(SD_logSRic_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m            1 17.419  17.419  45.023 5.01e-06 ***
    ## Site_type                          2  6.414   3.207   8.289  0.00339 ** 
    ## distance_to_ocean_min_m:Site_type  2  0.068   0.034   0.087  0.91689    
    ## Residuals                         16  6.190   0.387                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_lm_dist <- lm(log(row_sum) ~ distance_to_ocean_min_m + Site_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_lm_dist)
    ## W = 0.9606, p-value = 0.5015

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_lm_dist) ~ SR_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.1268 0.8816
    ##       19

``` r
# Check for outliers
outlierTest(SD_logSRic_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 1.941323           0.068983           NA

``` r
# Plot residuals
plot(SD_logSRic_lm_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-37.png)<!-- -->

``` r
# Test relationship
SD_logSRic_am_dist <- aov(log(row_sum) ~ distance_to_ocean_min_m + Site_type, data = SR_env[surveyed_sites,])

# Max depth
# Interaction
SD_logSRic_am_mxd <- aov(log(row_sum) ~ max_depth * Site_type, data = SR_env[surveyed_sites,])
summary(SD_logSRic_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth            1  1.336   1.336   3.599    0.076 .  
    ## Site_type            2 21.997  10.998  29.636 4.17e-06 ***
    ## max_depth:Site_type  2  0.820   0.410   1.105    0.355    
    ## Residuals           16  5.938   0.371                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_lm_mxd <- lm(log(row_sum) ~ max_depth + Site_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_lm_mxd)
    ## W = 0.96227, p-value = 0.5365

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_lm_mxd) ~ SR_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.6256 0.5456
    ##       19

``` r
# Check for outliers
outlierTest(SD_logSRic_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## TLN 3.003212          0.0079997      0.17599

``` r
# Plot residuals
plot(SD_logSRic_lm_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-38.png)<!-- -->

``` r
# Test relationship
SD_logSRic_am_mxd <- aov(log(row_sum) ~ max_depth + Site_type, data = SR_env[surveyed_sites,])

# Log area
# Interaction
SD_logSRic_am_lga <- aov(log(row_sum) ~ logArea * Site_type, data = SR_env[surveyed_sites,])
summary(SD_logSRic_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## logArea            1  1.167   1.167   3.179   0.0936 .  
    ## Site_type          2 21.837  10.918  29.740 4.08e-06 ***
    ## logArea:Site_type  2  1.213   0.607   1.652   0.2227    
    ## Residuals         16  5.874   0.367                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_lm_lga <- lm(log(row_sum) ~ logArea + Site_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_lm_lga)
    ## W = 0.96448, p-value = 0.5846

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_lm_lga) ~ SR_env[surveyed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.3441 0.7132
    ##       19

``` r
# Check for outliers
outlierTest(SD_logSRic_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## TLN 2.402145           0.028005       0.6161

``` r
# Plot residuals
plot(SD_logSRic_lm_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-39.png)<!-- -->

``` r
# Test relationship
SD_logSRic_am_lga <- aov(log(row_sum) ~ logArea + Site_type, data = SR_env[surveyed_sites,])

# Anova outputs
summary(SD_logSRic_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m  1 17.419  17.419  50.105 1.34e-06 ***
    ## Site_type                2  6.414   3.207   9.225  0.00175 ** 
    ## Residuals               18  6.258   0.348                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth    1  1.336   1.336   3.558   0.0755 .  
    ## Site_type    2 21.997  10.998  29.294 2.19e-06 ***
    ## Residuals   18  6.758   0.375                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_am_lga)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## logArea      1  1.167   1.167   2.964    0.102    
    ## Site_type    2 21.837  10.918  27.731 3.18e-06 ***
    ## Residuals   18  7.087   0.394                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
dist_p_values <- summary(SD_logSRic_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(SD_logSRic_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(SD_logSRic_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 8.024381e-06 1.048152e-02           NA 4.530155e-01 1.312572e-05
    ## [6]           NA 6.136847e-01 1.909931e-05           NA

``` r
## Mixed and stratified lakes
# Distance
# Interaction
SD_logSRic_MS_am_dist <- aov(log(row_sum) ~ distance_to_ocean_min_m * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m            1 12.814  12.814  31.294 0.000117 ***
    ## Site_type                          1  6.055   6.055  14.788 0.002329 ** 
    ## distance_to_ocean_min_m:Site_type  1  0.003   0.003   0.007 0.934494    
    ## Residuals                         12  4.914   0.409                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_MS_lm_dist <- lm(log(row_sum) ~ distance_to_ocean_min_m + Site_type, data = SR_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_MS_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_MS_lm_dist)
    ## W = 0.96462, p-value = 0.7457

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_MS_lm_dist) ~ SR_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0556  0.817
    ##       14

``` r
# Check for outliers
outlierTest(SD_logSRic_MS_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 1.906788            0.08077           NA

``` r
# Plot residuals
plot(SD_logSRic_MS_lm_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-40.png)<!-- -->

``` r
# Test relationship
SD_logSRic_MS_am_dist <- aov(log(row_sum) ~ distance_to_ocean_min_m + Site_type, data = SR_env[mixed_stratified_lakes,])

# Max depth
# Interaction
SD_logSRic_MS_am_mxd <- aov(log(row_sum) ~ max_depth * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth            1  1.817   1.817   4.215   0.0625 .  
    ## Site_type            1 16.023  16.023  37.164 5.37e-05 ***
    ## max_depth:Site_type  1  0.771   0.771   1.788   0.2059    
    ## Residuals           12  5.174   0.431                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_MS_lm_mxd <- lm(log(row_sum) ~ max_depth + Site_type, data = SR_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_MS_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_MS_lm_mxd)
    ## W = 0.94087, p-value = 0.3598

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_MS_lm_mxd) ~ SR_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.3162 0.5828
    ##       14

``` r
# Check for outliers
outlierTest(SD_logSRic_MS_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## TLN  2.75529           0.017431      0.27889

``` r
# Plot residuals
plot(SD_logSRic_MS_lm_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-41.png)<!-- -->

``` r
# Test relationship
SD_logSRic_MS_am_mxd <- aov(log(row_sum) ~ max_depth + Site_type, data = SR_env[mixed_stratified_lakes,])

# Log area
# Interaction
SD_logSRic_MS_am_lga <- aov(log(row_sum) ~ logArea * Site_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## logArea            1  0.004   0.004   0.010    0.924    
    ## Site_type          1 18.456  18.456  48.133 1.57e-05 ***
    ## logArea:Site_type  1  0.724   0.724   1.888    0.195    
    ## Residuals         12  4.601   0.383                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_MS_lm_lga <- lm(log(row_sum) ~ logArea + Site_type, data = SR_env[mixed_stratified_lakes,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_MS_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_MS_lm_lga)
    ## W = 0.93763, p-value = 0.321

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_MS_lm_lga) ~ SR_env[mixed_stratified_lakes,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1   1.588 0.2282
    ##       14

``` r
# Check for outliers
outlierTest(SD_logSRic_MS_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## TLN 3.109747          0.0090251       0.1444

``` r
# Plot residuals
plot(SD_logSRic_MS_lm_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-42.png)<!-- -->

``` r
# Test relationship
SD_logSRic_MS_am_lga <- aov(log(row_sum) ~ logArea + Site_type, data = SR_env[mixed_stratified_lakes,])

# Anova outputs
summary(SD_logSRic_MS_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m  1 12.814  12.814   33.88 5.97e-05 ***
    ## Site_type                1  6.055   6.055   16.01  0.00151 ** 
    ## Residuals               13  4.916   0.378                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_MS_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth    1  1.817   1.817   3.974   0.0676 .  
    ## Site_type    1 16.023  16.023  35.040 5.07e-05 ***
    ## Residuals   13  5.945   0.457                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_MS_am_lga)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## logArea      1  0.004   0.004   0.009    0.926    
    ## Site_type    1 18.456  18.456  45.056 1.44e-05 ***
    ## Residuals   13  5.325   0.410                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
dist_p_values <- summary(SD_logSRic_MS_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(SD_logSRic_MS_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(SD_logSRic_MS_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 3.579803e-04 9.049388e-03           NA 4.057944e-01 3.043065e-04
    ## [6]           NA 1.000000e+00 8.657848e-05           NA

``` r
## Ocean sites and mixed lakes
# Distance
# Interaction
SD_logSRic_OM_am_dist <- aov(log(row_sum) ~ distance_to_ocean_min_m * Site_type, data = SR_env[ocean_mixed_sites,])
summary(SD_logSRic_OM_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m            1  0.071  0.0712   0.190  0.672
    ## Site_type                          1  0.076  0.0763   0.203  0.662
    ## distance_to_ocean_min_m:Site_type  1  0.055  0.0554   0.147  0.709
    ## Residuals                         10  3.757  0.3757

``` r
# Linear model
SD_logSRic_OM_lm_dist <- lm(log(row_sum) ~ distance_to_ocean_min_m + Site_type, data = SR_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_OM_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_OM_lm_dist)
    ## W = 0.94615, p-value = 0.5028

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_OM_lm_dist) ~ SR_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.1052 0.7513
    ##       12

``` r
# Check for outliers
outlierTest(SD_logSRic_OM_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.127304            0.05929      0.83006

``` r
# Plot residuals
plot(SD_logSRic_OM_lm_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-43.png)<!-- -->

``` r
# Test relationship
SD_logSRic_OM_am_dist <- aov(log(row_sum) ~ distance_to_ocean_min_m + Site_type, data = SR_env[ocean_mixed_sites,])

# Max depth
# Interaction
SD_logSRic_OM_am_mxd <- aov(log(row_sum) ~ max_depth * Site_type, data = SR_env[ocean_mixed_sites,])
summary(SD_logSRic_OM_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value Pr(>F)  
    ## max_depth            1 1.6952  1.6952   7.879 0.0186 *
    ## Site_type            1 0.0010  0.0010   0.004 0.9483  
    ## max_depth:Site_type  1 0.1124  0.1124   0.522 0.4864  
    ## Residuals           10 2.1514  0.2151                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_OM_lm_mxd <- lm(log(row_sum) ~ max_depth + Site_type, data = SR_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_OM_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_OM_lm_mxd)
    ## W = 0.92331, p-value = 0.2453

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_OM_lm_mxd) ~ SR_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.2914 0.5992
    ##       12

``` r
# Check for outliers
outlierTest(SD_logSRic_OM_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 1.870151           0.090985           NA

``` r
# Plot residuals
plot(SD_logSRic_OM_lm_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-44.png)<!-- -->

``` r
# Test relationship
SD_logSRic_OM_am_mxd <- aov(log(row_sum) ~ max_depth + Site_type, data = SR_env[ocean_mixed_sites,])

# Log area
# Interaction
SD_logSRic_OM_am_lga <- aov(log(row_sum) ~ logArea * Site_type, data = SR_env[ocean_mixed_sites,])
summary(SD_logSRic_OM_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value Pr(>F)  
    ## logArea            1 0.6565  0.6565   3.147 0.1065  
    ## Site_type          1 0.1381  0.1381   0.662 0.4348  
    ## logArea:Site_type  1 1.0791  1.0791   5.173 0.0462 *
    ## Residuals         10 2.0862  0.2086                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_OM_lm_lga <- lm(log(row_sum) ~ logArea + Site_type, data = SR_env[ocean_mixed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_OM_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_OM_lm_lga)
    ## W = 0.88016, p-value = 0.05841

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_OM_lm_lga) ~ SR_env[ocean_mixed_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.0289 0.8679
    ##       12

``` r
# Check for outliers
outlierTest(SD_logSRic_OM_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OOO -2.705727           0.022095      0.30933

``` r
# Plot residuals
plot(SD_logSRic_OM_lm_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-45.png)<!-- -->

``` r
# Test relationship
SD_logSRic_OM_am_lga <- aov(log(row_sum) ~ logArea + Site_type, data = SR_env[ocean_mixed_sites,])

# Anova outputs
summary(SD_logSRic_OM_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)
    ## distance_to_ocean_min_m  1  0.071  0.0712   0.206  0.659
    ## Site_type                1  0.076  0.0763   0.220  0.648
    ## Residuals               11  3.812  0.3466

``` r
summary(SD_logSRic_OM_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)  
    ## max_depth    1  1.695  1.6952   8.237 0.0152 *
    ## Site_type    1  0.001  0.0010   0.005 0.9470  
    ## Residuals   11  2.264  0.2058                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OM_am_lga)
```

    ##             Df Sum Sq Mean Sq F value Pr(>F)
    ## logArea      1  0.656  0.6565   2.281  0.159
    ## Site_type    1  0.138  0.1381   0.480  0.503
    ## Residuals   11  3.165  0.2878

``` r
# p-values
dist_p_values <- summary(SD_logSRic_OM_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(SD_logSRic_OM_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(SD_logSRic_OM_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.00000000 1.00000000         NA 0.09144026 1.00000000         NA 0.95466089
    ## [8] 1.00000000         NA

``` r
## Stratified lakes and ocean sites
# Distance
# Interaction
SD_logSRic_SO_am_dist <- aov(log(row_sum) ~ distance_to_ocean_min_m * Site_type, data = SR_env[ocean_stratified_sites,])
summary(SD_logSRic_SO_am_dist)
```

    ##                                   Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m            1 14.817  14.817  39.938 8.68e-05 ***
    ## Site_type                          1  1.791   1.791   4.827   0.0527 .  
    ## distance_to_ocean_min_m:Site_type  1  0.066   0.066   0.177   0.6830    
    ## Residuals                         10  3.710   0.371                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_SO_lm_dist <- lm(log(row_sum) ~ distance_to_ocean_min_m + Site_type, data = SR_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_SO_lm_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_SO_lm_dist)
    ## W = 0.94979, p-value = 0.5574

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_SO_lm_dist) ~ SR_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.2878 0.6014
    ##       12

``` r
# Check for outliers
outlierTest(SD_logSRic_SO_lm_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OTM -2.041337           0.068494      0.95892

``` r
# Plot residuals
plot(SD_logSRic_SO_lm_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-46.png)<!-- -->

``` r
# Test relationship
SD_logSRic_SO_am_dist <- aov(log(row_sum) ~ distance_to_ocean_min_m + Site_type, data = SR_env[ocean_stratified_sites,])

# Max depth
# Interaction
SD_logSRic_SO_am_mxd <- aov(log(row_sum) ~ max_depth * Site_type, data = SR_env[ocean_stratified_sites,])
summary(SD_logSRic_SO_am_mxd)
```

    ##                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth            1  1.166   1.166   2.562 0.140513    
    ## Site_type            1 14.375  14.375  31.590 0.000221 ***
    ## max_depth:Site_type  1  0.291   0.291   0.641 0.442124    
    ## Residuals           10  4.550   0.455                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_SO_lm_mxd <- lm(log(row_sum) ~ max_depth + Site_type, data = SR_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_SO_lm_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_SO_lm_mxd)
    ## W = 0.93288, p-value = 0.3347

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_SO_lm_mxd) ~ SR_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.6516 0.4353
    ##       12

``` r
# Check for outliers
outlierTest(SD_logSRic_SO_lm_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## TLN 2.669414           0.023516      0.32923

``` r
# Plot residuals
plot(SD_logSRic_SO_lm_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-47.png)<!-- -->

``` r
# Test relationship
SD_logSRic_SO_am_mxd <- aov(log(row_sum) ~ max_depth + Site_type, data = SR_env[ocean_stratified_sites,])

# Log area
# Interaction
SD_logSRic_SO_am_lga <- aov(log(row_sum) ~ logArea * Site_type, data = SR_env[ocean_stratified_sites,])
summary(SD_logSRic_SO_am_lga)
```

    ##                   Df Sum Sq Mean Sq F value  Pr(>F)    
    ## logArea            1  2.487   2.487   4.914 0.05097 .  
    ## Site_type          1 12.831  12.831  25.354 0.00051 ***
    ## logArea:Site_type  1  0.005   0.005   0.010 0.92382    
    ## Residuals         10  5.061   0.506                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Linear model
SD_logSRic_SO_lm_lga <- lm(log(row_sum) ~ logArea + Site_type, data = SR_env[ocean_stratified_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_logSRic_SO_lm_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_logSRic_SO_lm_lga)
    ## W = 0.89998, p-value = 0.1127

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_logSRic_SO_lm_lga) ~ SR_env[ocean_stratified_sites,"Site_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  1  0.4666 0.5075
    ##       12

``` r
# Check for outliers
outlierTest(SD_logSRic_SO_lm_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## TLN 2.031971           0.069575      0.97405

``` r
# Plot residuals
plot(SD_logSRic_SO_lm_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANOVA-48.png)<!-- -->

``` r
# Test relationship
SD_logSRic_SO_am_lga <- aov(log(row_sum) ~ logArea + Site_type, data = SR_env[ocean_stratified_sites,])

# Anova outputs
summary(SD_logSRic_SO_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value   Pr(>F)    
    ## distance_to_ocean_min_m  1 14.817  14.817  43.168 4.02e-05 ***
    ## Site_type                1  1.791   1.791   5.218   0.0432 *  
    ## Residuals               11  3.776   0.343                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_SO_am_mxd)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth    1  1.166   1.166   2.649 0.131893    
    ## Site_type    1 14.375  14.375  32.658 0.000135 ***
    ## Residuals   11  4.842   0.440                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_SO_am_lga)
```

    ##             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## logArea      1  2.487   2.487   5.401 0.040297 *  
    ## Site_type    1 12.831  12.831  27.863 0.000261 ***
    ## Residuals   11  5.065   0.460                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# p-values
dist_p_values <- summary(SD_logSRic_SO_am_dist)[[1]][, "Pr(>F)"]
mxd_p_values <- summary(SD_logSRic_SO_am_mxd)[[1]][, "Pr(>F)"]
lga_p_values <- summary(SD_logSRic_SO_am_lga)[[1]][, "Pr(>F)"]
p_values <- c(dist_p_values[1:3], mxd_p_values[1:3], lga_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.0002413784 0.2592921788           NA 0.7913550208 0.0008113704
    ## [6]           NA 0.2417842531 0.0015652516           NA

``` r
## Ocean sites
SD_logSRic_geo_O_lm <- lm(log(row_sum) ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[ocean_sites,])
summary(SD_logSRic_geo_O_lm)
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
p_values <- summary(SD_logSRic_geo_O_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               0.5315567               1.0000000               1.0000000 
    ##                 logArea 
    ##               1.0000000

``` r
## Mixed lakes
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
Anova(SD_logSRic_geo_M_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: log(row_sum)
    ##                          Sum Sq Df F value Pr(>F)
    ## (Intercept)             0.11806  1  0.7369 0.4391
    ## distance_to_ocean_min_m 0.03884  1  0.2424 0.6483
    ## max_depth               0.01821  1  0.1137 0.7530
    ## logArea                 0.69363  1  4.3293 0.1059
    ## Residuals               0.64087  4

``` r
p_values <- summary(SD_logSRic_geo_M_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               1.0000000               1.0000000 
    ##                 logArea 
    ##               0.4237593

``` r
## Stratified lakes
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
Anova(SD_logSRic_geo_S_lm, type = 3)
```

    ## Anova Table (Type III tests)
    ## 
    ## Response: log(row_sum)
    ##                          Sum Sq Df F value Pr(>F)
    ## (Intercept)             0.05931  1  0.1014 0.7661
    ## distance_to_ocean_min_m 1.44031  1  2.4618 0.1917
    ## max_depth               0.00162  1  0.0028 0.9605
    ## logArea                 0.02762  1  0.0472 0.8386
    ## Residuals               2.34024  4

``` r
p_values <- summary(SD_logSRic_geo_S_lm)$coefficients[, "Pr(>|t|)"]
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ##             (Intercept) distance_to_ocean_min_m               max_depth 
    ##               1.0000000               0.7668993               1.0000000 
    ##                 logArea 
    ##               1.0000000

``` r
par(mfrow=c(1,1))
```

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
    ## 3 Stratified vs Reference  1 0.6261283 1.909357 0.2143091   0.126      0.756
    ## 4          Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.048      0.288
    ## 5      Mixed vs Reference  1 0.5979203 1.966332 0.2193017   0.116      0.696
    ## 6      Ocean vs Reference  1 0.5599217 1.628213 0.2456489   0.148      0.888
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
    ## 2 Stratified vs Ocean  1 1.3093089 3.913330 0.2459152   0.002      0.006   *
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.055      0.165

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
    ## 1 Stratified vs Mixed  1 1.3139866 4.158150 0.2289964   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.3841204 4.397984 0.2856207   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.6396174 2.135322 0.1625633   0.009      0.027   .

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
    ## 2 Stratified vs Ocean  1 1.3260066 4.147587 0.2931657   0.008      0.024   .
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
# Surveyed sites 
SD_beta_NMDS <- metaMDS(SD_beta_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
round(SD_beta_NMDS$stress, digits = 2)
SD_beta_rep_NMDS <- metaMDS(SD_beta_dist$Brepl, try = 1000, parallel = 4, trymax = 500, maxit = 500)
SD_beta_ric_NMDS <- metaMDS(SD_beta_dist$Brich, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed and stratified lakes
SD_beta_MS_NMDS <- metaMDS(SD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
SD_beta_OM_NMDS <- metaMDS(SD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
SD_beta_SO_NMDS <- metaMDS(SD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes
SD_beta_S_NMDS <- metaMDS(SD_beta_S_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)


### Environmental
# Surveyed sites 
SD_beta_env_NMDS <- metaMDS(SD_beta_env_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed and stratified lakes
SD_beta_env_MS_NMDS <- metaMDS(SD_beta_env_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
SD_beta_env_OM_NMDS <- metaMDS(SD_beta_env_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
SD_beta_env_SO_NMDS <- metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(SD_beta_env_SO_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
# Mixed lakes
SD_beta_env_M_NMDS <- metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(SD_beta_env_M_dist$Btotal, try = 1000, parallel = 4, :
    ## stress is (nearly) zero: you may have insufficient data

``` r
### Geographic
# Surveyed sites 
SD_beta_geo_NMDS <- metaMDS(SD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed and stratified lakes
SD_beta_geo_MS_NMDS <- metaMDS(SD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
SD_beta_geo_OM_NMDS <- metaMDS(SD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
SD_beta_geo_SO_NMDS <- metaMDS(SD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed lakes
SD_beta_geo_M_NMDS <- metaMDS(SD_beta_geo_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

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
#    see ''?p.adjust' for more options)
# n - optional, number of tests for which to correct; if not given, the number is
#    taken as the number of tests conducted by envfit function (for both vectors and factors).
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
    ## env[surveyed_sites_env, c(34)]  0.997830 -0.065891 0.7388  0.036 *
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
    ## env[surveyed_sites_env, c(34)]  0.997830 -0.065891 0.7388  0.036 *
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
    ## temperature_median -0.48786 -0.87292 0.0586  0.815  
    ## salinity_median     0.99783 -0.06589 0.7388  0.026 *
    ## oxygen_median       0.79309  0.60911 0.5362  0.820  
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
    ## temperature_median -0.48786 -0.87292 0.0586  1.000  
    ## salinity_median     0.99783 -0.06589 0.7388  0.078 .
    ## oxygen_median       0.79309  0.60911 0.5362  1.000  
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
    ## temperature_median -0.56513 -0.82500 0.0778  0.822  
    ## salinity_median     0.99173  0.12835 0.7458  0.034 *
    ## oxygen_median       0.94263  0.33384 0.3552  0.873  
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
    ## temperature_median -0.56513 -0.82500 0.0778  1.000
    ## salinity_median     0.99173  0.12835 0.7458  0.102
    ## oxygen_median       0.94263  0.33384 0.3552  1.000
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
    ## temperature_median -0.10614  0.99435 0.0705  0.702   
    ## salinity_median    -0.78731 -0.61655 0.7842  0.007 **
    ## oxygen_median      -0.96711  0.25435 0.7188  0.052 . 
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
    ## salinity_median    -0.78731 -0.61655 0.7842  0.021 *
    ## oxygen_median      -0.96711  0.25435 0.7188  0.156  
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
    ## temperature_median  0.00008632  1.00000000 0.2469  0.197
    ## salinity_median    -0.00064248 -1.00000000 0.5460  0.325
    ## oxygen_median      -0.00166877  1.00000000 0.7251  0.540
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
    ## temperature_median  0.00008632  1.00000000 0.2469  0.591
    ## salinity_median    -0.00064248 -1.00000000 0.5460  0.975
    ## oxygen_median      -0.00166877  1.00000000 0.7251  1.000
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
    ## temperature_median -0.0022389  1.0000000 0.0356  0.884  
    ## salinity_median     0.0025294 -1.0000000 0.6099  0.072 .
    ## oxygen_median       0.0081883 -0.9999700 0.6148  0.106  
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
    ## temperature_median -0.0022389  1.0000000 0.0356  1.000
    ## salinity_median     0.0025294 -1.0000000 0.6099  0.216
    ## oxygen_median       0.0081883 -0.9999700 0.6148  0.318
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
    ## temperature_median       -0.35881 -0.93341 0.1964  0.636
    ## salinity_median          -0.70036  0.71379 0.5991  0.111
    ## oxygen_median             0.98543 -0.17006 0.5185  0.169
    ## distance_to_ocean_mean_m  0.01326 -0.99991 0.1726  0.608
    ## max_depth                -0.49998 -0.86604 0.0383  0.866
    ## logArea                  -0.26513 -0.96421 0.2144  0.553
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
    ## salinity_median          -0.70036  0.71379 0.5991  0.666
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
    ## distance_to_ocean_min_m -0.99979  0.02055 0.6139  0.138
    ## max_depth               -0.17727 -0.98416 0.1185  0.777
    ## logArea                  0.26565 -0.96407 0.2080  0.129
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
    ## distance_to_ocean_min_m -0.99979  0.02055 0.6139  0.414
    ## max_depth               -0.17727 -0.98416 0.1185  1.000
    ## logArea                  0.26565 -0.96407 0.2080  0.387
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
    ## distance_to_ocean_min_m -0.99979  0.02055 0.6139  0.132
    ## max_depth               -0.17727 -0.98416 0.1185  0.782
    ## logArea                  0.26565 -0.96407 0.2080  0.107
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
    ## distance_to_ocean_min_m -0.99979  0.02055 0.6139  0.396
    ## max_depth               -0.17727 -0.98416 0.1185  1.000
    ## logArea                  0.26565 -0.96407 0.2080  0.321
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
    ## distance_to_ocean_min_m -0.94031  0.34032 0.5171  0.150
    ## max_depth               -0.21483 -0.97665 0.1850  0.655
    ## logArea                 -0.06535  0.99786 0.0118  0.957
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
    ## distance_to_ocean_min_m -0.94031  0.34032 0.5171   0.45
    ## max_depth               -0.21483 -0.97665 0.1850   1.00
    ## logArea                 -0.06535  0.99786 0.0118   1.00
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
    ## distance_to_ocean_min_m  0.86177 -0.50729 0.3198  0.056 .
    ## max_depth               -0.68923 -0.72454 0.5689  0.012 *
    ## logArea                 -0.84908  0.52826 0.2993  0.213  
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
    ## distance_to_ocean_min_m  0.86177 -0.50729 0.3198  0.168  
    ## max_depth               -0.68923 -0.72454 0.5689  0.036 *
    ## logArea                 -0.84908  0.52826 0.2993  0.639  
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
    ## distance_to_ocean_min_m  0.86963  0.49371 0.6793  0.196
    ## max_depth                0.68533  0.72823 0.0510  0.966
    ## logArea                 -0.52042  0.85391 0.2187  0.544
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
    ## distance_to_ocean_min_m  0.86963  0.49371 0.6793  0.588
    ## max_depth                0.68533  0.72823 0.0510  1.000
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
    ## distance_to_ocean_min_m -0.0016833 -1.0000000 0.3212  0.402  
    ## max_depth                0.0046668 -0.9999900 0.7062  0.059 .
    ## logArea                  0.0029950 -1.0000000 0.6575  0.067 .
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
    ## distance_to_ocean_min_m -0.0016833 -1.0000000 0.3212  1.000
    ## max_depth                0.0046668 -0.9999900 0.7062  0.177
    ## logArea                  0.0029950 -1.0000000 0.6575  0.201
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
    ## 0.296 0.328 0.350 0.364 
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
    ## 0.582 0.609 0.633 0.655 
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
    ##       Significance: 0.503 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.587 0.620 0.651 0.679 
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
    ##       Significance: 0.337 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.330 0.369 0.395 0.410 
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
    ## 0.558 0.593 0.630 0.667 
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
    ##       Significance: 0.595 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.472 0.526 0.559 0.597 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_MS_mant_pv <- rbind(SD_beta_env_MS_mant_t$signif, SD_beta_env_MS_mant_s$signif, SD_beta_env_MS_mant_o$signif)
SD_beta_env_MS_mant_pv <- SD_beta_env_MS_mant_pv[,1]
(SD_beta_env_MS_mant_pv <- p.adjust(SD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.033 1.000

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
    ##       Significance: 0.741 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.119 0.158 0.240 0.386 
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
    ##       Significance: 0.058 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.346 0.438 0.510 0.570 
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
    ##       Significance: 0.026 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.360 0.452 0.532 0.638 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_OM_mant_pv <- rbind(SD_beta_env_OM_mant_t$signif, SD_beta_env_OM_mant_s$signif, SD_beta_env_OM_mant_o$signif)
SD_beta_env_OM_mant_pv <- SD_beta_env_OM_mant_pv[,1]
(SD_beta_env_OM_mant_pv <- p.adjust(SD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.174 0.078

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
    ##       Significance: 0.19 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0230 0.0438 0.0645 0.0775 
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
    ##       Significance: 0.013 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.450 0.493 0.531 0.561 
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
    ##       Significance: 0.37 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.705 0.729 0.752 0.768 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_SO_mant_pv <- rbind(SD_beta_env_SO_mant_t$signif, SD_beta_env_SO_mant_s$signif, SD_beta_env_SO_mant_o$signif)
SD_beta_env_SO_mant_pv <- SD_beta_env_SO_mant_pv[,1]
(SD_beta_env_SO_mant_pv <- p.adjust(SD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 0.570 0.039 1.000

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
    ##       Significance: 0.76 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.253 0.321 0.427 0.690 
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
    ##       Significance: 0.158 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.259 0.326 0.365 0.405 
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
    ##       Significance: 0.086 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.237 0.417 0.472 0.579 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_M_mant_pv <- rbind(SD_beta_env_M_mant_t$signif, SD_beta_env_M_mant_s$signif, SD_beta_env_M_mant_o$signif)
SD_beta_env_M_mant_pv <- SD_beta_env_M_mant_pv[,1]
(SD_beta_env_M_mant_pv <- p.adjust(SD_beta_env_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.474 0.258

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
    ## 0.229 0.293 0.334 0.400 
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
    ##       Significance: 0.454 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.551 0.582 0.601 0.628 
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
    ## 0.199 0.225 0.242 0.272 
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
    ## 0.0857 0.0982 0.1090 0.1219 
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
    ##       Significance: 0.387 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.432 0.494 0.542 0.587 
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
    ##       Significance: 0.764 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.209 0.256 0.294 0.322 
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
    ## 0.0725 0.0948 0.1102 0.1317 
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
    ##       Significance: 0.599 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.241 0.281 0.308 0.338 
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
    ##       Significance: 0.058 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.164 0.249 0.290 0.416 
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
    ##       Significance: 0.116 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.207 0.252 0.298 0.333 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_OM_mant_pv <- rbind( SD_beta_geo_OM_mant_dmean$signif, SD_beta_geo_OM_mant_md$signif, SD_beta_geo_OM_mant_la$signif)
SD_beta_geo_OM_mant_pv <- SD_beta_geo_OM_mant_pv[,1]
(SD_beta_geo_OM_mant_pv <- p.adjust(SD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.174 0.348

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
    ##       Significance: 0.469 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.672 0.706 0.733 0.760 
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
    ##       Significance: 0.64 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.133 0.158 0.189 0.228 
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
    ##       Significance: 0.235 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.275 0.294 0.310 0.323 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_SO_mant_pv <- rbind( SD_beta_geo_SO_mant_dmean$signif, SD_beta_geo_SO_mant_md$signif, SD_beta_geo_SO_mant_la$signif)
SD_beta_geo_SO_mant_pv <- SD_beta_geo_SO_mant_pv[,1]
(SD_beta_geo_SO_mant_pv <- p.adjust(SD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 1.000 0.705

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
    ##       Significance: 0.471 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.221 0.381 0.470 0.733 
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
    ## 0.270 0.389 0.481 0.578 
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
    ##       Significance: 0.031 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.215 0.316 0.452 0.546 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_M_mant_pv <- rbind( SD_beta_geo_M_mant_dmean$signif, SD_beta_geo_M_mant_md$signif, SD_beta_geo_M_mant_la$signif)
SD_beta_geo_M_mant_pv <- SD_beta_geo_M_mant_pv[,1]
(SD_beta_geo_M_mant_pv <- p.adjust(SD_beta_geo_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.003 0.093

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
    ##  [1] MASS_7.3-65          ggvenn_0.1.10        pairwiseAdonis_0.4.1
    ##  [4] cluster_2.1.8        BAT_2.9.6            caret_7.0-1         
    ##  [7] ggrepel_0.9.6        ggplot2_3.5.1        picante_1.8.2       
    ## [10] nlme_3.1-167         vegan_2.6-10         lattice_0.22-6      
    ## [13] permute_0.9-7        car_3.1-3            carData_3.0-5       
    ## [16] tidyr_1.3.1          phytools_2.4-4       maps_3.4.2.1        
    ## [19] ape_5.8-1            reshape2_1.4.4       stringr_1.5.1       
    ## [22] dplyr_1.1.4          knitr_1.49          
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
    ##  [61] quadprog_1.5-8          checkmate_2.3.2         generics_0.1.3         
    ##  [64] recipes_1.1.1           gtable_0.3.6            class_7.3-23           
    ##  [67] data.table_1.17.0       hms_1.1.3               foreach_1.5.2          
    ##  [70] pillar_1.10.1           splines_4.4.3           survival_3.8-3         
    ##  [73] SparseM_1.84-2          ks_1.14.3               tidyselect_1.2.1       
    ##  [76] rms_7.0-0               gridExtra_2.3           stats4_4.4.3           
    ##  [79] xfun_0.51               expm_1.0-0              hardhat_1.4.1          
    ##  [82] timeDate_4041.110       proto_1.0.0             stringi_1.8.4          
    ##  [85] yaml_2.3.10             evaluate_1.0.3          codetools_0.2-20       
    ##  [88] tibble_3.2.1            cli_3.6.4               rpart_4.1.24           
    ##  [91] nls2_0.3-4              geometry_0.5.2          systemfonts_1.2.1      
    ##  [94] munsell_0.5.1           Rcpp_1.0.14             globals_0.16.3         
    ##  [97] coda_0.19-4.1           fastcluster_1.2.6       MatrixModels_0.5-3     
    ## [100] gower_1.0.2             prettyunits_1.2.0       mclust_6.1.1           
    ## [103] listenv_0.9.1           phangorn_2.12.1         mvtnorm_1.3-3          
    ## [106] ipred_0.9-15            scales_1.3.0            prodlim_2024.06.25     
    ## [109] e1071_1.7-16            purrr_1.0.4             crayon_1.5.3           
    ## [112] combinat_0.0-8          rlang_1.1.5             multcomp_1.4-28        
    ## [115] fastmatch_1.1-6         mnormt_2.1.1            hypervolume_3.1.5

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
