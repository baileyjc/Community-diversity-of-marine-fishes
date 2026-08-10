Species diversity analyses
================

### Load packages and files

``` r
# Load the knitr package if not already loaded
library(knitr)

# Source the R Markdown file
knit("/Users/bailey/Documents/research/Fish_Biodiversity/src/collection/load_collection_data.Rmd", output = "/Users/bailey/Documents/research/Fish_Biodiversity/src/collection/load_collection_data.md")
```

    ## 
    ## 
    ## processing file: /Users/bailey/Documents/research/Fish_Biodiversity/src/collection/load_collection_data.Rmd

    ##   |                  |          |   0%  |                  |          |   4%                                                                          |                  |.         |   8% [Bringing everything together load modifying files packages]             |                  |.         |  12%                                                                          |                  |..        |  17% [Bringing everything together load in modifying files]                   |                  |..        |  21%                                                                          |                  |..        |  25% [Check species names across files]                                       |                  |...       |  29%                                                                          |                  |...       |  33% [Modify environment data]                                                |                  |....      |  38%                                                                          |                  |....      |  42% [Modify incidence matrices]                                              |                  |.....     |  46%                                                                          |                  |.....     |  50% [Modify phylogeny]                                                       |                  |.....     |  54%                                                                          |                  |......    |  58% [Modify trait data]                                                      |                  |......    |  62%                                                                          |                  |.......   |  67% [Modify Location_type data frames]                                       |                  |.......   |  71%                                                                          |                  |........  |  75% [Modify Location_type trait data]                                        |                  |........  |  79%                                                                          |                  |........  |  83% [Location_type trait data tests]                                         |                  |......... |  88%                                                                          |                  |......... |  92% [Modify site trait data frames]                                          |                  |..........|  96%                                                                          |                  |..........| 100% [Bringing everything together load out modified files and session info]

    ## output file: /Users/bailey/Documents/research/Fish_Biodiversity/src/collection/load_collection_data.md

    ## [1] "/Users/bailey/Documents/research/Fish_Biodiversity/src/collection/load_collection_data.md"

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

    ## Warning: replacing previous import 'e1071::element' by 'ggplot2::element' when
    ## loading 'hypervolume'

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
library(emmeans)
```

    ## Welcome to emmeans.
    ## Caution: You lose important information if you filter this package's results.
    ## See '? untidy'

``` r
library(car)
```

    ## Loading required package: carData

    ## 
    ## Attaching package: 'car'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     recode

``` r
# Site vectors
surveyed_sites <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LCN", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OCM", "OCO", "OLO", "OTM", "RCA", "SLN", "TLN", "ULN")

surveyed_sites_wo_LCN <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OCM", "OCO", "OLO", "OTM", "RCA", "SLN", "TLN", "ULN")

surveyed_sites_wo_OCO <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LCN", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OCM", "OLO", "OTM", "RCA", "SLN", "TLN", "ULN")

surveyed_sites_wo_TLN_HLM <- c("BCM", "CLM", "FLK", "GLK", "HLO", "IBK", "LCN", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OCM", "OCO", "OLO", "OTM", "RCA", "SLN", "ULN")

surveyed_sites_env <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OTM", "RCA", "SLN", "TLN", "ULN")

ocean_mixed_sites <- c("FLK", "HLO", "IBK", "LCN", "LLN", "MLN", "NCN", "NLN", "NLU", "OCM", "OCO", "OLO", "RCA", "ULN")

ocean_mixed_sites_env <- c("FLK", "HLO", "IBK", "LLN", "MLN", "NCN", "NLN", "NLU", "OLO", "RCA", "ULN")

ocean_stratified_sites <- c("BCM", "CLM", "GLK", "HLM", "IBK", "LCN", "NCN", "NLK", "OCM", "OCO", "OTM", "RCA", "SLN", "TLN")

ocean_stratified_sites_env <- c("BCM", "CLM", "GLK", "HLM", "IBK", "NCN", "NLK", "OTM", "RCA", "SLN", "TLN")

mixed_stratified_lakes <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "LLN", "MLN", "NLK", "NLN", "NLU", "OLO", "OTM", "SLN", "TLN", "ULN")

ocean_sites <- c("IBK", "LCN", "NCN", "OCM", "OCO", "RCA")

ocean_sites_env <- c("IBK", "NCN", "RCA")

mixed_lakes <- c("FLK", "HLO", "LLN", "MLN", "NLN", "NLU", "OLO", "ULN")

stratified_lakes <- c("BCM", "CLM", "GLK", "HLM", "NLK", "OTM", "SLN", "TLN")

# Define your custom colors
custom_colors <- c("Reference" = "black", "Ocean" = "#EE6363", "Mixed" = "#87CEFA", "Stratified" = "#6E8B3D")

env_cont <- "#FFB90F"

set.seed(123)
```

### Edit input files

``` r
strat_fish <- stratified_fish[,-c(1,10)]
mix_fish <- mixed_fish[,-c(1,10)]
oc_fish <- ocean_fish[,-c(1,8)]

env <- env[order(env$X),]
presabs_lake <- presabs_lake[order(row.names(presabs_lake)),]
surveyed_sites_lake <- surveyed_sites_lake[order(row.names(surveyed_sites_lake)),]

Location_type_group_ref <- env[,"Location_type"]
Location_type_group <- env[surveyed_sites,"Location_type"]
```

# Species Diversity

## SD alpha diversity

### SD alpha bar plot

``` r
#Plot species richness with a bar plot
SR_env$Location_type <- factor(SR_env$Location_type, levels = c("Ocean", "Mixed", "Stratified"))

(SR_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = reorder(X, row_sum, decreasing = T), y = row_sum, color = Location_type, fill = Location_type, )) + 
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
  labs(x="Site", y="Species Richness", colour = "Location type:", fill = "Location type:"))
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20bar%20plot-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SR_plot.jpg", plot = SR_plot, width = 6, height = 6, units = "in")
```

### SD alpha venn diagram plot

``` r
Location_type_pres <- strat_presabs_lake[c(24:26),]
Location_type_pres <- t(Location_type_pres)
Location_type_pres <- Location_type_pres[which(rowSums(Location_type_pres) > 0),]

# Convert numerical values to logical
Location_type_pres_logical <- as.data.frame(Location_type_pres > 0)

# Check the structure of the transformed data
str(Location_type_pres_logical)
```

    ## 'data.frame':    249 obs. of  3 variables:
    ##  $ Ocean sites     : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ Mixed lakes     : logi  TRUE FALSE TRUE FALSE TRUE FALSE ...
    ##  $ Stratified lakes: logi  FALSE FALSE FALSE FALSE FALSE FALSE ...

``` r
# Create the Venn diagram
(venn_plot <- ggvenn(Location_type_pres_logical,
        show_percentage = F,
        fill_color = c("#EE6363", "#87CEFA", "#6E8B3D"),
        fill_alpha = 0.7,
        stroke_alpha = 0,
        stroke_size = 0.5, 
        set_name_size = 5,
        text_size = 5))
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## ℹ The deprecated feature was likely used in the ggvenn package.
    ##   Please report the issue to the authors.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](SD_analyses_files/figure-gfm/SD%20alpha%20venn%20diagram%20plot-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/venn_plot.jpg", plot = venn_plot, width = 4.25, height = 4, units = "in")
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

# Make Location type columns in each dataframe
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

# Keep only Location type
keep <- c("Ocean sites", "Mixed lakes", "Stratified lakes")
Location_type_pres_sample <- oc_mix_strat_merge[,keep]

# Replace NA with 0
Location_type_pres_sample[is.na(Location_type_pres_sample)] <- 0

## Remove all species not found in any locations
# Identifies which rows are greater than 0
Location_type_pres_sample <- Location_type_pres_sample[which(rowSums(Location_type_pres_sample) > 0),]

Location_type_pres_sample_logical <- as.data.frame(Location_type_pres_sample > 0)

venn_plot <- ggvenn(Location_type_pres_sample_logical,
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
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/venn_plot_s10.jpg", plot = venn_plot, width = 4.25, height = 4, units = "in")
```

### SD alpha outliers for each Location type

``` r
outlier_SD_alpha <- SR_env[surveyed_sites,] %>%
  group_by(Location_type) %>%
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
(outlier_SD_alpha_plot <- ggplot(outlier_SD_alpha, aes(x = Location_type, y = row_sum, fill = Location_type)) +
  geom_violin(alpha = 0.9, draw_quantiles = c(0.25, 0.5, 0.75), aes(fill = Location_type)) +
  geom_jitter(aes(color = is_outlier), width = 0.1, size = 2, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "purple")) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = outlier_SD_alpha, label = outlier_SD_alpha$X, size = 3, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme(text = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(color = "black", size = 12),
    axis.line = element_line(color = "black"),
    axis.title.y = element_text(margin = margin(t = 0)),  # reduce top margin (space from text)
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  scale_y_continuous(expand = c(0,2)) +
  guides(color = "none", fill = "none") + 
  labs(y = "SRic", x = "Location type", color = "Outlier:", tag = "a"))
```

    ## Warning: The `draw_quantiles` argument of `geom_violin()` is deprecated as of ggplot2
    ## 4.0.0.
    ## ℹ Please use the `quantiles.linetype` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](SD_analyses_files/figure-gfm/SD%20alpha%20outliers-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/outlier_SD_alpha.jpg", outlier_SD_alpha_plot, width = 3.25, height = 3.32, units = "in")
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


mod <- aov(temperature_median ~ Location_type, SR_env)
summary(mod)
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## Location_type  2   2.94   1.470   1.092  0.359
    ## Residuals     16  21.54   1.346               
    ## 4 observations deleted due to missingness

``` r
mod <- aov(row_sum ~ oxygen_median + Location_type, SR_env)
car::vif(mod)
```

    ##                   GVIF Df GVIF^(1/(2*Df))
    ## oxygen_median 3.374424  1        1.836961
    ## Location_type 3.374424  2        1.355345

``` r
mod <- aov(salinity_median ~ Location_type, SR_env)
summary(mod)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Location_type  2 147.34   73.67   13.61 0.000353 ***
    ## Residuals     16  86.62    5.41                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 4 observations deleted due to missingness

``` r
mod <- aov(row_sum ~ salinity_median + Location_type, SR_env)
car::vif(mod)
```

    ##                     GVIF Df GVIF^(1/(2*Df))
    ## salinity_median 2.701125  1        1.643510
    ## Location_type   2.701125  2        1.281995

``` r
mod <- aov(oxygen_median ~ Location_type, SR_env)
summary(mod)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Location_type  2 14.440    7.22      19 5.95e-05 ***
    ## Residuals     16  6.081    0.38                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 4 observations deleted due to missingness

``` r
mod <- aov(row_sum ~ oxygen_median + Location_type, SR_env)
car::vif(mod)
```

    ##                   GVIF Df GVIF^(1/(2*Df))
    ## oxygen_median 3.374424  1        1.836961
    ## Location_type 3.374424  2        1.355345

``` r
SR_ss_env <- SR_env[surveyed_sites_env,]
SR_ss_env$Location_type <- factor(SR_ss_env$Location_type, levels = c("Ocean", "Mixed", "Stratified"))

# Run linear discriminant analysis of environmental variables
LDA <- lda(SR_ss_env[,environment], SR_ss_env$Location_type)
spe.class <- predict(LDA)$class
spe.post <- predict(LDA)$posterior
(spe.table <- table(SR_ss_env$Location_type, spe.class))
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

mod <- aov(distance_to_ocean_min_m ~ Location_type, SR_env)
summary(mod)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Location_type  2  80935   40468   16.49 7.05e-05 ***
    ## Residuals     19  46637    2455                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 1 observation deleted due to missingness

``` r
mod <- aov(row_sum ~ distance_to_ocean_min_m + Location_type, SR_env)
car::vif(mod)
```

    ##                             GVIF Df GVIF^(1/(2*Df))
    ## distance_to_ocean_min_m 2.735421  1        1.653911
    ## Location_type           2.735421  2        1.286045

``` r
mod <- aov(max_depth ~ Location_type, SR_env)
summary(mod)
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## Location_type  2    535   267.5   2.235  0.134
    ## Residuals     19   2274   119.7               
    ## 1 observation deleted due to missingness

``` r
mod <- aov(row_sum ~ max_depth + Location_type, SR_env)
car::vif(mod)
```

    ##                   GVIF Df GVIF^(1/(2*Df))
    ## max_depth     1.235315  1        1.111447
    ## Location_type 1.235315  2        1.054252

``` r
mod <- aov(logArea ~ Location_type, SR_env)
summary(mod)
```

    ##               Df Sum Sq Mean Sq F value Pr(>F)
    ## Location_type  2  14.37   7.184   2.341  0.123
    ## Residuals     19  58.32   3.069               
    ## 1 observation deleted due to missingness

``` r
mod <- aov(row_sum ~ logArea + Location_type, SR_env)
car::vif(mod)
```

    ##                   GVIF Df GVIF^(1/(2*Df))
    ## logArea       1.246371  1        1.116410
    ## Location_type 1.246371  2        1.056603

``` r
SR_ss_geo <- SR_env[surveyed_sites,]
SR_ss_geo$Location_type <- factor(SR_ss_geo$Location_type, levels = c("Ocean", "Mixed", "Stratified"))

# Run linear discriminant analysis of geographical variables
LDA <- lda(SR_ss_geo[,geography], SR_ss_geo$Location_type)
spe.class <- predict(LDA)$class
spe.post <- predict(LDA)$posterior
(spe.table <- table(SR_ss_geo$Location_type, spe.class))
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
LDA <- lda(SR_ss_env[,envgeo], SR_ss_env$Location_type)
spe.class <- predict(LDA)$class
spe.post <- predict(LDA)$posterior
(spe.table <- table(SR_ss_env$Location_type, spe.class))
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

### SD alpha & Location type

``` r
### Run negative binomial regression
# All location types
SD_logSRic <- MASS::glm.nb(row_sum ~ Location_type, data = SR_env[surveyed_sites,])
summary(SD_logSRic)$deviance / summary(SD_logSRic)$df.residual
```

    ## [1] 1.20669

``` r
# Check for outliers
outlierTest(SD_logSRic)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## HLM 1.804773           0.087868           NA

``` r
# Plot residuals
plot(SD_logSRic)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-1.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-2.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-3.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-4.png)<!-- -->

``` r
# Summarize the Interaction results
summary(SD_logSRic)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ Location_type, data = SR_env[surveyed_sites, 
    ##     ], init.theta = 3.777842091, link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              3.96714    0.21742  18.246  < 2e-16 ***
    ## Location_typeMixed       0.01952    0.28754   0.068    0.946    
    ## Location_typeStratified -1.95224    0.31149  -6.267 3.67e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(3.7778) family taken to be 1)
    ## 
    ##     Null deviance: 70.758  on 21  degrees of freedom
    ## Residual deviance: 22.927  on 19  degrees of freedom
    ## AIC: 184.73
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  3.78 
    ##           Std. Err.:  1.34 
    ## 
    ##  2 x log-likelihood:  -176.734

``` r
# p-values
logSRic_drop1_results <- drop1(SD_logSRic, test = "Chisq")
logSRic_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ Location_type
    ##               Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>             22.927 182.73                     
    ## Location_type  2   70.758 226.56 47.831 4.108e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
SD_logSRic_p <- emmeans(SD_logSRic, pairwise ~ Location_type, adjust = "bonferroni")
SD_logSRic_p
```

    ## $emmeans
    ##  Location_type emmean    SE  df asymp.LCL asymp.UCL
    ##  Ocean           3.97 0.217 Inf      3.54      4.39
    ##  Mixed           3.99 0.188 Inf      3.62      4.36
    ##  Stratified      2.01 0.223 Inf      1.58      2.45
    ## 
    ## Results are given on the log (not the response) scale. 
    ## Confidence level used: 0.95 
    ## 
    ## $contrasts
    ##  contrast           estimate    SE  df z.ratio p.value
    ##  Ocean - Mixed       -0.0195 0.288 Inf  -0.068  1.0000
    ##  Ocean - Stratified   1.9522 0.311 Inf   6.267  <.0001
    ##  Mixed - Stratified   1.9718 0.292 Inf   6.757  <.0001
    ## 
    ## Results are given on the log (not the response) scale. 
    ## P value adjustment: bonferroni method for 3 tests

``` r
# Ocean vs. Mixed
SD_logSRic_om <- MASS::glm.nb(row_sum ~ Location_type, data = SR_env[ocean_mixed_sites,])
summary(SD_logSRic_om)$deviance / summary(SD_logSRic_om)$df.residual
```

    ## [1] 1.207842

``` r
# Check for outliers
outlierTest(SD_logSRic_om)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OLO -1.829756           0.094495           NA

``` r
# Plot residuals
plot(SD_logSRic_om)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-5.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-6.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-7.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-8.png)<!-- -->

``` r
# Summarize the Interaction results
summary(SD_logSRic_om)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ Location_type, data = SR_env[ocean_mixed_sites, 
    ##     ], init.theta = 4.392367269, link = log)
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)         3.96714    0.20273  19.569   <2e-16 ***
    ## Location_typeMixed  0.01952    0.26810   0.073    0.942    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(4.3924) family taken to be 1)
    ## 
    ##     Null deviance: 14.499  on 13  degrees of freedom
    ## Residual deviance: 14.494  on 12  degrees of freedom
    ## AIC: 135.23
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.39 
    ##           Std. Err.:  1.75 
    ## 
    ##  2 x log-likelihood:  -129.226

``` r
# p-values
logSRic_om_drop1_results <- drop1(SD_logSRic_om, test = "Chisq")
logSRic_om_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ Location_type
    ##               Df Deviance    AIC       LRT Pr(>Chi)
    ## <none>             14.494 133.23                   
    ## Location_type  1   14.499 131.23 0.0052991    0.942

``` r
# Ocean vs. Stratified
SD_logSRic_os <- MASS::glm.nb(row_sum ~ Location_type, data = SR_env[ocean_stratified_sites,])
summary(SD_logSRic_os)$deviance / summary(SD_logSRic_os)$df.residual
```

    ## [1] 1.217959

``` r
# Check for outliers
outlierTest(SD_logSRic_os)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## HLM 1.863906           0.089225           NA

``` r
# Plot residuals
plot(SD_logSRic_os)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-9.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-10.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-11.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-12.png)<!-- -->

``` r
# Summarize the Interaction results
summary(SD_logSRic_os)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ Location_type, data = SR_env[ocean_stratified_sites, 
    ##     ], init.theta = 3.809264115, link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               3.9671     0.2166  18.317  < 2e-16 ***
    ## Location_typeStratified  -1.9522     0.3105  -6.288 3.21e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(3.8093) family taken to be 1)
    ## 
    ##     Null deviance: 53.949  on 13  degrees of freedom
    ## Residual deviance: 14.616  on 12  degrees of freedom
    ## AIC: 107.81
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  3.81 
    ##           Std. Err.:  1.85 
    ## 
    ##  2 x log-likelihood:  -101.81

``` r
# p-values
logSRic_os_drop1_results <- drop1(SD_logSRic_os, test = "Chisq")
logSRic_os_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ Location_type
    ##               Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>             14.616 105.81                     
    ## Location_type  1   53.949 143.14 39.333 3.573e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Mixed vs. Stratified
SD_logSRic_ms <- MASS::glm.nb(row_sum ~ Location_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_ms)$deviance / summary(SD_logSRic_ms)$df.residual
```

    ## [1] 1.17698

``` r
# Check for outliers
outlierTest(SD_logSRic_ms)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## HLM 1.740995            0.10528           NA

``` r
# Plot residuals
plot(SD_logSRic_ms)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-13.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-14.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-15.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-16.png)<!-- -->

``` r
# Summarize the Interaction results
summary(SD_logSRic_ms)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ Location_type, data = SR_env[mixed_stratified_lakes, 
    ##     ], init.theta = 3.217246277, link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               3.9867     0.2029  19.647  < 2e-16 ***
    ## Location_typeStratified  -1.9718     0.3110  -6.341 2.28e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(3.2172) family taken to be 1)
    ## 
    ##     Null deviance: 54.011  on 15  degrees of freedom
    ## Residual deviance: 16.478  on 14  degrees of freedom
    ## AIC: 128.14
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  3.22 
    ##           Std. Err.:  1.34 
    ## 
    ##  2 x log-likelihood:  -122.142

``` r
# p-values
logSRic_ms_drop1_results <- drop1(SD_logSRic_ms, test = "Chisq")
logSRic_ms_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ Location_type
    ##               Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>             16.478 126.14                     
    ## Location_type  1   54.011 161.68 37.533 8.986e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
### Run ANOVA
anova_logSRic_result <- aov(log(row_sum) ~ Location_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(anova_logSRic_result))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(anova_logSRic_result)
    ## W = 0.94518, p-value = 0.2527

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(anova_logSRic_result) ~ SR_env[surveyed_sites,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  0.2139 0.8093
    ##       19

``` r
plot(anova_logSRic_result)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-17.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-18.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-19.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-20.png)<!-- -->

``` r
summary(anova_logSRic_result)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Location_type  2 22.346  11.173   27.41 2.51e-06 ***
    ## Residuals     19  7.745   0.408                     
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
    ## Fit: aov(formula = log(row_sum) ~ Location_type, data = SR_env[surveyed_sites, ])
    ## 
    ## $Location_type
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
Location_type <- tukey_logSRic_result$Location_type

# q-values
(OM_q <- (abs(Location_type[1,1]))/USE)
```

    ## [1] 0.1259392

``` r
(MS_q <- (abs(Location_type[3,1]))/BSE)
```

    ## [1] 9.222747

``` r
(SO_q <- (abs(Location_type[2,1]))/USE)
```

    ## [1] 8.664544

``` r
# Run ANOVA
anova_SRic_result <- aov(row_sum ~ Location_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(anova_SRic_result))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(anova_SRic_result)
    ## W = 0.94431, p-value = 0.2426

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(anova_SRic_result) ~ SR_env[surveyed_sites,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  4.3092 0.02863 *
    ##       19                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(anova_logSRic_result)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-21.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-22.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-23.png)<!-- -->![](SD_analyses_files/figure-gfm/SD%20alpha%20&%20Location%20type-24.png)<!-- -->

``` r
summary(anova_SRic_result)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## Location_type  2  10743    5371   11.01 0.000667 ***
    ## Residuals     19   9268     488                     
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
    ## Fit: aov(formula = row_sum ~ Location_type, data = SR_env[surveyed_sites, ])
    ## 
    ## $Location_type
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
Location_type <- tukey_SRic_result$Location_type

# q-values
(OM_q <- (abs(Location_type[1,1]))/USE)
```

    ## [1] 0.1235068

``` r
(MS_q <- (abs(Location_type[3,1]))/BSE)
```

    ## [1] 5.939085

``` r
(SO_q <- (abs(Location_type[2,1]))/USE)
```

    ## [1] 5.375017

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx
```

### SD alpha env and geo variance partitioning

``` r
envgeo_var <- env[surveyed_sites_env,c(environment, geography)]
surveyed_sites_lake_env <- surveyed_sites_lake[surveyed_sites_env,]
cca_var_all <- cca(surveyed_sites_lake_env ~ ., envgeo_var)
anova(cca_var_all)  # 0.001 *** - it is significant
```

    ## Permutation test for cca under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Model: cca(formula = surveyed_sites_lake_env ~ temperature_median + salinity_median + oxygen_median + distance_to_ocean_min_m + max_depth + logArea, data = envgeo_var)
    ##          Df ChiSquare     F Pr(>F)    
    ## Model     6    1.8862 1.465  0.001 ***
    ## Residual 12    2.5750                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
adjR2_cca <- RsquareAdj(cca_var_all)$adj.r.squared # 0.10 - adjusted R2 explained by the 6 variables is 10%
cca_var_0 <- cca(surveyed_sites_lake_env ~ 1, data = envgeo_var)

# Run forward selection
SD_alpha_os_R2 <- ordiR2step(cca_var_0, scope = formula (cca_var_all), R2scope = adjR2_cca, direction = 'forward', permutations = 999)
```

    ## Step: R2.adj= 0 
    ## Call: surveyed_sites_lake_env ~ 1 
    ##  
    ##                            R2.adjusted
    ## <All variables>            0.104921050
    ## + salinity_median          0.051869597
    ## + oxygen_median            0.051454373
    ## + distance_to_ocean_min_m  0.040597092
    ## + logArea                  0.015347056
    ## + max_depth                0.013469476
    ## <none>                     0.000000000
    ## + temperature_median      -0.000435785
    ## 
    ##                   Df    AIC      F Pr(>F)    
    ## + salinity_median  1 98.302 2.0853  0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.05179378 
    ## Call: surveyed_sites_lake_env ~ salinity_median 
    ##  
    ##                           R2.adjusted
    ## <All variables>            0.10492105
    ## + oxygen_median            0.07823307
    ## + logArea                  0.06912622
    ## + max_depth                0.06899039
    ## + distance_to_ocean_min_m  0.06809448
    ## + temperature_median       0.06068541
    ## <none>                     0.05179378
    ## 
    ##                 Df    AIC      F Pr(>F)   
    ## + oxygen_median  1 98.511 1.5813  0.003 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Step: R2.adj= 0.07859961 
    ## Call: surveyed_sites_lake_env ~ salinity_median + oxygen_median 
    ##  
    ##                           R2.adjusted
    ## <All variables>            0.10492105
    ## + temperature_median       0.08892883
    ## + logArea                  0.08794791
    ## + max_depth                0.08793223
    ## + distance_to_ocean_min_m  0.08254997
    ## <none>                     0.07859961
    ## 
    ##                      Df    AIC      F Pr(>F)
    ## + temperature_median  1 99.007 1.2351  0.135

``` r
SD_alpha_os_R2
```

    ## Call: cca(formula = surveyed_sites_lake_env ~ salinity_median +
    ## oxygen_median, data = envgeo_var)
    ## 
    ## -- Model Summary --
    ## 
    ##               Inertia Proportion Rank
    ## Total          4.4612     1.0000     
    ## Constrained    0.8449     0.1894    2
    ## Unconstrained  3.6164     0.8106   16
    ## 
    ## Inertia is scaled Chi-square
    ## 
    ## -- Note --
    ## 
    ## 18 species (variables) deleted due to missingness.
    ## 
    ## -- Eigenvalues --
    ## 
    ## Eigenvalues for constrained axes:
    ##   CCA1   CCA2 
    ## 0.5483 0.2965 
    ## 
    ## Eigenvalues for unconstrained axes:
    ##    CA1    CA2    CA3    CA4    CA5    CA6    CA7    CA8    CA9   CA10   CA11 
    ## 0.3771 0.3640 0.3300 0.2999 0.2814 0.2722 0.2502 0.2458 0.2207 0.2087 0.1997 
    ##   CA12   CA13   CA14   CA15   CA16 
    ## 0.1887 0.1581 0.1256 0.0601 0.0343

``` r
SD_alpha_os_R2$anova
```

    ##                     R2.adj Df    AIC      F Pr(>F)    
    ## + salinity_median 0.051794  1 98.302 2.0853  0.001 ***
    ## + oxygen_median   0.078600  1 98.511 1.5813  0.003 ** 
    ## <All variables>   0.104921                            
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
SD_alpha_os_R2_adj <- SD_alpha_os_R2
SD_alpha_os_R2_adj$anova$`Pr(>F)` <- p.adjust(SD_alpha_os_R2$anova$`Pr(>F)`, method = 'holm', n = ncol(envgeo_var))
SD_alpha_os_R2_adj$anova
```

    ##                     R2.adj Df    AIC      F Pr(>F)   
    ## + salinity_median 0.051794  1 98.302 2.0853  0.006 **
    ## + oxygen_median   0.078600  1 98.511 1.5813  0.015 * 
    ## <All variables>   0.104921                           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

### SD alpha env and geo plots

``` r
# Temperature
SD_alpha_T_plot <- ggplot(data = SR_env[surveyed_sites_env,], mapping = aes(y = log(row_sum), x = temperature_median, color = Location_type, fill = Location_type)) + 
  geom_point(stat = 'identity',
  size = 4,
  alpha = 1) + 
  geom_smooth(method = "lm", se = FALSE) +
  geom_text_repel(label = SR_env[surveyed_sites_env,1], size = 5, point.padding = 3) +
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
  labs(x="Temperature (ºC)", y="Log-SRic", colour = "Location type:", fill = "Location type:", tag = "a")
(SD_alpha_T_plot <- SD_alpha_T_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_alpha_T_plot.jpg", plot = SD_alpha_T_plot, width = 6.26, height = 6, units = "in")
```

    ## `geom_smooth()` using formula = 'y ~ x'

``` r
# Salinity
SD_alpha_S_plot <- ggplot(data = SR_env[surveyed_sites_env,], mapping = aes(y = log(row_sum), x = salinity_median, color = Location_type, fill = Location_type)) + 
  geom_point(stat = 'identity',
  size = 4,
  alpha = 1) + 
  geom_smooth(method = "lm", se = FALSE) +
  geom_text_repel(label = SR_env[surveyed_sites_env,1], size = 5, point.padding = 3) +
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
  labs(x="Salinity (ppt)", y="Log-SRic", colour = "Location type:", fill = "Location type:", tag = "b")
(SD_alpha_S_plot <- SD_alpha_S_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_alpha_S_plot.jpg", plot = SD_alpha_S_plot, width = 6.26, height = 6, units = "in")
```

    ## `geom_smooth()` using formula = 'y ~ x'

``` r
# Oxygen
SD_alpha_O_plot <- ggplot(data = SR_env[surveyed_sites_env,], mapping = aes(y = log(row_sum), x = oxygen_median, color = Location_type, fill = Location_type)) + 
  geom_point(stat = 'identity',
  size = 4,
  alpha = 1) + 
  geom_smooth(method = "lm", se = FALSE) +
  geom_text_repel(label = SR_env[surveyed_sites_env,1], size = 5, point.padding = 3) +
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
  labs(x="Oxygen (mg/L)", y="Log-SRic", colour = "Location type:", fill = "Location type:", tag = "c")
(SD_alpha_O_plot <- SD_alpha_O_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_alpha_O_plot.jpg", plot = SD_alpha_O_plot, width = 6.26, height = 6, units = "in")
```

    ## `geom_smooth()` using formula = 'y ~ x'

``` r
# Distance from the ocean mean
SD_alpha_D_plot <- ggplot(data = SR_env[surveyed_sites,], mapping = aes(y = log(row_sum), x = distance_to_ocean_min_m, color = Location_type, fill = Location_type)) + 
  geom_point(stat = 'identity',
  size = 4,
  alpha = 1) + 
  geom_smooth(method = "lm", se = FALSE) +
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
  labs(x="Isolation (m)", y="Log-SRic", colour = "Location type:", fill = "Location type:", tag = "a")
(SD_alpha_D_plot <- SD_alpha_D_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20plots-4.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_alpha_D_plot.jpg", plot = SD_alpha_D_plot, width = 6.26, height = 6, units = "in")
```

    ## `geom_smooth()` using formula = 'y ~ x'

``` r
# Max depth
SD_alpha_MD_plot <- ggplot(data = SR_env[surveyed_sites,], mapping = aes(y = log(row_sum), x = max_depth, color = Location_type, fill = Location_type)) + 
  geom_point(stat = 'identity',
  size = 4,
  alpha = 1) + 
  geom_smooth(method = "lm", se = FALSE) +
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
  labs(x="Age (m)", y="Log-SRic", colour = "Location type:", fill = "Location type:", tag = "b")
(SD_alpha_MD_plot <- SD_alpha_MD_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20plots-5.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_alpha_MD_plot.jpg", plot = SD_alpha_MD_plot, width = 6.26, height = 6, units = "in")
```

    ## `geom_smooth()` using formula = 'y ~ x'

``` r
# Log Area
SD_alpha_LA_plot <- ggplot(data = SR_env[surveyed_sites,], mapping = aes(y = log(row_sum), x = logArea, color = Location_type, fill = Location_type)) + 
  geom_point(stat = 'identity',
  size = 4,
  alpha = 1) + 
  geom_smooth(method = "lm", se = FALSE) +
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
  labs(x="Log Area (m"^"2"~")", y="Log-SRic", colour = "Location type:", fill = "Location type:", tag = "c")
(SD_alpha_LA_plot <- SD_alpha_LA_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20plots-6.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_alpha_LA_plot.jpg", plot = SD_alpha_LA_plot, width = 6.26, height = 6, units = "in")
```

    ## `geom_smooth()` using formula = 'y ~ x'

``` r
# 1. Create a linear model object
ocean_model <- lm(log(row_sum) ~ logArea, data = SR_env[ocean_sites,])
# 2. Extract the coefficients
ocean_coefficients <- coef(ocean_model)
# 3. Get the slope value (the second element)
ocean_slope <- ocean_coefficients[2]
# Print the slope
print(ocean_slope)
```

    ##    logArea 
    ## 0.02587099

``` r
# 1. Create a linear model object
mixed_model <- lm(log(row_sum) ~ logArea, data = SR_env[mixed_lakes,])
# 2. Extract the coefficients
mixed_coefficients <- coef(mixed_model)
# 3. Get the slope value (the second element)
mixed_slope <- mixed_coefficients[2]
# Print the slope
print(mixed_slope)
```

    ##   logArea 
    ## 0.3466022

``` r
# 1. Create a linear model object
stratified_model <- lm(log(row_sum) ~ logArea, data = SR_env[stratified_lakes,])
# 2. Extract the coefficients
stratified_coefficients <- coef(stratified_model)
# 3. Get the slope value (the second element)
stratified_slope <- stratified_coefficients[2]
# Print the slope
print(stratified_slope)
```

    ##      logArea 
    ## 0.0006553322

### SD alpha with env and geo linear models & ANCOVAs

``` r
par(mfrow=c(2,2)) 
####### logSRic GLMs
###### Environmental
##### Surveyed sites
#### Temperature
# Interaction
SD_logSRic_am_temp <- MASS::glm.nb(row_sum ~ temperature_median * Location_type, data = SR_env[surveyed_sites_env,])
summary(SD_logSRic_am_temp)$deviance / summary(SD_logSRic_am_temp)$df.residual
```

    ## [1] 1.550575

``` r
# Check for outliers
outlierTest(SD_logSRic_am_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -2.800064           0.016042       0.3048

``` r
# Plot residuals
plot(SD_logSRic_am_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-1.png)<!-- -->

``` r
# Summarize the Interaction results
summary(SD_logSRic_am_temp)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ temperature_median * Location_type, 
    ##     data = SR_env[surveyed_sites_env, ], init.theta = 4.647189241, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                                            Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)                                 -5.6462    12.1272  -0.466    0.642
    ## temperature_median                           0.3182     0.3942   0.807    0.420
    ## Location_typeMixed                          19.8085    15.1331   1.309    0.191
    ## Location_typeStratified                     10.4756    12.8753   0.814    0.416
    ## temperature_median:Location_typeMixed       -0.6536     0.4941  -1.323    0.186
    ## temperature_median:Location_typeStratified  -0.4085     0.4179  -0.978    0.328
    ## 
    ## (Dispersion parameter for Negative Binomial(4.6472) family taken to be 1)
    ## 
    ##     Null deviance: 80.600  on 18  degrees of freedom
    ## Residual deviance: 20.157  on 13  degrees of freedom
    ## AIC: 161.23
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.65 
    ##           Std. Err.:  1.94 
    ## 
    ##  2 x log-likelihood:  -147.234

``` r
# p-values
temp_drop1_results <- drop1(SD_logSRic_am_temp, test = "Chisq")
temp_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ temperature_median * Location_type
    ##                                  Df Deviance    AIC    LRT Pr(>Chi)
    ## <none>                                20.157 159.23                
    ## temperature_median:Location_type  2   22.094 157.17 1.9369   0.3797

``` r
# Additive model
SD_logSRic_am_temp <- MASS::glm.nb(row_sum ~ temperature_median + Location_type, data = SR_env[surveyed_sites_env,])
summary(SD_logSRic_am_temp)$deviance / summary(SD_logSRic_am_temp)$df.residual
```

    ## [1] 1.321835

``` r
# Check for outliers
outlierTest(SD_logSRic_am_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## SLN -2.16466           0.048182      0.91546

``` r
# Plot residuals
plot(SD_logSRic_am_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-2.png)<!-- -->

``` r
### Ocean & Mixed sites
# Interaction
SD_logSRic_OM_am_temp <- MASS::glm.nb(row_sum ~ temperature_median * Location_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_OM_am_temp)$deviance / summary(SD_logSRic_OM_am_temp)$df.residual
```

    ## [1] 1.633824

``` r
# Plot residuals
plot(SD_logSRic_OM_am_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-3.png)<!-- -->

``` r
# p-values
temp_drop1_results <- drop1(SD_logSRic_OM_am_temp, test = "Chisq")
temp_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ temperature_median * Location_type
    ##                                  Df Deviance    AIC    LRT Pr(>Chi)
    ## <none>                                11.437 107.31                
    ## temperature_median:Location_type  1   14.006 107.88 2.5692    0.109

``` r
# Additive model
SD_logSRic_OM_am_temp <- MASS::glm.nb(row_sum ~ temperature_median + Location_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_OM_am_temp)$deviance / summary(SD_logSRic_OM_am_temp)$df.residual
```

    ## [1] 1.424545

``` r
# Check for outliers
outlierTest(SD_logSRic_OM_am_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## NLU -1.877652            0.10252           NA

``` r
# Plot residuals
plot(SD_logSRic_OM_am_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-4.png)<!-- -->

``` r
### Ocean & Stratified sites
# Interaction
SD_logSRic_OS_am_temp <- MASS::glm.nb(row_sum ~ temperature_median * Location_type, data = SR_env[ocean_stratified_sites_env,])
summary(SD_logSRic_OS_am_temp)$deviance / summary(SD_logSRic_OS_am_temp)$df.residual
```

    ## [1] 1.723218

``` r
# Plot residuals
plot(SD_logSRic_OS_am_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-5.png)<!-- -->

``` r
# p-values
temp_drop1_results <- drop1(SD_logSRic_OS_am_temp, test = "Chisq")
temp_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ temperature_median * Location_type
    ##                                  Df Deviance    AIC    LRT Pr(>Chi)
    ## <none>                                12.062 81.616                
    ## temperature_median:Location_type  1   13.092 80.645 1.0294   0.3103

``` r
# Additive model
SD_logSRic_OS_am_temp <- MASS::glm.nb(row_sum ~ temperature_median + Location_type, data = SR_env[ocean_stratified_sites_env,])
summary(SD_logSRic_OS_am_temp)$deviance / summary(SD_logSRic_OS_am_temp)$df.residual
```

    ## [1] 1.460588

``` r
# Check for outliers
outlierTest(SD_logSRic_OS_am_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -2.040207           0.080688      0.88757

``` r
# Plot residuals
plot(SD_logSRic_OS_am_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-6.png)<!-- -->

``` r
### Mixed & Stratified lakes
SD_logSRic_MS_am_temp <- MASS::glm.nb(row_sum ~ temperature_median * Location_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_temp)$deviance / summary(SD_logSRic_MS_am_temp)$df.residual
```

    ## [1] 1.376018

``` r
# Plot residuals
plot(SD_logSRic_MS_am_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-7.png)<!-- -->

``` r
# p-values
temp_drop1_results <- drop1(SD_logSRic_MS_am_temp, test = "Chisq")
temp_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ temperature_median * Location_type
    ##                                  Df Deviance    AIC     LRT Pr(>Chi)
    ## <none>                                16.512 128.76                 
    ## temperature_median:Location_type  1   16.965 127.22 0.45276    0.501

``` r
# Additive model
SD_logSRic_MS_am_temp <- MASS::glm.nb(row_sum ~ temperature_median + Location_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_temp)$deviance / summary(SD_logSRic_MS_am_temp)$df.residual
```

    ## [1] 1.263246

``` r
# Check for outliers
outlierTest(SD_logSRic_MS_am_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## SLN -2.885009           0.013703      0.21924

``` r
# Plot residuals
plot(SD_logSRic_MS_am_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-8.png)<!-- -->

``` r
### Mixed lakes
SD_logSRic_M_am_temp <- MASS::glm.nb(row_sum ~ temperature_median, data = SR_env[mixed_lakes,])

### Stratified lakes
SD_logSRic_S_am_temp <- MASS::glm.nb(row_sum ~ temperature_median, data = SR_env[stratified_lakes,])


#### Salinity
# Interaction
SD_logSRic_am_sal <- MASS::glm.nb(row_sum ~ salinity_median * Location_type, data = SR_env[surveyed_sites_env,])
summary(SD_logSRic_am_sal)$deviance / summary(SD_logSRic_am_sal)$df.residual
```

    ## [1] 1.335344

``` r
# Check for outliers
outlierTest(SD_logSRic_am_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 1.983708           0.070644           NA

``` r
# Plot residuals
plot(SD_logSRic_am_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-9.png)<!-- -->

``` r
# p-values
sal_drop1_results <- drop1(SD_logSRic_am_sal, test = "Chisq")
sal_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ salinity_median * Location_type
    ##                               Df Deviance    AIC    LRT Pr(>Chi)
    ## <none>                             17.360 148.80                
    ## salinity_median:Location_type  2   20.221 147.66 2.8618   0.2391

``` r
# Additive model
SD_logSRic_am_sal <- MASS::glm.nb(row_sum ~ salinity_median + Location_type, data = SR_env[surveyed_sites_env,])
summary(SD_logSRic_am_sal)$deviance / summary(SD_logSRic_am_sal)$df.residual
```

    ## [1] 1.171295

``` r
# Check for outliers
outlierTest(SD_logSRic_am_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.090265            0.05532           NA

``` r
# Plot residuals
plot(SD_logSRic_am_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-10.png)<!-- -->

``` r
### Ocean & Mixed sites
# Interaction
SD_logSRic_OM_am_sal <- MASS::glm.nb(row_sum ~ salinity_median * Location_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_OM_am_sal)$deviance / summary(SD_logSRic_OM_am_sal)$df.residual
```

    ## [1] 1.603739

``` r
# Plot residuals
plot(SD_logSRic_OM_am_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-11.png)<!-- -->

``` r
# p-values
sal_drop1_results <- drop1(SD_logSRic_OM_am_sal, test = "Chisq")
sal_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ salinity_median * Location_type
    ##                               Df Deviance    AIC     LRT Pr(>Chi)
    ## <none>                             11.226 106.26                 
    ## salinity_median:Location_type  1   11.276 104.31 0.05009   0.8229

``` r
# Additive model
SD_logSRic_OM_am_sal <- MASS::glm.nb(row_sum ~ salinity_median + Location_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_OM_am_sal)$deviance / summary(SD_logSRic_OM_am_sal)$df.residual
```

    ## [1] 1.403741

``` r
# Check for outliers
outlierTest(SD_logSRic_OM_am_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 1.861137            0.10504           NA

``` r
# Plot residuals
plot(SD_logSRic_OM_am_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-12.png)<!-- -->

``` r
### Ocean & Stratified sites
# Interaction
SD_logSRic_OS_am_sal <- MASS::glm.nb(row_sum ~ salinity_median * Location_type, data = SR_env[ocean_stratified_sites_env,])
summary(SD_logSRic_OS_am_sal)$deviance / summary(SD_logSRic_OS_am_sal)$df.residual
```

    ## [1] 1.410686

``` r
# Plot residuals
plot(SD_logSRic_OS_am_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-13.png)<!-- -->

``` r
# p-values
sal_drop1_results <- drop1(SD_logSRic_OS_am_sal, test = "Chisq")
sal_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ salinity_median * Location_type
    ##                               Df Deviance    AIC      LRT Pr(>Chi)
    ## <none>                             9.8748 70.493                  
    ## salinity_median:Location_type  1   9.8834 68.502 0.008556   0.9263

``` r
# Additive model
SD_logSRic_OS_am_sal <- MASS::glm.nb(row_sum ~ salinity_median + Location_type, data = SR_env[ocean_stratified_sites_env,])
summary(SD_logSRic_OS_am_sal)$deviance / summary(SD_logSRic_OS_am_sal)$df.residual
```

    ## [1] 1.234366

``` r
# Check for outliers
outlierTest(SD_logSRic_OS_am_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## TLN 1.790903            0.11642           NA

``` r
# Plot residuals
plot(SD_logSRic_OS_am_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-14.png)<!-- -->

``` r
### Mixed & Stratified lakes
SD_logSRic_MS_am_sal <- MASS::glm.nb(row_sum ~ salinity_median * Location_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_sal)$deviance / summary(SD_logSRic_MS_am_sal)$df.residual
```

    ## [1] 1.188636

``` r
# Plot residuals
plot(SD_logSRic_MS_am_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-15.png)<!-- -->

``` r
# p-values
sal_drop1_results <- drop1(SD_logSRic_MS_am_sal, test = "Chisq")
sal_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ salinity_median * Location_type
    ##                               Df Deviance    AIC    LRT Pr(>Chi)
    ## <none>                             14.264 118.94                
    ## salinity_median:Location_type  1   16.767 119.44 2.5039   0.1136

``` r
# Additive model
SD_logSRic_MS_am_sal <- MASS::glm.nb(row_sum ~ salinity_median + Location_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_sal)$deviance / summary(SD_logSRic_MS_am_sal)$df.residual
```

    ## [1] 1.112888

``` r
# Check for outliers
outlierTest(SD_logSRic_MS_am_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 1.984553           0.070539           NA

``` r
# Plot residuals
plot(SD_logSRic_MS_am_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-16.png)<!-- -->

``` r
### Mixed lakes
SD_logSRic_M_am_sal <- MASS::glm.nb(row_sum ~ salinity_median, data = SR_env[mixed_lakes,])

### Stratified lakes
SD_logSRic_S_am_sal <- MASS::glm.nb(row_sum ~ salinity_median, data = SR_env[stratified_lakes,])

## Oxygen
# Interaction
SD_logSRic_am_oxy <- MASS::glm.nb(row_sum ~ oxygen_median * Location_type, data = SR_env[surveyed_sites_env,])
summary(SD_logSRic_am_oxy)$deviance / summary(SD_logSRic_am_oxy)$df.residual
```

    ## [1] 1.414637

``` r
# Check for outliers
outlierTest(SD_logSRic_am_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OLO -2.437004           0.031329      0.59525

``` r
# Plot residuals
plot(SD_logSRic_am_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-17.png)<!-- -->

``` r
# p-values
oxy_drop1_results <- drop1(SD_logSRic_am_oxy, test = "Chisq")
oxy_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ oxygen_median * Location_type
    ##                             Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>                           18.390 148.68                     
    ## oxygen_median:Location_type  2   34.928 161.22 16.537 0.0002564 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get the p-value for Temperature specifically for EACH group
slopes <- emtrends(SD_logSRic_am_oxy, ~ Location_type, var = "oxygen_median")
# view the p-values
test(slopes)
```

    ##  Location_type oxygen_median.trend    SE  df z.ratio p.value
    ##  Ocean                       0.837 0.739 Inf   1.133  0.2572
    ##  Mixed                       0.780 0.298 Inf   2.618  0.0088
    ##  Stratified                 -0.719 0.238 Inf  -3.021  0.0025

``` r
### Ocean & Mixed sites
# Interaction
SD_logSRic_OM_am_oxy <- MASS::glm.nb(row_sum ~ oxygen_median * Location_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_OM_am_oxy)$deviance / summary(SD_logSRic_OM_am_oxy)$df.residual
```

    ## [1] 1.615268

``` r
# Plot residuals
plot(SD_logSRic_OM_am_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-18.png)<!-- -->

``` r
# p-values
oxy_drop1_results <- drop1(SD_logSRic_OM_am_oxy, test = "Chisq")
oxy_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ oxygen_median * Location_type
    ##                             Df Deviance    AIC       LRT Pr(>Chi)
    ## <none>                           11.307 104.19                   
    ## oxygen_median:Location_type  1   11.312 102.20 0.0049625   0.9438

``` r
# Additive model
SD_logSRic_OM_am_oxy <- MASS::glm.nb(row_sum ~ oxygen_median + Location_type, data = SR_env[ocean_mixed_sites_env,])
summary(SD_logSRic_OM_am_oxy)$deviance / summary(SD_logSRic_OM_am_oxy)$df.residual
```

    ## [1] 1.413175

``` r
# Check for outliers
outlierTest(SD_logSRic_OM_am_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OLO -2.848788           0.024731      0.27204

``` r
# Plot residuals
plot(SD_logSRic_OM_am_oxy)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-19.png)<!-- -->

``` r
### Ocean & Stratified sites
# Interaction
SD_logSRic_OS_am_oxy <- MASS::glm.nb(row_sum ~ oxygen_median * Location_type, data = SR_env[ocean_stratified_sites_env,])
```

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

    ## Warning in theta.ml(Y, mu, sum(w), w, limit = control$maxit, trace =
    ## control$trace > : iteration limit reached

``` r
summary(SD_logSRic_OS_am_oxy)$deviance / summary(SD_logSRic_OS_am_oxy)$df.residual
```

    ## [1] 1.903384

``` r
# Plot residuals
plot(SD_logSRic_OS_am_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-20.png)<!-- -->

``` r
# p-values
oxy_drop1_results <- drop1(SD_logSRic_OS_am_oxy, test = "Chisq")
oxy_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ oxygen_median * Location_type
    ##                             Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>                           13.324 68.265                     
    ## oxygen_median:Location_type  1   40.480 93.422 27.157 1.876e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get the p-value for Temperature specifically for EACH group
slopes <- emtrends(SD_logSRic_OS_am_oxy, ~ Location_type, var = "oxygen_median")
# view the p-values
test(slopes)
```

    ##  Location_type oxygen_median.trend    SE  df z.ratio p.value
    ##  Ocean                       0.838 0.244 Inf   3.428  0.0006
    ##  Stratified                 -0.715 0.163 Inf  -4.373  <.0001

``` r
### Mixed & Stratified lakes
SD_logSRic_MS_am_oxy <- MASS::glm.nb(row_sum ~ oxygen_median * Location_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_oxy)$deviance / summary(SD_logSRic_MS_am_oxy)$df.residual
```

    ## [1] 1.249302

``` r
# Plot residuals
plot(SD_logSRic_MS_am_oxy)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-21.png)<!-- -->

``` r
# p-values
oxy_drop1_results <- drop1(SD_logSRic_MS_am_oxy, test = "Chisq")
oxy_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ oxygen_median * Location_type
    ##                             Df Deviance    AIC    LRT Pr(>Chi)    
    ## <none>                           14.992 120.03                    
    ## oxygen_median:Location_type  1   26.713 129.75 11.721 0.000618 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get the p-value for Temperature specifically for EACH group
slopes <- emtrends(SD_logSRic_MS_am_oxy, ~ Location_type, var = "oxygen_median")
# view the p-values
test(slopes)
```

    ##  Location_type oxygen_median.trend    SE  df z.ratio p.value
    ##  Mixed                       0.782 0.336 Inf   2.327  0.0200
    ##  Stratified                 -0.718 0.257 Inf  -2.797  0.0052

``` r
### Mixed lakes
SD_logSRic_M_am_oxy <- MASS::glm.nb(row_sum ~ oxygen_median, data = SR_env[mixed_lakes,])

### Stratified lakes
SD_logSRic_S_am_oxy <- MASS::glm.nb(row_sum ~ oxygen_median, data = SR_env[stratified_lakes,])

# Summarize GLM results and calculate pairwise comparisons
summary(SD_logSRic_am_temp)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ temperature_median + Location_type, 
    ##     data = SR_env[surveyed_sites_env, ], init.theta = 4.035920396, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              7.22177    3.88041   1.861   0.0627 .  
    ## temperature_median      -0.09887    0.12586  -0.786   0.4321    
    ## Location_typeMixed      -0.23937    0.35047  -0.683   0.4946    
    ## Location_typeStratified -2.12343    0.37233  -5.703 1.18e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(4.0359) family taken to be 1)
    ## 
    ##     Null deviance: 71.831  on 18  degrees of freedom
    ## Residual deviance: 19.828  on 15  degrees of freedom
    ## AIC: 159.04
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.04 
    ##           Std. Err.:  1.60 
    ## 
    ##  2 x log-likelihood:  -149.043

``` r
SD_logSRic_amp_temp <- emmeans(SD_logSRic_am_temp, pairwise ~ Location_type, adjust = "bonferroni")
SD_logSRic_amp_temp$contrasts
```

    ##  contrast           estimate    SE  df z.ratio p.value
    ##  Ocean - Mixed         0.239 0.350 Inf   0.683  1.0000
    ##  Ocean - Stratified    2.123 0.372 Inf   5.703  <.0001
    ##  Mixed - Stratified    1.884 0.301 Inf   6.251  <.0001
    ## 
    ## Results are given on the log (not the response) scale. 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(SD_logSRic_am_sal)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ salinity_median + Location_type, 
    ##     data = SR_env[surveyed_sites_env, ], init.theta = 6.679166694, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)             -2.35571    2.00440  -1.175  0.23989   
    ## salinity_median          0.19464    0.05943   3.275  0.00106 **
    ## Location_typeMixed      -0.08347    0.27798  -0.300  0.76397   
    ## Location_typeStratified -1.17562    0.42800  -2.747  0.00602 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(6.6792) family taken to be 1)
    ## 
    ##     Null deviance: 107.031  on 18  degrees of freedom
    ## Residual deviance:  17.569  on 15  degrees of freedom
    ## AIC: 149.45
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  6.68 
    ##           Std. Err.:  2.71 
    ## 
    ##  2 x log-likelihood:  -139.454

``` r
SD_logSRic_amp_sal <- emmeans(SD_logSRic_am_sal, pairwise ~ Location_type, adjust = "bonferroni")
SD_logSRic_amp_sal$contrasts
```

    ##  contrast           estimate    SE  df z.ratio p.value
    ##  Ocean - Mixed        0.0835 0.278 Inf   0.300  1.0000
    ##  Ocean - Stratified   1.1756 0.428 Inf   2.747  0.0181
    ##  Mixed - Stratified   1.0921 0.361 Inf   3.024  0.0075
    ## 
    ## Results are given on the log (not the response) scale. 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(SD_logSRic_am_oxy)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ oxygen_median * Location_type, 
    ##     data = SR_env[surveyed_sites_env, ], init.theta = 8.762466591, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                                       Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)                           -0.42181    4.03351  -0.105    0.917  
    ## oxygen_median                          0.83681    0.73863   1.133    0.257  
    ## Location_typeMixed                     0.72672    4.26943   0.170    0.865  
    ## Location_typeStratified                4.58007    4.09928   1.117    0.264  
    ## oxygen_median:Location_typeMixed      -0.05683    0.79644  -0.071    0.943  
    ## oxygen_median:Location_typeStratified -1.55569    0.77601  -2.005    0.045 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(8.7625) family taken to be 1)
    ## 
    ##     Null deviance: 130.60  on 18  degrees of freedom
    ## Residual deviance:  18.39  on 13  degrees of freedom
    ## AIC: 150.68
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  8.76 
    ##           Std. Err.:  3.95 
    ## 
    ##  2 x log-likelihood:  -136.683

``` r
SD_logSRic_amp_oxy <- emmeans(SD_logSRic_am_oxy, pairwise ~ Location_type, adjust = "bonferroni")
```

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
SD_logSRic_amp_oxy$contrasts
```

    ##  contrast           estimate    SE  df z.ratio p.value
    ##  Ocean - Mixed         -0.49 0.997 Inf  -0.492  1.0000
    ##  Ocean - Stratified     1.90 1.030 Inf   1.835  0.1993
    ##  Mixed - Stratified     2.39 0.392 Inf   6.089  <.0001
    ## 
    ## Results are given on the log (not the response) scale. 
    ## P value adjustment: bonferroni method for 3 tests

``` r
# p-values
temp_drop1_results <- drop1(SD_logSRic_am_temp, test = "Chisq")
temp_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ temperature_median + Location_type
    ##                    Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>                  19.828 157.04                     
    ## temperature_median  1   20.338 155.55  0.511    0.4748    
    ## Location_type       2   63.789 197.00 43.962 2.844e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
temp_p_values <- temp_drop1_results[["Pr(>Chi)"]]
sal_drop1_results <- drop1(SD_logSRic_am_sal, test = "Chisq")
sal_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ salinity_median + Location_type
    ##                 Df Deviance    AIC     LRT  Pr(>Chi)    
    ## <none>               17.569 147.45                      
    ## salinity_median  1   29.533 157.42 11.9635 0.0005425 ***
    ## Location_type    2   26.030 151.91  8.4611 0.0145444 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
sal_p_values <- sal_drop1_results[["Pr(>Chi)"]]
oxy_drop1_results <- drop1(SD_logSRic_am_oxy, test = "Chisq")
oxy_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ oxygen_median * Location_type
    ##                             Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>                           18.390 148.68                     
    ## oxygen_median:Location_type  2   34.928 161.22 16.537 0.0002564 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
oxy_p_values <- oxy_drop1_results[["Pr(>Chi)"]]
p_values <- c(temp_p_values[2:3], sal_p_values[2:3], oxy_p_values[2])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.000000e+00 1.421829e-09 2.712656e-03 7.272205e-02 1.282196e-03

``` r
# Summarize OM lm results
temp_drop1_OM_results <- drop1(SD_logSRic_OM_am_temp, test = "Chisq")
temp_drop1_OM_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ temperature_median + Location_type
    ##                    Df Deviance    AIC     LRT Pr(>Chi)
    ## <none>                  11.396 107.61                 
    ## temperature_median  1   11.671 105.89 0.27490   0.6001
    ## Location_type       1   11.931 106.15 0.53433   0.4648

``` r
summary(SD_logSRic_OM_am_temp)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ temperature_median + Location_type, 
    ##     data = SR_env[ocean_mixed_sites_env, ], init.theta = 5.015063015, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)          7.6203     7.0644   1.079    0.281
    ## temperature_median  -0.1118     0.2297  -0.487    0.627
    ## Location_typeMixed  -0.2471     0.3239  -0.763    0.446
    ## 
    ## (Dispersion parameter for Negative Binomial(5.0151) family taken to be 1)
    ## 
    ##     Null deviance: 11.997  on 10  degrees of freedom
    ## Residual deviance: 11.396  on  8  degrees of freedom
    ## AIC: 109.61
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  5.02 
    ##           Std. Err.:  2.30 
    ## 
    ##  2 x log-likelihood:  -101.613

``` r
sal_drop1_OM_results <- drop1(SD_logSRic_OM_am_sal, test = "Chisq")
sal_drop1_OM_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ salinity_median + Location_type
    ##                 Df Deviance    AIC    LRT Pr(>Chi)  
    ## <none>               11.230 104.31                  
    ## salinity_median  1   15.383 106.47 4.1534  0.04155 *
    ## Location_type    1   11.711 102.79 0.4814  0.48781  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OM_am_sal)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ salinity_median + Location_type, 
    ##     data = SR_env[ocean_mixed_sites_env, ], init.theta = 6.852629852, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)        -23.5262    12.1814  -1.931   0.0534 .
    ## salinity_median      0.8267     0.3636   2.274   0.0230 *
    ## Location_typeMixed   0.2442     0.3381   0.722   0.4702  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(6.8526) family taken to be 1)
    ## 
    ##     Null deviance: 15.815  on 10  degrees of freedom
    ## Residual deviance: 11.230  on  8  degrees of freedom
    ## AIC: 106.31
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  6.85 
    ##           Std. Err.:  3.25 
    ## 
    ##  2 x log-likelihood:  -98.313

``` r
oxy_drop1_OM_results <- drop1(SD_logSRic_OM_am_oxy, test = "Chisq")
oxy_drop1_OM_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ oxygen_median + Location_type
    ##               Df Deviance    AIC    LRT Pr(>Chi)   
    ## <none>             11.305 102.20                   
    ## oxygen_median  1   18.681 107.57 7.3754 0.006612 **
    ## Location_type  1   12.904 101.80 1.5982 0.206160   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OM_am_oxy)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ oxygen_median + Location_type, 
    ##     data = SR_env[ocean_mixed_sites_env, ], init.theta = 8.5999958, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)         -0.1581     1.5333  -0.103  0.91785   
    ## oxygen_median        0.7885     0.2785   2.831  0.00464 **
    ## Location_typeMixed   0.4233     0.3285   1.289  0.19748   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(8.6) family taken to be 1)
    ## 
    ##     Null deviance: 19.209  on 10  degrees of freedom
    ## Residual deviance: 11.305  on  8  degrees of freedom
    ## AIC: 104.2
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  8.60 
    ##           Std. Err.:  4.28 
    ## 
    ##  2 x log-likelihood:  -96.198

``` r
# p-values
temp_p_values <- temp_drop1_OM_results[["Pr(>Chi)"]]
sal_p_values <- sal_drop1_OM_results[["Pr(>Chi)"]]
oxy_p_values <- oxy_drop1_OM_results[["Pr(>Chi)"]]
p_values <- c(temp_p_values[2:3], sal_p_values[2:3], oxy_p_values[2:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.00000000 1.00000000 0.24931263 1.00000000 0.03967416 1.00000000

``` r
# Summarize OS lm results
temp_drop1_OS_results <- drop1(SD_logSRic_OS_am_temp, test = "Chisq")
temp_drop1_OS_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ temperature_median + Location_type
    ##                    Df Deviance     AIC    LRT  Pr(>Chi)    
    ## <none>                  11.685  80.568                     
    ## temperature_median  1   11.722  78.605  0.037    0.8466    
    ## Location_type       1   51.134 118.017 39.449 3.367e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OS_am_temp)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ temperature_median + Location_type, 
    ##     data = SR_env[ocean_stratified_sites_env, ], init.theta = 4.202979784, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              5.08309    4.21646   1.206    0.228    
    ## temperature_median      -0.02973    0.13683  -0.217    0.828    
    ## Location_typeStratified -2.14143    0.36806  -5.818 5.95e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(4.203) family taken to be 1)
    ## 
    ##     Null deviance: 52.244  on 10  degrees of freedom
    ## Residual deviance: 11.685  on  8  degrees of freedom
    ## AIC: 82.568
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.20 
    ##           Std. Err.:  2.63 
    ## 
    ##  2 x log-likelihood:  -74.568

``` r
sal_drop1_OS_results <- drop1(SD_logSRic_OS_am_sal, test = "Chisq")
sal_drop1_OS_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ salinity_median + Location_type
    ##                 Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>               9.8749 68.502                     
    ## salinity_median  1  27.3421 83.969 17.467 2.923e-05 ***
    ## Location_type    1  28.0551 84.682 18.180 2.010e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OS_am_sal)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ salinity_median + Location_type, 
    ##     data = SR_env[ocean_stratified_sites_env, ], init.theta = 25.63029503, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              -2.0421     1.5499  -1.318    0.188    
    ## salinity_median           0.1853     0.0461   4.019 5.84e-05 ***
    ## Location_typeStratified  -1.2226     0.2898  -4.219 2.46e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(25.6303) family taken to be 1)
    ## 
    ##     Null deviance: 161.2987  on 10  degrees of freedom
    ## Residual deviance:   9.8749  on  8  degrees of freedom
    ## AIC: 70.502
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  25.6 
    ##           Std. Err.:  24.9 
    ## 
    ##  2 x log-likelihood:  -62.502

``` r
oxy_drop1_OS_results <- drop1(SD_logSRic_OS_am_oxy, test = "Chisq")
oxy_drop1_OS_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ oxygen_median * Location_type
    ##                             Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>                           13.324 68.265                     
    ## oxygen_median:Location_type  1   40.480 93.422 27.157 1.876e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OS_am_oxy)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ oxygen_median * Location_type, 
    ##     data = SR_env[ocean_stratified_sites_env, ], init.theta = 239937.184, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                                       Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                            -0.4256     1.3496  -0.315 0.752477    
    ## oxygen_median                           0.8375     0.2443   3.428 0.000609 ***
    ## Location_typeStratified                 4.5714     1.4290   3.199 0.001379 ** 
    ## oxygen_median:Location_typeStratified  -1.5520     0.2939  -5.280 1.29e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(239937.2) family taken to be 1)
    ## 
    ##     Null deviance: 306.061  on 10  degrees of freedom
    ## Residual deviance:  13.324  on  7  degrees of freedom
    ## AIC: 70.265
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  239937 
    ##           Std. Err.:  5501508 
    ## Warning while fitting theta: iteration limit reached 
    ## 
    ##  2 x log-likelihood:  -60.265

``` r
# p-values
temp_p_values <- temp_drop1_OS_results[["Pr(>Chi)"]]
sal_p_values <- sal_drop1_OS_results[["Pr(>Chi)"]]
oxy_p_values <- oxy_drop1_OS_results[["Pr(>Chi)"]]
p_values <- c(temp_p_values[2:3], sal_p_values[2:3], oxy_p_values[2])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.000000e+00 1.683290e-09 1.461587e-04 1.004786e-04 9.380782e-07

``` r
# Summarize MS lm results
temp_drop1_MS_results <- drop1(SD_logSRic_MS_am_temp, test = "Chisq")
temp_drop1_MS_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ temperature_median + Location_type
    ##                    Df Deviance    AIC     LRT  Pr(>Chi)    
    ## <none>                  16.422 127.21                      
    ## temperature_median  1   17.383 126.17  0.9607     0.327    
    ## Location_type       1   43.130 151.91 26.7077 2.367e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_MS_am_temp)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ temperature_median + Location_type, 
    ##     data = SR_env[mixed_stratified_lakes, ], init.theta = 3.441261366, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               8.6994     4.2272   2.058   0.0396 *  
    ## temperature_median       -0.1555     0.1389  -1.119   0.2629    
    ## Location_typeStratified  -1.8325     0.3209  -5.710 1.13e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(3.4413) family taken to be 1)
    ## 
    ##     Null deviance: 57.156  on 15  degrees of freedom
    ## Residual deviance: 16.422  on 13  degrees of freedom
    ## AIC: 129.21
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  3.44 
    ##           Std. Err.:  1.46 
    ## 
    ##  2 x log-likelihood:  -121.208

``` r
sal_drop1_MS_results <- drop1(SD_logSRic_MS_am_sal, test = "Chisq")
sal_drop1_MS_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ salinity_median + Location_type
    ##                 Df Deviance    AIC     LRT Pr(>Chi)   
    ## <none>               14.467 119.25                    
    ## salinity_median  1   25.224 128.01 10.7562 0.001039 **
    ## Location_type    1   21.684 124.47  7.2164 0.007224 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_MS_am_sal)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ salinity_median + Location_type, 
    ##     data = SR_env[mixed_stratified_lakes, ], init.theta = 5.627390656, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)             -2.38800    2.05739  -1.161  0.24577   
    ## salinity_median          0.19308    0.06229   3.100  0.00194 **
    ## Location_typeStratified -1.09799    0.38443  -2.856  0.00429 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(5.6274) family taken to be 1)
    ## 
    ##     Null deviance: 84.907  on 15  degrees of freedom
    ## Residual deviance: 14.468  on 13  degrees of freedom
    ## AIC: 121.25
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  5.63 
    ##           Std. Err.:  2.45 
    ## 
    ##  2 x log-likelihood:  -113.252

``` r
oxy_drop1_MS_results <- drop1(SD_logSRic_MS_am_oxy, test = "Chisq")
oxy_drop1_MS_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ oxygen_median * Location_type
    ##                             Df Deviance    AIC    LRT Pr(>Chi)    
    ## <none>                           14.992 120.03                    
    ## oxygen_median:Location_type  1   26.713 129.75 11.721 0.000618 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_MS_am_oxy)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ oxygen_median * Location_type, 
    ##     data = SR_env[mixed_stratified_lakes, ], init.theta = 6.578880876, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                                       Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                             0.2973     1.5765   0.189  0.85042    
    ## oxygen_median                           0.7816     0.3359   2.327  0.01998 *  
    ## Location_typeStratified                 3.8590     1.7662   2.185  0.02890 *  
    ## oxygen_median:Location_typeStratified  -1.4998     0.4229  -3.547  0.00039 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(6.5789) family taken to be 1)
    ## 
    ##     Null deviance: 95.580  on 15  degrees of freedom
    ## Residual deviance: 14.992  on 12  degrees of freedom
    ## AIC: 122.03
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  6.58 
    ##           Std. Err.:  3.07 
    ## 
    ##  2 x log-likelihood:  -112.025

``` r
# p-values
temp_p_values <- temp_drop1_MS_results[["Pr(>Chi)"]]
sal_p_values <- sal_drop1_MS_results[["Pr(>Chi)"]]
oxy_p_values <- oxy_drop1_MS_results[["Pr(>Chi)"]]
p_values <- c(temp_p_values[2:3], sal_p_values[2:3], oxy_p_values[2])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.000000e+00 1.183385e-06 5.196396e-03 3.611930e-02 3.089819e-03

``` r
# Summarize M lm results
temp_drop1_M_results <- drop1(SD_logSRic_M_am_temp, test = "Chisq")
temp_drop1_M_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ temperature_median
    ##                    Df Deviance    AIC    LRT Pr(>Chi)
    ## <none>                  8.3080 77.599                
    ## temperature_median  1   9.7452 77.036 1.4372   0.2306

``` r
summary(SD_logSRic_M_am_temp)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ temperature_median, data = SR_env[mixed_lakes, 
    ##     ], init.theta = 4.444037406, link = log)
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)         14.1571     9.2409   1.532    0.126
    ## temperature_median  -0.3352     0.3040  -1.103    0.270
    ## 
    ## (Dispersion parameter for Negative Binomial(4.444) family taken to be 1)
    ## 
    ##     Null deviance: 9.7452  on 7  degrees of freedom
    ## Residual deviance: 8.3080  on 6  degrees of freedom
    ## AIC: 79.599
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.44 
    ##           Std. Err.:  2.37 
    ## 
    ##  2 x log-likelihood:  -73.599

``` r
sal_drop1_M_results <- drop1(SD_logSRic_M_am_sal, test = "Chisq")
sal_drop1_M_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ salinity_median
    ##                 Df Deviance    AIC    LRT Pr(>Chi)  
    ## <none>               8.1651 76.056                  
    ## salinity_median  1  11.5518 77.442 3.3867  0.06573 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_M_am_sal)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ salinity_median, data = SR_env[mixed_lakes, 
    ##     ], init.theta = 5.3688549, link = log)
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)     -24.0148    13.5432  -1.773   0.0762 .
    ## salinity_median   0.8490     0.4111   2.065   0.0389 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(5.3689) family taken to be 1)
    ## 
    ##     Null deviance: 11.5518  on 7  degrees of freedom
    ## Residual deviance:  8.1651  on 6  degrees of freedom
    ## AIC: 78.056
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  5.37 
    ##           Std. Err.:  2.90 
    ## 
    ##  2 x log-likelihood:  -72.056

``` r
oxy_drop1_M_results <- drop1(SD_logSRic_M_am_oxy, test = "Chisq")
oxy_drop1_M_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ oxygen_median
    ##               Df Deviance    AIC    LRT Pr(>Chi)  
    ## <none>             8.1869 75.378                  
    ## oxygen_median  1  12.5817 77.773 4.3947  0.03605 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_M_am_oxy)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ oxygen_median, data = SR_env[mixed_lakes, 
    ##     ], init.theta = 5.91189159, link = log)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)  
    ## (Intercept)     0.2949     1.6504   0.179   0.8582  
    ## oxygen_median   0.7821     0.3518   2.223   0.0262 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(5.9119) family taken to be 1)
    ## 
    ##     Null deviance: 12.5817  on 7  degrees of freedom
    ## Residual deviance:  8.1869  on 6  degrees of freedom
    ## AIC: 77.378
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  5.91 
    ##           Std. Err.:  3.26 
    ## 
    ##  2 x log-likelihood:  -71.378

``` r
# p-values
temp_p_values <- temp_drop1_M_results[["Pr(>Chi)"]]
sal_p_values <- sal_drop1_M_results[["Pr(>Chi)"]]
oxy_p_values <- oxy_drop1_M_results[["Pr(>Chi)"]]
p_values <- c(temp_p_values[2], sal_p_values[2], oxy_p_values[2])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.6917762 0.1971758 0.1081495

``` r
# Summarize S lm results
temp_drop1_S_results <- drop1(SD_logSRic_S_am_temp, test = "Chisq")
temp_drop1_S_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ temperature_median
    ##                    Df Deviance    AIC     LRT Pr(>Chi)
    ## <none>                  7.8306 50.793                 
    ## temperature_median  1   8.0576 49.020 0.22693   0.6338

``` r
summary(SD_logSRic_S_am_temp)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ temperature_median, data = SR_env[stratified_lakes, 
    ##     ], init.theta = 2.656583359, link = log)
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)         5.02556    5.26630   0.954    0.340
    ## temperature_median -0.09655    0.16857  -0.573    0.567
    ## 
    ## (Dispersion parameter for Negative Binomial(2.6566) family taken to be 1)
    ## 
    ##     Null deviance: 8.0576  on 7  degrees of freedom
    ## Residual deviance: 7.8306  on 6  degrees of freedom
    ## AIC: 52.793
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.66 
    ##           Std. Err.:  1.69 
    ## 
    ##  2 x log-likelihood:  -46.793

``` r
sal_drop1_S_results <- drop1(SD_logSRic_S_am_sal, test = "Chisq")
sal_drop1_S_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ salinity_median
    ##                 Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>               7.2449 41.682                     
    ## salinity_median  1  26.0806 58.518 18.836 1.425e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_S_am_sal)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ salinity_median, data = SR_env[stratified_lakes, 
    ##     ], init.theta = 36.79033432, link = log)
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     -3.27315    1.29696  -2.524   0.0116 *  
    ## salinity_median  0.18558    0.04423   4.195 2.72e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(36.7903) family taken to be 1)
    ## 
    ##     Null deviance: 26.0806  on 7  degrees of freedom
    ## Residual deviance:  7.2449  on 6  degrees of freedom
    ## AIC: 43.682
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  37 
    ##           Std. Err.:  107 
    ## 
    ##  2 x log-likelihood:  -37.682

``` r
oxy_drop1_S_results <- drop1(SD_logSRic_S_am_oxy, test = "Chisq")
oxy_drop1_S_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ oxygen_median
    ##               Df Deviance    AIC   LRT Pr(>Chi)   
    ## <none>             7.1886 44.482                  
    ## oxygen_median  1  17.2783 52.572 10.09 0.001491 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_S_am_oxy)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ oxygen_median, data = SR_env[stratified_lakes, 
    ##     ], init.theta = 9.56184012, link = log)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)     4.1586     0.7139   5.825  5.7e-09 ***
    ## oxygen_median  -0.7190     0.2329  -3.087  0.00202 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(9.5618) family taken to be 1)
    ## 
    ##     Null deviance: 17.2783  on 7  degrees of freedom
    ## Residual deviance:  7.1886  on 6  degrees of freedom
    ## AIC: 46.482
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  9.6 
    ##           Std. Err.:  10.7 
    ## 
    ##  2 x log-likelihood:  -40.482

``` r
# p-values
temp_p_values <- temp_drop1_S_results[["Pr(>Chi)"]]
sal_p_values <- sal_drop1_S_results[["Pr(>Chi)"]]
oxy_p_values <- oxy_drop1_S_results[["Pr(>Chi)"]]
p_values <- c(temp_p_values[2], sal_p_values[2], oxy_p_values[2])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.000000e+00 4.274267e-05 4.473182e-03

``` r
###### Geographical
##### Surveyed sites
#### Distance
# Interaction
SD_logSRic_am_dist <- MASS::glm.nb(row_sum ~ distance_to_ocean_min_m * Location_type, data = SR_env[surveyed_sites,])
summary(SD_logSRic_am_dist)$deviance / summary(SD_logSRic_am_dist)$df.residual
```

    ## [1] 1.345263

``` r
# Check for outliers
outlierTest(SD_logSRic_am_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.060057           0.057183           NA

``` r
# Plot residuals
plot(SD_logSRic_am_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   19

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-22.png)<!-- -->

``` r
# Summarize the Interaction results
summary(SD_logSRic_am_dist)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ distance_to_ocean_min_m * Location_type, 
    ##     data = SR_env[surveyed_sites, ], init.theta = 4.713439441, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                                                  Estimate Std. Error z value
    ## (Intercept)                                      3.962716   0.215021  18.429
    ## distance_to_ocean_min_m                          0.001142   0.022879   0.050
    ## Location_typeMixed                               0.188360   0.482286   0.391
    ## Location_typeStratified                         -0.969699   0.538409  -1.801
    ## distance_to_ocean_min_m:Location_typeMixed      -0.003582   0.023609  -0.152
    ## distance_to_ocean_min_m:Location_typeStratified -0.008268   0.023097  -0.358
    ##                                                 Pr(>|z|)    
    ## (Intercept)                                       <2e-16 ***
    ## distance_to_ocean_min_m                           0.9602    
    ## Location_typeMixed                                0.6961    
    ## Location_typeStratified                           0.0717 .  
    ## distance_to_ocean_min_m:Location_typeMixed        0.8794    
    ## distance_to_ocean_min_m:Location_typeStratified   0.7203    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(4.7134) family taken to be 1)
    ## 
    ##     Null deviance: 84.898  on 21  degrees of freedom
    ## Residual deviance: 21.524  on 16  degrees of freedom
    ## AIC: 185.34
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.71 
    ##           Std. Err.:  1.67 
    ## 
    ##  2 x log-likelihood:  -171.34

``` r
# p-values
dist_drop1_results <- drop1(SD_logSRic_am_dist, test = "Chisq")
dist_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ distance_to_ocean_min_m * Location_type
    ##                                       Df Deviance    AIC     LRT Pr(>Chi)
    ## <none>                                     21.524 183.34                 
    ## distance_to_ocean_min_m:Location_type  2   22.101 179.92 0.57658   0.7495

``` r
# Additive model
SD_logSRic_am_dist <- MASS::glm.nb(row_sum ~ distance_to_ocean_min_m + Location_type, data = SR_env[surveyed_sites,])
summary(SD_logSRic_am_dist)$deviance / summary(SD_logSRic_am_dist)$df.residual
```

    ## [1] 1.195106

``` r
# Check for outliers
outlierTest(SD_logSRic_am_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.345346           0.031402      0.69085

``` r
# Plot residuals
plot(SD_logSRic_am_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-23.png)<!-- -->

``` r
### Ocean & Mixed sites
# Interaction
SD_logSRic_OM_am_dist <- MASS::glm.nb(row_sum ~ distance_to_ocean_min_m * Location_type, data = SR_env[ocean_mixed_sites,])
summary(SD_logSRic_OM_am_dist)$deviance / summary(SD_logSRic_OM_am_dist)$df.residual
```

    ## [1] 1.448864

``` r
# Plot residuals
plot(SD_logSRic_OM_am_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   13

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-24.png)<!-- -->

``` r
# p-values
dist_drop1_results <- drop1(SD_logSRic_OM_am_dist, test = "Chisq")
dist_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ distance_to_ocean_min_m * Location_type
    ##                                       Df Deviance    AIC      LRT Pr(>Chi)
    ## <none>                                     14.489 137.07                  
    ## distance_to_ocean_min_m:Location_type  1   14.511 135.10 0.021964   0.8822

``` r
# Additive model
SD_logSRic_OM_am_dist <- MASS::glm.nb(row_sum ~ distance_to_ocean_min_m + Location_type, data = SR_env[ocean_mixed_sites,])
summary(SD_logSRic_OM_am_dist)$deviance / summary(SD_logSRic_OM_am_dist)$df.residual
```

    ## [1] 1.31723

``` r
# Check for outliers
outlierTest(SD_logSRic_OM_am_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 2.049628            0.06755       0.9457

``` r
# Plot residuals
plot(SD_logSRic_OM_am_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-25.png)<!-- -->

``` r
### Ocean & Stratified sites
# Interaction
SD_logSRic_OS_am_dist <- MASS::glm.nb(row_sum ~ distance_to_ocean_min_m * Location_type, data = SR_env[ocean_stratified_sites,])
summary(SD_logSRic_OS_am_dist)$deviance / summary(SD_logSRic_OS_am_dist)$df.residual
```

    ## [1] 1.339913

``` r
# Plot residuals
plot(SD_logSRic_OS_am_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   12

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-26.png)<!-- -->

``` r
# p-values
dist_drop1_results <- drop1(SD_logSRic_OS_am_dist, test = "Chisq")
dist_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ distance_to_ocean_min_m * Location_type
    ##                                       Df Deviance    AIC     LRT Pr(>Chi)
    ## <none>                                     13.399 104.18                 
    ## distance_to_ocean_min_m:Location_type  1   13.563 102.34 0.16356   0.6859

``` r
# Additive model
SD_logSRic_OS_am_dist <- MASS::glm.nb(row_sum ~ distance_to_ocean_min_m + Location_type, data = SR_env[ocean_stratified_sites,])
summary(SD_logSRic_OS_am_dist)$deviance / summary(SD_logSRic_OS_am_dist)$df.residual
```

    ## [1] 1.215853

``` r
# Check for outliers
outlierTest(SD_logSRic_OS_am_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OTM -1.785697            0.10446           NA

``` r
# Plot residuals
plot(SD_logSRic_OS_am_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-27.png)<!-- -->

``` r
### Mixed & Stratified lakes
SD_logSRic_MS_am_dist <- MASS::glm.nb(row_sum ~ distance_to_ocean_min_m * Location_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_dist)$deviance / summary(SD_logSRic_MS_am_dist)$df.residual
```

    ## [1] 1.275586

``` r
# Plot residuals
plot(SD_logSRic_MS_am_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-28.png)<!-- -->

``` r
# p-values
dist_drop1_results <- drop1(SD_logSRic_MS_am_dist, test = "Chisq")
dist_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ distance_to_ocean_min_m * Location_type
    ##                                       Df Deviance    AIC     LRT Pr(>Chi)
    ## <none>                                     15.307 125.18                 
    ## distance_to_ocean_min_m:Location_type  1   15.744 123.61 0.43658   0.5088

``` r
# Additive model
SD_logSRic_MS_am_dist <- MASS::glm.nb(row_sum ~ distance_to_ocean_min_m + Location_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_dist)$deviance / summary(SD_logSRic_MS_am_dist)$df.residual
```

    ## [1] 1.1762

``` r
# Check for outliers
outlierTest(SD_logSRic_MS_am_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN  2.36238           0.035896      0.57433

``` r
# Plot residuals
plot(SD_logSRic_MS_am_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-29.png)<!-- -->

``` r
### Mixed lakes
SD_logSRic_M_am_dist <- MASS::glm.nb(row_sum ~ distance_to_ocean_min_m, data = SR_env[mixed_lakes,])

### Stratified lakes
SD_logSRic_S_am_dist <- MASS::glm.nb(row_sum ~ distance_to_ocean_min_m, data = SR_env[stratified_lakes,])


#### Max depth
# Interaction
SD_logSRic_am_mxd <- MASS::glm.nb(row_sum ~ max_depth * Location_type, data = SR_env[surveyed_sites,])
summary(SD_logSRic_am_mxd)$deviance / summary(SD_logSRic_am_mxd)$df.residual
```

    ## [1] 1.47533

``` r
# Check for outliers
outlierTest(SD_logSRic_am_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## HLM 2.120647           0.051031           NA

``` r
# Plot residuals
plot(SD_logSRic_am_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-30.png)<!-- -->

``` r
# p-values
mxd_drop1_results <- drop1(SD_logSRic_am_mxd, test = "Chisq")
mxd_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ max_depth * Location_type
    ##                         Df Deviance    AIC    LRT Pr(>Chi)
    ## <none>                       23.605 183.04                
    ## max_depth:Location_type  2   26.553 181.98 2.9474   0.2291

``` r
# Additive model
SD_logSRic_am_mxd <- MASS::glm.nb(row_sum ~ max_depth + Location_type, data = SR_env[surveyed_sites,])
summary(SD_logSRic_am_mxd)$deviance / summary(SD_logSRic_am_mxd)$df.residual
```

    ## [1] 1.32018

``` r
# Check for outliers
outlierTest(SD_logSRic_am_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## TLN 2.646315            0.01697      0.37334

``` r
# Plot residuals
plot(SD_logSRic_am_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-31.png)<!-- -->

``` r
### Ocean & Mixed sites
# Interaction
SD_logSRic_OM_am_mxd <- MASS::glm.nb(row_sum ~ max_depth * Location_type, data = SR_env[ocean_mixed_sites,])
summary(SD_logSRic_OM_am_mxd)$deviance / summary(SD_logSRic_OM_am_mxd)$df.residual
```

    ## [1] 1.418821

``` r
# Plot residuals
plot(SD_logSRic_OM_am_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-32.png)<!-- -->

``` r
# p-values
mxd_drop1_results <- drop1(SD_logSRic_OM_am_mxd, test = "Chisq")
mxd_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ max_depth * Location_type
    ##                         Df Deviance    AIC     LRT Pr(>Chi)
    ## <none>                       14.188 130.23                 
    ## max_depth:Location_type  1   14.970 129.01 0.78224   0.3765

``` r
# Additive model
SD_logSRic_OM_am_mxd <- MASS::glm.nb(row_sum ~ max_depth + Location_type, data = SR_env[ocean_mixed_sites,])
summary(SD_logSRic_OM_am_mxd)$deviance / summary(SD_logSRic_OM_am_mxd)$df.residual
```

    ## [1] 1.30084

``` r
# Check for outliers
outlierTest(SD_logSRic_OM_am_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 1.892254           0.087736           NA

``` r
# Plot residuals
plot(SD_logSRic_OM_am_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-33.png)<!-- -->

``` r
### Ocean & Stratified sites
# Interaction
SD_logSRic_OS_am_mxd <- MASS::glm.nb(row_sum ~ max_depth * Location_type, data = SR_env[ocean_stratified_sites,])
summary(SD_logSRic_OS_am_mxd)$deviance / summary(SD_logSRic_OS_am_mxd)$df.residual
```

    ## [1] 1.507179

``` r
# Plot residuals
plot(SD_logSRic_OS_am_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-34.png)<!-- -->

``` r
# p-values
mxd_drop1_results <- drop1(SD_logSRic_OS_am_mxd, test = "Chisq")
mxd_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ max_depth * Location_type
    ##                         Df Deviance    AIC     LRT Pr(>Chi)
    ## <none>                       15.072 108.09                 
    ## max_depth:Location_type  1   16.040 107.06 0.96769   0.3253

``` r
# Additive model
SD_logSRic_OS_am_mxd <- MASS::glm.nb(row_sum ~ max_depth + Location_type, data = SR_env[ocean_stratified_sites,])
summary(SD_logSRic_OS_am_mxd)$deviance / summary(SD_logSRic_OS_am_mxd)$df.residual
```

    ## [1] 1.361541

``` r
# Check for outliers
outlierTest(SD_logSRic_OS_am_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## TLN 2.383627           0.038374      0.53724

``` r
# Plot residuals
plot(SD_logSRic_OS_am_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-35.png)<!-- -->

``` r
### Mixed & Stratified lakes
SD_logSRic_MS_am_mxd <- MASS::glm.nb(row_sum ~ max_depth * Location_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_mxd)$deviance / summary(SD_logSRic_MS_am_mxd)$df.residual
```

    ## [1] 1.401986

``` r
# Plot residuals
plot(SD_logSRic_MS_am_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-36.png)<!-- -->

``` r
# p-values
mxd_drop1_results <- drop1(SD_logSRic_MS_am_mxd, test = "Chisq")
mxd_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ max_depth * Location_type
    ##                         Df Deviance    AIC    LRT Pr(>Chi)
    ## <none>                       16.824 126.89                
    ## max_depth:Location_type  1   19.216 127.28 2.3923   0.1219

``` r
# Additive model
SD_logSRic_MS_am_mxd <- MASS::glm.nb(row_sum ~ max_depth + Location_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_mxd)$deviance / summary(SD_logSRic_MS_am_mxd)$df.residual
```

    ## [1] 1.293653

``` r
# Check for outliers
outlierTest(SD_logSRic_MS_am_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## TLN 2.464217           0.029807      0.47692

``` r
# Plot residuals
plot(SD_logSRic_MS_am_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-37.png)<!-- -->

``` r
### Mixed lakes
SD_logSRic_M_am_mxd <- MASS::glm.nb(row_sum ~ max_depth, data = SR_env[mixed_lakes,])

### Stratified lakes
SD_logSRic_S_am_mxd <- MASS::glm.nb(row_sum ~ max_depth, data = SR_env[stratified_lakes,])

## Log Area
# Interaction
SD_logSRic_am_lga <- MASS::glm.nb(row_sum ~ logArea * Location_type, data = SR_env[surveyed_sites,])
summary(SD_logSRic_am_lga)$deviance / summary(SD_logSRic_am_lga)$df.residual
```

    ## [1] 1.52311

``` r
# Check for outliers
outlierTest(SD_logSRic_am_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## HLM 2.283285           0.037412      0.82307

``` r
# Plot residuals
plot(SD_logSRic_am_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-38.png)<!-- -->

``` r
# p-values
lga_drop1_results <- drop1(SD_logSRic_am_lga, test = "Chisq")
lga_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ logArea * Location_type
    ##                       Df Deviance    AIC    LRT Pr(>Chi)  
    ## <none>                     24.370 180.58                  
    ## logArea:Location_type  2   30.694 182.91 6.3238  0.04235 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_am_lga)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ logArea * Location_type, data = SR_env[surveyed_sites, 
    ##     ], init.theta = 6.555462935, link = log)
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                      3.39757    0.86215   3.941 8.12e-05 ***
    ## logArea                          0.04820    0.07220   0.668   0.5044    
    ## Location_typeMixed              -2.89254    1.39149  -2.079   0.0376 *  
    ## Location_typeStratified         -1.09007    1.93698  -0.563   0.5736    
    ## logArea:Location_typeMixed       0.30069    0.13238   2.271   0.0231 *  
    ## logArea:Location_typeStratified -0.07677    0.18315  -0.419   0.6751    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(6.5555) family taken to be 1)
    ## 
    ##     Null deviance: 110.04  on 21  degrees of freedom
    ## Residual deviance:  24.37  on 16  degrees of freedom
    ## AIC: 182.58
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  6.56 
    ##           Std. Err.:  2.80 
    ## 
    ##  2 x log-likelihood:  -168.583

``` r
# Get the p-value for Log Area specifically for EACH group
slopes <- emtrends(SD_logSRic_am_lga, ~ Location_type, var = "logArea")
# view the p-values
test(slopes)
```

    ##  Location_type logArea.trend     SE  df z.ratio p.value
    ##  Ocean                0.0482 0.0722 Inf   0.668  0.5044
    ##  Mixed                0.3489 0.1110 Inf   3.144  0.0017
    ##  Stratified          -0.0286 0.1680 Inf  -0.170  0.8652

``` r
### Ocean & Mixed sites
# Interaction
SD_logSRic_OM_am_lga <- MASS::glm.nb(row_sum ~ logArea * Location_type, data = SR_env[ocean_mixed_sites,])
summary(SD_logSRic_OM_am_lga)$deviance / summary(SD_logSRic_OM_am_lga)$df.residual
```

    ## [1] 1.441194

``` r
# Plot residuals
plot(SD_logSRic_OM_am_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-39.png)<!-- -->

``` r
# p-values
lga_drop1_results <- drop1(SD_logSRic_OM_am_lga, test = "Chisq")
lga_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ logArea * Location_type
    ##                       Df Deviance    AIC    LRT Pr(>Chi)   
    ## <none>                     14.412 126.80                   
    ## logArea:Location_type  1   22.326 132.72 7.9141 0.004905 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OM_am_lga)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ logArea * Location_type, data = SR_env[ocean_mixed_sites, 
    ##     ], init.theta = 10.10193096, link = log)
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                 3.39604    0.71611   4.742 2.11e-06 ***
    ## logArea                     0.04833    0.05991   0.807  0.41981    
    ## Location_typeMixed         -2.91927    1.17151  -2.492  0.01271 *  
    ## logArea:Location_typeMixed  0.30339    0.11134   2.725  0.00643 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(10.1019) family taken to be 1)
    ## 
    ##     Null deviance: 29.882  on 13  degrees of freedom
    ## Residual deviance: 14.412  on 10  degrees of freedom
    ## AIC: 128.8
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  10.10 
    ##           Std. Err.:  4.67 
    ## 
    ##  2 x log-likelihood:  -118.803

``` r
# Get the p-value for Log Area specifically for EACH group
slopes <- emtrends(SD_logSRic_OM_am_lga, ~ Location_type, var = "logArea")
# view the p-values
test(slopes)
```

    ##  Location_type logArea.trend     SE  df z.ratio p.value
    ##  Ocean                0.0483 0.0599 Inf   0.807  0.4198
    ##  Mixed                0.3517 0.0939 Inf   3.748  0.0002

``` r
### Ocean & Stratified sites
# Interaction
SD_logSRic_OS_am_lga <- MASS::glm.nb(row_sum ~ logArea * Location_type, data = SR_env[ocean_stratified_sites,])
summary(SD_logSRic_OS_am_lga)$deviance / summary(SD_logSRic_OS_am_lga)$df.residual
```

    ## [1] 1.467081

``` r
# Plot residuals
plot(SD_logSRic_OS_am_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-40.png)<!-- -->

``` r
# p-values
lga_drop1_results <- drop1(SD_logSRic_OS_am_lga, test = "Chisq")
lga_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ logArea * Location_type
    ##                       Df Deviance    AIC     LRT Pr(>Chi)
    ## <none>                     14.671 109.50                 
    ## logArea:Location_type  1   14.797 107.63 0.12618   0.7224

``` r
# Additive model
SD_logSRic_OS_am_lga <- MASS::glm.nb(row_sum ~ logArea + Location_type, data = SR_env[ocean_stratified_sites,])
summary(SD_logSRic_OS_am_lga)$deviance / summary(SD_logSRic_OS_am_lga)$df.residual
```

    ## [1] 1.334658

``` r
# Check for outliers
outlierTest(SD_logSRic_OS_am_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## TLN 1.808691            0.10062           NA

``` r
# Plot residuals
plot(SD_logSRic_OS_am_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-41.png)<!-- -->

``` r
### Mixed & Stratified lakes
SD_logSRic_MS_am_lga <- MASS::glm.nb(row_sum ~ logArea * Location_type, data = SR_env[mixed_stratified_lakes,])
# Plot residuals
plot(SD_logSRic_MS_am_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-42.png)<!-- -->

``` r
# p-values
lga_drop1_results <- drop1(SD_logSRic_MS_am_lga, test = "Chisq")
lga_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ logArea * Location_type
    ##                       Df Deviance    AIC    LRT Pr(>Chi)  
    ## <none>                     18.434 122.96                  
    ## logArea:Location_type  1   22.114 124.64 3.6803  0.05506 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Additive model
SD_logSRic_MS_am_lga <- MASS::glm.nb(row_sum ~ logArea + Location_type, data = SR_env[mixed_stratified_lakes,])
summary(SD_logSRic_MS_am_lga)$deviance / summary(SD_logSRic_MS_am_lga)$df.residual
```

    ## [1] 1.383294

``` r
# Check for outliers
outlierTest(SD_logSRic_MS_am_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## TLN 3.024975           0.010565      0.16904

``` r
# Plot residuals
plot(SD_logSRic_MS_am_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-43.png)<!-- -->

``` r
### Mixed lakes
SD_logSRic_M_am_lga <- MASS::glm.nb(row_sum ~ logArea, data = SR_env[mixed_lakes,])

### Stratified lakes
SD_logSRic_S_am_lga <- MASS::glm.nb(row_sum ~ logArea, data = SR_env[stratified_lakes,])

# Summarize GLM results and calculate pairwise comparisons
summary(SD_logSRic_am_dist)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ distance_to_ocean_min_m + Location_type, 
    ##     data = SR_env[surveyed_sites, ], init.theta = 4.561278054, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              3.992157   0.199505  20.010   <2e-16 ***
    ## distance_to_ocean_min_m -0.006092   0.002762  -2.206   0.0274 *  
    ## Location_typeMixed       0.412158   0.317465   1.298   0.1942    
    ## Location_typeStratified -1.143719   0.483779  -2.364   0.0181 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(4.5613) family taken to be 1)
    ## 
    ##     Null deviance: 82.667  on 21  degrees of freedom
    ## Residual deviance: 21.512  on 18  degrees of freedom
    ## AIC: 181.91
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.56 
    ##           Std. Err.:  1.60 
    ## 
    ##  2 x log-likelihood:  -171.908

``` r
SD_logSRic_amp_dist <- emmeans(SD_logSRic_am_dist, pairwise ~ Location_type, adjust = "bonferroni")
SD_logSRic_amp_dist$contrasts
```

    ##  contrast           estimate    SE  df z.ratio p.value
    ##  Ocean - Mixed        -0.412 0.317 Inf  -1.298  0.5826
    ##  Ocean - Stratified    1.144 0.484 Inf   2.364  0.0542
    ##  Mixed - Stratified    1.556 0.345 Inf   4.506  <.0001
    ## 
    ## Results are given on the log (not the response) scale. 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(SD_logSRic_am_mxd)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ max_depth + Location_type, data = SR_env[surveyed_sites, 
    ##     ], init.theta = 4.647445581, link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              3.64979    0.25352  14.397  < 2e-16 ***
    ## max_depth                0.01995    0.01104   1.807   0.0708 .  
    ## Location_typeMixed       0.04241    0.26219   0.162   0.8715    
    ## Location_typeStratified -2.09908    0.30912  -6.791 1.12e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(4.6474) family taken to be 1)
    ## 
    ##     Null deviance: 83.933  on 21  degrees of freedom
    ## Residual deviance: 23.763  on 18  degrees of freedom
    ## AIC: 183.83
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.65 
    ##           Std. Err.:  1.79 
    ## 
    ##  2 x log-likelihood:  -173.828

``` r
SD_logSRic_amp_mxd <- emmeans(SD_logSRic_am_mxd, pairwise ~ Location_type, adjust = "bonferroni")
SD_logSRic_amp_mxd$contrasts
```

    ##  contrast           estimate    SE  df z.ratio p.value
    ##  Ocean - Mixed       -0.0424 0.262 Inf  -0.162  1.0000
    ##  Ocean - Stratified   2.0991 0.309 Inf   6.791  <.0001
    ##  Mixed - Stratified   2.1415 0.299 Inf   7.166  <.0001
    ## 
    ## Results are given on the log (not the response) scale. 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(SD_logSRic_am_lga)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ logArea * Location_type, data = SR_env[surveyed_sites, 
    ##     ], init.theta = 6.555462935, link = log)
    ## 
    ## Coefficients:
    ##                                 Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                      3.39757    0.86215   3.941 8.12e-05 ***
    ## logArea                          0.04820    0.07220   0.668   0.5044    
    ## Location_typeMixed              -2.89254    1.39149  -2.079   0.0376 *  
    ## Location_typeStratified         -1.09007    1.93698  -0.563   0.5736    
    ## logArea:Location_typeMixed       0.30069    0.13238   2.271   0.0231 *  
    ## logArea:Location_typeStratified -0.07677    0.18315  -0.419   0.6751    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(6.5555) family taken to be 1)
    ## 
    ##     Null deviance: 110.04  on 21  degrees of freedom
    ## Residual deviance:  24.37  on 16  degrees of freedom
    ## AIC: 182.58
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  6.56 
    ##           Std. Err.:  2.80 
    ## 
    ##  2 x log-likelihood:  -168.583

``` r
SD_logSRic_amp_lga <- emmeans(SD_logSRic_am_lga, pairwise ~ Location_type, adjust = "bonferroni")
```

    ## NOTE: Results may be misleading due to involvement in interactions

``` r
SD_logSRic_amp_lga$contrasts
```

    ##  contrast           estimate    SE  df z.ratio p.value
    ##  Ocean - Mixed        -0.243 0.254 Inf  -0.955  1.0000
    ##  Ocean - Stratified    1.891 0.272 Inf   6.956  <.0001
    ##  Mixed - Stratified    2.134 0.253 Inf   8.420  <.0001
    ## 
    ## Results are given on the log (not the response) scale. 
    ## P value adjustment: bonferroni method for 3 tests

``` r
# p-values
dist_drop1_results <- drop1(SD_logSRic_am_dist, test = "Chisq")
dist_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ distance_to_ocean_min_m + Location_type
    ##                         Df Deviance    AIC     LRT Pr(>Chi)    
    ## <none>                       21.512 179.91                     
    ## distance_to_ocean_min_m  1   26.625 183.02  5.1133  0.02374 *  
    ## Location_type            2   45.265 199.66 23.7536 6.95e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
dist_p_values <- dist_drop1_results[["Pr(>Chi)"]]
mxd_drop1_results <- drop1(SD_logSRic_am_mxd, test = "Chisq")
mxd_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ max_depth + Location_type
    ##               Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>             23.763 181.83                     
    ## max_depth      1   27.017 183.08  3.254   0.07125 .  
    ## Location_type  2   83.806 237.87 60.043 9.159e-14 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
mxd_p_values <- mxd_drop1_results[["Pr(>Chi)"]]
lga_drop1_results <- drop1(SD_logSRic_am_lga, test = "Chisq")
lga_drop1_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ logArea * Location_type
    ##                       Df Deviance    AIC    LRT Pr(>Chi)  
    ## <none>                     24.370 180.58                  
    ## logArea:Location_type  2   30.694 182.91 6.3238  0.04235 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
lga_p_values <- lga_drop1_results[["Pr(>Chi)"]]
p_values <- c(dist_p_values[2:3], mxd_p_values[2:3], lga_p_values[2])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.187125e-01 3.474931e-05 3.562572e-01 4.579724e-13 2.117284e-01

``` r
# Summarize OM lm results
dist_drop1_OM_results <- drop1(SD_logSRic_OM_am_dist, test = "Chisq")
dist_drop1_OM_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ distance_to_ocean_min_m + Location_type
    ##                         Df Deviance    AIC     LRT Pr(>Chi)
    ## <none>                       14.489 135.10                 
    ## distance_to_ocean_min_m  1   14.620 133.23 0.13006   0.7184
    ## Location_type            1   14.604 133.21 0.11402   0.7356

``` r
summary(SD_logSRic_OM_am_dist)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ distance_to_ocean_min_m + Location_type, 
    ##     data = SR_env[ocean_mixed_sites, ], init.theta = 4.434189514, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              3.975889   0.203058  19.580   <2e-16 ***
    ## distance_to_ocean_min_m -0.002197   0.005805  -0.378    0.705    
    ## Location_typeMixed       0.158634   0.459124   0.346    0.730    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(4.4342) family taken to be 1)
    ## 
    ##     Null deviance: 14.625  on 13  degrees of freedom
    ## Residual deviance: 14.490  on 11  degrees of freedom
    ## AIC: 137.1
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.43 
    ##           Std. Err.:  1.77 
    ## 
    ##  2 x log-likelihood:  -129.097

``` r
mxd_drop1_OM_results <- drop1(SD_logSRic_OM_am_mxd, test = "Chisq")
mxd_drop1_OM_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ max_depth + Location_type
    ##               Df Deviance    AIC    LRT Pr(>Chi)   
    ## <none>             14.309 128.99                   
    ## max_depth      1   22.076 134.76 7.7672  0.00532 **
    ## Location_type  1   14.356 127.04 0.0467  0.82896   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OM_am_mxd)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ max_depth + Location_type, data = SR_env[ocean_mixed_sites, 
    ##     ], init.theta = 7.053116644, link = log)
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)         3.46653    0.23203  14.940   <2e-16 ***
    ## max_depth           0.03281    0.01124   2.920   0.0035 ** 
    ## Location_typeMixed  0.04702    0.21801   0.216   0.8292    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(7.0531) family taken to be 1)
    ## 
    ##     Null deviance: 22.085  on 13  degrees of freedom
    ## Residual deviance: 14.309  on 11  degrees of freedom
    ## AIC: 130.99
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  7.05 
    ##           Std. Err.:  3.02 
    ## 
    ##  2 x log-likelihood:  -122.994

``` r
lga_drop1_OM_results <- drop1(SD_logSRic_OM_am_lga, test = "Chisq")
lga_drop1_OM_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ logArea * Location_type
    ##                       Df Deviance    AIC    LRT Pr(>Chi)   
    ## <none>                     14.412 126.80                   
    ## logArea:Location_type  1   22.326 132.72 7.9141 0.004905 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OM_am_lga)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ logArea * Location_type, data = SR_env[ocean_mixed_sites, 
    ##     ], init.theta = 10.10193096, link = log)
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                 3.39604    0.71611   4.742 2.11e-06 ***
    ## logArea                     0.04833    0.05991   0.807  0.41981    
    ## Location_typeMixed         -2.91927    1.17151  -2.492  0.01271 *  
    ## logArea:Location_typeMixed  0.30339    0.11134   2.725  0.00643 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(10.1019) family taken to be 1)
    ## 
    ##     Null deviance: 29.882  on 13  degrees of freedom
    ## Residual deviance: 14.412  on 10  degrees of freedom
    ## AIC: 128.8
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  10.10 
    ##           Std. Err.:  4.67 
    ## 
    ##  2 x log-likelihood:  -118.803

``` r
# p-values
dist_p_values <- dist_drop1_OM_results[["Pr(>Chi)"]]
mxd_p_values <- mxd_drop1_OM_results[["Pr(>Chi)"]]
lga_p_values <- lga_drop1_OM_results[["Pr(>Chi)"]]
p_values <- c(dist_p_values[2:3], mxd_p_values[2:3], lga_p_values[2])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.00000000 1.00000000 0.02660129 1.00000000 0.02452503

``` r
# Summarize OS lm results
dist_drop1_OS_results <- drop1(SD_logSRic_OS_am_dist, test = "Chisq")
dist_drop1_OS_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ distance_to_ocean_min_m + Location_type
    ##                         Df Deviance    AIC    LRT Pr(>Chi)  
    ## <none>                       13.374 102.34                  
    ## distance_to_ocean_min_m  1   19.567 106.54 6.1925  0.01283 *
    ## Location_type            1   17.894 104.86 4.5197  0.03351 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OS_am_dist)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ distance_to_ocean_min_m + Location_type, 
    ##     data = SR_env[ocean_stratified_sites, ], init.theta = 5.767550426, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              3.996011   0.179375  22.277   <2e-16 ***
    ## distance_to_ocean_min_m -0.007001   0.002961  -2.364   0.0181 *  
    ## Location_typeStratified -1.020421   0.483900  -2.109   0.0350 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(5.7676) family taken to be 1)
    ## 
    ##     Null deviance: 74.605  on 13  degrees of freedom
    ## Residual deviance: 13.374  on 11  degrees of freedom
    ## AIC: 104.34
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  5.77 
    ##           Std. Err.:  2.91 
    ## 
    ##  2 x log-likelihood:  -96.343

``` r
mxd_drop1_OS_results <- drop1(SD_logSRic_OS_am_mxd, test = "Chisq")
mxd_drop1_OS_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ max_depth + Location_type
    ##               Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>             14.977 107.03                     
    ## max_depth      1   15.812 105.86  0.835    0.3608    
    ## Location_type  1   58.260 148.31 43.283 4.738e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OS_am_mxd)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ max_depth + Location_type, data = SR_env[ocean_stratified_sites, 
    ##     ], init.theta = 4.243964504, link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              3.77383    0.28348  13.313  < 2e-16 ***
    ## max_depth                0.01181    0.01364   0.866    0.386    
    ## Location_typeStratified -2.03470    0.32693  -6.224 4.86e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(4.244) family taken to be 1)
    ## 
    ##     Null deviance: 58.845  on 13  degrees of freedom
    ## Residual deviance: 14.977  on 11  degrees of freedom
    ## AIC: 109.02
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.24 
    ##           Std. Err.:  2.19 
    ## 
    ##  2 x log-likelihood:  -101.025

``` r
lga_drop1_OS_results <- drop1(SD_logSRic_OS_am_lga, test = "Chisq")
lga_drop1_OS_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ logArea + Location_type
    ##               Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>             14.681 107.63                     
    ## logArea        1   14.866 105.81  0.184    0.6676    
    ## Location_type  1   48.873 139.82 34.192 4.993e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_OS_am_lga)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ logArea + Location_type, data = SR_env[ocean_stratified_sites, 
    ##     ], init.theta = 3.89825012, link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              3.55204    0.99493   3.570 0.000357 ***
    ## logArea                  0.03503    0.08310   0.422 0.673314    
    ## Location_typeStratified -1.89553    0.32941  -5.754  8.7e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(3.8983) family taken to be 1)
    ## 
    ##     Null deviance: 54.967  on 13  degrees of freedom
    ## Residual deviance: 14.681  on 11  degrees of freedom
    ## AIC: 109.63
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  3.90 
    ##           Std. Err.:  1.92 
    ## 
    ##  2 x log-likelihood:  -101.628

``` r
# p-values
dist_p_values <- dist_drop1_OS_results[["Pr(>Chi)"]]
mxd_p_values <- mxd_drop1_OS_results[["Pr(>Chi)"]]
lga_p_values <- lga_drop1_OS_results[["Pr(>Chi)"]]
p_values <- c(dist_p_values[2:3], mxd_p_values[2:3], lga_p_values[2:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 7.697408e-02 2.010450e-01 1.000000e+00 2.842727e-10 1.000000e+00
    ## [6] 2.995646e-08

``` r
# Summarize MS lm results
dist_drop1_MS_results <- drop1(SD_logSRic_MS_am_dist, test = "Chisq")
dist_drop1_MS_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ distance_to_ocean_min_m + Location_type
    ##                         Df Deviance    AIC     LRT  Pr(>Chi)    
    ## <none>                       15.291 123.61                      
    ## distance_to_ocean_min_m  1   20.232 126.55  4.9411   0.02623 *  
    ## Location_type            1   33.759 140.07 18.4681 1.728e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_MS_am_dist)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ distance_to_ocean_min_m + Location_type, 
    ##     data = SR_env[mixed_stratified_lakes, ], init.theta = 4.183207624, 
    ##     link = log)
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              4.412629   0.264721  16.669  < 2e-16 ***
    ## distance_to_ocean_min_m -0.006208   0.002864  -2.167   0.0302 *  
    ## Location_typeStratified -1.548477   0.357729  -4.329  1.5e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(4.1832) family taken to be 1)
    ## 
    ##     Null deviance: 67.141  on 15  degrees of freedom
    ## Residual deviance: 15.291  on 13  degrees of freedom
    ## AIC: 125.6
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.18 
    ##           Std. Err.:  1.76 
    ## 
    ##  2 x log-likelihood:  -117.605

``` r
mxd_drop1_MS_results <- drop1(SD_logSRic_MS_am_mxd, test = "Chisq")
mxd_drop1_MS_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ max_depth + Location_type
    ##               Df Deviance    AIC    LRT Pr(>Chi)    
    ## <none>             16.817 127.12                    
    ## max_depth      1   17.908 126.21  1.091   0.2963    
    ## Location_type  1   56.771 165.07 39.954  2.6e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_MS_am_mxd)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ max_depth + Location_type, data = SR_env[mixed_stratified_lakes, 
    ##     ], init.theta = 3.573806988, link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              3.74226    0.27592  13.563  < 2e-16 ***
    ## max_depth                0.01644    0.01520   1.081     0.28    
    ## Location_typeStratified -2.10825    0.34465  -6.117 9.53e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(3.5738) family taken to be 1)
    ## 
    ##     Null deviance: 58.987  on 15  degrees of freedom
    ## Residual deviance: 16.817  on 13  degrees of freedom
    ## AIC: 129.12
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  3.57 
    ##           Std. Err.:  1.56 
    ## 
    ##  2 x log-likelihood:  -121.116

``` r
lga_drop1_MS_results <- drop1(SD_logSRic_MS_am_lga, test = "Chisq")
lga_drop1_MS_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ logArea + Location_type
    ##               Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>             17.983 124.25                     
    ## logArea        1   22.974 127.25  4.991   0.02548 *  
    ## Location_type  1   72.939 177.21 54.956 1.233e-13 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_MS_am_lga)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ logArea + Location_type, data = SR_env[mixed_stratified_lakes, 
    ##     ], init.theta = 4.952902334, link = log)
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)               1.6405     1.0181   1.611   0.1071    
    ## logArea                   0.2333     0.1035   2.254   0.0242 *  
    ## Location_typeStratified  -2.0083     0.2731  -7.354 1.92e-13 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(4.9529) family taken to be 1)
    ## 
    ##     Null deviance: 76.863  on 15  degrees of freedom
    ## Residual deviance: 17.983  on 13  degrees of freedom
    ## AIC: 126.25
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  4.95 
    ##           Std. Err.:  2.56 
    ## 
    ##  2 x log-likelihood:  -118.254

``` r
# p-values
dist_p_values <- dist_drop1_MS_results[["Pr(>Chi)"]]
mxd_p_values <- mxd_drop1_MS_results[["Pr(>Chi)"]]
lga_p_values <- lga_drop1_MS_results[["Pr(>Chi)"]]
p_values <- c(dist_p_values[2:3], mxd_p_values[2:3], lga_p_values[2:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.573536e-01 1.036619e-04 1.000000e+00 1.560128e-09 1.528719e-01
    ## [6] 7.396315e-13

``` r
# Summarize M lm results
dist_drop1_M_results <- drop1(SD_logSRic_M_am_dist, test = "Chisq")
dist_drop1_M_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ distance_to_ocean_min_m
    ##                         Df Deviance    AIC     LRT Pr(>Chi)
    ## <none>                       8.3230 78.794                 
    ## distance_to_ocean_min_m  1   8.4529 76.924 0.12991   0.7185

``` r
summary(SD_logSRic_M_am_dist)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ distance_to_ocean_min_m, data = SR_env[mixed_lakes, 
    ##     ], init.theta = 3.803395785, link = log)
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              4.151359   0.476926   8.704   <2e-16 ***
    ## distance_to_ocean_min_m -0.002444   0.006432  -0.380    0.704    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(3.8034) family taken to be 1)
    ## 
    ##     Null deviance: 8.4529  on 7  degrees of freedom
    ## Residual deviance: 8.3230  on 6  degrees of freedom
    ## AIC: 80.794
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  3.80 
    ##           Std. Err.:  1.98 
    ## 
    ##  2 x log-likelihood:  -74.794

``` r
mxd_drop1_M_results <- drop1(SD_logSRic_M_am_mxd, test = "Chisq")
mxd_drop1_M_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ max_depth
    ##           Df Deviance    AIC   LRT Pr(>Chi)  
    ## <none>         8.0675 74.796                 
    ## max_depth  1  13.3144 78.043 5.247  0.02199 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_M_am_mxd)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ max_depth, data = SR_env[mixed_lakes, 
    ##     ], init.theta = 6.305525583, link = log)
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  3.33862    0.27637   12.08  < 2e-16 ***
    ## max_depth    0.04581    0.01755    2.61  0.00905 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(6.3055) family taken to be 1)
    ## 
    ##     Null deviance: 13.3144  on 7  degrees of freedom
    ## Residual deviance:  8.0675  on 6  degrees of freedom
    ## AIC: 76.796
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  6.31 
    ##           Std. Err.:  3.46 
    ## 
    ##  2 x log-likelihood:  -70.796

``` r
lga_drop1_M_results <- drop1(SD_logSRic_M_am_lga, test = "Chisq")
lga_drop1_M_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ logArea
    ##         Df Deviance    AIC    LRT  Pr(>Chi)    
    ## <none>        8.876 67.495                     
    ## logArea  1   36.477 93.096 27.601 1.491e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_M_am_lga)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ logArea, data = SR_env[mixed_lakes, 
    ##     ], init.theta = 22.83959653, link = log)
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   0.4010     0.7156   0.560    0.575    
    ## logArea       0.3593     0.0718   5.004  5.6e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(22.8396) family taken to be 1)
    ## 
    ##     Null deviance: 36.4772  on 7  degrees of freedom
    ## Residual deviance:  8.8761  on 6  degrees of freedom
    ## AIC: 69.495
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  22.8 
    ##           Std. Err.:  19.0 
    ## 
    ##  2 x log-likelihood:  -63.495

``` r
# p-values
dist_p_values <- dist_drop1_M_results[["Pr(>Chi)"]]
mxd_p_values <- mxd_drop1_M_results[["Pr(>Chi)"]]
lga_p_values <- lga_drop1_M_results[["Pr(>Chi)"]]
p_values <- c(dist_p_values[2], mxd_p_values[2], lga_p_values[2])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 1.000000e+00 6.595503e-02 4.472643e-07

``` r
# Summarize S lm results
dist_drop1_S_results <- drop1(SD_logSRic_S_am_dist, test = "Chisq")
dist_drop1_S_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ distance_to_ocean_min_m
    ##                         Df Deviance    AIC    LRT Pr(>Chi)  
    ## <none>                       7.3026 46.143                  
    ## distance_to_ocean_min_m  1  13.9174 50.757 6.6148  0.01011 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(SD_logSRic_S_am_dist)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ distance_to_ocean_min_m, data = SR_env[stratified_lakes, 
    ##     ], init.theta = 6.200945718, link = log)
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)              3.000445   0.453345   6.618 3.63e-11 ***
    ## distance_to_ocean_min_m -0.007181   0.002945  -2.438   0.0148 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(6.2009) family taken to be 1)
    ## 
    ##     Null deviance: 13.9174  on 7  degrees of freedom
    ## Residual deviance:  7.3026  on 6  degrees of freedom
    ## AIC: 48.143
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  6.20 
    ##           Std. Err.:  5.54 
    ## 
    ##  2 x log-likelihood:  -42.143

``` r
mxd_drop1_S_results <- drop1(SD_logSRic_S_am_mxd, test = "Chisq")
mxd_drop1_S_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ max_depth
    ##           Df Deviance    AIC       LRT Pr(>Chi)
    ## <none>         7.8876 51.017                   
    ## max_depth  1   7.8890 49.018 0.0013699   0.9705

``` r
summary(SD_logSRic_S_am_mxd)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ max_depth, data = SR_env[stratified_lakes, 
    ##     ], init.theta = 2.581247716, link = log)
    ## 
    ## Coefficients:
    ##               Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  2.0343083  0.5918560   3.437 0.000588 ***
    ## max_depth   -0.0008233  0.0226305  -0.036 0.970980    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(2.5812) family taken to be 1)
    ## 
    ##     Null deviance: 7.8890  on 7  degrees of freedom
    ## Residual deviance: 7.8876  on 6  degrees of freedom
    ## AIC: 53.017
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.58 
    ##           Std. Err.:  1.64 
    ## 
    ##  2 x log-likelihood:  -47.017

``` r
lga_drop1_S_results <- drop1(SD_logSRic_S_am_lga, test = "Chisq")
lga_drop1_S_results
```

    ## Single term deletions
    ## 
    ## Model:
    ## row_sum ~ logArea
    ##         Df Deviance    AIC      LRT Pr(>Chi)
    ## <none>       7.8857 51.002                  
    ## logArea  1   7.9014 49.018 0.015705   0.9003

``` r
summary(SD_logSRic_S_am_lga)
```

    ## 
    ## Call:
    ## MASS::glm.nb(formula = row_sum ~ logArea, data = SR_env[stratified_lakes, 
    ##     ], init.theta = 2.586773647, link = log)
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)
    ## (Intercept)  2.30829    2.34735   0.983    0.325
    ## logArea     -0.02864    0.22763  -0.126    0.900
    ## 
    ## (Dispersion parameter for Negative Binomial(2.5868) family taken to be 1)
    ## 
    ##     Null deviance: 7.9014  on 7  degrees of freedom
    ## Residual deviance: 7.8857  on 6  degrees of freedom
    ## AIC: 53.002
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  2.59 
    ##           Std. Err.:  1.64 
    ## 
    ##  2 x log-likelihood:  -47.002

``` r
# p-values
dist_p_values <- dist_drop1_S_results[["Pr(>Chi)"]]
mxd_p_values <- mxd_drop1_S_results[["Pr(>Chi)"]]
lga_p_values <- lga_drop1_S_results[["Pr(>Chi)"]]
p_values <- c(dist_p_values[2], mxd_p_values[2], lga_p_values[2])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.03034067 1.00000000 1.00000000

``` r
##### SRic ANCOVAs
#### Environmental
### Surveyed sites
## Temperature
# Interaction
SD_SRic_am_temp <- aov(row_sum ~ temperature_median * Location_type, data = SR_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_am_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_am_temp)
    ## W = 0.91852, p-value = 0.1063

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_am_temp) ~ SR_env[surveyed_sites_env,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  3.4144 0.05822 .
    ##       16                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_am_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## NLU -3.626467          0.0034731      0.06599

``` r
# Plot residuals
plot(SD_SRic_am_temp)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-44.png)<!-- -->

``` r
# Summarize the ANCOVA results
summary(SD_SRic_am_temp)
```

    ##                                  Df Sum Sq Mean Sq F value   Pr(>F)    
    ## temperature_median                1   1618    1618   4.099 0.063963 .  
    ## Location_type                     2   9895    4947  12.530 0.000928 ***
    ## temperature_median:Location_type  2   1878     939   2.378 0.131806    
    ## Residuals                        13   5133     395                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# ANCOVA
SD_SRic_am_temp <- aov(row_sum ~ temperature_median + Location_type, data = SR_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_am_temp))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_am_temp)
    ## W = 0.9529, p-value = 0.4421

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_am_temp) ~ SR_env[surveyed_sites_env,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  3.3867 0.05937 .
    ##       16                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_am_temp)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.313815          0.0051196     0.097272

``` r
# Plot residuals
plot(SD_SRic_am_temp)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-45.png)<!-- -->

``` r
## Salinity
# Interaction
SD_SRic_am_sal <- aov(row_sum ~ salinity_median * Location_type, data = SR_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_am_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_am_sal)
    ## W = 0.9208, p-value = 0.1171

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_am_sal) ~ SR_env[surveyed_sites_env,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  2.5635 0.1082
    ##       16

``` r
# Check for outliers
outlierTest(SD_SRic_am_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.517177          0.0042454     0.080663

``` r
# Plot residuals
plot(SD_SRic_am_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-46.png)<!-- -->

``` r
# Summarize the ANCOVA results
summary(SD_SRic_am_sal)
```

    ##                               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median                1   9173    9173  22.544 0.000381 ***
    ## Location_type                  2   2620    1310   3.220 0.073142 .  
    ## salinity_median:Location_type  2   1442     721   1.772 0.208663    
    ## Residuals                     13   5289     407                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# ANCOVA
SD_SRic_am_sal <- aov(row_sum ~ salinity_median + Location_type, data = SR_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_am_sal))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_am_sal)
    ## W = 0.93527, p-value = 0.2164

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_am_sal) ~ SR_env[surveyed_sites_env,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  4.3382 0.03124 *
    ##       16                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_am_sal)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.436669          0.0040095     0.076181

``` r
# Plot residuals
plot(SD_SRic_am_sal)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-47.png)<!-- -->

``` r
## Oxygen
# Interaction
SD_SRic_am_oxy <- aov(row_sum ~ oxygen_median * Location_type, data = SR_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_am_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_am_oxy)
    ## W = 0.88176, p-value = 0.02304

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_am_oxy) ~ SR_env[surveyed_sites_env,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value   Pr(>F)   
    ## group  2  6.8003 0.007287 **
    ##       16                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_am_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.514228          0.0042685     0.081102

``` r
# Plot residuals
plot(SD_SRic_am_oxy)
```

    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced
    ## Warning in sqrt(crit * p * (1 - hh)/hh): NaNs produced

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-48.png)<!-- -->

``` r
# Summarize the ANCOVA results
summary(SD_SRic_am_oxy)
```

    ##                             Df Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median                1   9605    9605  29.593 0.000113 ***
    ## Location_type                2   2160    1080   3.328 0.068075 .  
    ## oxygen_median:Location_type  2   2540    1270   3.913 0.046733 *  
    ## Residuals                   13   4219     325                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# ANCOVA
SD_SRic_am_oxy <- aov(row_sum ~ oxygen_median + Location_type, data = SR_env[surveyed_sites_env,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_am_oxy))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_am_oxy)
    ## W = 0.94821, p-value = 0.3682

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_am_oxy) ~ SR_env[surveyed_sites_env,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  2.2869 0.1338
    ##       16

``` r
# Check for outliers
outlierTest(SD_SRic_am_oxy)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.254152          0.0057649      0.10953

``` r
# Plot residuals
plot(SD_SRic_am_oxy)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-49.png)<!-- -->

``` r
# Summarize ANCOVA results and calculate pairwise comparisons
summary(SD_SRic_am_temp)
```

    ##                    Df Sum Sq Mean Sq F value  Pr(>F)   
    ## temperature_median  1   1618    1618   3.463 0.08249 . 
    ## Location_type       2   9895    4947  10.585 0.00136 **
    ## Residuals          15   7011     467                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
SD_SRic_amp_temp <- emmeans(SD_SRic_am_temp, pairwise ~ Location_type, adjust = "bonferroni")
SD_SRic_amp_temp$contrasts
```

    ##  contrast           estimate   SE df t.ratio p.value
    ##  Ocean - Mixed          11.0 14.7 15   0.746  1.0000
    ##  Ocean - Stratified     56.0 14.8 15   3.779  0.0055
    ##  Mixed - Stratified     45.1 11.5 15   3.912  0.0042
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(SD_SRic_am_sal)
```

    ##                 Df Sum Sq Mean Sq F value   Pr(>F)    
    ## salinity_median  1   9173    9173  20.440 0.000406 ***
    ## Location_type    2   2620    1310   2.919 0.084953 .  
    ## Residuals       15   6731     449                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
SD_SRic_amp_sal <- emmeans(SD_SRic_am_sal, pairwise ~ Location_type, adjust = "bonferroni")
SD_SRic_amp_sal$contrasts
```

    ##  contrast           estimate   SE df t.ratio p.value
    ##  Ocean - Mixed          9.36 14.4 15   0.650  1.0000
    ##  Ocean - Stratified    45.05 19.9 15   2.268  0.1156
    ##  Mixed - Stratified    35.69 16.4 15   2.182  0.1362
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(SD_SRic_am_oxy)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## oxygen_median  1   9605    9605  21.314 0.000336 ***
    ## Location_type  2   2160    1080   2.397 0.124957    
    ## Residuals     15   6759     451                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
SD_SRic_amp_oxy <- emmeans(SD_SRic_am_oxy, pairwise ~ Location_type, adjust = "bonferroni")
SD_SRic_amp_oxy$contrasts
```

    ##  contrast           estimate   SE df t.ratio p.value
    ##  Ocean - Mixed          4.87 15.9 15   0.306  1.0000
    ##  Ocean - Stratified    40.93 24.2 15   1.694  0.3327
    ##  Mixed - Stratified    36.06 16.5 15   2.189  0.1346
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
# p-values
temp_p_values <- summary(SD_SRic_am_temp)[[1]][, "Pr(>F)"]
sal_p_values <- summary(SD_SRic_am_sal)[[1]][, "Pr(>F)"]
oxy_p_values <- summary(SD_SRic_am_oxy)[[1]][, "Pr(>F)"]
p_values <- c(temp_p_values[1:3], sal_p_values[1:3], oxy_p_values[1:3])
(adjusted_p_values <- p.adjust(p_values, method = "bonferroni"))
```

    ## [1] 0.494932098 0.008150692          NA 0.002435311 0.509718354          NA
    ## [7] 0.002014076 0.749743578          NA

``` r
### Mixed lakes
SD_SRic_env_M_lm <- lm(row_sum ~ salinity_median + oxygen_median + temperature_median, SR_env[mixed_lakes,])
shapiro.test(residuals(SD_SRic_env_M_lm))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_env_M_lm)
    ## W = 0.9042, p-value = 0.315

``` r
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
car::Anova(SD_SRic_env_M_lm, type = 3)
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
### Stratified lakes
SD_SRic_env_S_lm <- lm(row_sum ~ salinity_median + oxygen_median + temperature_median, SR_env[stratified_lakes,])
shapiro.test(residuals(SD_SRic_env_S_lm))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_env_S_lm)
    ## W = 0.8672, p-value = 0.1415

``` r
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
car::Anova(SD_SRic_env_S_lm, type = 3)
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
### Surveyed sites
## Distance
# Interaction
SD_SRic_am_dist <- aov(row_sum ~ distance_to_ocean_min_m * Location_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_am_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_am_dist)
    ## W = 0.93911, p-value = 0.1898

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_am_dist) ~ SR_env[surveyed_sites,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  4.3497 0.02784 *
    ##       19                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_am_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.483922          0.0033317     0.069965

``` r
# Plot residuals
plot(SD_SRic_am_dist)
```

    ## Warning: not plotting observations with leverage one:
    ##   19

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-50.png)<!-- -->

``` r
# Summarize the ANCOVA results
summary(SD_SRic_am_dist)
```

    ##                                       Df Sum Sq Mean Sq F value  Pr(>F)   
    ## distance_to_ocean_min_m                1   6762    6762  11.938 0.00326 **
    ## Location_type                          2   4154    2077   3.667 0.04889 * 
    ## distance_to_ocean_min_m:Location_type  2     32      16   0.028 0.97262   
    ## Residuals                             16   9063     566                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# ANCOVA
SD_SRic_am_dist <- aov(row_sum ~ distance_to_ocean_min_m + Location_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_am_dist))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_am_dist)
    ## W = 0.94125, p-value = 0.2101

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_am_dist) ~ SR_env[surveyed_sites,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value  Pr(>F)  
    ## group  2  4.6736 0.02235 *
    ##       19                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Check for outliers
outlierTest(SD_SRic_am_dist)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.170152          0.0055954       0.1231

``` r
# Plot residuals
plot(SD_SRic_am_dist)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-51.png)<!-- -->

``` r
## Max Depth
# Interaction
SD_SRic_am_mxd <- aov(row_sum ~ max_depth * Location_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_am_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_am_mxd)
    ## W = 0.90823, p-value = 0.04356

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_am_mxd) ~ SR_env[surveyed_sites,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  2.2575 0.1319
    ##       19

``` r
# Check for outliers
outlierTest(SD_SRic_am_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.386996          0.0040648     0.089426

``` r
# Plot residuals
plot(SD_SRic_am_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-52.png)<!-- -->

``` r
# Summarize the ANCOVA results
summary(SD_SRic_am_mxd)
```

    ##                         Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth                1     33      33   0.093    0.764    
    ## Location_type            2  12617    6308  17.595 9.11e-05 ***
    ## max_depth:Location_type  2   1624     812   2.264    0.136    
    ## Residuals               16   5737     359                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# ANCOVA
SD_SRic_am_mxd <- aov(row_sum ~ max_depth + Location_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_am_mxd))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_am_mxd)
    ## W = 0.95124, p-value = 0.3342

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_am_mxd) ~ SR_env[surveyed_sites,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  1.0674 0.3636
    ##       19

``` r
# Check for outliers
outlierTest(SD_SRic_am_mxd)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##     rstudent unadjusted p-value Bonferroni p
    ## LLN 3.218891          0.0050385      0.11085

``` r
# Plot residuals
plot(SD_SRic_am_mxd)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-53.png)<!-- -->

``` r
## Log Area
# Interaction
SD_SRic_am_lga <- aov(row_sum ~ logArea * Location_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_am_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_am_lga)
    ## W = 0.92682, p-value = 0.1053

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_am_lga) ~ SR_env[surveyed_sites,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  2.0972 0.1503
    ##       19

``` r
# Check for outliers
outlierTest(SD_SRic_am_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OCO -3.113325          0.0071207      0.15666

``` r
# Plot residuals
plot(SD_SRic_am_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-54.png)<!-- -->

``` r
# Summarize the ANCOVA results
summary(SD_SRic_am_lga)
```

    ##                       Df Sum Sq Mean Sq F value   Pr(>F)    
    ## logArea                1   2166    2166   7.338   0.0155 *  
    ## Location_type          2  10588    5294  17.931 8.21e-05 ***
    ## logArea:Location_type  2   2532    1266   4.287   0.0323 *  
    ## Residuals             16   4724     295                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# ANCOVA
SD_SRic_am_lga <- aov(row_sum ~ logArea + Location_type, data = SR_env[surveyed_sites,])
# Shapiro-Wilk test on residuals
shapiro.test(residuals(SD_SRic_am_lga))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_am_lga)
    ## W = 0.97073, p-value = 0.7279

``` r
# Levene’s test for homogeneity of variance
car::leveneTest(residuals(SD_SRic_am_lga) ~ SR_env[surveyed_sites,"Location_type"])
```

    ## Levene's Test for Homogeneity of Variance (center = median)
    ##       Df F value Pr(>F)
    ## group  2  1.0393  0.373
    ##       19

``` r
# Check for outliers
outlierTest(SD_SRic_am_lga)
```

    ## No Studentized residuals with Bonferroni p < 0.05
    ## Largest |rstudent|:
    ##      rstudent unadjusted p-value Bonferroni p
    ## OCO -2.994945           0.008142      0.17912

``` r
# Plot residuals
plot(SD_SRic_am_lga)
```

![](SD_analyses_files/figure-gfm/SD%20alpha%20lm%20&%20ANCOVA-55.png)<!-- -->

``` r
# Summarize ANCOVA results and calculate pairwise comparisons
summary(SD_SRic_am_dist)
```

    ##                         Df Sum Sq Mean Sq F value Pr(>F)   
    ## distance_to_ocean_min_m  1   6762    6762  13.383 0.0018 **
    ## Location_type            2   4154    2077   4.111 0.0339 * 
    ## Residuals               18   9095     505                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
SD_SRic_amp_dist <- emmeans(SD_SRic_am_dist, pairwise ~ Location_type, adjust = "bonferroni")
SD_SRic_amp_dist$contrasts
```

    ##  contrast           estimate   SE df t.ratio p.value
    ##  Ocean - Mixed         -4.97 13.9 18  -0.358  1.0000
    ##  Ocean - Stratified    36.12 19.9 18   1.816  0.2580
    ##  Mixed - Stratified    41.09 14.4 18   2.849  0.0320
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(SD_SRic_am_mxd)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## max_depth      1     33      33   0.082 0.778361    
    ## Location_type  2  12617    6308  15.428 0.000125 ***
    ## Residuals     18   7360     409                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
SD_SRic_amp_mxd <- emmeans(SD_SRic_am_mxd, pairwise ~ Location_type, adjust = "bonferroni")
SD_SRic_amp_mxd$contrasts
```

    ##  contrast           estimate   SE df t.ratio p.value
    ##  Ocean - Mixed         -2.22 10.9 18  -0.203  1.0000
    ##  Ocean - Stratified    54.00 11.6 18   4.641  0.0006
    ##  Mixed - Stratified    56.22 11.1 18   5.069  0.0002
    ## 
    ## P value adjustment: bonferroni method for 3 tests

``` r
summary(SD_SRic_am_lga)
```

    ##               Df Sum Sq Mean Sq F value   Pr(>F)    
    ## logArea        1   2166    2166   5.375 0.032401 *  
    ## Location_type  2  10588    5294  13.134 0.000304 ***
    ## Residuals     18   7256     403                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
SD_SRic_amp_lga <- emmeans(SD_SRic_am_lga, pairwise ~ Location_type, adjust = "bonferroni")
SD_SRic_amp_lga$contrasts
```

    ##  contrast           estimate   SE df t.ratio p.value
    ##  Ocean - Mixed         -12.9 12.1 18  -1.070  0.8966
    ##  Ocean - Stratified     37.0 11.5 18   3.227  0.0140
    ##  Mixed - Stratified     49.9 10.2 18   4.912  0.0003
    ## 
    ## P value adjustment: bonferroni method for 3 tests

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
# Mixed lakes
SD_SRic_geo_M_lm <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[mixed_lakes,])
shapiro.test(residuals(SD_SRic_geo_M_lm))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_geo_M_lm)
    ## W = 0.97661, p-value = 0.9442

``` r
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
car::Anova(SD_SRic_geo_M_lm, type = 3)
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
### Stratified lakes
SD_SRic_geo_S_lm <- lm(row_sum ~ distance_to_ocean_min_m + max_depth + logArea, SR_env[stratified_lakes,])
shapiro.test(residuals(SD_SRic_geo_S_lm))
```

    ## 
    ##  Shapiro-Wilk normality test
    ## 
    ## data:  residuals(SD_SRic_geo_S_lm)
    ## W = 0.96411, p-value = 0.8482

``` r
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
car::Anova(SD_SRic_geo_S_lm, type = 3)
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
par(mfrow=c(1,1))
```

## SD beta diversity

### SD beta diversity distance calculations plus dispersion and PERMANOVA

``` r
### Regular
# Reference
SD_beta_ref_dist <- BAT::beta(presabs_lake, abund = F)
(SD_beta_ref_BD <- betadisper(SD_beta_ref_dist$Btotal, Location_type_group_ref))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_ref_dist$Btotal, group =
    ## Location_type_group_ref)
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
(SD_beta_ref_PM <- adonis2(SD_beta_ref_dist$Btotal ~ env[,"Location_type"], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = SD_beta_ref_dist$Btotal ~ env[, "Location_type"], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     3   2.6696 0.30291 2.7521  0.001 ***
    ## Residual 19   6.1435 0.69709                  
    ## Total    22   8.8131 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(SD_beta_ref_PM_pair <- pairwise.adonis(SD_beta_ref_dist$Btotal, env[,"Location_type"], p.adjust.m = "bonferroni", perm = 999))
```

    ## Set of permutations < 'minperm'. Generating entire set.

    ##                     pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted
    ## 1     Stratified vs Mixed  1 1.3139866 4.158150 0.2289964   0.001      0.006
    ## 2     Stratified vs Ocean  1 1.3093089 3.913330 0.2459152   0.001      0.006
    ## 3 Stratified vs Reference  1 0.6261283 1.909357 0.2143091   0.103      0.618
    ## 4          Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.054      0.324
    ## 5      Mixed vs Reference  1 0.5979203 1.966332 0.2193017   0.122      0.732
    ## 6      Ocean vs Reference  1 0.5599217 1.628213 0.2456489   0.139      0.834
    ##   sig
    ## 1   *
    ## 2   *
    ## 3    
    ## 4    
    ## 5    
    ## 6

``` r
# Surveyed sites total
SD_beta_dist <- BAT::beta(surveyed_sites_lake, abund = F)
(SD_beta_BD <- betadisper(SD_beta_dist$Btotal, Location_type_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_dist$Btotal, group = Location_type_group)
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
(SD_beta_PM <- adonis2(SD_beta_dist$Btotal ~ env[surveyed_sites,"Location_type"], permutations = 999, method = "euclidean"))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = SD_beta_dist$Btotal ~ env[surveyed_sites, "Location_type"], permutations = 999, method = "euclidean")
    ##          Df SumOfSqs      R2     F Pr(>F)    
    ## Model     2   2.1121 0.25583 3.266  0.001 ***
    ## Residual 19   6.1435 0.74417                 
    ## Total    21   8.2555 1.00000                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(SD_beta_PM_pair <- pairwise.adonis(SD_beta_dist$Btotal, env[surveyed_sites,"Location_type"], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.3139866 4.158150 0.2289964   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.3093089 3.913330 0.2459152   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.044      0.132

``` r
# Surveyed sites replacement
(SD_beta_rep_BD <- betadisper(SD_beta_dist$Brepl, Location_type_group))
```

    ## Warning in betadisper(SD_beta_dist$Brepl, Location_type_group): some squared
    ## distances are negative and changed to zero

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_dist$Brepl, group = Location_type_group)
    ## 
    ## No. of Positive Eigenvalues: 16
    ## No. of Negative Eigenvalues: 5
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.3052     0.2216     0.2754 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 21 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 0.6075 0.4576 0.4347 0.3451 0.3386 0.2780 0.2574 0.2153

``` r
(SD_beta_rep_AOV <- anova(SD_beta_rep_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Groups     2 0.02557 0.012787  0.4405 0.6501
    ## Residuals 19 0.55158 0.029031

``` r
(SD_beta_rep_THSD <- TukeyHSD(SD_beta_rep_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                         diff        lwr       upr     p adj
    ## Mixed-Ocean      -0.08355805 -0.3173249 0.1502088 0.6419253
    ## Stratified-Ocean -0.02982315 -0.2635900 0.2039437 0.9439116
    ## Stratified-Mixed  0.05373490 -0.1626912 0.2701610 0.8052064

``` r
(SD_beta_rep_PM <- adonis2(SD_beta_dist$Brepl ~ env[surveyed_sites,"Location_type"], permutations = 999, method = "euclidean"))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = SD_beta_dist$Brepl ~ env[surveyed_sites, "Location_type"], permutations = 999, method = "euclidean")
    ##          Df SumOfSqs       R2       F Pr(>F)
    ## Model     2 -0.31917 -0.18783 -1.5022  0.994
    ## Residual 19  2.01847  1.18783               
    ## Total    21  1.69930  1.00000

``` r
(SD_beta_rep_PM_pair <- pairwise.adonis(SD_beta_dist$Brepl, env[surveyed_sites,"Location_type"], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df  SumsOfSqs   F.Model         R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 -0.4165582 -4.325612 -0.4471200   0.994      1.000    
    ## 2 Stratified vs Ocean  1 -0.3967268 -3.415183 -0.3978166   0.997      1.000    
    ## 3      Mixed vs Ocean  1  0.3712364  3.440683  0.2228323   0.003      0.009   *

``` r
# Surveyed sites richness
(SD_beta_ric_BD <- betadisper(SD_beta_dist$Brich, Location_type_group))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_dist$Brich, group = Location_type_group)
    ## 
    ## No. of Positive Eigenvalues: 13
    ## No. of Negative Eigenvalues: 8
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.2127     0.2619     0.2856 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 21 eigenvalues)
    ##   PCoA1   PCoA2   PCoA3   PCoA4   PCoA5   PCoA6   PCoA7   PCoA8 
    ## 2.76461 1.02744 0.16424 0.11046 0.05150 0.04231 0.03669 0.02298

``` r
(SD_beta_ric_AOV <- anova(SD_beta_ric_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Groups     2 0.01851 0.009257  0.2602 0.7736
    ## Residuals 19 0.67594 0.035576

``` r
(SD_beta_ric_THSD <- TukeyHSD(SD_beta_ric_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                        diff        lwr       upr     p adj
    ## Mixed-Ocean      0.04925893 -0.2095210 0.3080388 0.8798615
    ## Stratified-Ocean 0.07289667 -0.1858832 0.3316766 0.7573587
    ## Stratified-Mixed 0.02363774 -0.2159459 0.2632214 0.9660225

``` r
(SD_beta_ric_PM <- adonis2(SD_beta_dist$Brich ~ env[surveyed_sites,"Location_type"], permutations = 999, method = "euclidean"))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = SD_beta_dist$Brich ~ env[surveyed_sites, "Location_type"], permutations = 999, method = "euclidean")
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2   2.0954 0.51574 10.118  0.001 ***
    ## Residual 19   1.9675 0.48426                  
    ## Total    21   4.0629 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(SD_beta_ric_PM_pair <- pairwise.adonis(SD_beta_dist$Brich, env[surveyed_sites,"Location_type"], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df   SumsOfSqs    F.Model          R2 p.value p.adjusted
    ## 1 Stratified vs Mixed  1  1.63476397 14.1534164  0.50272465   0.002      0.006
    ## 2 Stratified vs Ocean  1  1.45313620 14.0922202  0.54009280   0.003      0.009
    ## 3      Mixed vs Ocean  1 -0.02869442 -0.3186705 -0.02728032   1.000      1.000
    ##   sig
    ## 1   *
    ## 2   *
    ## 3

``` r
# Without LCN
SD_beta_wo_LCN_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_LCN,], abund = F)
(SD_beta_wo_LCN_BD <- betadisper(SD_beta_wo_LCN_dist$Btotal, Location_type_group[-8]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_wo_LCN_dist$Btotal, group =
    ## Location_type_group[-8])
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
    ## Mixed-Ocean      0.03896759 -0.05517574 0.1331109 0.5522541
    ## Stratified-Ocean 0.05828298 -0.03586035 0.1524263 0.2794657
    ## Stratified-Mixed 0.01931539 -0.06325377 0.1018846 0.8234561

``` r
(SD_beta_wo_LCN_PM <- adonis2(SD_beta_wo_LCN_dist$Btotal ~ env[surveyed_sites_wo_LCN,"Location_type"], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = SD_beta_wo_LCN_dist$Btotal ~ env[surveyed_sites_wo_LCN, "Location_type"], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2   2.2539 0.28733 3.6286  0.001 ***
    ## Residual 18   5.5904 0.71267                  
    ## Total    20   7.8444 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(SD_beta_wo_LCN_PM_pair <- pairwise.adonis(SD_beta_wo_LCN_dist$Btotal, env[surveyed_sites_wo_LCN,"Location_type"], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.3139866 4.158150 0.2289964   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.3841204 4.397984 0.2856207   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.6396174 2.135322 0.1625633   0.005      0.015   .

``` r
# Without OCO
SD_beta_wo_OCO_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_OCO,], abund = F)
(SD_beta_wo_OCO_BD <- betadisper(SD_beta_wo_OCO_dist$Btotal, Location_type_group[-8]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_wo_OCO_dist$Btotal, group =
    ## Location_type_group[-8])
    ## 
    ## No. of Positive Eigenvalues: 20
    ## No. of Negative Eigenvalues: 0
    ## 
    ## Average distance to median:
    ##      Ocean      Mixed Stratified 
    ##     0.4825     0.5409     0.5691 
    ## 
    ## Eigenvalues for PCoA axes:
    ## (Showing 8 of 20 eigenvalues)
    ##  PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
    ## 1.8371 1.0498 0.5591 0.4808 0.4579 0.4311 0.3607 0.3118

``` r
(SD_beta_wo_OCO_AOV <- anova(SD_beta_wo_OCO_BD))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Groups     2 0.023211 0.0116057  1.5833 0.2326
    ## Residuals 18 0.131937 0.0073298

``` r
(SD_beta_wo_OCO_THSD <- TukeyHSD(SD_beta_wo_OCO_BD))
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = distances ~ group, data = df)
    ## 
    ## $group
    ##                        diff         lwr       upr     p adj
    ## Mixed-Ocean      0.05841586 -0.06614958 0.1829813 0.4701006
    ## Stratified-Ocean 0.08661144 -0.03795400 0.2111769 0.2062131
    ## Stratified-Mixed 0.02819559 -0.08105553 0.1374467 0.7899115

``` r
(SD_beta_wo_OCO_PM <- adonis2(SD_beta_wo_OCO_dist$Btotal ~ env[surveyed_sites_wo_OCO,"Location_type"], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = SD_beta_wo_OCO_dist$Btotal ~ env[surveyed_sites_wo_OCO, "Location_type"], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2   2.0930 0.26865 3.3061  0.001 ***
    ## Residual 18   5.6978 0.73135                  
    ## Total    20   7.7908 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(SD_beta_wo_OCO_PM_pair <- pairwise.adonis(SD_beta_wo_OCO_dist$Btotal, env[surveyed_sites_wo_OCO,"Location_type"], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.3139866 4.15815 0.2289964   0.001      0.003   *
    ## 2 Stratified vs Ocean  1 1.2749118 3.92917 0.2631874   0.001      0.003   *
    ## 3      Mixed vs Ocean  1 0.4889043 1.58069 0.1256441   0.053      0.159

``` r
# Without TLN and HLM
SD_beta_wo_TLN_HLM_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites_wo_TLN_HLM,], abund = F)
(SD_beta_wo_TLN_HLM_BD <- betadisper(SD_beta_wo_TLN_HLM_dist$Btotal, Location_type_group[-c(5,21)]))
```

    ## 
    ##  Homogeneity of multivariate dispersions
    ## 
    ## Call: betadisper(d = SD_beta_wo_TLN_HLM_dist$Btotal, group =
    ## Location_type_group[-c(5, 21)])
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
(SD_beta_wo_TLN_HLM_PM <- adonis2(SD_beta_wo_TLN_HLM_dist$Btotal ~ env[surveyed_sites_wo_TLN_HLM,"Location_type"], permutations = 999))
```

    ## Permutation test for adonis under reduced model
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = SD_beta_wo_TLN_HLM_dist$Btotal ~ env[surveyed_sites_wo_TLN_HLM, "Location_type"], permutations = 999)
    ##          Df SumOfSqs      R2      F Pr(>F)    
    ## Model     2   2.1465 0.28727 3.4259  0.001 ***
    ## Residual 17   5.3256 0.71273                  
    ## Total    19   7.4721 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
(SD_beta_wo_TLN_HLM_PM_pair <- pairwise.adonis(SD_beta_wo_TLN_HLM_dist$Btotal, env[surveyed_sites_wo_TLN_HLM,"Location_type"], p.adjust.m = "bonferroni", perm = 999))
```

    ##                 pairs Df SumsOfSqs  F.Model        R2 p.value p.adjusted sig
    ## 1 Stratified vs Mixed  1 1.4219030 4.731564 0.2827927   0.002      0.006   *
    ## 2 Stratified vs Ocean  1 1.3260066 4.147587 0.2931657   0.004      0.012   .
    ## 3      Mixed vs Ocean  1 0.5079303 1.583987 0.1166069   0.051      0.153

``` r
# Mixed and stratified lakes
SD_beta_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes,], abund = F)
# Ocean sites and mixed lakes
SD_beta_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites,], abund = F)
# Stratified lakes and ocean sites
SD_beta_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites,], abund = F)
# Mixed lakes
SD_beta_M_dist<- BAT::beta(surveyed_sites_lake[mixed_lakes,], abund = F)
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


### Geographic
# Surveyed sites
SD_beta_geo_dist <- BAT::beta(surveyed_sites_lake[surveyed_sites,], abund = F)
# Mixed and stratified lakes
SD_beta_geo_MS_dist<- BAT::beta(surveyed_sites_lake[mixed_stratified_lakes,], abund = F)
# Ocean sites and mixed lakes
SD_beta_geo_OM_dist<- BAT::beta(surveyed_sites_lake[ocean_mixed_sites,], abund = F)
# Stratified lakes and ocean sites
SD_beta_geo_SO_dist<- BAT::beta(surveyed_sites_lake[ocean_stratified_sites,], abund = F)
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
(SO_q <- (abs(group[2,1]))/USE)
```

    ## [1] 0.2899708

``` r
(MS_q <- (abs(group[3,1]))/BSE)
```

    ## [1] 0.7104088

``` r
# Check values here https://www.socscistatistics.com/pvalues/qdistribution.aspx
```

### SD beta dispersions

- The dissimilarity between sites of the same Location type

``` r
# Dispersion within-groups
SD_beta_BD_dist <- SD_beta_BD$distances
SD_beta_BD_dist <- as.data.frame(SD_beta_BD_dist)
SD_beta_BD_dist$X <- row.names(SD_beta_BD_dist)
SD_beta_BD_dist_env <- merge(SD_beta_BD_dist, env[surveyed_sites,], by = "X", sort = F)

# Add levels to Location_type column
SD_beta_BD_dist_env$Location_type <- factor(SD_beta_BD_dist_env$Location_type, levels = c("Ocean", "Mixed", "Stratified"))

# Determine outliers based on dispersion of sites by Location type
outlier_SD_beta_BD_dist_env <- SD_beta_BD_dist_env %>%
  group_by(Location_type) %>%
  mutate(
    Q1 = quantile(SD_beta_BD_dist, 0.25),
    Q3 = quantile(SD_beta_BD_dist, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR,
    is_outlier = SD_beta_BD_dist < lower_bound | SD_beta_BD_dist > upper_bound
  )

outlier_SD_beta_BD_dist_env$is_outlier
```

    ##   25%   25%   25%   25%   25%   25%   25%   25%   25%   25%   25%   25%   25% 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
    ##   25%   25%   25%   25%   25%   25%   25%   25%   25% 
    ## FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

``` r
# Plot dispersion
(SD_beta_BD_dist_env_plot <- ggplot(SD_beta_BD_dist_env, aes(x = Location_type, y = SD_beta_BD_dist, color = Location_type, fill = Location_type)) +
  geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75), aes(color = Location_type, fill = Location_type)) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
    size = 4,
    alpha = 0.8,
    width = 0.1) +
  geom_text_repel(data = SD_beta_BD_dist_env, label = SD_beta_BD_dist_env$X, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme_bw() +
  theme(text = element_text(size = 17), legend.text = element_text(size = 17),
  panel.background = element_rect(fill='transparent'),
  plot.background = element_rect(fill='transparent', color=NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(), axis.text = element_text(size = 17, color = "black"),
  axis.line = element_line(color = "black")) +
  scale_y_continuous(expand = c(0,0.05)) +
  xlab("Location type") +
  ylab("Distance to Centroid") +
  labs(color = "Location_type", tag = "a"))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20dispersion%20boxplot-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_beta_BD_dist_plot.jpg", SD_beta_BD_dist_env_plot, width = 5.4, height = 6, units = "in")
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
# Without LCN
SD_beta_wo_LCN_NMDS <- metaMDS(SD_beta_wo_LCN_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Without TLN and HLM
SD_beta_wo_TLN_HLM_NMDS <- metaMDS(SD_beta_wo_TLN_HLM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed and stratified lakes
SD_beta_MS_NMDS <- metaMDS(SD_beta_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
SD_beta_OM_NMDS <- metaMDS(SD_beta_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
SD_beta_SO_NMDS <- metaMDS(SD_beta_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed lakes
SD_beta_M_NMDS <- metaMDS(SD_beta_M_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

    ## Warning in metaMDS(SD_beta_M_dist$Btotal, try = 1000, parallel = 4, trymax =
    ## 500, : stress is (nearly) zero: you may have insufficient data

``` r
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
### Geographic
# Surveyed sites 
SD_beta_geo_NMDS <- metaMDS(SD_beta_geo_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Mixed and stratified lakes
SD_beta_geo_MS_NMDS <- metaMDS(SD_beta_geo_MS_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Ocean sites and mixed lakes
SD_beta_geo_OM_NMDS <- metaMDS(SD_beta_geo_OM_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
# Stratified lakes and ocean sites
SD_beta_geo_SO_NMDS <- metaMDS(SD_beta_geo_SO_dist$Btotal, try = 1000, parallel = 4, trymax = 500, maxit = 500)
```

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

### SD beta varitation partitioning of env and geo variables

``` r
env_var <- env[surveyed_sites_env,environment]
geo_var <- env[surveyed_sites_env,geography]
SD_beta_varpart <- varpart(SD_beta_env_dist$Btotal, env_var, geo_var)
SD_beta_varpart$part
```

    ## No. of explanatory tables: 2 
    ## Total variation (SS): 7.0233 
    ## No. of observations: 19 
    ## 
    ## Partition table:
    ##                      Df R.squared Adj.R.squared Testable
    ## [a+c] = X1            3   0.34271       0.21125     TRUE
    ## [b+c] = X2            3   0.28407       0.14088     TRUE
    ## [a+b+c] = X1+X2       6   0.50065       0.25098     TRUE
    ## Individual fractions                                    
    ## [a] = X1|X2           3                 0.11009     TRUE
    ## [b] = X2|X1           3                 0.03973     TRUE
    ## [c]                   0                 0.10115    FALSE
    ## [d] = Residuals                         0.74902    FALSE
    ## ---
    ## Use function 'dbrda' to test significance of fractions of interest

``` r
# Open a jpg device
png("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_beta_varpart.jpg", width = 4.5, height = 4.5, units = "in", res = 300, type = "cairo")
# Plot the variation partitioning results
par(mar = c(3, 3, 3, 3) + 1)  # Increase bottom margin if needed
plot(SD_beta_varpart,
     Xnames = c("Env", "Geo"),
     bg = c("mediumpurple", "orange"), alpha = 80,
     digits = 1,
     asp = 1)
# Close the jpg device
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
plot(SD_beta_varpart,
     Xnames = c("Env", "Geo"),
     bg = c("mediumpurple", "orange"), alpha = 80,
     digits = 1,
     asp = 1)
```

![](SD_analyses_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

### SD beta with env and geo correlated variables using envfit

``` r
### Environmental
# Surveyed sites 
# For figure
(SD_beta_env_ef <- envfit(SD_beta_env_NMDS, env[surveyed_sites_env,"S"], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,"Location_type"]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                                  NMDS1     NMDS2     r2 Pr(>r)  
    ## env[surveyed_sites_env, "S"]  0.997810 -0.066189 0.7388  0.025 *
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
    ##                                  NMDS1     NMDS2     r2 Pr(>r)  
    ## env[surveyed_sites_env, "S"]  0.997810 -0.066189 0.7388  0.025 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by Location_type
(SD_beta_env_A_ef <- envfit(SD_beta_env_NMDS, env[surveyed_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,"Location_type"]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.48796 -0.87286 0.0586  0.828  
    ## salinity_median     0.99781 -0.06619 0.7388  0.026 *
    ## oxygen_median       0.79302  0.60919 0.5362  0.807  
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
    ## temperature_median -0.48796 -0.87286 0.0586  1.000  
    ## salinity_median     0.99781 -0.06619 0.7388  0.078 .
    ## oxygen_median       0.79302  0.60919 0.5362  1.000  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
(SD_beta_env_MS_ef <- envfit(SD_beta_env_MS_NMDS, env[mixed_stratified_lakes,environment], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,"Location_type"]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)  
    ## temperature_median -0.56525 -0.82492 0.0778  0.811  
    ## salinity_median     0.99173  0.12832 0.7458  0.044 *
    ## oxygen_median       0.94271  0.33361 0.3552  0.888  
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
    ## temperature_median -0.56525 -0.82492 0.0778  1.000
    ## salinity_median     0.99173  0.12832 0.7458  0.132
    ## oxygen_median       0.94271  0.33361 0.3552  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
(SD_beta_env_OM_ef <- envfit(SD_beta_env_OM_NMDS, env[ocean_mixed_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,"Location_type"]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                       NMDS1    NMDS2     r2 Pr(>r)   
    ## temperature_median  0.10614  0.99435 0.0705  0.707   
    ## salinity_median     0.78732 -0.61655 0.7842  0.005 **
    ## oxygen_median       0.96711  0.25436 0.7188  0.046 * 
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
    ## temperature_median  0.10614  0.99435 0.0705  1.000  
    ## salinity_median     0.78732 -0.61655 0.7842  0.015 *
    ## oxygen_median       0.96711  0.25436 0.7188  0.138  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
(SD_beta_env_SO_ef <- envfit(SD_beta_env_SO_NMDS, env[ocean_stratified_sites_env,environment], permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,"Location_type"]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                         NMDS1      NMDS2     r2 Pr(>r)
    ## temperature_median -0.0001648  1.0000000 0.0783  0.540
    ## salinity_median     0.0007670  1.0000000 0.5187  0.452
    ## oxygen_median       0.0036092  0.9999900 0.7137  0.820
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
    ## temperature_median -0.0001648  1.0000000 0.0783      1
    ## salinity_median     0.0007670  1.0000000 0.5187      1
    ## oxygen_median       0.0036092  0.9999900 0.7137      1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed lakes
(SD_beta_M_ef <- envfit(SD_beta_M_NMDS, env[mixed_lakes,c("temperature_median","salinity_median","oxygen_median","distance_to_ocean_min_m","max_depth","logArea")], permutations = 999, na.rm = TRUE))
```

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1       NMDS2     r2 Pr(>r)  
    ## temperature_median      -0.00008373 -1.00000000 0.2004  0.561  
    ## salinity_median          0.00049790  1.00000000 0.4751  0.199  
    ## oxygen_median            0.00083157  1.00000000 0.6341  0.088 .
    ## distance_to_ocean_min_m -0.00046802 -1.00000000 0.4010  0.297  
    ## max_depth                0.00081129  1.00000000 0.7707  0.035 *
    ## logArea                  0.00097962 -1.00000000 0.6937  0.042 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_M_efp <- p.adjust.envfit(SD_beta_M_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                               NMDS1       NMDS2     r2 Pr(>r)
    ## temperature_median      -0.00008373 -1.00000000 0.2004  1.000
    ## salinity_median          0.00049790  1.00000000 0.4751  1.000
    ## oxygen_median            0.00083157  1.00000000 0.6341  0.528
    ## distance_to_ocean_min_m -0.00046802 -1.00000000 0.4010  1.000
    ## max_depth                0.00081129  1.00000000 0.7707  0.210
    ## logArea                  0.00097962 -1.00000000 0.6937  0.252
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes
(SD_beta_S_ef <- envfit(SD_beta_S_NMDS, env[stratified_lakes,c("temperature_median","salinity_median","oxygen_median","distance_to_ocean_min_m","max_depth","logArea")], permutations = 999, na.rm = TRUE))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median      -0.35881 -0.93341 0.1964  0.609
    ## salinity_median         -0.70036  0.71379 0.5991  0.114
    ## oxygen_median            0.98543 -0.17006 0.5185  0.175
    ## distance_to_ocean_min_m  0.63373 -0.77355 0.2151  0.524
    ## max_depth               -0.49998 -0.86604 0.0383  0.899
    ## logArea                 -0.26513 -0.96421 0.2144  0.539
    ## Permutation: free
    ## Number of permutations: 999

``` r
(SD_beta_S_efp <- p.adjust.envfit(SD_beta_S_ef))
```

    ## Adjustment of significance by bonferroni method

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## temperature_median      -0.35881 -0.93341 0.1964  1.000
    ## salinity_median         -0.70036  0.71379 0.5991  0.684
    ## oxygen_median            0.98543 -0.17006 0.5185  1.000
    ## distance_to_ocean_min_m  0.63373 -0.77355 0.2151  1.000
    ## max_depth               -0.49998 -0.86604 0.0383  1.000
    ## logArea                 -0.26513 -0.96421 0.2144  1.000
    ## Permutation: free
    ## Number of permutations: 999

``` r
### Geographic
# Surveyed sites
# For figure
(SD_beta_geo_ef <- envfit(SD_beta_geo_NMDS, env[surveyed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,"Location_type"]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.99979  0.02028 0.6139  0.122
    ## max_depth               -0.17723 -0.98417 0.1185  0.780
    ## logArea                  0.26562 -0.96408 0.2081  0.129
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
    ## distance_to_ocean_min_m -0.99979  0.02028 0.6139  0.366
    ## max_depth               -0.17723 -0.98417 0.1185  1.000
    ## logArea                  0.26562 -0.96408 0.2081  0.387
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Surveyed sites by Location_type
(SD_beta_geo_A_ef <- envfit(SD_beta_geo_NMDS, env[surveyed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,"Location_type"]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.99979  0.02028 0.6139  0.147
    ## max_depth               -0.17723 -0.98417 0.1185  0.805
    ## logArea                  0.26562 -0.96408 0.2081  0.121
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
    ## distance_to_ocean_min_m -0.99979  0.02028 0.6139  0.441
    ## max_depth               -0.17723 -0.98417 0.1185  1.000
    ## logArea                  0.26562 -0.96408 0.2081  0.363
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Mixed and stratified lakes
(SD_beta_geo_MS_ef <- envfit(SD_beta_geo_MS_NMDS, env[mixed_stratified_lakes,geography], permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,"Location_type"]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.94031  0.34033 0.5171  0.147
    ## max_depth               -0.21483 -0.97665 0.1850  0.651
    ## logArea                 -0.06539  0.99786 0.0118  0.948
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
    ## distance_to_ocean_min_m -0.94031  0.34033 0.5171  0.441
    ## max_depth               -0.21483 -0.97665 0.1850  1.000
    ## logArea                 -0.06539  0.99786 0.0118  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Ocean sites and mixed lakes
(SD_beta_geo_OM_ef <- envfit(SD_beta_geo_OM_NMDS, env[ocean_mixed_sites,geography], permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,"Location_type"]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)  
    ## distance_to_ocean_min_m -0.86176  0.50732 0.3198  0.072 .
    ## max_depth                0.68924  0.72453 0.5689  0.021 *
    ## logArea                  0.84907 -0.52828 0.2993  0.204  
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
    ## distance_to_ocean_min_m -0.86176  0.50732 0.3198  0.216  
    ## max_depth                0.68924  0.72453 0.5689  0.063 .
    ## logArea                  0.84907 -0.52828 0.2993  0.612  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Stratified lakes and ocean sites
(SD_beta_geo_SO_ef <- envfit(SD_beta_geo_SO_NMDS, env[ocean_stratified_sites,geography], permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,"Location_type"]))
```

    ## 
    ## ***VECTORS
    ## 
    ##                            NMDS1    NMDS2     r2 Pr(>r)
    ## distance_to_ocean_min_m -0.86963 -0.49371 0.6793  0.188
    ## max_depth               -0.68534 -0.72822 0.0510  0.965
    ## logArea                  0.52042 -0.85391 0.2187  0.557
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
    ## distance_to_ocean_min_m -0.86963 -0.49371 0.6793  0.564
    ## max_depth               -0.68534 -0.72822 0.0510  1.000
    ## logArea                  0.52042 -0.85391 0.2187  1.000
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

### SD beta Mantel correlation tests

``` r
### Environmental
# Surveyed sites
env_dist_t <- dist(scaled_env[surveyed_sites_env,"temperature_median"], method = "euclidean")
(SD_beta_env_mant_t <- mantel(SD_beta_env_dist$Btotal, env_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_dist$Btotal, ydis = env_dist_t, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, "Location_type"],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.241 
    ##       Significance: 0.265 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.300 0.333 0.353 0.372 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_dist_s <- dist(scaled_env[surveyed_sites_env,"salinity_median"], method = "euclidean")
(SD_beta_env_mant_s <- mantel(SD_beta_env_dist$Btotal, env_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_dist$Btotal, ydis = env_dist_s, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, "Location_type"],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6697 
    ##       Significance: 0.007 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.566 0.608 0.639 0.658 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_dist_o <- dist(scaled_env[surveyed_sites_env,"oxygen_median"], method = "euclidean")
(SD_beta_env_mant_o <- mantel(SD_beta_env_dist$Btotal, env_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites_env,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_dist$Btotal, ydis = env_dist_o, method = "spearman",      permutations = 999, strata = env[surveyed_sites_env, "Location_type"],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r:  0.47 
    ##       Significance: 0.486 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.602 0.638 0.662 0.699 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_mant_pv <- rbind(SD_beta_env_mant_t$signif, SD_beta_env_mant_s$signif, SD_beta_env_mant_o$signif)
SD_beta_env_mant_pv <- SD_beta_env_mant_pv[,1]
(SD_beta_env_mant_pv <- p.adjust(SD_beta_env_mant_pv, method = "bonferroni"))
```

    ## [1] 0.795 0.021 1.000

``` r
# Mixed and stratified lakes
env_MS_dist_t <- dist(scaled_env[mixed_stratified_lakes,"temperature_median"], method = "euclidean")
(SD_beta_env_MS_mant_t <- mantel(SD_beta_env_MS_dist$Btotal, env_MS_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_MS_dist$Btotal, ydis = env_MS_dist_t,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2368 
    ##       Significance: 0.325 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.329 0.370 0.402 0.417 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_MS_dist_s <- dist(scaled_env[mixed_stratified_lakes,"salinity_median"], method = "euclidean")
(SD_beta_env_MS_mant_s <- mantel(SD_beta_env_MS_dist$Btotal, env_MS_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_MS_dist$Btotal, ydis = env_MS_dist_s,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6676 
    ##       Significance: 0.012 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.550 0.609 0.636 0.670 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_MS_dist_o <- dist(scaled_env[mixed_stratified_lakes,"oxygen_median"], method = "euclidean")
(SD_beta_env_MS_mant_o <- mantel(SD_beta_env_MS_dist$Btotal, env_MS_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_MS_dist$Btotal, ydis = env_MS_dist_o,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2979 
    ##       Significance: 0.602 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.483 0.521 0.564 0.603 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_MS_mant_pv <- rbind(SD_beta_env_MS_mant_t$signif, SD_beta_env_MS_mant_s$signif, SD_beta_env_MS_mant_o$signif)
SD_beta_env_MS_mant_pv <- SD_beta_env_MS_mant_pv[,1]
(SD_beta_env_MS_mant_pv <- p.adjust(SD_beta_env_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 0.975 0.036 1.000

``` r
# Ocean sites and mixed lakes
env_OM_dist_t <- dist(scaled_env[ocean_mixed_sites_env,"temperature_median"], method = "euclidean")
(SD_beta_env_OM_mant_t <- mantel(SD_beta_env_OM_dist$Btotal, env_OM_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_OM_dist$Btotal, ydis = env_OM_dist_t,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_env,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1214 
    ##       Significance: 0.749 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.122 0.168 0.228 0.405 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_OM_dist_s <- dist(scaled_env[ocean_mixed_sites_env,"salinity_median"], method = "euclidean")
(SD_beta_env_OM_mant_s <- mantel(SD_beta_env_OM_dist$Btotal, env_OM_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_OM_dist$Btotal, ydis = env_OM_dist_s,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_env,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4258 
    ##       Significance: 0.048 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.366 0.423 0.490 0.561 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_OM_dist_o <- dist(scaled_env[ocean_mixed_sites_env,"oxygen_median"], method = "euclidean")
(SD_beta_env_OM_mant_o <- mantel(SD_beta_env_OM_dist$Btotal, env_OM_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites_env,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_OM_dist$Btotal, ydis = env_OM_dist_o,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites_env,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5319 
    ##       Significance: 0.031 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.339 0.462 0.550 0.637 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_OM_mant_pv <- rbind(SD_beta_env_OM_mant_t$signif, SD_beta_env_OM_mant_s$signif, SD_beta_env_OM_mant_o$signif)
SD_beta_env_OM_mant_pv <- SD_beta_env_OM_mant_pv[,1]
(SD_beta_env_OM_mant_pv <- p.adjust(SD_beta_env_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.144 0.093

``` r
# Stratified lakes and ocean sites
env_SO_dist_t <- dist(scaled_env[ocean_stratified_sites_env,"temperature_median"], method = "euclidean")
(SD_beta_env_SO_mant_t <- mantel(SD_beta_env_SO_dist$Btotal, env_SO_dist_t, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_SO_dist$Btotal, ydis = env_SO_dist_t,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites_env,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.008188 
    ##       Significance: 0.223 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0223 0.0434 0.0633 0.0767 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_SO_dist_s <- dist(scaled_env[ocean_stratified_sites_env,"salinity_median"], method = "euclidean")
(SD_beta_env_SO_mant_s <- mantel(SD_beta_env_SO_dist$Btotal, env_SO_dist_s, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_SO_dist$Btotal, ydis = env_SO_dist_s,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites_env,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5528 
    ##       Significance: 0.013 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.456 0.497 0.521 0.563 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_SO_dist_o <- dist(scaled_env[ocean_stratified_sites_env,"oxygen_median"], method = "euclidean")
(SD_beta_env_SO_mant_o <- mantel(SD_beta_env_SO_dist$Btotal, env_SO_dist_o, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites_env,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_env_SO_dist$Btotal, ydis = env_SO_dist_o,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites_env,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6567 
    ##       Significance: 0.353 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.706 0.729 0.749 0.769 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_env_SO_mant_pv <- rbind(SD_beta_env_SO_mant_t$signif, SD_beta_env_SO_mant_s$signif, SD_beta_env_SO_mant_o$signif)
SD_beta_env_SO_mant_pv <- SD_beta_env_SO_mant_pv[,1]
(SD_beta_env_SO_mant_pv <- p.adjust(SD_beta_env_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 0.669 0.039 1.000

``` r
# Mixed lakes
env_M_dist_t <- dist(scaled_env[mixed_lakes,"temperature_median"], method = "euclidean")
(SD_beta_M_mant_t <- mantel(SD_beta_M_dist$Btotal, env_M_dist_t, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_M_dist$Btotal, ydis = env_M_dist_t, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1352 
    ##       Significance: 0.74 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.240 0.318 0.397 0.522 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_M_dist_s <- dist(scaled_env[mixed_lakes,"salinity_median"], method = "euclidean")
(SD_beta_M_mant_s <- mantel(SD_beta_M_dist$Btotal, env_M_dist_s, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_M_dist$Btotal, ydis = env_M_dist_s, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1877 
    ##       Significance: 0.157 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.253 0.319 0.371 0.453 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_M_dist_o <- dist(scaled_env[mixed_lakes,"oxygen_median"], method = "euclidean")
(SD_beta_M_mant_o <- mantel(SD_beta_M_dist$Btotal, env_M_dist_o, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_M_dist$Btotal, ydis = env_M_dist_o, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2781 
    ##       Significance: 0.084 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.233 0.377 0.452 0.581 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_M_dist_dm <- dist(scaled_env[mixed_lakes,"distance_to_ocean_min_m"], method = "euclidean")
(SD_beta_M_mant_dm <- mantel(SD_beta_M_dist$Btotal, geo_M_dist_dm, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_M_dist$Btotal, ydis = geo_M_dist_dm, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2877 
    ##       Significance: 0.04 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.212 0.265 0.343 0.728 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_M_dist_md <- dist(scaled_env[mixed_lakes,"max_depth"], method = "euclidean")
(SD_beta_M_mant_md <- mantel(SD_beta_M_dist$Btotal, geo_M_dist_md, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_M_dist$Btotal, ydis = geo_M_dist_md, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.6829 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.290 0.451 0.531 0.582 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_M_dist_la <- dist(scaled_env[mixed_lakes,"logArea"], method = "euclidean")
(SD_beta_M_mant_la <- mantel(SD_beta_M_dist$Btotal, geo_M_dist_la, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_M_dist$Btotal, ydis = geo_M_dist_la, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.416 
    ##       Significance: 0.037 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.236 0.345 0.460 0.531 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_M_mant_pv <- rbind(SD_beta_M_mant_t$signif, SD_beta_M_mant_s$signif, SD_beta_M_mant_o$signif, SD_beta_M_mant_dm$signif, SD_beta_M_mant_md$signif, SD_beta_M_mant_la$signif)
SD_beta_M_mant_pv <- SD_beta_M_mant_pv[,1]
(SD_beta_M_mant_pv <- p.adjust(SD_beta_M_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.942 0.504 0.240 0.006 0.222

``` r
# Stratified lakes
env_S_dist_t <- dist(scaled_env[stratified_lakes,"temperature_median"], method = "euclidean")
(SD_beta_S_mant_t <- mantel(SD_beta_S_dist$Btotal, env_S_dist_t, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_S_dist$Btotal, ydis = env_S_dist_t, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.1252 
    ##       Significance: 0.783 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.207 0.264 0.314 0.409 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_S_dist_s <- dist(scaled_env[stratified_lakes,"salinity_median"], method = "euclidean")
(SD_beta_S_mant_s <- mantel(SD_beta_S_dist$Btotal, env_S_dist_s, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_S_dist$Btotal, ydis = env_S_dist_s, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4001 
    ##       Significance: 0.02 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.241 0.317 0.379 0.425 
    ## Permutation: free
    ## Number of permutations: 999

``` r
env_S_dist_o <- dist(scaled_env[stratified_lakes,"oxygen_median"], method = "euclidean")
(SD_beta_S_mant_o <- mantel(SD_beta_S_dist$Btotal, env_S_dist_o, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_S_dist$Btotal, ydis = env_S_dist_o, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.307 
    ##       Significance: 0.048 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.235 0.301 0.369 0.435 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_S_dist_dm <- dist(scaled_env[stratified_lakes,"distance_to_ocean_min_m"], method = "euclidean")
(SD_beta_S_mant_dm <- mantel(SD_beta_S_dist$Btotal, geo_S_dist_dm, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_S_dist$Btotal, ydis = geo_S_dist_dm, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.06837 
    ##       Significance: 0.364 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.254 0.317 0.393 0.457 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_S_dist_md <- dist(scaled_env[stratified_lakes,"max_depth"], method = "euclidean")
(SD_beta_S_mant_md <- mantel(SD_beta_S_dist$Btotal, geo_S_dist_md, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_S_dist$Btotal, ydis = geo_S_dist_md, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.06552 
    ##       Significance: 0.593 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.250 0.325 0.419 0.462 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_S_dist_la <- dist(scaled_env[stratified_lakes,"logArea"], method = "euclidean")
(SD_beta_S_mant_la <- mantel(SD_beta_S_dist$Btotal, geo_S_dist_la, method = "spearman", permutations = 999, na.rm = TRUE))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_S_dist$Btotal, ydis = geo_S_dist_la, method = "spearman",      permutations = 999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.08347 
    ##       Significance: 0.286 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.189 0.273 0.402 0.503 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_S_mant_pv <- rbind(SD_beta_S_mant_t$signif, SD_beta_S_mant_s$signif, SD_beta_S_mant_o$signif, SD_beta_S_mant_dm$signif, SD_beta_S_mant_md$signif, SD_beta_S_mant_la$signif)
SD_beta_S_mant_pv <- SD_beta_S_mant_pv[,1]
(SD_beta_S_mant_pv <- p.adjust(SD_beta_S_mant_pv, method = "bonferroni"))
```

    ## [1] 1.000 0.120 0.288 1.000 1.000 1.000

``` r
### Geographic
# Surveyed sites
geo_dist_dm <- dist(scaled_env[surveyed_sites,"distance_to_ocean_min_m"], method = "euclidean")
(SD_beta_geo_mant_dm <- mantel(SD_beta_geo_dist$Btotal, geo_dist_dm, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_dist$Btotal, ydis = geo_dist_dm, method = "spearman",      permutations = 999, strata = env[surveyed_sites, "Location_type"],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.4904 
    ##       Significance: 0.084 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.480 0.510 0.537 0.559 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_dist_md <- dist(scaled_env[surveyed_sites,"max_depth"], method = "euclidean")
(SD_beta_geo_mant_md <- mantel(SD_beta_geo_dist$Btotal, geo_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_dist$Btotal, ydis = geo_dist_md, method = "spearman",      permutations = 999, strata = env[surveyed_sites, "Location_type"],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.06903 
    ##       Significance: 0.649 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.196 0.230 0.259 0.278 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_dist_la <- dist(scaled_env[surveyed_sites,"logArea"], method = "euclidean")
(SD_beta_geo_mant_la <- mantel(SD_beta_geo_dist$Btotal, geo_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[surveyed_sites,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_dist$Btotal, ydis = geo_dist_la, method = "spearman",      permutations = 999, strata = env[surveyed_sites, "Location_type"],      na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.04299 
    ##       Significance: 0.393 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0842 0.0971 0.1111 0.1236 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_mant_pv <- rbind(SD_beta_geo_mant_dm$signif, SD_beta_geo_mant_md$signif, SD_beta_geo_mant_la$signif)
SD_beta_geo_mant_pv <- SD_beta_geo_mant_pv[,1]
(SD_beta_geo_mant_pv <- p.adjust(SD_beta_geo_mant_pv, method = "bonferroni"))
```

    ## [1] 0.252 1.000 1.000

``` r
# Mixed and stratified lakes 
geo_MS_dist_dm <- dist(scaled_env[mixed_stratified_lakes,"distance_to_ocean_min_m"], method = "euclidean")
(SD_beta_geo_MS_mant_dm <- mantel(SD_beta_geo_MS_dist$Btotal, geo_MS_dist_dm, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_dm,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.3388 
    ##       Significance: 0.115 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.355 0.399 0.432 0.471 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_MS_dist_md <- dist(scaled_env[mixed_stratified_lakes,"max_depth"], method = "euclidean")
(SD_beta_geo_MS_mant_md <- mantel(SD_beta_geo_MS_dist$Btotal, geo_MS_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_md,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.007317 
    ##       Significance: 0.759 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.223 0.275 0.307 0.341 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_MS_dist_la <- dist(scaled_env[mixed_stratified_lakes,"logArea"], method = "euclidean")
(SD_beta_geo_MS_mant_la <- mantel(SD_beta_geo_MS_dist$Btotal, geo_MS_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[mixed_stratified_lakes,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_MS_dist$Btotal, ydis = geo_MS_dist_la,      method = "spearman", permutations = 999, strata = env[mixed_stratified_lakes,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.04736 
    ##       Significance: 0.773 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0707 0.0948 0.1183 0.1451 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_MS_mant_pv <- rbind(SD_beta_geo_MS_mant_dm$signif, SD_beta_geo_MS_mant_md$signif, SD_beta_geo_MS_mant_la$signif)
SD_beta_geo_MS_mant_pv <- SD_beta_geo_MS_mant_pv[,1]
(SD_beta_geo_MS_mant_pv <- p.adjust(SD_beta_geo_MS_mant_pv, method = "bonferroni"))
```

    ## [1] 0.345 1.000 1.000

``` r
# Ocean sites and mixed lakes
geo_OM_dist_dm <- dist(scaled_env[ocean_mixed_sites,"distance_to_ocean_min_m"], method = "euclidean")
(SD_beta_geo_OM_mant_dm <- mantel(SD_beta_geo_OM_dist$Btotal, geo_OM_dist_dm, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_dm,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2052 
    ##       Significance: 0.103 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.206 0.237 0.269 0.387 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_OM_dist_md <- dist(scaled_env[ocean_mixed_sites,"max_depth"], method = "euclidean")
(SD_beta_geo_OM_mant_md <- mantel(SD_beta_geo_OM_dist$Btotal, geo_OM_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_md,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2349 
    ##       Significance: 0.055 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.174 0.243 0.308 0.370 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_OM_dist_la <- dist(scaled_env[ocean_mixed_sites,"logArea"], method = "euclidean")
(SD_beta_geo_OM_mant_la <- mantel(SD_beta_geo_OM_dist$Btotal, geo_OM_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_mixed_sites,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_OM_dist$Btotal, ydis = geo_OM_dist_la,      method = "spearman", permutations = 999, strata = env[ocean_mixed_sites,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2007 
    ##       Significance: 0.113 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.207 0.256 0.287 0.321 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_OM_mant_pv <- rbind( SD_beta_geo_OM_mant_dm$signif, SD_beta_geo_OM_mant_md$signif, SD_beta_geo_OM_mant_la$signif)
SD_beta_geo_OM_mant_pv <- SD_beta_geo_OM_mant_pv[,1]
(SD_beta_geo_OM_mant_pv <- p.adjust(SD_beta_geo_OM_mant_pv, method = "bonferroni"))
```

    ## [1] 0.309 0.165 0.339

``` r
# Stratified lakes and ocean sites
geo_SO_dist_dm <- dist(scaled_env[ocean_stratified_sites,"distance_to_ocean_min_m"], method = "euclidean")
(SD_beta_geo_SO_mant_dm <- mantel(SD_beta_geo_SO_dist$Btotal, geo_SO_dist_dm, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_SO_dist$Btotal, ydis = geo_SO_dist_dm,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.5512 
    ##       Significance: 0.138 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.569 0.598 0.619 0.651 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_SO_dist_md <- dist(scaled_env[ocean_stratified_sites,"max_depth"], method = "euclidean")
(SD_beta_geo_SO_mant_md <- mantel(SD_beta_geo_SO_dist$Btotal, geo_SO_dist_md, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_SO_dist$Btotal, ydis = geo_SO_dist_md,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: -0.007031 
    ##       Significance: 0.618 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.147 0.175 0.196 0.248 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
geo_SO_dist_la <- dist(scaled_env[ocean_stratified_sites,"logArea"], method = "euclidean")
(SD_beta_geo_SO_mant_la <- mantel(SD_beta_geo_SO_dist$Btotal, geo_SO_dist_la, method = "spearman", permutations = 999, na.rm = TRUE, strata = env[ocean_stratified_sites,"Location_type"]))
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## mantel(xdis = SD_beta_geo_SO_dist$Btotal, ydis = geo_SO_dist_la,      method = "spearman", permutations = 999, strata = env[ocean_stratified_sites,          "Location_type"], na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.2461 
    ##       Significance: 0.224 
    ## 
    ## Upper quantiles of permutations (null model):
    ##   90%   95% 97.5%   99% 
    ## 0.276 0.293 0.308 0.320 
    ## Blocks:  strata 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Adjust p-values
SD_beta_geo_SO_mant_pv <- rbind( SD_beta_geo_SO_mant_dm$signif, SD_beta_geo_SO_mant_md$signif, SD_beta_geo_SO_mant_la$signif)
SD_beta_geo_SO_mant_pv <- SD_beta_geo_SO_mant_pv[,1]
(SD_beta_geo_SO_mant_pv <- p.adjust(SD_beta_geo_SO_mant_pv, method = "bonferroni"))
```

    ## [1] 0.414 1.000 0.672

### SD beta NMDS ordination plots

``` r
# SD beta total NMDS scores
SD_beta_NMDS_data.scores <- as.data.frame(scores(SD_beta_NMDS))
SD_beta_NMDS_data.scores$Location_type <- env[surveyed_sites,"Location_type"]
SD_beta_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
SD_beta_NMDS_data.scores$Location_type <- factor(SD_beta_NMDS_data.scores$Location_type, levels = c("Ocean", "Mixed", "Stratified"))

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
(SD_beta_ef_plot <- ggplot(data = SD_beta_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Location_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Location_type), type = "t", level = 0.95) +
  geom_point(data = SD_beta_NMDS_data.scores, aes(color = Location_type, fill = Location_type), size = 2, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
     data = SD_beta_env_ef_coord_cont, linewidth =1, alpha = 0.2, color = env_cont) +
  geom_text(data = SD_beta_env_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
    color = env_cont, label = row.names(SD_beta_env_ef_coord_cont), size = 3) +
  geom_text_repel(data = SD_beta_NMDS_data.scores, label = SD_beta_NMDS_data.scores$Lakes, 
      size = 3, force = 10, point.padding = 0, max.overlaps = 30, nudge_x = -0.05) +
  theme(text = element_text(size = 12), 
    legend.position = "bottom", 
    legend.title = element_text(size = 10), 
    legend.text = element_text(size = 10), 
    axis.title = element_text(size = 10),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.border = element_rect(fill = NA), 
    axis.text = element_text(color = "black", size = 10), 
    legend.key = element_blank()) +
  annotate("text", x = -0.5, y = 0.55, size = 3.5,
     label = paste("Stress: ", round(SD_beta_NMDS$stress, digits = 2))) +
  labs(colour = "Location type:  ", fill = "Location type:  ", tag = "a") +
  coord_fixed() +
  guides(color = "none", fill = "none") +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81)))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20NMDS%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_beta_NMDS.jpg", SD_beta_ef_plot, width = 3.25, height = 3.3, units = "in")


# Plot NMDS ordination of SD beta total with CI = 0.90
SD_beta_ef_plot_CI90 <- ggplot(data = SD_beta_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Location_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Location_type), type = "t", level = 0.90) +
  geom_point(data = SD_beta_NMDS_data.scores, aes(color = Location_type, fill = Location_type), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
     data = SD_beta_env_ef_coord_cont, linewidth =1, alpha = 0.2, color = env_cont) +
  geom_text(data = SD_beta_env_ef_coord_cont, aes(x = NMDS1, y = NMDS2),
    color = env_cont, label = row.names(SD_beta_env_ef_coord_cont), size = 5) +
  geom_text_repel(data = SD_beta_NMDS_data.scores, label = SD_beta_NMDS_data.scores$Lakes, 
      size = 5, force = 10, point.padding = 0, max.overlaps = 30, nudge_x = -0.05) +
  theme(text = element_text(size = 17), 
    legend.position = "bottom", 
    legend.title = element_text(size = 17), 
    legend.text = element_text(size = 17), 
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.border = element_rect(fill = NA), 
    axis.text = element_text(color = "black"), 
    legend.key = element_blank()) +
  annotate("text", x = -0.6, y = 0.5, size = 5,
     label = paste("Stress: ", round(SD_beta_NMDS$stress, digits = 2))) +
  labs(colour = "Location type:  ", fill = "Location type:  ", tag = "b") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81))
(SD_beta_ef_plot_CI90 <- SD_beta_ef_plot_CI90 + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20NMDS%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_beta_NMDS_CI90.jpg", SD_beta_ef_plot_CI90, width = 5.4, height = 6, units = "in")


# SD beta total NMDS scores without LCN
SD_beta_wo_LCN_NMDS_data.scores <- as.data.frame(scores(SD_beta_wo_LCN_NMDS))
SD_beta_wo_LCN_NMDS_data.scores$Location_type <- env[surveyed_sites_wo_LCN,"Location_type"]
SD_beta_wo_LCN_NMDS_data.scores$Lakes <- env[surveyed_sites_wo_LCN,1]
SD_beta_wo_LCN_NMDS_data.scores$Location_type <- factor(SD_beta_wo_LCN_NMDS_data.scores$Location_type, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of SD beta total without LCN with CI = 0.95
SD_beta_wo_LCN_plot <- ggplot(data = SD_beta_wo_LCN_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Location_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Location_type), type = "t", level = 0.95) +
  geom_point(data = SD_beta_wo_LCN_NMDS_data.scores, aes(color = Location_type, fill = Location_type), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = SD_beta_wo_LCN_NMDS_data.scores, label = SD_beta_wo_LCN_NMDS_data.scores$Lakes, 
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
  annotate("text", x = -0.6, y = 0.5, size = 5,
     label = paste("Stress: ", round(SD_beta_wo_LCN_NMDS$stress, digits = 2))) +
  labs(colour = "Location type:  ", fill = "Location type:  ", tag = "a") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81))
(SD_beta_wo_LCN_plot <- SD_beta_wo_LCN_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20NMDS%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_beta_NMDS_wo_LCN.jpg", SD_beta_wo_LCN_plot, width = 4.88, height = 6, units = "in")


# SD beta total NMDS scores without TLN and HLM
SD_beta_wo_TLN_HLM_NMDS_data.scores <- as.data.frame(scores(SD_beta_wo_TLN_HLM_NMDS))
SD_beta_wo_TLN_HLM_NMDS_data.scores$Location_type <- env[surveyed_sites_wo_TLN_HLM,"Location_type"]
SD_beta_wo_TLN_HLM_NMDS_data.scores$Lakes <- env[surveyed_sites_wo_TLN_HLM,1]
SD_beta_wo_TLN_HLM_NMDS_data.scores$Location_type <- factor(SD_beta_wo_TLN_HLM_NMDS_data.scores$Location_type, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of SD beta total without TLN and HLM with CI = 0.95
SD_beta_wo_TLN_HLM_plot <- ggplot(data = SD_beta_wo_TLN_HLM_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Location_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Location_type), type = "t", level = 0.95) +
  geom_point(data = SD_beta_wo_TLN_HLM_NMDS_data.scores, aes(color = Location_type, fill = Location_type), size = 4, alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = SD_beta_wo_TLN_HLM_NMDS_data.scores, label = SD_beta_wo_TLN_HLM_NMDS_data.scores$Lakes, 
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
  annotate("text", x = -0.6, y = 0.5, size = 5,
     label = paste("Stress: ", round(SD_beta_wo_TLN_HLM_NMDS$stress, digits = 2))) +
  labs(colour = "Location type:  ", fill = "Location type:  ", tag = "b") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-0.81,0.81))
(SD_beta_wo_TLN_HLM_plot <- SD_beta_wo_TLN_HLM_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20NMDS%20plots-4.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_beta_NMDS_wo_TLN_HLM.jpg", SD_beta_wo_TLN_HLM_plot, width = 4.88, height = 6, units = "in")


# SD beta replacement NMDS scores
SD_beta_rep_NMDS_data.scores <- as.data.frame(scores(SD_beta_rep_NMDS))
SD_beta_rep_NMDS_data.scores$Location_type <- env[surveyed_sites,"Location_type"]
SD_beta_rep_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
SD_beta_rep_NMDS_data.scores$Location_type <- factor(SD_beta_rep_NMDS_data.scores$Location_type, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of SD beta replacement with CI = 0.95
SD_beta_rep_plot <- ggplot(data = SD_beta_rep_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Location_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Location_type), type = "t", level = 0.95) +
  geom_point(data = SD_beta_rep_NMDS_data.scores, aes(color = Location_type, fill = Location_type), size = 4, alpha = 1) + 
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
  labs(colour = "Location type:  ", fill = "Location type:  ", tag = "a") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
(SD_beta_rep_plot <- SD_beta_rep_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20NMDS%20plots-5.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_beta_rep_NMDS.jpg", SD_beta_rep_plot, width = 4.88, height = 6, units = "in")


# SD beta richness NMDS scores
SD_beta_ric_NMDS_data.scores <- as.data.frame(scores(SD_beta_ric_NMDS))
SD_beta_ric_NMDS_data.scores$Location_type <- env[surveyed_sites,"Location_type"]
SD_beta_ric_NMDS_data.scores$Lakes <- env[surveyed_sites,1]
SD_beta_ric_NMDS_data.scores$Location_type <- factor(SD_beta_ric_NMDS_data.scores$Location_type, levels = c("Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of SD beta richness with CI = 0.95
SD_beta_ric_plot <- ggplot(data = SD_beta_ric_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Location_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Location_type), type = "t", level = 0.95) +
  geom_point(data = SD_beta_ric_NMDS_data.scores, aes(color = Location_type, fill = Location_type), size = 4, alpha = 1) + 
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
  labs(colour = "Location type:  ", fill = "Location type:  ", tag = "b") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
(SD_beta_ric_plot <- SD_beta_ric_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20NMDS%20plots-6.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_beta_ric_NMDS.jpg", SD_beta_ric_plot, width = 4.88, height = 6, units = "in")


# SD beta ref NMDS scores
SD_beta_ref_NMDS_data.scores <- as.data.frame(scores(SD_beta_ref_NMDS))
SD_beta_ref_NMDS_data.scores$Location_type <- env[,"Location_type"]
SD_beta_ref_NMDS_data.scores$Lakes <- env[,1]
SD_beta_ref_NMDS_data.scores$Location_type <- factor(SD_beta_ref_NMDS_data.scores$Location_type, levels = c("Reference", "Ocean", "Mixed", "Stratified"))

# Plot NMDS ordination of SD beta ref with CI = 0.95
SD_beta_ref_plot <- ggplot(data = SD_beta_ref_NMDS_data.scores, aes(x = NMDS1, y = NMDS2, color = Location_type)) + 
  stat_ellipse(aes(x = NMDS1, y = NMDS2, color = Location_type), type = "t", level = 0.95) +
  geom_point(data = SD_beta_ref_NMDS_data.scores, aes(color = Location_type, fill = Location_type), size = 4, alpha = 1) + 
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
  labs(colour = "Location type:  ", fill = "Location type:  ", tag = "a") +
  coord_fixed() +
  scale_x_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05)) +
  scale_y_continuous(breaks = c(-0.5,0,0.5), limits = c(-1.05,1.05))
(SD_beta_ref_plot <- SD_beta_ref_plot + guides(color = guide_legend(override.aes = list(label = ""))))
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_path()`).

![](SD_analyses_files/figure-gfm/SD%20beta%20NMDS%20plots-7.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_beta_ref_NMDS.jpg", SD_beta_ref_plot, width = 6.26, height = 6, units = "in")
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
png("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_beta_dendrogram.jpg", width = 6.5, height = 3, units = "in", res = 300, type = "cairo")

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
png("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_beta_dg_silhouette.jpg", width = 6.5, height = 3, units = "in", res = 300, type = "cairo")

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
png("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/SD_beta_dg_partitioning.jpg", width = 6.5, height = 3, units = "in", res = 300, type = "cairo")

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
  group_by(Location_type) %>%
  mutate(
  Q1 = quantile(NMDS1, 0.25),
  Q3 = quantile(NMDS1, 0.75),
  IQR = Q3 - Q1,
  lower_bound = Q1 - 1.5 * IQR,
  upper_bound = Q3 + 1.5 * IQR,
  is_outlier = NMDS1 < lower_bound | NMDS1 > upper_bound
  )

outlier_SD_beta_NMDS1 <- as.data.frame(outlier_SD_beta_NMDS1)

row.names(outlier_SD_beta_NMDS1) <- outlier_SD_beta_NMDS1$Lakes

# Create the plot
(outlier_SD_beta_NMDS1_plot <- ggplot(outlier_SD_beta_NMDS1, aes(x = Location_type, y = NMDS1, fill = Location_type)) +
  geom_violin(alpha = 0.9, draw_quantiles = c(0.25, 0.5, 0.75), aes(fill = Location_type)) +
  geom_jitter(aes(color = is_outlier), width = 0.1, size = 4, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "purple")) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = outlier_SD_beta_NMDS1, label = outlier_SD_beta_NMDS1$Lakes, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme(text = element_text(size = 16),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text = element_text(color = "black", size = 16)) +
  scale_y_continuous(expand = c(0,0.05)) +
  guides(fill = "none") + 
  labs(y = "NMDS1 Distances", x = "Location type", color = "Outlier:", tag = "a"))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20outliers-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/outlier_SD_beta_NMDS1.jpg", outlier_SD_beta_NMDS1_plot, width = 4.88, height = 6, units = "in")


# Determine NMDS2 outliers
ordered_SD_beta_NMDS2_data.scores <- SD_beta_NMDS_data.scores[order(SD_beta_NMDS_data.scores$NMDS2), ]

outlier_SD_beta_NMDS2 <- ordered_SD_beta_NMDS2_data.scores %>%
  group_by(Location_type) %>%
  mutate(
  Q1 = quantile(NMDS2, 0.25),
  Q3 = quantile(NMDS2, 0.75),
  IQR = Q3 - Q1,
  lower_bound = Q1 - 1.5 * IQR,
  upper_bound = Q3 + 1.5 * IQR,
  is_outlier = NMDS2 < lower_bound | NMDS2 > upper_bound
  )

outlier_SD_beta_NMDS2 <- as.data.frame(outlier_SD_beta_NMDS2)

row.names(outlier_SD_beta_NMDS2) <- outlier_SD_beta_NMDS2$Lakes

# Create the plot
(outlier_SD_beta_NMDS2_plot <- ggplot(outlier_SD_beta_NMDS2, aes(x = Location_type, y = NMDS2, fill = Location_type)) +
  geom_violin(alpha = 0.9, draw_quantiles = c(0.25, 0.5, 0.75), aes(fill = Location_type)) +
  geom_jitter(aes(color = is_outlier), width = 0.1, size = 4, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "purple")) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(data = outlier_SD_beta_NMDS2, label = outlier_SD_beta_NMDS2$Lakes, size = 5, point.padding = 7, max.overlaps = 30, nudge_x = -0.01) +
  theme(text = element_text(size = 16),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text = element_text(color = "black", size = 16)) +
  scale_y_continuous(expand = c(0,0.05)) +
  guides(fill = "none") + 
  labs(y = "NMDS2 Distances", x = "Location type", color = "Outlier:", tag = "b"))
```

![](SD_analyses_files/figure-gfm/SD%20beta%20outliers-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/SD/outlier_SD_beta_NMDS2.jpg", outlier_SD_beta_NMDS2_plot, width = 4.88, height = 6, units = "in")
```

### Package and version info

``` r
sessionInfo()
```

    ## R version 4.4.3 (2025-02-28)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS 26.0.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Chicago
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] MASS_7.3-65          car_3.1-3            carData_3.0-5       
    ##  [4] emmeans_1.11.0       ggvenn_0.1.10        pairwiseAdonis_0.4.1
    ##  [7] cluster_2.1.8.1      BAT_2.9.6            caret_7.0-1         
    ## [10] ggrepel_0.9.6        ggplot2_4.0.0        picante_1.8.2       
    ## [13] nlme_3.1-168         vegan_2.6-10         lattice_0.22-7      
    ## [16] permute_0.9-7        tidyr_1.3.1          phytools_2.4-4      
    ## [19] maps_3.4.2.1         ape_5.8-1            reshape2_1.4.4      
    ## [22] stringr_1.5.1        dplyr_1.1.4          knitr_1.50          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       magrittr_2.0.3         
    ##   [4] TH.data_1.1-3           estimability_1.5.1      farver_2.1.2           
    ##   [7] rmarkdown_2.29          ragg_1.4.0              vctrs_0.6.5            
    ##  [10] base64enc_0.1-3         terra_1.8-42            polspline_1.1.25       
    ##  [13] htmltools_0.5.8.1       progress_1.2.3          DEoptim_2.2-8          
    ##  [16] Formula_1.2-5           pROC_1.18.5             parallelly_1.43.0      
    ##  [19] pracma_2.4.4            KernSmooth_2.23-26      htmlwidgets_1.6.4      
    ##  [22] plyr_1.8.9              sandwich_3.1-1          palmerpenguins_0.1.1   
    ##  [25] zoo_1.8-14              lubridate_1.9.4         igraph_2.1.4           
    ##  [28] lifecycle_1.0.4         iterators_1.0.14        pkgconfig_2.0.3        
    ##  [31] Matrix_1.7-3            R6_2.6.1                fastmap_1.2.0          
    ##  [34] future_1.40.0           magic_1.6-1             digest_0.6.37          
    ##  [37] numDeriv_2016.8-1.1     colorspace_2.1-1        Hmisc_5.2-3            
    ##  [40] textshaping_1.0.0       pdist_1.2.1             labeling_0.4.3         
    ##  [43] clusterGeneration_1.3.8 timechange_0.3.0        abind_1.4-8            
    ##  [46] mgcv_1.9-3              compiler_4.4.3          proxy_0.4-27           
    ##  [49] withr_3.0.2             doParallel_1.0.17       backports_1.5.0        
    ##  [52] htmlTable_2.4.3         S7_0.2.0                optimParallel_1.0-2    
    ##  [55] quantreg_6.1            lava_1.8.1              scatterplot3d_0.3-44   
    ##  [58] ModelMetrics_1.2.2.2    tools_4.4.3             foreign_0.8-90         
    ##  [61] future.apply_1.11.3     nnet_7.3-20             glue_1.8.0             
    ##  [64] quadprog_1.5-8          checkmate_2.3.2         generics_0.1.3         
    ##  [67] recipes_1.3.0           gtable_0.3.6            class_7.3-23           
    ##  [70] data.table_1.17.0       hms_1.1.3               foreach_1.5.2          
    ##  [73] pillar_1.10.2           splines_4.4.3           survival_3.8-3         
    ##  [76] SparseM_1.84-2          ks_1.14.3               tidyselect_1.2.1       
    ##  [79] rms_8.0-0               gridExtra_2.3           stats4_4.4.3           
    ##  [82] xfun_0.52               expm_1.0-0              hardhat_1.4.1          
    ##  [85] timeDate_4041.110       proto_1.0.0             stringi_1.8.7          
    ##  [88] yaml_2.3.10             evaluate_1.0.3          codetools_0.2-20       
    ##  [91] tibble_3.2.1            cli_3.6.4               rpart_4.1.24           
    ##  [94] nls2_0.3-4              xtable_1.8-4            geometry_0.5.2         
    ##  [97] systemfonts_1.2.2       Rcpp_1.0.14             globals_0.17.0         
    ## [100] coda_0.19-4.1           fastcluster_1.2.6       MatrixModels_0.5-4     
    ## [103] gower_1.0.2             prettyunits_1.2.0       mclust_6.1.1           
    ## [106] listenv_0.9.1           phangorn_2.12.1         mvtnorm_1.3-3          
    ## [109] ipred_0.9-15            scales_1.4.0            prodlim_2024.06.25     
    ## [112] e1071_1.7-16            purrr_1.0.4             crayon_1.5.3           
    ## [115] combinat_0.0-8          rlang_1.1.6             fastmatch_1.1-6        
    ## [118] multcomp_1.4-28         mnormt_2.1.1            hypervolume_3.1.5

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
