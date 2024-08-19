Taxonomic diversity analyses
================

#### R Markdown

## Load packages

``` r
# Load the knitr package if not already loaded
library(knitr)
library(betapart)
library(ggplot2)
library(viridis)
```

    ## Loading required package: viridisLite

``` r
library(ggrepel)
library(ggvenn)
```

    ## Loading required package: dplyr

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## Loading required package: grid

``` r
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
strat_fish <- stratified_fish[,-c(1,10)]
mix_fish <- mixed_fish[,-c(1,10)]
oc_fish <- ocean_fish[,-c(1,8)]

# Define your custom colors
custom_colors <- c("Reference" = "black", "Ocean" = "#EE6363", "Mixed" = "#87CEFA", "Stratified" = "#6E8B3D")
```

## Plot species richness for all sites with more than 1 species

``` r
#Plot species richness
SR_env$Stratification <- factor(SR_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

SR_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = reorder(X, row_sum, decreasing = T), y = row_sum, color = Stratification, fill = Stratification)) + 
  geom_bar(stat = 'identity',
    alpha = 1) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.title = element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Site", y="Species Richness", colour = "Site type", fill = "Site type")
SR_plot
```

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_plot.png", plot = SR_plot, width = 6, height = 4, units = "in")

subset_SR_env <- subset(SR_env, Stratification == "Ocean")
ordered_subset_SR_env <- subset_SR_env[order(subset_SR_env$row_sum), ]
Q1 <- ordered_subset_SR_env[2,2]
Q3 <- ordered_subset_SR_env[5,2]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] -45.5

``` r
Q3 + 1.5*IQR1
```

    ## [1] 142.5

``` r
ordered_subset_SR_env$row_sum
```

    ## [1] 24 25 52 54 72 87

``` r
subset_SR_env <- subset(SR_env, Stratification == "Stratified")
ordered_subset_SR_env <- subset_SR_env[order(subset_SR_env$row_sum), ]
Q1 <- ordered_subset_SR_env[2,2]
Q3 <- ordered_subset_SR_env[5,2]
IQR1 <- Q3 - Q1
Q1 - 1.5*IQR1
```

    ## [1] 0

``` r
Q3 + 1.5*IQR1
```

    ## [1] 8

``` r
ordered_subset_SR_env$row_sum
```

    ## [1]  3  3  3  5  5  6 17 18

## Plot species richness venn diagram

``` r
site_type_pres <- strat_presabs_lake[c(24:26),]
site_type_pres <- t(site_type_pres)
site_type_pres <- site_type_pres[which(rowSums(site_type_pres) > 0),]

# Convert numerical values to logical
site_type_pres_logical <- as.data.frame(site_type_pres > 0)

# Check the structure of the transformed data
str(site_type_pres_logical)
```

    ## 'data.frame':    246 obs. of  3 variables:
    ##  $ Ocean sites     : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ Mixed lakes     : logi  TRUE FALSE TRUE FALSE TRUE FALSE ...
    ##  $ Stratified lakes: logi  FALSE FALSE FALSE FALSE FALSE FALSE ...

``` r
# Generate a color palette with three colors
# palette <- viridis_pal(alpha = 1, begin = 0.45, end = 0.75, option = "G")(3)
# palette

# Plot the colors
# barplot(rep(1, 3), col = palette, border = NA, axes = FALSE, main = "Viridis Palette")

# Create the Venn diagram
venn_plot <- ggvenn(site_type_pres_logical,
                    show_percentage = F,
                    fill_color = c("#EE6363", "#87CEFA", "#6E8B3D"),
                    fill_alpha = 0.7,
                    stroke_alpha = 0,
                    stroke_size = 0.5, 
                    set_name_size = 5,
                    text_size = 5)

venn_plot
```

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/venn_plot.png", plot = venn_plot, width = 6, height = 4, units = "in")
```

## Subsample species richness venn diagram

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

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/venn_plot_s10.png", plot = venn_plot, width = 6, height = 4, units = "in")
```

## SR complete island Linear Models

``` r
SR_model <- lm(log(row_sum) ~ surface_area_m2 + distance_to_ocean_mean_m + tidal_lag_time_minutes + max_depth + logArea, data = SR_env[c(1:2,4,5,12,17,21:22),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ surface_area_m2 + distance_to_ocean_mean_m + 
    ##     tidal_lag_time_minutes + max_depth + logArea, data = SR_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##     BCM     CLM     GLK     HLM     NLK     OTM     SLN     TLN 
    ## -0.2672 -0.1610  0.3630  1.0760 -0.1021 -0.8898 -0.2476  0.2285 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)               4.843e+00  1.057e+01   0.458    0.692
    ## surface_area_m2           3.686e-05  5.200e-05   0.709    0.552
    ## distance_to_ocean_mean_m -6.733e-03  6.842e-03  -0.984    0.429
    ## tidal_lag_time_minutes    1.374e-02  1.673e-02   0.821    0.498
    ## max_depth                -1.763e-02  7.675e-02  -0.230    0.840
    ## logArea                  -4.471e-01  1.244e+00  -0.359    0.754
    ## 
    ## Residual standard error: 1.073 on 2 degrees of freedom
    ## Multiple R-squared:  0.3921, Adjusted R-squared:  -1.128 
    ## F-statistic: 0.258 on 5 and 2 DF,  p-value: 0.9037

``` r
p_values <- summary(SR_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                          Df  Sum Sq Mean Sq F value Pr(>F)
    ## surface_area_m2           1 0.00262 0.00262  0.0023 0.9663
    ## distance_to_ocean_mean_m  1 0.60912 0.60912  0.5291 0.5426
    ## tidal_lag_time_minutes    1 0.54721 0.54721  0.4753 0.5618
    ## max_depth                 1 0.17773 0.17773  0.1544 0.7323
    ## logArea                   1 0.14864 0.14864  0.1291 0.7537
    ## Residuals                 2 2.30254 1.15127

``` r
# Get R-squared value
r_squared <- summary(SR_model)$r.squared

# Get number of observations
n <- nrow(SR_env[c(1:2,4,5,12,17,21:22),])

# Get number of predictor variables (excluding intercept)
k <- length(coef(SR_model)) - 1  # Subtract 1 for the intercept

# Calculate Cohen's f^2
f_squared <- r_squared / (1 - r_squared) * ((n - k - 1) / k)

# Print the result
print(f_squared)
```

    ## [1] 0.2580315

``` r
SR_model <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median + pH_median, data = SR_env[c(1:2,4,5,12,17,21:22),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = SR_env[c(1:2, 4, 5, 12, 
    ##     17, 21:22), ])
    ## 
    ## Residuals:
    ##     BCM     CLM     GLK     HLM     NLK     OTM     SLN     TLN 
    ##  0.7994 -1.1237  3.4120 -0.5078 -1.0232 -4.1916 -0.6347  3.2696 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)        -101.3096    67.4448  -1.502    0.230
    ## temperature_median   -1.2977     1.2020  -1.080    0.359
    ## salinity_median       0.2196     0.8161   0.269    0.805
    ## oxygen_median        -3.7251     2.9862  -1.247    0.301
    ## pH_median            20.7304    12.3223   1.682    0.191
    ## 
    ## Residual standard error: 3.808 on 3 degrees of freedom
    ## Multiple R-squared:  0.8424, Adjusted R-squared:  0.6322 
    ## F-statistic: 4.007 on 4 and 3 DF,  p-value: 0.1417

``` r
p_values <- summary(SR_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                    Df  Sum Sq Mean Sq F value Pr(>F)  
    ## temperature_median  1   4.710   4.710  0.3248 0.6086  
    ## salinity_median     1 176.532 176.532 12.1716 0.0398 *
    ## oxygen_median       1  10.197  10.197  0.7031 0.4633  
    ## pH_median           1  41.050  41.050  2.8303 0.1911  
    ## Residuals           3  43.511  14.504                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get R-squared value
r_squared <- summary(SR_model)$r.squared

# Get number of observations
n <- nrow(SR_env[c(1:2,4,5,12,17,21:22),])

# Get number of predictor variables (excluding intercept)
k <- length(coef(SR_model)) - 1  # Subtract 1 for the intercept

# Calculate Cohen's f^2
f_squared <- r_squared / (1 - r_squared) * ((n - k - 1) / k)

# Print the result
print(f_squared)
```

    ## [1] 4.007449

``` r
SR_model <- lm(row_sum ~ salinity_median, data = SR_env[c(1:2,4,5,12,17,21:22),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ salinity_median, data = SR_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -5.7935 -2.0573  0.5463  2.1924  5.4615 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)     -32.5070    11.8988  -2.732   0.0341 *
    ## salinity_median   1.4570     0.4303   3.386   0.0147 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 3.975 on 6 degrees of freedom
    ## Multiple R-squared:  0.6565, Adjusted R-squared:  0.5992 
    ## F-statistic: 11.46 on 1 and 6 DF,  p-value: 0.01475

``` r
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                 Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## salinity_median  1 181.181 181.181  11.465 0.01475 *
    ## Residuals        6  94.819  15.803                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(SR_env[c(1:2,4,5,12,17,21:22),2], SR_env[c(1:2,4,5,12,17,21:22),7])
```

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# Island biogeography common relationships
SR_model <- lm(log(row_sum) ~ logArea, data = SR_env[c(1:2,4,5,12,17,21:22),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ logArea, data = SR_env[c(1:2, 4, 
    ##     5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.6556 -0.6542 -0.1445  0.2983  1.1360 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)
    ## (Intercept) 1.7470332  2.5976226   0.673    0.526
    ## logArea     0.0006553  0.2516849   0.003    0.998
    ## 
    ## Residual standard error: 0.7946 on 6 degrees of freedom
    ## Multiple R-squared:  1.13e-06,   Adjusted R-squared:  -0.1667 
    ## F-statistic: 6.78e-06 on 1 and 6 DF,  p-value: 0.998

``` r
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##           Df Sum Sq Mean Sq F value Pr(>F)
    ## logArea    1 0.0000 0.00000       0  0.998
    ## Residuals  6 3.7879 0.63131

``` r
SR_model <- lm(log(row_sum) ~ log(distance_to_ocean_mean_m), data = SR_env[c(1:2,4,5,12,17,21:22),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ log(distance_to_ocean_mean_m), data = SR_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.7213 -0.4688 -0.1221  0.2440  1.2929 
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                     4.8788     3.9251   1.243    0.260
    ## log(distance_to_ocean_mean_m)  -0.5652     0.7083  -0.798    0.455
    ## 
    ## Residual standard error: 0.7555 on 6 degrees of freedom
    ## Multiple R-squared:  0.09595,    Adjusted R-squared:  -0.05472 
    ## F-statistic: 0.6368 on 1 and 6 DF,  p-value: 0.4553

``` r
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                               Df Sum Sq Mean Sq F value Pr(>F)
    ## log(distance_to_ocean_mean_m)  1 0.3635 0.36346  0.6368 0.4553
    ## Residuals                      6 3.4244 0.57073

``` r
SR_model <- lm(log(row_sum) ~ log(max_depth), data = SR_env[c(1:2,4,5,12,17,21:22),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ log(max_depth), data = SR_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.7527 -0.6017 -0.1620  0.3290  1.2110 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)      2.2403     1.8108   1.237    0.262
    ## log(max_depth)  -0.1565     0.5756  -0.272    0.795
    ## 
    ## Residual standard error: 0.7897 on 6 degrees of freedom
    ## Multiple R-squared:  0.01218,    Adjusted R-squared:  -0.1525 
    ## F-statistic: 0.07396 on 1 and 6 DF,  p-value: 0.7948

``` r
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                Df Sum Sq Mean Sq F value Pr(>F)
    ## log(max_depth)  1 0.0461 0.04612   0.074 0.7948
    ## Residuals       6 3.7417 0.62362

## SR habitat island Linear Models

``` r
SR_model <- lm(log(row_sum) ~ surface_area_m2 + distance_to_ocean_mean_m + tidal_lag_time_minutes + max_depth + logArea, data = SR_env[c(3,6,8,10,13:15,23),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ surface_area_m2 + distance_to_ocean_mean_m + 
    ##     tidal_lag_time_minutes + max_depth + logArea, data = SR_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##      FLK      HLO      MLN      NLN      NLU      OLO      ULN 
    ##  0.31825 -0.09661  0.22839 -0.16544 -0.08578 -0.40657  0.20775 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)               3.356e+00  4.508e+00   0.744    0.593
    ## surface_area_m2           2.975e-05  4.766e-05   0.624    0.645
    ## distance_to_ocean_mean_m -3.977e-03  1.018e-02  -0.391    0.763
    ## tidal_lag_time_minutes    1.951e-03  1.106e-02   0.176    0.889
    ## max_depth                 1.485e-02  6.045e-02   0.246    0.847
    ## logArea                  -3.815e-03  6.407e-01  -0.006    0.996
    ## 
    ## Residual standard error: 0.6372 on 1 degrees of freedom
    ##   (1 observation deleted due to missingness)
    ## Multiple R-squared:  0.7839, Adjusted R-squared:  -0.2966 
    ## F-statistic: 0.7255 on 5 and 1 DF,  p-value: 0.7068

``` r
p_values <- summary(SR_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                          Df  Sum Sq Mean Sq F value Pr(>F)
    ## surface_area_m2           1 1.35382 1.35382  3.3348 0.3189
    ## distance_to_ocean_mean_m  1 0.08058 0.08058  0.1985 0.7332
    ## tidal_lag_time_minutes    1 0.00002 0.00002  0.0001 0.9954
    ## max_depth                 1 0.03819 0.03819  0.0941 0.8105
    ## logArea                   1 0.00001 0.00001  0.0000 0.9962
    ## Residuals                 1 0.40597 0.40597

``` r
# Get R-squared value
r_squared <- summary(SR_model)$r.squared

# Get number of observations
n <- nrow(SR_env[c(3,6,8,10,13:15,23),])

# Get number of predictor variables (excluding intercept)
k <- length(coef(SR_model)) - 1  # Subtract 1 for the intercept

# Calculate Cohen's f^2
f_squared <- r_squared / (1 - r_squared) * ((n - k - 1) / k)

# Print the result
print(f_squared)
```

    ## [1] 1.450989

``` r
SR_model <- lm(row_sum ~ temperature_median + salinity_median + oxygen_median + pH_median, data = SR_env[c(3,6,8,10,13:15,23),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = SR_env[c(3, 6, 8, 10, 13:15, 
    ##     23), ])
    ## 
    ## Residuals:
    ##     FLK     HLO     LLN     MLN     NLN     NLU     OLO     ULN 
    ##   1.078  12.935   6.696 -17.610  12.044  -1.455 -17.313   3.624 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)        -680.961    892.125  -0.763   0.5008  
    ## temperature_median  -16.416     12.332  -1.331   0.2753  
    ## salinity_median      -6.374     21.275  -0.300   0.7840  
    ## oxygen_median       -17.272     24.035  -0.719   0.5243  
    ## pH_median           194.565     77.339   2.516   0.0865 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 18.11 on 3 degrees of freedom
    ## Multiple R-squared:  0.8364, Adjusted R-squared:  0.6184 
    ## F-statistic: 3.835 on 4 and 3 DF,  p-value: 0.1491

``` r
p_values <- summary(SR_model)$coefficients[, "Pr(>|t|)"]
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##                    Df  Sum Sq Mean Sq F value Pr(>F)  
    ## temperature_median  1 1209.80 1209.80  3.6904 0.1505  
    ## salinity_median     1  769.05  769.05  2.3459 0.2231  
    ## oxygen_median       1  975.73  975.73  2.9764 0.1830  
    ## pH_median           1 2074.82 2074.82  6.3290 0.0865 .
    ## Residuals           3  983.48  327.83                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Get R-squared value
r_squared <- summary(SR_model)$r.squared

# Get number of observations
n <- nrow(SR_env[c(3,6,8,10,13:15,23),])

# Get number of predictor variables (excluding intercept)
k <- length(coef(SR_model)) - 1  # Subtract 1 for the intercept

# Calculate Cohen's f^2
f_squared <- r_squared / (1 - r_squared) * ((n - k - 1) / k)

# Print the result
print(f_squared)
```

    ## [1] 3.835407

``` r
SR_model <- lm(row_sum ~ pH_median, data = SR_env[c(3,6,8,10,13:15,23),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = row_sum ~ pH_median, data = SR_env[c(3, 6, 8, 10, 
    ##     13:15, 23), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -24.192  -7.504   3.375   8.406  23.869 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept) -1157.43     320.03  -3.617  0.01114 * 
    ## pH_median     154.62      40.84   3.786  0.00912 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 17.2 on 6 degrees of freedom
    ## Multiple R-squared:  0.7049, Adjusted R-squared:  0.6557 
    ## F-statistic: 14.33 on 1 and 6 DF,  p-value: 0.00912

``` r
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: row_sum
    ##           Df Sum Sq Mean Sq F value  Pr(>F)   
    ## pH_median  1 4238.4  4238.4  14.331 0.00912 **
    ## Residuals  6 1774.5   295.7                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(SR_env[c(3,6,8,10,13:15,23),2], SR_env[c(3,6,8,10,13:15,23),11])
```

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
# Island biogeography common relationships
SR_model <- lm(log(row_sum) ~ logArea, data = SR_env[c(3,6,8,10,13:15,23),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ logArea, data = SR_env[c(3, 6, 8, 
    ##     10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.82627  0.03313  0.10586  0.15986  0.20691 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)  
    ## (Intercept)  0.48846    0.91576   0.533   0.6129  
    ## logArea      0.34660    0.09387   3.692   0.0102 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.3682 on 6 degrees of freedom
    ## Multiple R-squared:  0.6944, Adjusted R-squared:  0.6435 
    ## F-statistic: 13.63 on 1 and 6 DF,  p-value: 0.01018

``` r
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##           Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## logArea    1 1.84847 1.84847  13.634 0.01018 *
    ## Residuals  6 0.81344 0.13557                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
SR_model <- lm(log(row_sum) ~ log(distance_to_ocean_mean_m), data = SR_env[c(3,6,8,10,13:15,23),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ log(distance_to_ocean_mean_m), data = SR_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.9931 -0.4814  0.2609  0.3497  0.7527 
    ## 
    ## Coefficients:
    ##                               Estimate Std. Error t value Pr(>|t|)
    ## (Intercept)                     2.4274     3.0181   0.804    0.452
    ## log(distance_to_ocean_mean_m)   0.2924     0.6249   0.468    0.656
    ## 
    ## Residual standard error: 0.6542 on 6 degrees of freedom
    ## Multiple R-squared:  0.03521,    Adjusted R-squared:  -0.1256 
    ## F-statistic: 0.219 on 1 and 6 DF,  p-value: 0.6563

``` r
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                               Df  Sum Sq Mean Sq F value Pr(>F)
    ## log(distance_to_ocean_mean_m)  1 0.09373 0.09373   0.219 0.6563
    ## Residuals                      6 2.56818 0.42803

``` r
SR_model <- lm(log(row_sum) ~ log(max_depth), data = SR_env[c(3,6,8,10,13:15,23),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ log(max_depth), data = SR_env[c(3, 
    ##     6, 8, 10, 13:15, 23), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.8488 -0.1681  0.0071  0.2292  0.6257 
    ## 
    ## Coefficients:
    ##                Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)      2.5956     0.5944   4.367  0.00473 **
    ## log(max_depth)   0.5201     0.2382   2.184  0.07167 . 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4972 on 6 degrees of freedom
    ## Multiple R-squared:  0.4429, Adjusted R-squared:   0.35 
    ## F-statistic:  4.77 on 1 and 6 DF,  p-value: 0.07167

``` r
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                Df Sum Sq Mean Sq F value  Pr(>F)  
    ## log(max_depth)  1 1.1789 1.17891  4.7697 0.07167 .
    ## Residuals       6 1.4830 0.24717                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## SR vs log area, distance to the ocean, max depth

``` r
SR_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = logArea, y = log(row_sum), color = Stratification, fill = Stratification)) +
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = SR_env[-20,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.title = element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Log Area (m"^"2"~")", y= "Log SRic", colour = "Site type", fill = "Site type")
# SR_plot <- SR_plot + geom_smooth(mapping = aes(x = logArea, y = row_sum), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
SR_plot <- SR_plot + guides(color = guide_legend(override.aes = list(label = "")))
SR_plot
```

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_plot_area.png", SR_plot, width = 6, height = 4, units = "in")


SR_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = log(distance_to_ocean_mean_m), y = log(row_sum), color = Stratification, fill = Stratification)) +
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = SR_env[-20,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.title = element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Log Isolation", y= "Log SRic", colour = "Site type", fill = "Site type")
SR_plot <- SR_plot + guides(color = guide_legend(override.aes = list(label = "")))
SR_plot
```

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_plot_dist.png", SR_plot, width = 6, height = 4, units = "in")


SR_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = log(max_depth), y = log(row_sum), color = Stratification, fill = Stratification)) +
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = SR_env[-20,1], size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.title = element_text(size = 16), legend.text = element_text(size = 16),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Log Age", y= "Log SRic", colour = "Site type", fill = "Site type")
SR_plot <- SR_plot + guides(color = guide_legend(override.aes = list(label = "")))
SR_plot
```

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_plot_depth.png", SR_plot, width = 6, height = 4, units = "in")
```

## SR vs Oxygen Concentration

``` r
SR_oxygen_median_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = oxygen_median, y = row_sum, color = Stratification, fill = Stratification)) +
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 18),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Oxygen Concentration (mg/L)", y= "SRic", colour = "Site type", fill = "Site type") +
  ylim(c(0,110))
SR_oxygen_median_plot <- SR_oxygen_median_plot + geom_smooth(mapping = aes(x = oxygen_median, y = row_sum), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
SR_oxygen_median_plot
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 3 rows containing non-finite outside the scale range
    ## (`stat_smooth()`).

    ## Warning: Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_smooth()`).

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_oxygen_median_plot.png", SR_oxygen_median_plot, width = 6, height = 4, units = "in")
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 3 rows containing non-finite outside the scale range
    ## (`stat_smooth()`).

    ## Warning: Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values or values outside the scale range
    ## (`geom_smooth()`).

## SR vs tidal efficiency

``` r
SR_tidal_efficiency_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = tidal_efficiency, y = row_sum, color = Stratification, fill = Stratification)) +
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(text = element_text(size = 24),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Tidal Efficiency", y= "SRic", colour = "Site type", fill = "Site type") +
  ylim(c(0,110))
SR_tidal_efficiency_plot <- SR_tidal_efficiency_plot + geom_smooth(mapping = aes(x = tidal_efficiency, y = row_sum), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
SR_tidal_efficiency_plot
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1 row containing non-finite outside the scale range
    ## (`stat_smooth()`).

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_tidal_efficiency_plot.png", SR_tidal_efficiency_plot, width = 6, height = 4, units = "in")
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1 row containing non-finite outside the scale range (`stat_smooth()`).
    ## Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

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
    ## [1] parallel  grid      stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] tidyr_1.3.0       phytools_2.0-3    maps_3.4.1.1      ape_5.7-1        
    ##  [5] reshape2_1.4.4    stringr_1.5.1     ggvenn_0.1.10     dplyr_1.1.4      
    ##  [9] ggrepel_0.9.4     viridis_0.6.4     viridisLite_0.4.2 ggplot2_3.5.1    
    ## [13] betapart_1.6      knitr_1.45       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.0        farver_2.1.1            optimParallel_1.0-2    
    ##  [4] fastmap_1.1.1           combinat_0.0-8          digest_0.6.33          
    ##  [7] lifecycle_1.0.4         cluster_2.1.6           magrittr_2.0.3         
    ## [10] compiler_4.3.1          rlang_1.1.2             doSNOW_1.0.20          
    ## [13] tools_4.3.1             igraph_1.5.1            utf8_1.2.4             
    ## [16] yaml_2.3.7              geometry_0.4.7          phangorn_2.11.1        
    ## [19] clusterGeneration_1.3.8 labeling_0.4.3          mnormt_2.1.1           
    ## [22] scatterplot3d_0.3-44    plyr_1.8.9              abind_1.4-5            
    ## [25] expm_0.999-8            picante_1.8.2           purrr_1.0.2            
    ## [28] withr_2.5.2             numDeriv_2016.8-1.1     itertools_0.1-3        
    ## [31] fansi_1.0.6             colorspace_2.1-0        scales_1.3.0           
    ## [34] iterators_1.0.14        MASS_7.3-60             cli_3.6.1              
    ## [37] rmarkdown_2.25          vegan_2.6-4             ragg_1.2.6             
    ## [40] generics_0.1.3          rcdd_1.5-2              rstudioapi_0.15.0      
    ## [43] magic_1.6-1             splines_4.3.1           vctrs_0.6.5            
    ## [46] Matrix_1.6-4            minpack.lm_1.2-4        systemfonts_1.0.5      
    ## [49] foreach_1.5.2           snow_0.4-4              glue_1.6.2             
    ## [52] codetools_0.2-19        stringi_1.8.2           gtable_0.3.4           
    ## [55] quadprog_1.5-8          munsell_0.5.0           tibble_3.2.1           
    ## [58] pillar_1.9.0            htmltools_0.5.7         R6_2.5.1               
    ## [61] textshaping_0.3.7       doParallel_1.0.17       evaluate_0.23          
    ## [64] lattice_0.22-5          highr_0.10              Rcpp_1.0.11            
    ## [67] fastmatch_1.1-4         coda_0.19-4             gridExtra_2.3          
    ## [70] nlme_3.1-164            permute_0.9-7           mgcv_1.9-0             
    ## [73] xfun_0.41               pkgconfig_2.0.3

## Beta Diversity of SR

``` r
# Do the following to order lakes by the circulation pattern, then delete env data
####https://pedrohbraga.github.io/CommunityPhylogenetics-Workshop/CommunityPhylogenetics-Workshop.html#between-assemblage-phylogenetic-structure####
presabs_env <- merge(pres_abs_by_lake, env, by = 'X', sort = F)
presabs_env <- presabs_env[order(presabs_env$Stratification),]
row.names(presabs_env) <- presabs_env$X
presabs_env_ordered <- presabs_env[,-c(1, 1725:1770)]

SR_betap_jac <- beta.pair(presabs_env_ordered[-15,], index.family = "jaccard")

# Turnover matrix:
SR.jac.turn <- SR_betap_jac$beta.jtu %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("Sites") + ylab("Sites") + labs(fill = "Turnover") + # edit axis and legend titles
  scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels

# Nestedness matrix:
SR.jac.nest <- SR_betap_jac$beta.jne %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("Sites") + ylab("Sites") + labs(fill = "Nestedness") + # edit axis and legend titles
    scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels

# Nestedness matrix:
SR.jac <- SR_betap_jac$beta.jac %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("Sites") + ylab("Sites") + labs(fill = "Beta Diversity") + # edit axis and legend titles
    scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels
SR.jac

# plot both heatmaps next to each other
SR_beta_div_plot <- gridExtra::grid.arrange(SR.jac.turn, SR.jac.nest, ncol = 2)
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_beta_div_plot.png", SR_beta_div_plot, width = 6, height = 4.5)
```

## Determine what families of fish are most present in the lakes

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
ggsave("family.png", family, width = 6, height = 4.5)

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
#ggsave("order.png", order)

detach("package:rfishbase")
```
