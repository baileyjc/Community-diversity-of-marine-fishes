Trait plots
================

### Load packages

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
library(viridis)
```

    ## Loading required package: viridisLite

``` r
library(dunn.test)
library(car)
```

    ## Loading required package: carData

    ## 
    ## Attaching package: 'car'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     recode

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
# Define your custom colors
custom_shades <- c("R" = "black", "O" = "#EE6363", "M" = "#87CEFA", "S" = "#6E8B3D")

# Reference use Location_type_at_weighted
# Surveyed sites use Location_type_st_weighted
```

## Location trait data tests

``` r
# Create function
run_tests <- function(trait, data) {
    # Ensure the trait is numeric
    if (!is.numeric(data[[trait]])) {
        return(c(ANOVA_p = NA, Shapiro_p = NA, Levene_p = NA))}  # Return NA if not numeric
    formula <- as.formula(paste(trait, "~ Location_type"))  # Dynamically create formula
    anova_model <- aov(formula, data = data)  # Run ANOVA
    anova_p <- summary(anova_model)[[1]][["Pr(>F)"]][1]  # Extract p-value
    # Run Shapiro-Wilk test on residuals (normality assumption)
    shapiro_p <- tryCatch(
        shapiro.test(residuals(anova_model))$p.value,
        error = function(e) NA)  # Return NA if Shapiro fails
    # Run Levene's test for homogeneity of variance
    levene_p <- tryCatch(
        leveneTest(formula, data = data)$`Pr(>F)`[1],
        error = function(e) NA)  # Return NA if Levene fails
    return(c(ANOVA_p = anova_p, Shapiro_p = shapiro_p, Levene_p = levene_p))
}

# Run tests on numeric traits
traits <- colnames(Locations_at)[1:15]  # Extract trait names
results <- sapply(traits, run_tests, data = Locations_at)

# Convert to a tidy dataframe
results_df <- as.data.frame(t(results))
results_df$Trait <- rownames(results_df)
rownames(results_df) <- NULL

# Print results
print(results_df)
```

    ##         ANOVA_p    Shapiro_p    Levene_p            Trait
    ## 1            NA           NA          NA       BodyShapeI
    ## 2            NA           NA          NA      DemersPelag
    ## 3            NA           NA          NA OperculumPresent
    ## 4  1.562054e-01 8.280539e-36 0.079330384      MaxLengthTL
    ## 5  5.573640e-07 1.213909e-10 0.000287748            Troph
    ## 6  3.962927e-01 3.995155e-44 0.992085783         DepthMin
    ## 7  4.231113e-01 2.497339e-34 0.034674309         DepthMax
    ## 8  1.372934e-02 1.949571e-15 0.223943688      TempPrefMin
    ## 9  1.085826e-01 6.375706e-34 0.450032256      TempPrefMax
    ## 10           NA           NA          NA      FeedingPath
    ## 11           NA           NA          NA        RepGuild1
    ## 12           NA           NA          NA        RepGuild2
    ## 13           NA           NA          NA     ParentalCare
    ## 14           NA           NA          NA        WaterPref
    ## 15 1.846265e-02 3.016593e-19 0.001352179 DorsalSpinesMean

``` r
# Identify traits that are numeric
numerical <- c("MaxLengthTL", "Troph", "DepthMin", "DepthMax", "TempPrefMin", "TempPrefMax", "DorsalSpinesMean")

## MaxLengthTL
# Perform Dunn's test
dunn_result <- dunn.test(Locations_at$MaxLengthTL, Locations_at$Location_type, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 1.7554, df = 2, p-value = 0.42
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |      Mixed      Ocean
    ## ---------+----------------------
    ##    Ocean |   0.967279
    ##          |     0.5001
    ##          |
    ## Stratifi |   1.101793   0.569996
    ##          |     0.4058     0.8530
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 1.755429
    ## 
    ## $Z
    ## [1] 0.9672800 1.1017938 0.5699965
    ## 
    ## $P
    ## [1] 0.1667020 0.1352757 0.2843400
    ## 
    ## $P.adjusted
    ## [1] 0.5001061 0.4058270 0.8530201
    ## 
    ## $comparisons
    ## [1] "Mixed - Ocean"      "Mixed - Stratified" "Ocean - Stratified"

``` r
## Troph
# Perform Dunn's test
dunn_result <- dunn.test(Locations_at$Troph, Locations_at$Location_type, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 29.7428, df = 2, p-value = 0
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |      Mixed      Ocean
    ## ---------+----------------------
    ##    Ocean |   4.364086
    ##          |    0.0000*
    ##          |
    ## Stratifi |  -2.191762  -4.438696
    ##          |     0.0426    0.0000*
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 29.7428
    ## 
    ## $Z
    ## [1]  4.364086 -2.191762 -4.438696
    ## 
    ## $P
    ## [1] 6.382766e-06 1.419834e-02 4.525275e-06
    ## 
    ## $P.adjusted
    ## [1] 1.914830e-05 4.259501e-02 1.357583e-05
    ## 
    ## $comparisons
    ## [1] "Mixed - Ocean"      "Mixed - Stratified" "Ocean - Stratified"

``` r
## DepthMin
# Perform Dunn's test
dunn_result <- dunn.test(Locations_at$DepthMin, Locations_at$Location_type, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 41.2271, df = 2, p-value = 0
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |      Mixed      Ocean
    ## ---------+----------------------
    ##    Ocean |  -2.344649
    ##          |     0.0286
    ##          |
    ## Stratifi |   5.287029   6.406791
    ##          |    0.0000*    0.0000*
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 41.22713
    ## 
    ## $Z
    ## [1] -2.344649  5.287030  6.406792
    ## 
    ## $P
    ## [1] 9.522488e-03 6.215928e-08 7.430676e-11
    ## 
    ## $P.adjusted
    ## [1] 2.856746e-02 1.864779e-07 2.229203e-10
    ## 
    ## $comparisons
    ## [1] "Mixed - Ocean"      "Mixed - Stratified" "Ocean - Stratified"

``` r
## DepthMax
# Perform Dunn's test
dunn_result <- dunn.test(Locations_at$DepthMax, Locations_at$Location_type, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 2.0716, df = 2, p-value = 0.35
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |      Mixed      Ocean
    ## ---------+----------------------
    ##    Ocean |  -0.250346
    ##          |     1.0000
    ##          |
    ## Stratifi |   1.322206   1.426663
    ##          |     0.2791     0.2305
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 2.071623
    ## 
    ## $Z
    ## [1] -0.2503465  1.3222065  1.4266632
    ## 
    ## $P
    ## [1] 0.40115968 0.09304969 0.07683850
    ## 
    ## $P.adjusted
    ## [1] 1.0000000 0.2791491 0.2305155
    ## 
    ## $comparisons
    ## [1] "Mixed - Ocean"      "Mixed - Stratified" "Ocean - Stratified"

``` r
## TempPrefMin
# Perform Dunn's test
dunn_result <- dunn.test(Locations_at$TempPrefMin, Locations_at$Location_type, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 2.0253, df = 2, p-value = 0.36
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |      Mixed      Ocean
    ## ---------+----------------------
    ##    Ocean |  -0.592308
    ##          |     0.8305
    ##          |
    ## Stratifi |   1.125288   1.412638
    ##          |     0.3907     0.2366
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 2.02535
    ## 
    ## $Z
    ## [1] -0.5923085  1.1252884  1.4126388
    ## 
    ## $P
    ## [1] 0.27682202 0.13023341 0.07888097
    ## 
    ## $P.adjusted
    ## [1] 0.8304660 0.3907002 0.2366429
    ## 
    ## $comparisons
    ## [1] "Mixed - Ocean"      "Mixed - Stratified" "Ocean - Stratified"

``` r
## TempPrefMax
# Perform Dunn's test
dunn_result <- dunn.test(Locations_at$TempPrefMax, Locations_at$Location_type, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 1.8962, df = 2, p-value = 0.39
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |      Mixed      Ocean
    ## ---------+----------------------
    ##    Ocean |   1.202254
    ##          |     0.3439
    ##          |
    ## Stratifi |   0.927363   0.275785
    ##          |     0.5306     1.0000
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 1.896171
    ## 
    ## $Z
    ## [1] 1.2022541 0.9273638 0.2757858
    ## 
    ## $P
    ## [1] 0.1146326 0.1768688 0.3913563
    ## 
    ## $P.adjusted
    ## [1] 0.3438977 0.5306065 1.0000000
    ## 
    ## $comparisons
    ## [1] "Mixed - Ocean"      "Mixed - Stratified" "Ocean - Stratified"

``` r
## DorsalSpinesMean
# Perform Dunn's test
dunn_result <- dunn.test(Locations_at$DorsalSpinesMean, Locations_at$Location_type, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 16.838, df = 2, p-value = 0
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |      Mixed      Ocean
    ## ---------+----------------------
    ##    Ocean |  -1.323302
    ##          |     0.2786
    ##          |
    ## Stratifi |   3.481130   4.102539
    ##          |    0.0007*    0.0001*
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 16.83804
    ## 
    ## $Z
    ## [1] -1.323302  3.481131  4.102539
    ## 
    ## $P
    ## [1] 9.286743e-02 2.496508e-04 2.043203e-05
    ## 
    ## $P.adjusted
    ## [1] 2.786023e-01 7.489524e-04 6.129608e-05
    ## 
    ## $comparisons
    ## [1] "Mixed - Ocean"      "Mixed - Stratified" "Ocean - Stratified"

``` r
## # Identify traits that are factors
factor <- c("BodyShapeI", "WaterPref", "DemersPelag", "OperculumPresent", "FeedingPath", "RepGuild2", "ParentalCare")

## BodyShapeI
# Convert to contingency table
contingency_table <- table(Locations_at$Location_type, Locations_at$BodyShapeI)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]
# Perform Fisher's test
fisher_test_result <- fisher.test(contingency_table_clean, simulate.p.value=TRUE)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data with simulated p-value (based on
    ##  2000 replicates)
    ## 
    ## data:  contingency_table_clean
    ## p-value = 0.004998
    ## alternative hypothesis: two.sided

``` r
## WaterPref
# Convert to contingency table
contingency_table <- table(Locations_at$Location_type, Locations_at$WaterPref)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]
# Perform Fisher's test
fisher_test_result <- fisher.test(contingency_table_clean, simulate.p.value=TRUE)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data with simulated p-value (based on
    ##  2000 replicates)
    ## 
    ## data:  contingency_table_clean
    ## p-value = 0.0004998
    ## alternative hypothesis: two.sided

``` r
## DemersPelag
# Convert to contingency table
contingency_table <- table(Locations_at$Location_type, Locations_at$DemersPelag)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]
# Perform Fisher's test
fisher_test_result <- fisher.test(contingency_table_clean, simulate.p.value=TRUE)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data with simulated p-value (based on
    ##  2000 replicates)
    ## 
    ## data:  contingency_table_clean
    ## p-value = 0.0004998
    ## alternative hypothesis: two.sided

``` r
## OperculumPresent
# Convert to contingency table
contingency_table <- table(Locations_at$Location_type, Locations_at$OperculumPresent)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]
# Perform Fisher's test
fisher_test_result <- fisher.test(contingency_table_clean, simulate.p.value=TRUE)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data with simulated p-value (based on
    ##  2000 replicates)
    ## 
    ## data:  contingency_table_clean
    ## p-value = 0.05847
    ## alternative hypothesis: two.sided

``` r
## FeedingPath
# Convert to contingency table
contingency_table <- table(Locations_at$Location_type, Locations_at$FeedingPath)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]
# Perform Fisher's test
fisher_test_result <- fisher.test(contingency_table_clean, simulate.p.value=TRUE)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data with simulated p-value (based on
    ##  2000 replicates)
    ## 
    ## data:  contingency_table_clean
    ## p-value = 0.09595
    ## alternative hypothesis: two.sided

``` r
## RepGuild2
# Convert to contingency table
contingency_table <- table(Locations_at$Location_type, Locations_at$RepGuild2)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]
# Perform Fisher's test
fisher_test_result <- fisher.test(contingency_table_clean, simulate.p.value=TRUE)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data with simulated p-value (based on
    ##  2000 replicates)
    ## 
    ## data:  contingency_table_clean
    ## p-value = 0.03648
    ## alternative hypothesis: two.sided

``` r
## ParentalCare
# Convert to contingency table
contingency_table <- table(Locations_at$Location_type, Locations_at$ParentalCare)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]
# Perform Fisher's test
fisher_test_result <- fisher.test(contingency_table_clean, simulate.p.value=TRUE)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data with simulated p-value (based on
    ##  2000 replicates)
    ## 
    ## data:  contingency_table_clean
    ## p-value = 0.1059
    ## alternative hypothesis: two.sided

## Biotic traits

### Trophic level

``` r
Troph_plot_weighted <- ggplot(Location_type_st_weighted, mapping = aes(x= Strat, y= Troph, fill = Strat)) +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.6,
            width = 0.1) +  
  geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 1, aes(group = Strat, fill = Strat)) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  guides(fill = "none") +
  theme_bw() +
  theme(text = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylab("Trophic Level") +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "a")
```

    ## Warning: The `draw_quantiles` argument of `geom_violin()` is deprecated as of ggplot2
    ## 4.0.0.
    ## ℹ Please use the `quantiles.linetype` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
Troph_plot_weighted
```

![](trait_plots_files/figure-gfm/Trophic%20level-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/Troph_plot_weighted.jpg", Troph_plot_weighted, width = 3.25, height = 2.95, units = "in")


# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(Troph ~ Strat, data = Location_type_st_weighted)
# Print result
print(kruskal_result)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  Troph by Strat
    ## Kruskal-Wallis chi-squared = 29.743, df = 2, p-value = 3.479e-07

``` r
# Perform Dunn's test
dunn_result <- dunn.test(Location_type_st_weighted$Troph, Location_type_st_weighted$Strat, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 29.7428, df = 2, p-value = 0
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |          M          O
    ## ---------+----------------------
    ##        O |   4.364086
    ##          |    0.0000*
    ##          |
    ##        S |  -2.191762  -4.438696
    ##          |     0.0426    0.0000*
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 29.7428
    ## 
    ## $Z
    ## [1]  4.364086 -2.191762 -4.438696
    ## 
    ## $P
    ## [1] 6.382766e-06 1.419834e-02 4.525275e-06
    ## 
    ## $P.adjusted
    ## [1] 1.914830e-05 4.259501e-02 1.357583e-05
    ## 
    ## $comparisons
    ## [1] "M - O" "M - S" "O - S"

### Reproductive Guild2

``` r
Location_type_st_props_RepGuild2 <- Location_type_st %>%
  group_by(Strat, RepGuild2) %>%
  summarise(count = sum(LocationSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
Location_type_st_props_RepGuild2 <- na.omit(Location_type_st_props_RepGuild2)

RepGuild2_plot_weighted <- ggplot(Location_type_st_props_RepGuild2, mapping = aes(x= Strat, y= proportion, fill = Strat, shape = RepGuild2)) + 
  geom_point(position = "identity", size = 4, aes(group = Strat, fill = Strat)) +
  geom_line(position = "identity", aes(group = RepGuild2), linewidth = 1, linetype = "dotted") +
  scale_shape_manual(name = "Egg Strategy:", values = c(11, 21:25), labels = c('1ib' = 'bearers', '6s' = 'scatterers', '3n' = 'nesters', '5h' = 'hiders', '4t' = 'tenders', '2eb' = 'brooders')) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.title = element_text(size = 8),  # Adjust legend title size
    legend.text = element_text(size = 8),  # Adjust legend text size
    legend.key.size = unit(0, "cm"),  # Adjust legend key size
    legend.position = "bottom",
    legend.spacing = unit(0, "cm"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  guides(fill = "none", shape = guide_legend(nrow = 2, byrow = TRUE, title.position = "top")) +
  ylim(c(0,0.6)) +
  ylab("Proportion") +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "c")
RepGuild2_plot_weighted
```

![](trait_plots_files/figure-gfm/Reproductive%20Guild2-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/RepGuild2_plot_weighted.jpg", RepGuild2_plot_weighted, width = 3.25, height = 3.81, units = "in")

# Convert to contingency table
contingency_table <- table(Location_type_st_weighted$Strat, Location_type_st_weighted$RepGuild2)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]

fisher_test_result <- fisher.test(contingency_table_clean, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  contingency_table_clean
    ## p-value = 0.04161
    ## alternative hypothesis: two.sided

``` r
# Subset for groups O and M
pairwise_OM <- contingency_table_clean[c("O", "M"), ]
fisher_test_result_OM <- fisher.test(pairwise_OM, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OM)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OM
    ## p-value = 0.4749
    ## alternative hypothesis: two.sided

``` r
# Subset for groups O and S
pairwise_OS <- contingency_table_clean[c("O", "S"), ]
fisher_test_result_OS <- fisher.test(pairwise_OS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OS
    ## p-value = 0.01258
    ## alternative hypothesis: two.sided

``` r
# Subset for groups M and S
pairwise_MS <- contingency_table_clean[c("M", "S"), ]
fisher_test_result_MS <- fisher.test(pairwise_MS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_MS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_MS
    ## p-value = 0.0229
    ## alternative hypothesis: two.sided

``` r
# Adjust for multiple comparisons
p_values <- c(fisher_test_result_OM$p.value, fisher_test_result_OS$p.value, fisher_test_result_MS$p.value)
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
print(adjusted_p_values)
```

    ## [1] 1.00000000 0.03775467 0.06871058

### Dorsal Spines Mean

``` r
DorsalSpinesMean_plot_weighted <- ggplot(Location_type_st_weighted, mapping = aes(x= Strat, y= DorsalSpinesMean, fill = Strat)) +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.6,
            width = 0.1) +
  geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 1, aes(group = Strat, fill = Strat)) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  guides(fill = "none") +
  theme_bw() +
  theme(text = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylab("Dorsal Spines") +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "b")
DorsalSpinesMean_plot_weighted
```

![](trait_plots_files/figure-gfm/Dorsal%20Spines%20Mean-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 10 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 10 rows containing missing values (`geom_point()`). 
# Twelve values with NA
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/DorsalSpinesMean_plot_weighted.jpg", DorsalSpinesMean_plot_weighted, width = 3.25, height = 2.95, units = "in")


# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(DorsalSpinesMean ~ Strat, data = Location_type_st_weighted)
# Print result
print(kruskal_result)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  DorsalSpinesMean by Strat
    ## Kruskal-Wallis chi-squared = 16.838, df = 2, p-value = 0.0002206

``` r
# Perform Dunn's test
dunn_result <- dunn.test(Location_type_st_weighted$DorsalSpinesMean, Location_type_st_weighted$Strat, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 16.838, df = 2, p-value = 0
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |          M          O
    ## ---------+----------------------
    ##        O |  -1.323302
    ##          |     0.2786
    ##          |
    ##        S |   3.481130   4.102539
    ##          |    0.0007*    0.0001*
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 16.83804
    ## 
    ## $Z
    ## [1] -1.323302  3.481131  4.102539
    ## 
    ## $P
    ## [1] 9.286743e-02 2.496508e-04 2.043203e-05
    ## 
    ## $P.adjusted
    ## [1] 2.786023e-01 7.489524e-04 6.129608e-05
    ## 
    ## $comparisons
    ## [1] "M - O" "M - S" "O - S"

### Parental Care

``` r
Location_type_st_props_ParentalCare <- Location_type_st %>%
  group_by(Strat, ParentalCare) %>%
  summarise(count = sum(LocationSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
Location_type_st_props_ParentalCare <- na.omit(Location_type_st_props_ParentalCare)

ParentalCare_plot_weighted <- ggplot(Location_type_st_props_ParentalCare, mapping = aes(x= Strat, y= proportion, fill = Strat, shape = ParentalCare)) + 
  geom_point(position = "identity", size = 4, aes(group = Strat, fill = Strat)) +
  geom_line(position = "identity", aes(group = ParentalCare), linewidth = 1, linetype = "dotted") +
  scale_shape_manual(name = "Parental Care:", values = c(21:25), labels = c('4n' = 'none', '3p' = 'paternal', '2m' = 'maternal', '1b' = 'biparental')) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.title = element_text(size = 8),  # Adjust legend title size
    legend.text = element_text(size = 8),  # Adjust legend text size
    legend.key.size = unit(0, "cm"),  # Adjust legend key size
    legend.position = "bottom",
    legend.spacing = unit(0, "cm"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  guides(fill = "none", shape = guide_legend(nrow = 2, byrow = TRUE, title.position = "top")) +
  ylim(c(0,0.6)) +
  ylab("Proportion") +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "d")
ParentalCare_plot_weighted
```

![](trait_plots_files/figure-gfm/Parental%20Care-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/ParentalCare_plot_weighted.jpg", ParentalCare_plot_weighted, width = 3.25, height = 3.81, units = "in")


# Convert to contingency table
contingency_table <- table(Location_type_st_weighted$Strat, Location_type_st_weighted$ParentalCare)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]

fisher_test_result <- fisher.test(contingency_table_clean, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  contingency_table_clean
    ## p-value = 0.1106
    ## alternative hypothesis: two.sided

``` r
# Subset for groups O and M
pairwise_OM <- contingency_table_clean[c("O", "M"), ]
fisher_test_result_OM <- fisher.test(pairwise_OM, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OM)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OM
    ## p-value = 0.9447
    ## alternative hypothesis: two.sided

``` r
# Subset for groups O and S
pairwise_OS <- contingency_table_clean[c("O", "S"), ]
fisher_test_result_OS <- fisher.test(pairwise_OS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OS
    ## p-value = 0.01877
    ## alternative hypothesis: two.sided

``` r
# Subset for groups M and S
pairwise_MS <- contingency_table_clean[c("M", "S"), ]
fisher_test_result_MS <- fisher.test(pairwise_MS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_MS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_MS
    ## p-value = 0.02724
    ## alternative hypothesis: two.sided

``` r
# Adjust for multiple comparisons
p_values <- c(fisher_test_result_OM$p.value, fisher_test_result_OS$p.value, fisher_test_result_MS$p.value)
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
print(adjusted_p_values)
```

    ## [1] 1.00000000 0.05631347 0.08170764

## Abiotic traits

### Temperature Preference Minimum

``` r
TempPrefMin_plot_weighted <- ggplot(Location_type_st_weighted, mapping = aes(x= Strat, y= TempPrefMin, fill = Strat)) +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.6,
            width = 0.1) +
  geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 1, aes(group = Strat, fill = Strat)) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  guides(fill = "none") +
  theme_bw() +
  theme(text = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylab("Temp Min (Cº)") +
  scale_y_continuous(breaks = c(20,22,24,26,28)) +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "a")
TempPrefMin_plot_weighted
```

![](trait_plots_files/figure-gfm/Temperature%20Preference%20Minimum-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 9 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 9 rows containing missing values (`geom_point()`). 
# Two values with NA
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/TempPrefMin_plot_weighted.jpg", TempPrefMin_plot_weighted, width = 3.25, height = 2.95, units = "in")


# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(TempPrefMin ~ Strat, data = Location_type_st_weighted)
# Print result
print(kruskal_result)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  TempPrefMin by Strat
    ## Kruskal-Wallis chi-squared = 2.0253, df = 2, p-value = 0.3632

``` r
# Perform Dunn's test
dunn_result <- dunn.test(Location_type_st_weighted$TempPrefMin, Location_type_st_weighted$Strat, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 2.0253, df = 2, p-value = 0.36
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |          M          O
    ## ---------+----------------------
    ##        O |  -0.592308
    ##          |     0.8305
    ##          |
    ##        S |   1.125288   1.412638
    ##          |     0.3907     0.2366
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 2.02535
    ## 
    ## $Z
    ## [1] -0.5923085  1.1252884  1.4126388
    ## 
    ## $P
    ## [1] 0.27682202 0.13023341 0.07888097
    ## 
    ## $P.adjusted
    ## [1] 0.8304660 0.3907002 0.2366429
    ## 
    ## $comparisons
    ## [1] "M - O" "M - S" "O - S"

### Water

- Fresh/Brack/Salt

``` r
Location_type_st_props_Water <- Location_type_st %>%
  group_by(Strat, WaterPref) %>%
  summarise(count = sum(LocationSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
#Location_type_st_props_Water$Water <- factor(Location_type_st_props_Water$Water, levels = c("all", "fresh", "fresh-brack", "brack", "brack-salt", "salt"))

Water_plot_weighted <- ggplot(Location_type_st_props_Water, mapping = aes(x= Strat, y= proportion, fill = Strat, shape = WaterPref)) + 
  geom_point(position = "identity", size = 4, aes(group = Strat, fill = Strat)) +
  geom_line(position = "identity", aes(group = WaterPref), linewidth = 1, linetype = "dotted") +
  scale_shape_manual(name = "Water:", values = c(21:23,11,24:25), labels = c('3a' = 'all', '1s' = 'salt', '2bs' = 'brackish-salt', '4b' = 'brack', '5fb' = 'fresh-brackish', '6f' = 'fresh')) + 
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.title = element_text(size = 8),  # Adjust legend title size
    legend.text = element_text(size = 8),  # Adjust legend text size
    legend.key.size = unit(0, "cm"),  # Adjust legend key size
    legend.position = "bottom",
    legend.spacing = unit(0, "cm"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  guides(fill = "none", shape = guide_legend(nrow = 2, byrow = TRUE, title.position = "top")) +
  ylab("Proportion") +
  scale_y_continuous(breaks = c(0,0.5,1), limits = c(0,1)) +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "c")
Water_plot_weighted
```

![](trait_plots_files/figure-gfm/Water-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/Water_plot_weighted.jpg", Water_plot_weighted, width = 3.25, height = 3.81, units = "in")


# Convert to contingency table
contingency_table <- table(Location_type_st_weighted$Strat, Location_type_st_weighted$WaterPref)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]

fisher_test_result <- fisher.test(contingency_table_clean, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  contingency_table_clean
    ## p-value = 4.049e-14
    ## alternative hypothesis: two.sided

``` r
# Subset for groups O and M
pairwise_OM <- contingency_table_clean[c("O", "M"), ]
fisher_test_result_OM <- fisher.test(pairwise_OM, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OM)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OM
    ## p-value = 0.0001219
    ## alternative hypothesis: two.sided

``` r
# Subset for groups O and S
pairwise_OS <- contingency_table_clean[c("O", "S"), ]
fisher_test_result_OS <- fisher.test(pairwise_OS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OS
    ## p-value = 1.024e-15
    ## alternative hypothesis: two.sided

``` r
# Subset for groups M and S
pairwise_MS <- contingency_table_clean[c("M", "S"), ]
fisher_test_result_MS <- fisher.test(pairwise_MS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_MS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_MS
    ## p-value = 4.843e-08
    ## alternative hypothesis: two.sided

``` r
# Adjust for multiple comparisons
p_values <- c(fisher_test_result_OM$p.value, fisher_test_result_OS$p.value, fisher_test_result_MS$p.value)
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
print(adjusted_p_values)
```

    ## [1] 3.657151e-04 3.073329e-15 1.452993e-07

### Temperature Preference Maximum

``` r
TempPrefMax_plot_weighted <- ggplot(Location_type_st_weighted, mapping = aes(x= Strat, y= TempPrefMax, fill = Strat)) +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 0.6,
            width = 0.1) +
  geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 1, aes(group = Strat, fill = Strat)) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  guides(fill = "none") +
  theme_bw() +
  theme(text = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylab("Temp Max (Cº)") +
  scale_y_continuous(breaks = c(27,28,29,30)) +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "b")
TempPrefMax_plot_weighted
```

![](trait_plots_files/figure-gfm/Temperature%20Preference%20Maximum-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 11 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 11 rows containing missing values (`geom_point()`). 
# Two values with NA
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/TempPrefMax_plot_weighted.jpg", TempPrefMax_plot_weighted, width = 3.25, height = 2.95, units = "in")


# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(TempPrefMax ~ Strat, data = Location_type_st_weighted)
# Print result
print(kruskal_result)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  TempPrefMax by Strat
    ## Kruskal-Wallis chi-squared = 1.8962, df = 2, p-value = 0.3875

``` r
# Perform Dunn's test
dunn_result <- dunn.test(Location_type_st_weighted$TempPrefMax, Location_type_st_weighted$Strat, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 1.8962, df = 2, p-value = 0.39
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |          M          O
    ## ---------+----------------------
    ##        O |   1.202254
    ##          |     0.3439
    ##          |
    ##        S |   0.927363   0.275785
    ##          |     0.5306     1.0000
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 1.896171
    ## 
    ## $Z
    ## [1] 1.2022541 0.9273638 0.2757858
    ## 
    ## $P
    ## [1] 0.1146326 0.1768688 0.3913563
    ## 
    ## $P.adjusted
    ## [1] 0.3438977 0.5306065 1.0000000
    ## 
    ## $comparisons
    ## [1] "M - O" "M - S" "O - S"

## Supplementary traits

### Operculum Present

``` r
Location_type_st_props_OperculumPresent <- Location_type_st %>%
  group_by(Strat, OperculumPresent) %>%
  summarise(count = sum(LocationSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
OperculumPresent_plot_weighted <- ggplot(Location_type_st_props_OperculumPresent, mapping = aes(x= Strat, y= proportion, fill = Strat, shape = OperculumPresent)) + 
  geom_point(position = "identity", size = 10, aes(group = Strat, fill = Strat)) +
  geom_line(position = "identity", aes(group = OperculumPresent), linewidth = 2, linetype = "dotted") +
  scale_shape_manual(name = "Operculum:", values = c(21:22)) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  theme_bw() +
  theme(
    text = element_text(size = 30),
    legend.title = element_text(size = 30),  # Adjust legend title size
    legend.text = element_text(size = 30),  # Adjust legend text size
    axis.text = element_text(size = 30, color = "black"),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  guides(fill = "none") + # Remove the legends for Location type
  ylim(c(0,0.75)) +
  ylab("Proportion") +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "a")
OperculumPresent_plot_weighted
```

![](trait_plots_files/figure-gfm/Operculum%20Present-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/OperculumPresent_plot_weighted.jpg", OperculumPresent_plot_weighted,  width = 13, height = 10, units = "in")


# Convert to contingency table
contingency_table <- table(Location_type_st_weighted$Strat, Location_type_st_weighted$OperculumPresent)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]

fisher_test_result <- fisher.test(contingency_table_clean, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  contingency_table_clean
    ## p-value = 0.05839
    ## alternative hypothesis: two.sided

``` r
# Subset for groups O and M
pairwise_OM <- contingency_table_clean[c("O", "M"), ]
fisher_test_result_OM <- fisher.test(pairwise_OM, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OM)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OM
    ## p-value = 0.206
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.5877987 1.1259967
    ## sample estimates:
    ## odds ratio 
    ##  0.8145992

``` r
# Subset for groups O and S
pairwise_OS <- contingency_table_clean[c("O", "S"), ]
fisher_test_result_OS <- fisher.test(pairwise_OS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OS
    ## p-value = 0.0237
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.2831664 0.9461274
    ## sample estimates:
    ## odds ratio 
    ##  0.5161616

``` r
# Subset for groups M and S
pairwise_MS <- contingency_table_clean[c("M", "S"), ]
fisher_test_result_MS <- fisher.test(pairwise_MS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_MS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_MS
    ## p-value = 0.1132
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.354053 1.140195
    ## sample estimates:
    ## odds ratio 
    ##  0.6332218

``` r
# Adjust for multiple comparisons
p_values <- c(fisher_test_result_OM$p.value, fisher_test_result_OS$p.value, fisher_test_result_MS$p.value)
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
print(adjusted_p_values)
```

    ## [1] 0.61809404 0.07109036 0.33961596

### Max Length (TL)

``` r
MaxLengthTL_plot_weighted <- ggplot(Location_type_st_weighted, mapping = aes(x= Strat, y= MaxLengthTL, fill = Strat)) +
  geom_jitter(shape = 21,
            size = 4,
            alpha = 0.6,
            width = 0.1) +
  geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 2, aes(group = Strat, fill = Strat)) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  guides(fill = "none") +
  theme_bw() +
  theme(text = element_text(size = 30), legend.text = element_text(size = 30),
    axis.text = element_text(size = 30, color = "black"),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylab("Length (cm)") +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "b")
MaxLengthTL_plot_weighted
```

![](trait_plots_files/figure-gfm/Max%20Length%20(TL)-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 3 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 3 rows containing missing values (`geom_point()`). 
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/MaxLengthTL_plot_weighted.jpg", MaxLengthTL_plot_weighted,  width = 13, height = 10, units = "in")


# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(MaxLengthTL ~ Strat, data = Location_type_st_weighted)
# Print result
print(kruskal_result)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  MaxLengthTL by Strat
    ## Kruskal-Wallis chi-squared = 1.7554, df = 2, p-value = 0.4157

``` r
# Perform Dunn's test
dunn_result <- dunn.test(Location_type_st_weighted$MaxLengthTL, Location_type_st_weighted$Strat, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 1.7554, df = 2, p-value = 0.42
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |          M          O
    ## ---------+----------------------
    ##        O |   0.967279
    ##          |     0.5001
    ##          |
    ##        S |   1.101793   0.569996
    ##          |     0.4058     0.8530
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 1.755429
    ## 
    ## $Z
    ## [1] 0.9672800 1.1017938 0.5699965
    ## 
    ## $P
    ## [1] 0.1667020 0.1352757 0.2843400
    ## 
    ## $P.adjusted
    ## [1] 0.5001061 0.4058270 0.8530201
    ## 
    ## $comparisons
    ## [1] "M - O" "M - S" "O - S"

### Body Shape

``` r
Location_type_st_props_BodyShape <- Location_type_st %>%
  group_by(Strat, BodyShapeI) %>%
  summarise(count = sum(LocationSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
BodyShape_plot_weighted <- ggplot(Location_type_st_props_BodyShape, mapping = aes(x= Strat, y= proportion, fill = Strat, shape = BodyShapeI)) + 
  geom_point(position = "identity", size = 10, aes(group = Strat, fill = Strat)) +
  geom_line(position = "identity", aes(group = BodyShapeI), linewidth = 2, linetype = "dotted") +
  scale_shape_manual(name = "Body Shape:", values = c(11, 21:25), labels = c("1o" = "other", "2s" = "short deep", "3f" = "fusiform", "4e" = "elongated", "5l" = "eel-like")) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  theme_bw() +
  theme(
    text = element_text(size = 30),
    legend.title = element_text(size = 30),  # Adjust legend title size
    legend.text = element_text(size = 30),  # Adjust legend text size
    axis.text = element_text(size = 30, color = "black"),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  guides(fill = "none") + # Remove the legends for Location type
  ylim(c(0,0.6)) +
  ylab("Proportion") +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "c")
BodyShape_plot_weighted
```

![](trait_plots_files/figure-gfm/Body%20Shape-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/BodyShape_plot_weighted.jpg", BodyShape_plot_weighted,  width = 13, height = 10, units = "in")


# Convert to contingency table
contingency_table <- table(Location_type_st_weighted$Strat, Location_type_st_weighted$BodyShapeI)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]

fisher_test_result <- fisher.test(contingency_table_clean, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  contingency_table_clean
    ## p-value = 0.003607
    ## alternative hypothesis: two.sided

``` r
# Subset for groups O and M
pairwise_OM <- contingency_table_clean[c("O", "M"), ]
fisher_test_result_OM <- fisher.test(pairwise_OM, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OM)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OM
    ## p-value = 0.05253
    ## alternative hypothesis: two.sided

``` r
# Subset for groups O and S
pairwise_OS <- contingency_table_clean[c("O", "S"), ]
fisher_test_result_OS <- fisher.test(pairwise_OS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OS
    ## p-value = 0.001177
    ## alternative hypothesis: two.sided

``` r
# Subset for groups M and S
pairwise_MS <- contingency_table_clean[c("M", "S"), ]
fisher_test_result_MS <- fisher.test(pairwise_MS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_MS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_MS
    ## p-value = 0.06173
    ## alternative hypothesis: two.sided

``` r
# Adjust for multiple comparisons
p_values <- c(fisher_test_result_OM$p.value, fisher_test_result_OS$p.value, fisher_test_result_MS$p.value)
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
print(adjusted_p_values)
```

    ## [1] 0.157601430 0.003530602 0.185198232

### Depth Maximum

``` r
DepthMax_plot_weighted <- ggplot(Location_type_st_weighted, mapping = aes(x= Strat, y= DepthMax, fill = Strat)) +
  geom_jitter(shape = 21,
            size = 4,
            alpha = 0.6,
            width = 0.1) +
  geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 2, aes(group = Strat, fill = Strat)) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  guides(fill = "none") +
  theme_bw() +
  theme(text = element_text(size = 30), legend.text = element_text(size = 30),
    axis.text = element_text(size = 30, color = "black"),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylab("Depth Max (m)") +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "d")
DepthMax_plot_weighted
```

    ## Warning: Removed 6 rows containing non-finite outside the scale range
    ## (`stat_ydensity()`).

    ## Warning: Removed 6 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](trait_plots_files/figure-gfm/Depth%20Maximum-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 10 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 10 rows containing missing values (`geom_point()`). 
# One value with NA
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/DepthMax_plot_weighted.jpg", DepthMax_plot_weighted,  width = 13, height = 10, units = "in")
```

    ## Warning: Removed 6 rows containing non-finite outside the scale range
    ## (`stat_ydensity()`).
    ## Removed 6 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

``` r
# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(DepthMax ~ Strat, data = Location_type_st_weighted)
# Print result
print(kruskal_result)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  DepthMax by Strat
    ## Kruskal-Wallis chi-squared = 2.0716, df = 2, p-value = 0.3549

``` r
# Perform Dunn's test
dunn_result <- dunn.test(Location_type_st_weighted$DepthMax, Location_type_st_weighted$Strat, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 2.0716, df = 2, p-value = 0.35
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |          M          O
    ## ---------+----------------------
    ##        O |  -0.250346
    ##          |     1.0000
    ##          |
    ##        S |   1.322206   1.426663
    ##          |     0.2791     0.2305
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 2.071623
    ## 
    ## $Z
    ## [1] -0.2503465  1.3222065  1.4266632
    ## 
    ## $P
    ## [1] 0.40115968 0.09304969 0.07683850
    ## 
    ## $P.adjusted
    ## [1] 1.0000000 0.2791491 0.2305155
    ## 
    ## $comparisons
    ## [1] "M - O" "M - S" "O - S"

### Feeding Pathway

``` r
Location_type_st_props_FeedingPath <- Location_type_st %>%
  group_by(Strat, FeedingPath) %>%
  summarise(count = sum(LocationSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
Location_type_st_props_FeedingPath <- na.omit(Location_type_st_props_FeedingPath)

FeedingPath_plot_weighted <- ggplot(Location_type_st_props_FeedingPath, mapping = aes(x= Strat, y= proportion, fill = Strat, shape = FeedingPath)) + 
  geom_point(position = "identity", size = 10, aes(group = Strat, fill = Strat)) +
  geom_line(position = "identity", aes(group = FeedingPath), linewidth = 2, linetype = "dotted") +
  scale_shape_manual(name = "Diet Source:", values = c(21:22), labels = c("b" = "benthic", "p" = "pelagic")) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  theme_bw() +
  theme(
    text = element_text(size = 30),
    legend.title = element_text(size = 30),  # Adjust legend title size
    legend.text = element_text(size = 30),  # Adjust legend text size
    axis.text = element_text(size = 30, color = "black"),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  guides(fill = "none") + # Remove the legends for Location type
  ylim(c(0,1)) +
  ylab("Proportion") +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "e")
FeedingPath_plot_weighted
```

![](trait_plots_files/figure-gfm/Feeding%20Pathway-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/FeedingPath_plot_weighted.jpg", FeedingPath_plot_weighted,  width = 13, height = 10, units = "in")


# Convert to contingency table
contingency_table <- table(Location_type_st_weighted$Strat, Location_type_st_weighted$FeedingPath)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]

fisher_test_result <- fisher.test(contingency_table_clean, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  contingency_table_clean
    ## p-value = 0.08615
    ## alternative hypothesis: two.sided

``` r
# Subset for groups O and M
pairwise_OM <- contingency_table_clean[c("O", "M"), ]
fisher_test_result_OM <- fisher.test(pairwise_OM, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OM)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OM
    ## p-value = 0.09968
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.9376926 2.0182487
    ## sample estimates:
    ## odds ratio 
    ##   1.371045

``` r
# Subset for groups O and S
pairwise_OS <- contingency_table_clean[c("O", "S"), ]
fisher_test_result_OS <- fisher.test(pairwise_OS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OS
    ## p-value = 0.07332
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.914516 3.581409
    ## sample estimates:
    ## odds ratio 
    ##   1.839308

``` r
# Subset for groups M and S
pairwise_MS <- contingency_table_clean[c("M", "S"), ]
fisher_test_result_MS <- fisher.test(pairwise_MS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_MS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_MS
    ## p-value = 0.3326
    ## alternative hypothesis: true odds ratio is not equal to 1
    ## 95 percent confidence interval:
    ##  0.6860086 2.5310101
    ## sample estimates:
    ## odds ratio 
    ##   1.342524

``` r
# Adjust for multiple comparisons
p_values <- c(fisher_test_result_OM$p.value, fisher_test_result_OS$p.value, fisher_test_result_MS$p.value)
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
print(adjusted_p_values)
```

    ## [1] 0.2990474 0.2199597 0.9979309

### Depth Minimum

``` r
DepthMin_plot_weighted <- ggplot(Location_type_st_weighted, mapping = aes(x= Strat, y= DepthMin, fill = Strat)) +
  geom_jitter(shape = 21,
            size = 4,
            alpha = 0.6,
            width = 0.1) +
  geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 2, aes(group = Strat, fill = Strat)) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  guides(fill = "none") +
  theme_bw() +
  theme(text = element_text(size = 30), legend.text = element_text(size = 30),
    axis.text = element_text(size = 30, color = "black"),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylab("Depth Min (m)") +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "f")
DepthMin_plot_weighted
```

![](trait_plots_files/figure-gfm/Depth%20Minimum-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 8 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 8 rows containing missing values (`geom_point()`). 
# One value with NA
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/DepthMin_plot_weighted.jpg", DepthMin_plot_weighted,  width = 13, height = 10, units = "in")


# Perform Kruskal-Wallis test
kruskal_result <- kruskal.test(DepthMin ~ Strat, data = Location_type_st_weighted)
# Print result
print(kruskal_result)
```

    ## 
    ##  Kruskal-Wallis rank sum test
    ## 
    ## data:  DepthMin by Strat
    ## Kruskal-Wallis chi-squared = 41.227, df = 2, p-value = 1.116e-09

``` r
# Perform Dunn's test
dunn_result <- dunn.test(Location_type_st_weighted$DepthMin, Location_type_st_weighted$Strat, method = "Bonferroni")
```

    ##   Kruskal-Wallis rank sum test
    ## 
    ## data: x and group
    ## Kruskal-Wallis chi-squared = 41.2271, df = 2, p-value = 0
    ## 
    ## 
    ##                            Comparison of x by group                            
    ##                                  (Bonferroni)                                  
    ## Col Mean-|
    ## Row Mean |          M          O
    ## ---------+----------------------
    ##        O |  -2.344649
    ##          |     0.0286
    ##          |
    ##        S |   5.287029   6.406791
    ##          |    0.0000*    0.0000*
    ## 
    ## alpha = 0.05
    ## Reject Ho if p <= alpha/2

``` r
# Print Dunn's test result
print(dunn_result)
```

    ## $chi2
    ## [1] 41.22713
    ## 
    ## $Z
    ## [1] -2.344649  5.287030  6.406792
    ## 
    ## $P
    ## [1] 9.522488e-03 6.215928e-08 7.430676e-11
    ## 
    ## $P.adjusted
    ## [1] 2.856746e-02 1.864779e-07 2.229203e-10
    ## 
    ## $comparisons
    ## [1] "M - O" "M - S" "O - S"

### DemersPelag

``` r
Location_type_st_props_DemersPelag <- Location_type_st %>%
  group_by(Strat, DemersPelag) %>%
  summarise(count = sum(LocationSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
DemersPelag_plot_weighted <- ggplot(Location_type_st_props_DemersPelag, mapping = aes(x= Strat, y= proportion, fill = Strat, shape = DemersPelag)) + 
  geom_point(position = "identity", size = 10, aes(group = Strat, fill = Strat)) +
  geom_line(position = "identity", aes(group = DemersPelag), linewidth = 2, linetype = "dotted") +
  scale_shape_manual(name = "Demersal Pelagic:", values = c(21:25, 12:13), 
                     labels = c("1r" = "reef-associated", "2pn" = "pelagic-neritic", "3p" = "pelagic", "4po" = "pelagic-oceanic", "5d" = "demersal", '6bp' = 'benthopelagic', '7bd' = 'bathydemersal')) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  theme_bw() +
  theme(
    text = element_text(size = 30),
    legend.title = element_text(size = 30),  # Adjust legend title size
    legend.text = element_text(size = 30),  # Adjust legend text size
    axis.text = element_text(size = 30, color = "black"),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  guides(fill = "none") + # Remove the legends for Location type
  ylim(c(0,1)) +
  ylab("Proportion") +
  xlab("Location type") +
  labs(fill = "Location type:", tag = "g")
DemersPelag_plot_weighted
```

![](trait_plots_files/figure-gfm/DemersPelag-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/DemersPelag_plot_weighted.jpg", DemersPelag_plot_weighted,  width = 13, height = 10, units = "in")


# Convert to contingency table
contingency_table <- table(Location_type_st_weighted$Strat, Location_type_st_weighted$DemersPelag)
contingency_table_clean <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]

fisher_test_result <- fisher.test(contingency_table_clean, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  contingency_table_clean
    ## p-value = 3.554e-05
    ## alternative hypothesis: two.sided

``` r
# Subset for groups O and M
pairwise_OM <- contingency_table_clean[c("O", "M"), ]
fisher_test_result_OM <- fisher.test(pairwise_OM, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OM)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OM
    ## p-value = 0.4912
    ## alternative hypothesis: two.sided

``` r
# Subset for groups O and S
pairwise_OS <- contingency_table_clean[c("O", "S"), ]
fisher_test_result_OS <- fisher.test(pairwise_OS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_OS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_OS
    ## p-value = 6.411e-06
    ## alternative hypothesis: two.sided

``` r
# Subset for groups M and S
pairwise_MS <- contingency_table_clean[c("M", "S"), ]
fisher_test_result_MS <- fisher.test(pairwise_MS, workspace = 2e8)  # Adjust workspace size
print(fisher_test_result_MS)
```

    ## 
    ##  Fisher's Exact Test for Count Data
    ## 
    ## data:  pairwise_MS
    ## p-value = 0.0002654
    ## alternative hypothesis: two.sided

``` r
# Adjust for multiple comparisons
p_values <- c(fisher_test_result_OM$p.value, fisher_test_result_OS$p.value, fisher_test_result_MS$p.value)
adjusted_p_values <- p.adjust(p_values, method = "bonferroni")
print(adjusted_p_values)
```

    ## [1] 1.000000e+00 1.923188e-05 7.963498e-04

### Reproductive Guild1

``` r
Location_type_st_props_RepGuild1 <- Location_type_st %>%
  group_by(Strat, RepGuild1) %>%
  summarise(count = sum(LocationSums)) %>%
  group_by(Strat) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
Location_type_st_props_RepGuild1 <- na.omit(Location_type_st_props_RepGuild1)

RepGuild1_plot_weighted <- ggplot(Location_type_st_props_RepGuild1, mapping = aes(x = Strat, y = proportion, shape = RepGuild1)) + 
  geom_point(position = "identity", size = 10, aes(group = Strat, fill = Strat)) +
  geom_line(position = "identity", aes(group = RepGuild1), linewidth = 2, linetype = "dotted") +
  scale_shape_manual(name = "Egg Care:", values = c(21, 22, 23), labels = c('2g' = 'guarders', '1b' = 'bearers', '3n' = 'nonguarders')) +
  # scale_color_manual(values = custom_shades) +
  scale_fill_manual(values = custom_shades) +
  theme_bw() +
  theme(
    text = element_text(size = 30),
    legend.title = element_text(size = 30),  # Adjust legend title size
    legend.text = element_text(size = 30),  # Adjust legend text size
    axis.text = element_text(size = 30, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  guides(fill = "none") + # Remove the legends for Location type
  ylim(c(0, 0.6)) +
  ylab("Proportion") +
  xlab("Location type") +
  labs(tag = "none")

RepGuild1_plot_weighted
```

![](trait_plots_files/figure-gfm/Reproductive%20Guild1-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/Fish_Biodiversity/figures/FD/trait_plots/RepGuild1_plot_weighted.jpg", RepGuild1_plot_weighted,  width = 13, height = 10, units = "in")
```

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
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] tidyr_1.3.1       phytools_2.4-4    maps_3.4.2.1      ape_5.8-1        
    ##  [5] reshape2_1.4.4    stringr_1.5.1     knitr_1.50        car_3.1-3        
    ##  [9] carData_3.0-5     dunn.test_1.3.6   viridis_0.6.5     viridisLite_0.4.2
    ## [13] ggplot2_4.0.0     dplyr_1.1.4      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fastmatch_1.1-6         gtable_0.3.6            xfun_0.52              
    ##  [4] lattice_0.22-7          numDeriv_2016.8-1.1     quadprog_1.5-8         
    ##  [7] vctrs_0.6.5             tools_4.4.3             generics_0.1.3         
    ## [10] tibble_3.2.1            pkgconfig_2.0.3         Matrix_1.7-3           
    ## [13] RColorBrewer_1.1-3      S7_0.2.0                scatterplot3d_0.3-44   
    ## [16] lifecycle_1.0.4         compiler_4.4.3          farver_2.1.2           
    ## [19] textshaping_1.0.0       mnormt_2.1.1            combinat_0.0-8         
    ## [22] codetools_0.2-20        htmltools_0.5.8.1       yaml_2.3.10            
    ## [25] Formula_1.2-5           crayon_1.5.3            pillar_1.10.2          
    ## [28] MASS_7.3-65             clusterGeneration_1.3.8 iterators_1.0.14       
    ## [31] abind_1.4-8             foreach_1.5.2           nlme_3.1-168           
    ## [34] phangorn_2.12.1         tidyselect_1.2.1        digest_0.6.37          
    ## [37] stringi_1.8.7           purrr_1.0.4             labeling_0.4.3         
    ## [40] fastmap_1.2.0           grid_4.4.3              expm_1.0-0             
    ## [43] cli_3.6.4               magrittr_2.0.3          optimParallel_1.0-2    
    ## [46] withr_3.0.2             scales_1.4.0            DEoptim_2.2-8          
    ## [49] rmarkdown_2.29          igraph_2.1.4            gridExtra_2.3          
    ## [52] ragg_1.4.0              coda_0.19-4.1           evaluate_1.0.3         
    ## [55] doParallel_1.0.17       rlang_1.1.6             Rcpp_1.0.14            
    ## [58] glue_1.8.0              rstudioapi_0.17.1       R6_2.6.1               
    ## [61] plyr_1.8.9              systemfonts_1.2.2
