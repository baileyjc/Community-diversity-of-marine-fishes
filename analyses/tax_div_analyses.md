Taxonomic diversity analyses
================
Bailey Carlson

## R Markdown

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

# Source the R Markdown file
knit("/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.Rmd", output = "/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md")
```

    ## 
    ## 
    ## processing file: /Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.Rmd

    ##   |                  |          |   0%  |                  |          |   3%                                                                          |                  |.         |   6% [Bringing everything together load modifying files packages]             |                  |.         |   9%                                                                          |                  |.         |  12% [Bringing everything together load in modifying files]                   |                  |..        |  16%                                                                          |                  |..        |  19% [Check species names across files]                                       |                  |..        |  22%                                                                          |                  |..        |  25% [Modify environment data]                                                |                  |...       |  28%                                                                          |                  |...       |  31% [Modify incidence matrices]                                              |                  |...       |  34%                                                                          |                  |....      |  38% [Modify phylogeny]                                                       |                  |....      |  41%                                                                          |                  |....      |  44% [Modify trait data]                                                      |                  |.....     |  47%                                                                          |                  |.....     |  50% [Modify community data frames]                                           |                  |.....     |  53%                                                                          |                  |......    |  56% [Modify community trait data]                                            |                  |......    |  59%                                                                          |                  |......    |  62% [Community trait data tests]                                             |                  |.......   |  66%                                                                          |                  |.......   |  69% [Modify stratification data frames]                                      |                  |.......   |  72%                                                                          |                  |........  |  75% [Modify stratification trait data]                                       |                  |........  |  78%                                                                          |                  |........  |  81% [Stratification trait data tests]                                        |                  |........  |  84%                                                                          |                  |......... |  88% [Modify site trait data frames]                                          |                  |......... |  91%                                                                          |                  |......... |  94% [Site trait data tests]                                                  |                  |..........|  97%                                                                          |                  |..........| 100% [Bringing everything together load out modified files and session info]

    ## output file: /Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md

    ## [1] "/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md"

Plot species richness for all sites with more than 1 species

``` r
#Plot species density
SR_env$Stratification <- factor(SR_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

SR_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = reorder(X, row_sum, decreasing = T), y = row_sum, color = Stratification, fill = Stratification)) + 
  geom_bar(stat = 'identity',
    alpha = 1) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +  
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Site", y="Species Richness")
SR_plot
```

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_plot.png", plot = SR_plot, width = 6, height = 4, units = "in")
```

SR Linear Models

``` r
SR_model <- lm(log(row_sum) ~ -1 + volume_m3_w_chemocline + surface_area_m2 + distance_to_ocean_mean_m + tidal_lag_time_minutes + perimeter_fromSat + max_depth + logArea, data = SR_env[c(1:2,4,5,12,17,21:22),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ -1 + volume_m3_w_chemocline + surface_area_m2 + 
    ##     distance_to_ocean_mean_m + tidal_lag_time_minutes + perimeter_fromSat + 
    ##     max_depth + logArea, data = SR_env[c(1:2, 4, 5, 12, 17, 21:22), 
    ##     ])
    ## 
    ## Residuals:
    ##       BCM       CLM       GLK       HLM       NLK       OTM       SLN       TLN 
    ## -0.011296  0.016837  0.004355  0.003482 -0.011481  0.003345 -0.010996  0.004345 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)  
    ## volume_m3_w_chemocline    2.710e-06  4.623e-08  58.618   0.0109 *
    ## surface_area_m2          -6.767e-05  3.602e-06 -18.788   0.0339 *
    ## distance_to_ocean_mean_m -2.574e-03  2.071e-04 -12.429   0.0511 .
    ## tidal_lag_time_minutes    2.100e-02  5.235e-04  40.108   0.0159 *
    ## perimeter_fromSat         3.076e-03  1.751e-04  17.567   0.0362 *
    ## max_depth                -4.648e-03  4.353e-03  -1.068   0.4791  
    ## logArea                  -1.992e-01  2.157e-02  -9.239   0.0686 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02693 on 1 degrees of freedom
    ## Multiple R-squared:      1,  Adjusted R-squared:  0.9998 
    ## F-statistic:  5595 on 7 and 1 DF,  p-value: 0.01029

``` r
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                          Df  Sum Sq Mean Sq   F value   Pr(>F)   
    ## volume_m3_w_chemocline    1 16.5453 16.5453 22821.932 0.004214 **
    ## surface_area_m2           1  1.4372  1.4372  1982.434 0.014296 * 
    ## distance_to_ocean_mean_m  1  3.9781  3.9781  5487.262 0.008594 **
    ## tidal_lag_time_minutes    1  5.0166  5.0166  6919.729 0.007653 **
    ## perimeter_fromSat         1  1.1016  1.1016  1519.466 0.016328 * 
    ## max_depth                 1  0.2517  0.2517   347.184 0.034134 * 
    ## logArea                   1  0.0619  0.0619    85.354 0.068641 . 
    ## Residuals                 1  0.0007  0.0007                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
SR_model <- lm(log(row_sum) ~ -1 + temperature_median + salinity_median + oxygen_median + pH_median, data = SR_env[c(1:2,4,5,12,17,21:22),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ -1 + temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = SR_env[c(1:2, 4, 5, 12, 
    ##     17, 21:22), ])
    ## 
    ## Residuals:
    ##       BCM       CLM       GLK       HLM       NLK       OTM       SLN       TLN 
    ##  0.353551  0.007887  0.315659  0.296944 -0.039081 -0.878352 -0.406700  0.313130 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## temperature_median -0.04927    0.18228  -0.270    0.800
    ## salinity_median     0.07951    0.12381   0.642    0.556
    ## oxygen_median      -0.37828    0.44881  -0.843    0.447
    ## pH_median           0.31010    1.17333   0.264    0.805
    ## 
    ## Residual standard error: 0.5808 on 4 degrees of freedom
    ## Multiple R-squared:  0.9525, Adjusted R-squared:  0.905 
    ## F-statistic: 20.04 on 4 and 4 DF,  p-value: 0.006561

``` r
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                    Df  Sum Sq Mean Sq F value   Pr(>F)   
    ## temperature_median  1 24.5255 24.5255 72.7025 0.001038 **
    ## salinity_median     1  2.2360  2.2360  6.6283 0.061682 . 
    ## oxygen_median       1  0.2588  0.2588  0.7671 0.430561   
    ## pH_median           1  0.0236  0.0236  0.0699 0.804612   
    ## Residuals           4  1.3494  0.3373                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

SR vs log Area

``` r
SR_model <- lm(log(row_sum) ~ -1 + logArea, data = SR_env[c(1:2,4,5,12,17,21:22),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ -1 + logArea, data = SR_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.7653 -0.4049 -0.1650  0.2205  1.3378 
    ## 
    ## Coefficients:
    ##         Estimate Std. Error t value Pr(>|t|)    
    ## logArea  0.16893    0.02613   6.465 0.000345 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.7628 on 7 degrees of freedom
    ## Multiple R-squared:  0.8565, Adjusted R-squared:  0.836 
    ## F-statistic: 41.79 on 1 and 7 DF,  p-value: 0.0003454

``` r
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##           Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## logArea    1 24.3198 24.3198  41.792 0.0003454 ***
    ## Residuals  7  4.0734  0.5819                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
SR_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = logArea, y = log(row_sum), color = Stratification, fill = Stratification)) +
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = SR_env[-20,1], size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.5, begin = 0.45, end = 0.75, discrete = T, option = "G") +  
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Log Area (m"^"2"~")", y= "Log Species Richness")
# SR_plot <- SR_plot + geom_smooth(mapping = aes(x = logArea, y = row_sum), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
SR_plot <- SR_plot + guides(color = guide_legend(override.aes = list(label = "")))
SR_plot
```

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_plot_area.png", SR_plot, width = 6, height = 4, units = "in")

SR_model <- lm(log(row_sum) ~ -1 + distance_to_ocean_mean_m, data = SR_env[c(1:2,4,5,12,17,21:22),])
summary(SR_model)
```

    ## 
    ## Call:
    ## lm(formula = log(row_sum) ~ -1 + distance_to_ocean_mean_m, data = SR_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -1.18400 -0.27517  0.01573  0.69373  2.14391 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)   
    ## distance_to_ocean_mean_m 0.005650   0.001309   4.315   0.0035 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 1.053 on 7 degrees of freedom
    ## Multiple R-squared:  0.7268, Adjusted R-squared:  0.6877 
    ## F-statistic: 18.62 on 1 and 7 DF,  p-value: 0.003502

``` r
anova(SR_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(row_sum)
    ##                          Df Sum Sq Mean Sq F value   Pr(>F)   
    ## distance_to_ocean_mean_m  1 20.635 20.6352  18.619 0.003502 **
    ## Residuals                 7  7.758  1.1083                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
SR_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = distance_to_ocean_mean_m, y = log(row_sum), color = Stratification, fill = Stratification)) +
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = SR_env[-20,1], size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.5, begin = 0.45, end = 0.75, discrete = T, option = "G") +  
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x="Mean Distance to the Ocean (m)", y= "Log Species Richness")
# SR_plot <- SR_plot + geom_smooth(mapping = aes(x = logArea, y = row_sum), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
SR_plot <- SR_plot + guides(color = guide_legend(override.aes = list(label = "")))
SR_plot
```

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_plot_dist.png", SR_plot, width = 6, height = 4, units = "in")
```

SR vs Oxygen Concentration

``` r
SR_oxygen_median_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = oxygen_median, y = row_sum, color = Stratification, fill = Stratification)) +
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.5, begin = 0.45, end = 0.75, discrete = T, option = "G") +  
  theme_bw() +
  theme(text = element_text(size = 18),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Oxygen Concentration (mg/L)", y= "Species Richness") +
  ylim(c(0,110))
SR_oxygen_median_plot <- SR_oxygen_median_plot + geom_smooth(mapping = aes(x = oxygen_median, y = row_sum), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
SR_oxygen_median_plot
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 3 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values (`geom_smooth()`).

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_oxygen_median_plot.png", SR_oxygen_median_plot, width = 6, height = 4, units = "in")
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 3 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

    ## Warning: Removed 15 rows containing missing values (`geom_smooth()`).

SR vs tidal efficiency

``` r
SR_tidal_efficiency_plot <- ggplot(data = SR_env[-20,], mapping = aes(x = tidal_efficiency, y = row_sum, color = Stratification, fill = Stratification)) +
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.5, begin = 0.45, end = 0.75, discrete = T, option = "G") +  
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(text = element_text(size = 24),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Tidal Efficiency", y= "Species Richness") +
  ylim(c(0,110))
SR_tidal_efficiency_plot <- SR_tidal_efficiency_plot + geom_smooth(mapping = aes(x = tidal_efficiency, y = row_sum), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
SR_tidal_efficiency_plot
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1 rows containing non-finite values (`stat_smooth()`).

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_tidal_efficiency_plot.png", SR_tidal_efficiency_plot, width = 6, height = 4, units = "in")
```

    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning: Removed 1 rows containing non-finite values (`stat_smooth()`).
    ## Removed 1 rows containing missing values (`geom_point()`).

Beta Diversity of SR

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
```

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# plot both heatmaps next to each other
SR_beta_div_plot <- gridExtra::grid.arrange(SR.jac.turn, SR.jac.nest, ncol = 2)
```

![](tax_div_analyses_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/species_richness/SR_beta_div_plot.png", SR_beta_div_plot, width = 6, height = 4.5)
```

Determine what families of fish are most present in the lakes

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

``` r
sessionInfo()
```

    ## R version 4.3.1 (2023-06-16)
    ## Platform: aarch64-apple-darwin20 (64-bit)
    ## Running under: macOS Ventura 13.4.1
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
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] tidyr_1.3.0       phytools_2.0-3    maps_3.4.1.1      ape_5.7-1        
    ##  [5] reshape2_1.4.4    stringr_1.5.1     dplyr_1.1.4       ggrepel_0.9.4    
    ##  [9] viridis_0.6.4     viridisLite_0.4.2 ggplot2_3.4.4     betapart_1.6     
    ## [13] knitr_1.45       
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
    ## [25] expm_0.999-8            picante_1.8.2           withr_2.5.2            
    ## [28] purrr_1.0.2             numDeriv_2016.8-1.1     itertools_0.1-3        
    ## [31] grid_4.3.1              fansi_1.0.6             colorspace_2.1-0       
    ## [34] scales_1.3.0            iterators_1.0.14        MASS_7.3-60            
    ## [37] cli_3.6.1               rmarkdown_2.25          vegan_2.6-4            
    ## [40] ragg_1.2.6              generics_0.1.3          rcdd_1.5-2             
    ## [43] rstudioapi_0.15.0       magic_1.6-1             splines_4.3.1          
    ## [46] vctrs_0.6.5             Matrix_1.6-4            minpack.lm_1.2-4       
    ## [49] systemfonts_1.0.5       foreach_1.5.2           snow_0.4-4             
    ## [52] glue_1.6.2              codetools_0.2-19        stringi_1.8.2          
    ## [55] gtable_0.3.4            quadprog_1.5-8          munsell_0.5.0          
    ## [58] tibble_3.2.1            pillar_1.9.0            htmltools_0.5.7        
    ## [61] R6_2.5.1                textshaping_0.3.7       doParallel_1.0.17      
    ## [64] evaluate_0.23           lattice_0.22-5          highr_0.10             
    ## [67] Rcpp_1.0.11             fastmatch_1.1-4         coda_0.19-4            
    ## [70] gridExtra_2.3           nlme_3.1-164            permute_0.9-7          
    ## [73] mgcv_1.9-0              xfun_0.41               pkgconfig_2.0.3

``` r
citation()
```

    ## To cite R in publications use:
    ## 
    ##   R Core Team (2023). _R: A Language and Environment for Statistical
    ##   Computing_. R Foundation for Statistical Computing, Vienna, Austria.
    ##   <https://www.R-project.org/>.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {R: A Language and Environment for Statistical Computing},
    ##     author = {{R Core Team}},
    ##     organization = {R Foundation for Statistical Computing},
    ##     address = {Vienna, Austria},
    ##     year = {2023},
    ##     url = {https://www.R-project.org/},
    ##   }
    ## 
    ## We have invested a lot of time and effort in creating R, please cite it
    ## when using it for data analysis. See also 'citation("pkgname")' for
    ## citing R packages.
