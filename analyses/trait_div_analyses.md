Trait diversity analyses
================
Bailey Carlson

## R Markdown

``` r
# Load the knitr package if not already loaded
library(knitr)

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
library(picante)
```

    ## Loading required package: vegan

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.6-4

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
```

# Functional Diversity MNTD, MPD, and PD

``` r
# Trait distance calculation
# BodyShapeI, DemersPelag, OperculumPresent, Troph, DepthMin, RepGuild1, Habitat, DorsalSpinesMean
traits_dist <- gowdis(traits)
traits_clust <- hclust(traits_dist,"average")
traits_tree <- as.phylo(traits_clust)


#Mean pairwise differences
traits_sesmpd <- ses.mpd(presabs_lake, traits_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(traits_sesmpd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/traits_sesmpd.csv")


#Mean nearest taxon distance
traits_sesmntd <- ses.mntd(presabs_lake, traits_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(traits_sesmntd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/traits_sesmntd.csv")


#Faith's PD
traits_sespd <- ses.pd(presabs_lake, traits_tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)

write.csv(traits_sespd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/traits_sespd.csv")


# Trait distance calculation
straits_dist <- gowdis(straits)
straits_clust <- hclust(straits_dist,"average")
straits_tree <- as.phylo(straits_clust)


#Mean pairwise differences
straits_sesmpd <- ses.mpd(surveyed_sites_lake, straits_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(straits_sesmpd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/straits_sesmpd.csv")


#Mean nearest taxon distance
straits_sesmntd <- ses.mntd(surveyed_sites_lake, straits_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(straits_sesmntd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/straits_sesmntd.csv")


#Faith's PD
straits_sespd <- ses.pd(surveyed_sites_lake, straits_tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)

write.csv(straits_sespd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/straits_sespd.csv")
```

Subset of traits Functional Diversity MNTD, MPD, and PD

``` r
# Trait distance calculation
traits_subset <- traits[,c(1,3,5,13,15,16,17)]
traits_sub_dist <- gowdis(traits_subset)
traits_sub_clust <- hclust(traits_sub_dist,"average")
traits_sub_tree <- as.phylo(traits_sub_clust)

#Mean pairwise differences
traits_sub_sesmpd <- ses.mpd(presabs_lake, traits_sub_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

#Mean nearest taxon distance
traits_sub_sesmntd <- ses.mntd(presabs_lake, traits_sub_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

#Faith's PD
traits_sub_sespd <- ses.pd(presabs_lake, traits_sub_tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)


# Trait distance calculation
straits_subset <- straits[,c(1,3,5,13,15,16,17)]
straits_sub_dist <- gowdis(straits_subset)
straits_sub_clust <- hclust(straits_sub_dist,"average")
straits_sub_tree <- as.phylo(straits_sub_clust)

#Mean pairwise differences
straits_sub_sesmpd <- ses.mpd(surveyed_sites_lake, straits_sub_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

#Mean nearest taxon distance
straits_sub_sesmntd <- ses.mntd(surveyed_sites_lake, straits_sub_dist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

#Faith's PD
straits_sub_sespd <- ses.pd(surveyed_sites_lake, straits_sub_tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)
```

# Read in FD mpd, mntd, pd

Read in Functional Diversity files and combine with env data

``` r
traits_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/traits_sesmpd.csv")

traits_sesmpd_env <- merge(traits_sesmpd, env, by = "X", sort = F)
traits_sesmpd_env$measure <- "traits_sesmpd"
traits_sesmpd_env$Stratification <- factor(traits_sesmpd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))


traits_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/traits_sesmntd.csv")

traits_sesmntd_env <- merge(traits_sesmntd, env, by = "X", sort = F)
traits_sesmntd_env$measure <- "traits_sesmntd"
traits_sesmntd_env$Stratification <- factor(traits_sesmntd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))


traits_sespd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/traits_sespd.csv")

traits_sespd_env <- merge(traits_sespd, env, by = "X", sort = F)
traits_sespd_env$measure <- "traits_sespd"
traits_sespd_env$Stratification <- factor(traits_sespd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))


straits_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/straits_sesmpd.csv")

straits_sesmpd_env <- merge(straits_sesmpd, env, by = "X", sort = F)
straits_sesmpd_env$measure <- "straits_sesmpd"
straits_sesmpd_env$Stratification <- factor(straits_sesmpd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))


straits_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/straits_sesmntd.csv")

straits_sesmntd_env <- merge(straits_sesmntd, env, by = "X", sort = F)
straits_sesmntd_env$measure <- "straits_sesmntd"
straits_sesmntd_env$Stratification <- factor(straits_sesmntd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))


straits_sespd <-read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/straits_sespd.csv")

straits_sespd_env <- merge(straits_sespd, env, by = "X", sort = F)
straits_sespd_env$measure <- "straits_sespd"
straits_sespd_env$Stratification <- factor(straits_sespd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))
```

FD Linear Models

``` r
FD_model <- lm(pd.obs ~ -1 + volume_m3_w_chemocline + surface_area_m2 + distance_to_ocean_mean_m + tidal_lag_time_minutes + perimeter_fromSat + max_depth + logArea, data = traits_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(FD_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ -1 + volume_m3_w_chemocline + surface_area_m2 + 
    ##     distance_to_ocean_mean_m + tidal_lag_time_minutes + perimeter_fromSat + 
    ##     max_depth + logArea, data = traits_sespd_env[c(1:2, 4, 5, 
    ##     12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##         1         2         4         5        12        17        21        22 
    ## -0.032982  0.049161  0.012716  0.010167 -0.033523  0.009766 -0.032105  0.012687 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)  
    ## volume_m3_w_chemocline    1.559e-06  1.350e-07  11.546   0.0550 .
    ## surface_area_m2          -1.808e-05  1.052e-05  -1.719   0.3354  
    ## distance_to_ocean_mean_m -2.441e-03  6.046e-04  -4.037   0.1546  
    ## tidal_lag_time_minutes    1.027e-02  1.529e-03   6.719   0.0941 .
    ## perimeter_fromSat         8.845e-04  5.113e-04   1.730   0.3337  
    ## max_depth                -3.067e-02  1.271e-02  -2.413   0.2501  
    ## logArea                  -1.848e-02  6.297e-02  -0.294   0.8182  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.07862 on 1 degrees of freedom
    ## Multiple R-squared:  0.999,  Adjusted R-squared:  0.9923 
    ## F-statistic: 147.7 on 7 and 1 DF,  p-value: 0.06328

``` r
anova(FD_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##                          Df Sum Sq Mean Sq  F value  Pr(>F)  
    ## volume_m3_w_chemocline    1 3.8113  3.8113 616.6456 0.02562 *
    ## surface_area_m2           1 0.1324  0.1324  21.4251 0.13545  
    ## distance_to_ocean_mean_m  1 0.6600  0.6600 106.7896 0.06141 .
    ## tidal_lag_time_minutes    1 1.2481  1.2481 201.9445 0.04472 *
    ## perimeter_fromSat         1 0.3527  0.3527  57.0664 0.08379 .
    ## max_depth                 1 0.1837  0.1837  29.7195 0.11549  
    ## logArea                   1 0.0005  0.0005   0.0862 0.81822  
    ## Residuals                 1 0.0062  0.0062                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
FD_model <- lm(pd.obs ~ -1 + temperature_median + salinity_median + oxygen_median + pH_median, data = traits_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(FD_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ -1 + temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = traits_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##        1        2        4        5       12       17       21       22 
    ##  0.25021 -0.01287  0.19216  0.20076 -0.10330 -0.39603 -0.28330  0.13057 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)
    ## temperature_median -0.06690    0.09984  -0.670    0.540
    ## salinity_median     0.02767    0.06781   0.408    0.704
    ## oxygen_median      -0.30943    0.24583  -1.259    0.277
    ## pH_median           0.41550    0.64266   0.647    0.553
    ## 
    ## Residual standard error: 0.3181 on 4 degrees of freedom
    ## Multiple R-squared:  0.9367, Adjusted R-squared:  0.8734 
    ## F-statistic:  14.8 on 4 and 4 DF,  p-value: 0.01151

``` r
anova(FD_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##                    Df Sum Sq Mean Sq F value   Pr(>F)   
    ## temperature_median  1 4.9583  4.9583 48.9932 0.002193 **
    ## salinity_median     1 0.8673  0.8673  8.5696 0.042925 * 
    ## oxygen_median       1 0.1223  0.1223  1.2080 0.333437   
    ## pH_median           1 0.0423  0.0423  0.4180 0.553165   
    ## Residuals           4 0.4048  0.1012                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Plot functional richness for all sites with more than 1 species

``` r
fr_model <- lm(pd.obs ~ -1 + logArea, data = traits_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(fr_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ -1 + logArea, data = traits_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.44729 -0.23022 -0.14599  0.08025  0.76882 
    ## 
    ## Coefficients:
    ##         Estimate Std. Error t value Pr(>|t|)   
    ## logArea  0.07620    0.01558   4.892  0.00177 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.4547 on 7 degrees of freedom
    ## Multiple R-squared:  0.7737, Adjusted R-squared:  0.7413 
    ## F-statistic: 23.93 on 1 and 7 DF,  p-value: 0.00177

``` r
anova(fr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##           Df Sum Sq Mean Sq F value  Pr(>F)   
    ## logArea    1 4.9475  4.9475  23.927 0.00177 **
    ## Residuals  7 1.4474  0.2068                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
fr_plot <- ggplot(data = traits_sespd_env[-20,], mapping = aes(y = pd.obs, x = logArea, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = traits_sespd_env[-20,1], size = 7, point.padding = 3) +
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
  labs(x="Log Area (m"^"2"~")", y="Functional Richness")
fr_plot <- fr_plot + guides(color = guide_legend(override.aes = list(label = "")))
fr_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/fr_plot_area.png", plot = fr_plot, width = 6, height = 4, units = "in")

fr_model <- lm(pd.obs ~ -1 + distance_to_ocean_mean_m, data = traits_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(fr_model)
```

    ## 
    ## Call:
    ## lm(formula = pd.obs ~ -1 + distance_to_ocean_mean_m, data = traits_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.63608 -0.17152 -0.05091  0.30606  1.13571 
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error t value Pr(>|t|)   
    ## distance_to_ocean_mean_m 0.0025214  0.0007107   3.548  0.00937 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.5714 on 7 degrees of freedom
    ## Multiple R-squared:  0.6426, Adjusted R-squared:  0.5916 
    ## F-statistic: 12.59 on 1 and 7 DF,  p-value: 0.009369

``` r
anova(fr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: pd.obs
    ##                          Df Sum Sq Mean Sq F value   Pr(>F)   
    ## distance_to_ocean_mean_m  1 4.1096  4.1096  12.588 0.009369 **
    ## Residuals                 7 2.2853  0.3265                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
fr_plot <- ggplot(data = traits_sespd_env[-20,], mapping = aes(y = pd.obs, x = distance_to_ocean_mean_m, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = traits_sespd_env[-20,1], size = 7, point.padding = 3) +
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
  labs(x="Mean Distance to the Ocean (m)", y="Functional Richness")
fr_plot <- fr_plot + guides(color = guide_legend(override.aes = list(label = "")))
fr_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/fr_plot_dist.png", plot = fr_plot, width = 6, height = 4, units = "in")
```

Traits Z score vs p value plots

``` r
traits_sesmpd_plot <- ggplot(data = traits_sesmpd_env, mapping = aes(y = mpd.obs.p, x = mpd.obs.z, color = Stratification, fill = Stratification, label = X)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = traits_sesmpd_env[,1], size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x= "z-score", y= "p-value", colour = "Stratification", tag = "A") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
traits_sesmpd_plot <- traits_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
traits_sesmpd_plot
```

    ## Warning: ggrepel: 2 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](trait_div_analyses_files/figure-gfm/FD%20z%20score%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/traits_sesmpd_plot_z.png", traits_sesmpd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
traits_sesmntd_plot <- ggplot(data = traits_sesmntd_env, mapping = aes(y = mntd.obs.p, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = traits_sesmntd_env[,1], size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x= "z-score", y= "p-value", colour = "Stratification", tag = "B") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
traits_sesmntd_plot <- traits_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
traits_sesmntd_plot
```

    ## Warning: ggrepel: 10 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](trait_div_analyses_files/figure-gfm/FD%20z%20score%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/traits_sesmntd_plot_z.png", traits_sesmntd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 12 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
traits_sespd_plot <- ggplot(data = traits_sespd_env, mapping = aes(y = pd.obs.p, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = traits_sespd_env[,1], size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(x= "z-score", y= "p-value", colour = "Stratification", tag = "C") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) + 
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
traits_sespd_plot <- traits_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
traits_sespd_plot
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1 rows containing missing values (`geom_text_repel()`).

    ## Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](trait_div_analyses_files/figure-gfm/FD%20z%20score%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/traits_sespd_plot_z.png", traits_sespd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1 rows containing missing values (`geom_text_repel()`).

    ## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

Traits Z score vs obs value plots

``` r
traits_sesmpd_plot <- ggplot(data = traits_sesmpd_env, mapping = aes(y = mpd.obs, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed MPD", x= "z-score") +
  # ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
traits_sesmpd_plot <- traits_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
traits_sesmpd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/traits_sesmpd_plot_obs+z.png", traits_sesmpd_plot, width = 6, height = 4, units = "in")

traits_sesmntd_plot <- ggplot(data = traits_sesmntd_env, mapping = aes(y = mntd.obs, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed MNTD", x= "z-score") +
  # ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
traits_sesmntd_plot <- traits_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
traits_sesmntd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/traits_sesmntd_plot_obs+z.png", traits_sesmntd_plot, width = 6, height = 4, units = "in")

traits_sespd_plot <- ggplot(data = traits_sespd_env, mapping = aes(y = pd.obs, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed FD", x= "z-score") +
  ylim(c(0,6)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=6, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
traits_sespd_plot <- traits_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
traits_sespd_plot
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/traits_sespd_plot_obs+z.png", traits_sespd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

Straits Z score vs p value plots

``` r
straits_sesmpd_plot <- ggplot(data = straits_sesmpd_env, mapping = aes(y = mpd.obs.p, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = straits_sesmpd_env[,1], size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "p-value", x= "z-score", colour = "Stratification", tag = "D") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sesmpd_plot <- straits_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sesmpd_plot
```

    ## Warning: ggrepel: 7 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/straits_sesmpd_plot_z.png", straits_sesmpd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
straits_sesmntd_plot <- ggplot(data = straits_sesmntd_env, mapping = aes(y = mntd.obs.p, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = straits_sesmntd_env[,1], size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "p-value", x= "z-score", colour = "Stratification", tag = "E") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sesmntd_plot <- straits_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sesmntd_plot
```

    ## Warning: ggrepel: 2 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/straits_sesmntd_plot_z.png", straits_sesmntd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
straits_sespd_plot <- ggplot(data = straits_sespd_env, mapping = aes(y = pd.obs.p, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = straits_sespd_env[,1], size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "p-value", x= "z-score", colour = "Stratification", tag = "F") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sespd_plot <- straits_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sespd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/straits_sespd_plot_z.png", straits_sespd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 2 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

Straits Z score vs obs value plots

``` r
straits_sesmpd_plot <- ggplot(data = straits_sesmpd_env, mapping = aes(y = mpd.obs, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed MPD", x= "z-score") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sesmpd_plot <- straits_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sesmpd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/straits_sesmpd_plot_obs+z.png", straits_sesmpd_plot, width = 6, height = 4, units = "in")

straits_sesmntd_plot <- ggplot(data = straits_sesmntd_env, mapping = aes(y = mntd.obs, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed MNTD", x= "z-score") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sesmntd_plot <- straits_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sesmntd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/straits_sesmntd_plot_obs+z.png", straits_sesmntd_plot, width = 6, height = 4, units = "in")

straits_sespd_plot <- ggplot(data = straits_sespd_env, mapping = aes(y = pd.obs, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  labs(y= "Observed FD", x= "z-score") +
  ylim(c(0,6)) +
  scale_x_continuous(breaks = c(-6,-4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=6, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
straits_sespd_plot <- straits_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
straits_sespd_plot
```

![](trait_div_analyses_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/straits_sespd_plot_obs+z.png", straits_sespd_plot, width = 6, height = 4, units = "in")
```

Regression of MPD

``` r
#MPD
man_traits_sesmpd_env_chem <- manova(cbind(temperature_median, salinity_median, oxygen_median, pH_median) ~ mpd.obs.z, data = traits_sesmpd_env)
summary.aov(man_traits_mpd_env_chem)

# traits_sesmpd_env_chem_model <- glmulti("mpd.obs.z", c("temperature_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = traits_sesmpd_env, method = "h")

man_traits_sesmpd_env_phys <- manova(cbind(surface_area_m2, volume_m3_w_chemocline, distance_to_ocean_min_m, tidal_efficiency, depth) ~ mpd.obs.z, data = traits_sesmpd_env)
summary.aov(man_traits_mpd_env_phys)

# traits_sesmpd_env_phys_model <- glmulti("mpd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "depth", "logArea"), intercept = FALSE, level = 1, data = traits_sesmpd_env, method = "h")
# 
# traits_sesmpd_env_model <- glmulti("mpd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "tidal_efficiency", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = traits_sesmpd_env, method = "h")

traits_sesmpd_env_model <- lm(mpd.obs.z ~ -1 + distance_to_ocean_min_m + tidal_efficiency + salinity_median + oxygen_median, data = traits_sesmpd_env)
summary(traits_mpd_env_model)
anova(traits_mpd_env_model)
```

Plot MPD vs Oxygen Concentration

``` r
traits_sesmpd_oxygen_plot <- ggplot(data = traits_sesmpd_env, mapping = aes(x = oxygen_median, y = mpd.obs.z, color = Stratification, fill = Stratification, label = X)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 12),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Median Oxygen Concentration (mg/L)", y= "sesMPD") +
  ylim(c(-6,3))
traits_sesmpd_oxygen_plot <- traits_sesmpd_oxygen_plot + geom_smooth(mapping = aes(x = oxygen_median, y = mpd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("traits_sesmpd_oxygen_plot.png", traits_sesmpd_oxygen_plot, width = 2, height = 2)
```

Plot MPD vs Tidal Efficiency

``` r
traits_sesmpd_oxygen_plot <- ggplot(data = traits_sesmpd_env, mapping = aes(x = oxygen_median, y = mpd.obs.z, color = Stratification, fill = Stratification, label = X)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 12),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Median Oxygen Concentration (mg/L)", y= "sesMPD") +
  ylim(c(-6,3))
traits_sesmpd_oxygen_plot <- traits_sesmpd_oxygen_plot + geom_smooth(mapping = aes(x = oxygen_median, y = mpd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("traits_sesmpd_oxygen_plot.png", traits_sesmpd_oxygen_plot, width = 2, height = 2)
```

Regression of MNTD

``` r
#MNTD
man_traits_sesmntd_env_chem <- manova(cbind(temperature_median, conductivity_median, salinity_median, oxygen_median, pH_median) ~ mntd.obs.z, data = traits_sesmntd_env)
summary.aov(man_traits_mntd_env_chem)

man_traits_sesmntd_env_phys <- manova(cbind(surface_area_m2, volume_m3_w_chemocline, distance_to_ocean_min_m, tidal_efficiency, depth) ~ mntd.obs.z, data = traits_sesmntd_env)
summary.aov(man_traits_mntd_env_phys)

traits_sesmntd_env_phys_model <- glmulti("mntd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "depth", "logArea"), intercept = FALSE, level = 1, data = traits_sesmntd_env, method = "h")

traits_sesmntd_env_chem_model <- glmulti("mntd.obs.z", c("temperature_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = traits_sesmntd_env, method = "h")

traits_sesmntd_env_model <- glmulti("mntd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "tidal_efficiency", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = traits_sesmntd_env, method = "h")

traits_sesmntd_env_model <- lm(mntd.obs.z ~ -1 + distance_to_ocean_min_m + salinity_median + oxygen_median, data = traits_sesmntd_env)
summary(traits_mntd_env_model)
anova(traits_mntd_env_model)
```

Plot MNTD vs Oxygen Concentration

``` r
traits_sesmntd_oxygen_plot <- ggplot(data = traits_sesmntd_env, mapping = aes(x = oxygen_median, y = mntd.obs.z, color = Stratification, fill = Stratification, label = X)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 12),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Median Oxygen Concentration (mg/L)", y= "sesMNTD") +
  ylim(c(-6,2))
traits_sesmntd_oxygen_plot <- traits_sesmntd_oxygen_plot + geom_smooth(mapping = aes(x = oxygen_median, y = mntd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("traits_sesmntd_oxygen_plot.png", traits_sesmntd_oxygen_plot, width = 2, height = 2)
```

Plot MNTD vs Tidal Efficiency

``` r
traits_sesmntd_oxygen_plot <- ggplot(data = traits_sesmntd_env, mapping = aes(x = oxygen_median, y = mntd.obs.z, color = Stratification, fill = Stratification, label = X)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 12),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Median Oxygen Concentration (mg/L)", y= "sesMNTD") +
  ylim(c(-6,2))
traits_sesmntd_oxygen_plot <- traits_sesmntd_oxygen_plot + geom_smooth(mapping = aes(x = oxygen_median, y = mntd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("traits_sesmntd_oxygen_plot.png", traits_sesmntd_oxygen_plot, width = 2, height = 2)
```

Regression of PD

``` r
#PD
man_traits_sespd_env_chem <- manova(cbind(temperature_median, salinity_median, oxygen_median, pH_median) ~ pd.obs.z, data = traits_sespd_env)
summary.aov(man_traits_sespd_env_chem)

lm_traits_sespd_env_chem <- lm(pd.obs.z ~ oxygen_median, data = traits_sespd_env)
summary(lm_traits_sespd_env_chem)

traits_sespd_env_chem_model <- glmulti("pd.obs.z", c("temperature_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = traits_sespd_env, method = "h")

man_traits_sespd_env_phys <- manova(cbind(surface_area_m2, volume_m3_w_chemocline, distance_to_ocean_min_m, tidal_efficiency, depth, logArea) ~ pd.obs.z, data = traits_sespd_env)
summary.aov(man_traits_sespd_env_phys)

lm_traits_sespd_tidal <- lm(pd.obs.z ~ tidal_efficiency, data = traits_sespd_env)
summary(lm_traits_sespd_tidal)

traits_sespd_env_phys_model <- glmulti("pd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "depth", "logArea"), intercept = FALSE, level = 1, data = traits_sespd_env, method = "h")

traits_sespd_env_model <- lm(pd.obs.z ~ -1 + distance_to_ocean_min_m + tidal_efficiency + oxygen_median + pH_median, data = traits_sespd_env)
summary(traits_pd_env_model)
anova(traits_pd_env_model)
```

Plot PD vs Oxygen Concentration

``` r
traits_sespd_env_plot <- ggplot(data = traits_sespd_env, mapping = aes(x = oxygen_median, y = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 3,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.5, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 18),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Env Variable", y= "sesFPD")
traits_sespd_env_plot + geom_smooth(mapping = aes(x = env_variable, y = pd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("traits_sespd_env_plot.png", traits_sespd_env_plot, width = 6, height = 4, units = "in")
```

Plot PD vs Tidal Efficiency

``` r
traits_sespd_env_plot <- ggplot(data = traits_sespd_env, mapping = aes(x = oxygen_median, y = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 3,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.5, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 18),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Env Variable", y= "sesFPD")
traits_sespd_env_plot + geom_smooth(mapping = aes(x = env_variable, y = pd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("traits_sespd_env_plot.png", traits_sespd_env_plot, width = 6, height = 4, units = "in")
```

Plot z-scores and environmental variables

``` r
traits_phylo_ses_tidal_plot <- ggplot(data =  traits_phylo_ses_env, mapping = aes(x = tidal_efficiency, y = obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.5, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 24),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(x= "Tidal Efficiency", y= "z-score") +
  ylim(c(-10,10))
traits_phylo_ses_tidal_plot <- traits_phylo_ses_tidal_plot + geom_smooth(mapping = aes(x = tidal_efficiency, y = obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
traits_phylo_ses_tidal_plot <- traits_phylo_ses_tidal_plot + facet_wrap(ncol=3)
traits_phylo_ses_tidal_plot
ggsave("traits_phylo_ses_tidal_plot.png", traits_phylo_ses_tidal_plot, width = 2, height = 2)


traits_phylo_ses_oxygen_plot <- ggplot(data = traits_phylo_FD_all_env[-c(65:72,121:128,137:144,153:160),], mapping = aes(x = oxygen_median, y = obs.z, color = Environment)) + 
  geom_point(stat = 'identity',
    size = 8,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, end = .95, discrete = TRUE, option = "H") +
  theme_bw() +
  theme(text = element_text(size = 24),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Median Oxygen Concentration (mg/L)", y= "z-score") +
  ylim(c(-10,10))
traits_phylo_ses_oxygen_plot <- traits_phylo_ses_oxygen_plot + geom_smooth(mapping = aes(x = oxygen_median, y = obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
traits_phylo_ses_oxygen_plot <- traits_phylo_ses_oxygen_plot + facet_wrap(vars(Indice), ncol=2)
traits_phylo_ses_oxygen_plot
ggsave("traits_phylo_ses_oxygen_plot.png", traits_phylo_ses_oxygen_plot, width = 2, height = 2)
```

Beta Diversity of Traits

``` r
traits_dist <- gowdis(main_traits)
traits_dist2 <- cailliez(traits_dist)
traits_dist2_pco <- dudi.pco(traits_dist2, scannf = FALSE, full = TRUE)
traits_pco_round <- round(traits_dist2_pco$li, .Machine$double.exponent)
traits_pco_complete <- traits_pco_round[, 1:2]
qual.FRic <- sum(traits_dist2_pco$eig[1:2]) / sum(traits_dist2_pco$eig)

# Do the following to order lakes by the circulation pattern, then delete env data
####https://pedrohbraga.github.io/StratificationPhylogenetics-Workshop/StratificationPhylogenetics-Workshop.html#between-assemblage-phylogenetic-structure####
presabs_env <- merge(pres_abs_by_lake, env_byLake, by = 'X')
presabs_env <- presabs_env[order(presabs_env$Environment),]
row.names(presabs_env) <- presabs_env$X
presabs_env_ordered <- presabs_env[,-c(1, 1692:1722)]

trait_beta <- functional.betapart.core.pairwise(presabs_env_ordered[-23,], traits_dist, parallel = TRUE)


trait_betap_jac <- functional.beta.pair(presabs_env_ordered[-23,], main_traits[,-c(1,5)], index.family="jaccard")

# Turnover matrix:
trait.jac.turn <- trait_betap_jac$beta.jtu %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("sites") + ylab("sites") + labs(fill = "Turnover") + # edit axis and legend titles
  scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels

# Nestedness matrix:
trait.jac.nest <- trait_betap_jac$beta.jne %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("sites") + ylab("sites") + labs(fill = "Nestedness") + # edit axis and legend titles
    scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels

# Nestedness matrix:
trait.jac <- trait_betap_jac$beta.jac %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("sites") + ylab("sites") + labs(fill = "Beta Diversity") + # edit axis and legend titles
    scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels
trait.jac

# plot both heatmaps next to each other
gridExtra::grid.arrange(trait.jac.turn, trait.jac.nest, ncol = 2)
```

Phylogenetic Signal

``` r
phylo_tree <- tree

phylo_tree$node.label<-NULL

check_diff <- setdiff(phylo_tree$tip.label, row.names(traits))
check_diff

#Determine phylogenetic signal for each continuous trait
ps_TL<-setNames(traits$MaxLengthTL, row.names(traits))
phylosig(phylo_tree, ps_TL, method = "lambda", test = TRUE)

multiPhylosignal(traits, phylo_tree, checkdata = T)

phylosignal(x = ps_TL, phylo_tree)

match.phylo.data(phylo_tree, main_traits)

str(main_traits)
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
    ##  [1] ggrepel_0.9.4     viridis_0.6.4     viridisLite_0.4.2 ggplot2_3.4.4    
    ##  [5] picante_1.8.2     nlme_3.1-164      vegan_2.6-4       lattice_0.22-5   
    ##  [9] permute_0.9-7     tidyr_1.3.0       phytools_2.0-3    maps_3.4.1.1     
    ## [13] ape_5.7-1         reshape2_1.4.4    stringr_1.5.1     dplyr_1.1.4      
    ## [17] knitr_1.45       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fastmatch_1.1-4         gtable_0.3.4            xfun_0.41              
    ##  [4] numDeriv_2016.8-1.1     quadprog_1.5-8          vctrs_0.6.5            
    ##  [7] tools_4.3.1             generics_0.1.3          tibble_3.2.1           
    ## [10] fansi_1.0.6             highr_0.10              cluster_2.1.6          
    ## [13] pkgconfig_2.0.3         Matrix_1.6-4            scatterplot3d_0.3-44   
    ## [16] lifecycle_1.0.4         farver_2.1.1            compiler_4.3.1         
    ## [19] textshaping_0.3.7       munsell_0.5.0           mnormt_2.1.1           
    ## [22] combinat_0.0-8          codetools_0.2-19        htmltools_0.5.7        
    ## [25] yaml_2.3.7              pillar_1.9.0            MASS_7.3-60            
    ## [28] clusterGeneration_1.3.8 iterators_1.0.14        foreach_1.5.2          
    ## [31] phangorn_2.11.1         tidyselect_1.2.0        digest_0.6.33          
    ## [34] stringi_1.8.2           purrr_1.0.2             labeling_0.4.3         
    ## [37] splines_4.3.1           fastmap_1.1.1           grid_4.3.1             
    ## [40] colorspace_2.1-0        expm_0.999-8            cli_3.6.1              
    ## [43] magrittr_2.0.3          optimParallel_1.0-2     utf8_1.2.4             
    ## [46] withr_2.5.2             scales_1.3.0            rmarkdown_2.25         
    ## [49] igraph_1.5.1            gridExtra_2.3           ragg_1.2.6             
    ## [52] coda_0.19-4             evaluate_0.23           doParallel_1.0.17      
    ## [55] mgcv_1.9-0              rlang_1.1.2             Rcpp_1.0.11            
    ## [58] glue_1.6.2              rstudioapi_0.15.0       R6_2.5.1               
    ## [61] plyr_1.8.9              systemfonts_1.0.5
