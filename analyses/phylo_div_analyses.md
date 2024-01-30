Phylogenetic diversity analyses
================
Bailey Carlson

## R Markdown

# Load in files

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

# Phylogenetic diversity MPD, MNTD, PD for all

``` r
#Mean pairwise differences
phylo_sesmpd <- ses.mpd(presabs_lake, cophenetic(tree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(phylo_sesmpd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/phylo_sesmpd.csv")


#Mean nearest taxon distance
phylo_sesmntd <- ses.mntd(presabs_lake, cophenetic(tree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(phylo_sesmntd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/phylo_sesmntd.csv")


#Faiths PD
phylo_sespd <- ses.pd(presabs_lake, tree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)

write.csv(phylo_sespd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/phylo_sespd.csv")
```

# Phylogenetic diversity MPD, MNTD, PD for sites

``` r
#Mean pairwise differences
sphylo_sesmpd <- ses.mpd(surveyed_sites_lake, cophenetic(stree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(sphylo_sesmpd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/sphylo_sesmpd.csv")


#Mean nearest taxon distance
sphylo_sesmntd <- ses.mntd(surveyed_sites_lake, cophenetic(stree), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999, iterations = 1000)

write.csv(sphylo_sesmntd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/sphylo_sesmntd.csv")


#Faiths PD
sphylo_sespd <- ses.pd(surveyed_sites_lake, stree, null.model = "taxa.labels", runs = 999, iterations = 1000, include.root = TRUE)

write.csv(sphylo_sespd, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/sphylo_sespd.csv")
```

# Read in PD mpd, mntd, pd

Read in Phylogenetic Diversity files and combine with env data

``` r
phylo_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/phylo_sesmpd.csv")

phylo_sesmpd_env <- merge(phylo_sesmpd, env, by = "X", sort = F)
phylo_sesmpd_env$measure <- "phylo_sesmpd"
phylo_sesmpd_env$Stratification <- factor(phylo_sesmpd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))


phylo_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/phylo_sesmntd.csv")

phylo_sesmntd_env <- merge(phylo_sesmntd, env, by = "X", sort = F)
phylo_sesmntd_env$measure <- "phylo_sesmntd"
phylo_sesmntd_env$Stratification <- factor(phylo_sesmntd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))


phylo_sespd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/phylo_sespd.csv")

phylo_sespd_env <- merge(phylo_sespd, env, by = "X", sort = F)
phylo_sespd_env$measure <- "phylo_sespd"
phylo_sespd_env$Stratification <- factor(phylo_sespd_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))


sphylo_sesmpd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/sphylo_sesmpd.csv")

sphylo_sesmpd_env <- merge(sphylo_sesmpd, env[-20,], by = "X", sort = F)
sphylo_sesmpd_env$measure <- "sphylo_sesmpd"
sphylo_sesmpd_env$Stratification <- factor(sphylo_sesmpd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))


sphylo_sesmntd <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/sphylo_sesmntd.csv")

sphylo_sesmntd_env <- merge(sphylo_sesmntd, env[-20,], by = "X", sort = F)
sphylo_sesmntd_env$measure <- "sphylo_sesmntd"
sphylo_sesmntd_env$Stratification <- factor(sphylo_sesmntd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))


sphylo_sespd <-read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/phylo/sphylo_sespd.csv")

sphylo_sespd_env <- merge(sphylo_sespd, env[-20,], by = "X", sort = F)
sphylo_sespd_env$measure <- "sphylo_sespd"
sphylo_sespd_env$Stratification <- factor(sphylo_sespd_env$Stratification, levels = c("Ocean", "Mixed", "Stratified"))
```

PD Linear Models

``` r
PD_model <- lm(log(pd.obs) ~ -1 + volume_m3_w_chemocline + surface_area_m2 + distance_to_ocean_mean_m + tidal_lag_time_minutes + perimeter_fromSat + max_depth + logArea, data = phylo_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(PD_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ -1 + volume_m3_w_chemocline + surface_area_m2 + 
    ##     distance_to_ocean_mean_m + tidal_lag_time_minutes + perimeter_fromSat + 
    ##     max_depth + logArea, data = phylo_sespd_env[c(1:2, 4, 5, 
    ##     12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##        1        2        4        5       12       17       21       22 
    ##  0.11034 -0.16446 -0.04254 -0.03401  0.11215 -0.03267  0.10740 -0.04244 
    ## 
    ## Coefficients:
    ##                            Estimate Std. Error t value Pr(>|t|)
    ## volume_m3_w_chemocline    1.340e-06  4.516e-07   2.966    0.207
    ## surface_area_m2          -2.678e-06  3.518e-05  -0.076    0.952
    ## distance_to_ocean_mean_m -1.798e-03  2.022e-03  -0.889    0.537
    ## tidal_lag_time_minutes    3.531e-03  5.113e-03   0.690    0.615
    ## perimeter_fromSat        -1.345e-03  1.711e-03  -0.787    0.576
    ## max_depth                -8.870e-02  4.252e-02  -2.086    0.285
    ## logArea                   9.575e-01  2.107e-01   4.546    0.138
    ## 
    ## Residual standard error: 0.263 on 1 degrees of freedom
    ## Multiple R-squared:  0.9998, Adjusted R-squared:  0.9985 
    ## F-statistic: 775.9 on 7 and 1 DF,  p-value: 0.02764

``` r
anova(PD_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                          Df  Sum Sq Mean Sq  F value  Pr(>F)  
    ## volume_m3_w_chemocline    1 198.004 198.004 2862.617 0.01190 *
    ## surface_area_m2           1  75.143  75.143 1086.375 0.01931 *
    ## distance_to_ocean_mean_m  1  63.878  63.878  923.503 0.02094 *
    ## tidal_lag_time_minutes    1  29.166  29.166  421.665 0.03098 *
    ## perimeter_fromSat         1   7.057   7.057  102.024 0.06282 .
    ## max_depth                 1   1.017   1.017   14.698 0.16243  
    ## logArea                   1   1.429   1.429   20.662 0.13786  
    ## Residuals                 1   0.069   0.069                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
PD_model <- lm(log(pd.obs) ~ -1 + temperature_median + salinity_median + oxygen_median + pH_median, data = phylo_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(PD_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ -1 + temperature_median + salinity_median + 
    ##     oxygen_median + pH_median, data = phylo_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##        1        2        4        5       12       17       21       22 
    ##  0.15285 -0.14069  0.20773 -0.01537 -0.12652 -0.16591 -0.06790  0.15652 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)  
    ## temperature_median -0.02824    0.06261  -0.451   0.6753  
    ## salinity_median     0.03080    0.04252   0.724   0.5089  
    ## oxygen_median      -0.08952    0.15416  -0.581   0.5926  
    ## pH_median           0.95750    0.40300   2.376   0.0763 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.1995 on 4 degrees of freedom
    ## Multiple R-squared:  0.9996, Adjusted R-squared:  0.9992 
    ## F-statistic:  2359 on 4 and 4 DF,  p-value: 5.383e-07

``` r
anova(PD_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                    Df Sum Sq Mean Sq   F value   Pr(>F)    
    ## temperature_median  1 374.25  374.25 9403.9088 6.78e-08 ***
    ## salinity_median     1   1.07    1.07   26.9282 0.006564 ** 
    ## oxygen_median       1   0.06    0.06    1.4792 0.290755    
    ## pH_median           1   0.22    0.22    5.6449 0.076329 .  
    ## Residuals           4   0.16    0.04                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Plot phylogenetic richness for all sites with more than 1 species

``` r
pr_model <- lm(log(pd.obs) ~ -1 + logArea, data = phylo_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(pr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ -1 + logArea, data = phylo_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -0.6329 -0.5931 -0.1973  0.5237  1.3923 
    ## 
    ## Coefficients:
    ##         Estimate Std. Error t value Pr(>|t|)    
    ## logArea  0.65999    0.02769   23.84 5.81e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.8082 on 7 degrees of freedom
    ## Multiple R-squared:  0.9878, Adjusted R-squared:  0.9861 
    ## F-statistic: 568.3 on 1 and 7 DF,  p-value: 5.812e-08

``` r
anova(pr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##           Df Sum Sq Mean Sq F value    Pr(>F)    
    ## logArea    1 371.19  371.19  568.27 5.812e-08 ***
    ## Residuals  7   4.57    0.65                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
pr_plot <- ggplot(data = phylo_sespd_env[-20,], mapping = aes(y = log(pd.obs), x = logArea, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = phylo_sespd_env[-20,1], size = 7, point.padding = 3) +
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
  labs(x="Log Area (m"^"2"~")", y="Log Phylo Richness")
pr_plot <- pr_plot + guides(color = guide_legend(override.aes = list(label = "")))
pr_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/pr_plot_area.png", plot = pr_plot, width = 6, height = 4, units = "in")

pr_model <- lm(log(pd.obs) ~ -1 + distance_to_ocean_mean_m, data = phylo_sespd_env[c(1:2,4,5,12,17,21:22),])
summary(pr_model)
```

    ## 
    ## Call:
    ## lm(formula = log(pd.obs) ~ -1 + distance_to_ocean_mean_m, data = phylo_sespd_env[c(1:2, 
    ##     4, 5, 12, 17, 21:22), ])
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6262 -0.5679  0.8395  1.8450  4.4606 
    ## 
    ## Coefficients:
    ##                          Estimate Std. Error t value Pr(>|t|)    
    ## distance_to_ocean_mean_m 0.022738   0.003031   7.502 0.000137 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 2.437 on 7 degrees of freedom
    ## Multiple R-squared:  0.8894, Adjusted R-squared:  0.8736 
    ## F-statistic: 56.28 on 1 and 7 DF,  p-value: 0.0001371

``` r
anova(pr_model)
```

    ## Analysis of Variance Table
    ## 
    ## Response: log(pd.obs)
    ##                          Df Sum Sq Mean Sq F value    Pr(>F)    
    ## distance_to_ocean_mean_m  1 334.19  334.19  56.276 0.0001371 ***
    ## Residuals                 7  41.57    5.94                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
pr_plot <- ggplot(data = phylo_sespd_env[-20,], mapping = aes(y = log(pd.obs), x = distance_to_ocean_mean_m, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1) + 
  geom_text_repel(label = phylo_sespd_env[-20,1], size = 7, point.padding = 3) +
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
  labs(x="  Mean Distance to the Ocean (m)", y="Log Phylo Richness")
pr_plot <- pr_plot + guides(color = guide_legend(override.aes = list(label = "")))
pr_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/pr_plot_dist.png", plot = pr_plot, width = 6, height = 4, units = "in")
```

# Phylo Z score vs p value plots

``` r
phylo_sesmpd_plot <- ggplot(data = phylo_sesmpd_env, mapping = aes(y = mpd.obs.p, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = phylo_sesmpd_env[,1], size = 7, point.padding = 3) +
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
  labs(y= "p-value", x= "z-score", colour = "Stratification", tag = "A") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
phylo_sesmpd_plot <- phylo_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
phylo_sesmpd_plot
```

    ## Warning: ggrepel: 17 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](phylo_div_analyses_files/figure-gfm/create%20PD%20plots-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/phylo_sesmpd_plot_z.png", phylo_sesmpd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 18 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
phylo_sesmntd_plot <- ggplot(data = phylo_sesmntd_env, mapping = aes(y = mntd.obs.p, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = phylo_sesmntd_env[,1], size = 7, point.padding = 3) +
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
  labs(y= "p-value", x= "z-score", colour = "Stratification", tag = "B") +
  ylim(c(0,1))  +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
phylo_sesmntd_plot <- phylo_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
phylo_sesmntd_plot
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1 rows containing missing values (`geom_text_repel()`).

    ## Warning: ggrepel: 16 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](phylo_div_analyses_files/figure-gfm/create%20PD%20plots-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/phylo_sesmntd_plot_z.png", phylo_sesmntd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1 rows containing missing values (`geom_text_repel()`).

    ## Warning: ggrepel: 17 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
phylo_sespd_plot <- ggplot(data = phylo_sespd_env, mapping = aes(y = pd.obs.p, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = phylo_sespd_env[,1], size = 7, point.padding = 3) +
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
  labs(y= "p-value", x= "z-score", colour = "Stratification", tag = "C") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
phylo_sespd_plot <- phylo_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
phylo_sespd_plot
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1 rows containing missing values (`geom_text_repel()`).

    ## Warning: ggrepel: 15 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](phylo_div_analyses_files/figure-gfm/create%20PD%20plots-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/phylo_sespd_plot_z.png", phylo_sespd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1 rows containing missing values (`geom_text_repel()`).

    ## Warning: ggrepel: 17 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

Phylo Z score vs obs value plots

``` r
phylo_sesmpd_plot <- ggplot(data = phylo_sesmpd_env, mapping = aes(y = mpd.obs, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
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
  ylim(c(0,300)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=300, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
phylo_sesmpd_plot <- phylo_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
phylo_sesmpd_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/phylo_sesmpd_plot_obs+z.png", phylo_sesmpd_plot, width = 6, height = 4, units = "in")

phylo_sesmntd_plot <- ggplot(data = phylo_sesmntd_env, mapping = aes(y = mntd.obs, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
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
  ylim(c(0,250))  +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=250, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
phylo_sesmntd_plot <- phylo_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
phylo_sesmntd_plot
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/phylo_sesmntd_plot_obs+z.png", phylo_sesmntd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

``` r
phylo_sespd_plot <- ggplot(data = phylo_sespd_env, mapping = aes(y = pd.obs, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
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
  labs(y= "Observed PD", x = "z-score") +
  ylim(c(0,6000)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=6000, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
phylo_sespd_plot <- phylo_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
phylo_sespd_plot
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/phylo_sespd_plot_obs+z.png", phylo_sespd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 1 rows containing missing values (`geom_point()`).

Sphylo Z score vs p value plots

``` r
sphylo_sesmpd_plot <- ggplot(data = sphylo_sesmpd_env, mapping = aes(y = mpd.obs.p, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = sphylo_sesmpd_env[,1], size = 7, point.padding = 3) +
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
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=-0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
sphylo_sesmpd_plot <- sphylo_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
sphylo_sesmpd_plot
```

    ## Warning: ggrepel: 8 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/sphylo_sesmpd_plot_z.png", sphylo_sesmpd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
sphylo_sesmntd_plot <- ggplot(data = sphylo_sesmntd_env, mapping = aes(y = mntd.obs.p, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = sphylo_sesmntd_env[,1], size = 7, point.padding = 3) +
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
  ylim(c(0,1))  +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
sphylo_sesmntd_plot <- sphylo_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
sphylo_sesmntd_plot
```

    ## Warning: ggrepel: 11 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/sphylo_sesmntd_plot_z.png", sphylo_sesmntd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 12 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
sphylo_sespd_plot <- ggplot(data = sphylo_sespd_env, mapping = aes(y = pd.obs.p, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  geom_text_repel(label = sphylo_sespd_env[,1], size = 7, point.padding = 3) +
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
  labs(y= "p-value", x = "z-score", colour = "Stratification", tag = "F") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
sphylo_sespd_plot <- sphylo_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
sphylo_sespd_plot
```

    ## Warning: ggrepel: 12 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/sphylo_sespd_plot_z.png", sphylo_sespd_plot, width = 6, height = 4, units = "in")
```

    ## Warning: ggrepel: 13 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

Sphylo Z score vs obs value plots

``` r
sphylo_sesmpd_plot <- ggplot(data = sphylo_sesmpd_env, mapping = aes(y = mpd.obs, x = mpd.obs.z, color = Stratification, fill = Stratification)) + 
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
  ylim(c(0,250)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=250, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
sphylo_sesmpd_plot <- sphylo_sesmpd_plot + guides(color = guide_legend(override.aes = list(label = "")))
sphylo_sesmpd_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/sphylo_sesmpd_plot_obs+z.png", sphylo_sesmpd_plot, width = 6, height = 4, units = "in")

sphylo_sesmntd_plot <- ggplot(data = sphylo_sesmntd_env, mapping = aes(y = mntd.obs, x = mntd.obs.z, color = Stratification, fill = Stratification)) + 
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
  ylim(c(0,250))  +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=250, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
sphylo_sesmntd_plot <- sphylo_sesmntd_plot + guides(color = guide_legend(override.aes = list(label = "")))
sphylo_sesmntd_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/sphylo_sesmntd_plot_obs+z.png", sphylo_sesmntd_plot, width = 6, height = 4, units = "in")

sphylo_sespd_plot <- ggplot(data = sphylo_sespd_env, mapping = aes(y = pd.obs, x = pd.obs.z, color = Stratification, fill = Stratification)) + 
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
  labs(y= "Observed PD", x = "z-score") +
  ylim(c(0,6000)) +
  scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,3)) +
  annotate('rect', ymin=0, ymax=6000, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
sphylo_sespd_plot <- sphylo_sespd_plot + guides(color = guide_legend(override.aes = list(label = "")))
sphylo_sespd_plot
```

![](phylo_div_analyses_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/PD_plots/sphylo_sespd_plot_obs+z.png", sphylo_sespd_plot, width = 6, height = 4, units = "in")
```

Regression of MPD

``` r
#MPD
phylo_sesmpd_env_chem <- phylo_sesmpd_env[-c(9, 16, 18, 20),]

man_phylo_sesmpd_env_chem <- manova(cbind(temperature_median, salinity_median, oxygen_median, pH_median) ~ mpd.obs.z, data = phylo_sesmpd_env_chem)
summary.aov(man_phylo_sesmpd_env_chem)

lm_phylo_sesmpd_env_chem <- lm(mpd.obs.z ~ temperature_median * salinity_median * oxygen_median * pH_median, data = phylo_sesmpd_env_chem)
summary(lm_phylo_sesmpd_env_chem)

step( object = lm_phylo_sesmpd_env_chem,     # start at the full model
       direction = "backward"   # allow it remove predictors but not add them
)

phylo_sesmpd_env_chem_model <- glmulti("mpd.obs.z", c("temperature_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = phylo_sesmpd_env, method = "h")
```

Regression of MPD

``` r
phylo_sesmpd_env_phys <- phylo_sesmpd_env[-c(20),]

man_phylo_sesmpd_env_phys <- manova(cbind(surface_area_m2, volume_m3_w_chemocline, distance_to_ocean_min_m, tidal_efficiency, depth, logArea) ~ mpd.obs.z, data = phylo_sesmpd_env)
summary.aov(man_phylo_sesmpd_env_phys)

lm_phylo_sesmpd_env_phys <- lm(mpd.obs.z ~ surface_area_m2 + volume_m3_w_chemocline + distance_to_ocean_min_m + tidal_efficiency + depth + logArea, data = phylo_sesmpd_env_phys)
summary(lm_phylo_sesmpd_env_phys)

step( object = lm_phylo_sesmpd_env_phys,     # start at the full model
       direction = "both"   # allow it remove predictors but not add them
)


phylo_sesmpd_env_phys_model <- glmulti("mpd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "depth", "logArea"), intercept = FALSE, level = 1, data = phylo_sesmpd_env, method = "h")
```

Plot MPD vs Oxygen Concentration

``` r
lm_phylo_sesmpd_te <- lm(mpd.obs.z ~tidal_efficiency, data = phylo_sesmpd_env_phys)
summary(lm_phylo_sesmpd_te)

phylo_sesmpd_tidal_efficiency_plot <- ggplot(data = phylo_sesmpd_env, mapping = aes(x = tidal_efficiency, y = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 12),    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Tidal Efficiency", y= "sesMPD") +
  ylim(c(-6,2))
phylo_sesmpd_tidal_efficiency_plot <- phylo_sesmpd_tidal_efficiency_plot + geom_smooth(mapping = aes(x = tidal_efficiency, y = mpd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
phylo_sesmpd_tidal_efficiency_plot
ggsave("phylo_sesmpd_tidal_efficiency_plot.png", phylo_sesmpd_tidal_efficiency_plot, width = 6, height = 4, units = "in")
```

Plot MPD vs Tidal efficiency

``` r
lm_phylo_sesmpd_te <- lm(mpd.obs.z ~tidal_efficiency, data = phylo_sesmpd_env_phys)
summary(lm_phylo_sesmpd_te)

phylo_sesmpd_tidal_efficiency_plot <- ggplot(data = phylo_sesmpd_env, mapping = aes(x = tidal_efficiency, y = mpd.obs.z, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 5,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = T, option = "G") +
  theme_bw() +
  theme(text = element_text(size = 12),    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Tidal Efficiency", y= "sesMPD") +
  ylim(c(-6,2))
phylo_sesmpd_tidal_efficiency_plot <- phylo_sesmpd_tidal_efficiency_plot + geom_smooth(mapping = aes(x = tidal_efficiency, y = mpd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
phylo_sesmpd_tidal_efficiency_plot
ggsave("phylo_sesmpd_tidal_efficiency_plot.png", phylo_sesmpd_tidal_efficiency_plot, width = 6, height = 4, units = "in")
```

Regression of MNTD

``` r
#MNTD
man_phylo_sesmntd_env_chem <- manova(cbind(temperature_median, salinity_median, oxygen_median, pH_median) ~ mntd.obs.z, data = phylo_sesmntd_env)
summary.aov(man_phylo_sesmntd_env_chem)

phylo_sesmntd_env_chem_model <- glmulti("mntd.obs.z", c("temperature_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = phylo_sesmntd_env, method = "h")

man_phylo_sesmntd_env_phys <- manova(cbind(surface_area_m2, volume_m3_w_chemocline, distance_to_ocean_min_m, tidal_efficiency, depth) ~ mntd.obs.z, data = phylo_sesmntd_env)
summary.aov(man_phylo_sesmntd_env_phys)

phylo_sesmntd_env_phys_model <- glmulti("mntd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "depth", "logArea"), intercept = FALSE, level = 1, data = phylo_sesmntd_env, method = "h")

#phylo_sesmntd_env_model <- glmulti("mntd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "conductivity_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = phylo_sesmntd_env, method = "h")

phylo_mntd_env_model <- lm(mntd.obs.z ~ -1 + volume_m3_w_chemocline + distance_to_ocean_min_m + tidal_efficiency + salinity_median, data = phylo_sesmntd_env)
summary(phylo_mntd_env_model)
anova(phylo_mntd_env_model)
```

Plot MNTD vs Oxygen Concentration

``` r
phylo_sesmntd_oxygen_plot <- ggplot(data = phylo_sesmntd_env, mapping = aes(x = oxygen_median, y = mntd.obs.z, color = Environment)) + 
  geom_point(stat = 'identity',
    size = 12,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = .3, end = .95, discrete = TRUE, option = "H") +
  theme_bw() +
  theme(text = element_text(size = 36),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Median Oxygen Concentration (mg/L)", y= "sesMNTD") +
  ylim(c(-7,1))
phylo_sesmntd_oxygen_plot <- phylo_sesmntd_oxygen_plot + geom_smooth(mapping = aes(x = oxygen_median, y = mntd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("phylo_sesmntd_oxygen_plot.png", phylo_sesmntd_oxygen_plot, width = 2, height = 2)
```

Plot MNTD vs Tidal Efficiency

``` r
phylo_sesmntd_tidal_efficiency_plot <- ggplot(data = phylo_sesmntd_env, mapping = aes(x = tidal_efficiency, y = mntd.obs.z, color = Environment)) + 
  geom_point(stat = 'identity',
    size = 12,
    alpha = 1) + 
  scale_color_viridis(alpha = 1, begin = .3, end = .95, discrete = TRUE, option = "H") +
  theme_bw() +
  theme(text = element_text(size = 36),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Tidal Efficiency", y= "sesMNTD") +
  ylim(c(-7,1))
phylo_sesmntd_tidal_efficiency_plot <- phylo_sesmntd_tidal_efficiency_plot + geom_smooth(mapping = aes(x = tidal_efficiency, y = mntd.obs.z), method = lm, se = FALSE, inherit.aes = FALSE, color = 'black')
ggsave("phylo_sesmntd_tidal_efficiency_plot.png", phylo_sesmntd_tidal_efficiency_plot, width = 2, height = 2)
```

Regression of PD

``` r
#PD
man_phylo_sespd_env_chem <- manova(cbind(temperature_median, salinity_median, oxygen_median, pH_median) ~ pd.obs.z, data = phylo_sespd_env)
summary.aov(man_phylo_sespd_env_chem)

phylo_sespd_env_chem_model <- glmulti("pd.obs.z", c("temperature_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = phylo_sespd_env, method = "h")

man_phylo_sespd_env_phys <- manova(cbind(surface_area_m2, volume_m3_w_chemocline, distance_to_ocean_min_m, tidal_efficiency, depth, logArea) ~ pd.obs.z, data = phylo_sespd_env)
summary.aov(man_phylo_sespd_env_phys)

phylo_sespd_env_phys_model <- glmulti("pd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "depth", "logArea"), intercept = FALSE, level = 1, data = phylo_sespd_env, method = "h")

#phylo_sespd_env_model <- glmulti("pd.obs.z", c("surface_area_m2", "volume_m3_w_chemocline", "distance_to_ocean_min_m", "tidal_efficiency", "conductivity_median", "salinity_median", "oxygen_median", "pH_median"), intercept = FALSE, level = 1, data = phylo_sespd_env, method = "h")

phylo_sespd_env_model <- lm(pd.obs.z ~ -1 + volume_m3_w_chemocline + distance_to_ocean_min_m + tidal_efficiency + salinity_median, data = phylo_sespd_env)
summary(phylo_pd_env_model)
anova(phylo_pd_env_model)
```

Plot PD vs Oxygen Concentration

``` r
phylo_sespd_env_plot <- ggplot(data = phylo_sespd_env, mapping = aes(x = env_variable, y = pd.obs.z, color = Environment, shape = Environment, label = X)) + 
  geom_point(stat = 'identity',
    size = 12,
    alpha = 1) + 
  geom_text_repel(min.segment.length = 2, show.legend = FALSE) +
  scale_color_viridis(alpha = 1, begin = .3, end = .95, discrete = TRUE, option = "H") +
  theme_bw() +
  theme(text = element_text(size = 36),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Env Variable", y= "sesFPD")
phylo_sespd_env_plot
ggsave("phylo_sespd_env_plot.png", phylo_sespd_env_plot, width = 2, height = 2)
```

Plot PD vs Tidal Efficiency

``` r
phylo_sespd_env_plot <- ggplot(data = phylo_sespd_env, mapping = aes(x = env_variable, y = pd.obs.z, color = Environment, shape = Environment, label = X)) + 
  geom_point(stat = 'identity',
    size = 12,
    alpha = 1) + 
  geom_text_repel(min.segment.length = 2, show.legend = FALSE) +
  scale_color_viridis(alpha = 1, begin = .3, end = .95, discrete = TRUE, option = "H") +
  theme_bw() +
  theme(text = element_text(size = 36),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
  labs(x= "Env Variable", y= "sesFPD")
phylo_sespd_env_plot
ggsave("phylo_sespd_env_plot.png", phylo_sespd_env_plot, width = 2, height = 2)
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

Beta Diversity of phylogenetic relatedness

``` r
# Do the following to order lakes by the circulation pattern, then delete env data
####https://pedrohbraga.github.io/CommunityPhylogenetics-Workshop/CommunityPhylogenetics-Workshop.html#between-assemblage-phylogenetic-structure####
presabs_env <- merge(pres_abs_by_lake, env_byLake, by = 'X')
presabs_env <- presabs_env[order(presabs_env$Environment),]
row.names(presabs_env) <- presabs_env$X
presabs_env_ordered <- presabs_env[,-c(1, 1692:1722)]

phylo_betap_jac <- phylo.beta.pair(presabs_env_ordered[-23,], phylo_tree, index.family = "jaccard")

# Turnover matrix:
phylo.jac.turn <- phylo_betap_jac$phylo.beta.jtu %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("Sites") + ylab("Sites") + labs(fill = "Turnover") + # edit axis and legend titles
  scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels

# Nestedness matrix:
phylo.jac.nest <- phylo_betap_jac$phylo.beta.jne %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("Sites") + ylab("Sites") + labs(fill = "Nestedness") + # edit axis and legend titles
    scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels

# Beta Diversity matrix:
phylo.jac <- phylo_betap_jac$phylo.beta.jac %>% 
  as.matrix() %>% melt() %>% 
  ggplot() + geom_tile(aes(x = Var1, y = Var2, fill = value)) + # create the heatmap
  xlab("Sites") + ylab("Sites") + labs(fill = "Beta Diversity") + # edit axis and legend titles
    scale_fill_viridis_c(option = "H", limits = range(0,1)) +
  theme(axis.text.x = element_text(angle = 90)) # rotates x axis labels
phylo.jac

# plot both heatmaps next to each other
phylo_beta_div_plot <- gridExtra::grid.arrange(phylo.jac.turn, phylo.jac.nest, ncol = 2)
ggsave("phylo_beta_div_plot.png", phylo_beta_div_plot, width = 6, height = 4.5)
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
