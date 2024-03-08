Stratification abundance weighted mean trait plots
================

### R Markdown

## Load packages

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

## Body Shape

``` r
Stratification_at_props_BodyShape <- Stratification_at %>%
  group_by(Strat, BodyShapeI) %>%
  summarise(count = sum(SiteSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
BodyShape_plot_weighted <- ggplot(Stratification_at_props_BodyShape, mapping = aes(x= Strat, y= proportion, color = "black", fill = Strat, shape = BodyShapeI)) + 
  geom_point(position = "identity", size = 6, aes(group = Strat)) +
  geom_line(position = "identity", aes(group = BodyShapeI), linewidth = 0.5, linetype = "dotted") +
  scale_shape_manual(name = "Body Shape", values = c(11, 21:25), labels = c("1o" = "other", "2s" = "short deep", "3f" = "fusiform", "4e" = "elongated", "5l" = "eel-like")) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +  
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(0,0.6)) +
  ylab("Proportion") +
  xlab("Stratification") +
  labs(colour = "Stratification", tag = "A")
BodyShape_plot_weighted
```

![](stratification_trait_plots_files/figure-gfm/Body%20Shape-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/BodyShape_plot_weighted.png", BodyShape_plot_weighted, width = 6, height = 4, units = "in")


# BodyShape_plot_weighted <- ggplot(Stratification_at_props_BodyShape, aes(fill=Stratification, y=proportion, x=BodyShapeI)) + 
#   geom_bar(position='dodge', stStratification_at='identity') +
#   scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +  
#   guides(fill = "none", color = "none") +
#   theme_bw() +
#   theme(text = element_text(size = 24),
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("Proportion") +
#   xlab("Body Shape")
# BodyShape_plot_weighted
# ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/BodyShape_plot_weightedbar.png", BodyShape_plot_weighted, width = 6, height = 4, units = "in")


# # Identify how many individuals have one of the trait factors for each Stratification
# BodyShape_counts <- table(Stratification_at$Stratification, row.names = Stratification_at$BodyShapeI)
# 
# # Divide each count by the total number of rows to find the proportion
# BodyShape_props <- BodyShape_counts/(rowSums(BodyShape_counts))
# BodyShape_props
# 
# BodyShape_props <- as.dStratification_ata.frame(BodyShape_props)
# 
# BodyShape_props$Var1 <- factor(BodyShape_props$Var1, levels = c("Reference", "Ocean", "Holomictic", "Meromictic"))
# 
# BodyShape_plot <- ggplot(BodyShape_props, mapping = aes(x= Var1, y= Freq, color = Var1, shape = row.names)) + 
#   geom_point(position = "identity", size = 5, aes(group = Var1)) +
#   scale_color_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("Proportion") +
#   xlab("BodyShape")
# BodyShape_plot
```

## DemersPelag

``` r
Stratification_at_props_DemersPelag <- Stratification_at %>%
  group_by(Strat, DemersPelag) %>%
  summarise(count = sum(SiteSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
DemersPelag_plot_weighted <- ggplot(Stratification_at_props_DemersPelag, mapping = aes(x= Strat, y= proportion, color = "black", fill = Strat, shape = DemersPelag)) + 
  geom_point(position = "identity", size = 6, aes(group = Strat)) +
  geom_line(position = "identity", aes(group = DemersPelag), linewidth = 0.5, linetype = "dotted") +
  scale_shape_manual(name = "Demersal Pelagic", values = c(21:25, 12:13), 
                     labels = c("1r" = "reef-associated", "2pn" = "pelagic-neritic", "3p" = "pelagic", "4po" = "pelagic-oceanic", "5d" = "demersal", '6bp' = 'benthopelagic', '7bd' = 'bathydemersal')) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +  
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(0,1)) +
  ylab("Proportion") +
  xlab("Stratification") +
  labs(colour = "Stratification")
DemersPelag_plot_weighted
```

![](stratification_trait_plots_files/figure-gfm/DemersPelag-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/DemersPelag_plot_weighted.png", DemersPelag_plot_weighted, width = 6, height = 4, units = "in")
```

## Operculum Present

``` r
Stratification_at_props_OperculumPresent <- Stratification_at %>%
  group_by(Strat, OperculumPresent) %>%
  summarise(count = sum(SiteSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
OperculumPresent_plot_weighted <- ggplot(Stratification_at_props_OperculumPresent, mapping = aes(x= Strat, y= proportion, color = "black", fill = Strat, shape = OperculumPresent)) + 
  geom_point(position = "identity", size = 6, aes(group = Strat)) +
  geom_line(position = "identity", aes(group = OperculumPresent), linewidth = 0.5, linetype = "dotted") +
  scale_shape_manual(name = "Operculum", values = c(21:22)) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(0,0.75)) +
  ylab("Proportion") +
  xlab("Stratification") +
  labs(colour = "Stratification", tag = "B")
OperculumPresent_plot_weighted
```

![](stratification_trait_plots_files/figure-gfm/Operculum%20Present-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/OperculumPresent_plot_weighted.png", OperculumPresent_plot_weighted, width = 6, height = 4, units = "in")


# # Identify how many individuals have one of the trait factors for each Stratification
# OperculumPresent_counts <- table(Stratification_at$Stratification, row.names = Stratification_at$OperculumPresent)
# 
# # Divide each count by the total number of rows to find the proportion
# OperculumPresent_props <- OperculumPresent_counts/(rowSums(OperculumPresent_counts))
# OperculumPresent_props
# 
# OperculumPresent_props <- as.dStratification_ata.frame(OperculumPresent_props)
# 
# OperculumPresent_props$Var1 <- factor(OperculumPresent_props$Var1, levels = c("Reference", "Ocean", "Holomictic", "Meromictic"))
# 
# OperculumPresent_plot <- ggplot(OperculumPresent_props, mapping = aes(x= Var1, y= Freq, color = Var1, shape = row.names)) + 
#   geom_point(position = "identity", size = 5, aes(group = Var1)) +
#   scale_color_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("Proportion") +
#   xlab("OperculumPresent")
# OperculumPresent_plot
```

## Dorsal Spines Mean

``` r
DorsalSpinesMean_plot_weighted <- ggplot(Stratification_at_weighted, mapping = aes(x= Strat, y= DorsalSpinesMean, color = "black", fill = Strat)) +
  geom_violin(alpha = 1, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 1,
            width = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylab("Dorsal Spines") +
  xlab("Stratification") +
  labs(colour = "Stratification", tag = "C")
DorsalSpinesMean_plot_weighted
```

    ## Warning: Removed 10 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 10 rows containing missing values (`geom_point()`).

![](stratification_trait_plots_files/figure-gfm/Dorsal%20Spines%20Mean-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 10 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 10 rows containing missing values (`geom_point()`). 
# Twelve values with NA
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/DorsalSpinesMean_plot_weighted.png", DorsalSpinesMean_plot_weighted, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 10 rows containing non-finite values (`stat_ydensity()`).
    ## Removed 10 rows containing missing values (`geom_point()`).

``` r
# DorsalSpinesMax_plot <- ggplot(Stratification_at, mapping = aes(x= Strat, y= DorsalSpinesMax, fill = Strat)) +
#   geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75)) +
#   scale_fill_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   geom_jitter(shape = 21,
#             size = 3,
#             alpha = 0.75,
#             width = 0.1) +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("DorsalSpinesMax")
# DorsalSpinesMax_plot
```

## Max Length (TL)

``` r
MaxLengthTL_plot_weighted <- ggplot(Stratification_at_weighted, mapping = aes(x= Strat, y= MaxLengthTL, color = "black", fill = Strat)) +
  geom_violin(alpha = 1, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 1,
            width = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(0,450)) +
  ylab("Length (cm)") +
  xlab("Stratification") +
  labs(colour = "Stratification")
MaxLengthTL_plot_weighted
```

    ## Warning: Removed 3 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

![](stratification_trait_plots_files/figure-gfm/Max%20Length%20(TL)-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 3 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 3 rows containing missing values (`geom_point()`). 
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/MaxLengthTL_plot_weighted.png", MaxLengthTL_plot_weighted, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 3 rows containing non-finite values (`stat_ydensity()`).
    ## Removed 3 rows containing missing values (`geom_point()`).

``` r
# MaxLengthTL_plot <- ggplot(Stratification_at, mapping = aes(x= Strat, y= MaxLengthTL, fill = Strat)) +
#   geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75)) +
#   scale_fill_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   geom_jitter(shape = 21,
#             size = 3,
#             alpha = 0.75,
#             width = 0.1) +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("MaxLengthTL")
# MaxLengthTL_plot
```

## Trophic level

``` r
Troph_plot_weighted <- ggplot(Stratification_at_weighted, mapping = aes(x= Strat, y= Troph, color = "black", fill = Strat)) +
  geom_violin(alpha = 1, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 1,
            width = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylab("Trophic Level") +
  xlab("Stratification") +
  labs(colour = "Stratification", tag = "D")
Troph_plot_weighted
```

![](stratification_trait_plots_files/figure-gfm/Trophic%20level-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/Troph_plot_weighted.png", Troph_plot_weighted, width = 6, height = 4, units = "in")


# Troph_plot <- ggplot(Stratification_at, mapping = aes(x= Strat, y= Troph, fill = Strat)) +
#   geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75)) +
#   scale_fill_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   geom_jitter(shape = 21,
#             size = 3,
#             alpha = 0.75,
#             width = 0.1) +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("Troph")
# Troph_plot
```

## Depth Minimum

``` r
DepthMin_plot_weighted <- ggplot(Stratification_at_weighted, mapping = aes(x= Strat, y= DepthMin, color = "black", fill = Strat)) +
  geom_violin(alpha = 1, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 1,
            width = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(-1,201)) +
  ylab("Depth Min (m)") +
  xlab("Stratification") +
  labs(colour = "Stratification", tag = "A")
DepthMin_plot_weighted
```

    ## Warning: Removed 8 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 8 rows containing missing values (`geom_point()`).

![](stratification_trait_plots_files/figure-gfm/Depth%20Minimum-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 8 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 8 rows containing missing values (`geom_point()`). 
# One value with NA
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/DepthMin_plot_weighted.png", DepthMin_plot_weighted, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 8 rows containing non-finite values (`stat_ydensity()`).
    ## Removed 8 rows containing missing values (`geom_point()`).

``` r
# DepthMin_plot <- ggplot(Stratification_at, mapping = aes(x= Strat, y= DepthMin, fill = Strat)) +
#   geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75)) +
#   scale_fill_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   geom_jitter(shape = 21,
#             size = 3,
#             alpha = 0.75,
#             width = 0.1) +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("DepthMin")
# DepthMin_plot
```

## Depth Maximum

``` r
DepthMax_plot_weighted <- ggplot(Stratification_at_weighted, mapping = aes(x= Strat, y= DepthMax, color = "black", fill = Strat)) +
  geom_violin(alpha = 1, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 1,
            width = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(-1,1001)) +
  ylab("Depth Max (m)") +
  xlab("Stratification") +
  labs(colour = "Stratification")
DepthMax_plot_weighted
```

    ## Warning: Removed 10 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 10 rows containing missing values (`geom_point()`).

![](stratification_trait_plots_files/figure-gfm/Depth%20Maximum-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 10 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 10 rows containing missing values (`geom_point()`). 
# One value with NA
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/DepthMax_plot_weighted.png", DepthMax_plot_weighted, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 10 rows containing non-finite values (`stat_ydensity()`).
    ## Removed 10 rows containing missing values (`geom_point()`).

``` r
# DepthMax_plot <- ggplot(Stratification_at, mapping = aes(x= Strat, y= DepthMax, fill = Strat)) +
#   geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75)) +
#   scale_fill_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   geom_jitter(shape = 21,
#             size = 3,
#             alpha = 0.75,
#             width = 0.1) +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("DepthMax")
# DepthMax_plot
```

## Temperature Preference Minimum

``` r
TempPrefMin_plot_weighted <- ggplot(Stratification_at_weighted, mapping = aes(x= Strat, y= TempPrefMin, color = "black", fill = Strat)) +
  geom_violin(alpha = 1, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 1,
            width = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylab("Temp Min (Cº)") +
  xlab("Stratification") +
  labs(colour = "Stratification")
TempPrefMin_plot_weighted
```

    ## Warning: Removed 9 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 9 rows containing missing values (`geom_point()`).

![](stratification_trait_plots_files/figure-gfm/Temperature%20Preference%20Minimum-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 9 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 9 rows containing missing values (`geom_point()`). 
# Two values with NA
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/TempPrefMin_plot_weighted.png", TempPrefMin_plot_weighted, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 9 rows containing non-finite values (`stat_ydensity()`).
    ## Removed 9 rows containing missing values (`geom_point()`).

``` r
# TempPrefMin_plot <- ggplot(Stratification_at, mapping = aes(x= Strat, y= TempPrefMin, fill = Strat)) +
#   geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75)) +
#   scale_fill_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   geom_jitter(shape = 21,
#             size = 3,
#             alpha = 0.75,
#             width = 0.1) +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("TempPrefMin")
# TempPrefMin_plot
```

## Temperature Preference Maximum

``` r
TempPrefMax_plot_weighted <- ggplot(Stratification_at_weighted, mapping = aes(x= Strat, y= TempPrefMax, color = "black", fill = Strat)) +
  geom_violin(alpha = 1, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 1,
            width = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(17.5,30)) +
  ylab("Temp Max (Cº)") +
  xlab("Stratification") +
  labs(colour = "Stratification")
TempPrefMax_plot_weighted
```

    ## Warning: Removed 11 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 11 rows containing missing values (`geom_point()`).

![](stratification_trait_plots_files/figure-gfm/Temperature%20Preference%20Maximum-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 11 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 11 rows containing missing values (`geom_point()`). 
# Two values with NA
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/TempPrefMax_plot_weighted.png", TempPrefMax_plot_weighted, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 11 rows containing non-finite values (`stat_ydensity()`).
    ## Removed 11 rows containing missing values (`geom_point()`).

``` r
# TempPrefMax_plot <- ggplot(Stratification_at, mapping = aes(x= Strat, y= TempPrefMax, fill = Strat)) +
#   geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75)) +
#   scale_fill_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   geom_jitter(shape = 21,
#             size = 3,
#             alpha = 0.75,
#             width = 0.1) +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("TempPrefMax")
# TempPrefMax_plot
```

## Weight

``` r
Weight_plot_weighted <- ggplot(Stratification_at_weighted, mapping = aes(x= Strat, y= Weight, color = "black", fill = Strat)) +
  geom_violin(alpha = 1, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 1,
            width = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(0,1000000)) +
  ylab("Weight (g)") +
  xlab("Stratification")
Weight_plot_weighted
```

    ## Warning: Removed 3 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 3 rows containing missing values (`geom_point()`).

![](stratification_trait_plots_files/figure-gfm/Weight-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 3 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 3 rows containing missing values (`geom_point()`). 
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/Weight_plot_weighted.png", Weight_plot_weighted, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 3 rows containing non-finite values (`stat_ydensity()`).
    ## Removed 3 rows containing missing values (`geom_point()`).

``` r
# Weight_plot <- ggplot(Stratification_at, mapping = aes(x= Strat, y= Weight, fill = Strat)) +
#   geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75)) +
#   scale_fill_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   geom_jitter(shape = 21,
#             size = 3,
#             alpha = 0.75,
#             width = 0.1) +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("Weight")
# Weight_plot
```

## Caudal Fin Length

``` r
CaudalFinLength_plot_weighted <- ggplot(Stratification_at_weighted, mapping = aes(x= Strat, y= CaudalFinLength, color = "black", fill = Strat)) +
  geom_violin(alpha = 1, draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  geom_jitter(shape = 21,
            size = 2,
            alpha = 1,
            width = 0.1) +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(-1,100)) +
  ylab("Caudal Fin Length (cm)") +
  xlab("Stratification")
CaudalFinLength_plot_weighted
```

    ## Warning: Removed 13 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 13 rows containing missing values (`geom_point()`).

![](stratification_trait_plots_files/figure-gfm/Caudal%20Fin%20Length-1.png)<!-- -->

``` r
# Warning messages:
# 1: Removed 13 rows containing non-finite values (`stat_ydensity()`). 
# 2: Removed 13 rows containing missing values (`geom_point()`). 
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/CaudalFinLength_plot_weighted.png", CaudalFinLength_plot_weighted, width = 6, height = 4, units = "in")
```

    ## Warning: Removed 13 rows containing non-finite values (`stat_ydensity()`).
    ## Removed 13 rows containing missing values (`geom_point()`).

``` r
# CaudalFinLength_plot <- ggplot(Stratification_at, mapping = aes(x= Strat, y= CaudalFinLength, fill = Strat)) +
#   geom_violin(alpha = 0.75, draw_quantiles = c(0.25, 0.5, 0.75)) +
#   scale_fill_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   geom_jitter(shape = 21,
#             size = 3,
#             alpha = 0.75,
#             width = 0.1) +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("CaudalFinLength")
# CaudalFinLength_plot
```

## Feeding Pathway

``` r
Stratification_at_props_FeedingPath <- Stratification_at %>%
  group_by(Strat, FeedingPath) %>%
  summarise(count = sum(SiteSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
Stratification_at_props_FeedingPath <- na.omit(Stratification_at_props_FeedingPath)

FeedingPath_plot_weighted <- ggplot(Stratification_at_props_FeedingPath, mapping = aes(x= Strat, y= proportion, color = "black", fill = Strat, shape = FeedingPath)) + 
  geom_point(position = "identity", size = 6, aes(group = Strat)) +
  geom_line(position = "identity", aes(group = FeedingPath), linewidth = 0.5, linetype = "dotted") +
  scale_shape_manual(name = "Diet Source", values = c(21:22), labels = c("b" = "benthic", "p" = "pelagic")) +
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = 0.85, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(0,1)) +
  ylab("Proportion") +
  xlab("Stratification")
FeedingPath_plot_weighted
```

![](stratification_trait_plots_files/figure-gfm/Feeding%20Pathway-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/FeedingPath_plot_weighted.png", FeedingPath_plot_weighted, width = 6, height = 4, units = "in")


# # Identify how many individuals have one of the trait factors for each Stratification
# FeedingPath_counts <- table(Stratification_at$Stratification, row.names = Stratification_at$FeedingPStratification_ath)
# 
# # Divide each count by the total number of rows to find the proportion
# FeedingPath_props <- FeedingPath_counts/(rowSums(FeedingPath_counts))
# FeedingPath_props
# 
# FeedingPath_props <- as.dStratification_ata.frame(FeedingPath_props)
# 
# FeedingPath_props$Var1 <- factor(FeedingPath_props$Var1, levels = c("Reference", "Ocean", "Holomictic", "Meromictic"))
# 
# FeedingPath_plot <- ggplot(FeedingPath_props, mapping = aes(x= Var1, y= Freq, color = Var1, shape = row.names)) + 
#   geom_point(position = "identity", size = 5, aes(group = Var1)) +
#   scale_color_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("Proportion") +
#   xlab("Stratification")
# FeedingPath_plot
```

## Reproductive Guild1

``` r
Stratification_at_props_RepGuild1 <- Stratification_at %>%
  group_by(Strat, RepGuild1) %>%
  summarise(count = sum(SiteSums)) %>%
  group_by(Strat) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
Stratification_at_props_RepGuild1 <- na.omit(Stratification_at_props_RepGuild1)

RepGuild1_plot_weighted <- ggplot(Stratification_at_props_RepGuild1, mapping = aes(x= Strat, y= proportion, color = "black", fill = Strat, shape = RepGuild1)) + 
  geom_point(position = "identity", size = 6, aes(group = Strat)) +
  geom_line(position = "identity", aes(group = RepGuild1), linewidth = 0.5, linetype = "dotted") +
  scale_shape_manual(name = "Egg Care", values = c(21:23), labels = c('2g' = 'guarders', '1b' = 'bearers', '3n' =  'nonguarders')) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(0,0.6)) +
  ylab("Proportion") +
  xlab("Stratification") +
  labs(colour = "Stratification", tag = "E")
RepGuild1_plot_weighted
```

![](stratification_trait_plots_files/figure-gfm/Reproductive%20Guild1-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/RepGuild1_plot_weighted.png", RepGuild1_plot_weighted, width = 6, height = 4, units = "in")


# # Identify how many individuals have one of the trait factors for each Stratification
# RepGuild1_counts <- table(Stratification_at$Stratification, row.names = Stratification_at$RepGuild1)
# 
# # Divide each count by the total number of rows to find the proportion
# RepGuild1_props <- RepGuild1_counts/(rowSums(RepGuild1_counts))
# RepGuild1_props
# 
# RepGuild1_props <- as.dStratification_ata.frame(RepGuild1_props)
# 
# RepGuild1_props$Var1 <- factor(RepGuild1_props$Var1, levels = c("Reference", "Ocean", "Holomictic", "Meromictic"))
# 
# RepGuild1_plot <- ggplot(RepGuild1_props, mapping = aes(x= Var1, y= Freq, color = Var1, shape = row.names)) + 
#   geom_point(position = "identity", size = 5, aes(group = Var1)) +
#   scale_color_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("Proportion") +
#   xlab("Stratification")
# RepGuild1_plot
```

## Reproductive Guild2

``` r
Stratification_at_props_RepGuild2 <- Stratification_at %>%
  group_by(Strat, RepGuild2) %>%
  summarise(count = sum(SiteSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
Stratification_at_props_RepGuild2 <- na.omit(Stratification_at_props_RepGuild2)

RepGuild2_plot_weighted <- ggplot(Stratification_at_props_RepGuild2, mapping = aes(x= Strat, y= proportion, color = "black", fill = Strat, shape = RepGuild2)) + 
  geom_point(position = "identity", size = 6, aes(group = Strat)) +
  geom_line(position = "identity", aes(group = RepGuild2), linewidth = 0.5, linetype = "dotted") +
  scale_shape_manual(name = "Egg Strategy", values = c(11, 21:25), labels = c('1ib' = 'live bearers', '6s' = 'egg scatterers', '3n' = 'nesters', '5h' = 'brood hiders', '4t' = 'clutch tenders', '2eb' = 'external brooders')) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +  
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(0,0.6)) +
  ylab("Proportion") +
  xlab("Stratification")
RepGuild2_plot_weighted
```

![](stratification_trait_plots_files/figure-gfm/Reproductive%20Guild2-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/RepGuild2_plot_weighted.png", RepGuild2_plot_weighted, width = 6, height = 4, units = "in")


# # Identify how many individuals have one of the trait factors for each Stratification
# RepGuild2_counts <- table(Stratification_at$Stratification, row.names = Stratification_at$RepGuild2)
# 
# # Divide each count by the total number of rows to find the proportion
# RepGuild2_props <- RepGuild2_counts/(rowSums(RepGuild2_counts))
# RepGuild2_props
# 
# RepGuild2_props <- as.dStratification_ata.frame(RepGuild2_props)
# 
# RepGuild2_props$Var1 <- factor(RepGuild2_props$Var1, levels = c("Reference", "Ocean", "Holomictic", "Meromictic"))
# 
# RepGuild2_plot <- ggplot(RepGuild2_props, mapping = aes(x= Var1, y= Freq, color = Var1, shape = row.names)) + 
#   geom_point(position = "identity", size = 5, aes(group = Var1)) +
#   scale_color_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("Proportion") +
#   xlab("Stratification")
# RepGuild2_plot
```

## Parental Care

``` r
Stratification_at_props_ParentalCare <- Stratification_at %>%
  group_by(Strat, ParentalCare) %>%
  summarise(count = sum(SiteSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
Stratification_at_props_ParentalCare <- na.omit(Stratification_at_props_ParentalCare)

ParentalCare_plot_weighted <- ggplot(Stratification_at_props_ParentalCare, mapping = aes(x= Strat, y= proportion, color = "black", fill = Strat, shape = ParentalCare)) + 
  geom_point(position = "identity", size = 6, aes(group = Strat)) +
  geom_line(position = "identity", aes(group = ParentalCare), linewidth = 0.5, linetype = "dotted") +
  scale_shape_manual(name = "Parental Care", values = c(21:25), labels = c('4n' = 'none', '3p' = 'paternal', '2m' = 'maternal', '1b' = 'biparental')) +
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +  
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(0,0.6)) +
  ylab("Proportion") +
  xlab("Stratification") +
  labs(colour = "Stratification", tag = "F")
ParentalCare_plot_weighted
```

![](stratification_trait_plots_files/figure-gfm/Parental%20Care-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/ParentalCare_plot_weighted.png", ParentalCare_plot_weighted, width = 6, height = 4, units = "in")


# # Identify how many individuals have one of the trait factors for each Stratification
# ParentalCare_counts <- table(Stratification_at$Stratification, row.names = Stratification_at$ParentalCare)
# 
# # Divide each count by the total number of rows to find the proportion
# ParentalCare_props <- ParentalCare_counts/(rowSums(ParentalCare_counts))
# ParentalCare_props
# 
# ParentalCare_props <- as.dStratification_ata.frame(ParentalCare_props)
# 
# ParentalCare_props$Var1 <- factor(ParentalCare_props$Var1, levels = c("Reference", "Ocean", "Holomictic", "Meromictic"))
# 
# ParentalCare_plot <- ggplot(ParentalCare_props, mapping = aes(x= Var1, y= Freq, color = Var1, shape = row.names)) + 
#   geom_point(position = "identity", size = 5, aes(group = Var1)) +
#   scale_color_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("Proportion") +
#   xlab("Stratification")
# ParentalCare_plot
```

## Water

Fresh/Brack/Salt

``` r
Stratification_at_props_Water <- Stratification_at %>%
  group_by(Strat, WaterPref) %>%
  summarise(count = sum(SiteSums)) %>%
  mutate(proportion = count/sum(count))
```

    ## `summarise()` has grouped output by 'Strat'. You can override using the
    ## `.groups` argument.

``` r
#Stratification_at_props_Water$Water <- factor(Stratification_at_props_Water$Water, levels = c("all", "fresh", "fresh-brack", "brack", "brack-salt", "salt"))

Water_plot_weighted <- ggplot(Stratification_at_props_Water, mapping = aes(x= Strat, y= proportion, color = "black", fill = Strat, shape = WaterPref)) + 
  geom_point(position = "identity", size = 6, aes(group = Strat)) +
  geom_line(position = "identity", aes(group = WaterPref), linewidth = 0.5, linetype = "dotted") +
  scale_shape_manual(name = "Water", values = c(21:23,11,24:25), labels = c('3a' = 'all', '1s' = 'salt', '2bs' = 'brack-salt', '4b' = 'brack', '5fb' = 'fresh-brack', '6f' = 'fresh')) + 
  scale_color_viridis(alpha = 1, begin = 0, end = .85, discrete = T, option = "G") +
  scale_fill_viridis(alpha = 1, begin = 0.3, end = .85, discrete = T, option = "G") +  
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
  ylim(c(0,1)) +
  ylab("Proportion") +
  xlab("Stratification") +
  labs(colour = "Stratification", tag = "B")
Water_plot_weighted
```

![](stratification_trait_plots_files/figure-gfm/Water-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/stratification_trait_plots/Water_plot_weighted.png", Water_plot_weighted, width = 6, height = 4, units = "in")


# # Identify how many individuals have one of the trait factors for each Stratification
# Habitat_counts <- table(Stratification_at$Stratification, row.names = Stratification_at$Habitat)
# 
# # Divide each count by the total number of rows to find the proportion
# Habitat_props <- Habitat_counts/(rowSums(Habitat_counts))
# Habitat_props
# 
# Habitat_props <- as.dStratification_ata.frame(Habitat_props)
# 
# Habitat_props$Var1 <- factor(Habitat_props$Var1, levels = c("Reference", "Ocean", "Holomictic", "Meromictic"))
# 
# Habitat_plot <- ggplot(Habitat_props, mapping = aes(x= Var1, y= Freq, color = Var1, shape = row.names)) + 
#   geom_point(position = "identity", size = 5, aes(group = Var1)) +
#   scale_color_viridis(alpha = 0.5, end = 0.75, discrete = T, option = "G") +
#   theme_bw() +
#   theme(
#     plot.background = element_blank(),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     panel.border = element_blank(),
#     axis.line = element_line(color = "black")) +
#   ylab("Proportion") +
#   xlab("Habitat")
# Habitat_plot
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
    ##  [5] reshape2_1.4.4    stringr_1.5.1     knitr_1.45        viridis_0.6.4    
    ##  [9] viridisLite_0.4.2 ggplot2_3.4.4     dplyr_1.1.4      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fastmatch_1.1-4         gtable_0.3.4            xfun_0.41              
    ##  [4] lattice_0.22-5          numDeriv_2016.8-1.1     quadprog_1.5-8         
    ##  [7] vctrs_0.6.5             tools_4.3.1             generics_0.1.3         
    ## [10] tibble_3.2.1            fansi_1.0.6             highr_0.10             
    ## [13] pkgconfig_2.0.3         Matrix_1.6-4            scatterplot3d_0.3-44   
    ## [16] lifecycle_1.0.4         compiler_4.3.1          farver_2.1.1           
    ## [19] textshaping_0.3.7       munsell_0.5.0           mnormt_2.1.1           
    ## [22] combinat_0.0-8          codetools_0.2-19        htmltools_0.5.7        
    ## [25] yaml_2.3.7              pillar_1.9.0            MASS_7.3-60            
    ## [28] clusterGeneration_1.3.8 iterators_1.0.14        foreach_1.5.2          
    ## [31] nlme_3.1-164            phangorn_2.11.1         tidyselect_1.2.0       
    ## [34] digest_0.6.33           stringi_1.8.2           purrr_1.0.2            
    ## [37] labeling_0.4.3          fastmap_1.1.1           grid_4.3.1             
    ## [40] colorspace_2.1-0        expm_0.999-8            cli_3.6.1              
    ## [43] magrittr_2.0.3          optimParallel_1.0-2     utf8_1.2.4             
    ## [46] withr_2.5.2             scales_1.3.0            rmarkdown_2.25         
    ## [49] igraph_1.5.1            gridExtra_2.3           ragg_1.2.6             
    ## [52] coda_0.19-4             evaluate_0.23           doParallel_1.0.17      
    ## [55] rlang_1.1.2             Rcpp_1.0.11             glue_1.6.2             
    ## [58] rstudioapi_0.15.0       R6_2.5.1                plyr_1.8.9             
    ## [61] systemfonts_1.0.5
