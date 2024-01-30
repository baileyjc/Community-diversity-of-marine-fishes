TLN oxygen concentration
================
Bailey Carlson

## R Markdown

``` r
library(openxlsx)
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
oxy_con <- read.xlsx("/Users/bailey/Documents/research/fish_biodiversity/data/collection/environment/TLN/TLN_oxygen_concentrations.xlsx")

oxy_con$Year <- as.factor(oxy_con$Year)

oxy_con_plot <- ggplot(data = oxy_con, aes(y = Depth, x = Oxygen, color = Year, fill = Year)) + 
  geom_point(stat = 'identity', size = 4, alpha = 1) + 
  geom_line(aes(group = Year), orientation = "y", size = 1) +
  scale_color_viridis(alpha = 1, begin = 0.3, end = 0.75, discrete = TRUE, option = "A") +
  scale_fill_viridis(alpha = 0.1, begin = 0.3, end = 0.75, discrete = TRUE, option = "A") +
  scale_y_reverse() +
  theme_bw() +
  theme(
    text = element_text(size = 22),
    legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  labs(x = "Oxygen concentration (mg/L)", y = "Depth (m)", color = "Year")
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
oxy_con_plot
```

![](TLN_oxygen_concentration_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/oxy_con_plot.png", oxy_con_plot, width = 6, height = 4, units = "in")
```
