TLN oxygen concentration
================

#### R Markdown

## Load packages

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

## Read in file and create figure

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] viridis_0.6.4     viridisLite_0.4.2 ggplot2_3.5.1     dplyr_1.1.4      
    ## [5] openxlsx_4.2.5.2 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] gtable_0.3.4      compiler_4.3.1    tidyselect_1.2.0  Rcpp_1.0.11      
    ##  [5] zip_2.3.0         gridExtra_2.3     textshaping_0.3.7 systemfonts_1.0.5
    ##  [9] scales_1.3.0      yaml_2.3.7        fastmap_1.1.1     R6_2.5.1         
    ## [13] labeling_0.4.3    generics_0.1.3    knitr_1.45        tibble_3.2.1     
    ## [17] munsell_0.5.0     pillar_1.9.0      rlang_1.1.2       utf8_1.2.4       
    ## [21] stringi_1.8.2     xfun_0.41         cli_3.6.1         withr_2.5.2      
    ## [25] magrittr_2.0.3    digest_0.6.33     grid_4.3.1        rstudioapi_0.15.0
    ## [29] lifecycle_1.0.4   vctrs_0.6.5       evaluate_0.23     glue_1.6.2       
    ## [33] farver_2.1.1      ragg_1.2.6        fansi_1.0.6       colorspace_2.1-0 
    ## [37] rmarkdown_2.25    tools_4.3.1       pkgconfig_2.0.3   htmltools_0.5.7
