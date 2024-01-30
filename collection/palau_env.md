Palau lake environment
================
Bailey Carlson

## R Markdown

# Load in packages

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

# Load in and process files

Bring the environmental data into R. Missing environment data from Long
Lake and some of the ocean sites. We didn’t sample fish occurrences from
SLM or spooky lake. \#From trait-beta-div-processing.R

``` r
####################################
##### -- Environmental data -- #####
####################################
##### -- Load data -- #####
env_data <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/collection/environment/annual-environmental/annual-environmental-without_h2s_layer.csv", header = TRUE)

# Generate depth-by-lake environmental summaries
env_data_means <- env_data[, -c(2:4, 14)] %>%
  group_by(lake_code, depth) %>%
  summarise(across(.cols = everything(), ~mean(., na.rm = TRUE)))
```

    ## `summarise()` has grouped output by 'lake_code'. You can override using the
    ## `.groups` argument.

``` r
# Generate summaries of environmental variables
environment_by_lake <- env_data_means %>%
  dplyr::select(-depth) %>%
  group_by(lake_code) %>%
  summarise(across(.cols = everything(), list(median = ~median(., na.rm = TRUE), sd = ~sd(., na.rm = TRUE)), .names = "{.col}_{.fn}"))

# depth <- env_data_means %>%
#   dplyr::select(lake_code, depth) %>%
#   group_by(lake_code) %>%
#   summarise(across(.cols = depth, .fns = mean, .names = NULL, na.rm = TRUE))
# 
# depth[depth == "NaN"] <- NA
# 
# environment_byLake <- merge(environment_by_lake, depth, by = "lake_code", all = TRUE)

## Add additional lakes attributes
# Add measures of area and isolation
lake_attributes <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/collection/environment/lake_physical_attributes.csv", header = TRUE)

# Combine the physical and chemical factors for the lakes the by the site column
env_by_lake <- merge(environment_by_lake, lake_attributes, by = "lake_code", all = TRUE)

row.names(env_by_lake) <- env_by_lake$lake_code

env_by_lake <- env_by_lake[,-1]

# Reorder the columns based on the desired order
desired_order <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "LCN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOO", "OTM", "OOM", "RCA", "REF", "SLN", "TLN", "ULN")
env_by_lake <- env_by_lake[desired_order,]
```

# Load out files

``` r
# Contains environmental data for lakes sample from
write.csv(env_by_lake, "/Users/bailey/Documents/research/fish_biodiversity/data/collection/environment/env_by_lake.csv")

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] dplyr_1.1.4
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.33     utf8_1.2.4        R6_2.5.1          fastmap_1.1.1    
    ##  [5] tidyselect_1.2.0  xfun_0.41         magrittr_2.0.3    glue_1.6.2       
    ##  [9] tibble_3.2.1      knitr_1.45        pkgconfig_2.0.3   htmltools_0.5.7  
    ## [13] rmarkdown_2.25    generics_0.1.3    lifecycle_1.0.4   cli_3.6.1        
    ## [17] fansi_1.0.6       vctrs_0.6.5       withr_2.5.2       compiler_4.3.1   
    ## [21] rstudioapi_0.15.0 tools_4.3.1       pillar_1.9.0      evaluate_0.23    
    ## [25] yaml_2.3.7        rlang_1.1.2
