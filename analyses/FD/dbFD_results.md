---
title: "FD analyses results"
output: github_document
editor_options: 
  chunk_output_type: console
---
#### R Markdown

## Load packages and files

```r
# Load the knitr package if not already loaded
library(knitr)

# Source the R Markdown file
knit("/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.Rmd", output = "/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md")
```

```
## 
## 
## processing file: /Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.Rmd
```

```
## 
  |                
  |          |   0%
  |                
  |          |   3%                                                                        
  |                
  |.         |   6% [Bringing everything together load modifying files packages]           
  |                
  |.         |   9%                                                                        
  |                
  |.         |  12% [Bringing everything together load in modifying files]                 
  |                
  |..        |  16%                                                                        
  |                
  |..        |  19% [Check species names across files]                                     
  |                
  |..        |  22%                                                                        
  |                
  |..        |  25% [Modify environment data]                                              
  |                
  |...       |  28%                                                                        
  |                
  |...       |  31% [Modify incidence matrices]                                            
  |                
  |...       |  34%                                                                        
  |                
  |....      |  38% [Modify phylogeny]                                                     
  |                
  |....      |  41%                                                                        
  |                
  |....      |  44% [Modify trait data]                                                    
  |                
  |.....     |  47%                                                                        
  |                
  |.....     |  50% [Modify community data frames]                                         
  |                
  |.....     |  53%                                                                        
  |                
  |......    |  56% [Modify community trait data]                                          
  |                
  |......    |  59%                                                                        
  |                
  |......    |  62% [Community trait data tests]                                           
  |                
  |.......   |  66%                                                                        
  |                
  |.......   |  69% [Modify stratification data frames]                                    
  |                
  |.......   |  72%                                                                        
  |                
  |........  |  75% [Modify stratification trait data]                                     
  |                
  |........  |  78%                                                                        
  |                
  |........  |  81% [Stratification trait data tests]                                      
  |                
  |........  |  84%                                                                        
  |                
  |......... |  88% [Modify site trait data frames]                                        
  |                
  |......... |  91%                                                                        
  |                
  |......... |  94% [Site trait data tests]                                                
  |                
  |..........|  97%                                                                        
  |                
  |..........| 100% [Bringing everything together load out modified files and session info]
```

```
## output file: /Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md
```

```
## [1] "/Users/bailey/Documents/research/fish_biodiversity/src/collection/load_collection_data.md"
```

```r
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
```

```
## Loading required package: viridisLite
```

```
## 
## Attaching package: 'viridis'
```

```
## The following object is masked from 'package:maps':
## 
##     unemp
```

```r
library(ggrepel)

#Analyses
library(vegan)
```

```
## Loading required package: permute
```

```
## Loading required package: lattice
```

```
## This is vegan 2.6-4
```

```
## 
## Attaching package: 'vegan'
```

```
## The following object is masked from 'package:phytools':
## 
##     scores
```

```r
library(FD)
```

```
## Loading required package: ade4
```

```
## Loading required package: geometry
```

```r
surveyed_sites <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "LCN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOO", "OTM", "OOM", "RCA", "SLN", "TLN", "ULN")

# Define your custom colors
custom_colors <- c("Reference" = "black", "Ocean" = "#EE6363", "Mixed" = "#87CEFA", "Stratified" = "#6E8B3D")
```

## Functional Diversity FRic, FEve, FDiv, FDis, CWM

```r
# All sites including reference
# Using categorical excludes variance eitehr by mode or meidan whereas turning them into numerical allows the use of the mean which provides greater population variation measurement but same general outcome
FD_total <- dbFD(traits, presabs_lake, w.abun = T, corr = "cailliez")
FRic_total <- FD_total$FRic
FEve_total <- FD_total$FEve
FDiv_total <- FD_total$FDiv
FDis_total <- FD_total$FDis
CWM_total <- FD_total$CWM
FD_total <- cbind(FRic_total, FEve_total, FDiv_total, FDis_total, CWM_total)

write.csv(FD_total, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/FD_total.csv")

# All sites including reference but categorical traits are made numeric
ntraits <- traits

ntraits$BodyShapeI <- as.numeric(ntraits$BodyShapeI)
ntraits$DemersPelag <- as.numeric(ntraits$DemersPelag)
ntraits$OperculumPresent <- as.numeric(ntraits$OperculumPresent)
ntraits$FeedingPath <- as.numeric(ntraits$FeedingPath)
ntraits$RepGuild1 <- as.numeric(ntraits$RepGuild1)
ntraits$RepGuild2 <- as.numeric(ntraits$RepGuild2)
ntraits$ParentalCare <- as.numeric(ntraits$ParentalCare)
ntraits$WaterPref <- as.numeric(ntraits$WaterPref)

FD_total_an <- dbFD(ntraits, presabs_lake, w.abun = T, corr = "cailliez")
FRic_total_an <- FD_total_an$FRic
FEve_total_an <- FD_total_an$FEve
FDiv_total_an <- FD_total_an$FDiv
FDis_total_an <- FD_total_an$FDis
CWM_total_an <- FD_total_an$CWM
FD_total_an <- cbind(FRic_total_an, FEve_total_an, FDiv_total_an, FDis_total_an, CWM_total_an)

write.csv(FD_total_an, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/FD_total_an.csv")

# Only surveyed sites
FD_sites <- dbFD(straits, surveyed_sites_lake, corr = "cailliez", calc.CWM = F)
# CWM_sites <- functcomp(straits, surveyed_sites_lake)
FRic_sites <- FD_sites$FRic
FEve_sites <- FD_sites$FEve
FDiv_sites <- FD_sites$FDiv
FDis_sites <- FD_sites$FDis
# CWM_sites <- FD_sites$CWM
FD_sites <- cbind(FRic_sites, FEve_sites, FDiv_sites, FDis_sites)

write.csv(FD_sites, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/FD_sites.csv")

# Surveyed sites and community types without abundance
FD_total_strat <- dbFD(traits, strat_presabs_lake, corr = "cailliez")
# CWM_sites <- functcomp(straits, surveyed_sites_strat)
FRic_total_strat <- FD_total_strat$FRic
FEve_total_strat <- FD_total_strat$FEve
FDiv_total_strat <- FD_total_strat$FDiv
FDis_total_strat <- FD_total_strat$FDis
# CWM_total_strat <- FD_total_strat$CWM
FD_total_strat <- cbind(FRic_total_strat, FEve_total_strat, FDiv_total_strat, FDis_total_strat)

write.csv(FD_total_strat, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/FD_total_strat.csv")

# Surveyed sites and community types with abundance in each
FD_sites_strat <- dbFD(straits, surveyed_sites_strat, w.abun = T, corr = "cailliez", calc.CWM = F)
# CWM_sites <- functcomp(straits, surveyed_sites_sites_strat)
FRic_sites_strat <- FD_sites_strat$FRic
FEve_sites_strat <- FD_sites_strat$FEve
FDiv_sites_strat <- FD_sites_strat$FDiv
FDis_sites_strat <- FD_sites_strat$FDis
# CWM_sites_strat <- FD_sites_strat$CWM
FD_sites_strat <- cbind(FRic_sites_strat, FEve_sites_strat, FDiv_sites_strat, FDis_sites_strat)

write.csv(FD_sites_strat, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/FD_sites_strat.csv")
```

## Edit FD files
- Read in Functional Diversity file and combine with env data

```r
FD_total <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/FD_total.csv")

FD_total_env <- merge(FD_total, env[,-c(33:47)], by = "X", sort = F)
FD_total_env$Community <- factor(FD_total_env$Community, levels = c("Reference", "Ocean", "Holomictic", "Meromictic"))
FD_total_env$Stratification <- factor(FD_total_env$Stratification, levels = c("Reference", "Ocean", "Mixed", "Stratified"))
row.names(FD_total_env) <- FD_total_env$X
FD_total_env$B <- FD_total_env$BodyShapeI
FD_total_env$O <- FD_total_env$OperculumPresent
FD_total_env$DS <- FD_total_env$DorsalSpinesMean
FD_total_env$L <- FD_total_env$MaxLengthTL
FD_total_env$T <- FD_total_env$Troph
FD_total_env$DMin <- FD_total_env$DepthMin
FD_total_env$DMax <- FD_total_env$DepthMax
FD_total_env$TMin <- FD_total_env$TempPrefMin
FD_total_env$TMax <- FD_total_env$TempPrefMax
FD_total_env$FP <- FD_total_env$FeedingPath
FD_total_env$EC <- FD_total_env$RepGuild1
FD_total_env$ES <- FD_total_env$RepGuild2
FD_total_env$P <- FD_total_env$ParentalCare
FD_total_env$W <- FD_total_env$WaterPref

FD_total_env_numerical <- FD_total_env[, c("DS", "L", "T", "DMin", "DMax", "TMin", "TMax")]
FD_total_env_numerical$X <- row.names(FD_total_env_numerical)

FD_total_env_categorical <- as.data.frame(model.matrix(~ B + O + FP + EC + ES + P + W, data = FD_total_env))
FD_total_env_categorical$X <- row.names(FD_total_env_categorical)

FD_total_env_for_scale <- full_join(FD_total_env_numerical, FD_total_env_categorical, by = "X")
row.names(FD_total_env_for_scale) <- FD_total_env_for_scale$X
FD_total_env_for_scale <- FD_total_env_for_scale[,-c(8:9)]

# Standardize the variables
scaled_FD_total_env <- as.data.frame(scale(FD_total_env_for_scale))


FD_total_n <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/FD_total_n.csv")

FD_total_n_env <- merge(FD_total_n, FD_total_env[,-c(2:20)], by = "X", sort = F)
row.names(FD_total_n) <- FD_total_n$X
row.names(FD_total_n_env) <- FD_total_n_env$X 


FD_total_sub <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/FD_total_sub.csv")

FD_total_sub_env <- merge(FD_total_sub, FD_total_env[,-c(2:5)], by = "X", sort = F)
FD_total_sub <- merge(FD_total_sub, FD_total[,-c(2:5)], by = "X", sort = F)
row.names(FD_total_sub) <- FD_total_sub$X
row.names(FD_total_sub_env) <- FD_total_sub_env$X 


FD_sites <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/FD_sites.csv")

FD_sites_env <- merge(FD_sites, FD_total_env[,-c(2:5)], by = "X", sort = F)
FD_sites <- merge(FD_sites, FD_total[,-c(2:5)], by = "X", sort = F)
row.names(FD_sites) <- FD_sites$X
row.names(FD_sites_env) <- FD_sites_env$X 


FD_total_strat <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/FD_total_strat.csv")


FD_sites_strat <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/FD_sites_strat.csv")
```

## Biplot of FRic and FDis

```r
FRic_FDis_plot <- ggplot(data = FD_total_env, mapping = aes(y = FDis_total, x = FRic_total, color = Stratification, fill = Stratification)) + 
  geom_point(stat = 'identity',
    size = 4,
    alpha = 1, aes(group = Stratification, color = Stratification, fill = Stratification)) + 
  geom_text_repel(data = FD_total_env, label = FD_total_env$X, size = 5, point.padding = 3) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_bw() +
  theme(text = element_text(size = 22), 
    legend.title = element_text(size = 16),  # Adjust legend title size
    legend.text = element_text(size = 16),  # Adjust legend text size
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  labs(x= "FRic", y= "FDisp", colour = "Site type", fill = "Site type")
  # ylim(c(0,1)) +
  # scale_x_continuous(breaks = c(-6, -4, -2, 0, 2), limits = c(-7,2)) +
  # annotate('rect', ymin=0, ymax=1, xmin=-2, xmax=2, alpha = 0.3, fill='grey')
FRic_FDis_plot <- FRic_FDis_plot + guides(color = guide_legend(override.aes = list(label = "")))
FRic_FDis_plot
```


![](stratification_NMDS_envfit_analyses_files/figure-gfm/Biplot of FRic and FDis-1.png)<!-- -->

```r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/FD_plots/FRic_FDis_plot.png", FRic_FDis_plot, width = 6, height = 4, units = "in")
```


```r
sessionInfo()
```

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
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] FD_1.0-12.3       geometry_0.4.7    ade4_1.7-22       vegan_2.6-4      
##  [5] lattice_0.22-5    permute_0.9-7     ggrepel_0.9.4     viridis_0.6.4    
##  [9] viridisLite_0.4.2 ggplot2_3.5.1     tidyr_1.3.0       phytools_2.0-3   
## [13] maps_3.4.1.1      ape_5.7-1         reshape2_1.4.4    stringr_1.5.1    
## [17] dplyr_1.1.4       knitr_1.45       
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.2.0        farver_2.1.1            optimParallel_1.0-2    
##  [4] fastmap_1.1.1           combinat_0.0-8          digest_0.6.33          
##  [7] lifecycle_1.0.4         cluster_2.1.6           magrittr_2.0.3         
## [10] compiler_4.3.1          rlang_1.1.2             tools_4.3.1            
## [13] igraph_1.5.1            utf8_1.2.4              yaml_2.3.7             
## [16] phangorn_2.11.1         clusterGeneration_1.3.8 labeling_0.4.3         
## [19] mnormt_2.1.1            scatterplot3d_0.3-44    plyr_1.8.9             
## [22] abind_1.4-5             expm_0.999-8            withr_2.5.2            
## [25] purrr_1.0.2             numDeriv_2016.8-1.1     grid_4.3.1             
## [28] fansi_1.0.6             colorspace_2.1-0        scales_1.3.0           
## [31] iterators_1.0.14        MASS_7.3-60             cli_3.6.1              
## [34] rmarkdown_2.25          ragg_1.2.6              generics_0.1.3         
## [37] rstudioapi_0.15.0       magic_1.6-1             splines_4.3.1          
## [40] vctrs_0.6.5             Matrix_1.6-4            systemfonts_1.0.5      
## [43] foreach_1.5.2           glue_1.6.2              codetools_0.2-19       
## [46] stringi_1.8.2           gtable_0.3.4            quadprog_1.5-8         
## [49] munsell_0.5.0           tibble_3.2.1            pillar_1.9.0           
## [52] htmltools_0.5.7         R6_2.5.1                textshaping_0.3.7      
## [55] doParallel_1.0.17       evaluate_0.23           highr_0.10             
## [58] Rcpp_1.0.11             fastmatch_1.1-4         coda_0.19-4            
## [61] gridExtra_2.3           nlme_3.1-164            mgcv_1.9-0             
## [64] xfun_0.41               pkgconfig_2.0.3
```

## Null models of functional diversity

```r
#Try to bootstrap or null model community weighted mean
null_sites <- lapply(1:1000, function(x){
  randomizeMatrix(surveyed_sites_lake, null.model = "trialswap")
})
  
FD_null_sites <- lapply(1:1000, function(x){
dbFD(straits, null_sites[[x]], w.abun = T, corr = "cailliez")
})


null_sums_abund <- lapply(1:1000, function(x){
  randomizeMatrix(surveyed_sites_strat_abund, null.model = "trialswap")
})
  
FD_null_sums_abund <- lapply(1:1000, function(x){
dbFD(straits, null_sums_abund[[x]], w.abun = T, corr = "cailliez")
})

#Reference included
null <- lapply(1:1000, function(x){
  randomizeMatrix(presabs_lake, null.model = "trialswap")
})
  
FD_null <- lapply(1:1000, function(x){
dbFD(traits, null[[x]], w.abun = T, corr = "cailliez")
})
```

## Community weighted mean

```r
#Extract CWM dataframe
CWM_sums_abund <- FD_sums_abund[,c(1,6:21)]

#Extract CWM_null from list
CWM_null_sums_abund <- lapply(1:1000, function(x){
as.data.frame(FD_null_sums_abund[[x]]$CWM)
})

CWM_null_sums_abund <- do.call("rbind", CWM_null_sums_abund)


CWM_null_sums_abund_nvalues <- apply(CWM_null_sums_abund, 1, function(x){
          c(null_mean = mean(x), null_sd = sd(x))})

CWM_null_sums_abund_nvalues <- lapply(CWM_null_sums_abund, function(x){
  c(unlist(c(apply(x[, -c(1,2,4,12:16)], 1, function(y){
    null_mean = mean(y)}))),
    unlist(c(apply(x[, -c(1,2,4,12:16)], 1, function(z){
    null_sd = sd(z)})))
)})

CWM_null_sums_abund_summary <- lapply(FD_CWM_null_sums_abund, function(x){
  c(
  unlist(c(
    sapply(x[, c(1,2,4,12:16)], function(y){
      out_y <- table(y)
      ifelse(length(out_y) > 0, names(which.max(out_y)), NA)
    }))),    
  apply(x[, -c(1,2,4,12:16)], 2, mean, na.rm=T)
    )
})  
  

CWM_null_sas <- data.frame(do.call("rbind", FD_CWM_null_sums_abund))

CWM_null_values <- t(CWM_null__sums_abund_values)

CWM_null_values <- as.data.frame(CWM_null_values)

CWM_null_values$lake_code <- row.names(CWM_null_values)

CWM <- merge(CWM, CWM_null_values, by = 'lake_code')

CWM$z_score <- ((CWM$observed_mean - CWM$null_mean)/CWM$null_sd)

row.names(CWM) <- CWM$lake_code

CWM <- CWM[,-1]

write.csv(CWM, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/CWM.csv")
```

## Functional Richness

```r
#Extract FRic_sites dataframe
FRic_sites <- FD_sites[,c(1:2)]

#Extract FRic_sites_null from list
FRic_null_sites <- lapply(1:1000, function(x){
as.data.frame(FD_null_sites[[x]]$FRic)
})

FRic_null_sites <- do.call("cbind", FRic_null_sites)

FRic_null_sites_values <- apply(FRic_null_sites, 1, function(x){
  c(null_mean = mean(x), null_sd = sd(x))
})

FRic_null_sites_values <- as.data.frame(t(FRic_null_sites_values))

FRic_null_sites_values$X <- row.names(FRic_null_sites_values)

FRic_sites <- merge(FRic_sites, FRic_null_sites_values, by = 'X')

FRic_sites$z_score <- ((FRic_sites$FRic_sites - FRic_sites$null_mean)/FRic_sites$null_sd)

row.names(FRic_sites) <- FRic_sites$X

FRic_sites <- FRic_sites[,-1]

write.csv(FRic_sites, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/traits/FRic_sites.csv")
```

## Functional Dispersion

```r
#Extract FD_sites dataframe
FDis_sites <- FD_sites[,c(1,5)]

FDis_null_sites <- lapply(1:10, function(x){
as.data.frame(FD_null_sites[[x]]$FDis)
})

FDis_null_sites <- do.call("cbind", FDis_null_sites)

FDis_null_sites_values <- apply(FDis_null_sites, 1, function(x){
  c(null_mean = mean(x), null_sd = sd(x))
})

FDis_null_sites_values <- as.data.frame(t(FDis_null_sites_values))

FDis_null_sites_values$X <- row.names(FDis_null_sites_values)

FDis_sites <- merge(FDis_sites, FDis_null_sites_values, by = 'X')

FDis_sites$z_score <- ((FDis_sites$FDis_sites - FDis_sites$null_mean)/FDis_sites$null_sd)

row.names(FDis_sites) <- FDis_sites$X

FDis_sites <- FDis_sites[,-1]

write.csv(FDis_sites, "/Users/bailey/Documents/research/fish_biodiversity/data/analyses/FDis_sites.csv")
```
