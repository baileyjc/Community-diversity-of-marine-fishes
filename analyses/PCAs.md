Environment and trait PCAs
================
Bailey Carlson

## R Markdown

# Load PCA libaries and files

``` r
# Load the knitr package if not already loaded
library(knitr)
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

# All sites environment PCA

``` r
### Environment ###
## All sites
# Extract numerical
# all_sites_env_numerical_vars <- env[-c(8:9,16,18,20), c("temperature_median", "salinity_median", "oxygen_median", "pH_median", "volume_m3_w_chemocline", "surface_area_m2", "distance_to_ocean_mean_m", "tidal_efficiency", "max_depth")]
all_sites_env_numerical_vars <- env[-c(8:9,16,18,20), c("C", "S", "O", "pH", "V", "A", "meanD", "Te", "M")]

# Standardize the numerical variables
all_sites_env_scaled_data <- scale(all_sites_env_numerical_vars)

# Run PCA
all_sites_env_pca <- prcomp(all_sites_env_scaled_data, center = TRUE, scale. = TRUE)
summary(all_sites_env_pca)
```

    ## Importance of components:
    ##                           PC1    PC2    PC3     PC4     PC5     PC6     PC7
    ## Standard deviation     2.0193 1.5474 1.0004 0.75852 0.72518 0.55773 0.27873
    ## Proportion of Variance 0.4531 0.2660 0.1112 0.06393 0.05843 0.03456 0.00863
    ## Cumulative Proportion  0.4531 0.7191 0.8303 0.89423 0.95266 0.98723 0.99586
    ##                            PC8     PC9
    ## Standard deviation     0.19275 0.01068
    ## Proportion of Variance 0.00413 0.00001
    ## Cumulative Proportion  0.99999 1.00000

``` r
# Biplot
biplot(all_sites_env_pca)
```

![](PCAs_files/figure-gfm/All%20sites%20env%20PCA-1.png)<!-- -->

``` r
# Scree plot
plot(all_sites_env_pca)
```

![](PCAs_files/figure-gfm/All%20sites%20env%20PCA-2.png)<!-- -->

``` r
# Pull out principal component data frame
all_sites_env_coords <- all_sites_env_pca$x

# Combine these two data frames but we only want the stratification column from env
all_sites_env_coords <- as.data.frame(cbind(all_sites_env_coords, env[-c(8:9,16,18,20),19]))

# Renamme column
all_sites_env_coords <- all_sites_env_coords %>%
  rename(Stratification = V10)

# Make principal components numerical
all_sites_env_coords[,c(1:9)] <- sapply(all_sites_env_coords[,c(1:9)], as.numeric)

# Organize stratification column
all_sites_env_coords$Stratification <- factor(all_sites_env_coords$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# Loadings of variables on PC1
all_sites_env_loadings <- as.data.frame(all_sites_env_pca$rotation)

# Plot PC1 against PC2
all_sites_env_pca_plot <- ggplot(data = all_sites_env_coords, mapping = aes(x = PC1, y = PC2, color = Stratification)) + 
  stat_ellipse(data = all_sites_env_coords, aes(x = PC1, y = PC2), level = 0.5) +
  geom_point(data = all_sites_env_coords, stat = 'identity',
    size = 6,
    alpha = 1) + 
  geom_text_repel(data = all_sites_env_coords, label = row.names(all_sites_env_coords), size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  # geom_segment(data = all_sites_env_loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
  #              size =1, alpha = 0.2, color = "#CD661D") +
  # geom_text(data = all_sites_env_loadings, aes(x = PC1, y = PC2), 
  #           color = "#CD661D", label = row.names(all_sites_env_loadings), size = 6) + 
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
labs(x= "PC1 45%", y= "PC2 27%", color = "Stratification", tag = "A")

all_sites_env_pca_plot <- all_sites_env_pca_plot + guides(color = guide_legend(override.aes = list(label = "")))
all_sites_env_pca_plot
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values (`geom_path()`).

![](PCAs_files/figure-gfm/All%20sites%20env%20PCA-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/pca/all_sites_env_pca_plot.png", all_sites_env_pca_plot, width = 8.25, height = 4.13)
```

    ## Too few points to calculate an ellipse

    ## Warning: Removed 1 row containing missing values (`geom_path()`).

``` r
# Loadings of variables on PC1
all_sites_env_loadings_on_PC1 <- all_sites_env_pca$rotation[, 1]

# Sort variables by absolute loading on PC1
all_sites_env_sorted_vars_l1 <- names(all_sites_env_loadings_on_PC1[order(-abs(all_sites_env_loadings_on_PC1))])

# Print sorted variables
print(all_sites_env_sorted_vars_l1)
```

    ## [1] "Te"    "meanD" "S"     "O"     "pH"    "V"     "A"     "C"     "M"

``` r
# Loadings of variables on PC2
all_sites_env_loadings_on_PC2 <- all_sites_env_pca$rotation[, 2]

# Sort variables by absolute loading on PC2
all_sites_env_sorted_vars_l2 <- names(all_sites_env_loadings_on_PC2[order(-abs(all_sites_env_loadings_on_PC2))])

# Print sorted variables
print(all_sites_env_sorted_vars_l2)
```

    ## [1] "M"     "A"     "V"     "C"     "meanD" "S"     "Te"    "pH"    "O"

# Marine lakes environment PCA

``` r
## Only marine lakes
# Extract numerical
# marine_lakes_env_numerical_vars <- env[-c(7:9,11,16,18:20), c("temperature_median", "salinity_median", "oxygen_median", "pH_median", "volume_m3_w_chemocline", "surface_area_m2", "distance_to_ocean_mean_m", "tidal_efficiency", "max_depth")]
marine_lakes_env_numerical_vars <- env[-c(7:9,11,16,18:20), c("C", "S", "O", "pH", "V", "A", "meanD", "Te", "M")]

# Standardize the numerical variables
marine_lakes_env_scaled_data <- scale(marine_lakes_env_numerical_vars)

# Run PCA
marine_lakes_env_pca <- prcomp(marine_lakes_env_scaled_data, center = TRUE, scale. = TRUE)
summary(marine_lakes_env_pca)
```

    ## Importance of components:
    ##                           PC1    PC2     PC3     PC4     PC5     PC6     PC7
    ## Standard deviation     2.1915 1.4006 0.90238 0.81176 0.66078 0.41391 0.30557
    ## Proportion of Variance 0.5336 0.2180 0.09048 0.07322 0.04851 0.01904 0.01037
    ## Cumulative Proportion  0.5336 0.7516 0.84205 0.91527 0.96379 0.98282 0.99320
    ##                            PC8     PC9
    ## Standard deviation     0.23135 0.08777
    ## Proportion of Variance 0.00595 0.00086
    ## Cumulative Proportion  0.99914 1.00000

``` r
# Biplot
biplot(marine_lakes_env_pca)
```

![](PCAs_files/figure-gfm/Marine%20lakes%20env%20PCA-1.png)<!-- -->

``` r
# Scree plot
plot(marine_lakes_env_pca)
```

![](PCAs_files/figure-gfm/Marine%20lakes%20env%20PCA-2.png)<!-- -->

``` r
# Pull out principal component data frame
marine_lakes_env_coords <- marine_lakes_env_pca$x

# Combine these two data frames but we only want the stratification column from env
marine_lakes_env_coords <- as.data.frame(cbind(marine_lakes_env_coords, env[-c(7:9,11,16,18:20),19]))

# Renamme column
marine_lakes_env_coords <- marine_lakes_env_coords %>%
  rename(Stratification = V10)

# Make principal components numerical
marine_lakes_env_coords[,c(1:9)] <- sapply(marine_lakes_env_coords[,c(1:9)], as.numeric)

# Organize stratification column
marine_lakes_env_coords$Stratification <- factor(marine_lakes_env_coords$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# Loadings of variables
marine_lakes_env_loadings <- as.data.frame(marine_lakes_env_pca$rotation)

# Plot PC1 against PC2
marine_lakes_env_pca_plot <- ggplot(data = marine_lakes_env_coords, mapping = aes(x = PC1, y = PC2, color = Stratification)) + 
  stat_ellipse(data = marine_lakes_env_coords, aes(x = PC1, y = PC2, color = Stratification), level = 0.5) +
  geom_point(data = marine_lakes_env_coords, stat = 'identity',
    size = 6,
    alpha = 1) + 
  geom_text_repel(data = marine_lakes_env_coords, label = row.names(marine_lakes_env_coords), size = 7, point.padding = 3) +
  scale_color_viridis(alpha = 1, begin = 0.6, end = 0.75, discrete = T, option = "G") +
  # geom_segment(data = marine_lakes_env_loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
  #              size =1, alpha = 0.2, color = "#CD661D") +
  # geom_text(data = marine_lakes_env_loadings, aes(x = PC1, y = PC2), 
  #           color = "#CD661D", label = row.names(marine_lakes_env_loadings), size = 6) + 
  theme_bw() +
  theme(text = element_text(size = 22), legend.text = element_text(size = 22),
    axis.text = element_text(size = 22, color = "black"),
    axis.line = element_line(color = "black"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) + 
labs(x= "PC1 53%", y= "PC2 22%", color = "Stratification", tag = "B")

marine_lakes_env_pca_plot <- marine_lakes_env_pca_plot + guides(color = guide_legend(override.aes = list(label = "")))
marine_lakes_env_pca_plot
```

![](PCAs_files/figure-gfm/Marine%20lakes%20env%20PCA-3.png)<!-- -->

``` r
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/pca/marine_lakes_env_pca_plot.png", marine_lakes_env_pca_plot, width = 8.25, height = 4.13)

# Loadings of variables on PC1
marine_lakes_env_loadings_on_PC1 <- marine_lakes_env_pca$rotation[, 1]

# Sort variables by absolute loading on PC1
marine_lakes_env_sorted_vars_l1 <- names(marine_lakes_env_loadings_on_PC1[order(-abs(marine_lakes_env_loadings_on_PC1))])

# Print sorted variables
print(marine_lakes_env_sorted_vars_l1)
```

    ## [1] "meanD" "V"     "A"     "Te"    "M"     "S"     "C"     "O"     "pH"

``` r
# Loadings of variables on PC1
marine_lakes_env_loadings_on_PC2 <- marine_lakes_env_pca$rotation[, 2]

# Sort variables by absolute loading on PC1
marine_lakes_env_sorted_vars_l2 <- names(marine_lakes_env_loadings_on_PC2[order(-abs(marine_lakes_env_loadings_on_PC2))])

# Print sorted variables
print(marine_lakes_env_sorted_vars_l2)
```

    ## [1] "pH"    "S"     "O"     "M"     "A"     "V"     "Te"    "C"     "meanD"

# All sites trait PCA

``` r
St <- Sites_at[,c(1:18,36:37)] %>%
  rename(B = BodyShapeI) %>%
  rename(O = OperculumPresent) %>%
  rename(DS = DorsalSpinesMean) %>%
  rename(L = MaxLengthTL) %>%
  rename(T = Troph) %>%
  rename(DMin = DepthMin) %>%
  rename(DMax = DepthMax) %>%
  rename(TMin = TempPrefMin) %>%
  rename(TMax = TempPrefMax) %>%
  rename(F = FeedingPath) %>%
  rename(EC = RepGuild1) %>%
  rename(ES = RepGuild2) %>%
  rename(P = ParentalCare) %>%
  rename(W = WaterPref) 
St$X <- row.names(St)

### Traits ###
## All sites
# Extract numerical and encoded categorical variables
# all_sites_trait_numerical_vars <- Sites_at[, c("DorsalSpinesMean", "Troph", "DepthMin", "DepthMax", "TempPrefMin", "TempPrefMax")]
all_sites_trait_numerical_vars <- St[, c("DS", "L", "T", "DMin", "DMax", "TMin", "TMax")]
all_sites_trait_numerical_vars$X <- row.names(all_sites_trait_numerical_vars)

# One-hot encode categorical variables
# all_sites_trait_encoded <- as.data.frame(model.matrix(~ BodyShapeI + RepGuild1 + RepGuild2, data = Sites_at))
all_sites_trait_encoded <- as.data.frame(model.matrix(~ B + EC + ES + P, data = St))
all_sites_trait_encoded$X <- row.names(all_sites_trait_encoded)

# Assuming 'BinaryVar' is a binary categorical variable with 2 levels
# binary_encoded_OperculumPresent <- model.matrix(~ OperculumPresent - 1, data = Sites_at)
# binary_encoded_OperculumPresent <- as.data.frame(model.matrix(~ O - 1, data = St))
# binary_encoded_OperculumPresent$X <- row.names(binary_encoded_OperculumPresent)

binary_encoded_FeedingPath <- model.matrix(~ FeedingPath - 1, data = Sites_at)
binary_encoded_FeedingPath <- as.data.frame(model.matrix(~ F - 1, data = St))
binary_encoded_FeedingPath$X <- row.names(binary_encoded_FeedingPath)

# Combine numerical and encoded categorical variables
# data_for_pca <- cbind(numerical_vars, trait_encoded)
all_sites_trait_data_for_pca <- full_join(all_sites_trait_numerical_vars, all_sites_trait_encoded, by = "X")
# all_sites_trait_data_for_pca <- full_join(all_sites_trait_data_for_pca, binary_encoded_OperculumPresent, by = "X")
all_sites_trait_data_for_pca <- full_join(all_sites_trait_data_for_pca, binary_encoded_FeedingPath, by = "X")
row.names(all_sites_trait_data_for_pca) <- all_sites_trait_data_for_pca$X
all_sites_trait_data_for_pca <- all_sites_trait_data_for_pca[,-c(8:9)]

# Standardize the variables
all_sites_trait_scaled_data <- scale(all_sites_trait_data_for_pca)

# Bring in site information
all_sites_trait_scaled_data <- as.data.frame(cbind(all_sites_trait_scaled_data[,-c(17)], St$Site))

# Convert scaled columns to numeric
all_sites_trait_scaled_data[,c(1:22)] <- sapply(all_sites_trait_scaled_data[,c(1:22)], as.numeric)

# Convert the 'V19' column to a factor
all_sites_trait_scaled_data$V23 <- as.factor(all_sites_trait_scaled_data$V23)

# Get rid of any species with an NA
all_sites_trait_scaled_data <- na.omit(all_sites_trait_scaled_data)

# Assuming 'group_column' is your grouping variable
all_sites_trait_data_splits <- split(all_sites_trait_scaled_data, all_sites_trait_scaled_data$V23)

# Subset the scaled data by stratification group
ocean_sites_trait_scaled_data <- subset(all_sites_trait_scaled_data, V23 %in% c("IBK", "LLNC", "NCN", "OLOO", "OTMO", "RCA"))
mixed_sites_trait_scaled_data <- subset(all_sites_trait_scaled_data, V23 %in% c("FLK", "HLO", "LLN", "MLN", "NLN", "NLU", "OLO", "ULN"))
stratified_sites_trait_scaled_data <- subset(all_sites_trait_scaled_data, V23 %in% c("BCM", "CLM", "GLK", "HLM", "NLK", "OTM", "SLN", "TLN"))

# Run PCA assuming 'group_column' is your grouping variable
all_sites_trait_pca <- prcomp(all_sites_trait_scaled_data[, -which(names(all_sites_trait_scaled_data) == "V23")], center = TRUE)
summary(all_sites_trait_pca)

ocean_sites_trait_pca <- prcomp(ocean_sites_trait_scaled_data[, -which(names(ocean_sites_trait_scaled_data) == "V23")], center = TRUE)
mixed_sites_trait_pca <- prcomp(mixed_sites_trait_scaled_data[, -which(names(mixed_sites_trait_scaled_data) == "V23")], center = TRUE)
stratified_sites_trait_pca <- prcomp(stratified_sites_trait_scaled_data[, -which(names(stratified_sites_trait_scaled_data) == "V23")], center = TRUE)

# Biplot
biplot(all_sites_trait_pca, xlabs = all_sites_trait_scaled_data$V23)

biplot(ocean_sites_trait_pca, xlabs = ocean_sites_trait_scaled_data$V23)
biplot(mixed_sites_trait_pca, xlabs = mixed_sites_trait_scaled_data$V23)
biplot(stratified_sites_trait_pca, xlabs = stratified_sites_trait_scaled_data$V23)

# Scree plot
plot(all_sites_trait_pca)

# Pull out principal component data frame
all_sites_trait_coords <- as.data.frame(all_sites_trait_pca$x)

# Make species column since the next step will delete row names
all_sites_trait_coords$X <- row.names(all_sites_trait_coords)

# Add stratification column
all_sites_trait_coords <- merge(all_sites_trait_coords, St[,c(19,21)], by = "X")

# Re-add row names of species
row.names(all_sites_trait_coords) <- all_sites_trait_coords$X

# Delete species column
all_sites_trait_coords <- all_sites_trait_coords[,-1]

# Make principal components numerical
all_sites_trait_coords[,c(1:22)] <- sapply(all_sites_trait_coords[,c(1:22)], as.numeric)

# Organize stratification column
all_sites_trait_coords$Stratification <- factor(all_sites_trait_coords$Stratification, levels = c("Ocean", "Mixed", "Stratified"))

# Loadings of variables
all_sites_trait_loadings <- as.data.frame(all_sites_trait_pca$rotation)

ocean_sites_trait_loadings <- as.data.frame(ocean_sites_trait_pca$rotation)
mixed_sites_trait_loadings <- as.data.frame(mixed_sites_trait_pca$rotation)
stratified_sites_trait_loadings <- as.data.frame(stratified_sites_trait_pca$rotation)

# Plot PC1 against PC2
all_sites_trait_pca_plot <- ggplot(data = all_sites_trait_coords, mapping = aes(x = PC1, y = PC2, color = Stratification)) + 
  stat_ellipse(data = all_sites_trait_coords, aes(x = PC1, y = PC2, color = Stratification), level = 0.5) +
  geom_point(data = all_sites_trait_coords, stat = 'identity',
    size = 6,
    alpha = 1) + 
  # geom_text(data = all_sites_trait_coords, label = row.names(all_sites_trait_coords), size = 5, vjust = -0.5) +
  scale_color_viridis(alpha = 1, begin = 0.45, end = 0.75, discrete = T, option = "G") +
  # geom_segment(data = all_sites_trait_loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
  #              size =1, alpha = 0.2, color = "#698B22") +
  # geom_text(data = all_sites_trait_loadings, aes(x = PC1, y = PC2), 
  #           color = "#698B22", label = row.names(all_sites_trait_loadings), size = 6) + 
  theme_bw() +
  theme(text = element_text(size = 18),
  plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
labs(x= "PC1 30%", y= "PC2 12%")
all_sites_trait_pca_plot
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/pca/all_sites_trait_pca_plot.png", all_sites_trait_pca_plot, width = 6, height = 4)

# Loadings of variables on PC1
all_sites_trait_loadings_on_PC1 <- all_sites_trait_pca$rotation[, 1]

ocean_sites_trait_loadings_on_PC1 <- ocean_sites_trait_pca$rotation[, 1]
mixed_sites_trait_loadings_on_PC1 <- mixed_sites_trait_pca$rotation[, 1]
stratified_sites_trait_loadings_on_PC1 <- stratified_sites_trait_pca$rotation[, 1]

# Sort variables by absolute loading on PC1
all_sites_trait_sorted_vars_l1 <- names(all_sites_trait_loadings_on_PC1[order(-abs(all_sites_trait_loadings_on_PC1))])

ocean_sites_trait_sorted_vars_l1 <- names(ocean_sites_trait_loadings_on_PC1[order(-abs(ocean_sites_trait_loadings_on_PC1))])
mixed_sites_trait_sorted_vars_l1 <- names(mixed_sites_trait_loadings_on_PC1[order(-abs(mixed_sites_trait_loadings_on_PC1))])
stratified_sites_trait_sorted_vars_l1 <- names(stratified_sites_trait_loadings_on_PC1[order(-abs(stratified_sites_trait_loadings_on_PC1))])

# Print sorted variables
print(all_sites_trait_sorted_vars_l1)

print(ocean_sites_trait_sorted_vars_l1)
print(mixed_sites_trait_sorted_vars_l1)
print(stratified_sites_trait_sorted_vars_l1)
```

# Marine lakes trait PCA

``` r
Lt <- Lakes_at[,c(1:18,36:37)] %>%
  rename(B = BodyShapeI) %>%
  rename(O = OperculumPresent) %>%
  rename(DS = DorsalSpinesMean) %>%
  rename(L = MaxLengthTL) %>%
  rename(T = Troph) %>%
  rename(DMin = DepthMin) %>%
  rename(DMax = DepthMax) %>%
  rename(TMin = TempPrefMin) %>%
  rename(TMax = TempPrefMax) %>%
  rename(F = FeedingPath) %>%
  rename(EC = RepGuild1) %>%
  rename(ES = RepGuild2) %>%
  rename(P = ParentalCare) %>%
  rename(W = WaterPref) 
Lt$X <- row.names(Lt)

### Traits ### 
## Marine lake sites -c(7,9,11,16,18:20)
# One-hot encode categorical variables
# marine_lakes_trait_numerical_vars <- Lakes_at[, c("DorsalSpinesMean", "Troph", "DepthMin", "DepthMax", "TempPrefMin", "TempPrefMax")]
marine_lakes_trait_numerical_vars <- Lt[, c("DS", "L", "T", "DMax", "TMax")]
marine_lakes_trait_numerical_vars$X <- row.names(marine_lakes_trait_numerical_vars)

# One-hot encode categorical variables
# marine_lakes_trait_encoded <- as.data.frame(model.matrix(~ BodyShapeI + RepGuild1 + RepGuild2, data = Lakes_at))
marine_lakes_trait_encoded <- as.data.frame(model.matrix(~ B + EC + ES + P, data = Lt))
marine_lakes_trait_encoded$X <- row.names(marine_lakes_trait_encoded)

# Assuming 'BinaryVar' is a binary categorical variable with 2 levels
# binary_encoded_OperculumPresent <- model.matrix(~ OperculumPresent - 1, data = Lakes_at)
# binary_encoded_OperculumPresent <- as.data.frame(model.matrix(~ O - 1, data = Lt))
# binary_encoded_OperculumPresent$X <- row.names(binary_encoded_OperculumPresent)

# binary_encoded_FeedingPath <- model.matrix(~ FeedingPath - 1, data = Lakes_at)
# binary_encoded_FeedingPath <- as.data.frame(model.matrix(~ F - 1, data = Lt))
# binary_encoded_FeedingPath$X <- row.names(binary_encoded_FeedingPath)

# Combine numerical and encoded categorical variables
# data_for_pca <- cbind(numerical_vars, trait_encoded)
marine_lakes_trait_data_for_pca <- full_join(marine_lakes_trait_numerical_vars, marine_lakes_trait_encoded, by = "X")
# marine_lakes_trait_data_for_pca <- full_join(marine_lakes_trait_data_for_pca, binary_encoded_OperculumPresent, by = "X")
# marine_lakes_trait_data_for_pca <- full_join(marine_lakes_trait_data_for_pca, binary_encoded_FeedingPath, by = "X")
row.names(marine_lakes_trait_data_for_pca) <- marine_lakes_trait_data_for_pca$X
marine_lakes_trait_data_for_pca <- marine_lakes_trait_data_for_pca[,-c(6:7)]

# Standardize the numerical variables
marine_lakes_trait_scaled_data <- scale(marine_lakes_trait_data_for_pca)

# Bring in site information
marine_lakes_trait_scaled_data <- as.data.frame(cbind(marine_lakes_trait_scaled_data[,-c(15)], Lt$Site))

# Convert scaled columns to numeric
marine_lakes_trait_scaled_data[,c(1:18)] <- sapply(marine_lakes_trait_scaled_data[,c(1:18)], as.numeric)

# Convert the 'V19' column to a factor
marine_lakes_trait_scaled_data$V19 <- as.factor(marine_lakes_trait_scaled_data$V19)

# Get rid of any species with an NA
marine_lakes_trait_scaled_data <- na.omit(marine_lakes_trait_scaled_data)

# Assuming 'group_column' is your grouping variable
marine_lakes_trait_data_splits <- split(marine_lakes_trait_scaled_data, marine_lakes_trait_scaled_data$V19)

# Run PCA assuming 'group_column' is your grouping variable
marine_lakes_trait_pca <- prcomp(marine_lakes_trait_scaled_data[, -which(names(marine_lakes_trait_scaled_data) == "V19")], center = TRUE)
summary(marine_lakes_trait_pca)

# Biplot
biplot(marine_lakes_trait_pca, xlabs = marine_lakes_trait_scaled_data$V19)

# Scree plot
plot(marine_lakes_trait_pca)

# Pull out principal component data frame
marine_lakes_trait_coords <- as.data.frame(marine_lakes_trait_pca$x)

# Make species column since the next step will delete row names
marine_lakes_trait_coords$X <- row.names(marine_lakes_trait_coords)

# Add stratification column
marine_lakes_trait_coords <- merge(marine_lakes_trait_coords, Lt[,c(19,21)], by = "X")

# Re-add row names of species
row.names(marine_lakes_trait_coords) <- marine_lakes_trait_coords$X

# Delete species column
marine_lakes_trait_coords <- marine_lakes_trait_coords[,-1]

# Make prinicipal components numerical
marine_lakes_trait_coords[,c(1:18)] <- sapply(marine_lakes_trait_coords[,c(1:18)], as.numeric)

# Organize stratification column
marine_lakes_trait_coords$Stratification <- factor(marine_lakes_trait_coords$Stratification, levels = c("Mixed", "Stratified"))

# Loadings of variables on PC1
marine_lakes_trait_loadings <- as.data.frame(marine_lakes_trait_pca$rotation)

# Plot PC1 against PC2
marine_lakes_trait_pca_plot <- ggplot(data = marine_lakes_trait_coords, mapping = aes(x = PC1, y = PC2, color = Stratification)) +
  stat_ellipse(data = marine_lakes_trait_coords, aes(x = PC1, y = PC2, color = Stratification), level = 0.5) +
  geom_point(data = marine_lakes_trait_coords, stat = 'identity',
    size = 6,
    alpha = 1) + 
  # geom_text(data = marine_lakes_trait_coords, label = row.names(marine_lakes_trait_coords), size = 5, vjust = -0.5) +
  scale_color_viridis(alpha = 1, begin = 0.6, end = 0.75, discrete = T, option = "G") +
  # geom_segment(data = marine_lakes_trait_loadings, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
  #              size =1, alpha = 0.2, color = "#698B22") +
  # geom_text(data = marine_lakes_trait_loadings, aes(x = PC1, y = PC2), 
  #           color = "#698B22", label = row.names(marine_lakes_trait_loadings), size = 6) + 
  theme_bw() +
  theme(text = element_text(size = 18),
  plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black")) +
labs(x= "PC1 36%", y= "PC2 13%")
marine_lakes_trait_pca_plot
ggsave("/Users/bailey/Documents/research/fish_biodiversity/figures/pca/marine_lakes_trait_pca_plot.png", marine_lakes_trait_pca_plot, width = 6, height = 4)

# Loadings of variables on PC1
marine_lakes_trait_loadings_on_PC1 <- marine_lakes_trait_pca$rotation[, 1]

# Sort variables by absolute loading on PC1
marine_lakes_trait_sorted_vars_l1 <- names(marine_lakes_trait_loadings_on_PC1[order(-abs(marine_lakes_trait_loadings_on_PC1))])

# Print sorted variables
print(marine_lakes_trait_sorted_vars_l1)
```

``` r
marine_lakes_trait_numerical_vars <- Lt[, c("DS", "L", "T", "DMin", "DMax", "TMin", "TMax")]
marine_lakes_trait_numerical_vars$X <- row.names(marine_lakes_trait_numerical_vars)

marine_lakes_trait_encoded <- as.data.frame(model.matrix(~ B + EC + ES + P + H, data = Lt))
marine_lakes_trait_encoded$X <- row.names(marine_lakes_trait_encoded)

marine_lakes_trait_data_for_pca <- full_join(marine_lakes_trait_numerical_vars, marine_lakes_trait_encoded, by = "X")
row.names(marine_lakes_trait_data_for_pca) <- marine_lakes_trait_data_for_pca$X
marine_lakes_trait_data_for_pca <- marine_lakes_trait_data_for_pca[,-c(8:9)]

marine_lakes_trait_scaled_data <- scale(marine_lakes_trait_data_for_pca)

marine_lakes_trait_scaled_data <- as.data.frame(cbind(marine_lakes_trait_scaled_data[,-c(17,24:26)], Lt$Site))

marine_lakes_trait_scaled_data[,c(1:22)] <- sapply(marine_lakes_trait_scaled_data[,c(1:22)], as.numeric)

marine_lakes_trait_scaled_data$V23 <- as.factor(marine_lakes_trait_scaled_data$V23)

marine_lakes_trait_scaled_data <- na.omit(marine_lakes_trait_scaled_data)

marine_lakes_trait_scaled_data <- marine_lakes_trait_scaled_data %>%
  rename(Site = V23)

aggregated_data <- marine_lakes_trait_scaled_data %>%
  group_by(Site) %>%
  summarize_all(mean, na.rm = TRUE)

marine_lakes_trait_data_splits <- split(aggregated_data, aggregated_data$Site)

marine_lakes_trait_pca <- prcomp(aggregated_data[, -which(names(aggregated_data) == "Site")], center = TRUE)
summary(marine_lakes_trait_pca)

biplot(marine_lakes_trait_pca, xlabs = aggregated_data$Site)

marine_lakes_trait_loadings_on_PC1 <- marine_lakes_trait_pca$rotation[, 1]

marine_lakes_trait_sorted_vars_l1 <- names(marine_lakes_trait_loadings_on_PC1[order(-abs(marine_lakes_trait_loadings_on_PC1))])

print(marine_lakes_trait_sorted_vars_l1)
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
    ##  [5] reshape2_1.4.4    stringr_1.5.1     ggrepel_0.9.4     viridis_0.6.4    
    ##  [9] viridisLite_0.4.2 ggplot2_3.4.4     dplyr_1.1.4       knitr_1.45       
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

# Extra PCAs

``` r
#Plot eigenvalues and percentages of variation of an ordination object
# Kaiser rule and broken stick model
# Usage:evplot(ev)
# License: GPL-2 
# Author: Francois Gillet, 25 August 2012
evplot <- function(ev)
{
  # Broken stick model (MacArthur 1957)
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}

#### Environment ####
env_chem_pca <- prcomp(env[-c(9,16,18,20),c(1,3,5,7,9)], center = TRUE, scale. = TRUE)
summary(env_chem_pca)
loadings_chem <- env_chem_pca$rotation
correlations_chem <- t(loadings_chem)*env_chem_pca$sdev
correlations_chem
ev_chem <- env_chem_pca$sdev^2
evplot(ev_chem)

env_phys_pca <- prcomp(env[-c(8,20),c(22,24,25,29,31)], center = TRUE, scale. = TRUE)
summary(env_phys_pca)
loadings_phys <- env_phys_pca$rotation
correlations_phys <- t(loadings_phys)*env_phys_pca$sdev
correlations_phys
ev_phys <- env_phys_pca$sdev^2
evplot(ev_phys)

SR_env_phys_pca <- prcomp(SR_env[-c(8,20),c(22,24,25,29,31)], center = TRUE, scale. = TRUE)
summary(env_phys_pca)
loadings_phys <- env_phys_pca$rotation
correlations_phys <- t(loadings_phys)*env_phys_pca$sdev
correlations_phys
ev_phys <- env_phys_pca$sdev^2
evplot(ev_phys)

#### Traits #####
test <- traits[,-c(1:2,10,13:16)]
test <- na.omit(test)
traits_pca <- prcomp(test, center = TRUE, scale. = TRUE)
summary(traits_pca)
loadings_traits <- traits_pca$rotation
correlations_traits <- t(loadings_traits)*traits_pca$sdev
correlations_traits
ev_traits <- traits_pca$sdev^2
evplot(ev_traits)

env_phys_pca <- prcomp(env[-c(8,20),c(22,24,25,29,31)], center = TRUE, scale. = TRUE)
summary(env_phys_pca)
loadings_phys <- env_phys_pca$rotation
correlations_phys <- t(loadings_phys)*env_phys_pca$sdev
correlations_phys
ev_phys <- env_phys_pca$sdev^2
evplot(ev_phys)

SR_env_phys_pca <- prcomp(SR_env[-c(8,20),c(22,24,25,29,31)], center = TRUE, scale. = TRUE)
summary(env_phys_pca)
loadings_phys <- env_phys_pca$rotation
correlations_phys <- t(loadings_phys)*env_phys_pca$sdev
correlations_phys
ev_phys <- env_phys_pca$sdev^2
evplot(ev_phys)
```
