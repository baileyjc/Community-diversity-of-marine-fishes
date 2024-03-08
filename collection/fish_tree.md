Fish phylogeny
================

### R Markdown

## See tutorial to download initial tree

See Emily Jane McTavishes’ tutorial for how to download a date tree for
species in your list.
<https://github.com/McTavishLab/jupyter_OpenTree_tutorials/blob/master/notebooks/DatedTree_Bailey.ipynb>
Here is the output tree file form Emily Jane’s tutorial
labelled_dated_tree.tre Before writing in labelled_dated_tree.tre you
need to remove spaces between the genus and species of some of the
species in the file.

Below is a python code sample from the tutorial Emily Jane provides for
other trees which is on her github.

``` python
# #This is the python code from Emily Jane McTavish's OpenTree tutorial demo
# from opentree import OT
# 
# fi = open("../tutorial/main.csv").readlines()
# 
# ott_ids = set()
# 
# for lin in fi[1:]: #skip the header
#     lii = lin.split(',')#split on commas
#     ott_id = int(lii[2])#grab the opentree id
#     ott_ids.add(ott_id)#add to the set
# 
# 
# treefile = "GA_waterfowl.tre"
# #Get the synthetic tree from OpenTree
# output = OT.synth_induced_tree(ott_ids=list(ott_ids),  label_format='name')
# output.tree.write(path = treefile, schema = "newick")
# output.tree.print_plot(width=100)
```

## Load packages

``` r
library(phytools)
```

    ## Loading required package: ape

    ## Loading required package: maps

``` r
library(rfishbase)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:ape':
    ## 
    ##     where

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(stringr)
library(fuzzyjoin)
library(reticulate)
```

## Some old code

We won’t use this but this is another way of acquiring a tree from the
Open Tree of Life.

``` r
#### This is previous information from the other way I figured out how to download a phylo tree frm the Open Tree of life
### Tree help from Emily Jane, see https://github.com/OpenTreeOfLife/chronosynth
## Code used to generate tree with "chronosynth" is: ott278114
## python examples/synthpriordate.py --node_id 'ott212201' --method bladj --output_dir teleosts
## Teleostei
## python examples/synthpriordate.py --node_id 'ott278114' --output gnathostomata
## Gnasthostomata
# I had to manually delete apostrophes and put "Teleostei:0" after the last parentheses in the tree file before loading into R
# Get rid of additional information on tip.labels other than species name, will also get rid of subspecies names but FishBase does not use subspecies nomenclature yet based upon the several I have looked up.
#gnathostomata_tree$tip.label <- str_extract(gnathostomata_tree$tip.label, "[^_]+_[^_]+")
```

## Load and check tree

We are importing the tree we downloaded above using Emily Jane’s GitHub
tutorial. The species names on the Open Tree of Life differ slightly
from the FishBase names so we need to validate the species names on the
tree. There are a few species that FishBase does not recognize or fails
to correct when validating.

``` r
# Read in your tree to R and check the details
fish_tree <- read.tree("/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/phylo/labelled_dated_tree.tre")
fish_tree
```

    ## 
    ## Phylogenetic tree with 1631 tips and 2635 internal nodes.
    ## 
    ## Tip labels:
    ##   Plotosus_lineatus, Chanos_chanos, Amblygaster_sirm, Amblygaster_clupeoides, Amblygaster_leiogaster, Herklotsichthys_quadrimaculatus, ...
    ## Node labels:
    ##   ott278114, ott114656, ott114654, ott773483, ott285821, ott471203, ...
    ## 
    ## Rooted; includes branch lengths.

``` r
# Test if the tree is rooted, should come back as TRUE
is.rooted(fish_tree)
```

    ## [1] TRUE

``` r
# Test if the tree is ultrametric, should come back as FALSE
is.ultrametric(fish_tree)
```

    ## [1] FALSE

``` r
# Rename species that are missing the underscore or are at the subspecies level
fish_tree$tip.label[fish_tree$tip.label=="mrcaott896554ott896558"] = "Tylosurus_melanotus"
fish_tree$tip.label[fish_tree$tip.label=="Moolgarda_seheli"] = "Crenimugil_seheli"
fish_tree$tip.label[fish_tree$tip.label=="Yongeichthys_nebulosus"] = "Acentrogobius_nebulosus"
fish_tree$tip.label[fish_tree$tip.label=="Acentrogobius_chusanensis"] = "Ctenogobius_chusanensis"
fish_tree$tip.label[fish_tree$tip.label=="Antennarius analis"] = "Abantennarius_analis"
fish_tree$tip.label[fish_tree$tip.label=="Antennarius dorehensis"] = "Abantennarius_dorehensis"
fish_tree$tip.label[fish_tree$tip.label=="Diagramma_picta_picta"] = "Diagramma_pictum"
fish_tree$tip.label[fish_tree$tip.label=="Coranthus polyacanthus"] = "Amioides_polyacanthus"

# Check to determine if names were changed
fish_tree_labs <- fish_tree$tip.label
fish_tree_labs <- as.data.frame(fish_tree_labs)

# Need to make a new tree file to compared differences after FishBase validation
tree <- fish_tree

# Replace _ with a space
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Validate the tip label names that you can using FishBase
tree$tip.label <- validate_names(tree$tip.label)
```

    ## Joining with `by = join_by(Subfamily, GenCode, FamCode)`
    ## Joining with `by = join_by(FamCode)`
    ## Joining with `by = join_by(Order, Ordnum, Class, ClassNum)`
    ## Joining with `by = join_by(Class, ClassNum)`

``` r
# Replace space with an _ because later adding tips by genus doesn't seem to work without the _
tree$tip.label <- gsub(" ", "_", tree$tip.label)

# Find fish names that have been corrected by fish base located in tree file, should be the same number as fish_tree_species_names_only
tree_species_names_only <- setdiff(tree$tip.label, fish_tree$tip.label)
tree_species_names_only
```

    ##  [1] "Stegastes_lacrymatus"           "Plectroglyphidodon_fasciolatus"
    ##  [3] "Amphiprion_biaculeatus"         "Amblyglyphidodon_batunaorum"   
    ##  [5] "Pycnochromis_amboinensis"       "Pycnochromis_caudalis"         
    ##  [7] "Pycnochromis_agilis"            "Pycnochromis_retrofasciatus"   
    ##  [9] "Pycnochromis_atripes"           "Pycnochromis_margaritifer"     
    ## [11] "Pycnochromis_vanderbilti"       "Pycnochromis_acares"           
    ## [13] "Azurina_lepidolepis"            "Azurina_brevirostris"          
    ## [15] "Pycnochromis_delta"             "Pictichromis_paccagnellorum"   
    ## [17] "Discotrema_crinophilum"         "Osteomugil_engeli"             
    ## [19] "Osteomugil_perusii"             "Plicomugil_labiosus"           
    ## [21] "Doboatherina_duodecimalis"      "Ferdauia_orthogrammus"         
    ## [23] "Ostracion_cubicum"              "Ostracion_solorense"           
    ## [25] "Rhynchostracion_nasus"          "Abantennarius_analis"          
    ## [27] "Abantennarius_dorehensis"       "Abantennarius_rosaceus"        
    ## [29] "Abantennarius_coccineus"        "Abantennarius_nummifer"        
    ## [31] "Centropyge_loriculus"           "Leiognathus_equula"            
    ## [33] "Scorpaenopsis_oxycephalus"      "Chromileptes_altivelis"        
    ## [35] "Mirolabrichthys_pascalus"       "Pyronotanthias_lori"           
    ## [37] "Mirolabrichthys_tuka"           "Nemanthias_bartlettorum"       
    ## [39] "Nemanthias_bicolor"             "Pyronotanthias_flavoguttatus"  
    ## [41] "Pyronotanthias_smithvanizi"     "Nemanthias_dispar"             
    ## [43] "Pyronotanthias_parvirostris"    "Kyphosus_ocyurus"              
    ## [45] "Lophocampus_brevidorsalis"      "Lophocampus_retzii"            
    ## [47] "Microphis_brachyurus"           "Microphis_leiaspis"            
    ## [49] "Synchiropus_ocellatus"          "Synchiropus_moyeri"            
    ## [51] "Synchiropus_morrisoni"          "Pseudogobius_poicilosoma"      
    ## [53] "Fusigobius_gracilis"            "Vanderhorstia_phaeosticta"     
    ## [55] "Acentrogobius_signatus"         "Stonogobiops_xanthorhinicus"   
    ## [57] "Trimma_macrophthalmus"          "Trimma_hotsarihiense"          
    ## [59] "Eleotris_acanthopomus"          "Ostorhinchus_semilineatus"     
    ## [61] "Zoramia_leptacanthus"           "Cheilodipterus_isostigma"      
    ## [63] "Yarica_hyalosoma"               "Amioides_polyacanthus"         
    ## [65] "Encheliophis_boraborensis"      "Encheliophis_homei"            
    ## [67] "Lamnostoma_polyophthalmum"      "Lamnostoma_orientale"          
    ## [69] "Gymnothorax_rueppelliae"        "Rhina_ancylostomus"            
    ## [71] "Taeniurops_meyeni"              "Stegostoma_tigrinum"

``` r
# Find fish names that were changed or dropped by rfishbase in fish_tree, should be the same number as tree_species_names_only
fish_tree_species_names_only <- setdiff(fish_tree$tip.label, tree$tip.label)
fish_tree_species_names_only
```

    ##  [1] "Plectroglyphidodon_lacrymatus" "Stegastes_fasciolatus"        
    ##  [3] "Premnas_biaculeatus"           "Amblyglyphidodon_batunai"     
    ##  [5] "Chromis_amboinensis"           "Chromis_caudalis"             
    ##  [7] "Chromis_agilis"                "Chromis_retrofasciata"        
    ##  [9] "Chromis_atripes"               "Chromis_margaritifer"         
    ## [11] "Chromis_vanderbilti"           "Chromis_acares"               
    ## [13] "Chromis_lepidolepis"           "Chromis_brevirostris"         
    ## [15] "Chromis_delta"                 "Pictichromis_paccagnellae"    
    ## [17] "Discotrema_crinophila"         "Moolgarda_engeli"             
    ## [19] "Moolgarda_perusii"             "Oedalechilus_labiosus"        
    ## [21] "Atherinomorus_duodecimalis"    "Carangoides_orthogrammus"     
    ## [23] "Ostracion_cubicus"             "Ostracion_solorensis"         
    ## [25] "Ostracion_nasus"               "Antennarius_analis"           
    ## [27] "Antennarius_dorehensis"        "Antennatus_rosaceus"          
    ## [29] "Antennatus_coccineus"          "Antennatus_nummifer"          
    ## [31] "Centropyge_loricula"           "Leiognathus_equulus"          
    ## [33] "Scorpaenopsis_oxycephala"      "Cromileptes_altivelis"        
    ## [35] "Pseudanthias_pascalus"         "Pseudanthias_lori"            
    ## [37] "Pseudanthias_tuka"             "Pseudanthias_bartlettorum"    
    ## [39] "Pseudanthias_bicolor"          "Pseudanthias_flavoguttatus"   
    ## [41] "Pseudanthias_smithvanizi"      "Pseudanthias_dispar"          
    ## [43] "Pseudanthias_parvirostris"     "Sectator_ocyurus"             
    ## [45] "Microphis_brevidorsalis"       "Microphis_retzii"             
    ## [47] "Oostethus_brachyurus"          "Coelonotus_leiaspis"          
    ## [49] "Neosynchiropus_ocellatus"      "Neosynchiropus_moyeri"        
    ## [51] "Neosynchiropus_morrisoni"      "Pseudogobius_javanicus"       
    ## [53] "Coryphopterus_gracilis"        "Ctenogobiops_phaeostictus"    
    ## [55] "Amoya_signata"                 "Stonogobiops_xanthorhinica"   
    ## [57] "Trimma_macrophthalmum"         "Trimma_hotsarihiensis"        
    ## [59] "Eleotris_acanthopoma"          "Apogon_semilineatus"          
    ## [61] "Zoramia_leptacantha"           "Cheilodipterus_isostigmus"    
    ## [63] "Apogon_hyalosoma"              "Coranthus_polyacanthus"       
    ## [65] "Carapus_boraborensis"          "Carapus_homei"                
    ## [67] "Lamnostoma_polyophthalma"      "Lamnostoma_orientalis"        
    ## [69] "Gymnothorax_rueppellii"        "Rhina_ancylostoma"            
    ## [71] "Taeniura_meyeni"               "Stegostoma_fasciatum"

``` r
# Validate the names using FishBase, this allows you to see what the names were changed into and those that have NAs
fish_tree_species_names_only <- gsub("_", " ", fish_tree_species_names_only)
val_fish_tree_species_names_only <- validate_names(fish_tree_species_names_only)
```

    ## Joining with `by = join_by(Subfamily, GenCode, FamCode)`
    ## Joining with `by = join_by(FamCode)`
    ## Joining with `by = join_by(Order, Ordnum, Class, ClassNum)`
    ## Joining with `by = join_by(Class, ClassNum)`

``` r
fish_tree_species_names_only <- gsub(" ", "_", fish_tree_species_names_only)

# This will show us the names that were changed but also those that did not register on FishBase which we will edit
# use this to go back and edit the orignial tree file after loading it in. Anything with an NA needs to be looked up on FishBase and edited manual using the code above to alter species names.
compare_species_names <- as.data.frame(cbind(val_fish_tree_species_names_only, fish_tree_species_names_only))

# Check file to see names changed and are now part of the new tree file
species_names <- as.data.frame(cbind(tree$tip.label, fish_tree$tip.label))

# Check your tree details
tree
```

    ## 
    ## Phylogenetic tree with 1631 tips and 2635 internal nodes.
    ## 
    ## Tip labels:
    ##   Plotosus_lineatus, Chanos_chanos, Amblygaster_sirm, Amblygaster_clupeoides, Amblygaster_leiogaster, Herklotsichthys_quadrimaculatus, ...
    ## Node labels:
    ##   ott278114, ott114656, ott114654, ott773483, ott285821, ott471203, ...
    ## 
    ## Rooted; includes branch lengths.

``` r
# Should still fail to be ultrametric
is.ultrametric(tree)
```

    ## [1] FALSE

## Load in species list

We can use the species list to know which species were not found on the
Open Tree of Life and attach them to the tree.

``` r
# Load in list of species made up of reference pool and sites
# From fish_species_data.rmd
existing_species <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/species_lists/existing_species.csv")

# Get rid of extra # column
existing_species <- existing_species[,-1]

existing_species <- gsub(" ", "_", existing_species)

# existing_species[existing_species=="Platybelone_argalus"] <- "Platybelone_argalus_argalus"
# existing_species[existing_species=="Tylosurus_melanotus"] <- "Tylosurus_acus_melanotus"

diff_check <- setdiff(tree$tip.label, existing_species)
diff_check
```

    ## [1] "Cirrhilabrus_cyanopleura" "Mugilogobius_stigmaticus"
    ## [3] "Callogobius_liolepis"     "Gymnothorax_neglectus"

``` r
# These species names were changed after the tree was made "Cirrhilabrus cyanopleura" "Mugilogobius stigmaticus" "Callogobius liolepis" "Gymnothorax neglectus" so they will be dropped at the end before the final tree is made.

# Only use for looking up species name changes for those that are not validated by FishBase
# existing_species <- as.data.frame(existing_species)
```

## Old code

Previous code to validate species names that was helpful but a more
efficient method was developed.

``` r
# # Once you have found the matching synonyms in the tree in relation to the existing_species file you will load in this file. We will use this file to correct the names in the tree to match FishBase.
# conversionTable <- read.csv("/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/phylo/tax_con_corrected.csv")
# 
# ## Translate tree taxa to match FishBase taxonomy (involves renaming and dropping tips)
# # Determines whether species in column 1 are found in the tree. It should come back as FALSE
# all(conversionTable[,1] %in% fish_tree$tip.label)
# # We run this to change the names in the tree to match those in the "existing_species" file
# for (i in 1:nrow(conversionTable)) {
#   if (!is.na(conversionTable[i,1])) {
#     fish_tree$tip.label[which(fish_tree$tip.label == conversionTable[i,2])] <- conversionTable[i,1]
#   }
# }
# 
# tip_labs <- fish_tree$tip.label
# 
# # This should print out a shorter list than before. Here you are comparing the "existing_species" file with the tree. What remains are unmatched species still missing in the tree.
# um_DataOnly <- setdiff(existing_species, tip_labs)
# um_DataOnly
# 
# um_DataOnly <- setdiff(tip_labs, existing_species)
# um_DataOnly
```

## Jonathan Chang’s Fix \#3

This is code to make the tree ultrametric but it unfortunately still
fails
<https://jonathanchang.org/blog/three-ways-to-check-and-fix-ultrametric-phylogenies/>

``` r
# N <- Ntip(tree)
# root_node <- N + 1
# root_to_tip <- dist.nodes(tree)[1:N, root_node]
# 
# e1 <- tree$edge[, 1] # parent node
# e2 <- tree$edge[, 2] # child node
# EL <- tree$edge.length
# 
# ages <- numeric(N + tree$Nnode)
# 
# for (ii in seq_along(EL)) {
#      if (ages[e1[ii]] == 0) {
#          ages[e1[ii]] <- ages[e2[ii]] + EL[ii]
#      } else {
#          recorded_age <- ages[e1[ii]]
#          new_age <- ages[e2[ii]] + EL[ii]
#          if (recorded_age != new_age) {
#              cat(sprintf("node %i age %.6f != %.6f\n", e1[ii], recorded_age, new_age))
#              EL[ii] <- recorded_age - ages[e2[ii]]
#          }
#      }
#  }
# 
# tree$edge.length <- EL
# is.ultrametric(tree)
```

## Old code

Was used to make tree ultrametric but was having difficulties with the
code working properly. Now for the species in “existing_species” that
had a Genus name change and that Genus was also present in the tree or
for those species that did not have a synonym in the tree. Here we will
make the tree ultrametric to attach them. See Cadotte, M. W. 2015.
Phylogenetic diversity–ecosystem function relationships are insensitive
to phylogenetic edge lengths. Functional Ecology 29:718–723 \#From
traitDependent_functions.R written by Gio Rappacioulo

``` r
# check_and_fix_ultrametric <- function(phy){
#   
#   if (!is.ultrametric(phy)){
#       
#       vv <- vcv.phylo(phy)
#       dx <- diag(vv)
#       mxx <- max(dx) - dx
#       for (i in 1:length(mxx)){
#           phy$edge.length[phy$edge[,2] == i] <- phy$edge.length[phy$edge[,2] == i] + mxx[i]
#       }
#       if (!is.ultrametric(phy)){
#           stop("Ultrametric fix failed\n")
#       }   
#   }
#   
#   return(phy)
# }
```

## Make tree ultrametric and attach missing species

Here we use some commands from phytools to make our tree ultrametric.
Again the main goal for our analyses is that the tree topology is
conserved even if we lose the full branch lengths. The
add.species.to.genus command attaches species that share the same genus
already found in the tree to be attached there. We then find the node
closest to the tips we can attach species too whose genus is not already
in the tree which is around the subfamily and family taxonomic level.

``` r
# See which species are missing from the tree
missing_species_from_tree <- setdiff(existing_species, tree$tip.label)
missing_species_from_tree
```

    ##  [1] "Acanthurus_sp"                 "Ambassis_nalua"               
    ##  [3] "Amsichthys_knighti"            "Apogon_sp"                    
    ##  [5] "Apogon_susanae"                "Apogon_talboti"               
    ##  [7] "Asterropteryx_sp"              "Azurina_elerae"               
    ##  [9] "Bathygobius_sp"                "Bryaninops_sp"                
    ## [11] "Callogobius_sp"                "Callogobius_sp_B"             
    ## [13] "Calloplesiops_altivelis"       "Carcharhinus_amblyrhynchoides"
    ## [15] "Cercamia_sp"                   "Chromis_albomaculata"         
    ## [17] "Chromis_anadema"               "Chromis_scotochiloptera"      
    ## [19] "Chromis_sp"                    "Chrysiptera_ellenae"          
    ## [21] "Cirrhilabrus_ryukyuensis"      "Cirrhilabrus_sp"              
    ## [23] "Cryptocentrus_nanus"           "Cryptocentrus_sp"             
    ## [25] "Dischistodus_fasciatus"        "Drombus_sp"                   
    ## [27] "Dussumieria_sp"                "Epinephelus_sp"               
    ## [29] "Eviota_diyhritisma"            "Eviota_sp"                    
    ## [31] "Eviota_sp_2"                   "Fusigobius_sp"                
    ## [33] "Glossogobius_sandakanensis"    "Gobiodon_macrochir"           
    ## [35] "Gobiodon_sp"                   "Gobiopsis_aporia"             
    ## [37] "Gorgasia_sp"                   "Gymnocranius_microdon"        
    ## [39] "Liopropoma_sp"                 "Liopropoma_sp_B"              
    ## [41] "Lubbockichthys_sp"             "Manonichthys_polynemus"       
    ## [43] "Mugilogobius_platystoma"       "Myripristis_sp"               
    ## [45] "Neoglyphidodon_mitratus"       "Neopomacentrus_violascens"    
    ## [47] "Ophichthus_sp"                 "Opistognathus_solorensis"     
    ## [49] "Opistognathus_variabilis"      "Opistognathus_wassi"          
    ## [51] "Oxyurichthys_zeta"             "Pandaka_sp"                   
    ## [53] "Parioglossus_sp"               "Planiliza_melinoptera"        
    ## [55] "Platax_sp"                     "Plectranthias_sp"             
    ## [57] "Plesiops_corallicola"          "Plesiops_gracilis"            
    ## [59] "Plesiops_oxycephalus"          "Plesiops_verecundus"          
    ## [61] "Pomacentrus_albiaxillaris"     "Pomacentrus_flavoaxillaris"   
    ## [63] "Pomacentrus_micronesicus"      "Pomacentrus_simsiang"         
    ## [65] "Pomacentrus_sp"                "Pomacentrus_taeniometopon"    
    ## [67] "Pseudamia_amblyuroptera"       "Pseudamia_hayashii"           
    ## [69] "Pseudamia_zonata"              "Pseudochromis_marshallensis"  
    ## [71] "Pseudochromis_pylei"           "Pseudochromis_tapeinosoma"    
    ## [73] "Pseudojuloides_zeus"           "Pseudoplesiops_annae"         
    ## [75] "Pseudoplesiops_immaculatus"    "Pseudoplesiops_rosae"         
    ## [77] "Pseudoplesiops_typus"          "Pycnochromis_lineatus"        
    ## [79] "Redigobius_oyensi"             "Rhabdoblennius_sp"            
    ## [81] "Schindleria_sp"                "Sebastapistes_sp"             
    ## [83] "Silhouettea_sp"                "Soleichthys_sp"               
    ## [85] "Sphyraena_qenie"               "Steeneichthys_plesiopsus"     
    ## [87] "Stegastes_punctatus"           "Tetronarce_sp"                
    ## [89] "Tomiyamichthys_sp"             "Tosanoides_sp"                
    ## [91] "Trimma_sp"                     "Trimma_sp_2"                  
    ## [93] "Trimmatom_sp"                  "Valenciennea_yanoi"           
    ## [95] "Xenisthmus_sp"                 "Xenisthmus_sp_2"

``` r
#tree <- check_and_fix_ultrametric(tree)

# Use phytools to make tree ultrametric even though tree should be ultrametric it barely fails the test
tree <- force.ultrametric(tree, method = "extend")
```

    ## ***************************************************************
    ## *                          Note:                              *
    ## *    force.ultrametric does not include a formal method to    *
    ## *    ultrametricize a tree & should only be used to coerce    *
    ## *   a phylogeny that fails is.ultrametric due to rounding --  *
    ## *    not as a substitute for formal rate-smoothing methods.   *
    ## ***************************************************************

``` r
# Check that code above worked
is.ultrametric(tree)
```

    ## [1] TRUE

``` r
# Added species to the tree based upon genus, groups may not be monophyletic or species cannot be added so a warning will be issued.
for(i in 1:length(missing_species_from_tree)) tree<-add.species.to.genus(tree,missing_species_from_tree[i],where="root")
```

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Acanthurus may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Apogon may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Apogon may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Apogon may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Azurina may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Bryaninops may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Callogobius may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Callogobius may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Carcharhinus may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Chromis may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Chromis may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Chromis may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Chromis may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Chrysiptera may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Cirrhilabrus may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Cirrhilabrus may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Cryptocentrus may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Cryptocentrus may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Drombus may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Epinephelus may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Eviota may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Eviota may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Eviota may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Fusigobius may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Glossogobius may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Mugilogobius may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Ophichthus may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Oxyurichthys may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Pseudochromis may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Pseudochromis may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Pseudochromis may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Pycnochromis may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): could not match your species to a genus
    ##   check spelling, including case

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Tomiyamichthys may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Trimma may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Trimma may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

    ## Warning in add.species.to.genus(tree, missing_species_from_tree[i], where = "root"): Valenciennea may not be monophyletic
    ##   attaching to the most inclusive group containing members of this genus

``` r
# What remains should be species with no genus to attach to in the tree
missing_species_from_tree <- setdiff(existing_species, tree$tip.label)
missing_species_from_tree
```

    ##  [1] "Amsichthys_knighti"         "Calloplesiops_altivelis"   
    ##  [3] "Lubbockichthys_sp"          "Manonichthys_polynemus"    
    ##  [5] "Opistognathus_solorensis"   "Opistognathus_variabilis"  
    ##  [7] "Opistognathus_wassi"        "Pseudoplesiops_annae"      
    ##  [9] "Pseudoplesiops_immaculatus" "Pseudoplesiops_rosae"      
    ## [11] "Pseudoplesiops_typus"       "Rhabdoblennius_sp"         
    ## [13] "Schindleria_sp"             "Soleichthys_sp"            
    ## [15] "Steeneichthys_plesiopsus"   "Tetronarce_sp"

``` r
# Get rid of everything after the _ to get genera only
missing_genera_from_tree <- str_extract(missing_species_from_tree, "[^_]+")

# Create data frame with Genus and Species in separate columns
missing_species_from_tree <- as.data.frame(cbind(missing_species_from_tree, missing_genera_from_tree))
names(missing_species_from_tree)[1] <- "Species"
names(missing_species_from_tree)[2] <- "Genus"

# Download taxonomic information from FishBase because we will use this information to figure out where to attach the remaining species on the tree
missing_genera_from_tree <- rfishbase::load_taxa() %>% 
  filter(Genus %in% missing_genera_from_tree) %>%
  collect()
```

    ## Joining with `by = join_by(Subfamily, GenCode, FamCode)`
    ## Joining with `by = join_by(FamCode)`
    ## Joining with `by = join_by(Order, Ordnum, Class, ClassNum)`
    ## Joining with `by = join_by(Class, ClassNum)`

``` r
# Choose which classification levels you want to keep
missing_genera_from_tree <- missing_genera_from_tree[,-c(1:2,8)]

# Merge dataframes together
missing_from_tree <- merge(missing_genera_from_tree, missing_species_from_tree, by = "Genus", all.x = F, all.y = F)

# Remove duplicates
missing_from_tree <- missing_from_tree[!duplicated(missing_from_tree[c('Species')]), ]

# Place species column before Genus column
missing_from_tree <- missing_from_tree %>% relocate(Species, .before = Genus)

# You have to do this part manually to know where to attach each species
# For species still not present in the tree, use missing_from_tree to determine the family the genus is in or higher taxonomic levels. Using the Open Tree of Life (https://tree.opentreeoflife.org/opentree/argus/ottol@278114) you can find the node.label for the family or subfamily if available. Sometimes the subfamily or higher taxonomic level is present on Open Tree of Life but may not be present in our tree so you have to find the lowest taxonomic level that is present on Open Tree of Life and present in our tree. Then using the code below you can determine the node number:
which(tree$node.label=="ott195928")
```

    ## [1] 2566

``` r
# Genus_species, node id in synthetic tree from OTL, number of node.label of our tree
# Amsichthys_knighti, mrcaott26915ott295046, 379 #One step up from subfamily level
# Calloplesiops_altivelis, ott224901, 387 #Subfamily taxonomic level
# Lubbockichthys_sp, mrcaott5034ott26915, 379 #One step up from subfamily level
# Manonichthys_polynemus, mrcaott392ott5034, 94 #Three steps down from the order
# Opistognathus_solorensis, mrcaott5034ott259991, 278 #One step up from the family level since our tree lacks it
# Opistognathus_variabilis, mrcaott5034ott259991, 278 #One step up from the family level since our tree lacks it
# Opistognathus_wassi, mrcaott5034ott259991, 278 #One step up from the family level since our tree lacks it
# Pseudoplesiops_annae, mrcaott26915ott295046, 379 #One step up from the subfamily level
# Pseudoplesiops_immaculatus, mrcaott26915ott295046, 379 #One step up from the subfamily level
# Pseudoplesiops_rosae, mrcaott26915ott295046, 379 #One step up from the subfamily level
# Pseudoplesiops_typus, mrcaott26915ott295046, 379 #One step up from the subfamily level
# Rhabdoblennius_sp, mrcaott13624ott324193, 327 #Two steps up from genus level
# Schindleria_sp, mrcaott36673ott193341, 2146 #One step up from the family level
# Soleichthys_sp, mrcaott21417ott120430, 523 #Thre steps up from the genus level, two steps down from family
# Steeneichthys_plesiopsus, ott224901, 387 #Subfamily taxonomic level
# Tetronarce_sp, ott195928, 2566 #Four steps up from the genus level

# Using the node number you can attach the species to the tree at the lowest taxa classification as possible
tree <- bind.tip(tree, tip.label = "Amsichthys_knighti", where = 379)
tree <- bind.tip(tree, tip.label = "Calloplesiops_altivelis", where = 387)
tree <- bind.tip(tree, tip.label = "Lubbockichthys_sp", where = 379)
tree <- bind.tip(tree, tip.label = "Manonichthys_polynemus", where = 94)
tree <- bind.tip(tree, tip.label = "Opistognathus_solorensis", where = 278)
tree <- bind.tip(tree, tip.label = "Opistognathus_variabilis", where = 278)
tree <- bind.tip(tree, tip.label = "Opistognathus_wassi", where = 278)
tree <- bind.tip(tree, tip.label = "Pseudoplesiops_annae", where = 379)
tree <- bind.tip(tree, tip.label = "Pseudoplesiops_immaculatus", where = 379)
tree <- bind.tip(tree, tip.label = "Pseudoplesiops_rosae", where = 379)
tree <- bind.tip(tree, tip.label = "Pseudoplesiops_typus", where = 379)
tree <- bind.tip(tree, tip.label = "Rhabdoblennius_sp", where = 327)
tree <- bind.tip(tree, tip.label = "Schindleria_sp", where = 2146)
tree <- bind.tip(tree, tip.label = "Soleichthys_sp", where = 523)
tree <- bind.tip(tree, tip.label = "Steeneichthys_plesiopsus", where = 387)
tree <- bind.tip(tree, tip.label = "Tetronarce_sp", where = 2566)

# Do one last check to ensure all species from your list are present in the tree now
missing_species <- setdiff(existing_species, tree$tip.label)
missing_species
```

    ## character(0)

``` r
# Create tree with only species from Palau
palau_fish_tree <- drop.tip(tree, setdiff(tree$tip.label, existing_species), trim.internal = TRUE, rooted = is.rooted(tree))
```

## Load out files

``` r
# Change names in tree back to what we had, unfortunately you have to do it this way because these species are not validated correctly by FishBase unless the subspecies identifier is present
# palau_fish_tree$tip.label[palau_fish_tree$tip.label=="Platybelone_argalus_argalus"] <- "Platybelone_argalus"
# palau_fish_tree$tip.label[palau_fish_tree$tip.label=="Tylosurus_acus_melanotus"] <- "Tylosurus_melanotus"

# Check dataframe to ensure changes were added
pft_labs <- as.data.frame(palau_fish_tree$tip.label)

existing_species <- as.data.frame(existing_species)

## Load out files
# Phylogenetic tree of species present in marine lakes and ocean of Palau
write.tree(palau_fish_tree,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/phylo/palau_fish_tree.tre")

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
    ## [1] reticulate_1.34.0 fuzzyjoin_0.1.6   stringr_1.5.1     dplyr_1.1.4      
    ## [5] rfishbase_4.1.2   phytools_2.0-3    maps_3.4.1.1      ape_5.7-1        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fastmatch_1.1-4         xfun_0.41               lattice_0.22-5         
    ##  [4] tzdb_0.4.0              numDeriv_2016.8-1.1     quadprog_1.5-8         
    ##  [7] vctrs_0.6.5             tools_4.3.1             generics_0.1.3         
    ## [10] curl_5.2.0              parallel_4.3.1          tibble_3.2.1           
    ## [13] fansi_1.0.6             pkgconfig_2.0.3         contentid_0.0.18       
    ## [16] Matrix_1.6-4            dbplyr_2.4.0            scatterplot3d_0.3-44   
    ## [19] lifecycle_1.0.4         compiler_4.3.1          progress_1.2.3         
    ## [22] mnormt_2.1.1            combinat_0.0-8          codetools_0.2-19       
    ## [25] htmltools_0.5.7         yaml_2.3.7              pillar_1.9.0           
    ## [28] crayon_1.5.2            MASS_7.3-60             openssl_2.1.1          
    ## [31] cachem_1.0.8            clusterGeneration_1.3.8 iterators_1.0.14       
    ## [34] foreach_1.5.2           nlme_3.1-164            phangorn_2.11.1        
    ## [37] tidyselect_1.2.0        digest_0.6.33           stringi_1.8.2          
    ## [40] duckdb_0.9.2-1          purrr_1.0.2             fastmap_1.1.1          
    ## [43] grid_4.3.1              expm_0.999-8            cli_3.6.1              
    ## [46] magrittr_2.0.3          optimParallel_1.0-2     utf8_1.2.4             
    ## [49] withr_2.5.2             readr_2.1.4             prettyunits_1.2.0      
    ## [52] httr_1.4.7              rmarkdown_2.25          igraph_1.5.1           
    ## [55] askpass_1.2.0           png_0.1-8               hms_1.1.3              
    ## [58] memoise_2.0.1           coda_0.19-4             evaluate_0.23          
    ## [61] knitr_1.45              doParallel_1.0.17       rlang_1.1.2            
    ## [64] Rcpp_1.0.11             glue_1.6.2              DBI_1.1.3              
    ## [67] rstudioapi_0.15.0       jsonlite_1.8.8          R6_2.5.1               
    ## [70] fs_1.6.3
