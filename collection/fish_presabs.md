Fish incidence
================

#### R Markdown

## Packages

- These are the packages that are needed to extract data from online
  repositories and manipulate the data

``` r
library(openxlsx) # data manipulation
library(rfishbase) # data manipulation
library(tidyr) # data manipulation
library(dplyr) # data manipulation
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
library(phytools) # data manipulation
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     where

    ## Loading required package: maps

``` r
library(taxize)
```

    ## 
    ## Attaching package: 'taxize'

    ## The following object is masked from 'package:rfishbase':
    ## 
    ##     synonyms

## Reference

- Take the Micronesian fish file and extract the Palau species from it

``` r
# File from Rob containing the reference pool of all Micronesian fish with a column dedicated to the fish species in Palau.
micronesian_records <- read.xlsx("/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/species_surveys/micronesian_fish_corrected.xlsx", sheet = 1, startRow = 3, colNames = T, cols =  1:31)

# Extract the three columns of Palauan fish data with genus and species names and whether a fish is present or absent there
palau_fish <- micronesian_records[c(1:2282),c(5:6,8,24,29)]

# Some fish lack the information of whether they occur in Palau. Therefore we need to get rid of any fish species that we are not confident occur in Palau. A 1 means that the species occurs there and anything else is speculation currently.

# Convert all NAs (and a few remaining "x") to 0
palau_fish$X29[palau_fish$X29 == "x"] <- NA; palau_fish$X29[palau_fish$X29 == "o"] <- 1; palau_fish$X29[palau_fish$X29 == "1?"] <- 1; palau_fish$X29[palau_fish$X29 == "?"] <- NA; palau_fish$X29[palau_fish$X29 == "e"] <- NA; palau_fish$X29[palau_fish$X29 == "x?"] <- NA; palau_fish$X29[palau_fish$X29 == "2?"] <- 1; palau_fish$X29[palau_fish$X29 == " "] <- NA; palau_fish$X29[palau_fish$X29 == "3"] <- 1

# Convert occurrence column (apart from species_name) to numeric
palau_fish$X29 <- as.numeric(palau_fish$X29)

# Remove all species that were labeled with an NA which means currently they are not found in Palau.
fish <- subset(palau_fish, X29 > 0.5) 

## Remove species that live below 200m
# Paracaesio stonei
# Gymnothorax neglectus
# Create a logical condition to identify rows to remove
rows_to_remove <- fish$X8 %in% c("Paracaesio stonei", "Gymnothorax neglectus")

# Subset the data frame to exclude rows meeting the condition
fish <- fish[!rows_to_remove, ]

## Combine Genus and Species
# This will be the corrected combined Genus and Species column where as the other column X8 contains the uncorrected names from Rob's list
fish$Species <- paste(fish[, 1], fish[, 2])

## To preserve the original names Rob gave his species I kept the full Genus and species name uncorrected but corrected the Genus or species individually so I could compare the ID Rob gave versus the corrected name that FishBase gave his ID
# This will show original species names and those that have changed due to using FishBase command validate_names which will be executed later on down below
check <- setdiff(fish$Species, fish$X8)
check
```

    ##   [1] "Stegostoma tigrinum"            "Sphyrna mokarran"              
    ##   [3] "Rhina ancylostomus"             "Neotrygon kuhlii"              
    ##   [5] "Anguilla bicolor"               "Gymnothorax richardsonii"      
    ##   [7] "Gymnothorax rueppelliae"        "Gymnothorax pictus"            
    ##   [9] "Evips percinctus"               "Ichthyapus vulturis"           
    ##  [11] "Lamnostoma orientale"           "Lamnostoma polyophthalmum"     
    ##  [13] "Encrasicholina devisi"          "Thryssa baelama"               
    ##  [15] "Thryssa hamiltonii"             "Thryssa mystax"                
    ##  [17] "Abantennarius analis"           "Abantennarius coccineus"       
    ##  [19] "Abantennarius dorehensis"       "Abantennarius nummifer"        
    ##  [21] "Abantennarius rosaceus"         "Doboatherina duodecimalis"     
    ##  [23] "Stenatherina panatela"          "Hypoatherina temminckii"       
    ##  [25] "Platybelone argalus"            "Strongylura leiura"            
    ##  [27] "Strongylura urvillii"           "Tylosurus crocodilus"          
    ##  [29] "Cheilopogon spilonotopterus"    "Exocoetus monocirrhus"         
    ##  [31] "Photoblepharon palpebratum"     "Sargocentron diadema"          
    ##  [33] "Sargocentron ittodai"           "Sargocentron microstoma"       
    ##  [35] "Sargocentron punctatissimum"    "Corythoichthys intestinalis"   
    ##  [37] "Doryrhamphus janssi"            "Doryrhamphus excisus"          
    ##  [39] "Hippichthys cyanospilos"        "Microphis leiaspis"            
    ##  [41] "Lophocampus brevidorsalis"      "Lophocampus retzii"            
    ##  [43] "Microphis manadensis"           "Fistularia commersonii"        
    ##  [45] "Scorpaenopsis oxycephalus"      "Leiopotherapon plumbeus"       
    ##  [47] "Kyphosus ocyurus"               "Cirrhitus pinnulatus"          
    ##  [49] "Nemanthias bartlettorum"        "Nemanthias bicolor"            
    ##  [51] "Nemanthias dispar"              "Pyronotanthias flavoguttatus"  
    ##  [53] "Pyronotanthias lori"            "Pyronotanthias parvirostris"   
    ##  [55] "Mirolabrichthys pascalus"       "Pyronotanthias smithvanizi"    
    ##  [57] "Mirolabrichthys tuka"           "Cephalopholis boenak"          
    ##  [59] "Aethaloperca rogaa"             "Gracila albomarginata"         
    ##  [61] "Anyperodon leucogrammicus"      "Chromileptes altivelis"        
    ##  [63] "Epinephelus morrhua"            "Pictichromis paccagnellorum"   
    ##  [65] "Heteropriacanthus cruentatus"   "Cercamia eremia"               
    ##  [67] "Cheilodipterus isostigma"       "Cheilodipterus arabicus"       
    ##  [69] "Ostorhinchus neotes"            "Zoramia leptacanthus"          
    ##  [71] "Ferdauia orthogrammus"          "Trachinotus baillonii"         
    ##  [73] "Leiognathus equula"             "Paracaesio kusakarii"          
    ##  [75] "Symphorus nematophorus"         "Lutjanus ehrenbergii"          
    ##  [77] "Lutjanus russellii"             "Diagramma pictum"              
    ##  [79] "Gnathodentex aureolineatus"     "Nemipterus hexodon"            
    ##  [81] "Chaetodon auriga"               "Heniochus singularius"         
    ##  [83] "Centropyge loriculus"           "Centropyge fisheri"            
    ##  [85] "Paracentropyge multifasciata"   "Centropyge venusta"            
    ##  [87] "Pygoplites diacanthus"          "Chaetodontoplus poliourus"     
    ##  [89] "Stegastes lacrymatus"           "Plectroglyphidodon fasciolatus"
    ##  [91] "Stegastes punctatus"            "Pycnochromis acares"           
    ##  [93] "Pycnochromis agilis"            "Pycnochromis amboinensis"      
    ##  [95] "Pycnochromis atripes"           "Azurina brevirostris"          
    ##  [97] "Pycnochromis caudalis"          "Pycnochromis delta"            
    ##  [99] "Azurina elerae"                 "Azurina lepidolepis"           
    ## [101] "Pycnochromis lineatus"          "Pycnochromis margaritifer"     
    ## [103] "Pycnochromis retrofasciatus"    "Pycnochromis vanderbilti"      
    ## [105] "Amphiprion biaculeatus"         "Amblyglyphidodon batunaorum"   
    ## [107] "Chrysiptera brownriggii"        "Bodianus loxozonus"            
    ## [109] "Bodianus paraleucosticticus"    "Xiphocheilus typus"            
    ## [111] "Decodon pacificus"              "Cirrhilabrus roseafascia"      
    ## [113] "Novaculichthys taeniourus"      "Anampses neoguinaicus"         
    ## [115] "Halichoeres marginatus"         "Thalassoma jansenii"           
    ## [117] "Scarus ghobban"                 "Parapercis millepunctata"      
    ## [119] "Parapercis schauinslandii"      "Aspidontus taeniatus"          
    ## [121] "Meiacanthus atrodorsalis"       "Petroscirtes thepassii"        
    ## [123] "Plagiotremus laudandus"         "Plagiotremus rhinorhynchos"    
    ## [125] "Plagiotremus tapeinosoma"       "Crossosalarias macrospilus"    
    ## [127] "Nannosalarias nativitatis"      "Salarias ramosus"              
    ## [129] "Callionymus delicatulus"        "Callionymus simplicicornis"    
    ## [131] "Callionymus enneactis"          "Rhyacichthys aspro"            
    ## [133] "Bostrychus sinensis"            "Eleotris acanthopomus"         
    ## [135] "Giuris margaritaceus"           "Oxyeleotris lineolata"         
    ## [137] "Smilosicyopus fehlmanni"        "Sicyopus zosterophorus"        
    ## [139] "Lophogobius bleekeri"           "Mugilogobius platystoma"       
    ## [141] "Oxyurichthys microlepis"        "Oxyurichthys ophthalmonema"    
    ## [143] "Palutrus pruinosa"              "Pseudogobius poicilosoma"      
    ## [145] "Cryptocentrus maudae"           "Stonogobiops xanthorhinicus"   
    ## [147] "Gobius bontii"                  "Ctenogobius chusanensis"       
    ## [149] "Callogobius clitellus"          "Callogobius hasseltii"         
    ## [151] "Favonigobius melanobranchus"    "Glossogobius giuris"           
    ## [153] "Gobiopsis aporia"               "Istigobius goldmanni"          
    ## [155] "Istigobius nigroocellatus"      "Paragobiodon melanosoma"       
    ## [157] "Priolepis pallidicincta"        "Trimma gigantum"               
    ## [159] "Trimma hotsarihiense"           "Acentrogobius nebulosus"       
    ## [161] "Ptereleotris grammica"          "Platax boersii"                
    ## [163] "Siganus fuscescens"             "Siganus corallinus"            
    ## [165] "Acanthurus nigroris"            "Acanthurus triostegus"         
    ## [167] "Sphyraena putnamae"             "Ariomma brevimanus"            
    ## [169] "Pseudobalistes fuscus"          "Canthidermis maculata"         
    ## [171] "Balistoides viridescens"        "Rhinecanthus verrucosus"       
    ## [173] "Xanthichthys caeruleolineatus"  "Paramonacanthus curtorhynchos" 
    ## [175] "Ostracion cubicum"              "Rhynchostracion nasus"         
    ## [177] "Ostracion solorense"            "Sphoeroides pachygaster"       
    ## [179] "Canthigaster axiologus"

``` r
# Here you can see that most names that were corrected were due to minimal spelling mistakes
# Currently at 179

check2 <- setdiff(fish$X8, fish$Species)
check2
```

    ##   [1] "Stegastoma fasciatum"                       
    ##   [2] "Sphyrna mokorran"                           
    ##   [3] "Rhina anclystoma"                           
    ##   [4] "Neotrygon kuhli"                            
    ##   [5] "Anguilla pacifica"                          
    ##   [6] "Gymnothorax richardsoni"                    
    ##   [7] "Gymnothorax rueppellii"                     
    ##   [8] "Gymnothorax Siderea pictus"                 
    ##   [9] "Uropterygius alboguttatus"                  
    ##  [10] "Evipes percinctus"                          
    ##  [11] "Ichthyapus vulturus"                        
    ##  [12] "Lamnastoma orientalis"                      
    ##  [13] "Lamnastoma polyophthalmum"                  
    ##  [14] "Encrasicholina pseudoheteroloba"            
    ##  [15] "Thrissina baelama"                          
    ##  [16] "Thrissina hamiltoni"                        
    ##  [17] "Thrissina mystax"                           
    ##  [18] "Antennatus analis"                          
    ##  [19] "Antennatus coccineus"                       
    ##  [20] "Antennatus dorehensis"                      
    ##  [21] "Antennatus nummifer"                        
    ##  [22] "Antennatus rosaceus"                        
    ##  [23] "Atherinomorus duodecimalis"                 
    ##  [24] "Hypoatherina panatela"                      
    ##  [25] "Hypoatherina temmincki"                     
    ##  [26] "Platybelone argalus platyura"               
    ##  [27] "Strongylura leiura leiura"                  
    ##  [28] "Strongylura urvilli"                        
    ##  [29] "Tylosurus crocodilus crocodilus"            
    ##  [30] "Cheilopogon spilonotopteru"                 
    ##  [31] "Exocoetus monochirus"                       
    ##  [32] "Photoblepharon palpebratus"                 
    ##  [33] "Sargocentron Neoniphon diadema"             
    ##  [34] "Sargocentron Neoniphon ittodai"             
    ##  [35] "Sargocentron Neoniphon microstoma"          
    ##  [36] "Sargocentron Neoniphon punctatissimum"      
    ##  [37] "Corythoichthys waitei"                      
    ##  [38] "Doryrhamphus jansii"                        
    ##  [39] "Doryrhamphus melanopleura"                  
    ##  [40] "Hippichthys cyanospilus"                    
    ##  [41] "Microphis Oostethus leiaspis"               
    ##  [42] "Microphis brevidorsalis"                    
    ##  [43] "Microphis Lophocampus retzii"               
    ##  [44] "Microphis Oostethus manadensis"             
    ##  [45] "Fistularia commersoni"                      
    ##  [46] "Scorpaenopsis oxycephala"                   
    ##  [47] "Leiopotheropon plumbeus"                    
    ##  [48] "Kyphosus Sectator ocyurus"                  
    ##  [49] "Cirrhitus pinnulatus pinnulatus"            
    ##  [50] "Pseudanthias bartlettorum"                  
    ##  [51] "Pseudanthias bicolor"                       
    ##  [52] "Pseudanthias dispar"                        
    ##  [53] "Pseudanthias flavoguttatus"                 
    ##  [54] "Pseudanthias lori"                          
    ##  [55] "Pseudanthias parvirostris"                  
    ##  [56] "Pseudanthias pascalus"                      
    ##  [57] "Pseudanthias smithvanizi"                   
    ##  [58] "Pseudanthias tuka"                          
    ##  [59] "Cephalopholis boenack"                      
    ##  [60] "Cephalopholis Aethaloperca rogaa"           
    ##  [61] "Cephalopholis Gracila albomarginata"        
    ##  [62] "Epinephelus Anyperodon leucogrammicus"      
    ##  [63] "Epinephelus Cromileptes altivelis"          
    ##  [64] "Mycteroperca morrhua"                       
    ##  [65] "Pictichromis paccagnellae"                  
    ##  [66] "Heteropriacanthus carolinus"                
    ##  [67] "Cercamia eremeia"                           
    ##  [68] "Cheilodipterus isostigmus"                  
    ##  [69] "Cheilodipterus lineatus"                    
    ##  [70] "Ostorhinchus Brephamia neotes"              
    ##  [71] "Zoramia leptacantha"                        
    ##  [72] "Carangoides orthogrammus"                   
    ##  [73] "Trachinotus bailloni"                       
    ##  [74] "Leiognathus equulus"                        
    ##  [75] "Paracaesio kusakaraii"                      
    ##  [76] "Symphurus nematophorus"                     
    ##  [77] "Lutjanus ehrenbergi"                        
    ##  [78] "Lutjanus russelli"                          
    ##  [79] "Diagramma pictum pictum"                    
    ##  [80] "Gnathodentex aurolineatus"                  
    ##  [81] "Nemipterus hexadon"                         
    ##  [82] "Chaetodon auriga setifer"                   
    ##  [83] "Heniochus singularis"                       
    ##  [84] "Centropyge loricula"                        
    ##  [85] "Centropyge Xiphipos fisheri"                
    ##  [86] "Paracentropyge Centropyge multifasciata"    
    ##  [87] "Paracentropyge Centropyge venusta"          
    ##  [88] "Pygoplites diacanthus diacanthus"           
    ##  [89] "Chaetodontoplus polioururs"                 
    ##  [90] "Plectroglyphidodon lacrymatus"              
    ##  [91] "Stegastes fasciolatus"                      
    ##  [92] "Stegastes pr nsp"                           
    ##  [93] "Chromis acares"                             
    ##  [94] "Chromis agilis"                             
    ##  [95] "Chromis amboinensis"                        
    ##  [96] "Chromis atripes"                            
    ##  [97] "Chromis brevirostris"                       
    ##  [98] "Chromis caudalis"                           
    ##  [99] "Chromis delta"                              
    ## [100] "Chromis elerae"                             
    ## [101] "Chromis lepidolepis"                        
    ## [102] "Chromis lineata"                            
    ## [103] "Chromis margaritifer"                       
    ## [104] "Chromis retrofasciata"                      
    ## [105] "Chromis vanderbilti"                        
    ## [106] "Premnas biaculeatus"                        
    ## [107] "Amblyglyphidodon batunai"                   
    ## [108] "Chrysiptera leucopoma"                      
    ## [109] "Bodianus loxozonus loxozonus"               
    ## [110] "Bodianus paraleucostictus"                  
    ## [111] "Choerodon typus"                            
    ## [112] "Decadon pacificus"                          
    ## [113] "Oxycheilinus oxyrhinchus"                   
    ## [114] "Cirrhilabrus roseofascia"                   
    ## [115] "Novaculicthys taeniourus"                   
    ## [116] "Anampses neoguinaecus"                      
    ## [117] "Halichoeres annularis"                      
    ## [118] "Thalassoma janseni"                         
    ## [119] "Scarus pyrrostethus"                        
    ## [120] "Parapercis millipunctata"                   
    ## [121] "Parapercis schauinslandi"                   
    ## [122] "Aspidontus taeniatus taeniatus"             
    ## [123] "Meiacanthus atrodorsalis atrodorsalis"      
    ## [124] "Petroscirtes thepasi"                       
    ## [125] "Plagiotremus laudandus laudandus"           
    ## [126] "Plagiotremus rhinorhynchus"                 
    ## [127] "Plagiotremus tapienosoma"                   
    ## [128] "Cirripectes macrospilos"                    
    ## [129] "Nannosalarias nativitatus"                  
    ## [130] "Salarius ramosus"                           
    ## [131] "Callionymus (Calliurichthys) delicatulus"   
    ## [132] "Callionymus (Calliurichthys) simplicicornis"
    ## [133] "Callionymus (Paradiplogrammus) enneactis"   
    ## [134] "Rhyacichthys asporo"                        
    ## [135] "Bostrichthys sinensis"                      
    ## [136] "Eleotris acanthopoma"                       
    ## [137] "Giurus margaritaceus"                       
    ## [138] "Hypseleotris bipartita"                     
    ## [139] "Oxyeleotris lineolatus"                     
    ## [140] "Oxyeleotris pelewensis"                     
    ## [141] "Oxyeleotris pinguis"                        
    ## [142] "Sicyopus fehlmanni"                         
    ## [143] "Sicyopus zosterophorum"                     
    ## [144] "Lophiogobius bleekeri"                      
    ## [145] "Mugiligobius platystomus"                   
    ## [146] "Oxyurichthys longicauda"                    
    ## [147] "Oxyurichthys ophthalmonemus"                
    ## [148] "Palutris pruinosa"                          
    ## [149] "Pseudogobius javanicus"                     
    ## [150] "Cryptocentrus maudei"                       
    ## [151] "Stonogobiops xanthorhinica"                 
    ## [152] "Acentrogobius bonti"                        
    ## [153] "Acentrogobius chusanensis"                  
    ## [154] "Callogobius clittelus"                      
    ## [155] "Callogobius hasselti"                       
    ## [156] "Favonigobius Papillogobius melanobranchus"  
    ## [157] "Glossogobius guirus"                        
    ## [158] "Gobiopsis liolepis"                         
    ## [159] "Istigobius goldmanini"                      
    ## [160] "Istigobius nigrioocellatus"                 
    ## [161] "Paragobiodon melanosomus"                   
    ## [162] "Priolepis pallidocincta"                    
    ## [163] "Trimma giganteum"                           
    ## [164] "Trimma hotsarihiensis"                      
    ## [165] "Yongeichthys nebulosus"                     
    ## [166] "Ptereleotris grammica grammica"             
    ## [167] "Platax boersi"                              
    ## [168] "Siganus margaritiferus"                     
    ## [169] "Siganus studeri"                            
    ## [170] "Acanthurus nigros"                          
    ## [171] "Acanthurus triostegus triostegus"           
    ## [172] "Sphyraena putnamiae"                        
    ## [173] "Arioma brevimanum"                          
    ## [174] "Balistes fuscus"                            
    ## [175] "Canthidermis maculatus"                     
    ## [176] "Pseudobalistes viridescens"                 
    ## [177] "Rhinecanthus verrucosa"                     
    ## [178] "Xanthichthys careuleolineatus"              
    ## [179] "Paramonacanthus curtorhynchus"              
    ## [180] "Ostracion cubicus"                          
    ## [181] "Ostracion nasus"                            
    ## [182] "Ostracion solorensis"                       
    ## [183] "Sphaeroides pachygaster"                    
    ## [184] "Canthigaster axilogus"

``` r
# check2 includes species names that are obsolete and were changed to something else that was already in the list
compare_ref_species_names <- as.data.frame(cbind(check, check2))
```

    ## Warning in cbind(check, check2): number of rows of result is not a multiple of
    ## vector length (arg 1)

``` r
# Will be a different length when compared to "check" if there are duplicates for species names that are now something else that previously existed
# The warning is because the vectors are different lengths
write.csv(compare_ref_species_names, "/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/species_lists/compare_ref_species_names.csv")

# Change column name, it will duplicate the column with a new column name
fish$REF <- fish$X29

# Extract only the genus and species names now that we know that this file only contains species from Palau
fish <- fish[,c(6:7)]
```

## Surveys

- Use the community surveys and extract the species and incidence data

``` r
# A file with the presence and absence of species observed in the marine lakes. Names have been edited to reflect spelling in FishBase.
palau_records <- read.xlsx("/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/species_surveys/palau_marine_lake_fish_biodiversity_ordered&corrected.xlsx", sheet = 1, startRow = 3, colNames = T)

# Extract the lake columns with species names. The lake columns have 3 letter codes and if a species is present there will be a 1 in the column.
site_fish <- palau_records[-c(254:259), c(2:3,6,11,16,21,26,31,37,43:44,53,58,63,68,73:74,79,84,89,94,101,106,111)]

# Combine Genus and Species
site_fish$Species <- paste(site_fish[, 1], site_fish[, 2])

site_fish[is.na(site_fish)] <- 0; site_fish[site_fish == "x"] <- 0; site_fish[site_fish == "?"] <- 0

# Convert all fields (apart from species_name) to numeric
# Avoid the species column when using this command
site_fish[, -c(1,2,25)] <- lapply(site_fish[, -c(1,2,25)], function(x) as.numeric(x))

# Identify what species are found in the lakes but not the reference pool
site_data_only <- setdiff(site_fish$Species, fish$Species)
site_data_only
```

    ## [1] "Acanthurus sp"    "Atherinomorus sp" "Epinephelus sp"   "Myripristis sp"  
    ## [5] "Parioglossus sp"  "Platax sp"        "Pomacentrus sp 2" "Pomacentrus sp"  
    ## [9] "Silhouettea sp"

``` r
# [1] "Acanthurus sp" "Atherinomorus sp" "Epinephelus sp" "Myripristis sp" "Parioglossus sp"
# [6] "Platax sp" "Pomacentrus sp 2" "Pomacentrus sp" "Silhouettea sp"
```

## Rid unwanted species

``` r
# Bring together marine lake presence/absence surveys and palau survey
pres_abs <- merge(site_fish, fish, by = "Species", all = TRUE)

# Create the presence/absence dataframe
pres_abs <- pres_abs[, c(1, 4:26)]

# Convert all NAs (and a few remaining "x" and "?") to 0
pres_abs[is.na(pres_abs)] <- 0; pres_abs[pres_abs == "x"] <- 0; pres_abs[pres_abs == "?"] <- 0

# Convert all fields (apart from species_name) to numeric
# Avoid the species column when using this command
pres_abs[, -1] <- lapply(pres_abs[, -1], function(x) as.numeric(x))

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
pres_abs <- pres_abs[which(rowSums(pres_abs[, -1]) > 0),]

# Removes the species that are not present in any of the environments
pres_abs <- aggregate(pres_abs[,-1], list(names = pres_abs$Species), max, na.rm = TRUE)

pres_abs$REF[pres_abs$REF > "-1"] <- 1

# Extract only the species name
species_names <- pres_abs[, 1]
# 1724 fish species found in Palau and the lakes surveyed but 9 species from the site surveys were only identified to the genus
```

## FishBase validation

``` r
# Run to check for duplicates, you will delete these duplicates in the next line
duplicates_check <- which(duplicated(species_names))
duplicates_check 
```

    ## integer(0)

``` r
## Numbers may change from the above code identifying duplicates, so update each time
# Remove duplicates if need be
#species_names <- unique(species_names)

## Validate the species names against FishBase database
# This file is used to revise the original files to prevent having to update species names every time
val_species_names <- validate_names(species_names)
```

    ## Joining with `by = join_by(Subfamily, GenCode, FamCode)`
    ## Joining with `by = join_by(FamCode)`
    ## Joining with `by = join_by(Order, Ordnum, Class, ClassNum)`
    ## Joining with `by = join_by(Class, ClassNum)`

``` r
## Get rid of any species that are not valid by FishBase
# Be careful because a species may be misspelled or not updated or added by FishBase leading to an NA result
val_species_names <- na.omit(val_species_names)
# 1674 of the 1724 species names are validated by FishBase

# For some reason these species are not validate in FishBase even though they show up there.
Ci <- "Cheilodipterus isostigma"
Cr <- "Cirrhilabrus ryukyuensis"
Tm <- "Trimma macrophthalmus"
Zl <- "Zoramia leptacanthus"

# Add this species to the validated group because you can still download its traits.
val_species_names <- append(val_species_names, values = c(Ci, Cr, Tm, Zl))
# Now we have 1678 validated

# For some reason FishBase adds Genus Species to some of the species that are not in its database when you validate species names
blanks <- setdiff(val_species_names, species_names)
blanks
```

    ## [1] "Genus Species"

``` r
word_to_delete <- "Genus Species"
val_species_names <- val_species_names[val_species_names != word_to_delete]
# This shows us that only 1675 species were validated

# Create a file of species that were only identified to genus level. They will be lost by the validate_names command above and will need to be added back in
incomplete_to_add <- setdiff(species_names, val_species_names)
incomplete_to_add
```

    ##  [1] "Acanthurus sp"              "Apogon sp"                 
    ##  [3] "Asterropteryx sp"           "Atherinomorus sp"          
    ##  [5] "Bathygobius sp"             "Bryaninops sp"             
    ##  [7] "Callogobius sp"             "Callogobius sp B"          
    ##  [9] "Cercamia sp"                "Chromis anadema"           
    ## [11] "Chromis sp"                 "Chrysiptera ellenae"       
    ## [13] "Cirrhilabrus sp"            "Cryptocentrus sp"          
    ## [15] "Drombus sp"                 "Dussumieria sp"            
    ## [17] "Epinephelus sp"             "Eviota diyhritisma"        
    ## [19] "Eviota sp"                  "Eviota sp 2"               
    ## [21] "Fusigobius sp"              "Glossogobius sandakanensis"
    ## [23] "Gobiodon macrochir"         "Gobiodon sp"               
    ## [25] "Gorgasia sp"                "Liopropoma sp"             
    ## [27] "Liopropoma sp B"            "Lubbockichthys sp"         
    ## [29] "Myripristis sp"             "Ophichthus sp"             
    ## [31] "Opistognathus wassi"        "Pandaka sp"                
    ## [33] "Parioglossus sp"            "Platax sp"                 
    ## [35] "Plectranthias sp"           "Pomacentrus sp"            
    ## [37] "Rhabdoblennius sp"          "Schindleria sp"            
    ## [39] "Sebastapistes sp"           "Silhouettea sp"            
    ## [41] "Soleichthys sp"             "Tetronarce sp"             
    ## [43] "Tomiyamichthys sp"          "Tosanoides sp"             
    ## [45] "Trimma sp"                  "Trimma sp 2"               
    ## [47] "Trimmatom sp"               "Xenisthmus sp"             
    ## [49] "Xenisthmus sp 2"

``` r
# 48 species were not validated by FishBase. Some have Genus and Species but have yet to be added to FishBase. Others with an sp cannot be validated
write.csv(incomplete_to_add, "/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/species_lists/incomplete_to_add.csv")

# Add back in species that were lost by validate_names command
existing_species <- c(val_species_names, incomplete_to_add)

# Checking for duplicates and comparing files
duplicates_check <- which(duplicated(existing_species))
duplicates_check
```

    ## integer(0)

``` r
# Check differences
diff_check <- setdiff(existing_species, species_names)
diff_check
```

    ## character(0)

``` r
# Resort your vector since we added some species/changed names. Will come in handy later mainly
existing_species <- sort(existing_species)

# Create these files for when revising the original file so that you have all the correct names and the NAs that now replace lost species
# Once you have loaded the packages and loaded the file containing all species from survey data, you can download trait data from fishbase using fish_traits.Rmd. #From trait-beta-div-processing-fish-traits.R

write.csv(species_names, "/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/species_lists/species_names.csv")

write.csv(val_species_names, "/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/species_lists/val_species_names.csv")

write.csv(existing_species, "/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/species_lists/existing_species.csv")
```

## Aesthetics

- Edit occurrence data from the locations that were sampled. Isolate
  just the occurrence data from everything else that is in the file.
  This is for future when we compare it with the environmental data.
  Take out species not found in any location and edit the file for
  formatting purposes. \#From trait-beta-div-processing.R

``` r
# 1. BCM (S), 2. CLM (S), 3. FLK (M), 4. GLK (S), 5. HLM (S), 6. HLO (M), 7. IBK (O), 8. LLN (M), 9. LCN (O), 10. MLN (M), 11. NCN (O), 12. NLK (S), 13. NLN (M), 14. NLU (M), 15. OLO (M), 16. OOO (O), 17. OTM (S), 18. OOM (O), 19. RCA (O), 20. REF (R), 21. SLN (S), 22. TLN (S), 23. ULN (M)

# Sort names
pres_abs$names <- sort(pres_abs$names)

# This should hopefully return a file without duplicates unless things are named differently
pres_abs <- pres_abs[match(existing_species, pres_abs$names),]

# Replace the space in the species name with an underscore
pres_abs$names <- gsub(" ", "_", pres_abs$names)

duplicates_check <- which(duplicated(pres_abs))
duplicates_check 
```

    ## integer(0)

``` r
# Rename rows
row.names(pres_abs) <- pres_abs$names

# Order species names by alphabetical order
pres_abs <- pres_abs[order(row.names(pres_abs)),]

# Reorder lake names
# Define the desired column order
desired_order <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "LCN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOO", "OTM", "OOM", "RCA", "REF", "SLN", "TLN", "ULN")

# Reorder the columns based on the desired order
pres_abs <- pres_abs[, desired_order]

# Transpose matrix to generate lake-by-species matrix
pres_abs_t <- as.data.frame(t(pres_abs))
```

## Community type incidence matrices

- Create a presence absence matrix for each community/stratification
  type

``` r
keep <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "LCN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOO", "OTM", "OOM", "RCA", "SLN", "TLN", "ULN")
surveyed_site_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
surveyed_site_fish <- surveyed_site_fish[which(rowSums(surveyed_site_fish) > 0),]


keep <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "LLN", "MLN", "NLK", "NLN", "NLU", "OLO", "OTM", "SLN", "TLN", "ULN")
marine_lake_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
marine_lake_fish <- marine_lake_fish[which(rowSums(marine_lake_fish) > 0),]


keep <- c("IBK", "LCN", "NCN", "OOO", "OOM", "RCA")
ocean_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
ocean_fish <- ocean_fish[which(rowSums(ocean_fish) > 0),]

keep <- c("FLK", "HLM", "HLO", "LLN", "MLN", "NLN", "NLU", "OLO", "ULN")
holomictic_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
holomictic_fish <- holomictic_fish[which(rowSums(holomictic_fish) > 0),]

keep <- c("BCM", "CLM", "GLK", "NLK", "OTM", "SLN", "TLN")
meromictic_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
meromictic_fish <- meromictic_fish[which(rowSums(meromictic_fish) > 0),]

keep <- c("FLK", "HLO", "LLN", "MLN", "NLN", "NLU", "OLO", "ULN")
mixed_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
mixed_fish <- mixed_fish[which(rowSums(mixed_fish) > 0),]

keep <- c("BCM", "CLM", "GLK", "HLM", "NLK", "OTM", "SLN", "TLN")
stratified_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
stratified_fish <- stratified_fish[which(rowSums(stratified_fish) > 0),]


# Create extra column to be able to extract out lake specific incidences
pres_abs$zero <- 0

keep <- c("BCM", "zero")
BCM_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
BCM_fish <- BCM_fish[which(rowSums(BCM_fish) > 0),]
BCM_fish <- BCM_fish["BCM"]

keep <- c("CLM", "zero")
CLM_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
CLM_fish <- CLM_fish[which(rowSums(CLM_fish) > 0),]
CLM_fish <- CLM_fish["CLM"]

keep <- c("FLK", "zero")
FLK_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
FLK_fish <- FLK_fish[which(rowSums(FLK_fish) > 0),]
FLK_fish <- FLK_fish["FLK"]

keep <- c("GLK", "zero")
GLK_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
GLK_fish <- GLK_fish[which(rowSums(GLK_fish) > 0),]
GLK_fish <- GLK_fish["GLK"]

keep <- c("HLM", "zero")
HLM_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
HLM_fish <- HLM_fish[which(rowSums(HLM_fish) > 0),]
HLM_fish <- HLM_fish["HLM"]

keep <- c("HLO", "zero")
HLO_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
HLO_fish <- HLO_fish[which(rowSums(HLO_fish) > 0),]
HLO_fish <- HLO_fish["HLO"]

keep <- c("IBK", "zero")
IBK_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
IBK_fish <- IBK_fish[which(rowSums(IBK_fish) > 0),]
IBK_fish <- IBK_fish["IBK"]

keep <- c("LLN", "zero")
LLN_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
LLN_fish <- LLN_fish[which(rowSums(LLN_fish) > 0),]
LLN_fish <- LLN_fish["LLN"]

keep <- c("LCN", "zero")
LCN_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
LCN_fish <- LCN_fish[which(rowSums(LCN_fish) > 0),]
LCN_fish <- LCN_fish["LCN"]

keep <- c("MLN", "zero")
MLN_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
MLN_fish <- MLN_fish[which(rowSums(MLN_fish) > 0),]
MLN_fish <- MLN_fish["MLN"]

keep <- c("NCN", "zero")
NCN_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
NCN_fish <- NCN_fish[which(rowSums(NCN_fish) > 0),]
NCN_fish <- NCN_fish["NCN"]

keep <- c("NLK", "zero")
NLK_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
NLK_fish <- NLK_fish[which(rowSums(NLK_fish) > 0),]
NLK_fish <- NLK_fish["NLK"]

keep <- c("NLN", "zero")
NLN_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
NLN_fish <- NLN_fish[which(rowSums(NLN_fish) > 0),]
NLN_fish <- NLN_fish["NLN"]

keep <- c("NLU", "zero")
NLU_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
NLU_fish <- NLU_fish[which(rowSums(NLU_fish) > 0),]
NLU_fish <- NLU_fish["NLU"]

keep <- c("OLO", "zero")
OLO_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
OLO_fish <- OLO_fish[which(rowSums(OLO_fish) > 0),]
OLO_fish <- OLO_fish["OLO"]

keep <- c("OOO", "zero")
OOO_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
OOO_fish <- OOO_fish[which(rowSums(OOO_fish) > 0),]
OOO_fish <- OOO_fish["OOO"]

keep <- c("OTM", "zero")
OTM_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
OTM_fish <- OTM_fish[which(rowSums(OTM_fish) > 0),]
OTM_fish <- OTM_fish["OTM"]

keep <- c("OOM", "zero")
OOM_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
OOM_fish <- OOM_fish[which(rowSums(OOM_fish) > 0),]
OOM_fish <- OOM_fish["OOM"]

keep <- c("RCA", "zero")
RCA_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
RCA_fish <- RCA_fish[which(rowSums(RCA_fish) > 0),]
RCA_fish <- RCA_fish["RCA"]

keep <- c("SLN", "zero")
SLN_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
SLN_fish <- SLN_fish[which(rowSums(SLN_fish) > 0),]
SLN_fish <- SLN_fish["SLN"]

keep <- c("TLN", "zero")
TLN_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
TLN_fish <- TLN_fish[which(rowSums(TLN_fish) > 0),]
TLN_fish <- TLN_fish["TLN"]

keep <- c("ULN", "zero")
ULN_fish <- pres_abs[,keep]

## Remove all species not found in any locations and duplicates
# Identifies which rows are greater than 0
# Avoid the species column when using this command
ULN_fish <- ULN_fish[which(rowSums(ULN_fish) > 0),]
ULN_fish <- ULN_fish["ULN"]

keep <- c("BCM", "CLM", "FLK", "GLK", "HLM", "HLO", "IBK", "LLN", "LCN", "MLN", "NCN", "NLK", "NLN", "NLU", "OLO", "OOO", "OTM", "OOM", "RCA", "REF", "SLN", "TLN", "ULN")
pres_abs <- pres_abs[,keep]
```

## Load out files

``` r
# Fish presence by lake
write.csv(pres_abs_t, "/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/fish_presence_matrix_by_lake.csv")

# Fish presence by species
write.csv(pres_abs,  "/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/fish_presence_matrix_by_species.csv")

# Marine lake fish presence by species 
write.csv(surveyed_site_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/surveyed_site_presence_by_species.csv")

# Marine lake fish presence by species 
write.csv(marine_lake_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/marine_lake_presence_by_species.csv")

# Ocean site fish presence by species 
write.csv(ocean_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/ocean_site_presence_by_species.csv")

# Holomictic fish presence by species 
write.csv(holomictic_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/holomictic_presence_by_species.csv")

# Meromictic fish presence by species 
write.csv(meromictic_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/meromictic_presence_by_species.csv")

# Mixed fish presence by species 
write.csv(mixed_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/mixed_presence_by_species.csv")

# Stratified fish presence by species 
write.csv(stratified_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/stratified_presence_by_species.csv")

# Surveyed sites fish presence by species each in separate files
write.csv(BCM_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/BCM_fish.csv")
write.csv(CLM_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/CLM_fish.csv")
write.csv(FLK_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/FLK_fish.csv")
write.csv(GLK_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/GLK_fish.csv")
write.csv(HLM_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/HLM_fish.csv")
write.csv(HLO_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/HLO_fish.csv")
write.csv(IBK_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/IBK_fish.csv")
write.csv(LLN_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/LLN_fish.csv")
write.csv(LCN_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/LCN_fish.csv")
write.csv(MLN_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/MLN_fish.csv")
write.csv(NCN_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/NCN_fish.csv")
write.csv(NLK_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/NLK_fish.csv")
write.csv(NLN_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/NLN_fish.csv")
write.csv(NLU_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/NLU_fish.csv")
write.csv(OLO_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/OLO_fish.csv")
write.csv(OOO_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/OOO_fish.csv")
write.csv(OTM_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/OTM_fish.csv")
write.csv(OOM_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/OOM_fish.csv")
write.csv(RCA_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/RCA_fish.csv")
write.csv(SLN_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/SLN_fish.csv")
write.csv(TLN_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/TLN_fish.csv")
write.csv(ULN_fish,"/Users/bailey/Documents/research/fish_biodiversity/data/collection/fish/pres_abs/ULN_fish.csv")

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
    ## [1] taxize_0.9.100   phytools_2.0-3   maps_3.4.1.1     ape_5.7-1       
    ## [5] dplyr_1.1.4      tidyr_1.3.0      rfishbase_4.1.2  openxlsx_4.2.5.2
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.0        optimParallel_1.0-2     fastmap_1.1.1          
    ##  [4] combinat_0.0-8          duckdb_0.9.2-1          digest_0.6.33          
    ##  [7] lifecycle_1.0.4         magrittr_2.0.3          compiler_4.3.1         
    ## [10] rlang_1.1.2             progress_1.2.3          tools_4.3.1            
    ## [13] igraph_1.5.1            utf8_1.2.4              yaml_2.3.7             
    ## [16] data.table_1.14.10      knitr_1.45              phangorn_2.11.1        
    ## [19] clusterGeneration_1.3.8 conditionz_0.1.0        askpass_1.2.0          
    ## [22] prettyunits_1.2.0       mnormt_2.1.1            scatterplot3d_0.3-44   
    ## [25] contentid_0.0.18        curl_5.2.0              xml2_1.3.6             
    ## [28] httpcode_0.3.0          expm_0.999-8            withr_2.5.2            
    ## [31] purrr_1.0.2             numDeriv_2016.8-1.1     grid_4.3.1             
    ## [34] fansi_1.0.6             iterators_1.0.14        MASS_7.3-60            
    ## [37] crul_1.4.0              cli_3.6.1               rmarkdown_2.25         
    ## [40] crayon_1.5.2            generics_0.1.3          rstudioapi_0.15.0      
    ## [43] httr_1.4.7              tzdb_0.4.0              DBI_1.1.3              
    ## [46] cachem_1.0.8            stringr_1.5.1           parallel_4.3.1         
    ## [49] vctrs_0.6.5             Matrix_1.6-4            jsonlite_1.8.8         
    ## [52] hms_1.1.3               foreach_1.5.2           glue_1.6.2             
    ## [55] codetools_0.2-19        stringi_1.8.2           quadprog_1.5-8         
    ## [58] tibble_3.2.1            pillar_1.9.0            htmltools_0.5.7        
    ## [61] openssl_2.1.1           R6_2.5.1                dbplyr_2.4.0           
    ## [64] doParallel_1.0.17       bold_1.3.0              evaluate_0.23          
    ## [67] lattice_0.22-5          readr_2.1.4             memoise_2.0.1          
    ## [70] Rcpp_1.0.11             zip_2.3.0               uuid_1.1-1             
    ## [73] fastmatch_1.1-4         coda_0.19-4             nlme_3.1-164           
    ## [76] xfun_0.41               fs_1.6.3                zoo_1.8-12             
    ## [79] pkgconfig_2.0.3
