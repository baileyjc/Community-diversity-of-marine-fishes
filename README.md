# Community diversity of marine fishes in complete and habitat islands

A repository dedicated to sharing how I processed and analyzed ecological data in R from fish surveys in Palau. The results of this project have now been published in the _Frontiers of Biogeography_.

Description of the data and file structure

The files and scripts from data collection and processing through analysis are included in this Dryad submission. The necessary files to run all the scripts are provided in CSV format, except for the phylogenetic trees, which are in "tre" format. The first script to be run is fish_presabs.Rmd. The files you will use in this R Markdown file are palau_fish.csv and palau_surveyed_fish.csv. The output files you generate will serve as the input files for fish_tree.Rmd and fish_traits.Rmd. You will also need the labeled_dated_tree.tre file to run fish_tree.Rmd. From there, you can run palau_env.Rmd to generate the environmental and geographic data summaries. You will use the annual-environmental-without_h2s_layer.csv and lake_physical_attributes.csv files as input to create the file env_by_lake.csv. The environmental and geographical data differ from those in a previously published manuscript (Rapacciuolo et al. 2019) because we calculated environmental data that included only measures from stratified lakes above the chemocline. Fish cannot survive below the chemocline. The surface area also had to be recalculated because several locations included in our study did not have measurements and we did not want user differences to influence our measurements. After running those collection scripts, you now have the necessary input files to run load_collection_data.Rmd. Key input files have been attached here, such as fish_presence_matrix_by_lake.csv (a presence/absence matrix of Palau fish species), palau_fish_tree.tre (final Palau species tree), and final_traits_gnathostomata.csv (the final trait matrix after consulting Patzner 2008). Once you have run load_collection_data.Rmd, you will produce key output files, surveyed_fish_tree.tre, surveyed_sites_traits.csv, and env.csv; run the remaining R Markdown analyses files in no particular order, SD_analyses.Rmd, PD_analyses.Rmd, FD_analyses.Rmd, and trait_plots.Rmd. With these scripts and the data files, you can run all the scripts and generate the same figures found in the paper and the results to create the tables. Keep in mind that the paths to input and output files need to change to your local environment. Also, FishBase does occasionally update the data for fishes, which means results may vary slightly.

The supplementary file, Fish_Biodiversity_Supporting_Information_final.pdf, from the manuscript is also provided because it includes the details regarding the traits, environmental, and geographical data collected. Additionally, the names and locations of marine lakes and ocean locations that are identified by acronyms in the files and scripts listed above are provided in the supplementary file.

The TLN_oxygen_concentration.Rmd file is used to show that the oxygen concentration at the bottom of T Lake (TLN) varies from oxic to anoxic. The file TLN_oxygen_concentrations.xlsx is used in that R Markdown file.

Rapacciuolo Giovanni, Beman J. Michael, Schiebelhut Lauren M. and Dawson Michael N. (2019). Microbes and macro-invertebrates show parallel β-diversity but contrasting α-diversity patterns in a marine natural experimentProc. R. Soc. B.28620190999 http://doi.org/10.1098/rspb.2019.0999

Code/software

All versions packages used by the scripts can be  are identified at the bottom of each .md file using the sessionInfo() command in R. This includes fish_presabs.md, fish_tree.md, fish_traits.md, palau_env.md, load_collection_data.md, SD_analyses.md, PD_analyses.md, FD_analyses.md, trait_plots.md, and TLN_oxygen_concentration.md.

The script p.adjust.envfit.R is a function ‘p.adjust.envfit’, which can be found here: https://www.davidzeleny.net/anadat-r/doku.php/en:start

Access information

Other publicly accessible locations of the data:

https://github.com/baileyjc/Community-diversity-of-marine-fishes

Data were derived from the following sources:

Fish traits were downloaded from FishBase (Froese & Pauly 2022) using the R package rfishbase *v4.1.2 *(Boettiger et al. 2012) and fish family-level reproduction from Patzner (2008).

Boettiger C, Lang DT, Wainwright PC (2012) rfishbase: exploring, manipulating and visualizing FishBase data from R. Journal of Fish Biology 81: 2030–2039. https://doi.org/10.1111/j.1095-8649.2012.03464.x

Froese R, Pauly D [Eds] (2024) FishBase. World Wide Web Electronic Publication. www.fishbase.org [version (10/2024)]

Patzner RA (2008) Reproductive strategies of fish. In: Rocha MJ, Arukwe A, Kapoor BG (Eds) Fish reproduction (311–350). Taylor & Francis Group.

The fish phylogeny was downloaded from the Open Tree of Life version 13.4 (Redelings et al. 2019) using the Chronosynth tool (McTavish & Sanchez-Reyes 2022).

McTavish EJ, Sanchez-Reyes LL (2022) Chronosnyth. https://github.com/OpenTreeOfLife/chronosynth

P. Colin, L. Bell, and staff at Coral Reef Research Foundation (CRRF) provided support in Palau and environmental data. T. Stieglitz & H. Stibor provided lake bathymetry and tidal propagation information. This work was conducted under permits from the Ministry of Natural Resources, Environment, and Tourism (RE-13-11) and Koror State Government (#13-233), Palau.
