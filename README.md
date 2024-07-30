# Amphibian Skin 16S Nanopore Microbiome Against Woodhams Data Base of Bd-Inhibitory-taxa

Report for reproducibility purposes.

### Woodhams DB assignation as "putative inhibitory to Bd taxa" and comparison with samples infected and not infected with Bd: 
### First time done using 16S rRNA amplicons obtained by nanopore sequencing
In this repository you will find

a) Woodhams Data Base : Woodhams, Douglas & Alford, Ross & Antwis, Rachael & Archer, Holly & Becker, Matthew & Belden, Lisa & Bell, Sara & Bletz, Molly & Daskin, Joshua & Davis, Leyla & Flechas, Sandra & Lauer, Antje & González, Antonio & Harris, Reid & Holden, Whitney & Hughey, Myra & Nez, Roberto & Knight, Rob & Kueneman, Jordan & McKenzie, Valerie. (2015). Antifungal isolates database of amphibian skin-associated bacteria and function against emerging fungal pathogens Ecological Archives E096-059. Ecology. 96. 2015-595. doi:10.1890/14-1837.1.

b) Personal R script to analyze the proportion of the community (based on taxa counts) and the fraction of the community (based on relative abundance) with putative inhibition to Batrachochytrium dendrobatidis (Bd) growth that is found in Amphibian skin microbiome samples and their surrounding environment, using minimap2 to map nanopore sequences and finally, to generate a publication plot.

   - Script is used to extract sequences (reads) from PAF files from Spaghetti pipeline (taxonomic assignation results) and to create a merged table with: seq ids, alignment ids, relative abundances (from each sample) and taxonomic assignations (and many other characteristics of PAF file)
   - Extracted sequences (reads) were then used to be mapped against Woodhams Data Base of putative inhibitory to Bd growth taxa, using minimap2

   - New mapping output PAFs were then processed with a strict threshold of:  min alingment length of 500 pb + a minimum alignment of 80% to target sequences and a mapQ of 60 (best quality alingments) was implemented 
   - Final results presented in csv files: "results_woodham_60" and "results_woodham_60_counts_proportions_fractions"
   - Output table "results_woodham_60_counts_proportions_fractions" was finally used to later generate publication plots


### NOTE: Nanopore reads processing and taxonomic assignation was done previously using Spaghetti pipeline: https://github.com/adlape95/Spaghetti 
This R script should be used after Spaghetti pipeline is completed.
