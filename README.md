# Amphibian Skin 16S Nanopore Microbiome Against Woodhams Data Base of Bd-Inhibitory-taxa
Report for reproducibility purposes for analysis in paper (CITE WHEN PUBLISHED), where Woodhams data base was used to separate "putative inhibitory to Bd taxa" and compare community structure in samples od frogs infected and not infected with Batrachochytrium dendrobatidis (Bd), as well as the surrounding environmental water.

Analysis was performed using 16S rRNA amplicons obtained by nanopore sequencing. Nanopore reads processing and taxonomic assignation was done previously using Spaghetti pipeline: https://github.com/adlape95/Spaghetti . This R script should be used after Spaghetti pipeline is completed. 

### Repository contents

a) Woodhams Database : Woodhams, Douglas & Alford, Ross & Antwis, Rachael & Archer, Holly & Becker, Matthew & Belden, Lisa & Bell, Sara & Bletz, Molly & Daskin, Joshua & Davis, Leyla & Flechas, Sandra & Lauer, Antje & González, Antonio & Harris, Reid & Holden, Whitney & Hughey, Myra & Nez, Roberto & Knight, Rob & Kueneman, Jordan & McKenzie, Valerie. (2015). Antifungal isolates database of amphibian skin-associated bacteria and function against emerging fungal pathogens Ecological Archives E096-059. Ecology. 96. 2015-595. doi:10.1890/14-1837.1.

b) Personal R script: This script analyzes the microbial community with respect to Bd-inhibitory taxa, using two complementary measures:

   - Proportion of taxa (taxonomic proportion): The fraction of distinct taxa in each sample that belong to the Woodhams inhibitory category. Calculated as the number of inhibitory taxa divided by the total number of taxa detected in the sample. This metric reflects taxonomic representation independent of abundance.

   - Summed relative abundance of inhibitory taxa: The fraction of sequencing reads belonging to inhibitory taxa. First, relative abundance was calculated for each taxon (reads assigned to that taxon divided by total reads in the sample). Then, all inhibitory taxa were summed to obtain the total relative abundance of the inhibitory category. This metric reflects dominance of inhibitory taxa based on sequencing depth.
    
   - NOTE: These two metrics capture different ecological dimensions. A sample may have many rare inhibitory taxa (high proportion of taxa, low relative abundance) or a few dominant inhibitory taxa (low proportion of taxa, high relative abundance). Both are reported separately to avoid conflating “proportion” with “abundance.” This measurements were made using minimap2 to map nanopore sequences and finally, to generate a publication plot.


### Workflow summary

   1. Script is used to extract sequences (reads) from PAF files produced by the Spaghetti pipeline (taxonomic assignment results). And then, to create a merged table containing: sequence IDs, alignment IDs, relative abundances (per sample), taxonomic assignments, and other metadata from PAF files.

   2. Map extracted sequences against the Woodhams database using minimap2.

   3. Apply strict filtering criteria: minimum alignment length = 500 bp, minimum alignment identity to target sequences = 80%, mapping quality = 60 (highest mapQ). This strict criteria is suggested as we are only using 16S gene reads to infer a putative protective function, which is normally not suggested. But as Woodhams has made this data base and we are profiting from it, we advise against relaxing thresholds.

   4. Generate final results in CSV format. Files named "results_woodham_60" & "results_woodham_60_counts_proportions_fractions"

   5. The output table "results_woodham_60_counts_proportions_fractions" was used to produce the publication plots.
