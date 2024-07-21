############################################################################
##############                                                ##############
##############                                                ##############
##############            MARIETTE VILADOMAT JASSO            ##############
##############                                                ##############
##############                                                ##############
##############     PUTATIVE INHIBITORY OF Bd GROWTH - TAXA    ##############
##############              Amphibian Microbiome              ##############
##############               USING: Woodhams DB               ##############
##############                                                ##############
##############                 16S Nanopore                   ##############
##############                                                ##############
##############                                                ##############
############################################################################

###  Nanopore reads processing and taxonomic assignation was done previously using Spaghetti pipeline: 
###  https://github.com/adlape95/Spaghetti 

###  This script was used to extract sequences (reads) from PAF files from Spaghetti pipeline (taxonomic assignation results) and
###     to create a merged table with: seq ids, alignment ids, relative abundances (from each sample) and taxonomic assignations (and many other characteristics of PAF file)
###  Extracted sequences (reads) were then used to be mapped against Woodhams Data Base of putative inhibitory to Bd growth taxa, using minimap2
######
###  New mapping output PAFs were then processed:
###  A threshold of:  min alingment length of 500 pb + a minimum alignment of 80% to target sequences and a mapQ of 60 (best quality alingments) was implemented 
###  Final results presented in csv files: "results_woodham_60" and "results_woodham_60_counts_proportions_fractions"
###  Output table "results_woodham_60_counts_proportions_fractions" was finally used to generate publication plot

### 0) LOAD LIBRARIES ----

library(dplyr)
library(Biostrings)
library(ShortRead)
library(msa)
library(purrr)
library(seqinr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggforce)
library(reshape2)

### 1) IMPORT TAX AND OTU FILES -----

tax <- read.csv("physeq_taxonomy_amfibios_COLAPSED_genus.csv")
colnames(tax)[1] <- "otu_tax_id"

otus <- read.csv("otu_table_amfibios_COLAPSED_genus_rel.csv")
colnames(otus)[1]<- "otu_tax_id"

### 2) IMPORT READ DATA (PAF FILES) ------------

folder_path <- "All_data/filteredPAFs"

# List all .paf files in the folder, excluding those with "mock" in the filename
paf_files <- list.files(folder_path, pattern = "\\.paf$", full.names = TRUE)
paf_files <- paf_files[!grepl("mock", paf_files)]

# Initialize an empty data frame to store the read_ids and seq_ids
paf_data <- data.frame(read_ids = character(), seq_ids = character(), stringsAsFactors = FALSE)

# Process each PAF file individually
for (paf_file in paf_files) {
  # Read the file, assuming tab-delimited format
  data <- read.table(paf_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Extract columns 1 and 6 and file name
  read_ids <- data[[1]]
  seq_ids <- data[[6]]
  file_name <- basename(paf_file)
  
  # Combine the columns into a temporary data frame
  temp_data <- data.frame(read_ids, seq_ids, file_name, stringsAsFactors = FALSE)
  
  # Append to the combined data frame
  paf_data <- rbind(paf_data, temp_data)
}

#remove temp data from memory
rm(temp_data, data, read_ids, seq_ids, file_name)

#free RAM
gc()

#format file_name to match fastq names
paf_data$file_name <- gsub("-f.paf", "", paf_data$file_name)
paf_data$file_name <- gsub(".scrubb", "", paf_data$file_name)
paf_data$read_ids <- gsub("_.*", "", paf_data$read_ids) #format read ID name

# unique(paf_data$seq_ids)
# unique(paf_data$read_ids)



### 3) EXTRACT READS WITH MATCH IN PAFs --------

folder_path <- "All_data"
fastq_filenames <- unique(paf_data$file_name)


# filter paf data with useful otu_tax_id from tax file
paf_data <- paf_data[paf_data$seq_ids %in% tax$otu_tax_id, ]


# Initialize an empty df for results
paf_with_sequences <- data.frame()

for (fastq_filename in fastq_filenames){
  
  print(fastq_filename)
  
  #subset fastq name from paf data
  paf_data_subset <- paf_data[paf_data$file_name == fastq_filename, ]
  
  # Read the FASTQ file
  fastq_file <- file.path(folder_path, fastq_filename)
  sequences <- readFastq(fastq_file)

  # Convert to DNAStringSet
  dna_sequences <- DNAStringSet(sread(sequences))
  names(dna_sequences) <- id(sequences)
  
  # Remove everything after the first space in sequence names
  names(dna_sequences) <- sub(" .*", "", names(dna_sequences))
  
  # delete fastq file to free RAM
  rm(sequences)
  
  # Subset sequences that got a hit with minimap2
  subset_sequences <- dna_sequences[names(dna_sequences) %in% paf_data_subset$read_ids]
  
  # Convert subset_sequences to a data frame
  subset_sequences <- data.frame(
    read_ids = names(subset_sequences),
    seq = as.character(subset_sequences),
    width = width(subset_sequences),
    stringsAsFactors = FALSE,
    row.names = NULL)
  
  #merge paf subset with corresponding sequences
  merged_subset <- merge(paf_data_subset, subset_sequences, by =c("read_ids"))
  
  # concat iterations into a single df output
  paf_with_sequences <- rbind(paf_with_sequences, merged_subset)
  
}

#free RAM
rm(merged_subset, subset_sequences, dna_sequences, paf_data_subset)
rm(paf_data)
gc()



### 4) MERGE SEQS +  OTUS +  TAX // CHECKPOINT-----

colnames(paf_with_sequences)[2] <- "otu_tax_id" #change name from seq_id to otu_tax_id for easier merging

otu_tax <- merge(tax, otus, by = c("otu_tax_id"))

#final merge
otu_tax_paf_with_sequences <- merge(otu_tax, paf_with_sequences, by = c("otu_tax_id"))

#save final merge
saveRDS(otu_tax_paf_with_sequences, "otu_tax_paf_with_sequences.RDS" )



### 5) OUTPUT READS WITH READ_ID AS HEADER // CONVERT TO FASTA FILES------

otu_tax_paf_with_sequences <- readRDS("otu_tax_paf_with_sequences.RDS" )

colnames(otu_tax_paf_with_sequences) # we're gonna need read_ids and seq columns

# Create a DNAStringSet object from the read_ids and seq columns
fasta_sequences_final <- DNAStringSet(otu_tax_paf_with_sequences$seq)
names(fasta_sequences_final) <- otu_tax_paf_with_sequences$read_ids

# Write the DNAStringSet object to a FASTA file
writeXStringSet(fasta_sequences_final, filepath = "otu_tax_paf_sequences.fasta")


### 6) MAPPING W/MINIMAP2 ----

# download minimap2 from github page:

   # curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2
   # tar -jxvf - ./minimap2-2.28_x64-linux/minimap2


# set up data base to use with minimap2:

   # minimap2 -d Amphibian-skin_bacteria_16S_sequences.mmi Amphibian-skin_bacteria_16S_sequences.fna


# otu_tax_paf_sequences.fasta was used as input for minimap2 against the Woodhams db using Spaghetti's defaults:

   # minimap2-2.28_x64-linux/minimap2 -x map-ont -t 4 --secondary=no -K 10M Amphibian-skin_bacteria_16S_sequences.mmi otu_tax_paf_sequences.fasta > otu_tax_paf_sequences.fasta.paf


# The result was otu_tax_paf_sequences.fasta.paf ; flag " --secondary=no " is for no secondary alignments in the output file. But sometimes, we still get secondary alignments
# So we used script " filterPAF.py " from the Spaghetti pipeline in order to remove secondary alignments AND alingments below 500 pb:

   # python filterPAF.py -i otu_tax_paf_sequences.fasta.paf > otu_tax_paf_sequences.fasta_FILTERED.paf


# The result was otu_tax_paf_sequences.fasta_FILTERED.paf, which needs to be processed and merged with the final table


### 7) PAF FILE FORMATTING  /// SELECTING THRESHOLD ----

# Read PAF file characteristics (https://cran.r-project.org/web/packages/pafr/vignettes/Introduction_to_pafr.html for reference)
paf <- read.table("otu_tax_paf_sequences.fasta_FILTERED.paf", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

paf <- paf[,1:12] #subset useful columns
colnames(paf) <- c("qname", "qlen", "qstart", "qend",
                   "strand", "tname", "tlen", "tstart", 
                   "tend", "nmatch", "alen", "mapq")


# Exploratory Data Analysis

# 1. check distributionsof relevant columns
hist(paf$mapq, breaks = 50); median(paf$mapq)
hist(paf$alen, breaks = 100); median(paf$alen)
hist(paf$qlen, breaks = 100); median(paf$qlen)
hist(paf$tlen, breaks = 100); median(paf$tlen); min(paf$tlen)
hist(paf$nmatch, breaks = 100); median(paf$nmatch)

# 2. make columns with %length query to target and target to query and check distributions
paf$perc_qlen_aligned <- paf$alen / paf$qlen 
hist(paf$perc_qlen_aligned, breaks = 100)
paf$perc_tlen_aligned <- paf$alen / paf$tlen 
hist(paf$perc_tlen_aligned, breaks = 100)


# Criteria chosen: mapq > 60 and >50 and paf$perc_tlen_aligned > 0.8 (used tlen because it's the limiting factor in terms of length (db has as low as 552 pb reference seqs))
paf_filtered_60 <- paf[paf$mapq >= 60 & paf$perc_tlen_aligned > 0.8,]
paf_filtered_50 <- paf[paf$mapq >= 50 & paf$perc_tlen_aligned > 0.8,] 



### 8) SUBSETTING otu_tax_paf_with_sequences AND MERGING WITH NEW paf_filtered_* ----

otu_tax_paf_with_sequences <- readRDS("otu_tax_paf_with_sequences.RDS" )

#remove columns no longer needed, for clarity
otu_tax_paf_with_sequences <- otu_tax_paf_with_sequences %>%
  select(c(-seq, -file_name, -width))

otu_tax_id_read_id <- otu_tax_paf_with_sequences %>%
  select(otu_tax_id, read_ids)

#remove reads_ids no longer needed, for clarity
otu_tax_paf_with_sequences <- otu_tax_paf_with_sequences %>%
  select(-read_ids)

#remove redundancies
otu_tax_paf_with_sequences <- distinct(otu_tax_paf_with_sequences)

#check if relative abundances from all samples add up to 1 (or 100).
sums <- otu_tax_paf_with_sequences %>%
  select_if(is.numeric) %>%
  summarise(across(everything(), ~ round(sum(.x, na.rm = TRUE), 0)))

if (sum(sums) == 100*length(sums)){
  print("relative abundances are ok")
}else{
  print("fuck")
}


### 9) MERGE otu_tax_id TO WOODHAMS HITS (mapq = 60) AND PLOT FRACTIONS AND PROPORTIONS IN MICROBIOMES-----

paf_filtered_selection_60 <- paf_filtered_60 %>%
  select(qname, tname)

# Merge the expanded data frame with the selected columns from paf_filtered_60
merged_df_60 <- otu_tax_id_read_id %>%
  inner_join(paf_filtered_selection_60, by = c("read_ids" = "qname"))


#check if there are otu_tax_ids with multiple woodhams assignations: YES (doesn't matter for the results, just curious)
woodhams_results_60  <- merged_df_60 %>% # 
  group_by(otu_tax_id)%>%
  summarize(unique_woodhams_hit_count = length(unique(tname)),
            unique_woodhams_hits = list(unique(tname)))

#access hits names: doesn't matter, just curious
names(woodhams_results_60$unique_woodhams_hits) <- woodhams_results_60$otu_tax_id
woodhams_results_60$unique_woodhams_hits

#chack hit counts: doesn't matter, just curious
hist(woodhams_results_60$unique_woodhams_hit_count, breaks = 100)

# add woodham results to otu_tax_id
otu_tax_paf_with_sequences$woodham_result_mapq60 <- ifelse(otu_tax_paf_with_sequences$otu_tax_id %in% woodhams_results_60$otu_tax_id, "putative_inhibitory_taxa", "other")

# OUTPUT TAX, REL ABUNDNACES, WOODHAM TAGS
write.csv(otu_tax_paf_with_sequences, "results_woodham_60.csv", row.names = F)



# get proportion of inhibitory taxa (number) and fraction of the sample that is inhibitory (relative abundance)
results_60 <- otu_tax_paf_with_sequences %>%
  group_by(woodham_result_mapq60) %>%
  summarise(
    across(
      starts_with(c("sample", "water")), 
      list(
        non_zero_count = ~ sum(. != 0),
        non_zero_sum = ~ sum(., na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  )

results_60_long <- results_60 %>%
  pivot_longer(
    cols = starts_with(c("sample_", "water")),
    names_to = c(".value", "sample"),
    names_sep = "_non_zero_"
  )

results_60_long <- as.data.frame(results_60_long)

rownames(results_60_long) <- paste0(results_60_long$woodham_result_mapq60, "_", results_60_long$sample)

results_60_long<- results_60_long %>%
  select(c(-woodham_result_mapq60, -sample))

results_60_long <- as.data.frame(t(results_60_long))

results_60_long$putative_inhibitory_taxa_proportion <- results_60_long$putative_inhibitory_taxa_count / (results_60_long$putative_inhibitory_taxa_count + results_60_long$other_count)
results_60_long$other_proportion <- results_60_long$other_count / (results_60_long$putative_inhibitory_taxa_count + results_60_long$other_count)

# OUTPUT COUNTS, PROPORTIONS, ABUNDNACES
write.csv(results_60_long, "results_woodham_60_counts_proportions_fractions.csv", row.names = T)



### INPUT FOR PLOTS 
results_60_long$sample <- rownames(results_60_long)

# a) FRACTION
results_long_melted_fraction <- results_60_long %>%
  pivot_longer(
    cols = c(other_sum, putative_inhibitory_taxa_sum),
    names_to = "category",
    values_to = "relative_abundance"
  ) %>%
  # Reorder sample based on putative_inhibitory_taxa_sum
  mutate(sample = factor(sample, levels = results_60_long %>%
                           arrange(desc(putative_inhibitory_taxa_sum)) %>%
                           pull(sample))) %>%
  select(c(sample, category, relative_abundance))

# b) PROPORTION OF TAXA
results_long_melted_proportion <- results_60_long %>%
  pivot_longer(
    cols = c(other_proportion, putative_inhibitory_taxa_proportion),
    names_to = "category",
    values_to = "proportion_of_taxa"
  ) %>%
  # Reorder sample based on putative_inhibitory_taxa_proportion
  mutate(sample = factor(sample, levels = results_60_long %>%
                           arrange(desc(putative_inhibitory_taxa_proportion)) %>%
                           pull(sample))) %>%
  select(c(sample, category, proportion_of_taxa))


# Replace "_sum" and "_proportion" with an empty string in the category column
results_long_melted_proportion <- results_long_melted_proportion %>%
  mutate(category = str_replace(category, "_sum|_proportion", ""))

results_long_melted_fraction <- results_long_melted_fraction %>%
  mutate(category = str_replace(category, "_sum|_proportion", ""))


# 3) PROPORTION + FRACTION MERGED
results_long_melted_ALL <- merge(results_long_melted_fraction, results_long_melted_proportion, by = c("sample", "category"))

#results_long_melted_ALL$relative_abundance<- results_long_melted_ALL$relative_abundance /100


### PLOTS!

ggplot(results_long_melted_ALL, aes(x = sample, y = relative_abundance, fill = category)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "",
    x = "Sample",
    y = "Fraction of the sample (relative_abundance)",
    fill = ""
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  

ggplot(results_long_melted_ALL, aes(x = sample, y = proportion_of_taxa, fill = category)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "",
    x = "Sample",
    y = "Proportion of taxa",
    fill = ""
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  



# try 1

p <- ggplot(results_long_melted_ALL, aes(x = sample, y = relative_abundance, fill = category)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  theme_minimal() +
  scale_fill_manual(values = c("other" = "#E1BE6A", "putative_inhibitory_taxa" = "#40B0A6")) +
  labs(
    title = "",
    x = "Sample",
    y = "Relative Abundance (%)",
    fill = ""
  ) +theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7.5),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.background = element_blank(), 
           axis.line = element_line(colour = "white"))


# Add pie charts as annotations
p <-p + 
  geom_arc_bar(
    data = results_long_melted_ALL %>%
      group_by(sample) %>%
      mutate(
        start_angle = cumsum(proportion_of_taxa * 2 * pi) - proportion_of_taxa * 2 * pi,
        end_angle = cumsum(proportion_of_taxa * 2 * pi),
        radius = 0.5
      ),
    aes(
      x0 = as.numeric(as.factor(sample)), 
      y0 = max(relative_abundance) + 2.5,  # Position above the bars
      r0 = 0,
      r = radius,
      start = start_angle,
      end = end_angle,
      fill = category
    ),
    inherit.aes = FALSE,
    color = NA  # Remove the border by setting the color to NA
  ) +
  coord_cartesian(clip = "off")
p

ggsave("stacked_piescharts_mapq60.png", p, dpi = 300, height = 7, width = 7, bg = "white")



# try 2

results_long_melted_ALL$relative_abundance<- results_long_melted_ALL$relative_abundance /100

# Convert to long format if necessary
results_long_melted_ALL_melted <- melt(results_long_melted_ALL, id.vars = c("sample", "category"), 
                                       measure.vars = c("relative_abundance", "proportion_of_taxa"),
                                       variable.name = "measurement", value.name = "value")

q<- ggplot(results_long_melted_ALL_melted, aes(x = sample, y = value, fill = category)) +
  geom_bar(stat = "identity", alpha = 0.7) + #colour="black"
  facet_wrap(~ measurement, scales = "free_y", nrow = 2) +
  labs(x = "Sample", y = NULL, fill = "Category") +
  scale_fill_manual(values = c("other" = "#E1BE6A", "putative_inhibitory_taxa" = "#40B0A6")) +
  scale_color_manual(values = c("relative_abundance" = "black", "proportion" = "black")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "right",  # Adjust as needed (or use "none" to hide all legends)
    legend.title = element_blank(),  # Hides the legend title if needed
    strip.text = element_text(size = 12, hjust = 0),  # Align facet labels to the left
    strip.background = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "white")
  ) +
  scale_y_continuous(labels = scales::comma)

q

ggsave("stacked_twopanels_mapq60.png", q, dpi = 300, height = 8, width = 10, bg = "white")


