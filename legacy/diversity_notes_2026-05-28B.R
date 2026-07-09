#packages
library(here) #directory management
library(tidyverse) #general wrangling
library(ape) #phylodata manipulation
library(phangorn) #phylodata manipulation
library(pegas) # phylodata manipulation
library(Biostrings) #phylodata manipulation
library(msa) #phylodata manipulation
library(adegenet) #DAPC package

# OSF
library(osfr) # data management and archiving


##for spatial analyses and disaggregation
library(sp) #spatial data manipulation
library(spdep) #spatial data manipulation
library(terra) #spatial data manipulation
library(sf) #spatial data manipulation
library(geosphere) #distance measurements
library(vegan) #Mantel tests and Redundancy Analysis
library(ggforce) #help with pie chart plots
library(ggpubr) #help with pie chart plots
library(viridis) #color pallete
library(ggtree) #plot trees
library(patchwork) #combine plots
library(cowplot) #pull legend from ggplot
library(geodata) #import spatial covaraites
library(INLA) #model inference
library(Matrix) #diagnostics to check matrices


## Read Alignments and Data  

#sequence metadata: sequence accession names with associated location/date/host metadata
load(here("data/prep-for-seqs.Rda")) #data object of interest: seq_meta
rm(seqs.g23, seqs.csv1, genbank_return1, genbank_return2, seq_meta_orig, prov_reg, seq_meta_reg)

#clean trimmed VP1 alignments
#might need to revisit outlier removal process eventually, but these should work for now
alignment.SEA<-read.dna(here("data/vp1_O-SEA_trimmed.fasta"),
                        format="fasta", as.matrix=TRUE)

alignment.Cath<-read.dna(here("data/vp1_O-Cath_trimmed.fasta"),
                         format="fasta", as.matrix=TRUE)

alignment.MESA<-read.dna(here("data/vp1_O-MESA_trimmed.fasta"),
                         format="fasta", as.matrix=TRUE)

alignment.A<-read.dna(here("data/vp1_A_trimmed.fasta"),
                      format="fasta", as.matrix=TRUE)

#Vietnam shapefile: country shapefile with provinces in equal area projection
load(here("data/vn-spat.Rda")) #data object of interest: vn_provinces_eq
rm(vn_provinces, vn_spat1, vn_prov_reg1)

#equal area projection for reference
equal_area_proj<-"+proj=aea +lat_1=20 +lat_2=40 +lat_0=30 +lon_0=105 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"

metadata <- seq_meta

head(metadata)

target_accessions <- c("KT153164", "ON868620", "ON868617", "ON868609", "ON868605")
image(alignment.SEA[target_accessions,]) #whack O's


alignment_list <- list(
  Cath = alignment.Cath,
  SEA = alignment.SEA,
  MESA = alignment.MESA,
  A = alignment.A
)

# 1. Pull the raw sequences and convert them strictly to UPPERCASE
raw_sequences_list <- lapply(alignment_list, function(mat) {
  sapply(1:nrow(mat), function(i) {
    # as.character() creates lowercase; toupper() forces it to uppercase
    toupper(paste(as.character(mat[i, ]), collapse = ""))
  })
})

# 2. Combine into a single master character vector with names
all_raw_strings <- unlist(raw_sequences_list)
names(all_raw_strings) <- unlist(lapply(alignment_list, rownames))

# 3. Convert to a single DNAStringSet (Now properly formatted in uppercase)
master_stringset <- DNAStringSet(all_raw_strings)

# 4. Translate to Amino Acids (This will now run smoothly)
# we use if_incomplete="drop" to handle the trailing partial codon seamlessly
# master_aa <- translate(master_stringset, if_incomplete = "drop")

# 5. Run the progressive alignment on the clean protein frames
message("Aligning protein frames...")
global_nucleotide_msa <- msa(master_stringset, method = "ClustalW")

master_alignment_fixed <- as.DNAbin(global_nucleotide_msa)

# 3. Check the dimensions of your clean, global matrix
dim(master_alignment_fixed)

lineage_mapping <- stack(lapply(alignment_list, rownames)) %>%
  rename(accession = values, lineage = ind)

# 2. Join this back to your metadata
working_metadata <- metadata %>%
  inner_join(lineage_mapping, by = "accession")

combinations <- working_metadata %>%
  group_by(location, lineage) %>%
  tally() %>%
  filter(n >= 2) # Diversity requires at least 2 sequences


# Initialize a dataframe to capture the clean results
stratified_diversity_fixed <- data.frame()

# Run the calculation loop on the new master alignment
for(i in 1:nrow(combinations)) {
  prov <- combinations$location[i]
  lin  <- combinations$lineage[i]
  n_seq <- combinations$n[i]
  
  # Extract accessions for this specific lineage/province group
  target_accessions <- working_metadata %>%
    filter(location == prov, lineage == lin) %>%
    pull(accession)
  
  # Subset the fixed master alignment matrix
  sub_alignment <- master_alignment_fixed[target_accessions, , drop = FALSE]
  
  # Calculate Pegas nucleotide diversity (pi)
  pi_val <- nuc.div(sub_alignment)
  
  # Append to results
  stratified_diversity_fixed <- rbind(stratified_diversity_fixed, data.frame(
    Province = prov,
    Lineage = lin,
    Sample_Size = n_seq,
    Nucleotide_Diversity = pi_val
  ))
}

# View your corrected, finalized diversity breakdown
print(stratified_diversity_fixed)


diversity_wide <- stratified_diversity_fixed %>%
  select(Province, Lineage, Nucleotide_Diversity) %>%
  mutate(Lineage = paste0("pi_", Lineage)) %>% 
  pivot_wider(names_from = Lineage, values_from = Nucleotide_Diversity)

# 2. Check alignment of province names with your sf object
# (Replace 'Province_Name_Field' with the actual column name in vn_provinces_eq)
shapefile_names <- vn_provinces_eq$prov_eng

mismatches <- setdiff(diversity_wide$Province, shapefile_names)
if(length(mismatches) > 0) {
  cat("The following province names do not match your shapefile:\n")
  print(mismatches)
} else {
  cat("All province names align perfectly with the shapefile!\n")
}
