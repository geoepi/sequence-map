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






# Pool all sequence groups into one master DNAbin matrix
master_alignment <- do.call(rbind, alignment_list)








image(alignment.SEA[target_accessions,]) #whack O's



# 1. Correctly collapse the DNAbin matrix into 188 full-length strings
sea_char <- sapply(1:nrow(alignment.SEA), function(i) {
  paste(as.character(alignment.SEA[i, ]), collapse = "")
})

# Re-assign the sequence names so they aren't lost
names(sea_char) <- rownames(alignment.SEA)

# 2. Convert to DNAStringSet (Should read: "A DNAStringSet instance of length 188")
sea_stringset <- DNAStringSet(sea_char)
print(sea_stringset) 

# 3. Run ClustalW alignment (This will now run instantly)
sea_realignment <- msa(sea_stringset, method = "ClustalW")

# 4. Convert back to DNAbin
alignment.SEA_fixed <- as.DNAbin(sea_realignment)

# 5. Check the dimensions of your new alignment
dim(alignment.SEA_fixed)

# Slice out the core 553 positions
alignment.SEA_trimmed <- alignment.SEA_fixed[, 1:553]

# Quick sanity check on the new dimensions
dim(alignment.SEA_trimmed)
# Should return: [1] 188 553

target_accessions <- c("KT153164", "ON868620", "ON868617", "ON868609", "ON868605")

image(alignment.SEA_trimmed[target_accessions, ])

alignment_list$SEA <- alignment.SEA_trimmed

# 2. Re-pool the updated sequences into the master alignment matrix
master_alignment <- do.call(rbind, alignment_list)

# 3. Double check the total row count matches your original 430 total sequences
dim(master_alignment) 
# Should return: [1] 430 5
