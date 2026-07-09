# ref: serotype O https://www.ebi.ac.uk/ena/browser/view/DQ164925
# ref: serotype A https://www.ebi.ac.uk/ena/browser/view/HQ116312

library(here)
library(tidyverse)
# BiocManager::install("pwalign")
# BiocManager::install("DECIPHER")
library(Biostrings)


# Load sequences
sequences <- readDNAStringSet(here("data/sequences-aligned.fasta"))
mean(width(sequences)) # average number of nucleotides

# Load VP1 reference sequence
vp1_O_ref <- readDNAStringSet(here("data/vp1_DQ164925.1.fasta"))
reference_width <- width(vp1_O_ref) # number of nucleotides
reference_width

# Perform pairwise alignment for each sequence
O_alignments <- lapply(sequences, function(seq) {
  pwalign::pairwiseAlignment(pattern = vp1_O_ref, subject = seq, type = "global")
})

# Inspect alignment scores
alignment_scores <- sapply(O_alignments, function(align) score(align))

# Filter sequences with high alignment scores (threshold can be adjusted)
threshold <- quantile(alignment_scores, 0.5)  # Top 50% alignment scores

ggplot(data = data.frame(score = alignment_scores), aes(x = score)) +
  geom_histogram(bins = 20, fill = "gray60", alpha = 0.7, color = "black") + # Histogram
  geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 1) + # Vertical line
  labs(title = " ", x = "Alignment Score", y = "Frequency") +
  theme_minimal()

# Filter tothose above threshold
vp1_matches <- sequences[alignment_scores >= threshold]

all_widths <- mean(width(vp1_matches))
all_widths # these sequences have the VP1 based on the threshold, but have other genes as well


# remove everything other than VP1
trimmed_sequences <- subseq(vp1_matches, start = 1, end = reference_width)
widths <- width(trimmed_sequences )
summary(widths)

DECIPHER::browseSeqs(trimmed_sequences) 
View(as.data.frame(trimmed_sequences))

# Save the VP1-matching sequences
writeXStringSet(vp1_matches, here("data/vp1_sequences.fasta"))
