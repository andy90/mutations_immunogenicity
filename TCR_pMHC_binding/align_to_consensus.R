rm(list = ls())
library(tidyverse)
library(seqinr)
library("Biostrings")
data("BLOSUM62")

IEDB_11mer_consensus <- "RLLDLSGLLLL"
CDR3_19mer_consensus <- toupper(read.table("CDR3_19mer_consensus_1.txt")$V1)

##############
# function used to perform global alignment
##############
blosum_align <- function(pep1, pep2){
  res <- pairwiseAlignment(pep1, pep2, type = "global", substitutionMatrix = BLOSUM62,
                           gapOpening = 11, gapExtension = 10) # see wikipedia for the detail of the Smith-Waterman alignment
  res
}
##############
#############

input_file <- "VDJdb_human_TRB.tsv"
df_TRB_peptide <- read_tsv(input_file)

df_TRB_peptide_ClassI <-
  df_TRB_peptide %>%
  filter(`MHC class` == "MHCI") %>%
  mutate(lseq = nchar(CDR3)) %>% # make a cutoff at 19 is more than capable
  filter(lseq < 20) 

CDR3_aligned <- as.character(blosum_align(df_TRB_peptide_ClassI$CDR3, CDR3_19mer_consensus))

df_MHC_Epitope <-
  df_TRB_peptide_ClassI %>%
  distinct_at(vars( Epitope)) %>%
  mutate(lseq = nchar(Epitope)) %>%
  filter(lseq < 12)

Epitope_aligned <- as.character(blosum_align(df_MHC_Epitope$Epitope, IEDB_11mer_consensus))

df_MHC_Epitope_aligned <-
  df_MHC_Epitope %>%
  mutate(Epitope_aligned = Epitope_aligned)

df_CDR3_peptide_ClassI_aligned <- 
  df_TRB_peptide_ClassI %>%
  mutate(CDR3_aligned = CDR3_aligned) %>%
  filter(nchar(Epitope) < 12) %>%
  mutate(Epitope_aligned = sapply(Epitope, function(pep){
    ind <- which(df_MHC_Epitope_aligned$Epitope == pep)
    df_MHC_Epitope_aligned$Epitope_aligned[ind]
  })) %>%
  select(CDR3, CDR3_aligned, Epitope, Epitope_aligned)

saveRDS(df_CDR3_peptide_ClassI_aligned, "df_CDR3_peptide_ClassI_aligned")
