rm(list = ls())
library(tidyverse)
library(here)
library(seqinr)
library(magrittr)
library(msa)
library("Biostrings")
data("BLOSUM62")

##############
# function used to perform global alignment
##############
blosum_align <- function(pep1, pep2){
  res <- pairwiseAlignment(pep1, pep2, type = "global", substitutionMatrix = BLOSUM62,
                           gapOpening = 11, gapExtension = 10) # see wikipedia for the detail of the Smith-Waterman alignment
  res
}

input_file <- here("VDJdb_human_TRB.tsv")
df_TRB_peptide <- read_tsv(input_file)


df_TRB_peptide %>%
  group_by(`MHC class`) %>%
  summarise(n_CDR3 = n_distinct(CDR3),
            n_Epitope = n_distinct(Epitope),
            n_V = n_distinct(V),
            n_J = n_distinct(J),
            n_HLA = n_distinct(`MHC A`) + n_distinct(`MHC B`),
            n_MHCClass = n_distinct(`MHC class`)) # most of the data are for MHC Class I, which is good


df_TRB_peptide %>%
  group_by(`MHC A`, `MHC B`) %>%
  summarise(n_CDR3 = n_distinct(CDR3),
            n_Epitope = n_distinct(Epitope),
            n_V = n_distinct(V),
            n_J = n_distinct(J),
            n_HLA = n_distinct(`MHC A`) + n_distinct(`MHC B`),
            n_MHCClass = n_distinct(`MHC class`)) # most of the data are for MHC Class I, which is good

unique(df_TRB_peptide$`MHC B`) # MHC B are class II allele and the so called B2M allele (which might not even be an allele)
unique(df_TRB_peptide$`MHC A`) # MHC A has both A allele and B allele and Class II allele
unique(df_TRB_peptide$`MHC class`)


df_TRB_peptide %>%
  filter(`MHC class` == "MHCI") %>%
  group_by(`MHC A`) %>%
  summarise(n_CDR3 = n_distinct(CDR3),
            n_Epitope = n_distinct(Epitope),
            n_V = n_distinct(V),
            n_J = n_distinct(J)) %>%
  arrange(desc(n_CDR3)) # the TCRs most focus on several HLAs


df_TRB_peptide_ClassI <-
  df_TRB_peptide %>%
  filter(`MHC class` == "MHCI") %>%
  mutate(lseq = nchar(CDR3)) %>% # make a cutoff at 19 is more than capable
  filter(lseq < 20) 

ref_CDR3 <- "CASSVQALLAGDWADTQYF"

CDR3_aligned <- as.character(blosum_align(df_TRB_peptide_ClassI$CDR3, ref_CDR3))
hist(nchar(CDR3_aligned))
df_MHC_Epitope <-
  df_TRB_peptide_ClassI %>%
  distinct_at(vars( Epitope)) %>%
  mutate(lseq = nchar(Epitope)) %>%
  filter(lseq < 12)
ref_Epitope <- "EPLPQGQLTAY"
Epitope_aligned <- as.character(blosum_align(df_MHC_Epitope$Epitope, ref_Epitope))

df_MHC_Epitope_aligned <-
  df_MHC_Epitope %>%
  mutate(Epitope_aligned = Epitope_aligned)

df_CDR3_peptide_ClassI_aligned <- 
  df_TRB_peptide_ClassI %>%
  mutate(CDR3_aligned = CDR3_aligned) %>%
  mutate(Epitope_aligned = sapply(Epitope, function(pep){
    ind <- which(df_MHC_Epitope_aligned$Epitope == pep)
    df_MHC_Epitope_aligned$Epitope_aligned[ind]
  })) %>%
  select(CDR3, CDR3_aligned, Epitope, Epitope_aligned)
  
saveRDS(df_CDR3_peptide_ClassI_aligned, "df_CDR3_peptide_ClassI_aligned")

#write.fasta(sequences = as.list(df_TRB_peptide_ClassI$CDR3), names = 1:length(df_TRB_peptide_ClassI$CDR3), file.out = "TRB_ClassI_CDR3.fasta")
#write.fasta(sequences = as.list(df_MHC_Epitope$Epitope), names = 1:length(df_MHC_Epitope$Epitope), file.out = "TRB_ClassI_Epitope.fasta")

df_TRB_peptide_ClassI_19mer <- 
  df_TRB_peptide_ClassI %>%
  filter(lseq == 19)

df_TRB_peptide_ClassI_18mer <- 
  df_TRB_peptide_ClassI %>%
  filter(lseq == 18)

write.fasta(sequences = as.list(df_TRB_peptide_ClassI_18mer$CDR3), names = 1:length(df_TRB_peptide_ClassI_18mer$CDR3), file.out = "TRB_ClassI_CDR3_1mer.fasta")

write.fasta(sequences = as.list(df_TRB_peptide_ClassI_19mer$CDR3), names = 1:length(df_TRB_peptide_ClassI_19mer$CDR3), file.out = "TRB_ClassI_CDR3_19mer.fasta")
a <- msa(df_MHC_Epitope$Epitope, gapOpening = 20, gapExtension = 15, substitutionMatrix = BLOSUM62, type = "protein")
 