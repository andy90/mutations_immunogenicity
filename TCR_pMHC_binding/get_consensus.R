library(tidyverse)
library(seqinr)
df_IEDB_binders <- read_csv(file = "IEDB_HLAI_epitopes.csv", skip = 1)

df_IEDB_11mers <- 
  df_IEDB_binders %>%
  filter(nchar(Description) == 11)

write.fasta(as.list(df_IEDB_11mers$Description), names = 1:length(df_IEDB_11mers), file.out = "IEDB_11mers.fasta")

IEDB_11mer_consensus <- "RLLDLSGLLLL"

