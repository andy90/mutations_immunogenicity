rm(list = ls())
library(foreach)
library(doParallel)
registerDoParallel(cores=2)
df_CDR3_peptide_ClassI_aligned <- readRDS("df_CDR3_peptide_ClassI_aligned")

df_positive <- 
  df_CDR3_peptide_ClassI_aligned %>%
  select(CDR3_aligned, Epitope_aligned)

df_negative <-
  foreach(i = 1:(nrow(df_positive) + 20000), .combine = rbind) %dopar%{
    iCDR3 <- sample(df_positive$CDR3_aligned, 1)
    iEpitope <- sample(df_positive$Epitope_aligned, 1)
    
    c(iCDR3, iEpitope)
  } 

df_negative <- 
  data.frame(df_negative, stringsAsFactors = FALSE) %>%
  set_names(c("CDR3_aligned", "Epitope_aligned")) %>%
  remove_rownames()

duplicated <- do.call(paste0, df_negative) %in% do.call(paste0, df_positive)

df_negative_nooverlap <-
  df_negative[!duplicated,]

saveRDS(df_negative_nooverlap, file = "df_negative")
saveRDS(df_positive, file = "df_positive")
