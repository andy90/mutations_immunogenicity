# for this model, I just implement the simple hydrophobicity only criteria
rm(list = ls())
library(tidyverse)
library(foreach)
library(doParallel)
registerDoParallel(cores=2)
aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W", "-")
p_Sdandinski <- c(-0.2, 0.0, -1, -0.6, -1.2, -0.4, -0.5, -0.05, -0.8, 1.5, 0.0, 0.0, 0.1, -0.1, 0.05, 0.0, 0.2, -0.05, 0.4, 0.8, 1.2, 0)
p_wif <- c(0.81, 0.17, 0.99, 1.23, 2.02, 0.13, 0.14, 0.42, 0.58, -0.24, 0, 0.01, 0.45, 0.17, 0.07, -0.31, -0.56, -0.23, -1.13, -0.94, -1.85, 0)

to_hydro <- function(pep){
  aas <- strsplit(pep, "")[[1]]
  sapply(aas, function(aa){
    ind <- which(aa_name == aa)
    p_wif[ind]
  })
}

df_positive <- readRDS("df_positive")
df_negative <- readRDS("df_negative")


positive_CDR3_hydro <-
  foreach(pep = df_positive$CDR3_aligned, .combine = rbind) %dopar% {
    to_hydro(pep)
  }

row.names(positive_CDR3_hydro) <- c()
colnames(positive_CDR3_hydro) <- paste0("T", 1:19)

positive_Epitope_hydro <-
  foreach(pep = df_positive$Epitope_aligned, .combine = rbind) %dopar% {
    to_hydro(pep)
  }

row.names(positive_Epitope_hydro) <- c()
colnames(positive_Epitope_hydro) <- paste0("E", 1:11)


negative_CDR3_hydro <-
  foreach(pep = df_negative$CDR3_aligned, .combine = rbind) %dopar% {
    to_hydro(pep)
  }


row.names(negative_CDR3_hydro) <- c()
colnames(negative_CDR3_hydro) <- paste0("T", 1:19)

negative_Epitope_hydro <-
  foreach(pep = df_negative$Epitope_aligned, .combine = rbind) %dopar% {
    to_hydro(pep)
  }

row.names(negative_Epitope_hydro) <- c()
colnames(negative_Epitope_hydro) <- paste0("E", 1:11)

df_hydro_positive <- data.frame(positive_CDR3_hydro, positive_Epitope_hydro, score = 1)
df_hydro_negative <- data.frame(negative_CDR3_hydro, negative_Epitope_hydro, score = 0)

df_hydro <- rbind(df_hydro_positive, df_hydro_negative)

saveRDS(df_hydro, file = "df_hydro")  
