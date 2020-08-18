rm(list=ls()) # clean the memeory

ctl_A0201 <- read.csv("ctl_A0201_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele A0201
ctl_B5701 <- read.csv("ctl_B5701_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele B5701
ctl_A0301 <- read.csv("ctl_A0301_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele A0301
ctl_B0702 <- read.csv("ctl_B0702_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele B0702

# only work on those peptides that bind to the MHC
bind_criteria <- 500 # the binding IC50 is 500nM
ctl_A0201_binded <- ctl_A0201[ctl_A0201$ic50 < bind_criteria, ]
ctl_B5701_binded <- ctl_B5701[ctl_B5701$ic50 < bind_criteria, ]
ctl_A0301_binded <- ctl_A0301[ctl_A0301$ic50 < bind_criteria, ]
ctl_B0702_binded <- ctl_B0702[ctl_B0702$ic50 < bind_criteria, ]

ctl_binded_combinded <- do.call("rbind", list(ctl_A0201_binded, ctl_A0301_binded, ctl_B0702_binded, ctl_B5701_binded))

ctl_binded_9mer <- ctl_binded_combinded[nchar(ctl_binded_combinded$peptide)==9, ]
ctl_binded_9mer$Immunogenicity <- ctl_binded_9mer$Immunogenicity + 0.05
ctl_binded_9mer$copies <- round(ctl_binded_9mer$Immunogenicity * 20)

file_epitopes <- file("epitopes.fasta")
index <- 1
lines_to_write <- c()
epitopes_repeated <- c()
for (i_epi in 1:nrow(ctl_binded_9mer)){
  i_row <- ctl_binded_9mer[i_epi, ]
  for (i in 1:i_row$copies){
    #writeLines(c(paste(">",as.character(index),sep=""), i_row$peptide), file_epitopes)
    lines_to_write <- c(lines_to_write, paste(">",as.character(index),sep=""), i_row$peptide)
    epitopes_repeated <- c(epitopes_repeated, i_row$peptide)
    index <- index + 1
  }
}
writeLines(lines_to_write, file_epitopes)
close(file_epitopes)

energys <- read.table("epitopes-energy.txt")
epitopes_E <- data.frame("epitopes" = epitopes_repeated, "E" = energys)
epitopes_E_nodup <- epitopes_E[!duplicated(epitopes_E$epitopes), ]
epitopes_E_nodup$Immunogenicity <- ctl_binded_9mer$Immunogenicity
