rm(list=ls()) # clean the memeory
library("seqinr")
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

ctl_binded_combinded <- do.call("rbind", list(ctl_A0201_binded))
ctl_binded_9mer <- ctl_binded_combinded[nchar(ctl_binded_combinded$peptide)==9, ]

ctl_binded_9mer_immunogen <- ctl_binded_9mer[ctl_binded_9mer$Immunogenicity > 0.6, ]
write.fasta(as.list(ctl_binded_9mer_immunogen$peptide), names = 1:nrow(ctl_binded_9mer_immunogen), file.out = "epitopes_immunogen.fasta")

ctl_binded_9mer_nonimmunogen <- ctl_binded_9mer[ctl_binded_9mer$Immunogenicity < 0.4, ]
write.fasta(as.list(ctl_binded_9mer_nonimmunogen$peptide), names = 1:nrow(ctl_binded_9mer_immunogen), file.out = "epitopes_nonimmunogen.fasta")
