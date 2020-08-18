# this R file is used to calculate the MHC-peptide binding affinity using the algorithm provided by the IEDB website
rm(list=ls()) # clean the memeory
library(dplyr) # we will use the mutate_all function of this library
library(seqinr) # need the write.fasta function of this library

all_content <- readLines("ctl_variant_pruned.csv") # this is the pruned los alamos database, which contains HIV epitopes and its variants
skip_second <- all_content[-2] # remove the date line
ctl_variant_pruned <- read.csv(textConnection(skip_second), stringsAsFactors = FALSE) # it's important to not take strings as vectors


ctl_variant_pruned_A0201 <- ctl_variant_pruned[grepl("A2(,|\\.|$)|A\\*02", ctl_variant_pruned$HLA) & (!grepl("A\\*020[2-9]", ctl_variant_pruned$HLA)), ] # grab the epitopes bound to the A0201 HLA type
rownames(ctl_variant_pruned_A0201) <- 1:nrow(ctl_variant_pruned_A0201) # reorder the row numbers

ctl_variant_pruned_A0201_SeqImmu <- ctl_variant_pruned_A0201[, c(15,23)] # only take the sequence and the immunogenicity score as the necessary data

wt_epitopes <- names(table(ctl_variant_pruned_A0201[,14])) # take out all the wild type epitopes 
wt_immunogenicity <- rep(1, length(wt_epitopes)) # all the wild type epitopes are assumed to have immunogenicity score 1

ctl_wt_A0201_SeqImmu <- data.frame("Variant.Epitope" = wt_epitopes, "Mutant.Immunogenicity" = wt_immunogenicity, stringsAsFactors = FALSE) # construct the dataframe for the wild type epitopes

ctl_A0201_SeqImmu_combined <- rbind(ctl_variant_pruned_A0201_SeqImmu, ctl_wt_A0201_SeqImmu) # combine the dataframe of the wt and variant epitopes

ctl_A0201_SeqImmu_combined <- mutate_all(ctl_A0201_SeqImmu_combined, toupper) # change all the lowercase letters to uppercase
ctl_A0201_SeqImmu_combined$Mutant.Immunogenicity <- as.numeric(ctl_A0201_SeqImmu_combined$Mutant.Immunogenicity) # make sure the immunogenicity is of numeric type
ctl_A0201_SeqImmu_combined_reduced <- aggregate(Mutant.Immunogenicity ~ Variant.Epitope, data = ctl_A0201_SeqImmu_combined,  FUN = mean) # if some sequences appear multiple times, merge these sequences and take an average of the immunogenicity score

write.fasta(strsplit(ctl_A0201_SeqImmu_combined_reduced$Variant.Epitope, split = ""), names = 1:nrow(ctl_A0201_SeqImmu_combined_reduced), file.out = "ctl_A0201_Seq.fasta") # write out the sequence in fasta format 

A0201_IC50_results <- read.csv("result_A0201_netMHCpan.csv", stringsAsFactors = FALSE) # read in the MHC binding results of the peptides
A0201_IC50_results <- A0201_IC50_results[order(A0201_IC50_results$seq_num, A0201_IC50_results$rank),] 
A0201_IC50_results_reduced <- A0201_IC50_results[!duplicated(A0201_IC50_results$seq_num),] # only get the strongest binding resutlt for each peptide
A0201_IC50_results_reduced$Immunogenicity <- ctl_A0201_SeqImmu_combined_reduced$Mutant.Immunogenicity # combine IC50 with immunogenicity

A0201_results_outlier <- A0201_IC50_results_reduced[(A0201_IC50_results_reduced$ic50 > 1000) & (A0201_IC50_results_reduced$Immunogenicity > 0.25), ] # some peptides bind weakly to a MHC, but still able to elicit a response
write.csv(A0201_IC50_results_reduced[, c(1,6,7,8,9)], "ctl_A0201_SeqImmuIC50.csv", row.names = FALSE) # write out the results for all A0201 peptides

ctl_variant_pruned_B5701 <- ctl_variant_pruned[(grepl("B\\*57|B57", ctl_variant_pruned$HLA) & (!grepl("B\\*570[2-9]", ctl_variant_pruned$HLA))) | grepl("B\\*5701", ctl_variant_pruned$HLA), ]  # only look at allels which has B5701, notice that sometimes B5701 and B5703 exist simultaneously. That's why the logical sentence is long
rownames(ctl_variant_pruned_B5701) <- 1:nrow(ctl_variant_pruned_B5701) # reorder the row numbers

ctl_variant_pruned_B5701_SeqImmu <- ctl_variant_pruned_B5701[, c(15,23)] # only take the sequence and the immunogenicity score
wt_epitopes_B5701 <- names(table(ctl_variant_pruned_B5701[,14])) # take out all the wild type epitopes
wt_immunogenicity_B5701 <- rep(1, length(wt_epitopes_B5701)) # all the wild type epitopes are assumed to have immunogenicity score 1
ctl_wt_B5701_SeqImmu <- data.frame("Variant.Epitope" = wt_epitopes_B5701, "Mutant.Immunogenicity" = wt_immunogenicity_B5701, stringsAsFactors = FALSE)

ctl_B5701_SeqImmu_combined <- rbind(ctl_variant_pruned_B5701_SeqImmu, ctl_wt_B5701_SeqImmu) # combine the dataframe of mutant and wt epitopes

ctl_B5701_SeqImmu_combined <- mutate_all(ctl_B5701_SeqImmu_combined, toupper)
ctl_B5701_SeqImmu_combined$Mutant.Immunogenicity <- as.numeric(ctl_B5701_SeqImmu_combined$Mutant.Immunogenicity)
ctl_B5701_SeqImmu_combined_reduced <- aggregate(Mutant.Immunogenicity ~ Variant.Epitope, data = ctl_B5701_SeqImmu_combined, FUN = mean) # merge the sequences that appear multiple times
write.fasta(strsplit(ctl_B5701_SeqImmu_combined_reduced$Variant.Epitope, split = ""), names = 1:nrow(ctl_B5701_SeqImmu_combined_reduced), file.out = "ctl_B5701_Seq.fasta") # write out the sequence in fasta format 

B5701_IC50_results <- read.csv("result_B5701_netMHCpan.csv", stringsAsFactors = FALSE)
B5701_IC50_results <- B5701_IC50_results[order(B5701_IC50_results$seq_num, B5701_IC50_results$rank), ]
B5701_IC50_results_reduced <- B5701_IC50_results[!duplicated(B5701_IC50_results$seq_num),]
B5701_IC50_results_reduced$Immunogenicity <- ctl_B5701_SeqImmu_combined_reduced$Mutant.Immunogenicity

B5701_results_outlier <- B5701_IC50_results_reduced[(B5701_IC50_results_reduced$ic50 > 1000) & (B5701_IC50_results_reduced$Immunogenicity > 0.25), ]
write.csv(B5701_IC50_results_reduced[, c(1,6,7,8,9)], "ctl_B5701_SeqImmuIC50.csv", row.names = FALSE)


#ctl_variant_pruned_A0301 <- ctl_variant_pruned[grepl("A3(,|\\.|$)|A\\*03", ctl_variant_pruned$HLA), ]
#rownames(ctl_variant_pruned_A0301) <- 1:nrow(ctl_variant_pruned_A0301)

#ctl_variant_pruned_A0301_SeqImmu <- ctl_variant_pruned_A0301[, c(15,23)]
#wt_epitopes_A0301 <- names(table(ctl_variant_pruned_B[,14]))
#wt_immunogenicity_B5701 <- rep(1, length(wt_epitopes_B5701))