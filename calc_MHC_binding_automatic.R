# this R file is used to calculate the MHC-peptide binding affinity using the algorithm provided by the IEDB website, the process is automaticed for all the HLA allels
rm(list=ls()) # clean the memeory
library(dplyr) # we will use the mutate_all function of this library
library(seqinr) # need the write.fasta function of this library

all_content <- readLines("ctl_variant_pruned.csv") # this is the pruned los alamos database, which contains HIV epitopes and its variants
skip_second <- all_content[-2] # remove the date line
ctl_variant_pruned <- read.csv(textConnection(skip_second), stringsAsFactors = FALSE) # it's important to not take strings as vectors

#allele_index <- grepl("A3(,|\\.|$)|A\\*03", ctl_variant_pruned$HLA) # this is to get the index of the data entries for the specific allele
#allele_name <- "A0301"

#allele_index <- grepl("B8(,|\\.|$)|B\\*08", ctl_variant_pruned$HLA)
#allele_name <- "B0801"

#allele_index <- grepl("B27(,|\\.|$)|B\\*27", ctl_variant_pruned$HLA)
#allele_name <- "B2705"

#allele_index <- grepl("A24(,|\\.|$)|A\\*24", ctl_variant_pruned$HLA)
#allele_name <- "A2402"

#allele_index <- grepl("A2(,|\\.|$)|A\\*02", ctl_variant_pruned$HLA) & (!grepl("A\\*020[2-9]", ctl_variant_pruned$HLA))
#allele_name <- "A0201"

#allele_index <- (grepl("B\\*57|B57", ctl_variant_pruned$HLA) & (!grepl("B\\*570[2-9]", ctl_variant_pruned$HLA))) | grepl("B\\*5701", ctl_variant_pruned$HLA)
#allele_name <- "B5701"

allele_index <- grepl("B7(,|\\.|$| )|B\\*07", ctl_variant_pruned$HLA)
allele_name <- "B0702"
ctl_variant_pruned_HLA <- ctl_variant_pruned[allele_index, ]  # this is the data entries for the specific allele
rownames(ctl_variant_pruned_HLA) <- 1:nrow(ctl_variant_pruned_HLA) # reorder the row numbers

ctl_variant_pruned_HLA_SeqImmu <- select(ctl_variant_pruned_HLA, Variant.Epitope, Mutant.Immunogenicity) 

wt_epitopes <- names(table(select(ctl_variant_pruned_HLA, Epitope))) # take out all the non repeated epitopes
wt_immunogenicity <- rep(1, length(wt_epitopes)) # all the wild type epitopes are assumed to have immunogenicity score 1

ctl_wt_SeqImmu <- data.frame("Variant.Epitope" = wt_epitopes, "Mutant.Immunogenicity" = wt_immunogenicity, stringsAsFactors = FALSE) # construct the dataframe for the wild type epitopes
ctl_SeqImmu_combined <- rbind(ctl_variant_pruned_HLA_SeqImmu, ctl_wt_SeqImmu) # combine the dataframe of the wt and variant epitopes

ctl_SeqImmu_combined <- mutate_all(ctl_SeqImmu_combined, toupper) # change all the lowercase letters to uppercase
ctl_SeqImmu_combined$Mutant.Immunogenicity <- as.numeric(ctl_SeqImmu_combined$Mutant.Immunogenicity) # make sure the immunogenicity is of numeric type
ctl_SeqImmu_combined_reduced <- aggregate(Mutant.Immunogenicity ~ Variant.Epitope, data = ctl_SeqImmu_combined,  FUN = mean) # if some sequences appear multiple times, merge these sequences and take an average of the immunogenicity score

if (allele_name == "B2705"){
  ctl_SeqImmu_combined_reduced <- ctl_SeqImmu_combined_reduced[-c(21), ]
}
write.fasta(strsplit(ctl_SeqImmu_combined_reduced$Variant.Epitope, split = ""), names = 1:nrow(ctl_SeqImmu_combined_reduced), file.out = paste("ctl_",allele_name,"_Seq.fasta", sep="")) # write out the sequence in fasta format 

# now we need to upload the sequence to the iedb website to get the binding affinity
method_name <- "netMHCpan" # the method we are gonna choose
IC50_results <- read.csv(paste("result_",allele_name,"_",method_name, ".csv", sep=""), stringsAsFactors = FALSE) # read in the MHC binding results of the peptides
IC50_results <- IC50_results[order(IC50_results$seq_num, IC50_results$rank),]
IC50_results_reduced <- IC50_results[!duplicated(IC50_results$seq_num),] # only get the strongest binding resutlt for each peptide
IC50_results_reduced$Immunogenicity <- ctl_SeqImmu_combined_reduced$Mutant.Immunogenicity # combine IC50 with immunogenicity

IC50_results_reduced_reduced <- aggregate(cbind(ic50, rank, Immunogenicity) ~ peptide, data = IC50_results_reduced, FUN = mean) # if after processing some peptides again repeats, we want to merge them
IC50_results_reduced_reduced$allele = rep(allele_name, nrow(IC50_results_reduced_reduced))

results_outlier <- IC50_results_reduced_reduced[(IC50_results_reduced_reduced$ic50 > 1000) & (IC50_results_reduced_reduced$Immunogenicity > 0.25), ] # some peptides bind weakly to a MHC, but still able to elicit a response
write.csv(IC50_results_reduced_reduced, paste("ctl_", allele_name, "_SeqImmuIC50.csv", sep=""), row.names = FALSE) # write out the results for all A0201 peptides
