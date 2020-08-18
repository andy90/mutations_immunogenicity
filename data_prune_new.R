# this R file is used to prune out the unwanted data entries from the original epitope CTL response dataset
rm(list=ls()) # clean the memeory
all_content <- readLines("ctl_variant.csv") # this is the los alamos database, which contains HIV epitopes and its variants
skip_second <- all_content[-2] # remove the date line
ctl_variant <- read.csv(textConnection(skip_second), stringsAsFactors = FALSE) # it's important to not take strings as vectors

subset_methods <- grepl("IFN|Elispot|Chromium|Tetramer|cytokine", ctl_variant$Methods) # we only want the data obtained by IFNgamma-Elispot, Chromium release assay, Tetramer binding, Intracellular cytokine staining and Flow cytometric T-cell cytokine assay
subset_mutcode <- grepl("(^|,|\\s+)(S?)SF|(^|,|\\s+)(S?)NSF|(^|,|\\s+)E($|,)|(^|,|\\s+)LE($|,)||(^|,|\\s+)IE($|,)TCR|DR", ctl_variant$Mutation.Type.Code) # the mutation type code has to be one of the following: E, NSF, NSF-2, SF, SNSF, SSF, TCR, DR
subset_sequence <- !grepl("\\{|\\||\\*|[:punct:]|-|X|x|B|b|J|j|O|o|Z|z|U|u", ctl_variant$Variant.Epitope) # some epitope sequences contain special characters like {}, *, we want to take out those entries
subset_pruned <- subset_methods & subset_mutcode & subset_sequence # this is the data entries we want

ctl_variant_pruned <- ctl_variant[subset_pruned, ] # this is the pruned data set
rownames(ctl_variant_pruned) <- 1:nrow(ctl_variant_pruned) # reorder the row number from 1

function_code_to_immunogenicity <- function(row_of_ctl){ # this function is used to transfer the mutatant code to a immunogenicity score which measures the strength of the mutation epitope - CTL response
  codes <- strsplit(as.character(row_of_ctl[17]), split=",")
  total_score <- sum(grepl("(^|,|\\s+)(S?)SF", codes[[1]]) + 0.5 * grepl("DR", codes[[1]])) # total immune scores carried by the mutation type code, SF or SSF is taken as 1, DR is taken as 0.5
  total_code <- sum(grepl("(^|,|\\s+)(S?)SF|(^|,|\\s+)(S?)NSF|(^|,|\\s+)E($|,)|(^|,|\\s+)LE($|,)||(^|,|\\s+)IE($|,)TCR|DR", codes[[1]])) # the total number of codes related to the change of immunogenicity
  score <- total_score/total_code # this is the averaged immunogenecity score, for example, one epitope can elicit response in some patients and do not elicit in others, that's why the score for the epitope should be an average
}

immune_score <- apply(ctl_variant_pruned, 1, function_code_to_immunogenicity) # calculate the immune score for each data entry
ctl_variant_pruned["Mutant.Immunogenicity"] <- immune_score # add the immune score to the data frame

ctl_variant_immuscore <- data.frame("peptide" = toupper(ctl_variant_pruned$Variant.Epitope), "Immunogenicity" = ctl_variant_pruned$Mutant.Immunogenicity, stringsAsFactors = FALSE)

##############################
# now we read in all the wild-type epitopes 
#############################
all_content <- readLines("ctl_summary.csv") # this is the los alamos database, which contains HIV epitopes and its variants
skip_second <- all_content[-2] # remove the date line
ctl_wt <- read.csv(textConnection(skip_second), stringsAsFactors = FALSE) # it's important to not take strings as vectors
ctl_wt_pruned <- ctl_wt[subset_pruned, ]
ctl_wt_immuscore <- data.frame("peptide" = toupper(ctl_wt_pruned$Epitope), "Immunogenicity" = rep(1, nrow(ctl_wt_pruned)), stringsAsFactors = FALSE)



#######################
# construct a total data frame
#######################
ctl_wt_variant_immuscore <- rbind(ctl_variant_immuscore, ctl_wt_immuscore)

ctl_wt_variant_immuscore_reduced <- aggregate(Immunogenicity ~ peptide, data = ctl_wt_variant_immuscore,  FUN = mean) # if some sequences appear multiple times, merge these sequences and take an average of the immunogenicity score
# there are about 700 non repeated entries of 9 mers
sum((nchar(ctl_wt_variant_immuscore_reduced$peptide) == 9) & (ctl_wt_variant_immuscore_reduced$Immunogenicity > 0.7))
sum((nchar(ctl_wt_variant_immuscore_reduced$peptide) == 9) & (ctl_wt_variant_immuscore_reduced$Immunogenicity < 0.3)) # there are alot of non immunogenic peptides in this data set

write.csv(ctl_wt_variant_immuscore_reduced, file = "ctl_wt_variant_immu.csv")

