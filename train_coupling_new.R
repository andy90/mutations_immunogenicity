# this r script is used to train a model which uses a real number to mark each amino acid
rm(list=ls()) # clean the memeory


###################
# the first step is to read in the peptides
##################

ctl_A0201 <- read.csv("ctl_A0201_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele A0201
ctl_B5701 <- read.csv("ctl_B5701_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele B5701
ctl_A0301 <- read.csv("ctl_A0301_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele A0301
ctl_B0702 <- read.csv("ctl_B0702_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele B0702

# only work on those peptides that bind to the MHC
bind_criteria <- 50000000 # now we just read in all the peptides
ctl_A0201_binded <- ctl_A0201[ctl_A0201$ic50 < bind_criteria, ]
ctl_B5701_binded <- ctl_B5701[ctl_B5701$ic50 < bind_criteria, ]
ctl_A0301_binded <- ctl_A0301[ctl_A0301$ic50 < bind_criteria, ]
ctl_B0702_binded <- ctl_B0702[ctl_B0702$ic50 < bind_criteria, ]

ctl_binded_combinded <- do.call("rbind", list(ctl_A0201_binded, ctl_A0301_binded, ctl_B0702_binded, ctl_B5701_binded))

###############
# we could use a different database of peptides
######
ctl_wt_variant <- read.csv("ctl_wt_variant_immu.csv", stringsAsFactors = FALSE) # this database collapsed all the alleles

#ctl_binded_combinded <- do.call("rbind", list(ctl_A0201_binded, ctl_A0301_binded, ctl_B0702_binded, ctl_B5701_binded))
ctl_binded_combinded <- ctl_wt_variant

ctl_binded_9mer <- ctl_binded_combinded[nchar(ctl_binded_combinded$peptide)==9, ]
ctl_binded_9mer_separated <- ctl_binded_9mer[(ctl_binded_9mer$Immunogenicity > 0.7) | (ctl_binded_9mer$Immunogenicity < 0.3), ]
ctl_binded_9mer_separated$Immunogenicity <- round(ctl_binded_9mer_separated$Immunogenicity)
ctl_binded_9mer_separated <- ctl_binded_9mer_separated[-1]
#################
# the second step is to define the function to calculate the score of aa
#################

score_of_aa <- function(aa, aa_name, p){ # returns the p score of each amino acid
  score <- p[aa == aa_name]
}


total_score <- function(row_of_ctl, x_xc, aa_name){ # returns the score of each peptide, for a 9mer, x_xc is of length 46
  peptide <- row_of_ctl[1]
  l_pep <- 9 # only focus on the 9mers in this case

  x <- x_xc[ 1 : (l_pep + 1)] # this is the coefficient for the single term, including the constant term
  x_coupling <- x_xc[(l_pep + 2 ) : length(x_xc)] # this is the coefficient for the coupling term
  
  scores <- sapply(strsplit(peptide, "")[[1]], score_of_aa, aa_name, p) # this is the scores for each peptide
  
  coupling_scores_matrix <- scores %o% scores # this is all the coupling term between the scores
  t_score <- sum(c(1,scores)*x) + sum(coupling_scores_matrix[upper.tri(coupling_scores_matrix)] * x_coupling)
}

cost_function <- function(x_xc, ctl_binded_combinded, lambda_x, lambda_xc){ # this calculate the total cost
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
  l_pep <- 9 # only focus on the 9mers in this case
  t_scores <- apply(ctl_binded_combinded, 1, total_score, x_xc, aa_name) # this is the total energy part Ec + x*p for each peptide
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
  
  y <- ctl_binded_combinded$Immunogenicity
  m <- length(y)
  x <- x_xc[ 2 : (l_pep + 1)] # this is the coefficient for the single term, excluding the constant term
  x_coupling <- x_xc[(l_pep + 2 ) : length(x_xc)] # this is the coefficient for the coupling term

  total_cost <- -mean(y*log(sigmoid_tscores) + (1 - y)*log(1 - sigmoid_tscores)) + lambda_x/(2*m)*sum(x * x) + lambda_xc/(2*m)*sum(x_coupling * x_coupling)
}

evaluate <- function(x_xc, ctl_binded_combinded){ # this evaluate the trained model
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
  
  t_scores <- apply(ctl_binded_combinded, 1, total_score, x_xc, aa_name) # this is the total energy part Ec + x*p for each peptide
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
}


####################
#read in p scores and starts training
#####################
aa_new_scale <- read.table(file = "combined_aa_scale.csv", stringsAsFactors = FALSE)
#p <- aa_new_scale$scale
p_Sdandinski <- c(-0.2, 0.0, -1, -0.6, -1.2, -0.4, -0.5, -0.05, -0.8, 1.5, 0.0, 0.0, 0.1, -0.1, 0.05, 0.0, 0.2, -0.05, 0.4, 0.8, 1.2)
p<- p_Sdandinski[-11]
res <- optim(runif(46, -1, 1), cost_function, gr = NULL, ctl_binded_9mer_separated, 0, 2, method = "Nelder-Mead") # this currently gives the best result

scores <- evaluate(res$par, ctl_binded_9mer_separated)

ctl_binded_9mer_separated$pred <- round(scores)
