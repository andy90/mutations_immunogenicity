# this is used to train the model, which is plan A
rm(list=ls()) # clean the memeory

ctl_A0201 <- read.csv("ctl_A0201_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele A0201
ctl_B5701 <- read.csv("ctl_B5701_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele B5701
ctl_A0301 <- read.csv("ctl_A0301_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele A0301
ctl_B0702 <- read.csv("ctl_B0702_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele B0702

# only work on those peptides that bind to the MHC
bind_criteria <- 1000 # the binding IC50 is 1000nM
ctl_A0201_binded <- ctl_A0201[ctl_A0201$ic50 < bind_criteria, ]
ctl_B5701_binded <- ctl_B5701[ctl_B5701$ic50 < bind_criteria, ]
ctl_A0301_binded <- ctl_A0301[ctl_A0301$ic50 < bind_criteria, ]
ctl_B0702_binded <- ctl_B0702[ctl_B0702$ic50 < bind_criteria, ]

ctl_binded_combinded <- do.call("rbind", list(ctl_A0201_binded, ctl_A0301_binded)) # bind all the ctls into one
ctl_binded_9mer <- ctl_binded_combinded[nchar(ctl_binded_combinded$peptide) == 9, ]
ctl_binded_9mer$Immunogenicity <- round(ctl_binded_9mer$Immunogenicity)
# the maximum length of the peptide is 12, the majority of them are of length 9 and 10
# we will do the left alignment to get the weight of each position
score_of_aa <- function(aa, aa_name, p){ # returns the p score of each amino acid
  score <- p[aa == aa_name]
}

score_of_pep <- function(peptide, aa_name, p){ # returns the p score of a peptide
  scores <- sapply(strsplit(peptide, "")[[1]], score_of_aa, aa_name, p)
}

coupling_of_pep <- function(peptide, aa_name, p){ # calculate the coupling between the peptide scores
  scores <- score_of_pep(peptide, aa_name, p)
  lscores <- length(scores)
  scores[2:lscores] * scores[1:(lscores - 1)] # this is all the coupling term between the scores
}

total_score <- function(row_of_ctl, px, aa_name){ # returns the p score of each peptide
  peptide <- row_of_ctl[1]
  l_pep <- 9 # only focus on the 9mers in this case
  p=px[1:length(aa_name)] # extract p from px
  
  x <- px[(length(aa_name) + 1) : (l_pep + 1 + length(aa_name))] # this is the coefficient for the single term, including the constant term
  x_coupling <- px[(l_pep + 2 + length(aa_name)) : length(px)] # this is the coefficient for the coupling term
  
  scores <- sapply(strsplit(peptide, "")[[1]], score_of_aa, aa_name, p) # this is the scores for each peptide
  
  coupling_scores <- scores[2:l_pep] * scores[1:(l_pep - 1)] # this is all the coupling term between the scores
  t_score <- sum(c(1,scores)*x) + sum(coupling_scores * x_coupling)
}

score_all_aa <- function(aa, peptide_aa, x){ # count how many times each amino acide appreas on each postion, take summation of x
  sum(x[which(aa == peptide_aa)])
}

cost_function <- function(px, ctl_binded_combinded, lambda_p, lambda_x){ # this calculate the total cost
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p

  t_scores <- apply(ctl_binded_combinded, 1, total_score, px, aa_name) # this is the total energy part Ec + x*p for each peptide
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
  
  y <- ctl_binded_combinded$Immunogenicity
  m <- length(y)
  p <- px[1:length(aa_name)] # extract p from px
  x <- px[(length(aa_name) + 2): length(px)] # extract x from px, exclude the constant term, include all the other coefficient terms
  
  total_cost <- -mean(y*log(sigmoid_tscores) + (1 - y)*log(1 - sigmoid_tscores)) + lambda_p/(2*m)*sum(px[1:length(aa_name)] * px[1:length(aa_name)]) + lambda_x/(2*m)*sum(x * x)
}





evaluate <- function(px, ctl_binded_combinded){ # this evaluate the trained model
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p

  t_scores <- apply(ctl_binded_combinded, 1, total_score, px, aa_name) # this is the total energy part Ec + x*p for each peptide
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
}

px <- c(rep(1,21), rep(2,18))
res <- optim(runif(39, -0.1, 0.1), cost_function, gr = NULL, ctl_binded_9mer, 0, 0, method = "Nelder-Mead") # this currently gives the best result
