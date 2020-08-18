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

ctl_binded_combinded <- do.call("rbind", list(ctl_A0201_binded, ctl_B5701_binded, ctl_A0301_binded, ctl_B0702_binded)) # bind all the ctls into one

# the maximum length of the peptide is 12, the majority of them are of length 9 and 10
# we will do the left alignment to get the weight of each position
score_of_aa <- function(aa, aa_name, p){ # returns the p score of each amino acid
  score <- p[aa == aa_name]
}

score_of_pep <- function(peptide, aa_name, p){ 
  scores <- sapply(strsplit(peptide, "")[[1]], score_of_aa, aa_name, p)
}

total_score <- function(row_of_ctl, px, aa_name, allele_name){ # returns the p score of each peptide
  peptide <- row_of_ctl[1]
  
  p=px[1:length(aa_name)] # extract p from px
  
  i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
  x <- px[((i_allele-1)*13 + length(aa_name) + 1): (i_allele*13 + length(aa_name))]
  
  scores <- sapply(strsplit(peptide, "")[[1]], score_of_aa, aa_name, p)
  lscores <- length(scores)
  t_score <- sum(c(1,scores)*x[1:(lscores+1)])
}

score_all_aa <- function(aa, peptide_aa, x){ # count how many times each amino acide appreas on each postion, take summation of x
  sum(x[which(aa == peptide_aa)])
}

cost_function <- function(px, ctl_binded_combinded, lambda_p, lambda_x){ # this calculate the total cost
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
  allele_name <- c("A0201", "B5701", "A0301", "B0702") # this is the map between the  allele name the the weighit vector x
  
  t_scores <- apply(ctl_binded_combinded, 1, total_score, px, aa_name, allele_name) # this is the total energy part Ec + x*p for each peptide
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
  
  y <- ctl_binded_combinded$Immunogenicity
  m <- length(y)
  p <- px[1:length(aa_name)] # extract p from px
  x <- px[(length(aa_name) + 1): length(px)] # extract x from px
  index_constx <- (1:length(allele_name) - 1)*13  + 1 # these are the Ec
  total_cost <- -mean(y*log(sigmoid_tscores) + (1 - y)*log(1 - sigmoid_tscores)) + lambda_p/(2*m)*sum(px[1:length(aa_name)] * px[1:length(aa_name)]) + lambda_x/(2*m)*sum(x[-index_constx] * x[-index_constx])
}



dscore_of_pep <- function(row_of_ctl, aa_name, allele_name, px){ # this is the dE for a specific peptide
  dpx <- rep(0, length(px)) # derivative of E with respect to p and x
  p <- px[1:length(aa_name)] # extract p from px
  i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
  x_nnull <- px[((i_allele-1)*13 + length(aa_name) + 2): (i_allele*13 + length(aa_name))] # this does not include the const x
  
  peptide_aa <- strsplit(row_of_ctl[1], "")[[1]] # we disassociate the peptide string into array of characters
  
  scores <- sapply(peptide_aa, score_of_aa, aa_name, p)
  dpx[((i_allele-1)*13 + length(aa_name) + 2): ((i_allele-1)*13 + length(aa_name) + length(scores) + 1)] <- scores # dE/dx is dp, but it depends on the length of the peptides or scores
  
  dpx[(i_allele-1)*13 + length(aa_name) + 1] <- 1 # dE/dE_c is always 1
  
  weight_aa <- sapply(aa_name, score_all_aa, peptide_aa, x_nnull) # for each amino acide in the aa_name matrix calculate this weight
                      
  dpx[1:length(aa_name)] <- weight_aa # assign weight_aa to the p part of the dpx
  
  dpx # return dpx
  
}


dtotal_score <- function(px, ctl_binded_combinded, lambda_p, lambda_x){ # this calculate the total derivative
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
  allele_name <- c("A0201", "B5701", "A0301", "B0702") # this is the map between the  allele name the the weighit vector x
  
  t_scores <- apply(ctl_binded_combinded, 1, total_score, px, aa_name, allele_name) # this is the total energy part Ec + x*p
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
  # need some manipulations to proceed.
  y <- ctl_binded_combinded$Immunogenicity
  m <- length(y)
  d_scores <- apply(ctl_binded_combinded, 1, dscore_of_pep, aa_name, allele_name, px)
  dpx_noreg <- (d_scores %*% (sigmoid_tscores - y))/m
  
  # get the derivative part with respect to the regularization part
  dpx_reg<- c(rep(lambda_p/m,length(aa_name)), rep(lambda_x/m, length(px)-length(aa_name)))
  index_constx <- (1:length(allele_name) - 1)*13 + length(aa_name) + 1 # these are the index of Ec
  dpx_reg[index_constx] <- 0 # the constants are not subject to regularizatio
  dpx_reg <- dpx_reg*px
  
  dpx <- dpx_noreg + dpx_reg
}

evaluate <- function(px, ctl_binded_combinded){ # this evaluate the trained model
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
  allele_name <- c("A0201", "B5701", "A0301", "B0702") # this is the map between the  allele name the the weighit vector x
  
  t_scores <- apply(ctl_binded_combinded, 1, total_score, px, aa_name, allele_name) # this is the total energy part Ec + x*p for each peptide
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
}

px <- c(rep(1,21), rep(2,13), rep(3,13), rep(4,13), rep(5,13))
px_lower <- c(rep(-1,21), rep(-100, 4*13))
px_upper <- c(rep(1,21), rep(100, 4*13))
#res <- optim(0*px, cost_function, dtotal_score, ctl_binded_combinded, 0, 1, method = "L-BFGS-B", lower = px_lower, upper = px_upper)
res <- optim(runif(73, -0.1, 0.1), cost_function, gr = NULL, ctl_binded_combinded, 0, 0, method = "BFGS") # this currently gives the best result
