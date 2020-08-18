# in this file, we do the recursive training of x and p, still not sure whether the underlying cost function is convex of not
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

ctl_binded_combinded$Immunogenicity <- round(ctl_binded_combinded$Immunogenicity) # round the Immunogenicity to 0 or 1, 0.5 is round to 0

ctl_binded_combinded <- ctl_binded_combinded[nchar(ctl_binded_combinded$peptide)==9, ] # only look at the 9 mers for now, which accounts for 3/4 of the data

cost_function_fixp <- function(x_input, ctl_binded_combinded, lambda_x, p_input){ # this calculate the total cost with p fixed
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
  allele_name <- c("A0201", "B5701", "A0301", "B0702") # this is the map between the  allele name the the weighit vector x
  
  px <- c(p_input, x_input) # re assemble the px based on input from x and p, p_input is now a parameter
  lambda_p <- 0 # set lambda_p = 0 in this case, since we are not gonna vary it
  
  t_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){ # returns the p score of each peptide
    peptide <- row_of_ctl[1]
    
    p=px[1:length(aa_name)] # extract p from px
    
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    x <- px[((i_allele-1)*13 + length(aa_name) + 1): (i_allele*13 + length(aa_name))]
    
    scores <- sapply(strsplit(peptide, "")[[1]], function(aa){
      score <- p[aa == aa_name]
    })
    lscores <- length(scores)
    t_score <- sum(c(1,scores)*x[1:(lscores+1)])
  }) # this is the total energy part Ec + x*p for each peptide
  
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
  
  y <- ctl_binded_combinded$Immunogenicity
  m <- length(y)
  p <- px[1:length(aa_name)] # extract p from px
  x <- px[(length(aa_name) + 1): length(px)] # extract x from px
  index_constx <- (1:length(allele_name) - 1)*13 + 1 # these are the Ec
  total_cost <- -mean(y*log(sigmoid_tscores) + (1 - y)*log(1 - sigmoid_tscores)) + lambda_p/(2*m)*sum(px[1:length(aa_name)] * px[1:length(aa_name)]) + lambda_x/(2*m)*sum(x[-index_constx] * x[-index_constx])
}

dtotal_score_fixp <- function(x_input, ctl_binded_combinded, lambda_x, p_input){ # this calculate the total derivative with p fixed
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
  allele_name <- c("A0201", "B5701", "A0301", "B0702") # this is the map between the  allele name the the weighit vector x
  px <- c(p_input, x_input) # re assemble the px based on input from x and p, p_input is now a parameter
  lambda_p <- 0 # set lambda_p = 0 in this case, since we are not gonna vary it
  t_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){
    peptide <- row_of_ctl[1]
    
    p=px[1:length(aa_name)] # extract p from px
    
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    x <- px[((i_allele-1)*13 + length(aa_name) + 1): (i_allele*13 + length(aa_name))]
    
    scores <- sapply(strsplit(peptide, "")[[1]], function(aa){
      score <- p[aa == aa_name]
    })
    lscores <- length(scores)
    t_score <- sum(c(1,scores)*x[1:(lscores+1)])
  }) # this is the total energy part Ec + x*p
  
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
  # need some manipulations to proceed.
  y <- ctl_binded_combinded$Immunogenicity
  m <- length(y)
  d_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){ # this returns the derivative of E with respect to p and x
    dpx <- rep(0, length(px)) # derivative of E with respect to p and x
    p <- px[1:length(aa_name)] # extract p from px
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    x_nnull <- px[((i_allele-1)*13 + length(aa_name) + 2): (i_allele*13 + length(aa_name))] # this does not include the const x
    
    peptide_aa <- strsplit(row_of_ctl[1], "")[[1]] # we disassociate the peptide string into array of characters
    
    scores <- sapply(peptide_aa, function(aa){ # returns the p score of each amino acid
      score <- p[aa == aa_name]
    })
    
    dpx[((i_allele-1)*13 + length(aa_name) + 2): ((i_allele-1)*13 + length(aa_name) + length(scores) + 1)] <- scores # dE/dx is dp, but it depends on the length of the peptides or scores
    
    dpx[(i_allele-1)*13 + length(aa_name) + 1] <- 1 # dE/dE_c is always 1
    
    weight_aa <- sapply(aa_name, function(aa){ # count how many times each amino acide appreas on each postion, take summation of x
      sum(x_nnull[which(aa == peptide_aa)])
    }) # for each amino acide in the aa_name matrix calculate this weight
    
    dpx[1:length(aa_name)] <- weight_aa # assign weight_aa to the p part of the dpx
    
    dpx # return dpx
  })
  dpx_noreg <- (d_scores %*% (sigmoid_tscores - y))/m
  
  # get the derivative part with respect to the regularization part
  dpx_reg<- c(rep(lambda_p/m,length(aa_name)), rep(lambda_x/m, length(px)-length(aa_name)))
  index_constx <- (1:length(allele_name) - 1)*13 + length(aa_name) + 1 # these are the index of Ec
  dpx_reg[index_constx] <- 0 # the constants are not subject to regularization
  dpx_reg <- dpx_reg*px
  
  dpx <- dpx_noreg + dpx_reg
  dpx[(length(aa_name)+1):length(dpx)] # only return the x part
}

cost_function_fixx <- function(p_input, ctl_binded_combinded, lambda_p, x_input){ # this calculate the total cost with x fixed
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
  allele_name <- c("A0201", "B5701", "A0301", "B0702") # this is the map between the  allele name the the weighit vector x
  
  px <- c(p_input, x_input) # re assemble the px based on input from x and p, p_input is now a parameter
  lambda_x <- 0 # set lambda_p = 0 in this case, since we are not gonna vary it
  
  t_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){ # returns the p score of each peptide
    peptide <- row_of_ctl[1]
    
    p=px[1:length(aa_name)] # extract p from px
    
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    x <- px[((i_allele-1)*13 + length(aa_name) + 1): (i_allele*13 + length(aa_name))]
    
    scores <- sapply(strsplit(peptide, "")[[1]], function(aa){
      score <- p[aa == aa_name]
    })
    lscores <- length(scores)
    t_score <- sum(c(1,scores)*x[1:(lscores+1)])
  }) # this is the total energy part Ec + x*p for each peptide
  
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
  
  y <- ctl_binded_combinded$Immunogenicity
  m <- length(y)
  p <- px[1:length(aa_name)] # extract p from px
  x <- px[(length(aa_name) + 1): length(px)] # extract x from px
  index_constx <- (1:length(allele_name) - 1)*13 + length(aa_name) + 1 # these are the Ec
  total_cost <- -mean(y*log(sigmoid_tscores) + (1 - y)*log(1 - sigmoid_tscores)) + lambda_p/(2*m)*sum(px[1:length(aa_name)] * px[1:length(aa_name)]) + lambda_x/(2*m)*sum(x[-index_constx] * x[-index_constx])
}

dtotal_score_fixx <- function(p_input, ctl_binded_combinded, lambda_p, x_input){ # this calculate the total derivative with x fixed
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
  allele_name <- c("A0201", "B5701", "A0301", "B0702") # this is the map between the  allele name the the weighit vector x
  px <- c(p_input, x_input) # re assemble the px based on input from x and p, p_input is now a parameter
  lambda_x <- 0 # set lambda_p = 0 in this case, since we are not gonna vary it
  t_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){
    peptide <- row_of_ctl[1]
    
    p=px[1:length(aa_name)] # extract p from px
    
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    x <- px[((i_allele-1)*13 + length(aa_name) + 1): (i_allele*13 + length(aa_name))]
    
    scores <- sapply(strsplit(peptide, "")[[1]], function(aa){
      score <- p[aa == aa_name]
    })
    lscores <- length(scores)
    t_score <- sum(c(1,scores)*x[1:(lscores+1)])
  }) # this is the total energy part Ec + x*p
  
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
  # need some manipulations to proceed.
  y <- ctl_binded_combinded$Immunogenicity
  m <- length(y)
  d_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){ # this returns the derivative of E with respect to p and x
    dpx <- rep(0, length(px)) # derivative of E with respect to p and x
    p <- px[1:length(aa_name)] # extract p from px
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    x_nnull <- px[((i_allele-1)*13 + length(aa_name) + 2): (i_allele*13 + length(aa_name))] # this does not include the const x
    
    peptide_aa <- strsplit(row_of_ctl[1], "")[[1]] # we disassociate the peptide string into array of characters
    
    scores <- sapply(peptide_aa, function(aa){ # returns the p score of each amino acid
      score <- p[aa == aa_name]
    })
    
    dpx[((i_allele-1)*13 + length(aa_name) + 2): ((i_allele-1)*13 + length(aa_name) + length(scores) + 1)] <- scores # dE/dx is dp, but it depends on the length of the peptides or scores
    
    dpx[(i_allele-1)*13 + length(aa_name) + 1] <- 1 # dE/dE_c is always 1
    
    weight_aa <- sapply(aa_name, function(aa){ # count how many times each amino acide appreas on each postion, take summation of x
      sum(x_nnull[which(aa == peptide_aa)])
    }) # for each amino acide in the aa_name matrix calculate this weight
    
    dpx[1:length(aa_name)] <- weight_aa # assign weight_aa to the p part of the dpx
    
    dpx # return dpx
  })
  dpx_noreg <- (d_scores %*% (sigmoid_tscores - y))/m
  
  # get the derivative part with respect to the regularization part
  dpx_reg<- c(rep(lambda_p/m,length(aa_name)), rep(lambda_x/m, length(px)-length(aa_name)))
  index_constx <- (1:length(allele_name) - 1)*13 + length(aa_name) + 1 # these are the index of Ec
  dpx_reg[index_constx] <- 0 # the constants are not subject to regularization
  dpx_reg <- dpx_reg*px
  
  dpx <- dpx_noreg + dpx_reg
  dpx[1: length(aa_name)] # only return the x part
}
evaluate <- function(px, ctl_binded_combinded){ # this evaluate the trained model
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
  allele_name <- c("A0201", "B5701", "A0301", "B0702") # this is the map between the  allele name the the weighit vector x
  
  t_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){ # returns the p score of each peptide
    peptide <- row_of_ctl[1]
    
    p=px[1:length(aa_name)] # extract p from px
    
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    x <- px[((i_allele-1)*13 + length(aa_name) + 1): (i_allele*13 + length(aa_name))]
    
    scores <- sapply(strsplit(peptide, "")[[1]], function(aa){ # returns the p score of each amino acid
      score <- p[aa == aa_name]
    })
    lscores <- length(scores)
    t_score <- sum(c(1,scores)*x[1:(lscores+1)])
  }) # this is the total energy part Ec + x*p for each peptide
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
}

# Now we will start with a hydrophobic scale
p_wif <- c(0.81, 0.17, 0.99, 0.81, 2.02, 0.13, 0.14, 0.42, 0.58, -0.24, 0, 0.01, 0.45, 0.17, 0.07, -0.31, -0.56, -0.23, -1.13, -0.94, -1.85)
p_names <- c("R+", "H", "K+", "D-", "E-", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W")
plot(factor(p_names), -p_wif)

res <- optim(runif(73-21, -0.8, 0.8), cost_function_fixp, dtotal_score_fixp, ctl_binded_combinded, 0, p_wif, method = "BFGS") # this currently gives the best result
# the result is very irregular