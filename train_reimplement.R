# this file re-implements train.R in a more structured way
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

cost_function <- function(px, ctl_binded_combinded, lambda_p, lambda_x){ # this calculate the total cost
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
  allele_name <- c("A0201", "B5701", "A0301", "B0702") # this is the map between the  allele name the the weighit vector x
  
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

dtotal_score <- function(px, ctl_binded_combinded, lambda_p, lambda_x){ # this calculate the total derivative
  aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
  allele_name <- c("A0201", "B5701", "A0301", "B0702") # this is the map between the  allele name the the weighit vector x
  
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

# px <- c(rep(1,21), rep(2,13), rep(3,13), rep(4,13), rep(5,13))#cost function will not return the correct value with this intialization
# px <- runif(73, -0.1, 0.1) # this initialization has no problem
# res <- cost_function(px, ctl_binded_combinded[1:400, ], 5, 1)
# dres <- dtotal_score(px, ctl_binded_combinded[1:400, ], 5, 1)
# 
# dres3 <- sapply(1:73, function(i){
#   move <- rep(0,73)
#   move[i] <- 0.001
#   dres2 <- (cost_function(px+move, ctl_binded_combinded[1:400, ], 5, 1) -cost_function(px, ctl_binded_combinded[1:400, ], 5, 1))/0.001
# }) # by comparing dres and dres3, we can see whether our implementation of the derivative of costfunction is correct

px_lower <- c(rep(-1,21), rep(c(-100, rep(0, 12)), 4))
px_upper <- c(rep(1,21), rep(c(100, rep(100, 12)), 4))
res <- optim(runif(73, -0.5, 0.5), cost_function, dtotal_score, ctl_binded_combinded, 0, 0, lower = px_lower, upper = px_upper, method = "L-BFGS-B") # this currently gives the best result
# different initialization always gives different result

#results <- data.frame(XP=as.numeric(), Cost=double())
results <- list()
for (i in 1:20){
  res <- optim(runif(73, -0.5, 0.5), cost_function, dtotal_score, ctl_binded_combinded[1:150, ], 0, 0, method = "BFGS") # this currently gives the best result
  # different initialization always gives different result
  #new_res <- data.frame(XP=res$par, Cost=res$value)
  #results <- rbind(results, new_res)
  #results[nrow(results) + 1, ] <- list(res$par, res$value)
  results[[i]] <- c(res$par,res$value)
}
results <- do.call(rbind, results)

results <- list()
for (i in 1:20){
  res <- optim(c(runif(21, -0.5, 0.5), rep(c(runif(1, -10.0, 10), runif(12, 0, 1)), 4)), cost_function, dtotal_score, ctl_binded_combinded, 0, 0, lower = px_lower, upper = px_upper, method = "L-BFGS-B") # this currently gives the best result
  # different initialization always gives different result
  #new_res <- data.frame(XP=res$par, Cost=res$value)
  #results <- rbind(results, new_res)
  #results[nrow(results) + 1, ] <- list(res$par, res$value)
  results[[i]] <- c(res$par,res$value)
}
results <- do.call(rbind, results)
p_scaled <- results[which.min(results[, 74]), 1:21]/results[which.min(results[, 74]), 21]
aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W")
plot(factor(aa_name), p_scaled)

pred <- evaluate(results[which.min(results[, 74]), 1:73], ctl_binded_combinded) # see what prediction is made by the trained model
exper_pred <- cbind(ctl_binded_combinded$Immunogenicity, round(pred)) # this is easy for us to see which prediction is right, which prediction is wrong
length(which(exper_pred[,1] != exper_pred[,2])) # find out how many predictions are wrong
fn <- length(which((exper_pred[,1] == 1) & (exper_pred[,2] == 0))) # this is false negative
fp <- length(which((exper_pred[,1] == 0) & (exper_pred[,2] == 1))) # this is false positive
pos <- length(which(exper_pred[,1] == 1))
neg <- length(which(exper_pred[,1] == 0))

fp_rate <- fp/neg # this is the significance level of the test
fn_rate <- fn/pos #  this is the power of the test

