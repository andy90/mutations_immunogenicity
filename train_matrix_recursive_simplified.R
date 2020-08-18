# this is a simplified version of train_matrix_recursive.R. We just focus on a single allele A0201 here. And we don't implement regularization
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

ctl_A0201_binded_9mer <- ctl_A0201_binded[nchar(ctl_A0201_binded$peptide) == 9, ]
ctl_A0201_binded_9mer$Immunogenicity <- round(ctl_A0201_binded_9mer$Immunogenicity)
ctl_A0201_binded_9mer_reshuffle <- ctl_A0201_binded_9mer[sample(nrow(ctl_A0201_binded_9mer)), ]
index_train_set <- round(nrow(ctl_A0201_binded_9mer_reshuffle)*0.7)
train_set <- ctl_A0201_binded_9mer_reshuffle[1:index_train_set, ]
test_set <- ctl_A0201_binded_9mer_reshuffle[(index_train_set+1):nrow(ctl_A0201_binded_9mer_reshuffle), ]

aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is how the map between the amino acid and the corresponding score vector p
allele_name <- c("A0201", "B5701", "A0301", "B0702") # this is the map between the  allele name the the weighit vector x
len_pep <- 9 # this is the part of peptide length that we will consider

p_wif <- c(0.81, 0.17, 0.99, 1.23, 2.02, 0.13, 0.14, 0.42, 0.58, -0.24, 0, 0.01, 0.45, 0.17, 0.07, -0.31, -0.56, -0.23, -1.13, -0.94, -1.85)
p_woct <- c(1.81, 2.33, 2.80, 3.64, 3.63, 0.46, 0.25, 0.85, 0.77, -0.02, 0, 1.15, 0.14, 0.50, -0.46, -1.12, -1.25, -0.67, -1.71, -0.71, -2.09)
p_Sdandinski <- c(-0.2, 0.0, -1, -0.6, -1.2, -0.4, -0.5, -0.05, -0.8, 1.5, 0.0, 0.0, 0.1, -0.1, 0.05, 0.0, 0.2, -0.05, 0.4, 0.8, 1.2)

peptide_to_matrix <- function(peptide, len_pep, aa_name, i_allele, n_allele){ # this function transforms a peptide sequence into a matrix which documents the amino acid type at each position of the matrix
  # peptide is a character array, len_pep is the length of peptide that we are interested in. If len_pep is 9, and the peptide
  # length is 10, we will only look at the first 9 amino acid. If the peptide length is 8, we will put a null on the 9th position
  # aa_name is a character array which documents the name of each amino accid
  # i_allele is the index for the MHC allele
  # n_allele is the total number of alleles
  n_aa <- length(aa_name)
  rdim <- n_aa + 1
  cdim <- len_pep + 1
  A <- matrix(integer(cdim * rdim * n_allele * n_allele), nrow = n_allele * rdim)
  
  A[(i_allele - 1) * rdim + 1, (i_allele - 1) * cdim + 1] <- 1
  
  if (length(peptide) < len_pep){
    peptide <- c(peptide, rep("X", len_pep - length(peptide)))
  }
  
  A_sub <- sapply(peptide[1 : len_pep], function(aa){
    i_aa <- which(aa == aa_name)
    v <- rep(0, n_aa)
    v[i_aa] <- 1
    v
  })
  
  A[((i_allele - 1) * rdim + 2) : (i_allele * rdim), ((i_allele - 1) * cdim + 2) : (i_allele * cdim)] <- A_sub
  
  A
}

costfunction_fixx <- function(p, x, ctl_binded_combinded, lambda_p, len_pep, aa_name, allele_name){
  n_allele <- length(allele_name)
  n_aa <- length(aa_name)
  u_expand <- matrix(rep(diag(n_aa + 1), n_allele), nrow = (n_aa + 1))
  
  
  t_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    peptide_aa <- strsplit(row_of_ctl[1], "")[[1]] # we disassociate the peptide string into array of characters
    A <- peptide_to_matrix(peptide_aa, len_pep, aa_name, i_allele, n_allele)
    score <- p %*% u_expand %*% A %*% x
  })
  
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
  
  y <- ctl_binded_combinded$Immunogenicity
  m <- length(y)
  
  total_cost <- -mean(y*log(sigmoid_tscores) + (1 - y)*log(1 - sigmoid_tscores)) + lambda_p/(2*m)*sum(p[2:length(p)] * p[2:length(p)])
}

d_costfunction_fixx <- function(p, x, ctl_binded_combinded, lambda_p, len_pep, aa_name, allele_name){
  n_allele <- length(allele_name)
  n_aa <- length(aa_name)
  u_expand <- matrix(rep(diag(n_aa + 1), n_allele), nrow = (n_aa + 1))
  
  t_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    peptide_aa <- strsplit(row_of_ctl[1], "")[[1]] # we disassociate the peptide string into array of characters
    A <- peptide_to_matrix(peptide_aa, len_pep, aa_name, i_allele, n_allele)
    score <- p %*% u_expand %*% A %*% x
  })
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
  
  d_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    peptide_aa <- strsplit(row_of_ctl[1], "")[[1]] # we disassociate the peptide string into array of characters
    A <- peptide_to_matrix(peptide_aa, len_pep, aa_name, i_allele, n_allele)
    dscore <- u_expand %*% A %*% x
  })
  
  y <- ctl_binded_combinded$Immunogenicity
  m <- length(y)
  dpx_noreg <- (d_scores %*% (sigmoid_tscores - y))/m
  
  dpx_reg<- c(0, rep(lambda_p/m,length(aa_name)))
  
  dpx <- dpx_noreg + dpx_reg
}

costfunction_fixp <- function(x, p, ctl_binded_combinded, lambda_x, len_pep, aa_name, allele_name){
  n_allele <- length(allele_name)
  n_aa <- length(aa_name)
  u_expand <- matrix(rep(diag(n_aa + 1), n_allele), nrow = (n_aa + 1))
  
  
  t_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    peptide_aa <- strsplit(row_of_ctl[1], "")[[1]] # we disassociate the peptide string into array of characters
    A <- peptide_to_matrix(peptide_aa, len_pep, aa_name, i_allele, n_allele)
    score <- p %*% u_expand %*% A %*% x
  })
  
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
  
  y <- ctl_binded_combinded$Immunogenicity
  m <- length(y)
  index_constx <- (1:n_allele - 1)*(len_pep + 1) + 1 # these are the Ec
  total_cost <- -mean(y*log(sigmoid_tscores) + (1 - y)*log(1 - sigmoid_tscores)) + lambda_x/(2*m)*sum(x[-index_constx] * x[-index_constx])
}

d_costfunction_fixp <- function(x, p, ctl_binded_combinded, lambda_x, len_pep, aa_name, allele_name){
  n_allele <- length(allele_name)
  n_aa <- length(aa_name)
  u_expand <- matrix(rep(diag(n_aa + 1), n_allele), nrow = (n_aa + 1))
  
  t_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    peptide_aa <- strsplit(row_of_ctl[1], "")[[1]] # we disassociate the peptide string into array of characters
    A <- peptide_to_matrix(peptide_aa, len_pep, aa_name, i_allele, n_allele)
    score <- p %*% u_expand %*% A %*% x
  })
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
  
  d_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    peptide_aa <- strsplit(row_of_ctl[1], "")[[1]] # we disassociate the peptide string into array of characters
    A <- peptide_to_matrix(peptide_aa, len_pep, aa_name, i_allele, n_allele)
    dscore <- p %*% u_expand %*% A
  })
  
  y <- ctl_binded_combinded$Immunogenicity
  m <- length(y)
  dpx_noreg <- (d_scores %*% (sigmoid_tscores - y))/m
  
  dpx_reg<- rep(c(0, rep(lambda_x/m,len_pep)), 4)
  
  dpx <- dpx_noreg + dpx_reg
}

evaluate <- function(p, x, ctl_binded_combinded, lambda_p, len_pep, aa_name, allele_name){
  n_allele <- length(allele_name)
  n_aa <- length(aa_name)
  u_expand <- matrix(rep(diag(n_aa + 1), n_allele), nrow = (n_aa + 1))
  
  
  t_scores <- apply(ctl_binded_combinded, 1, function(row_of_ctl){
    i_allele <- which(row_of_ctl[5] == allele_name) # extract x which corresponds to the specifc allele
    peptide_aa <- strsplit(row_of_ctl[1], "")[[1]] # we disassociate the peptide string into array of characters
    A <- peptide_to_matrix(peptide_aa, len_pep, aa_name, i_allele, n_allele)
    score <- p %*% u_expand %*% A %*% x
  })
  
  sigmoid_tscores <- 1/(1+exp(-t_scores)) # this is the sigmoid function based on the scores
}

x <- rep(1, length(allele_name) * (len_pep + 1)) # this is the initial x
p <- c(1, p_Sdandinski) # this is the initial p


for (i in 1:10){
  res <- optim(x, costfunction_fixp, d_costfunction_fixp, p, train_set, 0, len_pep, aa_name, allele_name, method = "BFGS") # this currently gives the best result
  x<- res$par
  print(res$value)
  res <- optim(p, costfunction_fixx, d_costfunction_fixx, x, train_set, 0, len_pep, aa_name, allele_name, method = "BFGS")
  p <- res$par
  print(res$value)
}

pred <- evaluate(p, x, train_set, 0, len_pep, aa_name, allele_name) # see what prediction is made by the trained model
exper_pred <- cbind(round(train_set$Immunogenicity), round(pred)) # this is easy for us to see which prediction is right, which prediction is wrong
length(which(exper_pred[,1] != exper_pred[,2])) # find out how many predictions are wrong
fn <- length(which((exper_pred[,1] == 1) & (exper_pred[,2] == 0))) # this is false negative
fp <- length(which((exper_pred[,1] == 0) & (exper_pred[,2] == 1))) # this is false positive
pos <- length(which(exper_pred[,1] == 1))
neg <- length(which(exper_pred[,1] == 0))

plot(factor(aa_name), -p[2:length(p)])

pred_test <- evaluate(p, x, test_set, 0, len_pep, aa_name, allele_name) # see what prediction is made by the trained model
exper_pred_test <- cbind(round(test_set$Immunogenicity), round(pred_test)) # this is easy for us to see which prediction is right, which prediction is wrong
length(which(exper_pred_test[,1] != exper_pred_test[,2])) # find out how many predictions are wrong
fn_test <- length(which((exper_pred_test[,1] == 1) & (exper_pred_test[,2] == 0))) # this is false negative
fp_test <- length(which((exper_pred_test[,1] == 0) & (exper_pred_test[,2] == 1))) # this is false positive
pos_test <- length(which(exper_pred_test[,1] == 1))
neg_test <- length(which(exper_pred_test[,1] == 0))
# without regularization, we have 20% error rate for the training set and 30% or higher for the test set, however, the results fluctuates a lot. Different reshuffle might produce pretty different results and error rate
# most mis predictions are false positives