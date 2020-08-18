# this file is used to associate single mutations with the immunogenicity
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

ctl_binded_combinded <- do.call("rbind", list(ctl_A0201_binded, ctl_A0301_binded, ctl_B0702_binded, ctl_B5701_binded))

ctl_binded_9mer <- ctl_binded_combinded[nchar(ctl_binded_combinded$peptide)==9, ]
#ctl_binded_9mer$Immunogenicity <- round(ctl_binded_9mer$Immunogenicity)

find_singlemut_immuincrease <- function(r_ctl, row_ctl_data, A){ # find out the single mutations that increase immunogenicity
    ctl_pep_array <- strsplit(as.character(row_ctl_data[1]), "")[[1]]
    pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
    index_mut <- which(ctl_pep_array != pep_array)
    if (length(index_mut) == 1) {
      if ((r_ctl[4] == 0) & (row_ctl_data[4] == 1)){
        A[pep_array[index_mut], ctl_pep_array[index_mut]] = A[pep_array[index_mut], ctl_pep_array[index_mut]] + 1
      }
    }
  A
}

find_singlemut_immuincrease_noround <- function(r_ctl, row_ctl_data, A){ # find out the single mutations that increase immunogenicity, the immunogenicity is not rounded to 0 or 1
  ctl_pep_array <- strsplit(as.character(row_ctl_data[1]), "")[[1]]
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  index_mut <- which(ctl_pep_array != pep_array)
  if (length(index_mut) == 1) {
    if ((r_ctl[4] < 0.3) & (row_ctl_data[4] > 0.7)){
      A[pep_array[index_mut], ctl_pep_array[index_mut]] = A[pep_array[index_mut], ctl_pep_array[index_mut]] + 1
    }
  }
  A
}
sum_singlemut_immuincrease <- function(ctl_binded_9mer, mut_m){
  apply(ctl_binded_9mer, 1, function(r_ctl){
    Ai <- apply(ctl_binded_9mer, 1, find_singlemut_immuincrease, r_ctl, mut_m)
    Ai <- rowSums(Ai)
  })  
}

sum_singlemut_immuincrease_noround <- function(ctl_binded_9mer, mut_m){
  apply(ctl_binded_9mer, 1, function(r_ctl){
    Ai <- apply(ctl_binded_9mer, 1, find_singlemut_immuincrease_noround, r_ctl, mut_m)
    Ai <- rowSums(Ai)
  })  
}

find_singlemut_immu_keepone <- function(r_ctl, row_ctl_data, A){ # find out the single mutations that keeps immunogenicity to be 1
  ctl_pep_array <- strsplit(as.character(row_ctl_data[1]), "")[[1]]
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  index_mut <- which(ctl_pep_array != pep_array)
  if (length(index_mut) == 1) {
    if ((r_ctl[4] == 1) & (row_ctl_data[4] == 1)){
      A[pep_array[index_mut], ctl_pep_array[index_mut]] = A[pep_array[index_mut], ctl_pep_array[index_mut]] + 1
    }
  }
  A
}

sum_singlemut_immuu_keepone <- function(ctl_binded_9mer, mut_m){
  apply(ctl_binded_9mer, 1, function(r_ctl){
    Ai <- apply(ctl_binded_9mer, 1, find_singlemut_immu_keepone, r_ctl, mut_m)
    Ai <- rowSums(Ai)
  })  
}

find_singlemut_immu_keepzero <- function(r_ctl, row_ctl_data, A){ # find out the single mutations that keeps immunogenicity to be 0
  ctl_pep_array <- strsplit(as.character(row_ctl_data[1]), "")[[1]]
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  index_mut <- which(ctl_pep_array != pep_array)
  if (length(index_mut) == 1) {
    if ((r_ctl[4] == 0) & (row_ctl_data[4] == 0)){
      A[pep_array[index_mut], ctl_pep_array[index_mut]] = A[pep_array[index_mut], ctl_pep_array[index_mut]] + 1
    }
  }
  A
}

sum_singlemut_immuu_keepzero <- function(ctl_binded_9mer, mut_m){
  apply(ctl_binded_9mer, 1, function(r_ctl){
    Ai <- apply(ctl_binded_9mer, 1, find_singlemut_immu_keepzero, r_ctl, mut_m)
    Ai <- rowSums(Ai)
  })  
}

mut_m <- matrix(integer(21*21), nrow = 21, ncol = 21) # this is a template for the mutation matrix, used as zero point to be passed to find_singlemut
aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W")
colnames(mut_m) <- aa_name
rownames(mut_m) <- aa_name



A <- sum_singlemut_immuincrease_noround(ctl_binded_9mer, mut_m)
Am <- matrix(rowSums(A), nrow = 21) # 0 to 1 # there are total of 67 0 to 1 transfer

B <- sum_singlemut_immuu_keepone(ctl_binded_9mer, mut_m)
Bm <- matrix(rowSums(B), nrow = 21) # 1 to 1 # there are total of 93(186/2) 1 to 1 transfer

C <- sum_singlemut_immuu_keepzero(ctl_binded_9mer, mut_m)
Cm <- matrix(rowSums(C), nrow = 21) # 0 to 0 # there are total of 29(58/2) 0 to 0 transfer

# we want the rations below to be as small as possible
sum( Am & Cm )/sum( Am & 1 )
sum( Am & Cm )/sum( Cm & 1 )
sum( t(Am) & Bm )/sum( t(Am) & 1 ) # currently this is the largest one, which menas that we can not really distinguish 1 to 0 and 1 to 1 transition
sum( t(Am) & Bm )/sum( Bm & 1 )

# we want to know the probability of a mutation that will transform the sequence from 0 to 1, instead of from 1 to 0
Am / (Am + t(Am) + 0.000000001) # those non-zero elements measures the probability of a mutation that transforms the sequence from 0 to 1, instead of 1 to 0. The 0.5 is the least we want -- it means that a mutation can equally take the sequence from 0 to 1 or 1 to 0. this is counter intuitive
Am + t(Am)

# we want to know the probability of a mutation that will transform the sequence from 0 to 1, instead of 1 to 0, 0 to 0 or 1 to 1
p_increase <- Am / (Am + t(Am) + Bm + Cm + 0.000000001)
# probability that a mutation will transform the sequence from 1 to 0, instead of 0 to 1 or unchange
p_decrease <- t(Am) / (Am + t(Am) + Bm + Cm + 0.000000001)
# probability that a mutation will keep the sequence immunogenicity unchanged
p_keep <- (Bm + Cm) / (Am + t(Am) + Bm + Cm + 0.000000001)
# the entropy of each item
entro <- -p_increase * log(p_increase + 0.0001) -p_decrease * log(p_decrease + 0.0001) -p_keep * log(p_keep + 0.0001) # the maximum is about 1.1, minimums is about -0.0001, all the 0s are those mutations that does not happen
# there are 80 elements of entro that has lowest entropy. it might be enough to use these elements to rank all the amino acids
# we can try to use a single parameter to describe these 3 matrix
# do the logistic regression to find the ranking of all the amino acid
costfunction <- function(p, p_increase, total_m){
  # the first element of p is delta E, which is the constant correction
  p_noconst <- p[2:22] # take out the const correction
  pi <- matrix(rep(p_noconst,21), ncol = 21, byrow = FALSE)
  pj <- matrix(rep(p_noconst,21), ncol = 21, byrow = TRUE)
  scores_ij <- 1.0/(1.0 + exp(-(pj - pi - p[1])))
  total_cost <- sum(- total_m * (p_increase * log(scores_ij) + (1- p_increase) * log(1-scores_ij))) # the cost is weighted by the number of mutations that actually happen
}

costfunction2 <- function(p, p_increase, total_m){
  # the first element of p is delta E, which is the constant correction
  p_noconst <- p[2:22] # take out the const correction
  pi <- matrix(rep(p_noconst,21), ncol = 21, byrow = FALSE)
  pj <- matrix(rep(p_noconst,21), ncol = 21, byrow = TRUE)
  scores_ij <- 1.0/(1.0 + exp(-(pj - pi - p[1])))
  total_cost <- sum(- (total_m & 1) * (p_increase * log(scores_ij) + (1- p_increase) * log(1-scores_ij))) # the cost is weighted by the number of mutations that actually happen
}

test <- function(p, total_m){
  p_noconst <- p[2:22] # take out the const correction
  pi <- matrix(rep(p_noconst,21), ncol = 21, byrow = FALSE)
  pj <- matrix(rep(p_noconst,21), ncol = 21, byrow = TRUE)
  scores_ij <- 1.0/(1.0 + exp(-(pj - pi - p[1]))) # just want to get the scores for the costfunction
  (total_m & 1) * scores_ij
}


total_m <- Am + t(Am) + Bm + Cm # this records the total number of mutations
p_Sdandinski <- c(-0.2, 0.0, -1, -0.6, -1.2, -0.4, -0.5, -0.05, -0.8, 1.5, 0.0, 0.0, 0.1, -0.1, 0.05, 0.0, 0.2, -0.05, 0.4, 0.8, 1.2)
p_wif <- c(0.81, 0.17, 0.99, 1.23, 2.02, 0.13, 0.14, 0.42, 0.58, -0.24, 0, 0.01, 0.45, 0.17, 0.07, -0.31, -0.56, -0.23, -1.13, -0.94, -1.85)
cost <- costfunction(c(1,p_wif), p_increase, total_m)
res <- optim(c(1,p_Sdandinski), costfunction2, gr = NULL, p_increase, total_m)
scores <- test(res$par, total_m) # for the scores we obtained, non get close to 1, most close to 0.
# we can associate the mutations with positions and see whether the overlap is better
find_singlemut_immuincrease_loc <- function(r_ctl, row_ctl_data, A, loc){ # find out the single mutations that increase immunogenicity at specific location
  ctl_pep_array <- strsplit(as.character(row_ctl_data[1]), "")[[1]]
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  index_mut <- which(ctl_pep_array != pep_array)
  if (length(index_mut) == 1) {
    if (index_mut == loc){
      if ((r_ctl[4] == 0) & (row_ctl_data[4] == 1)){
        A[pep_array[index_mut], ctl_pep_array[index_mut]] = A[pep_array[index_mut], ctl_pep_array[index_mut]] + 1
      }
    }
  }
  A
}

sum_singlemut_immuincrease_loc <- function(ctl_binded_9mer, mut_m, loc){
  apply(ctl_binded_9mer, 1, function(r_ctl){ # this is essentially a double for loop
    Ai <- apply(ctl_binded_9mer, 1, find_singlemut_immuincrease_loc, r_ctl, mut_m, loc)
    Ai <- rowSums(Ai)
  })  
}

find_singlemut_immu_keepone_loc <- function(r_ctl, row_ctl_data, A, loc){ # find out the single mutations that keep immunogenicity to be 1 at specific location
  ctl_pep_array <- strsplit(as.character(row_ctl_data[1]), "")[[1]]
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  index_mut <- which(ctl_pep_array != pep_array)
  if (length(index_mut) == 1) {
    if (index_mut == loc){
      if ((r_ctl[4] == 1) & (row_ctl_data[4] == 1)){
        A[pep_array[index_mut], ctl_pep_array[index_mut]] = A[pep_array[index_mut], ctl_pep_array[index_mut]] + 1
      }
    }
  }
  A
}

sum_singlemut_immu_keepone_loc <- function(ctl_binded_9mer, mut_m, loc){
  apply(ctl_binded_9mer, 1, function(r_ctl){
    Ai <- apply(ctl_binded_9mer, 1, find_singlemut_immu_keepone_loc, r_ctl, mut_m, loc)
    Ai <- rowSums(Ai)
  })  
}

find_singlemut_immu_keepzero_loc <- function(r_ctl, row_ctl_data, A, loc){ # find out the single mutations that keep immunogenicity to be 0 at specific location
  ctl_pep_array <- strsplit(as.character(row_ctl_data[1]), "")[[1]]
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  index_mut <- which(ctl_pep_array != pep_array)
  if (length(index_mut) == 1) {
    if (index_mut == loc){
      if ((r_ctl[4] == 0) & (row_ctl_data[4] == 0)){
        A[pep_array[index_mut], ctl_pep_array[index_mut]] = A[pep_array[index_mut], ctl_pep_array[index_mut]] + 1
      }
    }
  }
  A
}

sum_singlemut_immu_keepzero_loc <- function(ctl_binded_9mer, mut_m, loc){
  apply(ctl_binded_9mer, 1, function(r_ctl){
    Ai <- apply(ctl_binded_9mer, 1, find_singlemut_immu_keepzero_loc, r_ctl, mut_m, loc)
    Ai <- rowSums(Ai)
  })  
}

A1 <- sum_singlemut_immuincrease_loc(ctl_binded_9mer, mut_m, 1)
Am1 <- matrix(rowSums(A1), nrow = 21)

B1 <- sum_singlemut_immu_keepone_loc(ctl_binded_9mer, mut_m, 1)
Bm1 <- matrix(rowSums(B1), nrow = 21)

C1 <- sum_singlemut_immu_keepzero_loc(ctl_binded_9mer, mut_m, 1)
Cm1 <- matrix(rowSums(C1), nrow = 21)

A9 <- sum_singlemut_immuincrease_loc(ctl_binded_9mer, mut_m, 9)
Am9 <- matrix(rowSums(A9), nrow = 21)

B9 <- sum_singlemut_immu_keepone_loc(ctl_binded_9mer, mut_m, 9)
Bm9 <- matrix(rowSums(B9), nrow = 21)

C9 <- sum_singlemut_immu_keepzero_loc(ctl_binded_9mer, mut_m, 9)
Cm9 <- matrix(rowSums(C9), nrow = 21)

sum(t(Am9) & Bm9) / sum(Am9 & 1) # the distinguishablibity between t(Am9) and Bm9 is pretty high, we probably need more information other than positions to distinguish the mutations that will change immnogenicity
sum(t(Am9) & Bm9) / sum(Bm9 & 1)
sum(Am9 & Cm9) / sum(Am9 & 1)
sum(Am9 & Cm9) / sum(Cm9 & 1)

get_ABC_loc <- function(ctl, m, loc){
  A <- sum_singlemut_immuincrease_loc(ctl, m, loc)
  Am <- matrix(rowSums(A), nrow = 21)
  
  B <- sum_singlemut_immu_keepone_loc(ctl, m, loc)
  Bm <- matrix(rowSums(B), nrow = 21)
  
  C <- sum_singlemut_immu_keepzero_loc(ctl, m, loc)
  Cm <- matrix(rowSums(C), nrow = 21)
  
  list("Am" = Am, "Bm" = Bm, "Cm" = Cm)
}

get_ABC <- function(ctl, m){
  A <- sum_singlemut_immuincrease(ctl, m)
  Am <- matrix(rowSums(A), nrow = 21)
  
  B <- sum_singlemut_immuu_keepone(ctl, m)
  Bm <- matrix(rowSums(B), nrow = 21)
  
  C <- sum_singlemut_immuu_keepzero(ctl, m)
  Cm <- matrix(rowSums(C), nrow = 21)
  
  list("Am" = Am, "Bm" = Bm, "Cm" = Cm)
}

get_distinguish <- function(ctl, m, loc){ # see whether certain mutations are distinguishable
  res <- get_ABC_loc(ctl, m, loc)
  Am <- res$Am
  Bm <- res$Bm
  Cm <- res$Cm
  d0_01 <- sum( Am & Cm )/sum( Am & 1 ) # number of 0 to 0 and 0 to 1 indistinguish compared to total number of 0 to 1 at the same position 
  d0_00 <- sum( Am & Cm )/sum( Cm & 1 )
  d1_10 <- sum( t(Am) & Bm )/sum( t(Am) & 1 ) 
  d1_11 <- sum( t(Am) & Bm )/sum( Bm & 1 )
  d01_10 <- sum( Am & t(Am) )/sum( Am & 1) # those mutations that brought 0 to 1 and 1 to 0, they should not overlap
  
  list("d0_01" = d0_01, "d0_00" = d0_00, "d1_10" = d1_10, "d1_11" = d1_11, "d01_10" = d01_10)
}

get_m0 <- function(loc, ctl, m){ # all the mutations from 0, see wheteher is will keep the immunogenicity unchanged or not. use a number from 0 to 1 to measure it
  
  A0 <- matrix(rowSums(sum_singlemut_immuincrease_loc(ctl, m, loc)), nrow = 21)
  C0 <- matrix(rowSums(sum_singlemut_immu_keepzero_loc(ctl, m, loc)), nrow = 21)
  
  ((!C0) & A0) + (C0/(C0 + A0 + 0.000000001)) - ((!A0) & (!C0))
}

get_m1 <- function(loc, ctl, m){ # all the mutations from 1, see wheteher is will keep the immunogenicity unchanged or not. use a number from 0 to 1 to measure it
  
  A0 <- matrix(rowSums(sum_singlemut_immuincrease_loc(ctl, m, loc)), nrow = 21)
  B0 <- matrix(rowSums(sum_singlemut_immu_keepone_loc(ctl, m, loc)), nrow = 21)
  
  ((!B0) & t(A0)) + (B0/(B0 + t(A0) + 0.000000001)) - ((!t(A0)) & (!B0))
}

get_n0 <- function(loc, ctl, m){ # count the number of mutations from 0
  A0 <- matrix(rowSums(sum_singlemut_immuincrease_loc(ctl, m, loc)), nrow = 21)
  C0 <- matrix(rowSums(sum_singlemut_immu_keepzero_loc(ctl, m, loc)), nrow = 21)
  
  A0 + C0
}

get_n1 <- function(loc, ctl, m){ # count the number of mutations from 1
  A0 <- matrix(rowSums(sum_singlemut_immuincrease_loc(ctl, m, loc)), nrow = 21)
  B0 <- matrix(rowSums(sum_singlemut_immu_keepone_loc(ctl, m, loc)), nrow = 21)
  
  t(A0) + B0
}


m0s <- sapply(1:9, get_m0, ctl_binded_9mer, mut_m)
n0s <- sapply(1:9, get_n0, ctl_binded_9mer, mut_m)
m0_av <- matrix(rowSums((m0s * n0s))/rowSums(n0s + 0.000000001) - (!rowSums(n0s)), nrow = 21)

m1s <- sapply(1:9, get_m1, ctl_binded_9mer, mut_m)
n1s <- sapply(1:9, get_n1, ctl_binded_9mer, mut_m)
m1_av <- matrix(rowSums((m1s * n1s))/rowSums(n1s + 0.000000001) - (!rowSums(n1s)), nrow = 21)

n0_tot <- matrix(rowSums(n0s), nrow = 21)
n1_tot <- matrix(rowSums(n1s), nrow = 21)
