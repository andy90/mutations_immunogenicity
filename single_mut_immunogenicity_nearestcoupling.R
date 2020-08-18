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
ctl_binded_9mer$Immunogenicity <- round(ctl_binded_9mer$Immunogenicity)

find_singlemut_immuincrease <- function(r_ctl, row_ctl_data, A, aa_name){ # find out the single mutations that increase immunogenicity
  ctl_pep_array <- strsplit(as.character(row_ctl_data[1]), "")[[1]]
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  
  ctl_pep_array_num <- sapply(ctl_pep_array, function(aa){ # convert the array to the number
    which(aa == aa_name)
  })
  
  pep_array_num <- sapply(pep_array, function(aa){ # convert the array to the number
    which(aa == aa_name)
  })
  
  index_mut <- which(ctl_pep_array != pep_array)
  if (length(index_mut) == 1) {
    if ((r_ctl[4] == 0) & (row_ctl_data[4] == 1)){
      if (index_mut == 1){
        A[length(aa_name) + 1, pep_array_num[2], pep_array_num[1], ctl_pep_array_num[1]] = A[length(aa_name)+1, pep_array_num[2], pep_array_num[1], ctl_pep_array_num[1]] + 1
      }else if(index_mut == 9){
        A[pep_array_num[8], length(aa_name) + 1, pep_array_num[9], ctl_pep_array_num[9]] = A[pep_array_num[8], length(aa_name) + 1, pep_array_num[9], ctl_pep_array_num[9]] + 1
      }else{
        A[pep_array_num[index_mut - 1], pep_array_num[index_mut + 1], pep_array_num[index_mut], ctl_pep_array_num[index_mut]] = A[pep_array_num[index_mut - 1], pep_array_num[index_mut + 1], pep_array_num[index_mut], ctl_pep_array_num[index_mut]] + 1
      }
    }
  }
  A
}

sum_singlemut_immuincrease <- function(ctl_binded_9mer, mut_m, aa_name){
  apply(ctl_binded_9mer, 1, function(r_ctl){
    Ai <- apply(ctl_binded_9mer, 1, find_singlemut_immuincrease, r_ctl, mut_m, aa_name)
    Ai <- rowSums(Ai)
  })  
}

find_singlemut_immukeepone <- function(r_ctl, row_ctl_data, A, aa_name){ # find out the single mutations that increase immunogenicity
  ctl_pep_array <- strsplit(as.character(row_ctl_data[1]), "")[[1]]
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  
  ctl_pep_array_num <- sapply(ctl_pep_array, function(aa){ # convert the array to the number
    which(aa == aa_name)
  })
  
  pep_array_num <- sapply(pep_array, function(aa){ # convert the array to the number
    which(aa == aa_name)
  })
  
  index_mut <- which(ctl_pep_array != pep_array)
  if (length(index_mut) == 1) {
    if ((r_ctl[4] == 1) & (row_ctl_data[4] == 1)){
      if (index_mut == 1){
        A[length(aa_name) + 1, pep_array_num[2], pep_array_num[1], ctl_pep_array_num[1]] = A[length(aa_name)+1, pep_array_num[2], pep_array_num[1], ctl_pep_array_num[1]] + 1
      }else if(index_mut == 9){
        A[pep_array_num[8], length(aa_name) + 1, pep_array_num[9], ctl_pep_array_num[9]] = A[pep_array_num[8], length(aa_name) + 1, pep_array_num[9], ctl_pep_array_num[9]] + 1
      }else{
        A[pep_array_num[index_mut - 1], pep_array_num[index_mut + 1], pep_array_num[index_mut], ctl_pep_array_num[index_mut]] = A[pep_array_num[index_mut - 1], pep_array_num[index_mut + 1], pep_array_num[index_mut], ctl_pep_array_num[index_mut]] + 1
      }
    }
  }
  A
}

sum_singlemut_immukeepone <- function(ctl_binded_9mer, mut_m, aa_name){
  apply(ctl_binded_9mer, 1, function(r_ctl){
    Ai <- apply(ctl_binded_9mer, 1, find_singlemut_immukeepone, r_ctl, mut_m, aa_name)
    Ai <- rowSums(Ai)
  })  
}

find_singlemut_immukeepzero <- function(r_ctl, row_ctl_data, A, aa_name){ # find out the single mutations that increase immunogenicity
  ctl_pep_array <- strsplit(as.character(row_ctl_data[1]), "")[[1]]
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  
  ctl_pep_array_num <- sapply(ctl_pep_array, function(aa){ # convert the array to the number
    which(aa == aa_name)
  })
  
  pep_array_num <- sapply(pep_array, function(aa){ # convert the array to the number
    which(aa == aa_name)
  })
  
  index_mut <- which(ctl_pep_array != pep_array)
  if (length(index_mut) == 1) {
    if ((r_ctl[4] == 0) & (row_ctl_data[4] == 0)){
      if (index_mut == 1){
        A[length(aa_name) + 1, pep_array_num[2], pep_array_num[1], ctl_pep_array_num[1]] = A[length(aa_name)+1, pep_array_num[2], pep_array_num[1], ctl_pep_array_num[1]] + 1
      }else if(index_mut == 9){
        A[pep_array_num[8], length(aa_name) + 1, pep_array_num[9], ctl_pep_array_num[9]] = A[pep_array_num[8], length(aa_name) + 1, pep_array_num[9], ctl_pep_array_num[9]] + 1
      }else{
        A[pep_array_num[index_mut - 1], pep_array_num[index_mut + 1], pep_array_num[index_mut], ctl_pep_array_num[index_mut]] = A[pep_array_num[index_mut - 1], pep_array_num[index_mut + 1], pep_array_num[index_mut], ctl_pep_array_num[index_mut]] + 1
      }
    }
  }
  A
}

sum_singlemut_immukeepzero <- function(ctl_binded_9mer, mut_m, aa_name){
  apply(ctl_binded_9mer, 1, function(r_ctl){
    Ai <- apply(ctl_binded_9mer, 1, find_singlemut_immukeepzero, r_ctl, mut_m, aa_name)
    Ai <- rowSums(Ai)
  })  
}

find_singlemut_immudecrease <- function(r_ctl, row_ctl_data, A, aa_name){ # find out the single mutations that increase immunogenicity
  ctl_pep_array <- strsplit(as.character(row_ctl_data[1]), "")[[1]]
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  
  ctl_pep_array_num <- sapply(ctl_pep_array, function(aa){ # convert the array to the number
    which(aa == aa_name)
  })
  
  pep_array_num <- sapply(pep_array, function(aa){ # convert the array to the number
    which(aa == aa_name)
  })
  
  index_mut <- which(ctl_pep_array != pep_array)
  if (length(index_mut) == 1) {
    if ((r_ctl[4] == 1) & (row_ctl_data[4] == 0)){
      if (index_mut == 1){
        A[length(aa_name) + 1, pep_array_num[2], pep_array_num[1], ctl_pep_array_num[1]] = A[length(aa_name)+1, pep_array_num[2], pep_array_num[1], ctl_pep_array_num[1]] + 1
      }else if(index_mut == 9){
        A[pep_array_num[8], length(aa_name) + 1, pep_array_num[9], ctl_pep_array_num[9]] = A[pep_array_num[8], length(aa_name) + 1, pep_array_num[9], ctl_pep_array_num[9]] + 1
      }else{
        A[pep_array_num[index_mut - 1], pep_array_num[index_mut + 1], pep_array_num[index_mut], ctl_pep_array_num[index_mut]] = A[pep_array_num[index_mut - 1], pep_array_num[index_mut + 1], pep_array_num[index_mut], ctl_pep_array_num[index_mut]] + 1
      }
    }
  }
  A
}

sum_singlemut_immudecrease <- function(ctl_binded_9mer, mut_m, aa_name){
  apply(ctl_binded_9mer, 1, function(r_ctl){
    Ai <- apply(ctl_binded_9mer, 1, find_singlemut_immudecrease, r_ctl, mut_m, aa_name)
    Ai <- rowSums(Ai)
  })  
}

aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W")
mut_m <- integer((length(aa_name) + 1) * (length(aa_name) + 1) * length(aa_name) * length(aa_name))
dim(mut_m) <- c((length(aa_name) + 1), (length(aa_name) + 1), length(aa_name), length(aa_name))
#A <- find_singlemut_immuincrease(ctl_binded_9mer[44, ], ctl_binded_9mer[43, ], mut_m, aa_name)
A <- sum_singlemut_immuincrease(ctl_binded_9mer, mut_m, aa_name)
A <- rowSums(A)
dim(A) <- c((length(aa_name) + 1), (length(aa_name) + 1), length(aa_name), length(aa_name))

B <- rowSums(sum_singlemut_immukeepone(ctl_binded_9mer, mut_m, aa_name))
dim(B) <- c((length(aa_name) + 1), (length(aa_name) + 1), length(aa_name), length(aa_name))

C <- rowSums(sum_singlemut_immukeepzero(ctl_binded_9mer, mut_m, aa_name))
dim(C) <- c((length(aa_name) + 1), (length(aa_name) + 1), length(aa_name), length(aa_name))

At <- rowSums(sum_singlemut_immudecrease(ctl_binded_9mer, mut_m, aa_name))
dim(At) <- c((length(aa_name) + 1), (length(aa_name) + 1), length(aa_name), length(aa_name))

sum(A)
sum(B)
sum(C)
sum(At)

sum(A & At) # this measures the overlap between 0 -> 1 and 1 -> 0 transition
sum(A & B) # this measures the overlap between 0 -> 1 and 1 -> 1 transition
sum(A & C) # this measures the overlap between 0 -> 1 and 0 -> 0 transition
# these are all supposed to be 0, and they are pretty low in this case. Remember that we mixed all the allels
# now the problem is how to do the training

# extract out all the mutations, that increase or keep, make sure they don't overlap, then try to match
extract_1d_4d <- function(indexes, l_aa){
  sapply(indexes, function(N){
    N <- N-1
    l <- as.integer(N / ((l_aa + 1) * (l_aa + 1) * l_aa))
    l_left <- N %% ((l_aa + 1) * (l_aa + 1) * l_aa)
    k <- as.integer(l_left / ((l_aa + 1) * (l_aa + 1)))
    k_left <- l_left %% ((l_aa + 1) * (l_aa + 1))
    j <- as.integer(k_left / (l_aa + 1))
    i <- k_left %% (l_aa + 1)
    c(i+1, j+1, k+1, l+1)
  })
}
i_increase <- t(extract_1d_4d(which(A >= 1), length(aa_name)))
i_nonincrease <- t(extract_1d_4d(which((A == 0) & ((A + At + B + C) >=1 )), length(aa_name)))

cost_function <- function(pJ, index_increase, index_nonincrease, l_aa){
  p <- pJ[1:l_aa]
  J <- matrix(pJ[(l_aa + 1) : (l_aa * l_aa + l_aa)], nrow = l_aa)
  E <- pJ[l_aa * l_aa + l_aa + 1]
  scores_increase <- apply(index_increase, 1, function(ind){
    dE <- p[ind[4]] - p[ind[3]]
    if (ind[1] == (l_aa + 1)){
      dE <- dE + J[ind[4], ind[2]] - J[ind[3], ind[2]]
    }else if(ind[2] == (l_aa + 1)){
      dE <- dE + J[ind[4], ind[1]] - J[ind[3], ind[1]]
    }else{
      dE <- dE + J[ind[4], ind[2]] - J[ind[3], ind[2]] + J[ind[4], ind[1]] - J[ind[3], ind[1]] 
    }
    dE
  })
  
  scores_nonincrease <- apply(index_nonincrease, 1, function(ind){
    dE <- p[ind[4]] - p[ind[3]]
    if (ind[1] == (l_aa + 1)){
      dE <- dE + J[ind[4], ind[2]] - J[ind[3], ind[2]]
    }else if(ind[2] == (l_aa + 1)){
      dE <- dE + J[ind[4], ind[1]] - J[ind[3], ind[1]]
    }else{
      dE <- dE + J[ind[4], ind[2]] - J[ind[3], ind[2]] + J[ind[4], ind[1]] - J[ind[3], ind[1]] 
    }
    dE
  })
  
  scores_increase_sigmoid <- 1.0/(1.0 + exp(-(scores_increase - E)))
  scores_nonincrease_sigmoid <- 1.0/(1.0 + exp(-(scores_nonincrease - E)))
  cost <- -sum(log(scores_increase_sigmoid)) - sum(log(1-scores_nonincrease_sigmoid))
  #cost <- length(which((scores_increase - E) < 0)) + length(which((scores_nonincrease - E) > 0))
}

pJ <- rnorm(length(aa_name) * length(aa_name) + length(aa_name) + 1)
optim(pJ, cost_function, gr = NULL, i_increase, i_nonincrease, length(aa_name))
i