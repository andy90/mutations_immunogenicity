# use monte carlo algorithm (I just mean random algorithm) to find out best ways to coarse grain the amino acids
rm(list=ls()) # clean the memeory

partition_mc <- function(aas){ # read in an array of amino acids, randomly partition them into two sets
  N <- length(aas)
  k <- sample(1:(N-1), 1) # we don't accept empty set, the length of aas has to be larger than 1
  set1 <- sample(aas, k)
  set2 <- setdiff(aas, set1)
  list(set1, set2)
}

aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W")
sets <- partition_mc(aa_name)
subsets <- partition_mc(sets[[1]])

coarse_grain_mc <- function(aas, k){ # k is the number of subsets we want
  sets <- aas
  N <- length(sets)
  subsets <- partition_mc(sets)
  while((length(subsets[[1]]) == 1) | (length(subsets[[1]]) == (N - 1))){
    subsets <- partition_mc(sets)
  }
  sets <- subsets
  if (k > 2){
    for (i in 1:(k-2)){
      i_selected <- sample(1:2, 1)
      while (length(sets[[i_selected]]) < 2){
        i_selected <- sample(1:2, 1)
      }
      subsets <- partition_mc(sets[[i_selected]])
      sets <- c(sets[-i_selected], subsets)
    }
  }
 sets 
}

coarse_grain_mc(aa_name, 10)
