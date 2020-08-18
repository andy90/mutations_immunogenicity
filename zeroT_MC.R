rm(list=ls()) # clean the memeory
library("seqinr")
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

aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W")

coarse_grain <- function(aa){
  g1 <- c(1:2, 14)
  g2 <- c(3:6, 15)
  g3 <- c(7:13, 16)
  g4 <- c(17:21)
  
  aa_index <- which(aa == aa_name)
  if (aa_index %in% g1){
    aa_coarse <- "R"
  }else if (aa_index %in% g2){
    aa_coarse <- "K"
  }else if (aa_index %in% g3){
    aa_coarse <- "T"
  }else{
    aa_coarse <- "L"
  }
  aa_coarse
}

cg_seq <-function(r_ctl){
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  pep_array_cg <- sapply(pep_array, coarse_grain)
  pep_array_cg <- paste(pep_array_cg, collapse = "")
  pep_array_cg
}



ctl_binded_9mer_immunogen <- ctl_binded_9mer[ctl_binded_9mer$Immunogenicity > 0.7, ] # when I am looking at the 0 -> 1 transitions, I should only look at the transition from very non immunogenic to very immunogenic
ctl_binded_9mer_nonimmunogen <- ctl_binded_9mer[ctl_binded_9mer$Immunogenicity < 0.3, ]

cg_immunogen <- apply(ctl_binded_9mer_immunogen, 1, cg_seq)
cg_nonimmunogen <- apply(ctl_binded_9mer_nonimmunogen, 1, cg_seq)

#################
# in the above part, we read in all the sequences and coarse grain them according to our methods
################

################
# in the part below we implement the function to calculate the energy of sequences
################
aa_cg_name <- c("R", "K", "T", "L") # the coarse grained aa names

seq_to_array <- function(pep_seq){ # convert character sequences into array of numbers
  pep_seq_array <- strsplit(pep_seq, "")[[1]]
  pep_seq_numarray <- sapply(pep_seq_array, function(aa){
    ind <- which(aa == aa_cg_name)
    n_aa <- length(aa_cg_name)
    aa_numarray <- rep(0, n_aa)
    aa_numarray[ind] <- 1
    aa_numarray
  })
  as.vector(pep_seq_numarray)
}

array_to_seq <- function(pep_num_array){ 
  pep_matrix <- matrix(pep_num_array, nrow = n_letters)
  pep_seq <- apply(pep_matrix, 2, function(sub_array){
    ind <- which(sub_array == 1)
    aa_cg_name[ind]
  })
  pep_seq
}


library("matrixcalc")
library("lattice")

n_letters <- length(aa_cg_name) # number of coarse grained groups
n_sites <- 9 # number of sites
index_map <- {}
for (i in 1:(n_sites - 1)){
  for (j in (i+1):n_sites){
    line_num <- (i-1)*(n_sites-2) - (i-1)*(i-2)/2 + j - 1
    index_map <- rbind(index_map, c(i, j, line_num))
  }
} 

getJ <- function(J_reshaped){ # reshape J into a square matrix
  J <- matrix(0, nrow = n_letters * n_sites, ncol = n_letters * n_sites)
  for (line_num in 1:nrow(J_reshaped)){
    ind <- which(line_num == index_map[ ,3])
    i <- index_map[ind, 1]
    j <- index_map[ind, 2]
    
    J <- set.submatrix(J, matrix(J_reshaped[line_num,], byrow = TRUE, ncol = n_letters), (n_letters*(i-1)+1), (n_letters*(j-1)+1))
    J <- set.submatrix(J, matrix(J_reshaped[line_num,], byrow = TRUE, ncol = n_letters), (n_letters*(j-1)+1), (n_letters*(i-1)+1)) 
  }
  J
}

J_raw_immuno <- data.matrix(read.table("J_ZS_immu.txt"))
J_immuno <- getJ(J_raw_immuno)

h_raw_immuno <- data.matrix(read.table("h_ZS_immu.txt"))
h_immuno <- as.vector(t(h_raw_immuno)) # convert h to an array of 36 length

seq_energy <- function(pep_num_array){
  E <- h_immuno %*% pep_num_array + 0.5 * pep_num_array %*% J_immuno %*% pep_num_array
  E
}

seq_energy_single <- function(pep_num_array){ # this is the energy coming from the single field
  E <- h_immuno %*% pep_num_array 
  E
}

seq_energy_double <- function(pep_num_array){
  E <- 0.5 * pep_num_array %*% J_immuno %*% pep_num_array
  E
}

#################
# in the part below we will try to construct the zero temperatue monte carlo move
################

change_AA <- function(pep_num_array){
  isite <- sample(1:n_sites, 1) # determine which site to move
  iAA <- sample(1:n_letters, 1) # determine the new amino acid (coarse graind)
  sub_array <- rep(0, n_letters)
  sub_array[iAA] <- 1 
  new_numarray <- pep_num_array
  new_numarray[((isite-1)*n_letters+1) : (isite*n_letters)] <- sub_array
  new_numarray
}

change_AA_fixSiteAA <- function(pep_num_array, isite, iAA){
  sub_array <- rep(0, n_letters)
  sub_array[iAA] <- 1 
  new_numarray <- pep_num_array
  new_numarray[((isite-1)*n_letters+1) : (isite*n_letters)] <- sub_array
  new_numarray
}

change_AA_descent <- function(pep_num_array){ # this implements the gradient descent

  E1 <- seq_energy(pep_num_array)
  pep_descent <- pep_num_array
  for (isite in 1:n_sites){
    for (iAA in 1:n_letters){
      pep_new <- change_AA_fixSiteAA(pep_num_array, isite, iAA)
      E2 <- seq_energy(pep_new)
      if (E2 < E1){
        E1 <- E2
        pep_descent <- pep_new
      }
    }
  }
  pep_descent
}

mc_zero <- function(pep_cg_init, nsteps){ # given an initial array, find the local minimum using zero temperature monte carlo
  pep_num <- seq_to_array(pep_cg_init)
  E0 <- seq_energy(pep_num)
  for (i in 1:nsteps){
    pep_new <- change_AA_descent(pep_num)
    E1 <- seq_energy(pep_new)
    if (E1 < E0){
      pep_num <- pep_new
      E0 <- E1
    }
    print(list(E0, array_to_seq(pep_num)))
  }
  array_to_seq(pep_num)
}

mc_zero(cg_immunogen[130],2000)

