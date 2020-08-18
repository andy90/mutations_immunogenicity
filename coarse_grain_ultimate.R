rm(list=ls()) # clean the memeory


coarse_grain <- function(aa, aa_ordering, seps){ # this function can coarse grain the amino acids into arbitrary groups
  
  aa_index <- which(aa == aa_ordering$aa_name) # aa_ordering order the amino acids according to their similarities
  seps_ordered <- sort(seps)
  i_interval <- sum(aa_index > seps_ordered) + 1 # determine which interval should the amino acid fall into. seps define left open right closed interval (sep1, sep2]  
  
  aa_coarse <- LETTERS[i_interval]
  aa_coarse
}

cg_seq <-function(r_ctl, aa_ordering, seps){
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  pep_array_cg <- sapply(pep_array, coarse_grain, aa_ordering, seps)
  pep_array_cg <- paste(pep_array_cg, collapse = "")
  pep_array_cg
}

intersect_immu_nonimmu <- function(immunogen, nonimmunogen, aa_ordering, seps){
  cg_immunogen <- apply(immunogen, 1, cg_seq, aa_ordering, seps)
  cg_nonimmunogen <- apply(nonimmunogen, 1, cg_seq, aa_ordering, seps)
  
  list(sum(cg_nonimmunogen %in% cg_immunogen), sum(cg_immunogen %in% cg_nonimmunogen))
}

single_prob <- function(seq, pos, nletters){# this is used to calculate the probability of each coarse grained amino acid at each position, nletters is the number of coarse grained letters we used 
  aa_coarse <- strsplit(seq, "")[[1]][pos]
  aa_ind <- which(aa_coarse == LETTERS)
  p <- rep(0, nletters)
  p[aa_ind] <- 1
  p
}


double_prob_gauge <- function(seq, pos1, pos2, nletters){
  seq_array <- strsplit(seq, "")[[1]]
  aa_coarse1 <- seq_array[pos1]
  aa_coarse2 <- seq_array[pos2]
  coarse_letters <- LETTERS[1 : nletters]
  coarse_letters <- coarse_letters[-gauge_site] # remove the guage size, this applies to whichever gauge site we chose
  aa_ind1 <- which(aa_coarse1 == coarse_letters)
  aa_ind2 <- which(aa_coarse2 == coarse_letters)
  l_letters <- length(coarse_letters)
  p <- matrix(rep(0, l_letters*l_letters), nrow = l_letters)
  p[aa_ind1, aa_ind2] <- 1
  p
}


get_single_probs <- function(seqs, nsites, nletters){
  sapply(1:nsites, function(pos){
    rowSums(sapply(seqs, single_prob, pos, nletters))/length(seqs)
  })
}

get_double_probs <- function(seqs, nsites, nletters){
  lc_letters <- nletters
  pij_total <- {}
  for (pos1 in 1:(nsites-1)){
    for (pos2 in (pos1+1):nsites){
      pij <- matrix(rowSums(sapply(seqs, double_prob_gauge, pos1, pos2, nletters)), (lc_letters - 1))/length(seqs)
      pij_total <- c(pij_total, as.vector(t(pij)))
    }
  }
  pij_total <- matrix(pij_total, byrow = TRUE, ncol = (lc_letters - 1)*(lc_letters - 1))
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

sum_singlemut_immuincrease_noround <- function(ctl_binded_9mer, mut_m){
  apply(ctl_binded_9mer, 1, function(r_ctl){
    Ai <- apply(ctl_binded_9mer, 1, find_singlemut_immuincrease_noround, r_ctl, mut_m)
    Ai <- rowSums(Ai)
  })  
}

mut_m <- matrix(integer(20*20), nrow = 20, ncol = 20) # this is a template for the mutation matrix, used as zero point to be passed to find_singlemut
aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # I have deleted amino acid U from the list
colnames(mut_m) <- aa_name
rownames(mut_m) <- aa_name


##############################
# read in data
##############################

ctl_A0201 <- read.csv("ctl_A0201_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele A0201
ctl_B5701 <- read.csv("ctl_B5701_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele B5701
ctl_A0301 <- read.csv("ctl_A0301_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele A0301
ctl_B0702 <- read.csv("ctl_B0702_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele B0702

# only work on those peptides that bind to the MHC
bind_criteria <- 50000000 # the binding IC50 is 500nM
ctl_A0201_binded <- ctl_A0201[ctl_A0201$ic50 < bind_criteria, ]
ctl_B5701_binded <- ctl_B5701[ctl_B5701$ic50 < bind_criteria, ]
ctl_A0301_binded <- ctl_A0301[ctl_A0301$ic50 < bind_criteria, ]
ctl_B0702_binded <- ctl_B0702[ctl_B0702$ic50 < bind_criteria, ]

ctl_binded_combinded <- do.call("rbind", list(ctl_A0201_binded, ctl_A0301_binded, ctl_B0702_binded, ctl_B5701_binded))

ctl_binded_9mer <- ctl_binded_combinded[nchar(ctl_binded_combinded$peptide)==9, ]

A <- sum_singlemut_immuincrease_noround(ctl_binded_9mer, mut_m)
Am <- matrix(rowSums(A), nrow = 20) # 0 to 1 # there are total of 67 0 to 1 transfer
effective_mutations <- which(Am !=0, arr.ind = T) # this is a matrix, for each row, the first element is the starting aa, the second element is the end aa



##############################
# now set the functions to do the MC simulation clustering the amino acids
#############################
calc_distance <- function(AA_xyz){ # AA_xyz is the coordinate of the amino acids, it is a 20 by 2 matrix
  apply(AA_xyz, 1, function(ai_xyz){
    apply(AA_xyz, 1, function(aj_xyz){
      dist(rbind(ai_xyz, aj_xyz))
    })
  })
} 

mc_move <- function(AA_xyz, Am, nsteps){
  AA_dis <- calc_distance(AA_xyz)
  rcut_rep <- 0.1
  rcut_attract <- 0.05
  E <- sum((AA_dis < rcut_rep) * Am - (AA_dis < rcut_attract)/2)
  for (i in 1:nsteps){
    AA_xyz_new <- AA_xyz
    ia_change <- sample(1:nrow(AA_xyz), 1)
    AA_xyz_new[ia_change, ] <- runif(2) # we confine to the interval from 0 to 1 
    AA_dis_new <- calc_distance(AA_xyz_new)
    E_new <- sum((AA_dis_new < rcut_rep) * Am - (AA_dis_new < rcut_attract)/2)
    if (E_new <= E){
      AA_xyz <- AA_xyz_new
      E_new <- E
    }
  }
  AA_xyz
}

AA_xyz_initial <- matrix(runif(40), ncol = 2)
AA_xyz_final <- mc_move(AA_xyz_initial, Am, 1000)
calc_distance(AA_xyz_initial)

###################
# one possible way of separation amino acids, with very minimal overlap
####################
aa_manual_ordering <- aa_name[c(5,10,13,17,20, 1,2,4,12,14, 8,9,11,15,18, 3,6,7,16,19)]
aa_manual_seps <- c(5, 10, 15)
aa_manual <- data.frame(aa_name = aa_manual_ordering)

coarse_grain("Y", aa_ordering = aa_manual, seps = aa_manual_seps)

#################
# now do the coarse graining
#################
ctl_binded_9mer_immunogen <- ctl_binded_9mer[ctl_binded_9mer$Immunogenicity > 0.7, ] # when I am looking at the 0 -> 1 transitions, I should only look at the transition from very non immunogenic to very immunogenic
ctl_binded_9mer_nonimmunogen <- ctl_binded_9mer[ctl_binded_9mer$Immunogenicity < 0.3, ]

cg_immunogen <- apply(ctl_binded_9mer_immunogen, 1, cg_seq, aa_manual, aa_manual_seps)
cg_nonimmunogen <- apply(ctl_binded_9mer_nonimmunogen, 1, cg_seq, aa_manual, aa_manual_seps)
list(sum(cg_nonimmunogen %in% cg_immunogen), sum(cg_immunogen %in% cg_nonimmunogen)) # the intersect function only lists unique overlapping elements, that's why I do this function, which has information both about the total number of nonimmunogens in immunogen and immunogens in nonimmunogen
# with this coarse graining scheme, there is no overlap between the immunogens and nonimmunogens
library("seqinr")
write.fasta(as.list(gsub("B", "E", cg_immunogen)), names = 1:length(cg_immunogen), file.out = "epitopes_cgimmunogen.fasta")
write.fasta(as.list(gsub("B", "E", cg_nonimmunogen)), names = 1:length(cg_nonimmunogen), file.out = "epitopes_cgnonimmunogen.fasta")

######################
# write out the single probability and double probabilities
#######################
gauge_site <- 4 # this site will be used as gauge, the h's and J's of this site will be set to 0
single_probs <- get_single_probs(cg_immunogen, 9, 4)
double_probs <- get_double_probs(cg_immunogen, 9, 4)

single_probs_nonimmu <- get_single_probs(cg_nonimmunogen, 9, 4)
double_probs_nonimmu <- get_double_probs(cg_nonimmunogen, 9, 4)

write.table(t(single_probs)[, -gauge_site], file = "epitopes_cgimmunogen.p", append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# remember to delete the last line from the single probabilities since it is not needed

write.table(double_probs, file = "epitopes_cgimmunogen.p", append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(t(single_probs_nonimmu)[, -gauge_site], file = "epitopes_cgnonimmunogen.p", append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(double_probs_nonimmu, file = "epitopes_cgnonimmunogen.p", append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)


#######################
# now we start calculate the energy of each peptide
########################

library("matrixcalc")
library("lattice")

n_letters <- 4 # number of coarse grained groups
n_sites <- 9 # number of sites

seq_to_array <- function(pep_seq){ # convert character sequences into array of numbers
  pep_seq_array <- strsplit(pep_seq, "")[[1]]
  pep_seq_numarray <- sapply(pep_seq_array, function(aa){
    ind <- which(aa == LETTERS)
    n_aa <- n_letters
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
    LETTERS[ind]
  })
  pep_seq
}


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

sapply(cg_immunogen, function(pep){
  seq_energy(seq_to_array(pep))
})

sapply(cg_nonimmunogen, function(pep){
  seq_energy(seq_to_array(pep))
})
