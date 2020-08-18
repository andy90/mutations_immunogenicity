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

ctl_binded_combinded <- do.call("rbind", list(ctl_A0201_binded))
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

coarse_grain2 <- function(aa){
  g1 <- c(1:2, 14)
  g2 <- c(3:6, 15)
  g3 <- c(7:13)
  g4 <- c(16:17)
  g5 <- c(18:21)
  
  aa_index <- which(aa == aa_name)
  if (aa_index %in% g1){
    aa_coarse <- "R"
  }else if (aa_index %in% g2){
    aa_coarse <- "K"
  }else if (aa_index %in% g3){
    aa_coarse <- "T"
  }else if (aa_index %in% g4){
    aa_coarse <- "I"
  }else{
    aa_coarse <- "M"
  }
  aa_coarse
}

cg_seq2 <-function(r_ctl){
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  pep_array_cg <- sapply(pep_array, coarse_grain2)
  pep_array_cg <- paste(pep_array_cg, collapse = "")
  pep_array_cg
}

ctl_binded_9mer_immunogen <- ctl_binded_9mer[ctl_binded_9mer$Immunogenicity > 0.7, ] # when I am looking at the 0 -> 1 transitions, I should only look at the transition from very non immunogenic to very immunogenic
ctl_binded_9mer_nonimmunogen <- ctl_binded_9mer[ctl_binded_9mer$Immunogenicity < 0.3, ]

ratio <- 0.7
ctl_binded_9mer_immunogen_train <- ctl_binded_9mer_immunogen[sample(1:(ratio*nrow(ctl_binded_9mer_immunogen))), ]
ctl_binded_9mer_nonimmunogen_train <- ctl_binded_9mer_nonimmunogen[sample(1:(ratio*nrow(ctl_binded_9mer_nonimmunogen))), ]

cg_immunogen <- apply(ctl_binded_9mer_immunogen_train, 1, cg_seq)
cg_nonimmunogen <- apply(ctl_binded_9mer_nonimmunogen_train, 1, cg_seq)
intersect(cg_immunogen, cg_nonimmunogen) # the number of intersections between the immunogens and nonimmunogens are surprisingly low after coarse graining 

write.fasta(as.list(cg_immunogen), names = 1:length(cg_immunogen), file.out = "epitopes_cgimmunogen_train.fasta")
write.fasta(as.list(cg_nonimmunogen), names = 1:length(cg_nonimmunogen), file.out = "epitopes_cgnonimmunogen_train.fasta")

single_prob <- function(seq, pos){# this is used to calculate the probability of each coarse grained amino acid at each position
  aa_coarse <- strsplit(seq, "")[[1]][pos]
  coarse_letters <- c("R", "K", "T", "L")
  aa_ind <- which(aa_coarse == coarse_letters)
  p <- rep(0, length(coarse_letters))
  p[aa_ind] <- 1
  p
}


double_prob <- function(seq, pos1, pos2){
  seq_array <- strsplit(seq, "")[[1]]
  aa_coarse1 <- seq_array[pos1]
  aa_coarse2 <- seq_array[pos2]
  coarse_letters <- c("R", "K", "T", "L")
  aa_ind1 <- which(aa_coarse1 == coarse_letters)
  aa_ind2 <- which(aa_coarse2 == coarse_letters)
  l_letters <- length(coarse_letters)
  p <- matrix(rep(0, l_letters*l_letters), nrow = l_letters)
  p[aa_ind1, aa_ind2] <- 1
  p
}

double_prob_neglectlast <- function(seq, pos1, pos2){
  seq_array <- strsplit(seq, "")[[1]]
  aa_coarse1 <- seq_array[pos1]
  aa_coarse2 <- seq_array[pos2]
  coarse_letters <- c("R", "K", "T", "L")
  aa_ind1 <- which(aa_coarse1 == coarse_letters)
  aa_ind2 <- which(aa_coarse2 == coarse_letters)
  l_letters <- length(coarse_letters)
  p <- matrix(rep(0, (l_letters - 1)*(l_letters-1)), nrow = (l_letters - 1))
  if ((aa_ind1 < l_letters) & (aa_ind2 < l_letters)){
    p[aa_ind1, aa_ind2] <- 1
  }
  p
}

double_prob_gauge <- function(seq, pos1, pos2){
  seq_array <- strsplit(seq, "")[[1]]
  aa_coarse1 <- seq_array[pos1]
  aa_coarse2 <- seq_array[pos2]
  coarse_letters <- c("R", "K", "T", "L")
  coarse_letters <- coarse_letters[-gauge_site] # remove the guage size
  aa_ind1 <- which(aa_coarse1 == coarse_letters)
  aa_ind2 <- which(aa_coarse2 == coarse_letters)
  l_letters <- length(coarse_letters)
  p <- matrix(rep(0, l_letters*l_letters), nrow = l_letters)
  p[aa_ind1, aa_ind2] <- 1
  p
}

get_double_probs <- function(seqs, nsites){
  c_letters <- c("R", "K", "T", "L")
  lc_letters <- length(c_letters)
  pij_total <- {}
  for (pos1 in 1:(nsites-1)){
    for (pos2 in (pos1+1):nsites){
      pij <- matrix(rowSums(sapply(seqs, double_prob_gauge, pos1, pos2)), (lc_letters - 1))/length(seqs)
      pij_total <- c(pij_total, as.vector(t(pij)))
    }
  }
  pij_total <- matrix(pij_total, byrow = TRUE, ncol = (lc_letters - 1)*(lc_letters - 1))
}


gauge_site <- 4 # this site will be used as gauge, the h's and J's of this site will be set to 0
single_probs <- sapply(1:9, function(pos){
  rowSums(sapply(cg_immunogen, single_prob, pos))/length(cg_immunogen)
})
double_probs <- get_double_probs(cg_immunogen, 9)
write.table(t(single_probs)[, -gauge_site], file = "epitopes_cgimmunogen_train.p", append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
# remember to delete the last line from the single probabilities since it is not needed

write.table(double_probs, file = "epitopes_cgimmunogen_train.p", append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)

single_probs_nonimmu <- sapply(1:9, function(pos){
  rowSums(sapply(cg_nonimmunogen, single_prob, pos))/length(cg_nonimmunogen)
})
double_probs_nonimmu <- get_double_probs(cg_nonimmunogen, 9)
write.table(t(single_probs_nonimmu)[, -gauge_site], file = "epitopes_cgnonimmunogen_train.p", append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(double_probs_nonimmu, file = "epitopes_cgnonimmunogen_train.p", append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
