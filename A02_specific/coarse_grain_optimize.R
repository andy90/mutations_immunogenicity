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



coarse_grain_4group <- function(aa, aa_ordering, seps){
  sep1 <- seps[1]
  sep2 <- seps[2]
  sep3 <- seps[3]
  g1 <- 1:sep1
  g2 <- (sep1+1):sep2
  g3 <- (sep2+1):sep3
  
  aa_index <- which(aa == aa_ordering$aa_name)
  if (aa_index %in% g1){
    aa_coarse <- "A"
  }else if (aa_index %in% g2){
    aa_coarse <- "B"
  }else if (aa_index %in% g3){
    aa_coarse <- "C"
  }else{
    aa_coarse <- "D"
  }
  aa_coarse
}

cg_seq_4group <-function(r_ctl, aa_ordering, seps){
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  pep_array_cg <- sapply(pep_array, coarse_grain_4group, aa_ordering, seps)
  pep_array_cg <- paste(pep_array_cg, collapse = "")
  pep_array_cg
}

intersect_immu_nonimmu <- function(immunogen, nonimmunogen, aa_ordering, seps){
  cg_immunogen <- apply(immunogen, 1, cg_seq_4group, aa_ordering, seps)
  cg_nonimmunogen <- apply(nonimmunogen, 1, cg_seq_4group, aa_ordering, seps)
  
  #length(intersect(unique(cg_immunogen), unique(cg_nonimmunogen))) # the number of intersections between the immunogens and nonimmunogens
  list(sum(cg_nonimmunogen %in% cg_immunogen), sum(cg_immunogen %in% cg_nonimmunogen))
}


ctl_binded_9mer_immunogen <- ctl_binded_9mer[ctl_binded_9mer$Immunogenicity > 0.7, ] # when I am looking at the 0 -> 1 transitions, I should only look at the transition from very non immunogenic to very immunogenic
ctl_binded_9mer_nonimmunogen <- ctl_binded_9mer[ctl_binded_9mer$Immunogenicity < 0.3, ]

aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W") # this is the aa sequence ordered according to the AA chart on wikipedia
p_woct <- c(1.81, 0.11, 2.80, 3.64, 3.63, 0.46, 0.25, 0.85, 0.77, -0.02, 0, 1.15, 0.14, 0.50, -0.46, -1.12, -1.25, -0.67, -1.71, -0.71, -2.09) # this is the Wimley-White octale scale corresponds to the sequence of the aa_name
woct <- data.frame(aa_name, p_woct)
woct_ordering <- woct[order(woct$p_woct), ]

seps_intersets_woct <- data.frame(sep1 = numeric(), spe2 = numeric(), sep3 = numeric(), n_intersect1 = integer(), n_intersect2 = integer())
for (i in 1:6000){
  seps <- sort(sample(1:21, 3))
  n_intersect <- intersect_immu_nonimmu(ctl_binded_9mer_immunogen, ctl_binded_9mer_nonimmunogen, woct_ordering, seps)
  seps_intersets_woct[nrow(seps_intersets_woct)+1, ] <- c(seps, n_intersect[1], n_intersect[2])
}
seps_intersets_woct <- seps_intersets_woct[order(seps_intersets_woct$n_intersect2), ]

seps_woct_4group_1 <- c(4,7,13) # this is one of the best separatation scheme
seps_woct_4group_2 <- c(5,7,13) # this is another good separation scheme
seps_woct_4group_3 <- c(4, 10, 16) # this is another good spearation scheme

seps_intersets_woct_3group <- data.frame(sep1 = numeric(), spe2 = numeric(), sep3 = numeric(), n_intersect1 = integer(), n_intersect2 = integer())
for (i in 1:1000){
  seps <- sort(sample(1:20, 3))
  seps[3] <- 21
  n_intersect <- intersect_immu_nonimmu(ctl_binded_9mer_immunogen, ctl_binded_9mer_nonimmunogen, woct_ordering, seps)
  seps_intersets_woct_3group[nrow(seps_intersets_woct_3group)+1, ] <- c(seps, n_intersect[1], n_intersect[2])
}
seps_intersets_woct_3group <- seps_intersets_woct_3group[order(seps_intersets_woct_3group$n_intersect2), ]

seps_woct_3group_1 <- c(4, 13, 21) #these two are the good separation schemes if we separate them to 3 groups
seps_woct_3group_2 <- c(6, 13, 21) # there are about 40 repeated sequences

seps_intersets_woct_2group <- data.frame(sep1 = numeric(), spe2 = numeric(), sep3 = numeric(), n_intersect1 = integer(), n_intersect2 = integer())
for (i in 1:100){
  seps <- sort(sample(1:20, 3))
  seps[3] <- 22
  seps[2] <- 21
  n_intersect <- intersect_immu_nonimmu(ctl_binded_9mer_immunogen, ctl_binded_9mer_nonimmunogen, woct_ordering, seps)
  seps_intersets_woct_2group[nrow(seps_intersets_woct_2group)+1, ] <- c(seps, n_intersect[1], n_intersect[2])
}
seps_intersets_woct_2group <- seps_intersets_woct_2group[order(seps_intersets_woct_2group$n_intersect2), ]

seps_woct_2group <- c(6, 21, 22) # if we just divide it into two groups, better divide around 6
#########################
#########################

p_wif <- c(0.81, 0.17, 0.99, 1.23, 2.02, 0.13, 0.14, 0.42, 0.58, -0.24, 0, 0.01, 0.45, 0.17, 0.07, -0.31, -0.56, -0.23, -1.13, -0.94, -1.85)
wif <- data.frame(aa_name, p_wif)
wif_ordering <- wif[order(wif$p_wif), ]

seps_intersets_wif_4group <- data.frame(sep1 = numeric(), spe2 = numeric(), sep3 = numeric(), n_intersect1 = integer(), n_intersect2 = integer())
for (i in 1:6000){
  seps <- sort(sample(1:20, 3))
  n_intersect <- intersect_immu_nonimmu(ctl_binded_9mer_immunogen, ctl_binded_9mer_nonimmunogen, wif_ordering, seps)
  seps_intersets_wif_4group[nrow(seps_intersets_wif_4group)+1, ] <- c(seps, n_intersect[1], n_intersect[2])
}
seps_intersets_wif_4group <- seps_intersets_wif_4group[order(seps_intersets_wif_4group$n_intersect2), ]
seps_wif_4group <- c(6, 10, 13) # there about 25 sharing sequences

seps_intersets_wif_3group <- data.frame(sep1 = numeric(), spe2 = numeric(), sep3 = numeric(), n_intersect1 = integer(), n_intersect2 = integer())
for (i in 1:1000){
  seps <- sort(sample(1:20, 3))
  seps[3] <- 21
  n_intersect <- intersect_immu_nonimmu(ctl_binded_9mer_immunogen, ctl_binded_9mer_nonimmunogen, wif_ordering, seps)
  seps_intersets_wif_3group[nrow(seps_intersets_wif_3group)+1, ] <- c(seps, n_intersect[1], n_intersect[2])
}
seps_intersets_wif_3group <- seps_intersets_wif_3group[order(seps_intersets_wif_3group$n_intersect2), ]
seps_wif_3group <- c(6, 11, 21) # there about 38 repeated sequences

##################
##################
martini <- data.frame(aa_name = character(), polarity = numeric(), volume = numeric(), stringsAsFactors=FALSE)
martini[nrow(martini)+1, ] <- c("A", 0, 0)
martini[nrow(martini)+1, ] <- c("C", 2, 1)
martini[nrow(martini)+1, ] <- c("H", 9, 3)
martini[nrow(martini)+1, ] <- c("M", 2, 1)
martini[nrow(martini)+1, ] <- c("T", 4, 1)
martini[nrow(martini)+1, ] <- c("R", 9, 2)
martini[nrow(martini)+1, ] <- c("Q", 5, 1)
martini[nrow(martini)+1, ] <- c("I", 1, 1)
martini[nrow(martini)+1, ] <- c("F", 3, 3)
martini[nrow(martini)+1, ] <- c("W", 7, 4)
martini[nrow(martini)+1, ] <- c("N", 5, 1)
martini[nrow(martini)+1, ] <- c("E", 6, 1)
martini[nrow(martini)+1, ] <- c("L", 1, 1)
martini[nrow(martini)+1, ] <- c("P", 1, 1)
martini[nrow(martini)+1, ] <- c("Y", 6, 3)
martini[nrow(martini)+1, ] <- c("D", 6, 1)
martini[nrow(martini)+1, ] <- c("G", 0, 0)
martini[nrow(martini)+1, ] <- c("K", 7, 2)
martini[nrow(martini)+1, ] <- c("S", 4, 1)
martini[nrow(martini)+1, ] <- c("V", 1, 1)

martini_polar_ordering <- martini[order(martini$polarity), ]
martini_volume_ordering <- martini[order(martini$volume, decreasing = TRUE), ]

##################
###################
seps_intersets_mpolarity_4group <- data.frame(sep1 = numeric(), spe2 = numeric(), sep3 = numeric(), n_intersect1 = integer(), n_intersect2 = integer())
for (i in 1:1000){
  seps <- sort(sample(1:20, 3))
  seps[1] <- 2
  n_intersect <- intersect_immu_nonimmu(ctl_binded_9mer_immunogen, ctl_binded_9mer_nonimmunogen, martini_polar_ordering, seps)
  seps_intersets_mpolarity_4group[nrow(seps_intersets_mpolarity_4group)+1, ] <- c(seps, n_intersect[1], n_intersect[2])
}
seps_intersets_mpolarity_4group <- seps_intersets_mpolarity_4group[order(seps_intersets_mpolarity_4group$n_intersect2), ]
seps_mpolarity_4group <- c() 

seps_intersets_mpolarity_3group <- data.frame(sep1 = numeric(), spe2 = numeric(), sep3 = numeric(), n_intersect1 = integer(), n_intersect2 = integer())
for (i in 1:1000){
  seps <- sort(sample(1:20, 3))
  seps[3] <- 20
  n_intersect <- intersect_immu_nonimmu(ctl_binded_9mer_immunogen, ctl_binded_9mer_nonimmunogen, martini_polar_ordering, seps)
  seps_intersets_mpolarity_3group[nrow(seps_intersets_mpolarity_3group)+1, ] <- c(seps, n_intersect[1], n_intersect[2])
}
seps_intersets_mpolarity_3group <- seps_intersets_mpolarity_3group[order(seps_intersets_mpolarity_3group$n_intersect2), ]
seps_mpolarity_3group <- c()

# there are two much degeneracy from martini force field
###################
###################
AA_volume_ordering  <- data.frame(aa_name = character(), volume = integer(), stringsAsFactors=FALSE)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("G", 1)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("A", 2)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("S", 3)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("C", 4)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("D", 5)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("P", 6)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("N", 7)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("T", 8)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("E", 9)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("V", 10)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("Q", 11)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("H", 12)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("M", 13)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("I", 14)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("L", 15)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("K", 16)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("R", 17)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("F", 18)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("Y", 19)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("W", 20)

AA_volume_ordering <- AA_volume_ordering[order(AA_volume_ordering$volume, decreasing = TRUE), ]
seps_intersets_volume_4group <- data.frame(sep1 = numeric(), spe2 = numeric(), sep3 = numeric(), n_intersect1 = integer(), n_intersect2 = integer())
for (i in 1:6000){
  seps <- sort(sample(1:20, 3))
  n_intersect <- intersect_immu_nonimmu(ctl_binded_9mer_immunogen, ctl_binded_9mer_nonimmunogen, AA_volume_ordering, seps)
  seps_intersets_volume_4group[nrow(seps_intersets_volume_4group)+1, ] <- c(seps, n_intersect[1], n_intersect[2])
}
seps_intersets_volume_4group <- seps_intersets_volume_4group[order(seps_intersets_volume_4group$n_intersect2), ]

seps_volume_4group_1 <- c(9, 11, 18)
seps_volume_4group_2 <- c(9, 12, 18)
seps_volume_4group_3 <- c(7, 11, 18)
seps_volume_4group_4 <- c(7, 12, 18) # there are about 20 overlaps for the 4 different ways of separating
seps_volume_4group_5 <- c(7, 12, 17)
intersect_immu_nonimmu(ctl_binded_9mer_immunogen, ctl_binded_9mer_nonimmunogen, AA_volume_ordering, c(4,9,16)) # there are 35 overlapped sequences

#################
################
# for the polarity, what we obtained is very coarse grained
AA_polarity_ordering  <- data.frame(aa_name = character(), polarity = integer(), stringsAsFactors=FALSE)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("R", 1)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("H", 1)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("K", 1)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("D", 1)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("E", 1)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("S", 2)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("T", 2)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("N", 2)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("Q", 2)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("Y", 2)
#AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("A", 3)
#AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("G", 3)
#AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("P", 3)
#AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("C", 3)
#AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("M", 3)
#AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("V", 4)
#AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("I", 4)
#AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("L", 4)
#AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("F", 4)
#AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("W", 4)

AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("A", 3)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("G", 3)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("V", 3)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("I", 3)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("L", 3)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("P", 4)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("C", 4)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("M", 4)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("F", 4)
AA_polarity_ordering[nrow(AA_polarity_ordering)+1, ] <- list("W", 4)


AA_polarity_ordering <- AA_polarity_ordering[order(AA_polarity_ordering$polarity, decreasing = TRUE), ]
seps_polarity_4group <- c(5, 10, 15)
n_intersect <- intersect_immu_nonimmu(ctl_binded_9mer_immunogen, ctl_binded_9mer_nonimmunogen, AA_polarity_ordering, seps_polarity_4group) # there are 35 overlapped sequences


#####################
#####################
test_overlap <- function(seqs1, seqs2){
  overlap_status <- sapply(seqs1, function(seq){
    seq %in% seqs2
  })
}
cg_woct_immunogen <- apply(ctl_binded_9mer_immunogen, 1, cg_seq_4group, woct_ordering, c(4, 10, 16))
cg_woct_nonimmunogen <- apply(ctl_binded_9mer_nonimmunogen, 1, cg_seq_4group, woct_ordering, c(4, 10, 16))

cg_volume_immunogen <- apply(ctl_binded_9mer_immunogen, 1, cg_seq_4group, AA_volume_ordering, c(7, 11, 16))
cg_volume_nonimmunogen <- apply(ctl_binded_9mer_nonimmunogen, 1, cg_seq_4group, AA_volume_ordering, c(7, 11, 16))

cg_polarity_immunogen <- apply(ctl_binded_9mer_immunogen, 1, cg_seq_4group, AA_polarity_ordering, seps_polarity_4group)
cg_polarity_nonimmunogen <- apply(ctl_binded_9mer_nonimmunogen, 1, cg_seq_4group, AA_polarity_ordering, seps_polarity_4group)

immunogen_distinguishability <- data.frame(seq = ctl_binded_9mer_immunogen$peptide, hydrophobic_cg = test_overlap(cg_woct_immunogen, cg_woct_nonimmunogen), volume_cg = test_overlap(cg_volume_immunogen, cg_volume_nonimmunogen), polarity_cg = test_overlap(cg_polarity_immunogen, cg_polarity_nonimmunogen))
nonimmunogen_distinguishability <- data.frame(seq = ctl_binded_9mer_nonimmunogen$peptide, hydrophobic_cg = test_overlap(cg_woct_nonimmunogen, cg_woct_immunogen), volume_cg = test_overlap(cg_volume_nonimmunogen, cg_volume_immunogen), polarity_cg = test_overlap(cg_polarity_nonimmunogen, cg_polarity_immunogen))

which((nonimmunogen_distinguishability$hydrophobic_cg + nonimmunogen_distinguishability$volume_cg +nonimmunogen_distinguishability$polarity_cg) == 3)
which((immunogen_distinguishability$hydrophobic_cg + immunogen_distinguishability$volume_cg +immunogen_distinguishability$polarity_cg) == 3)

######################
######################

single_prob <- function(seq, pos){# this is used to calculate the probability of each coarse grained amino acid at each position
  aa_coarse <- strsplit(seq, "")[[1]][pos]
  coarse_letters <- c("A", "B", "C", "D")
  aa_ind <- which(aa_coarse == coarse_letters)
  p <- rep(0, length(coarse_letters))
  p[aa_ind] <- 1
  p
}


double_prob <- function(seq, pos1, pos2){
  seq_array <- strsplit(seq, "")[[1]]
  aa_coarse1 <- seq_array[pos1]
  aa_coarse2 <- seq_array[pos2]
  coarse_letters <- c("A", "B", "C", "D")
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
  coarse_letters <- c("A", "B", "C", "D")
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
  coarse_letters <- c("A", "B", "C", "D")
  coarse_letters <- coarse_letters[-gauge_site] # remove the guage size
  aa_ind1 <- which(aa_coarse1 == coarse_letters)
  aa_ind2 <- which(aa_coarse2 == coarse_letters)
  l_letters <- length(coarse_letters)
  p <- matrix(rep(0, l_letters*l_letters), nrow = l_letters)
  p[aa_ind1, aa_ind2] <- 1
  p
}

get_single_probs <- function(seqs, nsites){
  sapply(1:nsites, function(pos){
    rowSums(sapply(seqs, single_prob, pos))/length(seqs)
  })
}

get_double_probs <- function(seqs, nsites){
  c_letters <- c("A", "B", "C", "D")
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

write_probs <- function(gauge_site, seqs, file_name){
  single_probs <- get_single_probs(seqs, 9)
  double_probs <- get_double_probs(seqs, 9)
  write.table(t(single_probs)[, -gauge_site], file = file_name, append = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  # remember to delete the last line from the single probabilities since it is not needed
  
  write.table(double_probs, file = file_name, append = TRUE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

gauge_site <- 4 # this site will be used as gauge, the h's and J's of this site will be set to 0


write_probs(gauge_site, cg_woct_immunogen, "cg_immu_woct_full.p")
write_probs(gauge_site, cg_woct_nonimmunogen, "cg_nonimmu_woct_full.p")

write_probs(gauge_site, cg_volume_immunogen, "cg_immu_vol_full.p")
write_probs(gauge_site, cg_volume_nonimmunogen, "cg_nonimmu_vol_full.p")

write_probs(gauge_site, cg_polarity_immunogen, "cg_immu_pol_full.p")
write_probs(gauge_site, cg_polarity_nonimmunogen, "cg_nonimmu_pol_full.p")

write_probs(gauge_site, sample(cg_woct_immunogen, 48), "cg_immu_woct_half.p")
write_probs(gauge_site, sample(cg_woct_nonimmunogen, 19), "cg_nonimmu_woct_half.p")

write_probs(gauge_site, sample(cg_volume_immunogen, 48), "cg_immu_vol_half.p")
write_probs(gauge_site, sample(cg_volume_nonimmunogen, 19), "cg_nonimmu_vol_half.p")

write_probs(gauge_site, sample(cg_polarity_immunogen, 48), "cg_immu_pol_half.p")
write_probs(gauge_site, sample(cg_polarity_nonimmunogen, 19), "cg_nonimmu_pol_half.p")


