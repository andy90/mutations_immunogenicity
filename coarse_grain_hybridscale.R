# this r script is used to find out a coarse grain scale that can best separate the amino acids
rm(list=ls()) # clean the memeory


###################
# the first step is to produce the matrix that documents the transition that will cause immunogenicity from 0 to 1
##################

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

mut_m <- matrix(integer(21*21), nrow = 21, ncol = 21) # this is a template for the mutation matrix, used as zero point to be passed to find_singlemut
aa_name <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W")
colnames(mut_m) <- aa_name
rownames(mut_m) <- aa_name



A <- sum_singlemut_immuincrease_noround(ctl_binded_9mer, mut_m)
Am <- matrix(rowSums(A), nrow = 21) # 0 to 1 # there are total of 67 0 to 1 transfer
effective_mutations <- which(Am !=0, arr.ind = T) # this is a matrix, for each row, the first element is the starting aa, the second element is the end aa

###################
# define several scales for amino acid
###################
p_woct <- c(1.81, 0.11, 2.80, 3.64, 3.63, 0.46, 0.25, 0.85, 0.77, -0.02, 0, 1.15, 0.14, 0.50, -0.46, -1.12, -1.25, -0.67, -1.71, -0.71, -2.09) # this is the Wimley-White octale scale corresponds to the sequence of the aa_name
woct <- data.frame(aa_name, p_woct)
woct$p_woct <- woct$p_woct/max(woct$p_woct)

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
AA_volume_ordering$volume <- AA_volume_ordering$volume/max(AA_volume_ordering$volume)

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
AA_polarity_ordering$polarity <- AA_polarity_ordering$polarity/max(AA_polarity_ordering$polarity)
####################
# define a function to calculate the distance between all the mutations
####################
aa_scale <- function(a1, a2, a3, aa){ # a1, a2 and a3 are the coefficients for the aa scale
  a1 * woct[aa_name == aa, 2] + a2 * AA_volume_ordering[aa_name == aa, 2] + a3 * AA_polarity_ordering[aa_name == aa, 2]
}

aa_scale_fixa1 <- function(a2, a3, aa){ #  a2 and a3 are the coefficients for the aa scale
   0*woct[woct$aa_name == aa, 2] + a2 * AA_volume_ordering[AA_volume_ordering$aa_name == aa, 2] + a3 * AA_polarity_ordering[AA_polarity_ordering$aa_name == aa, 2]
}

total_distance <- function(par){
  a2 <- par[1]
  a3 <- par[2]
  sum(apply(effective_mutations, 1, function(aas){
    aa_before <- aa_name[aas[1]]
    aa_after <- aa_name[aas[2]]
    abs(aa_scale_fixa1( a2, a3, aa_before) - aa_scale_fixa1(a2, a3, aa_after))
  }))
}

average_distance <- function(par){
  a2 <- par[1]
  a3 <- par[2]
  aa_name_noU <- aa_name[aa_name != 'U']
  sum(sapply(aa_name_noU, function(aa_before){
    sapply(aa_name_noU, function(aa_after){
      abs(aa_scale_fixa1(a2, a3, aa_before) - aa_scale_fixa1(a2, a3, aa_after))
    })
  }))/(length(aa_name_noU) * length(aa_name_noU))
}



negative_total_distance <- function(par){  # assume a1 is 1
  -total_distance(par)/average_distance(par)
}
res <- optim(c(1,1), negative_total_distance)

new_scale <- numeric()
j <- 1
for (i in aa_name[aa_name != 'U']){
  new_scale[j] <- aa_scale_fixa1(res$par[1], res$par[2], i)
  j <- j+1
}
aa_new_scale <- data.frame(aa_name = aa_name[aa_name != 'U'], scale = new_scale)
aa_new_scale_ordered <- aa_new_scale[order(aa_new_scale$scale), ]
aa_new_scale$scale <- aa_new_scale$scale/max(aa_new_scale$scale)
write.table(aa_new_scale, "combined_aa_scale.csv")
##########################
# now coarse grain the amino acids and test the overlap between the immuno and non-immunogens
##########################

coarse_grain_5group <- function(aa, aa_ordering, seps){
  sep1 <- seps[1]
  sep2 <- seps[2]
  sep3 <- seps[3]
  sep4 <- seps[4]
  
  g1 <- 1:sep1
  g2 <- (sep1+1):sep2
  g3 <- (sep2+1):sep3
  g4 <- (sep3+1):sep4
  
  aa_index <- which(aa == aa_ordering$aa_name)
  if (aa_index %in% g1){
    aa_coarse <- "A"
  }else if (aa_index %in% g2){
    aa_coarse <- "B"
  }else if (aa_index %in% g3){
    aa_coarse <- "C"
  }else if (aa_index %in% g4){
    aa_coarse <- "D"
  }else{
    aa_coarse <- "E"
  }
  aa_coarse
}

cg_seq_5group <-function(r_ctl, aa_ordering, seps){
  pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
  pep_array_cg <- sapply(pep_array, coarse_grain_5group, aa_ordering, seps)
  pep_array_cg <- paste(pep_array_cg, collapse = "")
  pep_array_cg
}

intersect_immu_nonimmu <- function(immunogen, nonimmunogen, aa_ordering, seps){
  cg_immunogen <- apply(immunogen, 1, cg_seq_5group, aa_ordering, seps)
  cg_nonimmunogen <- apply(nonimmunogen, 1, cg_seq_5group, aa_ordering, seps)
  
  list(sum(cg_nonimmunogen %in% cg_immunogen), sum(cg_immunogen %in% cg_nonimmunogen))
}


ctl_binded_9mer_immunogen <- ctl_binded_9mer[ctl_binded_9mer$Immunogenicity > 0.7, ] # when I am looking at the 0 -> 1 transitions, I should only look at the transition from very non immunogenic to very immunogenic
ctl_binded_9mer_nonimmunogen <- ctl_binded_9mer[ctl_binded_9mer$Immunogenicity < 0.3, ]

seps_intersets_5group <- data.frame(sep1 = numeric(), spe2 = numeric(), sep3 = numeric(), sep4 = numeric(), n_intersect1 = integer(), n_intersect2 = integer())
for (i in 1:6000){
  seps <- sort(sample(1:20, 4))
  n_intersect <- intersect_immu_nonimmu(ctl_binded_9mer_immunogen, ctl_binded_9mer_nonimmunogen, aa_new_scale_ordered, seps)
  seps_intersets_5group[nrow(seps_intersets_5group)+1, ] <- c(seps, n_intersect[1], n_intersect[2])
}
seps_intersets_5group <- seps_intersets_5group[order(seps_intersets_5group$n_intersect2), ]


####################
# we can count the occurance frequency of each amino acid in the peptides
#####################
count_freq <- function(ctl){
  rowSums(apply(ctl, 1, function(r_ctl){
    pep_array <- strsplit(as.character(r_ctl[1]), "")[[1]]
    aa_name %in% pep_array
  }))
}

aa_freq <- data.frame(name = aa_name, freq = count_freq(ctl_binded_9mer))
