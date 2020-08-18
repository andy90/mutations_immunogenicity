rm(list=ls()) # clean the memeory

compare_pep <- function(pep1, pep2){
  pep1_array <- strsplit(as.character(pep1), "")[[1]]
  pep2_array <- strsplit(as.character(pep2), "")[[1]]
  sum(pep1_array != pep2_array)
}

get_mdiff <- function(ctl_base, ctl_mut){
  m_diff <- sapply(ctl_mut$peptide, function(pep1){
    which.min(sapply(ctl_base$peptide, function(pep2){
      compare_pep(pep1, pep2)
    }))
  })
}

get_ndiff <- function(ctl_base, ctl_mut){
  n_diff <- sapply(ctl_mut$peptide, function(pep1){
    min(sapply(ctl_base$peptide, function(pep2){
      compare_pep(pep1, pep2)
    }))
  })
}

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

ctl_base <- ctl_binded_9mer[ctl_binded_9mer$peptide == "SLYNTVATL", ]

n_diff <- get_ndiff(ctl_base, ctl_binded_9mer)

ctl_local_to_base <- ctl_binded_9mer[(n_diff < 4) , ]

ctl_local_to_base <- ctl_local_to_base[(ctl_local_to_base$Immunogenicity > 0.7) | (ctl_local_to_base$Immunogenicity < 0.3), ]


##################################
# construct amino acid scales
##################################
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
AA_polarity_ordering$polarity <- (AA_polarity_ordering$polarity - mean(AA_polarity_ordering$polarity))/(max(AA_polarity_ordering$polarity) - min(AA_polarity_ordering$polarity))

AA_volume_ordering  <- data.frame(aa_name = character(), volume = integer(), stringsAsFactors=FALSE)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("G", 60.1)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("A", 88.6)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("S", 89.0)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("C", 108.5)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("D", 111.1)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("P", 112.7)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("N", 114.1)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("T", 116.1)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("E", 138.4)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("V", 140.0)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("Q", 143.8)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("H", 153.2)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("M", 162.9)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("I", 166.7)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("L", 166.7)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("K", 168.6)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("R", 173.4)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("F", 189.9)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("Y", 193.6)
# AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("W", 227.8)

AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("G", 1)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("A", 1)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("S", 1)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("C", 2)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("D", 2)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("P", 2)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("N", 3)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("T", 3)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("E", 3)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("V", 3)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("Q", 4)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("H", 4)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("M", 5)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("I", 5)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("L", 5)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("K", 5)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("R", 5)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("F", 6)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("Y", 6)
AA_volume_ordering[nrow(AA_volume_ordering)+1, ] <- list("W", 6)
AA_volume_ordering$volume <- (AA_volume_ordering$volume - mean(AA_volume_ordering$volume))/(max(AA_volume_ordering$volume) - min(AA_volume_ordering$volume))

aa_scales <- merge(AA_volume_ordering, AA_polarity_ordering)


##############################
# convert the peptide to an array
###############################
peptide_to_array <- function(pep){
  pep_separated <- strsplit(as.character(pep), "")[[1]]
  res <- sapply(pep_separated, function(aa){
    i_index <- which(aa == aa_scales$aa_name)
    x <- rep(0, length(aa_scales$aa_name))
    x[i_index] <- 1
    x
  })
  as.vector(res)
}

mut <- t(sapply(ctl_local_to_base$peptide, peptide_to_array))
wt <- peptide_to_array("SLYNTVATL")
diff <- sweep(mut, 2, wt)

#########################
# now we can define the matrix related 
#########################
p <- aa_scales$volume


p_subfill <- numeric(0)
for (i in 1:9){
  x <- rep(0, 9*length(p))
  x[((i-1)*length(p) + 1) : (i*length(p))] <- p
  p_subfill <- cbind(p_subfill, x)
}

diff_vol_scores <- diff %*% p_subfill # this is the same as what we get using the seq_to_scores function

#########################
# now we can define the cost function based on diff and p_subfill. I will first do that to learn the position dependent parameters. Then I will learn p with fixed x
#########################
