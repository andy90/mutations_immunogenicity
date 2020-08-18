# this r script is used to train a model which uses a real number to mark each amino acid, and also used the neural network
rm(list=ls()) # clean the memeory
###################
# the first step is to read in the peptides
##################

ctl_A0201 <- read.csv("ctl_A0201_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele A0201
ctl_B5701 <- read.csv("ctl_B5701_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele B5701
ctl_A0301 <- read.csv("ctl_A0301_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele A0301
ctl_B0702 <- read.csv("ctl_B0702_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele B0702

# only work on those peptides that bind to the MHC
bind_criteria <- 50000000 # now we just read in all the peptides
ctl_A0201_binded <- ctl_A0201[ctl_A0201$ic50 < bind_criteria, ]
ctl_B5701_binded <- ctl_B5701[ctl_B5701$ic50 < bind_criteria, ]
ctl_A0301_binded <- ctl_A0301[ctl_A0301$ic50 < bind_criteria, ]
ctl_B0702_binded <- ctl_B0702[ctl_B0702$ic50 < bind_criteria, ]

###############
# we could use a different database of peptides
######
ctl_wt_variant <- read.csv("ctl_wt_variant_immu.csv", stringsAsFactors = FALSE) # this database collapsed all the alleles

#ctl_binded_combinded <- do.call("rbind", list(ctl_A0201_binded, ctl_A0301_binded, ctl_B0702_binded, ctl_B5701_binded))
ctl_binded_combinded <- ctl_wt_variant
ctl_binded_9mer <- ctl_binded_combinded[nchar(ctl_binded_combinded$peptide)==9, ]
ctl_binded_9mer_separated <- ctl_binded_9mer[(ctl_binded_9mer$Immunogenicity > 0.7) | (ctl_binded_9mer$Immunogenicity < 0.3), ]
ctl_binded_9mer_separated$Immunogenicity <- round(ctl_binded_9mer_separated$Immunogenicity)

####################
# define several measures for the amino acid
####################
aa_scales_unnormed <- read.csv("aminoacid_scales.csv", stringsAsFactors = FALSE)
maxs <- apply(aa_scales_unnormed[, 2:4], 2, max)
mins <- apply(aa_scales_unnormed[, 2:4], 2, min)
aa_scales <- data.frame("AminoAcid" = aa_scales_unnormed$AminoAcid, scale(aa_scales_unnormed[, 2:4], center = mins, scale = maxs - mins)) 
score_of_aa <- function(aa, i_metric){
  aa_scales[which(aa_scales$AminoAcid == aa), i_metric]
}


seqs_to_scores <- function(seqs, i_metric){
  t(unname(sapply(seqs, function(seq){
    sapply(strsplit(seq, "")[[1]], score_of_aa, i_metric)
  })))
}

hydro_scores_binded_9mer_separated <- seqs_to_scores(ctl_binded_9mer_separated$peptide, 2)
vol_scores_binded_9mer_separated <- seqs_to_scores(ctl_binded_9mer_separated$peptide, 3)
pol_score_binded_9mer_separated <- seqs_to_scores(ctl_binded_9mer_separated$peptide, 4)

#scores_binded_9mer_separated <- data.frame(cbind(vol_scores_binded_9mer_separated, pol_score_binded_9mer_separated))
scores_binded_9mer_separated <- data.frame(cbind(vol_scores_binded_9mer_separated))
scores_binded_9mer_separated$Immunogenicity <- ctl_binded_9mer_separated$Immunogenicity
##########################
# starts training the model
#########################
library(neuralnet)

get_excluded_weights <- function(n_sites, n_scores, n_hidden_each_score){
  exclude_weigh <- numeric(0)
  for (i in 1:n_scores){
    for (j in 1:n_scores){
      if (j != i){
        for (k in 1:n_hidden_each_score){
          temp <- cbind(rep(1, n_sites), ((i-1)*n_sites+2):(i*n_sites+1), rep((j-1)*n_hidden_each_score+k, n_sites))
          exclude_weigh <- rbind(exclude_weigh, temp)
        }
      }
    }
  }
  exclude_weigh
}

n_scores <- 1
n_hidden_per_score <- 3
excluded <- get_excluded_weights(9, n_scores, n_hidden_per_score)

index_positive <- sample(which(scores_binded_9mer_separated$Immunogenicity > 0.7), round(0.25*nrow(scores_binded_9mer_separated)))
index_negative <- sample(which(scores_binded_9mer_separated$Immunogenicity < 0.3), round(0.25*nrow(scores_binded_9mer_separated)))
index <- c(index_positive, index_negative)
#index <- sample(1:nrow(scores_binded_9mer_separated), round(0.75*nrow(scores_binded_9mer_separated)))
#ind_p_ok <- read.table("positive_index_ok.txt")
#ind_n_ok <- read.table("negative_index_ok.txt")
#index <- c(ind_p_ok$x, ind_n_ok$x)

train_ <- scores_binded_9mer_separated[index, ]
test_ <- scores_binded_9mer_separated[-index, ]
n <- names(train_)
f <- as.formula(paste("Immunogenicity ~", paste(n[!n %in% "Immunogenicity"], collapse = " + ")))
#nn <- neuralnet(f, data = train_, hidden = n_scores*n_hidden_per_score, linear.output = FALSE, exclude = excluded)
nn <- neuralnet(f, data = train_, hidden = n_scores*n_hidden_per_score, linear.output = FALSE)
#nn$weights[[1]][[1]][is.na(nn$weights[[1]][[1]])]  <- 0
pr_train_nn <- predict(nn, train_)
train_accuracy <- data.frame("exper" = train_$Immunogenicity, "pred" = round(pr_train_nn))

pr_test_nn <- predict(nn, test_)
test_accuracy <- data.frame("exper" = test_$Immunogenicity, "pred" = round(pr_test_nn))


#########################
# define a function to calculate the false negative rate and false positive rate 
########################
fp_fn <- function(data_exp_pred){
  fp <- sum((data_exp_pred$exper == 0) & (data_exp_pred$pred == 1))/sum(data_exp_pred$exper == 0)
  fn <- sum((data_exp_pred$exper == 1) & (data_exp_pred$pred == 0))/sum(data_exp_pred$exper == 1)
  list("false_positive_rate" = fp, "false_negative_rate" = fn)
}
fp_fn(train_accuracy) # the false positive rate is very high, even higher than 50%, which means that it can not distinguish the escape mutations. false positive means the negatives that are wrongly characterized as positive
fp_fn(test_accuracy) # of course the false positve rate is even higher for the test set

