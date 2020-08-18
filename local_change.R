rm(list=ls()) # clean the memeory


#####
# read in the data and separate the immuno and non immuno group
#####
ctl_wt_variant <- read.csv("ctl_wt_variant_immu.csv", stringsAsFactors = FALSE) # this database collapsed all the alleles

ctl_9mer <- ctl_wt_variant[nchar(ctl_wt_variant$peptide)==9, ]
ctl_9mer_immuno <- ctl_9mer[(ctl_9mer$Immunogenicity > 0.7) ,]
ctl_9mer_nonimmuno <- ctl_9mer[(ctl_9mer$Immunogenicity < 0.3) ,]


#####
# read in the amino acid scales
#####
aa_scales_unnormed <- read.csv("aminoacid_scales.csv", stringsAsFactors = FALSE)
maxs <- apply(aa_scales_unnormed[, 2:4], 2, max)
mins <- apply(aa_scales_unnormed[, 2:4], 2, min)
#aa_scales <- data.frame("AminoAcid" = aa_scales_unnormed$AminoAcid, scale(aa_scales_unnormed[, 2:4], center = (maxs + mins)/2, scale = (maxs - mins)/2))
aa_scales <- data.frame("AminoAcid" = aa_scales_unnormed$AminoAcid, scale(aa_scales_unnormed[, 2:4], center = mins, scale = (maxs - mins)))
score_of_aa <- function(aa, i_metric){
  aa_scales[which(aa_scales$AminoAcid == aa), i_metric]
}


seqs_to_scores <- function(seqs, i_metric){
  t(unname(sapply(seqs, function(seq){
    sapply(strsplit(seq, "")[[1]], score_of_aa, i_metric)
  })))
}


##################
# the next step is to transform each amino acid to a score, get both the single score and coupling score
##################
hydro_scores_immuno <- seqs_to_scores(ctl_9mer_immuno$peptide, 2)
vol_scores_immuno <- seqs_to_scores(ctl_9mer_immuno$peptide, 3)
pol_scores_immuno <- seqs_to_scores(ctl_9mer_immuno$peptide, 4)

hydro_scores_nonimmuno <- seqs_to_scores(ctl_9mer_nonimmuno$peptide, 2)
vol_scores_nonimmuno <- seqs_to_scores(ctl_9mer_nonimmuno$peptide, 3)
pol_scores_nonimmuno <- seqs_to_scores(ctl_9mer_nonimmuno$peptide, 4)

get_scores_errorbar <- function(scores){
  scores <- scores[sample(1:nrow(scores)), ] # reshuffle the rows
  scores_subgroups <-  numeric(0)
  n = 4
  group_size <- as.integer(ceiling(nrow(scores)/n))
  sub_indexes <- ceiling((1:nrow(scores))/group_size)
  for (i in 1:n){
    scores_subgroups <- rbind(scores_subgroups, colMeans(scores[i == sub_indexes, ]))
  }
  apply(scores_subgroups, 2, sd)
}

get_coupling_scores <- function(scores){ # read in single scores, return the coupling scores
  nsites <- ncol(scores)
  nsamples <- nrow(scores)
  
  coupling <- matrix(rep(0, nsamples*nsites*(nsites-1)/2), nrow = nsamples)
  counter <- 1
  for (i in 1:(nsites - 1)){
    for (j in (i+1):nsites){
      coupling[, counter] <- scores[, i] * scores[, j]
      counter <- counter + 1
    }
  }
  coupling
} 

hydro_coupling_scores_immuno <- get_coupling_scores(hydro_scores_immuno)
hydro_coupling_scores_nonimmuno <- get_coupling_scores(hydro_scores_nonimmuno)

vol_coupling_scores_immuno <- get_coupling_scores(vol_scores_immuno)
vol_coupling_scores_nonimmuno <- get_coupling_scores(vol_scores_nonimmuno)

pol_coupling_scores_immuno <- get_coupling_scores(pol_scores_immuno)
pol_coupling_scores_nonimmuno <- get_coupling_scores(pol_scores_nonimmuno)


nsites <- 9
pos <- {}
for (i in 1:nsites){
  pos <- c(pos, as.character(i))
}
pos_pair <- {}
for (i in 1:(nsites - 1)){
  for (j in (i+1):nsites){
    pos_pair <- c(pos_pair, paste(i,j,sep = "_"))
  }
}

df_pep_scores_immuno <- setNames(data.frame(ctl_9mer_immuno$peptide, 1, vol_scores_immuno, vol_coupling_scores_immuno), c("peptide", "ImmuType", paste("vol", pos, sep = "_"), paste("vol", pos_pair, sep = "_")))
df_pep_scores_nonimmuno <- setNames(data.frame(ctl_9mer_nonimmuno$peptide, 0, vol_scores_nonimmuno, vol_coupling_scores_nonimmuno), c("peptide", "ImmuType", paste("vol", pos, sep = "_"), paste("vol", pos_pair, sep = "_")))

df_pep_scores <- rbind(df_pep_scores_immuno, df_pep_scores_nonimmuno)
#####################
# now I construct the same database for the A list only
######################
all_content <- readLines("optimal_ctl_summary.csv") # this is the los alamos database, which contains HIV epitopes and its variants
skip_second <- all_content[-2] # remove the date line
ctl_wt <- read.csv(textConnection(skip_second), stringsAsFactors = FALSE) # it's important to not take strings as vectors
ctl_wt_immuscore <- data.frame("peptide" = toupper(ctl_wt$Epitope), "Immunogenicity" = rep(1, nrow(ctl_wt)), stringsAsFactors = FALSE)

ctl_9mer_wt <- ctl_wt_immuscore[nchar(ctl_wt_immuscore$peptide)==9, ]
hydro_scores_wt <- seqs_to_scores(ctl_9mer_wt$peptide, 2)
vol_scores_wt <- seqs_to_scores(ctl_9mer_wt$peptide, 3)
pol_scores_wt <- seqs_to_scores(ctl_9mer_wt$peptide, 4)

hydro_coupling_scores_wt <- get_coupling_scores(hydro_scores_wt)
vol_coupling_scores_wt <- get_coupling_scores(vol_scores_wt)
pol_coupling_scores_wt <- get_coupling_scores(pol_scores_wt)

df_pep_scores_wt <- setNames(data.frame(ctl_9mer_wt$peptide, "WT", vol_scores_wt, vol_coupling_scores_wt), c("peptide", "ImmuType", paste("vol", pos, sep = "_"), paste("vol", pos_pair, sep = "_")))


##############
#now I define a function to compare the difference between peptides
###############
compare_pep <- function(pep1, pep2){
  pep1_array <- strsplit(as.character(pep1), "")[[1]]
  pep2_array <- strsplit(as.character(pep2), "")[[1]]
  sum(pep1_array != pep2_array)
}

m_diff <- sapply(df_pep_scores$peptide, function(pep1){
  which.min(sapply(df_pep_scores_wt$peptide, function(pep2){
    compare_pep(pep1, pep2)
  }))
})

n_diff <- sapply(df_pep_scores$peptide, function(pep1){
  min(sapply(df_pep_scores_wt$peptide, function(pep2){
    compare_pep(pep1, pep2)
  }))
})

diff_score <- data.frame()
for (i in 1:nrow(df_pep_scores)){
  diff_temp <- df_pep_scores[i, 3:47]-df_pep_scores_wt[m_diff[i], 3:47]
  #print(class(as.vector(diff_temp)))
  diff_score <- rbind(diff_score, diff_temp)
}

diff_score$peptide <- df_pep_scores$peptide
diff_score$ImmuType <- df_pep_scores$ImmuType
diff_score$ndiff <- n_diff


################
# now we do the Tsne plot for the data
################

tsne <- Rtsne(diff_score[, 1:45], dims = 2,check_duplicates = FALSE) 
plot(tsne$Y, t = 'n', main = 'tsne')
train_label <- as.factor(diff_score$ndiff)
#train_label <- as.factor(diff_score$ImmuType)
#levels(train_label) <- list(A = "Immuno", B = "NonImmuno")
colors = rainbow(length(unique(train_label)))
names(colors) = unique(train_label)
text(tsne$Y, labels = train_label, col = colors[train_label]) # the tsne result is not so good, similar to the case when we use all the single probabilities 

############
# do logistic regression for the diff score data
###########
library(neuralnet)
diff_score_nonzero <- diff_score[(diff_score$ndiff > 0) & (diff_score$ndiff < 3), ]
index <- sample(1:nrow(diff_score_nonzero), round(0.5*nrow(diff_score_nonzero)))
diff_score_nonzero <- diff_score_nonzero[, c(1:9, 47)]
train_ <- diff_score_nonzero[index, ]
test_ <- diff_score_nonzero[-index, ]
n <- names(train_)
f <- as.formula(paste("ImmuType ~", paste(n[!n %in% c("ImmuType", "ndiff", "peptide")], collapse = " + ")))
#nn <- neuralnet(f, data = train_, hidden = n_scores*n_hidden_per_score, linear.output = FALSE, exclude = excluded)
nn <- neuralnet(f, data = train_, hidden = 1, linear.output = FALSE)
#nn$weights[[1]][[1]][is.na(nn$weights[[1]][[1]])]  <- 0
pr_train_nn <- predict(nn, train_)
train_accuracy <- data.frame("exper" = train_$ImmuType, "pred" = round(pr_train_nn))

pr_test_nn <- predict(nn, test_)
test_accuracy <- data.frame("exper" = test_$ImmuType, "pred" = round(pr_test_nn))


#########################
# define a function to calculate the false negative rate and false positive rate 
########################
fp_fn <- function(data_exp_pred){
  fp <- sum((data_exp_pred$exper == 0) & (data_exp_pred$pred == 1))/sum(data_exp_pred$exper == 0)
  fn <- sum((data_exp_pred$exper == 1) & (data_exp_pred$pred == 0))/sum(data_exp_pred$exper == 1)
  list("false_positive_rate" = fp, "false_negative_rate" = fn)
}
fp_fn(train_accuracy) # the false positive rate is very high, even higher than 50%, which means that it can not distinguish the escape mutations
fp_fn(test_accuracy) # of course the false positve rate is even higher for the test set
