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
aa_scales <- read.csv("aminoacid_mutscales.csv", stringsAsFactors = FALSE)
#maxs <- apply(aa_scales_unnormed[, 2], 2, max)
#mins <- apply(aa_scales_unnormed[, 2], 2, min)
#aa_scales <- data.frame("AminoAcid" = aa_scales_unnormed$AminoAcid, scale(aa_scales_unnormed[, 2:4], center = (maxs + mins)/2, scale = (maxs - mins)/2))
#aa_scales <- data.frame("AminoAcid" = aa_scales_unnormed$AminoAcid, scale(aa_scales_unnormed[, 2:4], center = mins, scale = (maxs - mins)))
score_of_aa <- function(aa, i_metric){
  aa_scales[which(aa_scales$AminoAcid == aa), i_metric]
}


seqs_to_scores <- function(seqs, i_metric){
  t(unname(sapply(seqs, function(seq){
    sapply(strsplit(seq, "")[[1]], score_of_aa, i_metric)
  })))
}

#########
# define the function to get the frequency of each amino acid in the dataset
#########

get_AAfreq <- function(seqs, nsites){ # read in a vector if peptides, nsites is the length of the peptides
  sapply(aa_scales$AminoAcid, function(aa){
    sum(sapply(seqs, function(seq){
      sum(strsplit(seq, "")[[1]] == aa)
    }))
  })/(length(seqs)*nsites)
}

get_AAfreq_errorbar <- function(seqs, nsites){
  seqs <- sample(seqs)
  AA_freq_subgroups <- numeric(0)
  n = 4
  group_size <- as.integer(ceiling(length(seqs)/n))
  seqs_subs <- split(seqs, ceiling(seq_along(seqs)/group_size))
  for (iseqs in seqs_subs){
    AA_freq_subgroups <- rbind(AA_freq_subgroups, get_AAfreq(iseqs, nsites))
  }
  apply(AA_freq_subgroups, 2, sd)
}

aa_freq_immuno <- get_AAfreq(ctl_9mer_immuno$peptide, 9)
aa_freq_immuno_eb <- get_AAfreq_errorbar(ctl_9mer_immuno$peptide, 9)
aa_freq_nonimmuno <- get_AAfreq(ctl_9mer_nonimmuno$peptide, 9)
aa_freq_nonimmuno_eb <- get_AAfreq_errorbar(ctl_9mer_nonimmuno$peptide, 9)

df_aa_freq_immuno <- data.frame("ImmuType" = "Immuno", "AminoAcid" = aa_scales$AminoAcid, "mean"= aa_freq_immuno, "sd" = aa_freq_immuno_eb)
df_aa_freq_nonimmuno <- data.frame("ImmuType" = "Non-Immuno", "AminoAcid" = aa_scales$AminoAcid, "mean"= aa_freq_nonimmuno, "sd" = aa_freq_nonimmuno_eb)
df_aa_freq <- rbind(df_aa_freq_immuno, df_aa_freq_nonimmuno)

library(ggplot2)
ggplot(df_aa_freq,aes(x=AminoAcid,y=mean,fill=ImmuType))+
  geom_bar(stat="identity",position="dodge")+
  xlab("Amino Acid")+ylab("Mean Percentage")+geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9))

##################
# the next step is to compare the hydrophobicity scores at each position
##################

mut_scores_immuno <- seqs_to_scores(ctl_9mer_immuno$peptide, 2)
mut_scores_nonimmuno <- seqs_to_scores(ctl_9mer_nonimmuno$peptide, 2)

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

df_scores_immuno <-  data.frame("ImmuType" = "Immuno", "pos" = 1:9, "mut_mean" = colMeans(mut_scores_immuno), "mut_sd" = get_scores_errorbar(mut_scores_immuno))

df_scores_nonimmuno <-  data.frame("ImmuType" = "Non-Immuno", "pos" = 1:9, "mut_mean" = colMeans(mut_scores_nonimmuno), "mut_sd" = get_scores_errorbar(mut_scores_nonimmuno))

df_scores <- rbind(df_scores_immuno, df_scores_nonimmuno)

ggplot(df_scores,aes(x=as.factor(pos),y=mut_mean,fill=ImmuType))+
  geom_bar(stat="identity",position="dodge")+
  xlab("Residue position")+ylab("Mean mutphobicity")+geom_errorbar(aes(ymin=mut_mean-mut_sd, ymax=mut_mean+mut_sd), width=.2,position=position_dodge(.9))

#######################
# do the PCA analysis for the data
#######################
df_scores_eachsample_immuno <- data.frame("ImmuType" = "Immuno", mut_scores_immuno)
df_scores_eachsample_nonimmuno <- data.frame("ImmuType" = "Non-Immuno", mut_scores_nonimmuno)
df_scores_eachsample <- rbind(df_scores_eachsample_immuno, df_scores_eachsample_nonimmuno)
pca_mut_score <-prcomp(df_scores_eachsample[, 2:10], center = TRUE,scale. = TRUE)
summary(pca_mut_score)



# we need a lot of principle components to capture the variation


##########################
# do the tsne plot for the data
###########################
library(Rtsne)
tsne <- Rtsne(df_scores_eachsample[, 2:10], dims = 2,check_duplicates = FALSE)
plot(tsne$Y, t = 'n', main = 'tsne')
train_label <- as.factor(df_scores_eachsample$ImmuType)
levels(train_label) <- list(A = "Immuno", B = "Non-Immuno")
colors = rainbow(length(unique(train_label)))
names(colors) = unique(train_label)
text(tsne$Y, labels = train_label, col = colors[train_label])
