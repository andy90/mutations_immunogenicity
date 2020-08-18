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

#hydro_scores_binded_9mer_separated <- seqs_to_scores(ctl_binded_9mer_separated$peptide, 2)
#vol_scores_binded_9mer_separated <- seqs_to_scores(ctl_binded_9mer_separated$peptide, 3)
#pol_score_binded_9mer_separated <- seqs_to_scores(ctl_binded_9mer_separated$peptide, 4)



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

df_scores_immuno <-  data.frame("ImmuType" = "Immuno", "pos" = 1:9, "hydro_mean" = colMeans(hydro_scores_immuno), "hydro_sd" = get_scores_errorbar(hydro_scores_immuno),
                                "vol_mean" = colMeans(vol_scores_immuno), "vol_sd" = get_scores_errorbar(vol_scores_immuno),
                                "pol_mean" = colMeans(pol_scores_immuno), "pol_sd" = get_scores_errorbar(pol_scores_immuno))

df_scores_nonimmuno <-  data.frame("ImmuType" = "Non-Immuno", "pos" = 1:9, "hydro_mean" = colMeans(hydro_scores_nonimmuno), "hydro_sd" = get_scores_errorbar(hydro_scores_nonimmuno),
                                   "vol_mean" = colMeans(vol_scores_nonimmuno), "vol_sd" = get_scores_errorbar(vol_scores_nonimmuno),
                                   "pol_mean" = colMeans(pol_scores_nonimmuno), "pol_sd" = get_scores_errorbar(pol_scores_nonimmuno))

df_scores <- rbind(df_scores_immuno, df_scores_nonimmuno)

ggplot(df_scores,aes(x=as.factor(pos),y=hydro_mean,fill=ImmuType))+
  geom_bar(stat="identity",position="dodge")+
  xlab("Residue position")+ylab("Mean Hydrophobicity")+geom_errorbar(aes(ymin=hydro_mean-hydro_sd, ymax=hydro_mean+hydro_sd), width=.2,position=position_dodge(.9))

ggplot(df_scores,aes(x=as.factor(pos),y=vol_mean,fill=ImmuType))+
  geom_bar(stat="identity",position="dodge")+
  xlab("Residue position")+ylab("Mean Volume")+geom_errorbar(aes(ymin=vol_mean-vol_sd, ymax=vol_mean+vol_sd), width=.2,position=position_dodge(.9))
# the difference in mean volume at position 3 might also be important

ggplot(df_scores,aes(x=as.factor(pos),y=pol_mean,fill=ImmuType))+
  geom_bar(stat="identity",position="dodge")+
  xlab("Residue position")+ylab("Mean Polarity")+geom_errorbar(aes(ymin=pol_mean-pol_sd, ymax=pol_mean+pol_sd), width=.2,position=position_dodge(.9))

# no consistent enrichment of hydrophobicity, volume or polarity is found


#######################
# do the PCA analysis for the data
#######################
df_scores_eachsample_immuno <- data.frame("ImmuType" = "Immuno", hydro_scores_immuno, vol_scores_immuno, pol_scores_immuno)
df_scores_eachsample_nonimmuno <- data.frame("ImmuType" = "Non-Immuno", hydro_scores_nonimmuno, vol_scores_nonimmuno, pol_scores_nonimmuno)
df_scores_eachsample <- rbind(df_scores_eachsample_immuno, df_scores_eachsample_nonimmuno)
pca_hydro_score <-prcomp(df_scores_eachsample[, 2:10], center = TRUE,scale. = TRUE)
summary(pca_hydro_score)

pca_vol_score <- prcomp(df_scores_eachsample[, 11:19], center = TRUE,scale. = TRUE)
summary(pca_vol_score)

pca_pol_score <- prcomp(df_scores_eachsample[, 20:28], center = TRUE,scale. = TRUE)
summary(pca_pol_score)

pca_total_score <- prcomp(df_scores_eachsample[, 2:28], center = TRUE,scale. = TRUE)
summary(pca_total_score)

# we need a lot of principle components to capture the variation


##########################
# do the tsne plot for the data
###########################
library(Rtsne)
tsne <- Rtsne(df_scores_eachsample[, 11:19], dims = 2, check_duplicate = FALSE)
plot(tsne$Y, t = 'n', main = 'tsne')
train_label <- as.factor(df_scores_eachsample$ImmuType)
levels(train_label) <- list(A = "Immuno", B = "Non-Immuno")
colors = rainbow(length(unique(train_label)))
names(colors) = unique(train_label)
text(tsne$Y, labels = train_label, col = colors[train_label])


#######################
# now we will check out the coupling zi*zj, where i and j are two arbitrary positions
########################

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
pos_pair <- {}
for (i in 1:(nsites - 1)){
  for (j in (i+1):nsites){
    pos_pair <- c(pos_pair, paste(i,j,sep = ","))
  }
}

df_coupling_scores_immuno <-  data.frame("ImmuType" = "Immuno", "pos" = pos_pair, "hydro_coup_mean" = colMeans(hydro_coupling_scores_immuno), "hydro_coup_sd" = get_scores_errorbar(hydro_coupling_scores_immuno),
                                "vol_coup_mean" = colMeans(vol_coupling_scores_immuno), "vol_coup_sd" = get_scores_errorbar(vol_coupling_scores_immuno),
                                "pol_coup_mean" = colMeans(pol_coupling_scores_immuno), "pol_coup_sd" = get_scores_errorbar(pol_coupling_scores_immuno))

df_coupling_scores_nonimmuno <-  data.frame("ImmuType" = "Non-Immuno", "pos" = pos_pair, "hydro_coup_mean" = colMeans(hydro_coupling_scores_nonimmuno), "hydro_coup_sd" = get_scores_errorbar(hydro_coupling_scores_nonimmuno),
                                            "vol_coup_mean" = colMeans(vol_coupling_scores_nonimmuno), "vol_coup_sd" = get_scores_errorbar(vol_coupling_scores_nonimmuno),
                                            "pol_coup_mean" = colMeans(pol_coupling_scores_nonimmuno), "pol_coup_sd" = get_scores_errorbar(pol_coupling_scores_nonimmuno))
df_coupling_scores <- rbind(df_coupling_scores_immuno, df_coupling_scores_nonimmuno)

ggplot(df_coupling_scores,aes(x=as.factor(pos),y=hydro_coup_mean,fill=ImmuType))+
  geom_bar(stat="identity",position="dodge")+
  xlab("Residue pair position")+ylab("Mean Coupling Hydrophobicity")+geom_errorbar(aes(ymin=hydro_coup_mean-hydro_coup_sd, ymax=hydro_coup_mean+hydro_coup_sd), width=.2,position=position_dodge(.9))
# 5,8  coupling might be good features to distinguish the immuno and non immuno types

ggplot(df_coupling_scores,aes(x=as.factor(pos),y=vol_coup_mean,fill=ImmuType))+
  geom_bar(stat="identity",position="dodge")+
  xlab("Residue pair position")+ylab("Mean Coupling Volume")+geom_errorbar(aes(ymin=vol_coup_mean-vol_coup_sd, ymax=vol_coup_mean+vol_coup_sd), width=.2,position=position_dodge(.9))
# 1,7  2,3  2,7 3,4  3,5 3,7 5,7 4,7 7,9 might be good features to distinguish the immuno and non immuno types. 1,7 and 3,7 are the best

ggplot(df_coupling_scores,aes(x=as.factor(pos),y=pol_coup_mean,fill=ImmuType))+
  geom_bar(stat="identity",position="dodge")+
  xlab("Residue pair position")+ylab("Mean Coupling Polarity")+geom_errorbar(aes(ymin=pol_coup_mean-pol_coup_sd, ymax=pol_coup_mean+pol_coup_sd), width=.2,position=position_dodge(.9))
# 1,4 6,8 might be good

########################
# look at the probability distribution of the dataset
########################
colnames(hydro_coupling_scores_immuno) <- pos_pair
colnames(hydro_coupling_scores_nonimmuno) <- pos_pair

colnames(vol_coupling_scores_immuno) <- pos_pair
colnames(vol_coupling_scores_nonimmuno) <- pos_pair

colnames(pol_coupling_scores_immuno) <- pos_pair
colnames(pol_coupling_scores_nonimmuno) <- pos_pair
plot(density(vol_coupling_scores_immuno[, "3,7"]))
plot(density(vol_coupling_scores_nonimmuno[, "3,7"]))

plot(density(vol_coupling_scores_immuno[, "1,7"]))
plot(density(vol_coupling_scores_nonimmuno[, "1,7"]))

mean(vol_coupling_scores_immuno[, "3,7"])
mean(vol_coupling_scores_nonimmuno[, "3,7"])

df_coupling_scores_eachsample_immuno <- setNames(data.frame("Immuno", hydro_coupling_scores_immuno, vol_coupling_scores_immuno, pol_coupling_scores_immuno), c("ImmuType", paste("hydro", pos_pair, sep = ","), paste("vol", pos_pair, sep = ","), paste("pol", pos_pair, sep = ",")))
df_coupling_scores_eachsample_nonimmuno <- setNames(data.frame("Non-Immuno", hydro_coupling_scores_nonimmuno, vol_coupling_scores_nonimmuno, pol_coupling_scores_nonimmuno), c("ImmuType", paste("hydro", pos_pair, sep = ","), paste("vol", pos_pair, sep = ","), paste("pol", pos_pair, sep = ",")))
df_coupling_scores_eachsample <- rbind(df_coupling_scores_eachsample_immuno, df_coupling_scores_eachsample_nonimmuno)

#tsne <- Rtsne(df_coupling_scores_eachsample[, c("vol,1,7", "vol,3,7")], dims = 2,check_duplicates = FALSE) 
#tsne <- Rtsne(df_coupling_scores_eachsample[, c("hydro,5,8", "vol,1,7", "vol,3,7", "pol,1,4", "pol,6,8")], dims = 2,check_duplicates = FALSE) 
tsne <- Rtsne(df_coupling_scores_eachsample[, 2:ncol(df_coupling_scores_eachsample)], dims = 2,check_duplicates = FALSE) 
plot(tsne$Y, t = 'n', main = 'tsne')
train_label <- as.factor(df_scores_eachsample$ImmuType)
levels(train_label) <- list(A = "Immuno", B = "Non-Immuno")
colors = rainbow(length(unique(train_label)))
names(colors) = unique(train_label)
text(tsne$Y, labels = train_label, col = colors[train_label]) # the tsne result is not so good, similar to the case when we use all the single probabilities 
