rm(list=ls()) # clean the memeory
library(dendextend)
set.seed(1)
rstr <- function(n,k){   # vector of n random char(k) strings
  sapply(1:n,function(i){do.call(paste0,as.list(sample(letters,k,replace=T)))})
}

str<- c(paste0("aa",rstr(10,3)),paste0("bb",rstr(10,3)),paste0("cc",rstr(10,3)))
# Levenshtein Distance
d  <- adist(str)
rownames(d) <- str
hc <- hclust(as.dist(d))
hc_tocolor <- as.dendrogram(hc) # this creates a new object, which could be used to color the different labels
colors_to_use <- c(rep(1,15), rep(2,15)) # the labels that corresponds the the str
colors_to_use <- colors_to_use[order.dendrogram(hc_tocolor)] # reorder the color labels to match the order in the dendrogram
labels_colors(hc_tocolor)  <- colors_to_use # assign the colors
plot(hc_tocolor) # plot
rect.hclust(hc_tocolor,k=2)
df <- data.frame(str,cutree(hc,k=3))


rm(list=ls()) # clean the memeory
######################################
# above are some sample code to do the clustering plot
#####################################
ctl_wt_variant <- read.csv("ctl_wt_variant_immu.csv", stringsAsFactors = FALSE) # this database collapsed all the alleles

ctl_9mer <- ctl_wt_variant[nchar(ctl_wt_variant$peptide)==9, ]
ctl_9mer_immuno <- ctl_9mer[(ctl_9mer$Immunogenicity > 0.7) ,]
ctl_9mer_nonimmuno <- ctl_9mer[(ctl_9mer$Immunogenicity < 0.3) ,]

ctl_9mer_immuno$Immunogenicity <- round(ctl_9mer_immuno$Immunogenicity)
ctl_9mer_nonimmuno$Immunogenicity <- round(ctl_9mer_nonimmuno$Immunogenicity)
ctl_9mer_combined <- rbind(ctl_9mer_immuno[, 2:3], ctl_9mer_nonimmuno[, 2:3])

d <- adist(ctl_9mer_combined$peptide)
rownames(d) <- ctl_9mer_combined$peptide
hc <- hclust(as.dist(d))
hc_tocolor <- as.dendrogram(hc) # this creates a new object, which could be used to color the different labels
colors_to_use <- ctl_9mer_combined$Immunogenicity
colors_to_use <- colors_to_use[order.dendrogram(hc_tocolor)] # reorder the color labels to match the order in the dendrogram
labels_colors(hc_tocolor)  <- colors_to_use # assign the colors
pdf("hclust.pdf", width=40, height=15)
plot(hc_tocolor) # plot
dev.off()

library(ape)
plot_pcoa <- pcoa(d)
biplot(plot_pcoa)
####################
# now try to use much less data
#####################
rm(list=ls()) # clean the memeory
ctl_A0201 <- read.csv("ctl_A0201_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele A0201
ctl_B5701 <- read.csv("ctl_B5701_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele B5701
ctl_A0301 <- read.csv("ctl_A0301_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele A0301
ctl_B0702 <- read.csv("ctl_B0702_SeqImmuIC50.csv", stringsAsFactors = FALSE) # read in allele B0702

# only work on those peptides that bind to the MHC
bind_criteria <- 500 # now we just read in all the peptides
ctl_A0201_binded <- ctl_A0201[ctl_A0201$ic50 < bind_criteria, ]
ctl_B5701_binded <- ctl_B5701[ctl_B5701$ic50 < bind_criteria, ]
ctl_A0301_binded <- ctl_A0301[ctl_A0301$ic50 < bind_criteria, ]
ctl_B0702_binded <- ctl_B0702[ctl_B0702$ic50 < bind_criteria, ]

