rm(list=ls()) # clean the memeory
n_letters <- 4 # number of coarse grained groups
n_sites <- 9 # number of sites

library("matrixcalc")
library("lattice")


index_map <- {}
for (i in 1:(n_sites - 1)){
  for (j in (i+1):n_sites){
    line_num <- (i-1)*(n_sites-2) - (i-1)*(i-2)/2 + j - 1
    index_map <- rbind(index_map, c(i, j, line_num))
  }
}  

J_raw_immuno <- data.matrix(read.table("J_ZS.txt")) # this is the J matrix obtained from John's method, we need to fill it into the J we want

getJ <- function(J_reshaped){
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

J_immuno <- getJ(J_raw_immuno)
image(t(J_immuno[nrow(J_immuno):1, ]), zlim = c(-1, 1))
image(abs(t(J_immuno[nrow(J_immuno):1, ])), zlim = c(0, 1)) # this is looking at the strength of coupling

J_raw_nonimmuno <- data.matrix(read.table("J_ZS_noimmuno.txt"))
J_nonimmuno <- getJ(J_raw_nonimmuno)
image(abs(t(J_nonimmuno[nrow(J_nonimmuno):1, ])), zlim = c(0, 1)) # this is looking at the strength of coupling

levelplot(abs(t(J_immuno[nrow(J_immuno):1, ])), asp = "iso", col.regions = heat.colors, xlab = NULL, ylab = NULL)
levelplot(abs(t(J_nonimmuno[nrow(J_nonimmuno):1, ])), asp = "iso", col.regions = heat.colors, xlab = NULL, ylab = NULL)

#filled.contour(abs(t(J_immuno[nrow(J_immuno):1, ])), asp = "iso", color = heat.colors)
#heatmap(abs(t(J_immuno[nrow(J_immuno):1, ])))

get_coarseJ <- function(J_reshaped){
  J <- matrix(0, nrow = n_sites, ncol = n_sites)
  for (line_num in 1:nrow(J_reshaped)){
    ind <- which(line_num == index_map[ ,3])
    i <- index_map[ind, 1]
    j <- index_map[ind, 2]
    
    J[i, j] <- sum(abs(J_reshaped[line_num,]))
    J[j, i] <- sum(abs(J_reshaped[line_num,]))
  }
  J
}

J_immu_coarse <- get_coarseJ(J_raw_immuno)
J_nonimmu_coarse <- get_coarseJ(J_raw_nonimmuno)
levelplot(t(J_immu_coarse[nrow(J_immu_coarse):1,]), asp = "iso", col.regions = heat.colors, xlab = NULL, ylab = NULL)
levelplot(t(J_nonimmu_coarse[nrow(J_nonimmu_coarse):1,]), asp = "iso", col.regions = heat.colors, xlab = NULL, ylab = NULL)
