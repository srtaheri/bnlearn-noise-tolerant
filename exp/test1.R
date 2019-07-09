library(igraph)
library(bnlearn)
library(stringr)
library(dplyr)
library(parallel)
addDiffNoiseLevelsToGraphData <- function(data, noiseType = 'gaussian',
                                          noiseLevel = seq(0,1,1/(ncol(data)-1))) {
  #noiseLevel is the standard deviation of the noise of variables
  stopifnot(noiseType == 'gaussian')
  for (ci in 1:ncol(data)) {
    if(noiseLevel[ci] != 0) {
      data[, ci] = data[, ci] + rnorm(nrow(data), 0, noiseLevel[ci])
    }
  }
  return(data)
}
ham_dist_btw_matrices <- function(matrix_one, matrix_two) {
  if(nrow(matrix_one) != nrow(matrix_two)) {
    stop("The number of rows should be the same.")
  }
  if(ncol(matrix_one) != ncol(matrix_two)) {
    stop("The number of columns should be the same.")
  }
  sum = 0
  for(i in 1:nrow(matrix_one)) {
    for(j in 1:ncol(matrix_one)) {
      sum = sum + abs(matrix_one[i,j] - matrix_two[i,j])
    }
  }
  return(sum)
}
###################### Experiment 1 - Tree with 16 nodes ###################################################
# Add different noise levels to nodes in the network
set.seed(1)
setwd("C:/Users/Ehsan/Documents/GitHub/bnlearn_forked/exp/")
data_t16 <- readRDS("data_t16.RData")
g_t16 <- readRDS("g_t16.RData")
mat_true <- as.matrix(as_adj(g_t16)) + t(as.matrix(as_adj(g_t16))) # true network adjacency matrix
num_iterations = 20
t16_noise_allNodes <- data.frame("beforeCancel" = rep(-1,num_iterations),
                                 "afterCancel" = rep(-1,num_iterations))
for (iteration in 1:num_iterations) {
  noiseLevels <- sample(seq(0,0.9,0.1),ncol(data_t16),replace = T)
  gData_noisy <- addDiffNoiseLevelsToGraphData(data_t16,noiseLevel = noiseLevels)
  noiseLevels <- data.frame(t(noiseLevels))
  colnames(noiseLevels) <- colnames(gData_noisy)
  learned_noisy <- si.hiton.pc(data.frame(gData_noisy))
  adj_mat_learned_noisy <- amat(learned_noisy) + t(amat(learned_noisy))
  adj_mat_learned_noisy[adj_mat_learned_noisy == 2] <- 1
  t16_noise_allNodes[iteration,"beforeCancel"] <- ham_dist_btw_matrices(matrix_one = mat_true,
                                                                        matrix_two = adj_mat_learned_noisy)/2
  # Cancel the noise
  learned_noise_cancel <- si.hiton.pc(data.frame(gData_noisy), noise.levels = noiseLevels)
  adj_mat_noise_cancel <- amat(learned_noise_cancel) + t(amat(learned_noise_cancel))
  adj_mat_noise_cancel[adj_mat_noise_cancel == 2] <- 1
  t16_noise_allNodes[iteration,"afterCancel"] <- ham_dist_btw_matrices(matrix_one = mat_true,
                                                                       matrix_two = adj_mat_noise_cancel)/2
}

print(t16_noise_allNodes)
################### Experiment 2 - graph with 16 nodes#############################################

###################Experiment 3 - CHDI data with around 300 nodes#############################################
CHDI <- read.csv("prunedData.csv")
#only get the first half of the data which is the gene expression data
CHDI_gexpr <- CHDI[,which(vapply(names(CHDI),function(x) grepl("_gexpr",x),FUN.VALUE = logical(1)) == TRUE)]
#only get the second half of the data which is the protein
CHDI_prot <- CHDI[,which(vapply(names(CHDI),function(x) grepl("_prot",x),FUN.VALUE = logical(1)) == TRUE)]

### learn the structure of network from all the available variables
cl = makeCluster(3, type = "SOCK")
pdag = mmpc(CHDI[,seq(4,ncol(CHDI))], cluster = cl)
TruePositives <- c()

for (geneVar in colnames(CHDI_gexpr)) {
  protVar <- str_replace(geneVar,"_gexpr","_prot")
  if(protVar %in% ((pdag$nodes)[geneVar][[1]][1][[1]]))
    TruePositives <- c(TruePositives, c(geneVar,protVar))
}

#calculating noise level of each gene
#we have 147 gene expreseeion
sigma_gex <- list()
for (mid in 1:12) {
  new_data <- CHDI[which(CHDI[,"MID"] == mid),]
  gex_pridx <- grep("gexpr",colnames(new_data))
  gex_data <- new_data[,gex_pridx]
  sigma_gex[[mid]] <- apply(gex_data,2,sd)
}

#compute the pooled average for each gene expression
sum = 0
sum_length_mid= 0
geneVar_noise_level <- c()
for (geneVar in colnames(CHDI_gexpr)) {
  for (mid in 1:12) {
    len = length(which(CHDI[,"MID"] == mid))-1
    sum <- sum + sigma_gex[[mid]][[geneVar]]*(len)
    sum_length_mid <- sum_length_mid + len 
  }
  geneVar_noise_level <- c(geneVar_noise_level,sum/sum_length_mid)
  sum = 0
  sum_length_mid = 0
}
# create noise levels
geneVar_noise_level <- data.frame(t(rbind(data.frame("nLevel" = geneVar_noise_level),
                               data.frame("nLevel" = rep(0,147)))))
colnames(geneVar_noise_level) <- colnames(CHDI)[4:ncol(CHDI)]

data <- CHDI[,seq(4,ncol(CHDI))]
start_time <- Sys.time()
pdag_n_cancel <- mmpc(data,
                      noise.levels = geneVar_noise_level,
                      cluster = cl 
                           )
end_time <- Sys.time()
end_time - start_time

TruePositives_n_cancelled <- c()

for (geneVar in colnames(CHDI_gexpr)) {
  protVar <- str_replace(geneVar,"_gexpr","_prot")
  if(protVar %in% ((pdag_n_cancel$nodes)[geneVar][[1]][1][[1]]))
    TruePositives_n_cancelled <- c(TruePositives_n_cancelled, c(geneVar,protVar))
}
