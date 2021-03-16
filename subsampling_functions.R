library(lme4)
library(lmerTest)
library(dplyr)
library(pldist)
library(MiRKAT)
library(vegan)
source("./mp_analysis_functions.R")

##
## SUBSAMPLING AT TIME INTERVALS
## 

subsample <- function(file.path, n.subj, n.time, interval) {
  sim.dataset <- as.matrix(read.table(file = file.path))
  
  n.otu <- nrow(sim.dataset)
  n.interval <- n.time %/% interval
  
  ss.dataset <- matrix(nrow = n.otu, ncol = 0)
  for (i in 1:(n.interval)) {
    sim.start <- ((i-1) * interval) * n.subj + 1
    sim.end <- sim.start + n.subj - 1

    ss.dataset <- cbind(ss.dataset, sim.dataset[ , sim.start:sim.end])
  }
  return(ss.dataset)
}



##
## DETECTING DIFFERENCE IN QUALITATIVE VOLATILITY
##

qualitative.test <- function(file.path, n.subj, direction = "disappearance") {
  sim.dataset <- as.matrix(read.table(file = file.path))
  
  input_1 <- list() # group 1: otutab and subject ID
  input_2 <- list() # group 2
  for (i in 1:n.subj) {
    subj.id <- paste("SUBJ_", i, sep = "")
    sim.otutab_i <- as.data.frame(sim.dataset) %>%
      dplyr::select(contains(paste(subj.id, ".", sep = "")))
    if(i <= (n.subj/2)) {
      input_1[[i]] <- list(sim.otutab_i, subj.id)
    }
    else {
      input_2[[i-(n.subj/2)]] <- list(sim.otutab_i, subj.id)
    } 
  }
  
  if(direction == "disappearance") {
    big.disapp.tab_1 <- disapp_prop_table(input_1) %>%    # table of possible/actual disappearances
      mutate(group2 = 0)                                # indicator variable for membership in group 2
    big.disapp.tab_2 <- disapp_prop_table(input_2) %>%
      mutate(group2 = 1)
    big.disapp.tab <- rbind(big.disapp.tab_1, big.disapp.tab_2)
    
    glmm <- glmer(prop.disapp ~ 1 + group2 + (1 | subj.id) + (1 | otu.id),   # mixed effects logistic regression
                  data = big.disapp.tab,
                  family = binomial,
                  weights = n.presences)
  }
  else {
    big.reapp.tab_1 <- reapp_prop_table(input_1) %>%    # table of actual / possible reappearances
      mutate(group2 = 0)                              # indicator variable for membership in group 2
    big.reapp.tab_2 <- reapp_prop_table(input_2) %>%
      mutate(group2 = 1)
    big.reapp.tab <- rbind(big.reapp.tab_1, big.reapp.tab_2)
    
    glmm <- glmer(prop.reapp ~ 1 + group2 + (1 | subj.id) + (1 | otu.id),   # mixed effects logistic regression
                  data = big.reapp.tab,
                  family = binomial, 
                  weights = n.absences)
  }
  glmm.summary <- summary(glmm)
  return(glmm.summary$coefficients["group2", ])
}

doManyQualTest <- function(n.set, folder.path, n.subj, direction = "disappearance") {
  registerDoParallel(cores = detectCores()-1)
  res <- foreach(i = 1:n.set) %dopar% qualitative.test(file.path = paste(folder.path, "/set", sprintf("%04d", i), ".txt", sep = ""), n.subj, direction)
  res <- matrix(unlist(res), ncol = 4, byrow = TRUE)
  return(res)
}

##
## DETECTING DIFFERENCE IN QUALITATIVE VOLATILITY
##

quantitative.test <- function(file.path, n.subj, n.time, subsample = FALSE, interval = NULL) {
  if (subsample) {
    stopifnot(!is.null(interval))
    sim.dataset <- subsample(file.path, n.subj, n.time, interval)
    n.fc <- n.time %/% interval
  } else {
    sim.dataset <- as.matrix(read.table(file = file.path))
    n.fc <- n.time - 1
  }
  
  big.logFC.tab <- data.frame(logFC = numeric(n.subj * n.fc),
                              group2 = numeric(n.subj * n.fc),
                              subj.id = character(n.subj * n.fc),
                              stringsAsFactors = FALSE)
  
  sim.dataset <- as.matrix(read.table(file = file.path))
  for (i in 1:n.subj) {
    subj.id_i <- paste("SUBJ_", i, sep = "")
    sim.otutab_i <- as.data.frame(sim.dataset) %>%
      dplyr::select(contains(paste(subj.id_i, ".", sep = "")))
    sim.logFC.tab <- fold_difference_table(sim.otutab_i, logRatio = TRUE)
    # within-sample squared L2 norm of log fold-changes excluding reappearances and continued absences
    logFC_i <- apply(sim.logFC.tab, 2, FUN = function(x) sum((x[is.finite(x)])**2, na.rm = TRUE))
    
    index.start <- ((i-1) * n.fc) + 1
    index.end <- i * n.fc
    
    big.logFC.tab[index.start:index.end, 1] <- logFC_i # log fold-change values
    big.logFC.tab[index.start:index.end, 3] <- rep(subj.id_i, n.fc) # subj.id
    if(i > (n.subj / 2)) {
      big.logFC.tab[index.start:index.end, 2] <- rep(1, n.fc) # group2 = 1
    }
  }
  lmm <- lmer(logFC ~ 1 + group2 + (1 | subj.id), data = big.logFC.tab)
  lmm.summary <- summary(lmm)
  return(lmm.summary$coefficients["group2", ])
}

doManyQuantTest <- function(n.set, folder.path, n.subj, n.time) {
  registerDoParallel(cores = detectCores()-1)
  res <- foreach(i = 1:n.set) %dopar% quantitative.test(file.path = paste(folder.path, "/set", sprintf("%04d", i), ".txt", sep = ""), n.subj, n.time)
  res <- matrix(unlist(res), ncol = 5, byrow = TRUE)
  return(res)
}

## 
## Microbiome Regression-based Association Tests
##

# pldist(): returns a NxN (N = number of subjects) distance matrix
# otus: rows are samples, columns are taxa
# metadata: subjID, sampID, time
# paired = FALSE
# binary: qualitative (TRUE) or quantitative (FALSE)
# method: braycurtis, jaccard, kulczynski, gower, unifrac (requires phylo tree)

mirkat.test <- function(file.path, n.subj, n.time, response, method = "bray", subsample = FALSE, interval = NULL) {
  if (subsample) {
    stopifnot(!is.null(interval))
    sim.dataset <- subsample(file.path, n.subj, n.time, interval)
  } else {
    sim.dataset <- as.matrix(read.table(file = file.path))
  }
  
  otus <- t(sim.dataset)
  sampID <- rownames(otus)
  subjID <- sapply(strsplit(sampID, split = "[.]"), "[[", 1)
  time <- sapply(strsplit(sampID, split = "[.]"), "[[", 2)
  metadata <- data.frame(subjID, sampID, time)
  
  Ds <- pldist_all(otu = otus, metadata = metadata, method = c("bray"))
  Ks <- lapply(Ds, FUN = function(d) D2K(d))
  res <- MiRKAT(y = response, Ks = Ks)
  return(unlist(res))
}

doManyMirkatTest <- function(n.set, folder.path, n.subj, n.time, response, method = "bray", subsample = FALSE, interval = NULL) {
  registerDoParallel(cores = detectCores()-1)
  res.list <- foreach(i = 1:n.set) %dopar% mirkat.test(file.path = paste(folder.path, "/set", sprintf("%04d", i), ".txt", sep = ""), n.subj, n.time, response, binary, method)
  res <- matrix(unlist(res.list), nrow = length(res.list), byrow = TRUE)
  return(res)
}

adhoc.test <- function(file.path, n.subj, method = "bray", binary = FALSE, subsample = FALSE, n.time = NULL, interval = NULL) {
  if (subsample) {
    stopifnot(!is.null(interval), !is.null(n.time))
    sim.dataset <- subsample(file.path, n.subj, n.time, interval)
  } else {
    sim.dataset <- as.matrix(read.table(file = file.path))
  }
  
  ## create data dataset with: distance, subj.id, group2
  big.dist.tab <- data.frame(distance = numeric(),
                             group2 = numeric(),
                             subj.id = character(),
                             stringsAsFactors = FALSE)
  for (i in 1:n.subj) {
    subj.id_i <- paste("SUBJ_", i, sep = "")
    sim.otutab_i <- as.data.frame(sim.dataset) %>%
      dplyr::select(contains(paste(subj.id_i, ".", sep = "")))
    
    n.pairs <- ncol(sim.otutab_i) - 1
    
    if (i <= (n.subj / 2)) {
      group2 <- rep(0, n.pairs)
    } else {
      group2 <- rep(1, n.pairs)
    }
    
    dist.mat <- as.matrix(vegdist(x = t(sim.otutab_i), method = method, binary = binary))
    
    dist.tab_i <- data.frame(distance = diag(dist.mat[2:nrow(dist.mat), 1:(nrow(dist.mat)-1)]),
                             group2 = group2,
                             subj.id = rep(subj.id_i, n.pairs))
    ## rbind with data dataset
    big.dist.tab <- rbind(big.dist.tab, dist.tab_i)
  }
  
  # LMM: distance ~ 1 + group2 + (1|subj.id)
  lmm <- lmer(distance ~ 1 + group2 + (1 | subj.id), data = big.dist.tab)
  lmm.summary <- summary(lmm)
  return(lmm.summary$coefficients["group2", ])
}








