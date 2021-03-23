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

subsample <- function(file.path, n.subj, n.time, n.interval, sub.subj = n.subj / 2, rand.interval = FALSE) {
  # assumed: samples are columns S1.1 S2.1 ... SN.1 ... S1.T S2.T ... SN.T
  # assumed: first n.subj / 2 subjects are in group 1, the rest are in group 2
  sim.dataset <- as.matrix(read.table(file = file.path)) # read simulated OTU table
  n.otu <- nrow(sim.dataset) # taxa are rows
  ss.dataset <- matrix(nrow = n.otu, ncol = 0) # subsampled data
  
  if (rand.interval) {
    # extract first (sub)group
    for (i in 1:sub.subj) {
      t <- sort(sample(x = n.time, size = n.interval)) # random time points
      samp.ind <- (t-1) * n.subj + i # indices for this subject and time points in sim.dataset
      ss.dataset <- cbind(ss.dataset, sim.dataset[ ,samp.ind]) # bind to result
    }
    # extract second (sub)group
    group2.start <- (n.subj / 2) + 1 # first subjID in group2
    group2.end <- (n.subj / 2) + sub.subj # last subjID in (sub)group2
    for (i in group2.start:group2.end) {
      t <- sort(sample(x = n.time, size = n.interval)) # random time points
      samp.ind <- (t-1) * n.subj + i # indices for this subject and time points in sim.dataset
      ss.dataset <- cbind(ss.dataset, sim.dataset[ , samp.ind]) # bind to result
    }
  } else {
    interval <- n.time %/% n.interval # length of evenly spaced intervals
    for (i in 1:(n.interval)) {
      # select subset from group1
      sim.start <- ((i-1) * interval) * n.subj + 1
      sim.end <- sim.start + sub.subj - 1
      ss.dataset <- cbind(ss.dataset, sim.dataset[ , sim.start:sim.end])
      # select subset from group2
      sim.start <- ((i-1) * interval) * n.subj + (n.subj/2) + 1
      sim.end <- sim.start + sub.subj - 1
      ss.dataset <- cbind(ss.dataset, sim.dataset[ , sim.start:sim.end])
    }
  }
  
  return(ss.dataset)
}

##
## DETECTING DIFFERENCE IN QUALITATIVE VOLATILITY
##

qualitative.test <- function(file.path, n.subj, n.time, n.interval = n.time, sub.subj = n.subj / 2, rand.interval = FALSE, direction = "both") {
  sim.dataset <- subsample(file.path, n.subj, n.time, n.interval, sub.subj, rand.interval)   # subsampled OTU tab
  otu.avg.abnd <- apply(sim.dataset, 1, mean)       # global OTU average abundances
  otu.quintiles <- data.frame(otu.avg.abnd) %>%     # abundance-based OTU quintiles
    rownames_to_column(var = "otu.id") %>%
    mutate(quintile = ntile(otu.avg.abnd, n = 5)) 
  
  n.totalsubj <- sub.subj * 2 # total number of subsampled subjects
  subj.ind <- c(1:sub.subj, (n.subj/2 + 1):(n.subj/2+sub.subj)) # subjIDs
  
  input_1 <- list() # group 1: otutab and subject ID
  input_2 <- list() # group 2
  for (i in 1:n.totalsubj) {
    subj.id <- paste("SUBJ_", subj.ind[i], sep = "")
    sim.otutab_i <- as.data.frame(sim.dataset) %>%
      dplyr::select(contains(paste(subj.id, ".", sep = "")))
    if(i <= sub.subj) {
      input_1[[i]] <- list(sim.otutab_i, subj.id)
    }
    else {
      input_2[[i-sub.subj]] <- list(sim.otutab_i, subj.id)
    } 
  }
  
  if(direction == "both") {
    big.statechange.tab_1 <- statechange_prop_table(input_1) %>%
      mutate(group2 = 0)
    
    big.statechange.tab_2 <- statechange_prop_table(input_2) %>%
      mutate(group2 = 1)
    
    big.statechange.tab <- rbind(big.statechange.tab_1, big.statechange.tab_2)
    big.statechange.tab <- left_join(x = big.statechange.tab, y = otu.quintiles, by = "otu.id")
    big.statechange.tab$quintile <- as.factor(big.statechange.tab$quintile)
    
    glmm <- glmer(prop.statechange ~ 1 + group2 + quintile + (1 | subj.id) + (1 | otu.id),
                  data = big.statechange.tab, 
                  family = binomial, 
                  weights = timepoints)
    
  } else if(direction == "disappearance") {
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

doManyQualTest <- function(n.set, folder.path, n.subj, n.time, n.interval = n.time, sub.subj = n.subj / 2, rand.interval = FALSE, direction = "both") {
  registerDoParallel(cores = detectCores()-1)
  res <- foreach(i = 1:n.set) %dopar% qualitative.test(file.path = paste(folder.path, "/set", sprintf("%04d", i), ".txt", sep = ""), 
                                                       n.subj, n.time, n.interval, sub.subj, rand.interval, direction)
  res <- matrix(unlist(res), nrow = length(res), byrow = TRUE)
  return(res)
}

##
## DETECTING DIFFERENCE IN QUALITATIVE VOLATILITY
##

quantitative.test <- function(file.path, n.subj, n.time, n.interval = n.time, sub.subj = n.subj / 2, rand.interval = FALSE) {
  sim.dataset <- subsample(file.path, n.subj, n.time, n.interval, sub.subj, rand.interval)
  n.fc <- n.interval - 1 # number of fold-changes
  n.totalsubj <- sub.subj * 2 # total number of subsampled subjects
  big.logFC.tab <- data.frame(logFC = numeric(n.totalsubj * n.fc),
                              group2 = numeric(n.totalsubj * n.fc),
                              subj.id = character(n.totalsubj * n.fc),
                              stringsAsFactors = FALSE)
  subj.ind <- c(1:sub.subj, (n.subj/2 + 1):(n.subj/2+sub.subj)) # subjIDs
  for (i in 1:n.totalsubj) {
    subj.id_i <- paste("SUBJ_", subj.ind[i], sep = "")
    sim.otutab_i <- as.data.frame(sim.dataset) %>%
      dplyr::select(contains(paste(subj.id_i, ".", sep = "")))
    sim.logFC.tab <- fold_difference_table(sim.otutab_i, logRatio = TRUE)
    # within-sample squared L2 norm of log fold-changes excluding dis/reappearances and continued absences
    logFC_i <- apply(sim.logFC.tab, 2, FUN = function(x) mean((x[is.finite(x)])**2, na.rm = TRUE))
    
    index.start <- ((i-1) * n.fc) + 1
    index.end <- i * n.fc
    
    big.logFC.tab[index.start:index.end, 1] <- logFC_i # log fold-change values
    big.logFC.tab[index.start:index.end, 3] <- rep(subj.id_i, n.fc) # subj.id
    if(i > sub.subj) {
      big.logFC.tab[index.start:index.end, 2] <- rep(1, n.fc) # group2 = 1
    }
  }
  lmm <- lmer(logFC ~ 1 + group2 + (1 | subj.id), data = big.logFC.tab)
  lmm.summary <- summary(lmm)
  return(lmm.summary$coefficients["group2", ])
}

doManyQuantTest <- function(n.set, folder.path, n.subj, n.time, n.interval = n.time, sub.subj = n.subj / 2, rand.interval = FALSE) {
  registerDoParallel(cores = detectCores()-1)
  res.list <- foreach(i = 1:n.set) %dopar% quantitative.test(file.path = paste(folder.path, "/set", sprintf("%04d", i), ".txt", sep = ""), 
                                                             n.subj, n.time, n.interval, sub.subj, rand.interval)
  res <- matrix(unlist(res.list), nrow = length(res.list), byrow = TRUE)
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

mirkat.test <- function(file.path, n.subj, n.time, n.interval = n.time, sub.subj = n.subj / 2, rand.interval = FALSE, method = "bray") {
  sim.dataset <- subsample(file.path, n.subj, n.time, n.interval, sub.subj, rand.interval)
  
  otus <- t(sim.dataset)
  sampID <- rownames(otus)
  subjID <- sapply(strsplit(sampID, split = "[.]"), "[[", 1)
  time <- sapply(strsplit(sampID, split = "[.]"), "[[", 2)
  metadata <- data.frame(subjID, sampID, time)
  
  y <- c(rep(0, sub.subj), rep(1, sub.subj))
  Ds <- pldist_all(otus = otus, metadata = metadata, method = method)
  Ks <- lapply(Ds, FUN = function(d) D2K(d))
  res <- MiRKAT(y = y, Ks = Ks)
  return(unlist(res))
}

doManyMirkatTest <- function(n.set, folder.path, n.subj, n.time, n.interval = n.time, sub.subj = n.subj / 2, rand.interval = FALSE, method = "bray") {
  registerDoParallel(cores = detectCores()-1)
  res.list <- foreach(i = 1:n.set) %dopar% mirkat.test(file.path = paste(folder.path, "/set", sprintf("%04d", i), ".txt", sep = ""), 
                                                       n.subj, n.time, n.interval, sub.subj, rand.interval, method)
  res <- matrix(unlist(res.list), nrow = length(res.list), byrow = TRUE)
  return(res)
}

adhoc.test <- function(file.path, n.subj, n.time, n.interval = n.time, sub.subj = n.subj / 2, rand.interval = FALSE, method = "bray") {
  sim.dataset <- subsample(file.path, n.subj, n.time, n.interval, sub.subj, rand.interval)
  n.totalsubj <- sub.subj * 2 # total number of subsampled subjects
  subj.ind <- c(1:sub.subj, (n.subj/2 + 1):(n.subj/2 + sub.subj)) # subjIDs
  
  ## create data dataset with: distance, subj.id, group2
  ## binary (qualitative) and quantitative distances
  big.dist.tab.qual <- data.frame(distance = numeric(),
                                  group2 = numeric(),
                                  subj.id = character(),
                                  stringsAsFactors = FALSE)
  big.dist.tab.quant <- data.frame(distance = numeric(),
                                   group2 = numeric(),
                                   subj.id = character(),
                                   stringsAsFactors = FALSE)
  for (i in 1:n.totalsubj) {
    subj.id_i <- paste("SUBJ_", subj.ind[i], sep = "")
    sim.otutab_i <- as.data.frame(sim.dataset) %>%
      dplyr::select(contains(paste(subj.id_i, ".", sep = "")))
    
    n.pairs <- ncol(sim.otutab_i) - 1 # number of consecutive pairs of samples
    
    if (i <= sub.subj) {
      group2 <- rep(0, n.pairs) # not in group2 if subject in first half
    } else {
      group2 <- rep(1, n.pairs) # in group2 if subject in second half
    }
    
    dist.mat.qual <- as.matrix(vegdist(x = t(sim.otutab_i), method = method, binary = TRUE)) # sample-sample distance matrix
    dist.tab.qual_i <- data.frame(distance = diag(dist.mat.qual[2:nrow(dist.mat.qual), 1:(nrow(dist.mat.qual)-1)]),
                                  group2 = group2,
                                  subj.id = rep(subj.id_i, n.pairs))
    dist.mat.quant <- as.matrix(vegdist(x = t(sim.otutab_i), method = method, binary = FALSE)) # sample-sample distance matrix
    dist.tab.quant_i <- data.frame(distance = diag(dist.mat.quant[2:nrow(dist.mat.quant), 1:(nrow(dist.mat.quant)-1)]),
                                   group2 = group2,
                                   subj.id = rep(subj.id_i, n.pairs))
    ## rbind with data dataset
    big.dist.tab.qual <- rbind(big.dist.tab.qual, dist.tab.qual_i)
    big.dist.tab.quant <- rbind(big.dist.tab.quant, dist.tab.quant_i)
  }
  
  # LMM: distance ~ 1 + group2 + (1|subj.id)
  lmm.qual <- lmer(distance ~ 1 + group2 + (1 | subj.id), data = big.dist.tab.qual)
  lmm.qual.summary <- summary(lmm.qual)
  lmm.qual.coefs <- lmm.qual.summary$coefficients
  
  lmm.quant <- lmer(distance ~ 1 + group2 + (1 | subj.id), data = big.dist.tab.quant)
  lmm.quant.summary <- summary(lmm.quant)
  lmm.quant.coefs <- lmm.quant.summary$coefficients
  
  lmm.qual.pval <- lmm.qual.coefs["group2", ncol(lmm.qual.coefs)]
  lmm.quant.pval <- lmm.quant.coefs["group2", ncol(lmm.quant.coefs)]
  res <- matrix(c(lmm.qual.pval, lmm.quant.pval), nrow = 1, ncol = 2, dimnames = list(NULL, c("adhoc_qual_pval", "adhoc_quant_pval")))
  return(res)
}

doManyAdhocTest <- function(n.set, folder.path,	n.subj, n.time, n.interval = n.time, sub.subj = n.subj / 2, rand.interval = FALSE, method = "bray") {
  registerDoParallel(cores = detectCores())
  res <- foreach(i = 1:n.set) %dopar% adhoc.test(file.path = paste(folder.path, "/set", sprintf("%04d", i), ".txt", sep = ""),
                                                 n.subj, n.time, n.interval, sub.subj, rand.interval, method)
  res <- matrix(unlist(res), nrow = length(res), byrow = TRUE)
  return(res)
}






