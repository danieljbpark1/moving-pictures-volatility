library(microbiome)
library(phyloseq)
library(dplyr)
library(tidyr)
library(dirmult)
library(doParallel)
source("./simulations/sim_functions.R")
load("./simulations/dethlefsen_relman_for_sim.Rdata") 

## DIFFERENCE IN QUALITATIVE VOLATILITY BETWEEN 2 GROUPS
## SUBJECTS 1-50: GROUP 1
## SUBJECTS 51-100: GROUP 2 (HIGHER VOLATILITY)

n.otu <- length(otu.ids) # number of OTUs
n.subj <- 100 # number of subjects
subj.ids <- character(n.subj) # vector of subject IDs: SUBJ_[i]
for (i in 1:n.subj) {
  subj.ids[i] <- paste("SUBJ_", i, sep="")
}

# input for later functions
input <- list(list(rel.otutab.D, "D"),
              list(rel.otutab.E, "E"),
              list(rel.otutab.F, "F"))
##
## SIMULATE OTU TABLES AT T=1
##
load("DirMultOutput_DR.Rda") 

set.seed(0) 
n.set = 500 # number of simulated datasets

# required: a pre-existing folder
simDM <- function(set, n.subj, dd, folder.path) {
  # simulate from Dirichlet-multinomial 
  S = simPop(J = n.subj, n = 10000, pi = dd$pi, theta = dd$theta) # matrix of OTU counts, taxa are columns
  S.otutab <- otu_table(t(S$data), taxa_are_rows = TRUE) # transpose matrix so that taxa are rows, subjects are columns
  S.rel.otutab <- transform(S.otutab, transform = "compositional") # transform into relative abundances
  
  # save in pre-existing folder, labeled as set0001, set0002, etc. 
  write.table(S.rel.otutab, file = paste(folder.path, "/set", sprintf("%04d", set), ".txt", sep = ""), 
              sep = "\t", col.names = T, row.names = T)
}

registerDoParallel(cores=2) # parallel processing with 2 clusters
foreach(set = 1:n.set) %dopar% simDM(set, n.subj, dd.DR, "./SimSets_DR_v02_n100")

## 
## SIMULATE LONGITUDINAL OTUs (multiple time points -- balanced) 
## 

n.time = 120 # number of time points
samp.ids <- character(n.subj*n.time) # sample IDs: SUBJ_[i].[t] 
for (t in 1:n.time) {
  for (i in 1:n.subj) {
    samp.ids[(t-1)*n.subj + i] <- paste(subj.ids[i], t, sep=".")
  }
}

## DISAPPEARANCE PROBABILITIES
disapp.logistic <- fit.disapp.model(input)
prob.disapp_1 <- predict.prob.disapp(disapp.logistic, otu.ids, n.subj = 50, mean.group = 0, sd.group = 0.1) # group 1
prob.disapp_2 <- predict.prob.disapp(disapp.logistic, otu.ids, n.subj = 50, mean.group = 2, sd.group = 0.1) # group 2
prob.disapp <- cbind(prob.disapp_1, prob.disapp_2) # prob. disapp. for each OTU
## REAPPEARANCE PROBABILITIES
reapp.logistic <- fit.reapp.model(input)
prob.reapp_1 <- predict.prob.reapp(reapp.logistic, otu.ids, n.subj = 50, mean.group = 0, sd.group = 0.1) # group 1
prob.reapp_2 <- predict.prob.reapp(reapp.logistic, otu.ids, n.subj = 50, mean.group = 2, sd.group = 0.1) # group 2
prob.reapp <- cbind(prob.reapp_1, prob.reapp_2) # prob. reapp. for each OTU
## SIMULATE REAPPEARANCES
otu.rarity.newdata <- otu.rarity.tab(otutab.D, otu.ids, subj.ids) # OTU average abundance with new subject ID levels
reapp.beta <- fit.reapp.beta(input) # model for reappeared relative abundance
beta.params.df <- data.frame(mu = predict(reapp.beta, newdata = otu.rarity.newdata, allow.new.levels = TRUE, type = "response"),
                             phi = predict(reapp.beta, newdata = otu.rarity.newdata, allow.new.levels = TRUE, type = "disp"))
## CALCULATE WITHIN-SUBJECT QUINTILES OF AVERAGE RELATIVE ABUNDANCE
otu.quintiles <- otu.quintiles.tab(list(list(rel.otutab.D, "D")))

# requires pre-existing folder
simLong <- function(set, n.otu, n.subj, n.time, otu.ids, samp.ids,
                    prob.disapp, prob.reapp, reapp.data, reapp.glmm, 
                    avg.pos.rel.abnd, otu.quintiles, read.path, write.path) {
  # set up matrix and fill in time 1 
  this.otus <- matrix(nrow = n.otu, ncol = n.subj*n.time, dimnames = list(otu.ids, samp.ids))
  this.otus[, c(1:n.subj)] <- as.matrix(read.table(paste(read.path, "/set", sprintf("%04d", set), ".txt", sep = "")))
  
  for (tt in 2:n.time) {
    prev.start <- (tt - 2)*n.subj + 1 
    prev.end <- (tt - 1)*n.subj
    prev.otus <- this.otus[, c(prev.start:prev.end)]
    
    # does taxon disappear? 
    indic.disapp <- sapply(prob.disapp, FUN = function(p) rbinom(1, size = 1, prob = (1-p)))  # indicator that cell does *not* disappear
    mat.disapp <- matrix(indic.disapp, ncol = n.subj) 
    
    # if not, how much change from t1 to t2? 
    perturb.lval <- sim.log.foldchange(otu.quintiles, n.subj, base.sd = 0.75, scaling.factors = c(.5, .75, 1, 1.25, 1.5))
    mat.perturb <- exp(perturb.lval) # exponentiate the log fold-changes
    
    # final perturbation 
    current.otus <- prev.otus * mat.disapp * mat.perturb # disappeared OTUs go to 0
    
    # does taxon reappear?
    indic.reapp <- sapply(prob.reapp, FUN = function(p) rbinom(1, size = 1, prob = p)) # indicator that cell *does* reappear
    mat.reapp <- matrix(indic.reapp, ncol = n.subj)
    
    # at what relative abundance does it reappear?
    sim.reapp.abnd <- sim.reapp(reapp.means, reapp.disps, n.otu, n.subj)
    mat.sim.reapp <- matrix(sim.reapp.abnd, ncol = n.subj, byrow = TRUE)
    
    # if prev.otu == 0 and mat.reapp == 1 then set to new reappeared relative abundance
    prev.absent <- prev.otus == 0 # matrix indices of previous absences
    mat.sim.reapp <- mat.sim.reapp * mat.reapp # non-reappeared OTUs go to 0
    current.otus[prev.absent] <- mat.sim.reapp[prev.absent] # fill in previously absent, reappeared OTUs
    
    # re-normalize
    this.start <- (tt - 1)*n.subj + 1 
    this.end <- tt*n.subj 
    this.otus[, this.start:this.end] <- transform(current.otus, transform = 'compositional')
  }
  
  # again write to pre-created folder 
  write.table(this.otus, 
              file = paste("./SimSets_DR_n", n.subj, "_t", n.time, "/set", sprintf("%04d", set), ".txt", sep = ""), 
              sep = "\t", col.names = T, row.names = T)
}

foreach(set = 1:15) %dopar% simLong(set, n.otu, n.subj, n.time, 
                                    prob.disapp, prob.reapp, reapp.data, reapp.glmm, 
                                    avg.pos.rel.abnd, otu.quintiles, read.path = "./SimSets_DR_v02_n100", write.path = "./SimSets_DR_v02_n100_t120")


