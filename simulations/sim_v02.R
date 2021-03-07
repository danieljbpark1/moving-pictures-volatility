args = commandArgs(trailingOnly=TRUE)

library(microbiome)
library(dplyr)
library(tidyr)
library(dirmult)
library(doParallel)
source("./sim_functions.R")
load("./simulations/dethlefsen_relman_for_sim.Rdata") 

## DIFFERENCE IN QUALITATIVE VOLATILITY BETWEEN 2 GROUPS
## SUBJECTS 1-50: GROUP 1
## SUBJECTS 51-100: GROUP 2 (HIGHER VOLATILITY)

n.otu <- length(otu.ids) # number of OTUs
n.subj <- as.numeric(args[1]) # number of subjects
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
load("./DirMultOutput_DR.Rda") # load pre-trained Dirichlet-Multinomial object named dd.DR

set.seed(0) 
n.set = as.numeric(args[2]) # number of simulated datasets

# only need to run once
foreach(set = 1:n.set) %dopar% simDM(set, n.subj, dd.DR, "./SimSets_DR_v02_n100") 

## 
## SIMULATE LONGITUDINAL OTUs (multiple time points -- balanced) 
## 

n.time = as.numeric(args[3]) # number of time points
samp.ids <- character(n.subj*n.time) # sample IDs: SUBJ_[i].[t] 
for (t in 1:n.time) {
  for (i in 1:n.subj) {
    samp.ids[(t-1)*n.subj + i] <- paste(subj.ids[i], t, sep=".")
  }
}

qual.delta <- as.numeric(args[4])
## DISAPPEARANCE PROBABILITIES
disapp.logistic <- fit.disapp.model(input)
## REAPPEARANCE PROBABILITIES
reapp.logistic <- fit.reapp.model(input)
## SIMULATE REAPPEARANCES
otu.rarity.newdata <- otu.rarity.tab(otutab.D, otu.ids, subj.ids) # OTU average abundance with new subject ID levels
reapp.beta <- fit.reapp.beta(input) # model for reappeared relative abundance
beta.shapes.df <- get.beta.shapes(reapp.beta, otu.rarity.newdata) # predicted beta distributions
## CALCULATE WITHIN-SUBJECT QUINTILES OF AVERAGE RELATIVE ABUNDANCE
otu.quintiles <- otu.quintiles.tab(list(list(rel.otutab.D, "D")))
scaled.sd <- get.scaled.sd(otu.quintiles, n.subj, base.sd = 0.75, scaling.factors = c(.5, .75, 1, 1.25, 1.5))

# requires existing directory at write.path argument
foreach(set = 1:n.set) %dopar% {
  prob.disapp_1 <- predict.logistic(disapp.logistic, otu.ids, n.subj = (n.subj/2), mean.group = 0, sd.group = 0.1) # group 1
  prob.disapp_2 <- predict.logistic(disapp.logistic, otu.ids, n.subj = (n.subj/2), mean.group = qual.delta, sd.group = 0.1) # group 2
  prob.disapp <- cbind(prob.disapp_1, prob.disapp_2) # prob. disapp. for each OTU
  
  prob.reapp_1 <- predict.logistic(reapp.logistic, otu.ids, n.subj = (n.subj/2), mean.group = 0, sd.group = 0.1) # group 1
  prob.reapp_2 <- predict.logistic(reapp.logistic, otu.ids, n.subj = (n.subj/2), mean.group = qual.delta, sd.group = 0.1) # group 2
  prob.reapp <- cbind(prob.reapp_1, prob.reapp_2) # prob. reapp. for each OTU
  
  simLong(set, n.otu, n.subj, n.time, otu.ids, samp.ids, prob.disapp, prob.reapp, beta.shapes.df, scaled.sd, 
          read.path = "./SimSets_DR_v02_n100", write.path = as.character(args[5]))
}

n <- 2
foreach(i = 1:10) %dopar% {
  print(n)
}

