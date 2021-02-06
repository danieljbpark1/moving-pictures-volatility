args = commandArgs(trailingOnly=TRUE)
# args[1] : number of subjects
# args[2] : number of simulated datsets
# args[3] : number of time points

library(microbiome)
library(dplyr)
library(tidyr)
library(dirmult)
library(doParallel)
source("./sim_functions.R") # for simulation functions
load("./dethlefsen_relman_for_sim.Rdata") # OTU tables and IDs from Dethlefsen-Relman study

## NULL MODEL 
## NO DIFFERENCE IN VOLATILITY BETWEEN GROUPS

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
load("./DirMultOutput_DR.Rda") # load Dirichlet-Multinomial object

set.seed(0) 
n.set = as.numeric(args[2]) # number of simulated datasets

foreach(set = 1:n.set) %dopar% simDM(set, n.subj, dd.DR, "./SimSets_DR_v01_n100")

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

## DISAPPEARANCE PROBABILITIES
disapp.logistic <- fit.disapp.model(input)
prob.disapp <- predict.logistic(disapp.logistic, otu.ids, n.subj, mean.group = 0, sd.group = 0.1) # prob. disapp. for each OTU
## REAPPEARANCE PROBABILITIES
reapp.logistic <- fit.reapp.model(input)
prob.reapp <- predict.logistic(reapp.logistic, otu.ids, n.subj, mean.group = 0, sd.group = 0.1) # prob. reapp. for each OTU
## SIMULATE REAPPEARANCES
otu.rarity.newdata <- otu.rarity.tab(otutab.D, otu.ids, subj.ids) # OTU average abundance with new subject ID levels
reapp.beta <- fit.reapp.beta(input) # model for reappeared relative abundance
beta.shapes.df <- get.beta.shapes(reapp.beta, otu.rarity.newdata) # predicted beta distributions
## CALCULATE WITHIN-SUBJECT QUINTILES OF AVERAGE RELATIVE ABUNDANCE
otu.quintiles <- otu.quintiles.tab(list(list(rel.otutab.D, "D"))) # subject D is baseline for OTU quintiles
scaled.sd <- get.scaled.sd(otu.quintiles, n.subj, base.sd = 0.75, scaling.factors = c(.5, .75, 1, 1.25, 1.5))

foreach(set = 1:n.set) %dopar% simLong(set, n.otu, n.subj, n.time, otu.ids, samp.ids, prob.disapp, prob.reapp, beta.shapes.df, scaled.sd, 
                                       read.path = "./SimSets_DR_v01_n100", write.path = "./SimSets_DR_v01_n100_t120")
