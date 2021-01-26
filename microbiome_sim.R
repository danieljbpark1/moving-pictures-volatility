## SIMULATING LONGITUDINAL MICROBIOME COMPOSITION 
library(microbiome)
library(phyloseq)
library(dplyr)
library(tidyr)
library(dirmult)
library(doParallel)
source("./sim_functions.R")
# load("mp_F4_data.Rdata") # Moving Pictures subject F4
# load("mp_M3_data.Rdata") # Moving Pictures subject M3
load("dethlefsen_relman.Rdata") # Dethlefsen Relman subjects D, E, F

# start by filtering out singletons/doubletons
# join tables within studies but not across studies
otutab.D <- dr.D.otutab[apply(dr.D.otutab > 0, 1, sum) > 2, ]
otutab.E <- dr.E.otutab[apply(dr.E.otutab > 0, 1, sum) > 2, ]
otutab.F <- dr.F.otutab[apply(dr.F.otutab > 0, 1, sum) > 2, ]

dim(otutab.D); dim(otutab.E); dim(otutab.F)

otu.ids.D <- rownames(otutab.D)
otu.ids.E <- rownames(otutab.E)
otu.ids.F <- rownames(otutab.F)
# OTU IDs found in at least one subject
otu.ids <- unique(c(otu.ids.D,
                    otu.ids.E,
                    otu.ids.F))
n.otu <- length(otu.ids) # number of OTUs


otutab.D <- as.data.frame(otutab.D) %>%
  rownames_to_column(var = "otu.id")
otutab.D <- left_join(data.frame(otu.id = otu.ids),
                      otutab.D,
                      by = "otu.id") %>%
  mutate_all(list(~replace_na(.,0))) %>% 
  column_to_rownames("otu.id")

otutab.E <- as.data.frame(otutab.E) %>%
  rownames_to_column(var = "otu.id")
otutab.E <- left_join(data.frame(otu.id = otu.ids),
                      otutab.E,
                      by = "otu.id") %>%
  mutate_all(list(~replace_na(.,0))) %>% 
  column_to_rownames("otu.id")

otutab.F <- as.data.frame(otutab.F) %>%
  rownames_to_column(var = "otu.id")
otutab.F <- left_join(data.frame(otu.id = otu.ids),
                      otutab.F,
                      by = "otu.id") %>%
  mutate_all(list(~replace_na(.,0))) %>% 
  column_to_rownames("otu.id")

# combined otutab
otutab <- cbind(otutab.D,
                otutab.E,
                otutab.F)

# relative abundance OTU tables
rel.otutab.D <- transform(otutab.D, transform = 'compositional')
rel.otutab.E <- transform(otutab.E, transform = 'compositional')
rel.otutab.F <- transform(otutab.F, transform = 'compositional')

# input for later functions
input <- list(list(rel.otutab.D, "D"),
              list(rel.otutab.E, "E"),
              list(rel.otutab.F, "F"))

subj.ids <- c("D", "E", "F") 
n.subj <- length(subj.ids) # number of subjects

# dataframe of all otu.id / subj.id pairs
full.otu.subj.data <- data.frame(otu.id = rep(otu.ids, each = n.subj),
                                 subj.id = rep(subj.ids, n.otu))

# average within-subject present relative abundance for each OTU
avg.pos.rel.abnd <- matrix(nrow = n.otu, ncol = n.subj, dimnames = list(otu.ids, subj.ids)) 
for (i in 1:n.subj) {
  subj.otutab <- input[[i]][[1]]
  # avg. rel. abnd. presence for each OTU
  otu.avg.presence <- apply(subj.otutab, 1, function(a) mean(a[a > 0], na.rm = TRUE)) 
  otu.avg.presence[is.na(otu.avg.presence)] <- 0 # 0 if OTU never appears in subject
  avg.pos.rel.abnd[ ,i] <- otu.avg.presence
}

##
## SIMULATE OTU TABLES AT T=1
##

# estimate Dirichlet-multinomial parameters 
dd.DR = dirmult(t(otutab)) # taxa are columns

# save(dd.DR, file = "DirMultOutput_DR.Rda")
# load("DirMultOutput_DR.Rda") 

# Function parameters: 
#   J = number of subpopulations sampled = "number of samples"
#   K = number of different OTUs = length of pi
#   n = number of OTUs sampled in each subpopulation (i.e., read count)
#   pi = vector of OTU probabilities
#   theta = theta value used in dirmult() function = overdispersion parameter (due to intra-class correlation)
# Creates S datasets and saves them in (pre-existing) SimSets[n] folder
# Each simset is a JxK matrix with OTU counts. 

set.seed(0) 
n.set = 100 # number of simulated datasets

# required: a pre-existing folder
simDM <- function(set, n.subj, dd, folder.path) {
  # simulate from Dirichlet-multinomial 
  S = simPop(J = n.subj, n = 10000, pi = dd$pi, theta = dd$theta) # matrix of OTU counts, taxa are columns
  S.otutab <- otu_table(t(S$data), taxa_are_rows = TRUE) # transpose matrix so taxa are rows, subjects are columns
  S.rel.otutab <- transform(S.otutab, transform = "compositional") # transform into relative abundances
  
  # save in pre-existing folder, labeled as set0001, set0002, etc. 
  write.table(S.rel.otutab, file = paste(folder.path, "/set", sprintf("%04d", set), ".txt", sep = ""), 
              sep = "\t", col.names = T, row.names = T)
}

registerDoParallel(cores=2) # parallel processing with 2 clusters
foreach(set = 1:n.set) %dopar% simDM(set, n.subj, dd.DR, "./SimSets_DR3")

## 
## SIMULATE LONGITUDINAL OTUs (multiple time points -- balanced) 
## 

n.time = 120 # number of time points
samp.ids <- character() # sample identifiers
for (t in 1:n.time) {
  for (s in subj.ids) {
    samp.ids <- c(samp.ids, paste(s,t,sep = "."))
  }
}

## DISAPPEARANCE PROBABILITIES
prob.disapp <- predict.prob.disapp(input, full.otu.subj.data) # prob. disapp. for each OTU
## REAPPEARANCE PROBABILITIES
prob.reapp <- predict.prob.reapp(input, full.otu.subj.data) # prob. reapp. for each OTU
## SIMULATE REAPPEARANCES
reapp.data <- asin.reapp.tab(input) # training data for reapp. models
reapp.glmm <- fit.reapp.glmm(reapp.data) # model for reapp. rel. abnd.
## CALCULATE WITHIN-SUBJECT QUINTILES OF AVERAGE RELATIVE ABUNDANCE
otu.rarity.data <- rarity.quintiles.tab(input)

# requires pre-existing folder
simLong <- function(set, n.otu, n.subj, n.time, 
                    prob.disapp, prob.reapp, reapp.data, reapp.glmm, 
                    avg.pos.rel.abnd, otu.rarity.data, full.otu.subj.data) {
  # set up matrix and fill in time 1 
  this.otus <- matrix(nrow = n.otu, ncol = n.subj*n.time, dimnames = list(otu.ids, samp.ids))
  this.otus[, c(1:n.subj)] <- as.matrix(read.table(paste("./SimSets_DR", n.subj, "/set", sprintf("%04d", set), ".txt", sep = "")))
  
  for (tt in 2:n.time) {
    prev.start <- (tt - 2)*n.subj + 1 
    prev.end <- (tt - 1)*n.subj
    prev.otus <- this.otus[, c(prev.start:prev.end)]
    
    # does taxon disappear? 
    indic.disapp <- sapply(prob.disapp, FUN = function(p) rbinom(1, size = 1, prob = (1-p)))  # indicator that cell does *not* disappear
    mat.disapp <- matrix(indic.disapp, ncol = n.subj, byrow = TRUE) 
    
    # if not, how much change from t1 to t2? 
    perturb.lval <- sim.log.foldchange(otu.rarity.data, full.otu.subj.data, 
                                       base.sd = 0.5, scaling.factors = c(.4, .7, 1, 1.4, 2))
    perturb.val <- exp(perturb.lval) 
    mat.perturb <- matrix(perturb.val, ncol = n.subj, byrow = TRUE)
    
    # final perturbation 
    current.otus <- prev.otus * mat.disapp * mat.perturb # disappeared OTUs go to 0

    # does taxon reappear?
    indic.reapp <- sapply(prob.reapp, FUN = function(p) rbinom(1, size = 1, prob = p)) # indicator that cell does reappear
    mat.reapp <- matrix(indic.reapp, ncol = n.subj, byrow = TRUE)
    
    # at what relative abundance does it reappear?
    sim.reapp.abnd <- sim.reapp(reapp.data = reapp.data, 
                                reapp.glmmTMB = reapp.glmm, 
                                avg.pos.rel.abnd = avg.pos.rel.abnd, 
                                otu.subj.data = full.otu.subj.data)
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
                                       avg.pos.rel.abnd, otu.rarity.data, full.otu.subj.data)
