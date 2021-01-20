## SIMULATING LONGITUDINAL MICROBIOME COMPOSITION 
library(microbiome)
library(nlme)
library(lme4)
library(glmmTMB)
library(dplyr)
library(tidyr)
library(dirmult)
library(doParallel)
source("mp_analysis_functions.R")
load("mp_F4_data.Rdata") # Moving Pictures subject F4
load("mp_M3_data.Rdata") # Moving Pictures subject M3
load("dethlefsen_relman.Rdata") # Dethlefsen Relman subjects D, E, F

# start by filtering out singletons/doubletons
# join tables within studies but not across studies
otutab.f4 <- otu.counts.f4[apply(otu.counts.f4 > 0, 1, sum) > 2, ]
otutab.m3 <- otu.counts.m3[apply(otu.counts.m3 > 0, 1, sum) > 2, ]

otu.ids.f4 <- rownames(otutab.f4)
otu.ids.m3 <- rownames(otutab.m3)
# OTU IDs found in at least one subject
otu.ids <- unique(c(otu.ids.f4,
                    otu.ids.m3))
n.otu <- length(otu.ids) # number of OTUs

otutab.f4 <- as.data.frame(otutab.f4) %>%
  rownames_to_column(var = "otu.id")
otutab.f4 <- left_join(data.frame(otu.id = otu.ids),
                       otutab.f4,
                       by = "otu.id") %>%
  mutate_all(list(~replace_na(.,0)))
otutab.f4 <- otutab.f4 %>% column_to_rownames("otu.id")

otutab.m3 <- as.data.frame(otutab.m3) %>%
  rownames_to_column(var = "otu.id")
otutab.m3 <- left_join(data.frame(otu.id = otu.ids),
                       otutab.m3,
                       by = "otu.id") %>%
  mutate_all(list(~replace_na(.,0)))
otutab.m3 <- otutab.m3 %>% column_to_rownames("otu.id")
# combined otutab
otutab <- cbind(otutab.f4,
                otutab.m3)

subj.ids <- c("F4", "M3") 
n.subj <- length(subj.ids) # number of subjects
# dataframe of all otu.id / subj.id pairs
full.otu.subj.data <- data.frame(otu.id = rep(otu.ids, each = n.subj),
                                 subj.id = rep(subj.ids, n.otu))

# relative abundance OTU tables
rel.otutab.f4 <- transform(otutab.f4, transform = 'compositional')
rel.otutab.m3 <- transform(otutab.m3, transform = 'compositional')

# input for later functions
input <- list(list(rel.otutab.f4, "F4"),
              list(rel.otutab.m3, "M3"))

# average within-subject present relative abundance for each OTU
avg.pos.rel.abnd <- matrix(nrow = n.otu, ncol = n.subj, dimnames = list(otu.ids, subj.ids)) 
for (i in 1:n.subj) {
  subj.otutab <- input[[i]][[1]]
  # avg. rel. abnd. presence for each OTU
  avg.pos.rel.abnd[ ,i] <- apply(subj.otutab, 1, function(a) mean(a[a > 0], na.rm = TRUE)) 
}
avg.pos.rel.abnd[is.na(avg.pos.rel.abnd)] <- 0 # 0 if OTU never appears in subject

##
## SIMULATE OTU TABLES AT T=1
##

# estimate Dirichlet-multinomial parameters 
# dd.mp = dirmult(t(otutab)) # taxa are columns

# save(dd.mp, file = "DirMultOutput_MP.Rda")  
# load("DirMultOutput_MP.Rda") 

# Function parameters: 
#   J = number of subpopulations sampled = "number of samples"
#   K = number of different OTUs = length of pi
#   n = number of OTUs sampled in each subpopulation (i.e., read count)
#   pi = vector of OTU probabilities
#   theta = theta value used in dirmult() function = overdispersion parameter (due to intra-class correlation)
# Creates S datasets and saves them in (pre-existing) SimSets[n] folder
# Each simset is a JxK matrix with OTU counts. 

set.seed(0) 
n.set = 500 # number of simulated datasets

simDM <- function(set, n.subj, dd.mp) {
  # track progress 
  if (set %% 100 == 0) print(set)
  
  # simulate from Dirichlet-multinomial 
  S = simPop(J = n.subj, n = 1000, pi = dd.mp$pi, theta = dd.mp$theta) # matrix of counts, taxa are columns
  S.otutab <- otu_table(t(S$data), taxa_are_rows = TRUE) # transpose matrix so taxa are rows, subjects are columns
  S.rel.otutab <- transform(S.otutab, transform = "compositional") # transform into relative abnd. table
  
  # save in pre-existing folder, labeled as set0001, set0002, etc. 
  write.table(S.rel.otutab, file = paste("./SimSets_MP", n.subj, "/set", sprintf("%04d", set), ".txt", sep = ""), 
              sep = "\t", col.names = T, row.names = T)
}

cl <- makeCluster(2, type="FORK") # parallel processing with 2 clusters
registerDoParallel(cl)
foreach(set = 1:n.set) %dopar% simDM(set, n.subj, dd.mp)

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

# returns dataframe of all post-presence instances 
asin.disapp.tab <- function(input) {
  big.disapp.tab <- data.frame(otu.id = character(),
                               disappeared = numeric(),
                               otu.rarity = numeric(),
                               subj.id = character())
  for (subj in input) {
    subj.otutab <- subj[[1]]
    subj.id <- subj[[2]]
    
    subj.otu.rarity <- avg.asin.abnd(subj.otutab)
    disapp.tab <- format_disappearance_data(otutab = subj.otutab,
                                            otu.rarity = subj.otu.rarity, 
                                            subj.id = subj.id)
    big.disapp.tab <- rbind(big.disapp.tab, disapp.tab)
  }
  return(big.disapp.tab)
}

# returns vector of pred. prob. of disapp. sorted by otu.ids / subj.id
predict.prob.disapp <- function(input, otu.subj.data) {
  big.disapp.tab <- asin.disapp.tab(input) 
  # mixed effects logistic regression
  disapp.glmm <- glmer(disappeared ~ 1 + subj.id + (1 | otu.id),
                        data = big.disapp.tab,
                        family = binomial)
  # predict for all subj.id / otu.id
  pred.prob.disapp <- predict(disapp.glmm, 
                              re.form = ~ (1 | otu.id), 
                              type = "response")
  # bind to training dataset
  big.disapp.tab$pred.prob.disapp <- pred.prob.disapp
  # choose with uniform probability one prediction for each subj.id / otu.id
  pred.disapp <- big.disapp.tab %>%
    group_by(subj.id, otu.id) %>%
    sample_n(size = 1) %>%
    dplyr::select(-disappeared)
  # join to otu.ids so that predictions are ordered by otu.id / subj.id
  res <- left_join(otu.subj.data,
                   pred.disapp,
                   by = c("otu.id", "subj.id"))
  # if there's no prediction for an OTU, it should always disappear
  # because it shouldn't be present in that subject
  res$pred.prob.disapp[is.na(res$pred.prob.disapp)] <- 1
  return(res$pred.prob.disapp)
}

# returns vector of simulated log fold-changes
sim.log.foldchange <- function(n) {
  return(rnorm(n = n, mean = 0, sd = 1))
}

# returns dataframe of all post-absence instances
asin.reapp.tab <- function(input) {
  big.reapp.tab <- data.frame(otu.id = character(),
                              reappeared = numeric(),
                              reapp.abnd = numeric(),
                              otu.rarity = numeric(),
                              subj.id = character())
  for (subj in input) {
    subj.otutab <- subj[[1]]
    subj.id <- subj[[2]]
    
    subj.otu.rarity <- avg.asin.abnd(subj.otutab)
    reapp.tab <- format_reappearance_data(subj.otutab,
                                          subj.otu.rarity,
                                          subj.id)
    big.reapp.tab <- rbind(big.reapp.tab, reapp.tab)
  }
  return(big.reapp.tab)
}

# returns vector of prob. of reapp. forsorted by otu.ids / subj.id
predict.prob.reapp <- function(reapp.data, otu.subj.data) {
  # mixed effects logistic regression
  reapp.glmm <- glmer(reappeared ~ 1 + subj.id + (1 | otu.id), 
                      data = reapp.data,
                      family = binomial)
  # predict for all rows in training data
  pred.prob.reapp <- predict(reapp.glmm, 
                             re.form = ~ (1 | otu.id), 
                             type = "response")
  
  big.reapp.tab <- reapp.data
  big.reapp.tab$pred.prob.reapp <- pred.prob.reapp
  # choose with uniform probability one prediction for each subj.id / otu.id in training data
  pred.reapp <- big.reapp.tab %>%
    group_by(subj.id, otu.id) %>%
    sample_n(size = 1) %>%
    dplyr::select(otu.id, subj.id, pred.prob.reapp)
  
  # NA predictions for otu.id / subj.id that were never absent (always present)
  # so fill in NA with 1
  res <- left_join(otu.subj.data,
                   pred.reapp,
                   by = c("otu.id", "subj.id")) %>%
    replace_na(replace = list(pred.prob.reapp = 1))
  
  return(res$pred.prob.reapp)
}

# returns fitted glmmTMB mixed effects beta regression model
# reappearance relative abundance ~ otu.rarity + (otu.rarity|subj.id)
fit.reapp.glmm <- function(reapp.data) {
  reapp.abnd.df <- reapp.data %>%
    filter(reappeared == 1) %>%
    mutate(subj.id = as.factor(subj.id))
  
  reapp.glmmTMB <- glmmTMB(reapp.abnd ~ otu.rarity + (otu.rarity|subj.id), 
                           data = reapp.abnd.df, 
                           family = beta_family)
  return(reapp.glmmTMB)
}

# returns dataframe of simulated reappeared relative abundances
# if otu.id / subj.id is not simulated by reapp.glmmTMB
# fill in with its average within-subject presence
sim.reapp <- function(reapp.data, reapp.glmmTMB, avg.pos.rel.abnd, otu.subj.data) {
  # all reappearance relative abundance instances
  reapp.abnd.df <- reapp.data %>%
    filter(reappeared == 1)
  # simulated relative abundances 
  sim.reapp <- simulate(reapp.glmmTMB)
  reapp.abnd.df$sim <- sim.reapp[ ,1]
  
  # choose with uniform probability one prediction for each subj.id / otu.id in training data
  sim.reapp.data <- reapp.abnd.df %>%
    group_by(subj.id, otu.id) %>%
    sample_n(size = 1)
  
  # NA simulations for otu.id / subj.id that were never absent (always present)
  # or are never present in a subject
  res <- left_join(otu.subj.data,
                   sim.reapp.data,
                   by = c("otu.id", "subj.id"))
  # average within-subject present relative abundance for each OTU
  avg.presence <- melt(data.table(t(avg.pos.rel.abnd), keep.rownames = TRUE), id.vars = 'rn')
  # fill in NA simulations with the OTU's average within-subject presence
  res$sim[is.na(res$sim)] <- avg.presence$value[is.na(res$sim)]
  
  return(res$sim)
}

## DISAPPEARANCE PROBABILITIES
prob.disapp <- predict.prob.disapp(input = input, otu.subj.data = full.otu.subj.data) # prob. disapp. for each OTU
## REAPPEARANCE PROBABILITIES
reapp.data <- asin.reapp.tab(input = input) # training data for reapp. models
prob.reapp <- predict.prob.reapp(reapp.data, full.otu.subj.data) # prob. reapp. for each OTU
## SIMULATE REAPPEARANCES
reapp.glmm <- fit.reapp.glmm(reapp.data = reapp.data) # model for reapp. rel. abnd.

simLong <- function(set, n.subj, n.otu, n.time, 
                    prob.disapp, prob.reapp, reapp.data, reapp.glmm, 
                    avg.pos.rel.abnd, full.otu.subj.data) {
  # track progress 
  if (set %% 100 == 0) print(set) 
  
  # set up matrix and fill in time 1 
  this.otus <- matrix(nrow = n.otu, ncol = n.subj*n.time, dimnames = list(otu.ids, samp.ids))
  this.otus[, c(1:n.subj)] <- as.matrix(read.table(paste("./SimSets_MP", n.subj, "/set", sprintf("%04d", set), ".txt", sep = "")))
  
  for (tt in 2:n.time) {
    prev.start <- (tt - 2)*n.subj + 1 
    prev.end <- (tt - 1)*n.subj
    prev.otus <- this.otus[, c(prev.start:prev.end)]
    
    # does taxon disappear? 
    indic.disapp <- sapply(prob.disapp, FUN = function(p) rbinom(1, size = 1, prob = (1-p)))  # indicator that cell does *not* disappear
    mat.disapp <- matrix(indic.disapp, ncol = n.subj, byrow = TRUE) 
    
    # if not, how much change from t1 to t2? 
    perturb.lval <- sim.log.foldchange(n = n.subj*n.otu)
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
    current.absent <- prev.otus == 0 # matrix indices of current absences
    mat.sim.reapp <- mat.sim.reapp * mat.reapp # non-reappeared OTUs go to 0
    current.otus[current.absent] <- mat.sim.reapp[current.absent] # fill in currently absent, reappeared OTUs
    
    # re-normalize
    this.start <- (tt - 1)*n.subj + 1 
    this.end <- tt*n.subj 
    this.otus[, this.start:this.end] <- transform(current.otus, transform = 'compositional')
  }
  
  # again write to pre-created folder 
  write.table(this.otus, 
              file = paste("./SimSets_MP_n", n.subj, "_t", n.time, "/set", sprintf("%04d", set), ".txt", sep = ""), 
              sep = "\t", col.names = T, row.names = T)
}

foreach(set = 1:n.set) %dopar% simLong(set, n.subj, n.otu, n.time, 
                                       prob.disapp, prob.reapp, reapp.data, reapp.glmm, 
                                       avg.pos.rel.abnd, full.otu.subj.data)
