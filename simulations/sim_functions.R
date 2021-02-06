library(nlme)
library(lme4)
library(glmmTMB)
library(microbiome)
library(phyloseq)
library(dplyr)
library(tidyr)
library(boot)
library(simstudy)
source("./mp_analysis_functions.R")

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

# fits a mixed effects logistic regression to observed data in input
fit.disapp.model <- function(input) {
  big.disapp.tab <- asin.disapp.tab(input) # table of possible/actual disappearances
  disapp.glmm <- glmer(disappeared ~ 1 + subj.id + (1 | otu.id),   # mixed effects logistic regression
                       data = big.disapp.tab,
                       family = binomial)
  return(disapp.glmm)
}

# returns matrix of simulated probabilities of dis/(re)appearance
# taxa are rows, subjects are columns
# uses fitted coefficients to simulate new probabilities where simulated subjects intercepts ~ N(mean.group, sd.group)
predict.logistic <- function(logistic.glmm, otu.ids, n.subj, mean.group, sd.group) {
  baseline.intercepts <- coef(logistic.glmm)[[1]]
  baseline.intercepts <- baseline.intercepts[otu.ids, 1] # baseline subject intercept + OTU random intercepts in order
  
  subj.intercepts <- rnorm(n = n.subj, mean = mean.group, sd = sd.group) # simulated subject fixed intercepts
  
  res <- matrix(nrow = length(otu.ids), ncol = n.subj) # return matrix with taxa as rows, subjects as columns
  for (i in 1:n.subj) {
    res[ ,i] <- inv.logit(baseline.intercepts + subj.intercepts[i]) # predict probabilities for all otu.id / subj.id
  }
  return(res)
}

# returns (KxN) x 4 dataframe of OTUs and their within-subject average relative abundance quintiles
otu.quintiles.tab <- function(input) {
  res <- data.frame(otu.id = character(),
                    subj.id = character(),
                    avg.rel.abnd = numeric(),
                    quintile = numeric())
  for (subj in input) {
    subj.otutab <- subj[[1]] # OTU relative abundance table
    subj.id <- subj[[2]] # current subject ID
    subj.res <- data.frame(otu.id = rownames(subj.otutab)) %>%
      mutate(subj.id = subj.id,
             avg.rel.abnd = apply(subj.otutab, 1, mean),
             quintile = ntile(avg.rel.abnd, n=5))
    res <- rbind(res, subj.res)
  }  
  return(res)
}

# return K x N matrix of scaled standard deviations
# base s.d. scaled by factor according to an OTU's quintiles
get.scaled.sd <- function(quintiles.tab, n.subj, base.sd = 1, scaling.factors = c(0.8, 0.9, 1, 1.1, 1.2)) {
  scaled.sd <- base.sd * scaling.factors[quintiles.tab$quintile] # multiply base standard deviation by scaling factor according to rarity quintile
  scaled.sd <- matrix(rep(scaled.sd, n.subj), ncol = n.subj)
  return(scaled.sd)
}

# returns K x N matrix of simulated log fold-changes
# intput is matrix of scaled standard deviations 
sim.log.foldchange <- function(scaled.sd) {
  res <- sapply(scaled.sd, FUN = function(s) rnorm(n = 1, mean = 0, sd = s))
  res <- matrix(res, ncol = ncol(scaled.sd)) # taxa are rows, subjects are columns
  return(res)
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

# fits a mixed effects logistic regression to observed data in input
fit.reapp.model <- function(input) {
  big.reapp.tab <- asin.reapp.tab(input) # table of possible/actual reappearances
  reapp.glmm <- glmer(reappeared ~ 1 + subj.id + (1 | otu.id),   # mixed effects logistic regression
                      data = big.reapp.tab,
                      family = binomial)
 return(reapp.glmm) 
}

# returns dataframe of otu.id, otu.rarity, subj.id
otu.rarity.tab <- function(otutab, otu.ids, subj.ids) {
  n.subj <- length(subj.ids)
  n.otu <- length(otu.ids)
  otu.rarity <- avg.asin.abnd(otutab)
  res <- data.frame(otu.id = rep(otu.ids, n.subj),
                    otu.rarity = rep(otu.rarity, n.subj),
                    subj.id = rep(subj.ids, each = n.otu))
  return(res)
}

# returns fitted glmmTMB mixed effects beta regression model
# reappearance relative abundance ~ otu.rarity + (otu.rarity|subj.id)
fit.reapp.beta <- function(input) { 
  reapp.abnd.df <- asin.reapp.tab(input) %>%    # table of all reappearances
    filter(reappeared == 1) %>% 
    mutate(subj.id = as.factor(subj.id))
  
  reapp.glmmTMB <- glmmTMB(reapp.abnd ~ otu.rarity + (otu.rarity|subj.id), 
                           data = reapp.abnd.df, 
                           family = beta_family,
                           control = glmmTMBControl(optimizer=optim, 
                                                    optArgs=list(method="BFGS")))
  return(reapp.glmmTMB)
}

# returns (KxN) x 2 dataframe of beta distribution shape1, shape2
# input is dataframe of beta distribution mean, precision
get.beta.shapes <- function(reapp.model, newdata) {
  beta.params.df <- data.frame(mu = predict(reapp.model, newdata = newdata, allow.new.levels = TRUE, type = "response"),
                               phi = predict(reapp.model, newdata = newdata, allow.new.levels = TRUE, type = "disp"))
  
  shapes.list <- apply(beta.params.df, 1, function(x) betaGetShapes(x[1], x[2]))
  shapes.df <- as.data.frame(do.call(rbind, shapes.list))
  return(shapes.df)
}

# returns K x N dataframe of simulated reappeared relative abundances
# input is dataframe of beta distribution shape1, shape2 parameters
sim.reapp <- function(shapes.df, n.subj) {
  # predict reappearances for all OTUs with new subject ID levels 
  res <- rbeta(n = nrow(shapes.df), shape1 = unlist(shapes.df[ ,1]), shape2 = unlist(shapes.df[ ,2]))
  res <- matrix(res, ncol = n.subj)
  return(res)
}

# required: a pre-existing folder
simDM <- function(set, n.subj, dd, write.path) {
  # simulate from Dirichlet-multinomial 
  S = simPop(J = n.subj, n = 10000, pi = dd$pi, theta = dd$theta) # matrix of OTU counts, taxa are columns
  S.otutab <- otu_table(t(S$data), taxa_are_rows = TRUE) # transpose matrix so that taxa are rows, subjects are columns
  S.rel.otutab <- transform(S.otutab, transform = "compositional") # transform into relative abundances
  
  # save in pre-existing folder, labeled as set0001, set0002, etc. 
  write.table(S.rel.otutab, file = paste(write.path, "/set", sprintf("%04d", set), ".txt", sep = ""), 
              sep = "\t", col.names = T, row.names = T)
}

# requires pre-existing folder
simLong <- function(set, n.otu, n.subj, n.time, otu.ids, samp.ids, prob.disapp, prob.reapp, beta.shapes.df, scaled.sd, read.path, write.path) {
  # set up matrix and fill in time 1 
  this.otus <- matrix(nrow = n.otu, ncol = n.subj*n.time, dimnames = list(otu.ids, samp.ids))
  this.otus[ ,c(1:n.subj)] <- as.matrix(read.table(paste(read.path, "/set", sprintf("%04d", set), ".txt", sep = "")))
  
  for (tt in 2:n.time) {
    prev.start <- (tt - 2)*n.subj + 1 
    prev.end <- (tt - 1)*n.subj
    prev.otus <- this.otus[ ,c(prev.start:prev.end)]
    
    # does taxon disappear? 
    indic.disapp <- sapply(prob.disapp, FUN = function(p) rbinom(1, size = 1, prob = (1-p)))  # indicator that cell does *not* disappear
    mat.disapp <- matrix(indic.disapp, ncol = n.subj) 
    
    # if not, how much change from t1 to t2? 
    perturb.lval <- sim.log.foldchange(scaled.sd)
    mat.perturb <- exp(perturb.lval) # exponentiate the log fold-changes
    
    # final perturbation 
    current.otus <- prev.otus * mat.disapp * mat.perturb # disappeared OTUs go to 0
    
    # does taxon reappear?
    indic.reapp <- sapply(prob.reapp, FUN = function(p) rbinom(1, size = 1, prob = p)) # indicator that cell *does* reappear
    mat.reapp <- matrix(indic.reapp, ncol = n.subj)
    
    # at what relative abundance does it reappear?
    mat.sim.reapp <- sim.reapp(beta.shapes.df, n.subj)
    
    # if prev.otu == 0 and mat.reapp == 1 then set to new reappeared relative abundance
    prev.absent <- prev.otus == 0 # matrix indices of previous absences
    mat.sim.reapp <- mat.sim.reapp * mat.reapp # non-reappeared OTUs go to 0
    current.otus[prev.absent] <- mat.sim.reapp[prev.absent] # fill in previously absent, reappeared OTUs
    
    # re-normalize
    this.start <- (tt - 1)*n.subj + 1 
    this.end <- tt*n.subj 
    this.otus[ ,this.start:this.end] <- transform(current.otus, transform = 'compositional')
  }
  
  # again write to pre-created folder 
  write.table(this.otus, 
              file = paste(write.path, "/set", sprintf("%04d", set), ".txt", sep = ""), 
              sep = "\t", col.names = T, row.names = T)
}
