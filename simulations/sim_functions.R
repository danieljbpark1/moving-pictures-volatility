library(nlme)
library(lme4)
library(glmmTMB)
library(microbiome)
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

# returns matrix of simulated probabilities of disappearance
# taxa are rows, subjects are columns
# uses fitted coefficients to simulate new probabilities where simulated subjects intercepts ~ N(mean.group, sd.group)
predict.prob.disapp <- function(disapp.glmm, otu.ids, n.subj, mean.group, sd.group) {
  baseline.intercepts <- coef(disapp.glmm)[[1]]
  baseline.intercepts <- baseline.intercepts[otu.ids, 1] # baseline subject intercept + OTU random intercepts in order
  
  subj.intercepts <- rnorm(n = n.subj, mean = mean.group, sd = sd.group) # simulated subject fixed intercepts
  
  res <- matrix(nrow = length(baseline.intercepts), ncol = n.subj) # return matrix with taxa as rows, subjects as columns
  for (i in 1:n.subj) {
    res[ ,i] <- inv.logit(baseline.intercepts + subj.intercepts[i]) # predict probabilities for all otu.id / subj.id
  }
  return(res)
}

# returns dataframe of OTUs and their within-subject average relative abundance quintiles
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

# returns matrix of simulated log fold-changes
sim.log.foldchange <- function(quintiles.tab, n.subj, base.sd = 1, scaling.factors = c(0.8, 0.9, 1, 1.1, 1.2)) {
  res <- matrix(nrow = nrow(quintiles.tab), ncol = n.subj) # taxa are rows, subjects are columns
  scaled.sd <- base.sd * scaling.factors[quintiles.tab$quintile] # multiply base standard deviation by scaling factor according to rarity quintile
  for (i in 1:n.subj) {
    res[ ,i] <- sapply(scaled.sd, FUN = function(s) rnorm(n = 1, mean = 0, sd = s)) 
  }
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

# returns matrix of simulated probabilities of reappearance
# taxa are rows, subjects are columns
# uses fitted coefficients to simulate new probabilities where simulated subjects intercepts ~ N(mean.group, sd.group)
predict.prob.reapp <- function(reapp.glmm, otu.ids, n.subj, mean.group, sd.group) {
  baseline.intercepts <- coef(reapp.glmm)[[1]]
  baseline.intercepts <- baseline.intercepts[otu.ids, 1] # baseline subject intercept + OTU random intercepts in order
  
  subj.intercepts <- rnorm(n = n.subj, mean = mean.group, sd = sd.group) # simulated subject fixed intercepts
  
  res <- matrix(nrow = length(baseline.intercepts), ncol = n.subj) # return matrix with taxa as rows, subjects as columns
  for (i in 1:n.subj) {
    res[ ,i] <- inv.logit(baseline.intercepts + subj.intercepts[i]) # predict probabilities for all otu.id / subj.id
  }
  return(res)
  # NA predictions for otu.id / subj.id that were always present across all observed subjects
  # so fill in NA with 1
}

# returns dataframe of otu.id, otu.rarity, subj.id
otu.rarity.tab <- function(otutab, otu.ids, subj.ids) {
  otu.rarity <- avg.asin.abnd(otutab)
  res <- data.frame(otu.id = rep(otu.ids, length(subj.ids)),
                    otu.rarity = rep(otu.rarity, length(subj.ids)),
                    subj.id = rep(subj.ids, each = length(otu.ids)))
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

# returns dataframe of beta distribution shape1, shape2
# input is dataframe of beta distribution mean, precision
get.beta.params <- function(df) {
  shapes.list <- apply(df, 1, function(x) betaGetShapes(x[1], x[2]))
}

# returns dataframe of simulated reappeared relative abundances
# if otu.id / subj.id is not simulated by reapp.glmmTMB
# fill in with its average within-subject presence
sim.reapp <- function(beta.means, beta.disps, n.otu, n.subj) {
  # predict reappearances for all OTUs with new subject ID levels 
  res <- numeric(n.otu*n.subj)
  for (i in 1:length(beta.means)) {
    shapes <- betaGetShapes(beta.means[i], beta.disps[i])
    res[i] <- rbeta(n = 1, shapes$shape1, shapes$shape2)
  }
  res <- matrix(res, nrow = n.otu, ncol = n.subj)
  return(res)
}
