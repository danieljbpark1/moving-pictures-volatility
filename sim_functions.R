library(nlme)
library(lme4)
library(glmmTMB)
library(microbiome)
library(dplyr)
library(tidyr)
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

# returns vector of pred. prob. of disapp. sorted by otu.ids / subj.id
predict.prob.disapp <- function(input, otu.subj.data) {
  big.disapp.tab <- asin.disapp.tab(input) 
  # mixed effects logistic regression
  disapp.glmm <- glmer(disappeared ~ 1 + subj.id + (1 | otu.id),
                       data = big.disapp.tab,
                       family = binomial)
  # predict for all otu.id / subj.id
  # bind to training dataset
  big.disapp.tab$pred.prob.disapp <- predict(disapp.glmm, 
                                             type = "response")
  # choose with uniform probability one prediction for each subj.id / otu.id
  pred.disapp <- big.disapp.tab %>%
    group_by(subj.id, otu.id) %>%
    sample_n(size = 1) %>%
    dplyr::select(-disappeared)
  # join to otu.ids so that predictions are ordered by otu.id / subj.id
  res <- left_join(otu.subj.data,
                   pred.disapp,
                   by = c("otu.id", "subj.id"))
  # if prediction is NA for an OTU, it was never present in that subject
  # fill in NAs with draws from Unif(0.5, 1)
  res.pred <- res$pred.prob.disapp
  res.pred[is.na(res.pred)] <- runif(n = sum(is.na(res.pred)), min = 0.5, max = 1)
  return(res.pred)
}

# returns dataframe of OTUs and their within-subject average relative abundance quintiles
rarity.quintiles.tab <- function(input) {
  res <- data.frame(otu.id = character(),
                    subj.id = character(),
                    avg.rel.abnd = numeric(),
                    quintile = numeric())
  for (subj in input) {
    subj.otutab <- subj[[1]] # OTU relative abundance table
    subj.id <- subj[[2]] 
    subj.res <- data.frame(otu.id = rownames(subj.otutab)) %>%
      mutate(subj.id = subj.id,
             avg.rel.abnd = apply(subj.otutab, 1, mean),
             quintile = ntile(avg.rel.abnd, n=5))
    res <- rbind(res, subj.res)
  }  
  return(res)
}

# returns vector of simulated log fold-changes
sim.log.foldchange <- function(quintiles.tab, otu.subj.data, base.sd = 1, scaling.factors = c(0.8, 0.9, 1, 1.1, 1.2)) {
  quintiles.tab$scaled.sd <- base.sd * scaling.factors[quintiles.tab$quintile]
  quintiles.tab$sim.log.foldchange <- sapply(quintiles.tab$scaled.sd, FUN = function(s) rnorm(n=1, mean = 0, sd = s))
  res <- left_join(x=otu.subj.data,
                   y=quintiles.tab,
                   by=c("otu.id","subj.id"))
  return(res$sim.log.foldchange)
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
predict.prob.reapp <- function(input, otu.subj.data) {
  big.reapp.tab <- asin.reapp.tab(input)
  # mixed effects logistic regression
  reapp.glmm <- glmer(reappeared ~ 1 + subj.id + (1 | otu.id), 
                      data = big.reapp.tab,
                      family = binomial)
  # predict for all otu.id / subj.id
  # bind to training dataset
  big.reapp.tab$pred.prob.reapp <- predict(reapp.glmm, 
                                           type = "response")
  # choose with uniform probability one prediction for each otu.id / subj.id in training data
  pred.reapp <- big.reapp.tab %>%
    group_by(subj.id, otu.id) %>%
    sample_n(size = 1) %>%
    dplyr::select(otu.id, subj.id, pred.prob.reapp)
  # join to otu.ids so that predictions are ordered by otu.id / subj.id
  res <- left_join(otu.subj.data,
                   pred.reapp,
                   by = c("otu.id", "subj.id"))
  # NA predictions for otu.id / subj.id that were always present
  # so fill in NA with 1
  res.pred <- res$pred.prob.reapp
  res.pred[is.na(res.pred)] <- 1
  return(res.pred)
}

# returns fitted glmmTMB mixed effects beta regression model
# reappearance relative abundance ~ otu.rarity + (otu.rarity|subj.id)
fit.reapp.glmm <- function(reapp.data) {
  reapp.abnd.df <- reapp.data %>%
    filter(reappeared == 1) %>%
    mutate(subj.id = as.factor(subj.id))
  
  reapp.glmmTMB <- glmmTMB(reapp.abnd ~ otu.rarity + (otu.rarity|subj.id), 
                           data = reapp.abnd.df, 
                           family = beta_family,
                           control = glmmTMBControl(optimizer=optim, 
                                                    optArgs=list(method="BFGS")))
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
  
  # NA simulations for otu.id / subj.id that were always present
  res <- left_join(otu.subj.data,
                   sim.reapp.data,
                   by = c("otu.id", "subj.id"))
  # average within-subject present relative abundance for each OTU
  avg.presence <- melt(data.table(t(avg.pos.rel.abnd), keep.rownames = TRUE), id.vars = 'rn')
  # fill in NA simulations with the OTU's average within-subject presence
  reapp.abnd <- res$sim
  reapp.abnd[is.na(reapp.abnd)] <- avg.presence$value[is.na(reapp.abnd)]
  return(reapp.abnd)
}
