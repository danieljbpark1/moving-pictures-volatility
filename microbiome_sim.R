## SIMULATING LONGITUDINAL MICROBIOME COMPOSITION 
library(microbiome)
library(nlme)
library(lme4)
library(dplyr)
library(dirmult)
source("mp_analysis_functions.R")
load("mp_F4_data.Rdata")
load("mp_M3_data.Rdata")
load("dethlefsen_relman.Rdata")

otutab # OTU table of raw abundances, taxa are rows
n.subj <- 5 # number of subjects

# select only taxa that appear in at least 10% of samples
otutab <- otutab[apply(otutab > 0, 1, mean) > 0.1, ]
otutab.f4 <- otu.counts.f4[apply(otu.counts.f4 > 0, 1, mean) > 0.1, ]
otutab.m3 <- otu.counts.m3[apply(otu.counts.m3 > 0, 1, mean) > 0.1, ]
otutab.D <- dr.D.otutab[apply(dr.D.otutab > 0, 1, mean) > 0.1, ]
otutab.E <- dr.E.otutab[apply(dr.E.otutab > 0, 1, mean) > 0.1, ]
otutab.F <- dr.F.otutab[apply(dr.F.otutab > 0, 1, mean) > 0.1, ]

##
## SIMULATE OTU TABLES AT T=1
##

# estimate Dirichlet-multinomial parameters 
dd = dirmult(otutab) 
save(dd, file = "DirMultOutput.Rda")  
# load("DirMultOutput.Rda") 

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
for (set in 1:n.set) {
  # track progress 
  if (set %% 100 == 0) print(set)
  
  # simulate from Dirichlet-multinomial 
  S = simPop(J = n.subj, n = 1000, pi = dd$pi, theta = dd$theta)
  
  # save in pre-existing folder, labeled as set0001, set0002, etc. 
  write.table(S$data, file = paste("./SimSets", n.subj, "/set", sprintf("%04d", set), ".txt", sep = ""), 
              sep = "\t", col.names = F, row.names = F)
}

## 
## SIMULATE LONGITUDINAL OTUs (multiple time points -- balanced) 
## 

n.time = 4 # number of time points
n.otu # number of OTUs
prob.disapp <- predict.prob.disapp(input = input) # predict with SSasymp
prob.reapp <- predict.prob.reapp(input = input)
for (set in 1:n.set) {
  # track progress 
  if (set %% 100 == 0) print(set) 
  
  # set up matrix and fill in time 1 
  this.otus <- matrix(nrow = n.subj*n.time, ncol = n.otu)
  this.otus[c(1:n.subj), ] <- as.matrix(read.table(paste("./SimSets", n.subj, "/set", sprintf("%04d", set), ".txt", sep = "")))
  
  for (tt in 2:n.time) {
    prev.start <- (tt - 2)*n.subj + 1 
    prev.end <- (tt - 1)*n.subj
    prev.otus <- this.otus[c(prev.start:prev.end), ]
    
    # does taxon disappear? 
    indic.disapp <- sapply(prob.disapp, FUN = function(p) rbinom(1, size = 1, prob = (1-p)))  # indicator that cell does *not* disappear
    mat.disapp <- matrix(indic.disapp, nrow = n.subj) 
    
    # if not, how much change from t1 to t2? 
    perturb.lval <- sim.log.foldchange(n.subj*n.otu)
    perturb.val <- exp(perturb.lval) 
    mat.perturb <- matrix(perturb.val, nrow = n.subj)
    
    # final perturbation 
    this.start <- (tt - 1)*n.subj + 1 
    this.end <- tt*n.subj 
    this.otus[this.start:this.end, ] <- round(prev.otus * mat.disapp * mat.perturb) 
    
    # does taxon reappear?
    indic.reapp <- sapply(prob.reapp, FUN = function(p) rbinom(1, size = 1, prob = p)) # indicator that cell does reappear
    mat.reapp <- matrix(indic.reapp, nrow = n.subj)
  }
  
  # again write to pre-created folder 
  write.table(this.otus, 
              file = paste("./SimSets_n", n.subj, "_t", n.time, "/set", sprintf("%04d", set), ".txt", sep = ""), 
              sep = "\t", col.names = F, row.names = F)
}

input <- list(list(otutab.D, "D"),
              list(otutab.E, "E"),
              list(otutab.F, "F"))
# input: list of otutab, subj.id pairs
predict.prob.disapp <- function(input) {
  big.disapp.tab <- data.frame(prob.disappear = numeric(),
                               avg.asin = numeric(),
                               subj.id = character())
  for (subj in input) {
    subj.otutab <- subj[[1]]
    subj.id <- subj[[2]]
    big.disapp.tab <- rbind(big.disapp.tab, 
                            asin.disapp.tab(subj.otutab, subj.id))
  }
  
  big.disapp.tab <- big.disapp.tab %>%
    mutate(subj.id = as.factor(subj.id))
  
  grouped.disapp.tab <- groupedData(prob.disappear ~ avg.asin | subj.id, 
                                    data = big.disapp.tab)
  # grouped.disapp.tab <- grouped.disapp.tab %>% 
  #   drop_na()
  
  asymp.fixed <- nlsList(SSasymp, data = grouped.disapp.tab)
  asymp.mixed <- nlme(asymp.fixed, random = lrc ~ 1)
  
  return(predict(asymp.mixed))
}

asin.disapp.tab <- function(otutab, subj.id) {
  disapp.prob <- disappearance_probabilities(otutab)
  avg.asin <- avg.asin.abnd(otutab)
  disapp.tab <- data.frame(prob.disappear = disapp.prob,
                           avg.asin = avg.asin) %>%
    mutate(subj.id = subj.id)
  return(disapp.tab)
}

avg.asin.abnd <- function(otutab) {
  asin.otutab <- asin(transform(otutab, 
                                transform = 'hellinger', 
                                target = "OTU"))
  avg.asin <- apply(asin.otutab, 1, mean)
  return(avg.asin)
}

sim.log.foldchange <- function(n) {
  return(rnorm(n, 0, sd = 1))
}

predict.prob.reapp <- function(input) {
  big.reapp.tab <- data.frame(y = numeric(),
                              otu.rarity = numeric(),
                              otu.id = character(),
                              subj.id = character())
  for (subj in input) {
    subj.otutab <- subj[[1]]
    subj.id <- subj[[2]]
    big.reapp.tab <- rbind(big.reapp.tab, 
                           asin.reapp.tab(subj.otutab, subj.id))
  }
  big.reapp.tab <- big.reapp.tab %>%
    mutate(reappeared = if_else(y > 0, 1, 0),
           subj.id = as.factor(subj.id))
  
  reapp.glmm <- glmer(reappeared ~ 1 + (1 | subj.id / otu.id), 
                      data = big.reapp.tab,
                      family = binomial)
  
  pred.prob.reapp <- predict(reapp.glmm, 
                             re.form = ~ (1 | subj.id / otu.id), 
                             type = "response")
  return(pred.prob.reapp)
  # pred.reapp.df <- data.frame(pred.prob.reapp = pred.prob.reapp) %>%
  #   mutate(otu.id = big.reapp.tab$otu.id,
  #          subj.id = subj.id)
}

asin.reapp.tab <- function(otutab, subj.id) {
  avg.asin <- avg.asin.abnd(otutab)
  reapp.tab <- format_reappearance_data(otutab, 
                                        avg.asin, 
                                        subj.id)
  return(reapp.tab)
}

sim.reapp <- function() {
  
}






