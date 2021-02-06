library(boot)
library(data.table)
library(tibble)
library(microbiome)

# returns vector of avg. asin Hellinger-transformed abnd. for each OTU
avg.asin.abnd <- function(otutab) {
  asin.otutab <- asin(transform(otutab, 
                                transform = 'hellinger', 
                                target = "OTU"))
  avg.asin <- apply(asin.otutab, 1, mean)
  return(avg.asin)
}

# returns table of fold-changes between consecutive samples
# otutab : phyloseq OTU table
# taxa_are_rows : are taxa rows or columns?
# logRatio : option to return log fold-changes instead
fold_difference_table <- function(otutab, taxa_are_rows = TRUE, logRatio = FALSE) {
  if (taxa_are_rows) {
    num.samples <- ncol(otutab)
    fold.difference.table <- otutab[ , 2:num.samples] / otutab[ , 1:(num.samples-1)]  
  }
  else {
    num.samples <- nrow(otutab)
    fold.difference.table <- otutab[2:num.samples, ] / otutab[1:(num.samples-1), ]
  }
  
  if (logRatio) {
    fold.difference.table <- log(fold.difference.table)
  }
  
  return(fold.difference.table)
}

change.unweighted <- function(otutab, otu.rarity, subj.id, change.direction, taxa_are_rows = TRUE) {
  if (!taxa_are_rows) {
    otutab = t(otutab)
  }
  
  otu.rarity.df <- data.frame(otu.rarity = otu.rarity) %>%
    rownames_to_column(var = "otu.id")
  
  fold.diff.tab <- fold_difference_table(otutab) 
  fold.diff.dt <- setDT(data.frame(fold.diff.tab), 
                        keep.rownames = TRUE)
  fold.diff.melted <- melt(data = fold.diff.dt,
                           id.vars = "rn") %>%
    rename(otu.id = rn)
  
  if(change.direction == "disappearance") {
    disapp.df <- fold.diff.melted %>%
      filter(is.finite(value)) %>%
      mutate(disappeared = if_else(value == 0, 1, 0)) %>%
      dplyr::select(otu.id, disappeared)
    
    res.df <- left_join(x = disapp.df, 
                        y = otu.rarity.df,
                        by = "otu.id") %>%
      mutate(subj.id = subj.id)
  }
  else if(change.direction == "reappearance") {
    reapp.abnd.tab <- setDT(data.frame(otutab[ ,-1]),
                            keep.rownames = TRUE)
    reapp.abnd.melted <- melt(data = reapp.abnd.tab,
                              id.vars = "rn") %>%
      rename(otu.id = rn,
             reapp.abnd = value)
    
    fold.diff.melted$reapp.abnd <- reapp.abnd.melted$reapp.abnd
    
    reapp.df <- fold.diff.melted %>%
      filter(is.infinite(value) | is.nan(value)) %>%
      mutate(reappeared = if_else(is.infinite(value), 1, 0)) %>%
      dplyr::select(otu.id, reappeared, reapp.abnd)
    
    res.df <- left_join(x = reapp.df,
                        y = otu.rarity.df, 
                        by = "otu.id") %>%
      mutate(subj.id = subj.id)
  }
  else {return("change.direction argument must be either 'disappearance' or 'reappearance'.")}
  
  return(res.df)
}

# format 1 for disappeared, 0 for continued presence
format_disappearance_data <- function(otutab, otu.rarity, subj.id, taxa_are_rows = TRUE) {
  return(change.unweighted(otutab, otu.rarity, subj.id,
                           change.direction = "disappearance"))
  # if (!taxa_are_rows) {
  #   otutab = t(otutab)
  # }
  # 
  # num.samples <- ncol(otutab)
  # num.otus <- nrow(otutab)
  # 
  # fold.diff.tab <- fold_difference_table(otu.table = otutab) 
  # candidate.indices <- is.finite(fold.diff.tab) 
  # delta.time <- t.sample[2:num.samples] - t.sample[1:(num.samples-1)]
  # read.depths <- read.depths[2:num.samples]
  # 
  # disappeared <- numeric(); median.abundance <- numeric(); sample.read.depth <- numeric(); time.diff <- numeric(); otu_ids <- character()
  # for (i in 1:num.otus) {
  #   otu.candidate.indices <- candidate.indices[i, ]
  #   if (!any(otu.candidate.indices)) {
  #     next
  #   }
  #   
  #   otu.candidate.samples <- fold.diff.tab[i, otu.candidate.indices]
  #   
  #   disappeared <- append(disappeared, ifelse(otu.candidate.samples == 0, 1, 0))
  #   median.abundance <- append(median.abundance, rep(otu.median.abs[i], length(otu.candidate.samples)))
  #   sample.read.depth <- append(sample.read.depth, read.depths[otu.candidate.indices])
  #   time.diff <- append(time.diff, delta.time[otu.candidate.indices])
  #   otu_ids <- append(otu_ids, rep(otu_id[i], length(otu.candidate.samples)))
  # }
  # 
  # return(data.frame(disappeared, median.abundance, sample.read.depth, time.diff, otu_ids, row.names = NULL))
}

# function for formatting reappearance data for regression
# inputs: 
#   OTU table with samples ordered by time and OTUs in order
#   otu rarity metric (either raw, relative, or clr-transformed) of OTUs in order
# returns: data frame with response y being the magnitude of a reappeared count, including zeroes
format_reappearance_data <- function(otutab, otu.rarity, subj.id, taxa_are_rows = TRUE) {
  return(change.unweighted(otutab, otu.rarity, subj.id, 
                           change.direction = "reappearance"))
  # if (!taxaAreRows) {
  #   otutab = t(otutab)
  # }
  # 
  # num.samples <- ncol(otutab)
  # num.otus <- nrow(otutab)
  # otuID <- rownames(otutab)
  # 
  # # samples from 2 to N
  # reappearance.candidate.samples <- otutab[ , 2:num.samples]
  # fold.diff.tab <- fold_difference_table(otutab) 
  # # matrix indicating reappearance by OTU and sample
  # reappearance.indices <- is.infinite(fold.diff.tab) | is.nan(fold.diff.tab)
  # 
  # y <- numeric()
  # otu.rarity <- numeric()
  # otu.id <- character()
  # for (i in 1:num.otus) {
  #   otu.reappearance.samples <- reappearance.indices[i, ]
  #   # if OTU never reappeared, move on
  #   if (!any(otu.reappearance.samples)) {
  #     next
  #   }
  #   # extract the abundances after absences
  #   otu.reappearance.magnitudes <- reappearance.candidate.samples[i, which(otu.reappearance.samples==TRUE)]
  #   otu.reappearance.magnitudes <- as.numeric(otu.reappearance.magnitudes)
  #   
  #   y <- append(y, otu.reappearance.magnitudes)
  #   otu.rarity <- append(otu.rarity, rep(otuRarity[i], length(otu.reappearance.magnitudes)))
  #   otu.id <- append(otu.id, rep(otuID[i], length(otu.reappearance.magnitudes)))
  # }
  # 
  # res <- data.frame(y, otu.rarity, otu.id)
  # res$subj.id <- rep(subjID, nrow(res))
  # 
  # return(res)
}

prob.disapp.tab <- function(otutab, subj.id) {
  disapp.tab <- data.frame(otu.id = rownames(otutab),
                           prob.disappear = disappearance_probabilities(otutab),
                           avg.asin = avg.asin.abnd(otutab)) %>%
    mutate(subj.id = subj.id)
  return(disapp.tab)
}

prob.reapp.tab <- function(otutab, subj.id) {
  reapp.tab <- data.frame(otu.id = rownames(otutab),
                          prob.reappear = reappearance_probabilities(otutab),
                          avg.asin = avg.asin.abnd(otutab)) %>%
    mutate(subj.id = subj.id)
  return(reapp.tab)
}

# function for predicting reappearance counts from negative binomial model
predict_reappearance <- function(X_count, X_zero=NULL, negbin_model, zero_inflated=FALSE, seed = 0){
  set.seed(seed)
  n <- length(X_count[[1]])
  # calculate predicted probability of Y = 0
  if (zero_inflated) {
    negbin_coeffs <- negbin_model$coefficients$count
    if (!is.na(negbin_model$coefficients$zero[2])) {
      lambda <- negbin_model$coefficients$zero[1]
      for (i in 1:length(X_zero)) {
        lambda <- lambda + negbin_model$coefficients$zero[i+1]*X_zero[[i]]
      }
      pr.zero <- inv.logit(lambda)
    }
    else {
      pr.zero <- inv.logit(negbin_model$coefficients$zero[1])
    }
    # simulate n Bernoulli trials using p = pr.zero
    pred.zero <- rbinom(n = n, size = 1, prob = pr.zero)
  }
  else {
    negbin_coeffs <- negbin_model$coefficients
    pred.zero <- 0
  }
  
  # predicted NegBin distribution means are linked by log
  pred.nbinom.means <- exp(negbin_coeffs[1])
  for (i in 1:length(X_count)) {
    pred.nbinom.means <- pred.nbinom.means * exp(negbin_coeffs[i+1]*X_count[[i]])
  }

  dispersion <- negbin_model$theta
  pred.y <- numeric()
  for (i in 1:n) {
    mu <- pred.nbinom.means[i]
    p <- dispersion / (mu + dispersion)
    pred.y[i] <- rnbinom(n = 1, size = dispersion, prob = p)
  }
  
  pred.y[pred.zero == 1] <- 0
  return(pred.y)
}

# function for calculating probability of reappearance for each OTU
# inputs: OTU table with samples ordered by time 
# returns: vector of reappearance probabilities for each OTU
reappearance_probabilities <- function(otutab, taxa_are_rows = TRUE){
  fold.difference.table <- fold_difference_table(otutab, taxa_are_rows = taxa_are_rows)
  
  if (taxa_are_rows) {
    dimension <- 1
  } 
  else {
    dimension <- 2
  }
  
  # number of times where taxa reappeared
  reappearances <- apply(fold.difference.table, dimension, function(a) sum(is.infinite(a)))
  # number of times where taxa were absent at time t-1
  absences <- apply(fold.difference.table, dimension, function(a) sum(!is.finite(a)))
  return(as.numeric(reappearances / absences))
}

# function for simulating inflated beta distribution 
# inputs: coefficients for Beta mean and dispersion, coefficients for probability of one-inflation
# returns: draws from one-inflated Beta distribution with probability of one, Beta mean and dispersion regressed on predictors
sim.zoib <- function(zoib.model, beta.mean.coefs, beta.disp.coefs, infl.coefs, infl = 0, seed = 0){
  set.seed(seed)
  n <- nrow(zoib.model$Xb)
  
  # calculate predicted probability of Y = infl
  # probability of inflation uses logit link function
  if (infl == 0) {
    pr.infl <- inv.logit(as.matrix(zoib.model$Xb0) %*% infl.coefs)
  }
  else {
    pr.infl <- inv.logit(as.matrix(zoib.model$Xb1) %*% infl.coefs)
  }
  
  # simulate n Bernoulli trials using p = pr.infl
  pred.infl <- rbinom(n = n, size = 1, prob = pr.infl)
  
  # predicted Beta distribution means are linked by logit link function
  pred.beta.mean <- inv.logit(as.matrix(zoib.model$Xb) %*% beta.mean.coefs)
  
  # predicted Beta distribution dispersions are linked by log
  pred.beta.disp <- exp(as.matrix(zoib.model$Xd) %*% beta.disp.coefs)
  
  pred.alpha <- pred.beta.mean * pred.beta.disp
  pred.beta <- pred.beta.disp * (1 - pred.beta.mean)
  
  pred.y <- numeric()
  for (i in 1:n) {
    pred.y[i] <- rbeta(n = 1, shape1 = pred.alpha[i], shape2 = pred.beta[i])
  }
  # set inflated values
  pred.y[pred.infl == 1] <- infl
  return(pred.y)
}

# inputs: OTU table with samples ordered by time 
# returns: vector of disappearance probabilities named by OTU
disappearance_probabilities <- function(otutab, taxa_are_rows = TRUE) {
  fold.difference.table <- fold_difference_table(otutab, taxa_are_rows = taxa_are_rows)
  
  if (taxa_are_rows) {
    dimension <- 1
  } 
  else {
    dimension <- 2
  }
  
  # a count of times at which the taxa disappeared
  disappearances <- apply(fold.difference.table, dimension, function(a) sum(a == 0, na.rm = TRUE))
  # a count of times at which the taxa were present at time t-1
  presences <- apply(fold.difference.table, dimension, function(a) sum(is.finite(a), na.rm = TRUE))
  return(as.numeric(disappearances / presences))
}

# function for predicting an asymptotic regression
predict_AsympReg <- function(arr, NLSstAsymptotic_params) {
  b0 <- NLSstAsymptotic_params[1]
  b1 <- NLSstAsymptotic_params[2]
  lrc <- NLSstAsymptotic_params[3]
  return(b0 + b1*(1-exp(-exp(lrc)*arr)))
}

# function to format the time series dataframe
format_fold_change_data <- function(otu.tab) {
  # calculate consecutive fold changes
  fold.change.tab <- fold_difference_table(otu.tab)
  num.samples <- ncol(fold.change.tab)
  num.otus <- nrow(fold.change.tab)
  otu.names <- rownames(fold.change.tab)
  
  # appends rows (OTUs) one after another
  flat.fold.change <- as.vector(t(fold.change.tab))
  flat.otu.names <- rep(otu.names, each = num.samples)
  
  formatted.fold.change.df <- data.frame(otu.id = flat.otu.names, fold.change = flat.fold.change)
  
  # select only the rows containing finite and positive fold changes
  # no disappearances, no reappearances, no absence to absences
  formatted.fold.change.df <- subset(formatted.fold.change.df, is.finite(fold.change) & fold.change > 0)
  
  # create column for log fold change
  formatted.fold.change.df$log.fold.change <- log(formatted.fold.change.df$fold.change)
  return(formatted.fold.change.df)
}

# function to forecast an ARMA(1,1)
# older lags come first in the formula
predictArma <- function(arma.model, n.predict) {
  num.ar <- arma.model$arma[1]
  num.ma <- arma.model$arma[2]
  
  x.lags = as.vector(tail(arma.model$fitted, num.ar))
  e.lags = as.vector(tail(arma.model$residuals, num.ma))
  
  arma.coefs <- arma.model$coef
  
  intercept = arma.coefs["intercept"]
  ar.coefs = rev(arma.coefs[grep("ar", names(arma.coefs))])
  ma.coefs = rev(arma.coefs[grep("ma", names(arma.coefs))])
  
  start.t <- num.ar + 1
  i <- num.ma + 1
  x.pred <- x.lags
  e.pred <- e.lags
  
  for (t in start.t:(start.t + n.predict - 1) ) {
    e.pred[i] <- rnorm(n=1, mean=0, sd=sqrt(arma.model$sigma2))
    x.pred[t] <- intercept + as.numeric(ar.coefs %*% x.pred[(t-num.ar):(t-1)]) + as.numeric(ma.coefs %*% e.pred[(i-num.ma):(i-1)]) + e.pred[i]
    i <- i+1
  }
  
  return(x.pred[start.t:length(x.pred)])
}

# function to forecast a GARCH(1,1)
predictGarch <- function(garch.model, n.predict) {
  cond.dist <- garch.model@fit$params$cond.dist
  if (cond.dist == "norm") {
    rDist <- function() { return(rnorm(n = 1, mean = 0, sd = 1)) }
  }
  else if (cond.dist == "std") {
    cond.dist.shape <- garch.model@fit$params$shape
    rDist <- function() { return(rt(n = 1, df = cond.dist.shape)) }
  }

  model.coefs <- garch.model@fit$coef
  
  omega = model.coefs["omega"]
  alpha.coefs = rev(model.coefs[grep("alpha", names(model.coefs))])
  beta.coefs = rev(model.coefs[grep("beta", names(model.coefs))])
  
  num.alpha <- length(alpha.coefs)
  num.beta <- length(beta.coefs)
  
  start.t <- num.alpha + 1
  i <- num.beta + 1
  
  e.lags = as.vector(tail(garch.model@residuals, num.alpha))
  sigma.lags = as.vector(tail(garch.model@sigma.t, num.beta))
  
  sigma.pred <- sigma.lags
  e.pred <- e.lags
  for (t in start.t:(start.t + n.predict - 1) ) {
    sigma.pred[i] <- sqrt( omega + as.numeric(alpha.coefs %*% (e.pred[(t-num.alpha):(t-1)])^2) + as.numeric(beta.coefs %*% (sigma.pred[(i-num.beta):(i-1)])^2))
    e.pred[t] <- sigma.pred[i] * rDist()
    i <- i+1
  }
  
  return(e.pred[start.t:length(e.pred)])
}

# function to forecast a ARMA(1,1) / GARCH(1,1)
predictArmaGarch <- function(arma.garch, n.predict) {
  model.coefs <- arma.garch@fit$coef
  mu = model.coefs["mu"]
  ar.coefs = rev(model.coefs[grep("ar", names(model.coefs))])
  ma.coefs = rev(model.coefs[grep("ma", names(model.coefs))])
  
  num.ar <- length(ar.coefs)
  num.ma <- length(ma.coefs)
  
  x.lags = as.vector(tail(arma.garch@fitted, num.ar))
  e.lags = as.vector(tail(arma.garch@residuals, num.ma))
  
  start.t <- num.ar + 1
  i <- num.ma + 1
  x.pred <- x.lags
  e.pred <- c(e.lags, predictGarch(arma.garch, n.predict))
  
  for (t in start.t:(start.t + n.predict - 1) ) {
    x.pred[t] <- mu + as.numeric(ar.coefs %*% x.pred[(t-num.ar):(t-1)]) + as.numeric(ma.coefs %*% e.pred[(i-num.ma):(i-1)]) + e.pred[i]
    i <- i+1
  }
  
  return(x.pred[start.t:length(x.pred)])
}

forecastArmaGarch <- function(arma.model, garch.model, n.predict) {
  num.ar <- arma.model$arma[1]
  num.ma <- arma.model$arma[2]
  
  arma.coefs <- arma.model$coef
  intercept = arma.coefs["intercept"]
  ar.coefs = rev(arma.coefs[grep("ar", names(arma.coefs))])
  ma.coefs = rev(arma.coefs[grep("ma", names(arma.coefs))])
  
  x.lags = as.vector(tail(arma.model$fitted, num.ar))
  e.lags = as.vector(tail(arma.model$residuals, num.ma))
  
  start.t <- num.ar + 1
  i <- num.ma + 1
  x.pred <- x.lags
  e.pred <- c(e.lags, predictGarch(garch.model, n.predict))
  
  for (t in start.t:(start.t + n.predict - 1) ) {
    x.pred[t] <- intercept + as.numeric(ar.coefs %*% x.pred[(t-num.ar):(t-1)]) + as.numeric(ma.coefs %*% e.pred[(i-num.ma):(i-1)]) + e.pred[i]
    i <- i+1
  }
  
  return(x.pred[start.t:length(x.pred)])
}