library(fGarch)
library(tseries)
library(forecast)
library(dplyr)
library(tidyverse)
library(car)
source("mp_analysis_functions.R")

# <subjID>.metadata : metadata sorted by time
# otu.counts.<subjID> : OTU counts table
# otu.relabs.<subjID> : OTU relative abundances table
# otu.clrabs.<subjID> : OTU clr-transformed abundances
# median.otu.relabs.<subjID> : OTU's median relative abundance across all samples
# median.otu.clrabs.<subjID> : OTU's median clr-transformed abundance across all samples
# coverage.<subjID> : the total number of reads by sample
load("~/Downloads/mp_F4_data.Rdata")
load("~/Downloads/mp_M3_data.Rdata")

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

# function to calculate MSE
mse <- function(v.1, v.2) {
  n <- length(v.1)
  return(sum((v.1-v.2)^2)/n)
}

# function to calculate MSE of time series model on validation set
evaluate_model <- function(data, model) {
  fitted.model <- model(data)
  if (class(fitted.model) == "fGARCH") {
    return(mse(data, fitted.model@fitted))
  }
  else {
    return(mse(data, fitted.model$fitted))
  }
}

validate_model <- function(ts.data, model.function) {
  n.obs <- length(ts.data)
  # split train : validation by 2 : 1
  n.train <- n.obs %/% 3 * 2
  n.validation <- n.obs - n.train
  
  ts.train <- ts.data[1:n.train]
  ts.validation <- ts.data[(n.train+1):n.obs]
  
  trained.model <- model.function(ts.train)
  
  if (class(trained.model) == "list") {
    trained.arma <- trained.model[[1]]
    trained.garch <- trained.model[[2]]
    predictions <- forecastArmaGarch(trained.arma, trained.garch, n.validation)
  }
  else {
    predictions <- predictArma(arma.model = trained.model, n.predict = n.validation)
  }
  
  return(mse(ts.validation, predictions))
}

# format dataframe of consecutive fold changes for subject F4 and M3
f4.fold.change.df <- format_fold_change_data(otu.relabs.f4)
m3.fold.change.df <- format_fold_change_data(otu.relabs.m3)

fit.model.1 <- function(time.series) { return(garchFit(formula = ~ arma(1, 1) + garch(1, 1), data = time.series, cond.dist = "std", include.shape = TRUE, trace = FALSE)) }
fit.model.2 <- function(time.series) { return(Arima(time.series, order = c(1,0,1))) }
fit.model.3 <- function(time.series) { return(Arima(time.series, order = c(1,0,0))) }

# TEST VALIDATION MSE
model.function.1 <- function(time.series) {
  arma.model <- Arima(time.series, order = c(1,0,1), include.mean = TRUE, method="ML")
  garch.model <- garchFit(formula = ~ garch(1,1), data = arma.model$residuals, cond.dist = "std", include.shape = TRUE, trace = FALSE)
  return(list(arma.model, garch.model))
}

otu.list <- f4.fold.change.df %>% distinct(otu.id)

for (i in 1:nrow(otu.list)) {
  otu <- otu.list[i,1]
  print(otu)
  
  df <- subset(f4.fold.change.df, otu.id == otu)
  print(validate_model(df$log.fold.change))
}

model.1.validation <- as.vector(unlist(by(f4.fold.change.df, f4.fold.change.df$otu.id, function(otu.data) validate_model(otu.data$log.fold.change, model.function.1))))

# dataframes of MSE across models
f4.mse.df <- data.frame(otu.id <- as.vector(f4.fold.change.df %>% distinct(otu.id)),
                        mse.model.1 = as.vector(unlist(by(f4.fold.change.df, f4.fold.change.df$otu.id, function(otu.data) evaluate_model(otu.data$log.fold.change, fit.model.1)))),
                        mse.model.2 = as.vector(unlist(by(f4.fold.change.df, f4.fold.change.df$otu.id, function(otu.data) evaluate_model(otu.data$log.fold.change, fit.model.2)))),
                        mse.model.3 = as.vector(unlist(by(f4.fold.change.df, f4.fold.change.df$otu.id, function(otu.data) evaluate_model(otu.data$log.fold.change, fit.model.3))))
                        )
save(f4.mse.df, "f4_mse_df.Rdata")

m3.mse.df <- data.frame(otu.id <- as.vector(m3.fold.change.df %>% distinct(otu.id)),
                        mse.model.1 = as.vector(unlist(by(m3.fold.change.df, m3.fold.change.df$otu.id, function(otu.data) evaluate_model(otu.data$log.fold.change, fit.model.1)))),
                        mse.model.2 = as.vector(unlist(by(m3.fold.change.df, m3.fold.change.df$otu.id, function(otu.data) evaluate_model(otu.data$log.fold.change, fit.model.2)))),
                        mse.model.3 = as.vector(unlist(by(m3.fold.change.df, m3.fold.change.df$otu.id, function(otu.data) evaluate_model(otu.data$log.fold.change, fit.model.3))))
                        )
save(m3.mse.df, "m3_mse_df.Rdata")

# Model 1 and 2 have lower median and mean training MSE than Model 3
summary(f4.mse.df)

# same trend for subject M3
summary(m3.mse.df)

f4.mse.long <- f4.mse.df %>% 
  pivot_longer(., cols = c("mse.model.1", "mse.model.2", "mse.model.3"),
               names_to = "model", values_to = "training.mse")

f4.mse.long %>%
  ggplot(aes(x = model, y = training.mse, fill = model)) +
  geom_boxplot() +
  labs(title = "Subject F4 time series models", y = "Training MSE for each OTU") +
  scale_fill_discrete(name = "model", labels = c("ARMA-GARCH(1,1,1,1)", "ARMA(1,1)", "AR(1)"))

m3.mse.long <- m3.mse.df %>% 
  pivot_longer(., cols = c("mse.model.1", "mse.model.2", "mse.model.3"),
               names_to = "model", values_to = "training.mse") 

m3.mse.long %>%
  ggplot(aes(x = model, y = training.mse, fill = model)) +
  geom_boxplot() +
  labs(title = "Subject M3 time series models", y = "Training MSE for each OTU") +
  scale_fill_discrete(name = "model", labels = c("ARMA-GARCH(1,1,1,1)", "ARMA(1,1)", "AR(1)"))


res.aov.f4 <- aov(training.mse ~ model, data = f4.mse.long)
# no pattern between residuals and fitted values, but there are outliers
plot(res.aov.f4, 1)
# p-value is not significant, therefore we can assume homogeneity of group variances
leveneTest(training.mse ~ as.factor(model), data = f4.mse.long)
# normality of residuals not met; heavily right skewed
plot(res.aov.f4, 2)
# ANOVA test for F4 suggests there is a difference in models' mean training MSE
summary(res.aov.f4)

res.aov.m3 <- aov(training.mse ~ model, data = m3.mse.long)
# no pattern between residuals and fitted values, but there are outliers
plot(res.aov.m3, 1)
# p-value is not significant, therefore we can assume homogeneity of group variances
leveneTest(training.mse ~ as.factor(model), data = m3.mse.long)
# normality of residuals not met; heavily right skewed
plot(res.aov.m3, 2)
# no significant group differences for M3
summary(res.aov.m3)

### LOG transformation of MSE
log.res.aov.f4 <- aov(log(training.mse) ~ model, data = f4.mse.long)
# no pattern between residuals and fitted values
plot(log.res.aov.f4, 1)
# homogeneity of group variances
leveneTest(log(training.mse) ~ as.factor(model), data = f4.mse.long)
# more symmetric with log transformation, but only approximately normal
plot(log.res.aov.f4, 2)
# significant ANOVA p-value
summary(log.res.aov.f4)

# Model 3 - Model 1 and Model 3 - Model 2 are significant differences
# this suggests the AR(1) model has significantly higher mean training MSE than the other models
TukeyHSD(log.res.aov.f4)

log.res.aov.m3 <- aov(log(training.mse) ~ model, data = m3.mse.long)
# no pattern between residuals and fitted values
plot(log.res.aov.m3, 1)
# homogeneity of group variances
leveneTest(log(training.mse) ~ as.factor(model), data = m3.mse.long)
# still right-skewed residuals
plot(log.res.aov.m3, 2)
# significant ANOVA p-value
summary(log.res.aov.m3)

### NON-PARAMETRIC TESTS

# non-parametric Kruskal-Wallis rank sum test produces significant p-value
# which means there still is a difference in the group means without ANOVA assumptions
kruskal.test(training.mse ~ model, data = f4.mse.long)

kruskal.test(training.mse ~ model, data = m3.mse.long)

# Only Model 1 - Model 3 and Model 2 - Model 3 are significantly different
pairwise.wilcox.test(f4.mse.long$training.mse, f4.mse.long$model)

# Only Model 1 - Model 3 and Model 2 - Model 3 are significantly different
pairwise.wilcox.test(m3.mse.long$training.mse, m3.mse.long$model)
