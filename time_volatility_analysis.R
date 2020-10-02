# consider samples at a given time interval
# calculate log(fold-change) of OTUs between pairs of consecutive samples
# calculate the quantiles of each pair
load("mp_F4_data.Rdata")
load("mp_M3_data.Rdata")
source("mp_analysis_functions.R")
library(ggplot2)
library(gridExtra)

# are any OTUs present in 2 or fewer samples?
any(rowSums(otu.relabs.f4 != 0) <= 2)
any(rowSums(otu.relabs.m3 != 0) <= 2)

# a function to calculate the quantiles of log fold changes 
# between samples that are lag days apart
# metadata.tab: metadata containing sampID and time ordered by time ascending
# otu.tab: OTU table with taxa as rows and samples as columns ordered by time 
# lag.intervals: a vector of lags to use
# returns a dataframe with quantiles of log fold changes, 
# start time of interval, and lag of interval
lagged.foldchange <- function(metadata.tab, otu.tab, lag.intervals) {
  result.df <- data.frame(start.time =  numeric(),
                          lag =  numeric(),
                          quantile_0.1 = numeric(),
                          quantile_0.25 = numeric(),
                          quantile_0.5 = numeric(),
                          quantile_0.75 = numeric(),
                          quantile_0.9 = numeric()
                          )
  # the times at which the samples were observed
  time <- metadata.tab$days_since_experiment_start
  
  for (lag in lag.intervals) {
    # the sequence of times with lag in between 
    lag.sequence <- seq(min(time), max(time), by = lag)
    
    # subsetting the metadata where time is in the lag sequence
    lag.metadata <- subset(metadata.tab, 
                           days_since_experiment_start %in% lag.sequence)
    
    # selecting samples from OTU table where time is in the lag sequence
    lag.otutab <- otu.tab[ , rownames(lag.metadata)]
    
    # log fold changes in OTU (relative) abundance between consecutive samples
    log.foldchange.tab <- fold_difference_table(lag.otutab, 
                                                taxa_are_rows = TRUE, 
                                                logRatio = TRUE)
    
    # exclude log fold changes between samples that are consecutive but not exactly lag days apart
    num.samples <- nrow(lag.metadata)
    time.gap <- lag.metadata$days_since_experiment_start[2:num.samples] - lag.metadata$days_since_experiment_start[1:(num.samples-1)]
    true.pairs <- time.gap == lag
    
    log.foldchange.tab <- log.foldchange.tab[ , true.pairs]
    # quantiles within intervals of log fold changes of only OTUs that were present at start and end
    quantile.df <- t(apply(log.foldchange.tab, 2, 
                           function(x) quantile(x[is.finite(x)], 
                                                probs = c(0.1, 0.25, 0.5, 0.75, 0.9), 
                                                na.rm = TRUE)
                           )
                     )
    
    # the start times for pairs of consecutive samples exactly lag days apart
    lag.metadata <- lag.metadata[c(true.pairs, FALSE), "days_since_experiment_start"]
    lag.metadata$lag <- rep(lag, nrow(lag.metadata))
    
    lag.df <- cbind(lag.metadata, quantile.df)
    colnames(lag.df) <- c("start.time", 
                          "lag", 
                          "quantile_0.1",
                          "quantile_0.25",
                          "quantile_0.5",
                          "quantile_0.75",
                          "quantile_0.9"
                          )
    
    result.df <- rbind(result.df, lag.df)
  }
  return(result.df)
}

lag.quantile.f4 <- lagged.foldchange(f4.metadata, otu.relabs.f4, lag.intervals = c(3, 7, 14, 30))
lag.quantile.m3 <- lagged.foldchange(m3.metadata, otu.relabs.m3, lag.intervals = c(3, 7, 14, 30))

summary(lag.quantile.f4)
summary(lag.quantile.m3)

p.f4.1 <- ggplot(data = lag.quantile.f4, aes(x=lag, y=quantile_0.1, group=lag)) +
  geom_boxplot()
p.f4.2 <- ggplot(data = lag.quantile.f4, aes(x=lag, y=quantile_0.25, group=lag)) +
  geom_boxplot()
p.f4.3 <- ggplot(data = lag.quantile.f4, aes(x=lag, y=quantile_0.5, group=lag)) +
  geom_boxplot()
p.f4.4 <- ggplot(data = lag.quantile.f4, aes(x=lag, y=quantile_0.75, group=lag)) +
  geom_boxplot()
p.f4.5 <- ggplot(data = lag.quantile.f4, aes(x=lag, y=quantile_0.9, group=lag)) +
  geom_boxplot()
grid.arrange(p.f4.1, p.f4.2, p.f4.3, p.f4.4, p.f4.5, nrow=2,
             top = "Subject F4")

p.m3.1 <- ggplot(data = lag.quantile.m3, aes(x=lag, y=quantile_0.1, group=lag)) +
  geom_boxplot()
p.m3.2 <- ggplot(data = lag.quantile.m3, aes(x=lag, y=quantile_0.25, group=lag)) +
  geom_boxplot()
p.m3.3 <- ggplot(data = lag.quantile.m3, aes(x=lag, y=quantile_0.5, group=lag)) +
  geom_boxplot()
p.m3.4 <- ggplot(data = lag.quantile.m3, aes(x=lag, y=quantile_0.75, group=lag)) +
  geom_boxplot()
p.m3.5 <- ggplot(data = lag.quantile.m3, aes(x=lag, y=quantile_0.9, group=lag)) +
  geom_boxplot()
grid.arrange(p.m3.1, p.m3.2, p.m3.3, p.m3.4, p.m3.5, nrow=2,
             top = "Subject M3")

paired.wilcox <- function(quantile.df, lag.values, quantile.names){
  result.df <- data.frame(lag.1 = numeric(),
                          lag.2 = numeric(),
                          wilcox.pval = numeric(),
                          quantile.name = character()
                          )
  lag.pairs <- t(combn(lag.values, 2))
  n.pairs <- nrow(lag.pairs)
  for (q in quantile.names) {
    quantile.name <- rep(q, n.pairs)
    wilcox.pval <- numeric()
    for (i in 1:n.pairs) {
      wilcox.pval[i] <- wilcox.test(x = quantile.df[quantile.df$lag == lag.pairs[i,1], q],
                                    y = quantile.df[quantile.df$lag == lag.pairs[i,2], q])$p.value
    }
    wilcox.df <- data.frame(cbind(lag.pairs, wilcox.pval, quantile.name))
    colnames(wilcox.df) <- c("lag.1", "lag.2", "wilcox.pval", "quantile.name")
    result.df <- rbind(result.df, wilcox.df)
  }
  return(as.data.frame(result.df))
}

wilcox.tests.f4 <- paired.wilcox(lag.quantile.f4, 
                                 lag.values = unique(lag.quantile.f4$lag), 
                                 quantile.names = colnames(lag.quantile.f4)[c(-1,-2)])
wilcox.tests.m3 <- paired.wilcox(lag.quantile.m3,
                                lag.values = unique(lag.quantile.m3$lag),
                                quantile.names = colnames(lag.quantile.m3)[c(-1,-2)])

p.adjust(p = wilcox.tests.f4$wilcox.pval, method = "fdr")
p.adjust(p = wilcox.tests.m3$wilcox.pval, method = "fdr")

