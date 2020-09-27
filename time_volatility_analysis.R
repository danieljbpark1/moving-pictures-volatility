# consider samples at a given time interval
# calculate log(fold-change) of OTUs between pairs of consecutive samples
# calculate the quantiles of each pair
load("mp_F4_data.Rdata")
source("mp_analysis_functions.R")

# are any OTUs present in 2 or fewer samples?
any(rowSums(otu.relabs.f4 != 0) <= 2)

# the day on which the sample was observed
f4.time <- f4.metadata$days_since_experiment_start

# the times with lag between 
lag <- 7
f4.lags <- seq(min(f4.time), max(f4.time), by = lag)

# subsetting the metadata where time is in the lags
lag.metadata <- subset(f4.metadata, 
                       days_since_experiment_start %in% f4.lags)

# selecting samples from OTU table where sample time is in the lags
lag.otutab <- otu.relabs.f4[ , rownames(lag.metadata)]

# log fold changes in OTU relative abundance between consecutive samples
log.foldchange.tab <- fold_difference_table(lag.otutab, 
                                            taxa_are_rows = TRUE, 
                                            logRatio = TRUE)

big.df <- data.frame(start.time =  numeric(),
                        lag =  numeric(),
                        quantile_0.1 = numeric(),
                        quantile_0.25 = numeric(),
                        quantile_0.5 = numeric(),
                        quantile_0.75 = numeric(),
                        quantile_0.9 = numeric()
                        )

# exclude log fold changes between samples that are consecutive but not exactly lag days apart
num.samples <- nrow(lag.metadata)
time.gap <- lag.metadata$days_since_experiment_start[2:num.samples] - lag.metadata$days_since_experiment_start[1:(num.samples-1)]
true.lagpairs <- time.gap == lag

log.foldchange.tab <- log.foldchange.tab[ , true.lagpairs]
# quantiles within intervals of log fold changes of only OTUs that were present at start and end
quantile.df <- t(apply(log.foldchange.tab, 2, 
                       function(x) quantile(x[is.finite(x)], 
                                            probs = c(0.1, 0.25, 0.5, 0.75, 0.9), 
                                            na.rm = TRUE)
                       )
                 )

# the start times for pairs of consecutive samples exactly lag days apart
lag.metadata <- lag.metadata[c(true.lagpairs, FALSE), "days_since_experiment_start"]
lag.metadata$lag <- rep(lag, nrow(lag.metadata))

result.df <- cbind(lag.metadata, quantile.df)
colnames(result.df) <- c("start.time", 
                         "lag", 
                         "quantile_0.1",
                         "quantile_0.25",
                         "quantile_0.5",
                         "quantile_0.75",
                         "quantile_0.9"
                         )

big.df <- rbind(big.df, result.df)
