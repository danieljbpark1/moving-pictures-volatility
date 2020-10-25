library(dplyr)
library(zoib)
library(doParallel)

source("mp_analysis_functions.R")

# <subjID>.metadata : metadata sorted by time
# otu.counts.<subjID> : OTU counts table
# otu.relabs.<subjID> : OTU relative abundances table
# otu.clrabs.<subjID> : OTU clr-transformed abundances
# median.otu.relabs.<subjID> : OTU's median relative abundance across all samples
# median.otu.clrabs.<subjID> : OTU's median clr-transformed abundance across all samples
# coverage.<subjID> : the total number of reads by sample
load("mp_F4_data.Rdata")
load("mp_M3_data.Rdata")

# format the reappearance dataframe where y = relative abundance after an absence
reapp.relabs.f4 <- format_reappearance_data(otutab = otu.relabs.f4, 
                                           otu.median.abs = median.otu.clrabs.f4,
                                           read.depths = coverage.f4,
                                           otu_id = names(median.otu.clrabs.f4))

# remove outlier reappeared relative abundances from subject F4 data
new.reapp.relabs.f4 <- subset(reapp.relabs.f4, y < 0.04)

reapp.relabs.m3 <- format_reappearance_data(otutab = otu.relabs.m3,
                                            otu.median.abs = median.otu.clrabs.m3,
                                            read.depths = coverage.m3,
                                            otu_id = names(median.otu.clrabs.m3))

df.list <- list(new.reapp.relabs.f4, reapp.relabs.m3)

# mean ~ median.abundance + log(sample.read.depth)
# dispersion ~ median.abundance 
# Pr(Y = 0) ~ 1
fit.zoib <- function(d) {
  return(
    zoib(model = y ~ median.abundance + log(sample.read.depth) | median.abundance | 1, 
         data = d, joint = FALSE, n.response = 1, 
         zero.inflation = TRUE, one.inflation = FALSE, 
         EUID = 1:nrow(d),
         n.iter = 2000)
  )
}

# PARALLEL PROCESSING
registerDoParallel(cores=2)  

# returns list of zero-inflated beta regression models 
zoib.list <- foreach(df = df.list) %dopar% fit.zoib(df)

f4.reapp.relabs.ZIBR <- zoib.list[[1]]
m3.reapp.relabs.ZIBR <- zoib.list[[2]]
save(f4.reapp.relabs.ZIBR, file = "f4_reapp_rel_ab_ZIBR.Rda")
save(m3.reapp.relabs.ZIBR, file = "m3_reapp_rel_ab_ZIBR.Rda")


