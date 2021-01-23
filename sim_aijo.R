library(dplyr)
library(utilityf)
load("./otu_tables_MP.Rdata")
load("mp_M3_data.Rdata") # Moving Pictures subject M3

## 1. calculate proportion estimates of taxa at individual time points under multinomial model
##    for subject M3 this is stored in rel.otutab.m3 

plot(1:nrow(m3.metadata), m3.metadata$days_since_experiment_start,
     pch = 19, xlab = "sample index", ylab = "days since experiment start",
     main = "Moving Pictures Subject M3 Sample Days")

i.start <- which(m3.metadata$days_since_experiment_start == 50) # start at the sample taken on t = 50
i.end <- i.start + 199 # stop at 200th subsequent sample

## 2. filter proportion estimates using running median filter of length 15
##    filtered estimates are re-normalized to sum to 1 at each time point
filt.otutab.m3 <- matrix(nrow = nrow(rel.otutab.m3), ncol = 200)
for (i in i.start:i.end) {
  t.now <- m3.metadata$days_since_experiment_start[i] # what day is it currently
  filt.metadata <- filter(data.frame(m3.metadata), dplyr::between(days_since_experiment_start, t.now-7, t.now+7))
  filt.median <- apply(rel.otutab.m3[ ,filt.metadata$X.SampleID], 1, median) # running median 
  filt.otutab.m3[ ,i-(i.start-1)] <- filt.median / sum(filt.median) # renormalized
}

## 3. discard those taxa that are lowly abundant (average proportion is less than 5 × 10^-5) 
##    followed by a re-normalization step 
filt.otutab.m3 <- filt.otutab.m3[apply(filt.otutab.m3, 1, mean) > 5e-5, ]
filt.otutab.m3 <- apply(filt.otutab.m3, 2, function(x) x/sum(x))

## 4. transform the simplex-valued estimates to real space using the inverse softmax function
real.otutab.m3 <- apply(filt.otutab.m3, 2, function(x) reverseSoftmax(x)$par)

## 5. multiplicative Gaussian distributed noise 
##    log fold-change noise with zero-mean and standard deviation (SD) σn = 0.5
mat.noise <- matrix(exp(rnorm(n = prod(dim(real.otutab.m3)), mean = 0, sd = 0.5)), ncol = 200)
real.otutab.m3 <- real.otutab.m3 * mat.noise

##    impose zeros... HOW?

## 6. obtained noisy relative abundances by projecting the real values onto the simplex using the softmax function
soft.otutab.m3 <- apply(real.otutab.m3, 2, function(x) softmax(x))

## 7. generate noisy count data from multinomial distribution using noisy relative abundances and 
##    N_t sampled from Poisson distribution with rate lambda = 10,000
noisy.otutab.m3 <- matrix(nrow = nrow(soft.otutab.m3), ncol = 200)
for (i in 1:200) {
  N <- rpois(n = 1, lambda = 10000)
  noisy.otutab.m3[ ,i] <- rmultinom(n = 1, size = N, prob = soft.otutab.m3[ ,i])
}

# iterate
# required: a pre-existing folder at "./SimSets<n.subj>"
for (set in 1:500) {
  # save in pre-existing folder, labeled as set0001, set0002, etc. 
  write.table(noisy.otutab.m3, file = paste("./SimSets", n.subj, "/set", sprintf("%04d", set), ".txt", sep = ""), 
              sep = "\t", col.names = T, row.names = T)
}









