library(dplyr)
library(lme4)
library(lmerTest)
library(pldist)
library(MiRKAT)
library(vegan)
source("./mp_analysis_functions.R")


halfvarson.metadata <- read.table(file = "Halfvarson 2017 IBS data/Halfvarson_2017_Metadata.txt", 
                                  header = T, sep = "\t")

# this is raw abundances table
halfvarson.otutab <- read.table(file = "Halfvarson 2017 IBS data/Halfvarson_2017_OtuTable.txt",
                                header = T, sep = "\t")

# choose subjects with 4 or more timepoints
subjID <- halfvarson.metadata %>%
  group_by(patientnumber) %>%
  summarise(n = n()) %>%
  filter(n >= 4) %>%
  pull(patientnumber)

# the new subjects' metadata
new.metadata <- halfvarson.metadata %>%
  filter(patientnumber %in% subjID) 

# the new subjects' samples and only OTUs that appear 2 or more times globally
new.otutab <- halfvarson.otutab[apply(halfvarson.otutab, 1, function(x) sum(x > 0)) >= 2, new.metadata$sample_name]

# filter samples with fewer than 10,000 reads first
new.otutab <- new.otutab[ , apply(new.otutab, 2, sum) > 10000]
# relative abundances
rel.otutab <- transform(otu_table(new.otutab, taxa_are_rows = TRUE), transform = "compositional")

# update metadata
new.metadata <- new.metadata %>%
  filter(sample_name %in% colnames(new.otutab))

# QUAL LMVoltest

otu.avg.abnd <- apply(rel.otutab, 1, mean)       # global OTU average relative abundances
otu.quintiles <- data.frame(otu.avg.abnd) %>%     # abundance-based OTU quintiles
  rownames_to_column(var = "otu.id") %>%
  mutate(quintile = ntile(otu.avg.abnd, n = 5)) 

n.totalsubj <- length(subjID) # total number of subsampled subjects

input_1 <- list() # group 1: otutab and subject ID
input_2 <- list() # group 2
for (i in 1:n.totalsubj) {
  subj.id <- subjID[i]
  subj.metadata <- new.metadata %>%
    filter(patientnumber == subj.id) %>%
    mutate(date = as.Date(collection_timestamp, format = "%m/%d/%Y")) %>%
    arrange(date)
    
  subj.otutab <- as.data.frame(rel.otutab[ , subj.metadata$sample_name])
  if("HC" %in% subj.metadata$ibd_subtype) {
    input_1[[length(input_1)+1]] <- list(subj.otutab, subj.id)
  }
  else {
    input_2[[length(input_2)+1]] <- list(subj.otutab, subj.id)
  } 
}

big.statechange.tab_1 <- statechange_prop_table(input_1) %>%
  mutate(group2 = 0)

big.statechange.tab_2 <- statechange_prop_table(input_2) %>%
  mutate(group2 = 1)

big.statechange.tab <- rbind(big.statechange.tab_1, big.statechange.tab_2)
big.statechange.tab <- left_join(x = big.statechange.tab, y = otu.quintiles, by = "otu.id")
big.statechange.tab$quintile <- as.factor(big.statechange.tab$quintile)

glmm <- glmer(prop.statechange ~ 1 + group2 + quintile + (1 | subj.id) + (1 | otu.id),
              data = big.statechange.tab, 
              family = binomial, 
              weights = timepoints)
glmm.summary <- summary(glmm)
glmm.summary$coefficients

ss <- getME(glmm, c("theta", "fixef"))
glmm2 <- update(glmm, start = ss, control = glmerControl(optCtrl = list(maxfun = 2e4)))


# QUANT LMVoltest
big.logFC.tab <- data.frame(logFC = numeric(),
                            group2 = numeric(),
                            subj.id = character())
for (i in 1:n.totalsubj) {
  subj.id <- subjID[i]
  subj.metadata <- new.metadata %>%
    filter(patientnumber == subj.id) %>%
    arrange(timepoint)
  subj.otutab <- as.data.frame(new.otutab[ , subj.metadata$sample_name])
  sim.logFC.tab <- fold_difference_table(subj.otutab, logRatio = TRUE)
  # within-sample squared L2 norm of log fold-changes excluding dis/reappearances and continued absences
  logFC_i <- apply(sim.logFC.tab, 2, FUN = function(x) mean((x[is.finite(x)])**2, na.rm = TRUE))
  
  if("HC" %in% subj.metadata$ibd_subtype) {
    subj.df <- data.frame(logFC = logFC_i,
                          group2 = rep(0, length(logFC_i)),
                          subj.id = rep(subj.id, length(logFC_i)))
  }
  else {
    subj.df <- data.frame(logFC = logFC_i,
                          group2 = rep(1, length(logFC_i)),
                          subj.id = rep(subj.id, length(logFC_i)))
    
  }
  big.logFC.tab <- rbind(big.logFC.tab, subj.df)
}

lmm <- lmer(logFC ~ 1 + group2 + (1 | subj.id), data = big.logFC.tab)
sink("halfvarson-lmvoltest-quant.txt")
print(summary(lmm))
sink()

## PLDIST
otus <- t(rel.otutab)
metadata <- data.frame(subjID = new.metadata$patientnumber, 
                       sampID = new.metadata$sample_name, 
                       time = new.metadata$timepoint)

dp <- data_prep(otus = otus, metadata = metadata, paired = FALSE)
pt <- pltransform(otu.data = dp, paired = FALSE)

Ds <- pldist_all(otus = otus, metadata = metadata, paired = FALSE, method = "bray")
Ks <- lapply(Ds, FUN = function(d) D2K(d))

res <- MiRKAT(y = y, Ks = Ks)

## NAIVE Bray 

