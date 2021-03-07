## GIVEN SUBJ 
## CALC LOG FOLD-CHANGES FOR EACH OTU
## IGNORE REAPPS AND DISAPPS
## STACK LOG FOLD-CHANGES AS COL
## LABEL OTU BY REL ABND BIN (5 BINS)
## PLOT HISTOGRAM COLOR BY BIN 
library(tidyverse)
library(data.table)
source("mp_analysis_functions.R")
source("simulations/sim_functions.R") # for simulation functions
load("mp_F4_data.Rdata")
load("mp_M3_data.Rdata")
load("simulations/dethlefsen_relman_for_sim.Rdata")
load("../otu_tables_MP.Rdata")


# consider only OTUs present in at least 2 samples

# set <- sample(x = 1:10, size = 1)
#sim.dataset <- as.matrix(read.table(file = paste("./SimSets_DR_n3_t120", "/set", sprintf("%04d", set), ".txt", sep = "")))
sim.dataset <- as.matrix(read.table(file = "simset_DR_v03_0500.txt"))
dim(sim.dataset)

sim.otutab_33 <- as.data.frame(sim.dataset) %>%
  dplyr::select(contains("SUBJ_33"))
sim.otutab_77 <- as.data.frame(sim.dataset) %>%
  dplyr::select(contains("SUBJ_77"))

dim(sim.otutab_33); dim(sim.otutab_77)

summary(apply(otutab.D, 2, function(x) mean(x==0)))
summary(apply(sim.otutab_33, 2, function(x) mean(x==0)))
summary(apply(sim.otutab_77, 2, function(x) mean(x==0)))

otu.quintiles <- otu.quintiles.tab(list(list(rel.otutab.D, "D"))) # subject D is baseline for OTU quintiles
otu.quintiles <- otu.quintiles %>%
  dplyr::select(-subj.id)
head(otu.quintiles)

## DATAFRAME OF LOG FOLDCHANGES WITHIN ALL SUBJECTS 
## otu.id : the OTU id
## log.foldchange : the log foldchange of an OTU within a subject
## avg.otu.rel.abnd : average OTU relative abundance within-subject
## rel.abnd.quintile : OTU quintile based on within-subject average rel. abnd.
## subj.id : the subject id
log.foldchange.data <- data.frame(otu.id = character(),
                                  log.foldchange = numeric(),
                                  otu.avg.relabs = numeric(),
                                  rel.abnd.quintile = numeric(),
                                  subj.id = character()
                                  )

## DATAFRAME OF PROP. NON-ZERO SAMPLES
prop.positive.data <- data.frame(otu.id = character(),
                                 prop.nonzero = numeric(),
                                 otu.avg.relabs = numeric(),
                                 subj.id = character())

## iterate thru all the data for all subjects
data.list <- list(list(sim.otutab_33, "SUBJ_33"),
                  list(sim.otutab_77, "SUBJ_77"))
# data.list <- list(list(sim.otutab.m3, "M3"))
# data.list <- list(list(rel.otutab.f4, "F4"),
#                   list(rel.otutab.m3, "M3"))


for (subj in data.list) {
  rel.otutab <- subj[[1]] # OTU relative abundance table
  otu.avg.relabs <- apply(rel.otutab, 1, mean) # OTU average rel. abnd. 
  subj.id <- subj[[2]]
  
  # log foldchanges table
  log.foldchange.tab <- fold_difference_table(otutab = rel.otutab,
                                              taxa_are_rows = TRUE,
                                              logRatio = TRUE)
  # make otu.id a column
  log.foldchange.dt <- setDT(data.frame(log.foldchange.tab), 
                             keep.rownames = TRUE)
  # make long-format table by stacking row data into a column
  log.foldchange.melted <- melt(data = log.foldchange.dt,
                                id.vars = "rn")
  # only consider presence --> presence log foldchanges
  subj.data <- log.foldchange.melted %>%
    filter(is.finite(value)) %>%
    dplyr::select(rn, value) %>%
    rename(otu.id = rn,
           log.foldchange = value)
  # bin OTUs into quintiles based on avg. within-subject rel. abnd.
  # otu.rel.abnd.data <- data.frame(otu.avg.relabs = otu.avg.relabs) %>%
  #   rownames_to_column(var = "otu.id") %>%
  #   mutate(rel.abnd.quintile = ntile(otu.avg.relabs, 
  #                                    n = 5)) 
  
  # join the two dataframes
  res.data <- left_join(x = subj.data,
                        y = otu.quintiles,
                        by = "otu.id") %>%
    mutate(subj.id = subj.id)
  # bind to master dataframe
  log.foldchange.data <- rbind(log.foldchange.data, res.data)
  
  ## PROP. OF NON-ZERO SAMPLES 
  rare.otu.data <- otu.quintiles %>% 
    filter(quintile == 1)
  
  prop.data <- apply(X = rel.otutab > 0, MARGIN = 1, mean)
  prop.data <- data.frame(prop.nonzero = prop.data)
  prop.data <- setDT(prop.data, keep.rownames = TRUE) %>%
    rename(otu.id = rn)
  
  res.prop.data <- inner_join(x = prop.data,
                              y = rare.otu.data, 
                              by = "otu.id") %>%
    dplyr::select(-quintile) %>%
    mutate(subj.id = subj.id)
  
  prop.positive.data <- rbind(prop.positive.data, res.prop.data)
}

log.foldchange.data <- log.foldchange.data %>%
  mutate(quintile = as.factor(quintile))

dim(log.foldchange.data)
head(log.foldchange.data, 10)
tail(log.foldchange.data, 10)

log.foldchange.plot <- ggplot(data = log.foldchange.data, aes(x=log.foldchange, fill=quintile)) +
  geom_density(aes(color=quintile), alpha=0.1) +
  facet_wrap(vars(subj.id), scales = "free") +
  labs(title = "Distribution of OTU log fold-changes",
       #subtitle = "Observed Dethlefsen-Relman dataset",
       subtitle = paste("v03 Simulated Dethlefsen-Relman dataset set", sprintf("%04d", 500), sep = ""),
       x = "log fold-changes",
       color = "within-subject \n average \n relative \n abundance \n quintile",
       fill = "within-subject \n average \n relative \n abundance \n quintile")

log.foldchange.plot

log.foldchange.data %>%
  group_by(subj.id, quintile) %>%
  summarise(sd = sd(log.foldchange))

ggsave("sim_log_foldchanges_MP.png")

summary(prop.positive.data)

prop.positive.plot <- ggplot(data = prop.positive.data,
                             aes(x=prop.nonzero)) +
  geom_histogram(bins = 15) +
  facet_wrap(vars(subj.id), scales = "free") +
  labs(title = "Distributions of proportion of non-zero abundances across samples",
       subtitle = "OTUs in 1st quintile for rarity within each subject",
       x = "proportion of non-zero abundances across all samples within subject")

prop.positive.plot
ggsave("prop_nonzero_rare_otus.png")



