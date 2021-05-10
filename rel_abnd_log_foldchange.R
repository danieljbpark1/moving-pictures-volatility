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
load("./otu_tables_MP.Rdata")


# consider only OTUs present in at least 2 samples
rel.otutab.D <- rel.otutab.D[apply(rel.otutab.D, 1, function(x) sum(x > 0)) >= 2, ]
rel.otutab.E <- rel.otutab.E[apply(rel.otutab.E, 1, function(x) sum(x > 0)) >= 2, ]
rel.otutab.F <- rel.otutab.F[apply(rel.otutab.F, 1, function(x) sum(x > 0)) >= 2, ]
rel.otutab.F4 <- rel.otutab.f4[apply(rel.otutab.f4, 1, function(x) sum(x > 0)) >= 2, ]
rel.otutab.M3 <- rel.otutab.m3[apply(rel.otutab.m3, 1, function(x) sum(x > 0)) >= 2, ]
## DATAFRAME OF LOG FOLDCHANGES WITHIN ALL SUBECTS 
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
data.list <- list(list(rel.otutab.F4, "F4"),
                  list(rel.otutab.M3, "M3"),
                  list(rel.otutab.D, "D"),
                  list(rel.otutab.E, "E"),
                  list(rel.otutab.F, "F"))


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
  quintile.data <- data.frame(otu.avg.relabs = otu.avg.relabs) %>%
    rownames_to_column(var = "otu.id") %>%
    mutate(quintile = ntile(otu.avg.relabs,
                                     n = 5))
  
  # join the two dataframes
  res.data <- left_join(x = subj.data,
                        y = quintile.data,
                        by = "otu.id") %>%
    mutate(subj.id = subj.id)
  # bind to master dataframe
  log.foldchange.data <- rbind(log.foldchange.data, res.data)
  
  ## PROP. OF NON-ZERO SAMPLES 
  # rare.otu.data <- otu.quintiles %>% 
  #   filter(quintile == 1)
  # 
  # prop.data <- apply(X = rel.otutab > 0, MARGIN = 1, mean)
  # prop.data <- data.frame(prop.nonzero = prop.data)
  # prop.data <- setDT(prop.data, keep.rownames = TRUE) %>%
  #   rename(otu.id = rn)
  # 
  # res.prop.data <- inner_join(x = prop.data,
  #                             y = rare.otu.data, 
  #                             by = "otu.id") %>%
  #   dplyr::select(-quintile) %>%
  #   mutate(subj.id = subj.id)
  # 
  # prop.positive.data <- rbind(prop.positive.data, res.prop.data)
}

log.foldchange.data <- log.foldchange.data %>%
  mutate(quintile = as.factor(quintile))

dim(log.foldchange.data)
head(log.foldchange.data, 10)
tail(log.foldchange.data, 10)

log.foldchange.plot <- ggplot(data = log.foldchange.data, aes(x=log.foldchange, fill=quintile)) +
  geom_density(aes(color=quintile), alpha=0.1) +
  facet_wrap(vars(subj.id), scales = "free", ncol = 2) +
  labs(title = "Distribution of OTU log fold-changes",
       x = "Log fold-changes",
       color = "Within-subject average relative abundance quintile",
       fill = "Within-subject average relative abundance quintile") +
  theme(legend.position = "bottom")

log.foldchange.plot

log.foldchange.data %>%
  group_by(subj.id, quintile) %>%
  summarise(sd = sd(log.foldchange)) %>%
  print(n = Inf)

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



