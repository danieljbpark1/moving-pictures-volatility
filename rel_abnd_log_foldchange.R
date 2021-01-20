## GIVEN SUBJ 
## CALC LOG FOLD-CHANGES FOR EACH OTU
## IGNORE REAPPS AND DISAPPS
## STACK LOG FOLD-CHANGES AS COL
## LABEL OTU BY REL ABND BIN (5 BINS)
## PLOT HISTOGRAM COLOR BY BIN 
library(tidyverse)
library(data.table)
source("mp_analysis_functions.R")
load("mp_F4_data.Rdata")
load("mp_M3_data.Rdata")
load("dethlefsen_relman.Rdata")
load("../otu_tables_MP.Rdata")

# consider only OTUs present in at least 2 samples
otu.relabs.D <- dr.D.rel.otutab[rowSums(dr.D.rel.otutab != 0) > 2, ]
otu.relabs.E <- dr.E.rel.otutab[rowSums(dr.E.rel.otutab != 0) > 2, ]
otu.relabs.F <- dr.F.rel.otutab[rowSums(dr.F.rel.otutab != 0) > 2, ]

set <- 358
sim.dataset <- as.matrix(read.table(file = paste("./SimSets_MP_n2_t120", "/set", sprintf("%04d", set), ".txt", sep = "")))

sim.otutab.f4 <- as.data.frame(sim.dataset) %>%
  dplyr::select(contains("F4"))
sim.otutab.m3 <- as.data.frame(sim.dataset) %>%
  dplyr::select(contains("M3"))

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
data.list <- list(list(sim.otutab.f4, "F4"),
                  list(sim.otutab.m3, "M3"))

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
  otu.rel.abnd.data <- data.frame(otu.avg.relabs = otu.avg.relabs) %>%
    rownames_to_column(var = "otu.id") %>%
    mutate(rel.abnd.quintile = ntile(otu.avg.relabs, 
                                     n = 5)) 
  
  # join the two dataframes
  res.data <- left_join(x = subj.data,
                        y = otu.rel.abnd.data,
                        by = "otu.id") %>%
    mutate(subj.id = subj.id)
  # bind to master dataframe
  log.foldchange.data <- rbind(log.foldchange.data, res.data)
  
  ## PROP. OF NON-ZERO SAMPLES 
  rare.otu.data <- otu.rel.abnd.data %>% 
    filter(rel.abnd.quintile == 1)
  
  prop.data <- apply(X = rel.otutab > 0, MARGIN = 1, mean)
  prop.data <- data.frame(prop.nonzero = prop.data)
  prop.data <- setDT(prop.data, keep.rownames = TRUE) %>%
    rename(otu.id = rn)
  
  res.prop.data <- inner_join(x = prop.data,
                              y = rare.otu.data, 
                              by = "otu.id") %>%
    dplyr::select(-rel.abnd.quintile) %>%
    mutate(subj.id = subj.id)
  
  prop.positive.data <- rbind(prop.positive.data, res.prop.data)
}

log.foldchange.data <- log.foldchange.data %>%
  mutate(rel.abnd.quintile = as.factor(rel.abnd.quintile))

dim(log.foldchange.data)
head(log.foldchange.data, 10)
tail(log.foldchange.data, 10)

log.foldchange.plot <- ggplot(data = log.foldchange.data,
                              aes(x=log.foldchange,
                                  fill=rel.abnd.quintile)) +
  geom_density(aes(color=rel.abnd.quintile, 
                   ), alpha=0.1) +
  facet_wrap(vars(subj.id), scales = "free") +
  xlim(-8,8) +
  labs(title = "Distribution of OTU log fold-changes",
       # subtitle = "Observed Moving Pictures dataset",
       subtitle = paste("Simulated Moving Pictures dataset set", sprintf("%04d", set), sep = ""),
       x = "log fold-changes",
       color = "within-subject \n average \n relative \n abundance \n quintile",
       fill = "within-subject \n average \n relative \n abundance \n quintile")

log.foldchange.plot
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



