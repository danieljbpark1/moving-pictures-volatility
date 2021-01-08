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

# are any OTUs present in 2 or fewer samples?
any(rowSums(otu.relabs.f4 != 0) <= 2)
any(rowSums(otu.relabs.m3 != 0) <= 2)
any(rowSums(dr.D.rel.otutab != 0) <= 2)
any(rowSums(dr.E.rel.otutab != 0) <= 2)
any(rowSums(dr.F.rel.otutab != 0) <= 2)

# consider only OTUs present in at least 2 samples
otu.relabs.D <- dr.D.rel.otutab[rowSums(dr.D.rel.otutab != 0) > 2, ]
otu.relabs.E <- dr.E.rel.otutab[rowSums(dr.E.rel.otutab != 0) > 2, ]
otu.relabs.F <- dr.F.rel.otutab[rowSums(dr.F.rel.otutab != 0) > 2, ]

# calculate each OTU's average rel. abnd. within its subject
avg.otu.relabs.f4 <- apply(otu.relabs.f4, 1, mean)
avg.otu.relabs.m3 <- apply(otu.relabs.m3, 1, mean)
avg.otu.relabs.D <- apply(otu.relabs.D, 1, mean)
avg.otu.relabs.E <- apply(otu.relabs.E, 1, mean)
avg.otu.relabs.F <- apply(otu.relabs.F, 1, mean)

## DATAFRAME OF LOG FOLDCHANGES WITHIN ALL SUBJECTS 
## otu.id : the OTU id
## log.foldchange : the log foldchange of an OTU within a subject
## avg.otu.rel.abnd : average OTU relative abundance within-subject
## rel.abnd.quintile : OTU quintile based on within-subject average rel. abnd.
## subj.id : the subject id
log.foldchange.data <- data.frame(otu.id = character(),
                                  log.foldchange = numeric(),
                                  avg.otu.rel.abnd = numeric(),
                                  rel.abnd.quintile = numeric(),
                                  subj.id = character()
                                  )

## iterate thru all the data for all subjects
data.list <- list(list(otu.relabs.f4, avg.otu.relabs.f4, "F4"),
                  list(otu.relabs.m3, avg.otu.relabs.m3, "M3"),
                  list(otu.relabs.D, avg.otu.relabs.D, "D"),
                  list(otu.relabs.E, avg.otu.relabs.E, "E"),
                  list(otu.relabs.F, avg.otu.relabs.F, "F"))

for (subj in data.list) {
  # OTU rel. abnd. table
  otu.relabs.tab <- subj[[1]]
  # average OTU rel. abnd. 
  avg.otu.relabs <- subj[[2]]
  subj.id <- subj[[3]]
  
  # log foldchanges table
  log.foldchange.tab <- fold_difference_table(otu.table = otu.relabs.tab,
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
    select(rn, value) %>%
    rename(otu.id = rn,
           log.foldchange = value)
  # bin OTUs into quintiles based on avg. within-subject rel. abnd.
  otu.rel.abnd.data <- data.frame(avg.otu.rel.abnd = avg.otu.relabs)
  otu.rel.abnd.data <- otu.rel.abnd.data %>%
    rownames_to_column(var = "otu.id") %>%
    mutate(rel.abnd.quintile = ntile(avg.otu.rel.abnd, 
                                     n = 5)) 
  
  print(length(unique(subj.data$otu.id)))
  print(length(unique(otu.rel.abnd.data$otu.id)))
  # join the two dataframes
  res.data <- left_join(x = subj.data,
                        y = otu.rel.abnd.data,
                        by = "otu.id") %>%
    mutate(subj.id = subj.id)
  # bind to master dataframe
  log.foldchange.data <- rbind(log.foldchange.data, res.data)
}

log.foldchange.data <- log.foldchange.data %>%
  mutate(rel.abnd.quintile = as.factor(rel.abnd.quintile))

dim(log.foldchange.data)
head(log.foldchange.data, 10)
tail(log.foldchange.data, 10)

log.foldchange.plot <- ggplot(data = log.foldchange.data,
                              aes(x=log.foldchange,
                                  fill=rel.abnd.quintile)) +
  geom_density(alpha=0.4) +
  facet_wrap(vars(subj.id), scales = "free") +
  labs(title = "Distribution of log foldchanges",
       subtitle = "OTUs binned by average relative abundance within-subject quintiles",
       x = "log foldchanges",
       fill = "within-subject \n avg. rel. abnd. \n quintile")

log.foldchange.plot
ggsave("log_foldchanges_rel_abnd.png")







