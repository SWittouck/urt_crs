library(tidyverse)
library(tidyamplicons)

load("data/urt.rda")

lgc_genera <- c(
  "Lactobacillus", "Pediococcus", "Leuconostoc", "Weissella", "Oenococcus",
  "Fructobacillus", "Convivina"
)

urt$taxa[urt$taxa$genus %in% lgc_genera, "family"] <- "LGC"

urt_fam <- 
  urt %>%
  filter_samples(location %in% c("N", "NF")) %>%
  aggregate_taxa(rank = "family")

save(urt_fam, file = "results/parsed/urt_fam.rda")
