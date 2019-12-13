library(tidyverse)
library(tidyamplicons)

dout <- "results/coda"
dout_paper <- "results/paper_ambr2"
if (! dir.exists(dout)) dir.create(dout)
if (! dir.exists(dout_paper)) dir.create(dout_paper)

source("src/family-level/functions.R")

load("results/parsed/urt_fam.rda")

# rarefying (or not)

# urt_fam <-
#   urt_fam %>%
#   add_lib_size() %>%
#   filter_samples(lib_size >= 1000) %>%
#   rarefy(n = 1000)

# differential abundance (tidyamplicons)

urt_fam$samples %>% count(condition, location)

toplot_nose <- 
  urt_fam %>% 
  filter_samples(location == "N") %>%
  add_logratios(max_taxa = 30) %>%
  add_codifab(condition = condition)
toplot_nose %>%
  codifab_plot(CON_vs_CRS)
ggsave(
  paste0(dout, "/con_crs_nose.png"), 
  units = "cm", width = 16, height = 12
)
file.copy(
  paste0(dout, "/con_crs_nose.png"), 
  paste0(dout_paper, "/con_vs_crs_codifab_nose.png")
)

toplot_npx <- 
  urt_fam %>% 
  filter_samples(location == "NF") %>%
  add_logratios(max_taxa = 30) %>%
  add_codifab(condition = condition)
toplot_npx %>%
  codifab_plot(CON_vs_CRS)
ggsave(
  paste0(dout, "/con_crs_npx.png"), 
  units = "cm", width = 16, height = 12
)
file.copy(
  paste0(dout, "/con_crs_npx.png"), 
  paste0(dout_paper, "/con_crs_codifab_npx.png")
)

urt_fam_polyps <-
  urt_fam %>% 
  filter_samples(condition == "CRS", ! is.na(polyps)) %>%
  mutate_samples(polyps = if_else(polyps == 1, "CRSwNP", "CRSsNP"))

urt_fam_polyps$samples %>% count(location, polyps)

urt_fam_polyps %>% 
  filter_samples(location == "N") %>%
  add_logratios(max_taxa = 30) %>%
  add_codifab(condition = polyps) %>%
  codifab_plot(CRSsNP_vs_CRSwNP)
ggsave(
  paste0(dout, "/polyps_nose.png"), 
  units = "cm", width = 16, height = 12
)

urt_fam_polyps %>% 
  filter_samples(location == "NF") %>%
  add_logratios(max_taxa = 50) %>%
  add_codifab(condition = polyps) %>%
  codifab_plot(CRSsNP_vs_CRSwNP)
ggsave(
  paste0(dout, "/polyps_npx.png"), 
  units = "cm", width = 16, height = 12
)

# differential abundance (ANCOM)

ancom_nose <- 
  urt_fam %>%
  filter_samples(location == "N") %>%
  perform_ancom(condition = "condition") %>%
  left_join(urt_fam$taxa) %>%
  select(family, W_stat_nose = W_stat, W_perc_nose = W_perc) 

ancom_npx <- 
  urt_fam %>%
  filter_samples(location == "NF") %>%
  perform_ancom(condition = "condition") %>%
  left_join(urt_fam$taxa) %>%
  select(family, W_stat_npx = W_stat, W_perc_npx = W_perc) 

ancom <- left_join(ancom_nose, ancom_npx, by = "family")

write_csv(ancom, path = paste0(dout, "/families_con_crs_ancom.csv"))

# compositional pca (tidyamplicons)

urt_fam2 <- 
  urt_fam %>% 
  filter_samples(location == "N") %>% 
  add_comp_pca()

urt_fam2$samples %>%
  ggplot(aes(x = pca_1, y = pca_2, col = condition)) +
  geom_point() + 
  scale_color_brewer(palette = "Paired")
