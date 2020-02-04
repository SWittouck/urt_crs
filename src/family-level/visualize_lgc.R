library(tidyverse)
library(tidyamplicons)
library(ggforce)
library(ggpubr)

dout <- "results/LGC"
dout_paper <- "results/paper_ambr2"
if (! dir.exists(dout)) dir.create(dout)
if (! dir.exists(dout_paper)) dir.create(dout_paper)

load("results/parsed/urt_fam.rda")

# rarefying (or not)

# urt_fam <-
#   urt_fam %>%
#   add_lib_size() %>%
#   filter_samples(lib_size >= 1000) %>%
#   rarefy(n = 1000)

# CON vs CRS

samples <- 
  urt_fam %>%
  mutate_samples(
    location = recode(location, "N" = "Anterior nares", "NF" = "Nasopharynx")
  ) %>%
  add_lib_size() %>%
  filter_taxa(family == "LGC") %>%
  {left_join(.$samples, .$abundances)} %>%
  replace_na(replace = list(abundance = 0)) %>%
  mutate(
    present = if_else(abundance == 0, "no", "yes"),
    abundance = abundance + 1, 
    rel_abundance = abundance / (lib_size + 1)
  ) %>%
  mutate_at("present", factor, levels = c("yes", "no"))

(p1 <- samples %>%
  ggplot(aes(x = condition, y = rel_abundance)) +
  geom_violin(draw_quantiles = c(0.5), col = "darkgrey", fill = "#a6cee3") +
  geom_point(size = 0.5) + 
  facet_wrap(~ location, scales = "free_y") +
  theme_classic())
ggsave(
  paste0(dout, "/con_vs_crs_violin.png"), width = 8, height = 8, units = "cm"
)

(p2 <- samples %>%
  ggplot(aes(x = condition, y = rel_abundance)) +
  geom_sina(size = 0.2, col = "grey") +
  stat_summary(fun.y = "median", geom = "point", shape = 24) + 
  facet_wrap(~ location, scales = "free_y") +
  theme_classic())
ggsave(
  paste0(dout, "/con_vs_crs_sina.png"), width = 8, height = 8, units = "cm"
)

(p3 <- samples %>%
  ggplot(aes(x = condition, y = rel_abundance)) +
  geom_violin(draw_quantiles = c(0.5), fill = "#a6cee3", size = 0.3) +
  geom_point(size = 0.5) + 
  facet_wrap(~ location, scales = "free_y") +
  scale_y_log10() +
  theme_classic())
ggsave(
  paste0(dout, "/con_vs_crs_log_violin.png"), 
  width = 8, height = 8, units = "cm"
)

(p4 <- samples %>% 
  mutate_at("present", factor, levels = c("no", "yes")) %>%
  ggplot(aes(x = condition, y = rel_abundance)) +
  geom_boxplot(col = "grey", outlier.alpha = 0) +
  geom_sina(aes(col = present), size = 0.3) + 
  stat_compare_means(
    method = "wilcox.test", paired = F, 
    comparisons = list(c(1, 2)), size = 2
  ) +
  facet_wrap(~ location, scales = "free_y") +
  ylab("relative abundance") +
  scale_y_log10() +
  scale_color_manual(values = c("yes" = "#1f78b4", "no" = "#a6cee3")) +
  theme_classic()) +
  theme(
    legend.position = "bottom",
    text = element_text(size = 8)
  )
ggsave(
  paste0(dout, "/con_vs_crs_log_sina.png"), 
  width = 8, height = 8, units = "cm"
)
file.copy(
  paste0(dout, "/con_vs_crs_log_sina.png"), 
  paste0(dout_paper, "/con_vs_crs_lgc.png")
)

ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D"))
ggsave(
  paste0(dout, "/con_vs_crs.png"), width = 16, height = 16, units = "cm"
)

# CON vs CRS - presence/absence 

samples %>% 
  mutate_at("present", factor, levels = c("no", "yes")) %>%
  count(condition, location, present, name = "count") %>%
  group_by(condition, location) %>%
  mutate(percentage = count / sum(count)) %>%
  ungroup() %>%
  ggplot(aes(x = condition, y = percentage, fill = present)) +
  geom_col() + 
  facet_wrap(~ location, scales = "free_y") +
  scale_fill_manual(values = c("yes" = "#1f78b4", "no" = "#a6cee3")) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    text = element_text(size = 8)
  )
ggsave(
  paste0(dout_paper, "/con_vs_crs_lgc_presabs.png"), 
  width = 8, height = 8, units = "cm"
)

# CON vs CRSwNP vs CRSsNP

samples_polyps <-
  samples %>%
  filter(! (is.na(polyps) & condition == "CRS")) %>%
  mutate(condition = case_when(
    condition == "CON" ~ "CON",
    polyps == 1 ~ "CRSwNP",
    polyps == 0 ~ "CRSsNP"
  )) 

samples_polyps %>%
  ggplot(aes(x = condition, y = rel_abundance)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = "#a6cee3", size = 0.3) +
  geom_point(size = 0.5) + 
  facet_wrap(~ location, scales = "free_y") +
  scale_y_log10() +
  theme_classic()
ggsave(
  paste0(dout, "/con_crswnp_crssnp_log_violin.png"), 
  width = 16, height = 16, units = "cm"
)

samples_polyps %>%
  filter(! is.na(condition)) %>%
  ggplot(aes(x = condition, y = rel_abundance)) +
  geom_boxplot(col = "grey", outlier.alpha = 0) +
  geom_sina(size = 0.3, col = "black") +
  facet_wrap(~ location, scales = "free_y") +
  scale_y_log10() +
  theme_classic() +
  ylab("relative abundance (log)")
ggsave(
  paste0(dout, "/con_crswnp_crssnp_log_sina.png"), 
  width = 16, height = 16, units = "cm"
)
