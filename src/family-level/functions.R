source("src/family-level/ancom/ANCOM_updated_code.R")

#' Perform ANCOM analysis
#'
#' This function performs an "analysis of composition of microbes" (ANCOM), as
#' introduced by Mandal et al. (2015).
#' 
#' Beware! The ANCOM method adds pseudocounts of one to all abundances!! In my
#' opinion, this indroduces a lot of noise in the data, unless you rarefy (but
#' then you lose information).
#' 
#' W is defined as the number of taxa with respect to which a given taxon is
#' differentially abundant between the conditions. For example, if taxon_32 has
#' a W of 40, this means that it is significantly more/less abundant than 40
#' other taxa, when compared between the conditions.
#' 
#' Detected_0.9 means that W was equal or greater than 90% of the number of taxa
#' included in the analysis (taxa with too many zero counts are excluded from
#' the analysis). Similar for detected_0.8, detected_0.7 and detected_0.6. 
#'
#' @param ta A tidyamplicons object
#' @param condition A categorical variable in the sample table; should be quoted
#' @param covariates A character vector of covariates in the sample table
#'
#' @return A tibble with the ANCOM results for all taxa
#'
#' @references \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450248/}
#'   \url{https://sites.google.com/site/siddharthamandal1985/research}
perform_ancom <- function(ta, condition, covariates = NULL) {
  
  adjusted <- ! is.null(covariates)
  if (adjusted) {
    adj_formula <- str_c(covariates, collapse = " + ")
  } else {
    adj_formula <- NULL
  }
  
  taxon_table <- 
    as_abundances_matrix(ta$abundances) %>%
    as_tibble(rownames = "Sample.ID")
  
  metadata <-
    tibble(sample_id = taxon_table$Sample.ID) %>%
    left_join(ta$samples) %>%
    rename(Sample.ID = sample_id)
  
  comparison_test <- ANCOM.main(
    OTUdat = taxon_table, 
    Vardat = metadata[c("Sample.ID", condition, covariates)],
    adjusted = adjusted, repeated = F, main.var = "condition", 
    adj.formula = adj_formula, repeat.var = NULL, longitudinal = F, 
    random.formula = NULL, multcorr = 2, sig = 0.05, prev.cut = 0.90
  )
  
  n_taxa_included <-
    ta$abundances %>%
    count(taxon_id) %>%
    {sum(.$n >= 0.10 * nrow(ta$samples))}
  
  comparison_test$W.taxa %>%
    select(taxon_id = otu.names, W_stat) %>%
    mutate(W_perc = W_stat / n_taxa_included)
  
}

#' Add logratios
#'
#' This function computes pairwise logratio values between all taxa and adds
#' these to the tidyamplicons object in the form of a table called logratios.
#' 
#' IMPORTANT: this function add pseudocounts of one to all abundances. 
#'
#' @param ta A tidyamplicons object
#' @param condition A binary variable in the sample table (unquoted)
#'
#' @return A tidyamplicons object with an extra table logratios
add_logratios <- function(ta, max_taxa = 30) {
  
  if (nrow(ta$taxa) > max_taxa) {
    
    ta <- ta %>% add_occurrences() 
    
    ta$taxa <-
      ta$taxa %>%
      arrange(desc(occurrence)) %>%
      mutate(keep = F) %>%
      {.$keep[1:max_taxa] <- T; .}
    
    ta <- ta %>% filter_taxa(keep)
    
  }
  
  abundances_complete <- 
    ta$abundances %>%
    complete(sample_id, taxon_id, fill = list(abundance = 0))
  
  ta$logratios <-
    full_join(
      abundances_complete %>% 
        select(sample_id, taxon_id, abundance),
      abundances_complete %>% 
        select(sample_id, ref_taxon_id = taxon_id, ref_abundance = abundance),
      by = "sample_id"
    ) %>%
    mutate(
      taxon_ids = str_c(taxon_id, ref_taxon_id, sep = "_"),
      logratio = log10((abundance + 1) / (ref_abundance + 1))
    ) %>%
    select(- abundance, - ref_abundance)
  
  ta
  
}

#' Perform compositional differential abundance analysis
#'
#' This function performs a differential abundance test for all pairwise ratios
#' between taxa. 
#' 
#' A table called taxon_pairs will be added to the tidyamplicons object, with
#' for each pair of a taxon and a reference taxon, the differential abundance of
#' the taxon between the two conditions (with respect to the reference taxon).
#' The test that is performed is a Wilcoxon rank sum test and the test statistic
#' that is reported is the two-sample Hodges–Lehmann estimator (the median of
#' all pairwise differences between the samples).
#'
#' @param ta A tidyamplicons object
#' @param condition A binary variable in the sample table (unquoted)
#'
#' @return A tidyamplicons object with an extra table taxon_pairs
add_codifab <- function(ta, condition) {
  
  condition <- rlang::enexpr(condition)
  
  # if logratios not present: add and remove on exit
  if (! "logratios" %in% names(ta)) {
    ta <- add_logratios(ta)
  }
  
  samples <- ta$samples %>% mutate(condition = !! condition)
  
  conditions <- unique(samples$condition)
  
  a_vs_b <- paste0(conditions[1], "_vs_", conditions[2])
  
  ta$taxon_pairs <-
    ta$logratios %>%
    filter(taxon_id != ref_taxon_id) %>%
    left_join(samples) %>%
    group_by(taxon_ids, taxon_id, ref_taxon_id) %>%
    summarize(
      wilcox = list(wilcox.test(
        logratio ~ condition, data = tibble(logratio, condition), conf.int = T, 
        exact = F
      )),
      a_vs_b = map_dbl(wilcox, ~ .[["estimate"]]),
      wilcox_p = map_dbl(wilcox, ~ .[["p.value"]])
    ) %>%
    ungroup() %>%
    mutate(a_vs_b = 10 ^ a_vs_b) %>%
    rename(!! a_vs_b := a_vs_b)
  
  ta
  
}

add_comp_pca <- function(ta) {
  
  # if logratios not present: add and remove on exit
  if (! "logratios" %in% names(ta)) {
    ta <- add_logratios(ta)
    on.exit(ta$logratios <- NULL)
  }
  
  logratio_matrix <-
    ta$logratios %>%
    select(taxon_ids, sample_id, logratio) %>%
    spread(key = taxon_ids, value = logratio) %>%
    {
      m <- as.matrix(.[, -1])
      row.names(m) <- .$sample_id
      m
    }
  
  pca <- prcomp(logratio_matrix[, colSums(logratio_matrix) != 0], scale. = T)
  samples_pca <- tibble(
    sample_id = rownames(pca$x),
    pca_1 = unname(pca$x[, 1]),
    pca_2 = unname(pca$x[, 2])
  )
  
  # add PCoA dimensions to sample table
  ta$samples <- ta$samples %>%
    left_join(samples_pca)
  
  ta
  
}

#' Generate a compositional differential abundance plot
#'
#' This function returns a plot to visualize differential abundance of taxa
#' between conditions, compared to all other taxa as references. These
#' differential abundances should already have been calculated with
#' [add_codifab].
#' 
#' Significance of tests is determined by capping the false discovery rate at
#' 10%, using the method of Benjamini and Yekutieli, which is developed for
#' non-independent tests. See [p.adjust].
#'
#' @param ta A tidyamplicons object
#' @param diffabun_var The variable with differential abundances in the
#'   taxon_pair table
#'
#' @return A ggplot object
codifab_plot <- function(ta, diffabun_var) {
  
  diffabun_var <- rlang::enexpr(diffabun_var)
  
  taxon_pairs <-
    ta$taxon_pairs %>%
    mutate(wilcox_p = p.adjust(wilcox_p, "BY")) %>%
    left_join(
      ta$taxa %>% select(taxon_id = taxon_id, family = family)
    ) %>%
    left_join(
      ta$taxa %>% select(ref_taxon_id = taxon_id, ref_family = family)
    ) %>%
    mutate(
      direction = if_else(!! diffabun_var > 1, "+", "-"),
      sign = wilcox_p < 0.10
    )
  
  taxa_ordered <-
    taxon_pairs %>%
    group_by(family) %>%
    summarize(median_diffabun = median(!! diffabun_var)) %>%
    arrange(median_diffabun) %>%
    pull(family)
  
  taxon_pairs %>%
    mutate_at(c("family", "ref_family"), factor, levels = taxa_ordered) %>%
    ggplot(aes(x = ref_family, y = family, fill = !! diffabun_var)) +
    geom_tile() +
    geom_text(
      aes(label = if_else(sign, direction, ""), col = direction), size = 2
    ) +
    scale_color_manual(values = c("+" = "black", "-" = "white"), guide = F) +
    scale_fill_continuous(trans = "log10") +
    xlab("reference family") +
    theme_minimal() + 
    theme(
      text = element_text(size = 8),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = element_blank()
    ) 
  
}
