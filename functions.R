#' Add pairwise comparisons to a ggplot plot
#'
#' This function adds brackets and labels that represent pairwise comparisons
#' between groups of observations to a ggplot plot.
#' 
#' The comparison table should contain the following variables: 
#' 
#' * group_1
#' * group_2
#' * label
#' 
#' The groups in the group_1 and group_2 variables should correspond to groups
#' that are plotted on the x-axis in the plot.
#' 
#' The vscale parameter should be expressed as a fraction of the range of the y
#' variable
#'
#' @param plot A ggplot object
#' @param comparisons A comparison table
#' @param group_levels The group levels in x-axis order
#' @param facet_vars Variables used to subset the data (character vector)
#' @param vscale Vertical distance between brackets
#' @param size_brackets Size of the brackets
#' @param size_labels Size of the labels
#' 
#' @return A ggplot object
#' 
#' @export
add_comparisons <- function(
    plot, comparisons, group_levels, facet_vars = character(0), vscale = 0.05, 
    size_brackets = 0.5, size_labels = 4
  ) {

  facets <-
    plot$data %>%
    mutate(y = !! plot$mapping$y) %>%
    group_by_at(.vars = vars(one_of(facet_vars))) %>%
    summarize(ymax = max(y), yrange = max(y) - min(y))
  
  comparisons <-
    comparisons %>%
    mutate_at(c("group_1", "group_2"), factor, levels = group_levels) %>% 
    mutate_at(c("group_1", "group_2"), as.numeric) %>%
    mutate(
      x_middle = map2_dbl(group_1, group_2, ~ min(.x, .y) + abs(.x - .y) / 2)
    ) %>%
    arrange(group_1, group_2) %>%
    {
      if (length(facet_vars) == 0) 
        mutate(., ymax = facets$ymax, yrange = facets$yrange)
      else left_join(., facets)
    } %>%
    group_by_at(.vars = vars(one_of(facet_vars))) %>%
    mutate(y = ymax + 1:n() * (yrange * !! vscale)) %>% 
    ungroup()
  
  hsegments <- 
    comparisons %>%
    mutate(x = group_1, y = y, xend = group_2, yend = y)
  
  vsegments_left <-
    comparisons %>%
    mutate(
      x = group_1, y = y, xend = group_1, yend = y - (yrange * !! vscale) / 3
    )
  
  vsegments_right <-
    comparisons %>%
    mutate(
      x = group_2, y = y, xend = group_2, yend = y - (yrange * !! vscale) / 3
    )
  
  segments <-
    hsegments %>%
    bind_rows(vsegments_left) %>%
    bind_rows(vsegments_right)
  
  plot + 
    geom_segment(
      data = segments, aes(x = x, y = y, xend = xend, yend = yend), 
      size = size_brackets
    ) +
    geom_text(
      data = comparisons, 
      aes(x = x_middle, y = y + (yrange * !! vscale) / 3, label = label),
      size = size_labels
    )
  
}

#' Compute pairwise t-tests on a table
#'
#' This function computates pairwise t-test on a variable, with groups defined
#' by a second variable.
#'
#' @param table A tibble or data frame
#' @param value Variable containing the values (bare name)
#' @param group Variable containing the groups (bare name)
#' 
#' @return A tibble with pairwise comparisons
#' 
#' @export
pairwise_t_tests <- function(table, value, group) {
  
  value <- rlang::enexpr(value)
  group <- rlang::enexpr(group) 
  
  table <-
    table %>%
    mutate(value = !! value, group = !! group)
  
  levels <- as.character(unique(table$group))
  
  tibble(group_1 = levels, group_2 = levels) %>%
    expand(group_1, group_2) %>%
    filter(group_1 < group_2) %>%
    mutate(t_test = map2(group_1, group_2, function(grp1, grp2) {
      grp1 <- table$value[table$group == grp1]
      grp2 <- table$value[table$group == grp2]
      t.test(grp1, grp2)
    })) %>%
    mutate(
      p_value = map_dbl(t_test, "p.value"),
      t = map_dbl(t_test, "statistic")
    ) %>%
    select(- t_test)
  
}

#' Add letter to a plot
#'
#' This function adds a given letter to a plot, which is useful when
#' constructing figure panels.
#'
#' @param plot A ggplot plot
#' @param letter A single letter
#' 
#' @return A ggplot plot
#' 
#' @export
give_letter <- function(plot, letter) {
  
  g <- ggplotGrob(plot + ggtitle(letter))
  g$layout$l[g$layout$name == "title"] <- 1
  
  g
  
}