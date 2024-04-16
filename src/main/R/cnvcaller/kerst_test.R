#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)
library(class)

# Declare constants

# Declare function definitions
function(a, b) {

}

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Process input
  probabilities <- tibble(fread(
    "out.naive_clustering.updated_probabilities.csv.gz",
    data.table=F, header=T)) %>%
    pivot_longer(-c(`Sample ID`), names_to = "cnv", values_to = "prob") %>%
    mutate(prob = exp(prob)) %>%
    group_by(`Sample ID`) %>%
    mutate(prob = prob / sum(prob, na.rm=T),
           max = max(prob, na.rm=T)) %>%
    ungroup() %>%
    filter(max==prob) %>%
    mutate(cnv=ordered(cnv))

  assays_curated <- fread(
    "/groups/umcg-fg/tmp01/projects/pgx-passport/tools/PGx-passport-pilot/data/cyp2d6/configuration_bed_files/variants_curated_delta.bed",
    col.names = c("Chromosome", "Start", "End", "Name"))

  variants_curated <- assays_curated %>%
    group_by(Chromosome, Start, End) %>%
    arrange(Name) %>% slice_head(n=1)

  for (variant_name in variants_curated$Name) {
    message(variant_name)

    intensities <- tibble(fread(
      sprintf("out.ab_intensity_dataset.%s.csv.gz", variant_name), data.table=F))

    intensities_pivotted <- pivot_longer(intensities, -c(`Sample ID`, subspace)) %>%
      pivot_wider(names_from = subspace, values_from = value) %>%
      inner_join(probabilities, by = c("Sample ID"))

    intensities_df <- intensities_pivotted %>%
      complete(name, `Sample ID`) %>%
      group_by(`Sample ID`) %>%
      filter(all(!is.na(A)), all(!is.na(B)), prob > 0.90) %>%
      ungroup() %>%
      pivot_wider(id_cols = c("Sample ID", "cnv"),
                  names_from="name",
                  values_from = c("A", "B"))

    nearest_neighbours <- knn(train = intensities_df %>% select(where(is_double)),
        test = intensities_df %>% select(where(is_double)),
        k=5, l=4, cl=intensities_df %>% pull(cnv))

    intensities_nn <- intensities_df %>% mutate(nn = nearest_neighbours) %>%
      select(`Sample ID`, nn) %>%
      inner_join(intensities_pivotted, by = c("Sample ID")) %>%
      mutate(nn = ordered(nn)) %>%
      filter(nn == cnv)

    ggplot(intensities_nn,
           aes(A, B, color = cnv)) +
      geom_point() +
      facet_wrap(~ cnv + name)

    ggsave(sprintf("out.kerst_test.nn.%s.png", variant_name))

    medians <- intensities_nn %>% filter(cnv == "-2") %>%
      summarise(A = median(A), B = median(B))

    theta_nn <- intensities_nn %>%
      mutate(A2 = A - medians$A, B2 = B - medians$B) %>%
      mutate(theta = atan(A2/B2) * 180 / pi,
             theta = if_else(theta < -45, 180 + theta, theta))
    # pivot_wider(id_cols = c("Sample ID", "cnv"),
    #             names_from="name",
    #             values_from = c("A", "B", "theta"))
    #
    # d <- dist(theta_nn %>% filter(cnv == "1") %>% select(`Sample ID`, starts_with("theta")), method = "euclidean")
    # hc1 <- hclust(d, method = "complete")

    ggplot(theta_nn,
           aes(theta, color = cnv)) +
      geom_histogram() +
      facet_wrap(~ cnv + name)

    ggsave(sprintf("out.kerst_test.nn.theta.%s.png", variant_name))

  }


  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}