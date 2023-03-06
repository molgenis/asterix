#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)
library(class)
library(ggrepel)

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

  probabilities_processed <- tibble(fread(
    "out.cnv_probabilities.fitted.csv.gz",
    data.table=F, header=T)) %>%
    rename("Sample_ID" = "Sample ID") %>%
    mutate(Cnv = A + B, Probability_Log=log(Probability + .Machine$double.eps))

  cnv_per_variant_probabilities <- probabilities_processed %>%
    group_by(Sample_ID, Variant, Cnv) %>%
    summarise(Cnv_Probability_Log = sum(Probability))

  cnv_per_variant_dosages <- cnv_per_variant_probabilities %>%
    group_by(Sample_ID, Variant) %>%
    summarise(Cnv_Dosage = sum(Cnv * Cnv_Probability_Log)) %>%
    ungroup() %>%
    inner_join(variants_curated, by = c("Variant" = "Name")) %>%
    arrange(Start) %>%
    mutate(Variant = ordered(Variant, levels = unique(Variant)))

  library(heatmap3)

  cnv_matrix <- cnv_per_variant_dosages %>% pivot_wider(id_cols = Sample_ID, names_from = "Variant", values_from = "Cnv_Dosage") %>% select(-Sample_ID) %>% as.matrix


  png("out.feb_test.dosages.clustered.png")

  heatmap(cnv_matrix, Colv = NA, scale = "none")

  dev.off()

  ggplot(cnv_per_variant_dosages, aes(Variant, Cnv_Dosage, group = Sample_ID)) +
    geom_point(size=0.5) +
    geom_line(size=0.5)

  ggsave(sprintf("out.feb_test.dosages.png", variant_name))

  cnv_probabilities <- cnv_per_variant_probabilities%>%
    group_by(Sample_ID, Cnv) %>%
    summarise(Cnv_Probability_Total = exp(sum(log(Cnv_Probability_Log + .Machine$double.eps), na.rm=T))) %>%
    group_by(Sample_ID) %>%
    mutate(
      Cnv_Probability_Adjusted = Cnv_Probability_Total / sum(Cnv_Probability_Total, na.rm=T))

  probabilities_adjusted <- cnv_probabilities %>%
    ungroup() %>%
    inner_join(probabilities_processed, by = c("Sample_ID", "Cnv")) %>%
    group_by(Sample_ID, Variant, Cnv) %>%
    mutate(
      Cnv_Variant_Probability = sum(Probability, na.rm=T),
      Adjusted_Probability = Probability / Cnv_Variant_Probability * Cnv_Probability_Adjusted) %>%
    group_by(Sample_ID, Variant) %>%
    slice_max(Adjusted_Probability) %>%
    mutate(Probability = Adjusted_Probability)
  #
  # probabilities_fitted <- tibble(fread(
  #   "out.cnv_probabilities.fitted.csv.gz",
  #   data.table=F, header=T)) %>% rename("Sample_ID" = "Sample ID") %>%
  #   group_by(Variant, Sample_ID) %>%
  #   slice_max(Probability)

  predicted_genotypes <- probabilities_adjusted %>%
    mutate(Genotype = paste0(
      case_when(A == 0 & B == 0 ~ "0", TRUE ~ ""),
      case_when(A == 0 ~ "", A == 1 ~ "A", A == 2 ~ "AA", A > 2 ~ paste0(as.character(A), "A")),
      case_when(B == 0 ~ "", B == 1 ~ "B", B == 2 ~ "BB", B > 2 ~ paste0(as.character(B), "B"))),
          CNV = A + B,
          Genotype_Theta = log2(A / B)) %>%
    rename("A_Dosage" = "A", "B_Dosage" = "B")

  for (variant_name in variants_curated$Name) {
    message(variant_name)

    intensities <- tibble(fread(
      sprintf("out.ab_intensity_dataset.%s.csv.gz", variant_name), data.table=F))

    intensities_pivotted <- pivot_longer(intensities, -c(`Sample ID`, subspace)) %>%
      pivot_wider(names_from = subspace, values_from = value) %>%
      rename("Sample_ID" = "Sample ID") %>%
      inner_join(predicted_genotypes %>% filter(Variant == variant_name), by = c("Sample_ID")) %>%
      group_by(Genotype) %>%
      mutate(mean_diff_euclidian = sqrt((A - mean(A, na.rm=T))^2 + (B - mean(B, na.rm=T))^2)) %>%
      mutate(Label=if_else(min(mean_diff_euclidian) == mean_diff_euclidian, Genotype, ""))

    ggplot(intensities_pivotted,
           aes(A, B, color=Genotype, label=Label, shape = Probability > 0.95)) +
      geom_point(size=0.2, shape=16) +
      geom_text_repel(box.padding = 1,
                      show.legend = FALSE) +
      coord_fixed() +
      theme(aspect.ratio = 1) +
      facet_wrap(~ name)

    ggsave(sprintf("out.feb_test.nn.%s.png", variant_name))

    # medians <- intensities_nn %>% filter(cnv == "-2") %>%
    #   summarise(A = median(A), B = median(B))
    #
    # theta_nn <- intensities_nn %>%
    #   mutate(A2 = A - medians$A, B2 = B - medians$B) %>%
    #   mutate(theta = atan(A2/B2) * 180 / pi,
    #          theta = if_else(theta < -45, 180 + theta, theta))
    # # pivot_wider(id_cols = c("Sample ID", "cnv"),
    # #             names_from="name",
    # #             values_from = c("A", "B", "theta"))
    # #
    # # d <- dist(theta_nn %>% filter(cnv == "1") %>% select(`Sample ID`, starts_with("theta")), method = "euclidean")
    # # hc1 <- hclust(d, method = "complete")
    #
    # ggplot(theta_nn,
    #        aes(theta, color = cnv)) +
    #   geom_histogram() +
    #   facet_wrap(~ cnv + name)
    #
    # ggsave(sprintf("out.kerst_test.nn.theta.%s.png", variant_name))

  }




  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}