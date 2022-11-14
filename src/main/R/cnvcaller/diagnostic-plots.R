#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(argparse)
library(yaml)
library(ggrastr)
library(ggh4x)
library(extrafont)

loadNamespace("MASS")

# Declare constants
loadfonts()
old <- theme_set(theme_classic())
theme_update(
  line = element_line(
    colour = "black", size = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  panel.grid.minor = element_blank(),
  text = element_text(family="Lato"),
  title = element_text(colour = "#595A5C", face = "bold")
)

caw_color_vector <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000")

parser <- ArgumentParser(description = 'Visualize cnv calling output.')
parser$add_argument('--path_ref', metavar = 'path', type = 'character',
                    help = 'Path to a directory with output files, used as a reference to compare plot a new dataset to')
parser$add_argument('--path_new', metavar = 'path', type = 'character', required = TRUE,
                    help = 'Path to a directory with output files.')

# Declare function definitions

caw_pal <- function(
  primary = "#D55E00",
  other = "#009E73",
  direction = 1
) {
  stopifnot(primary %in% caw_color_vector)

  function(n) {
    if (n > 8) warning("CAW Color Palette only has 8 colors.")

    if (n == 2) {
      color_list <- c(other, primary)
    } else {
      color_list <- caw_color_vector[1:n]
    }

    if (direction >= 0) color_list else rev(color_list)
  }
}

scale_colour_caw <- function(
  primary = "#D55E00",
  other = "#009E73",
  direction = 1,
  ...
) {
  ggplot2::discrete_scale(
    "colour", "caw",
    caw_pal(primary, other, direction),
    ...
  )
}

scale_color_caw <- scale_colour_caw


plot_intensities <- function() {

  ggplot(
    projected_updated_data_frame %>% filter(dataset == "fit"),
    aes(factor(Start), R_intensity)) +
    rasterize(geom_jitter(
      alpha = 0.18, shape = 16, size = 0.75,
      aes(colour = R_intensity)), dpi = 300) +
    # scale_colour_manual(values = brewer.pal(n = 7, name = 'Set2')) +
    scale_colour_viridis_c(
      option = "plasma", name = "Intensity", na.value = "grey50") +
    rasterize(geom_jitter(
      data = projected_updated_data_frame %>% filter(Sample_ID %in% mom_samples_mistakes),
      shape = 21, colour = "black",
      aes(fill = CNV_status)), dpi = 300) +
    scale_fill_viridis_c(
      limits = c(0, NA_integer_),
      option = "plasma", name = "CNV,\nstatus", na.value = "grey50") +
    scale_colour_viridis_c(
      limits = c(0, NA_integer_),
      option = "plasma", name = "CNV,\nstatus", na.value = "grey50") +
    theme(legend.position = "bottom", axis.line = element_blank())

  dev.off()
}

add_highlight_sample <- function() {
  # Highlight a sample
  geom_bin2d() + scale_fill_viridis_c(
    limits = c(0, NA_integer_),
    option = "plasma", name = "CNV,\nstatus", na.value = "grey50")
}

plot_predicted_dosages_per_feature <- function() {
  ggplot(
    projected_updated_data_frame %>% filter(dataset == "fit"),
    aes(index, dosage)) +
    rasterize(geom_jitter(
      alpha = 0.18, shape = 16, size = 0.75,
      aes(colour = R_intensity)), dpi = 300) +
    # scale_colour_manual(values = brewer.pal(n = 7, name = 'Set2')) +
    scale_colour_viridis_c(
      option = "plasma", name = "Intensity", na.value = "grey50") +
    rasterize(geom_jitter(
      data = projected_updated_data_frame %>% filter(Sample_ID %in% mom_samples_mistakes),
      shape = 21, colour = "black",
      aes(fill = CNV_status)), dpi = 300) +
    scale_fill_viridis_c(
      limits = c(0, NA_integer_),
      option = "plasma", name = "CNV,\nstatus", na.value = "grey50") +
    scale_colour_viridis_c(
      limits = c(0, NA_integer_),
      option = "plasma", name = "CNV,\nstatus", na.value = "grey50") +
    theme(legend.position = "bottom", axis.line = element_blank())

  dev.off()
}

load_batch_effects <- function(path) {
  return(read_csv(
  paste(path, "intensity_correction.batcheffects.csv.gz", sep = "."),
  col_types = cols("Sample ID" = col_character())) %>%
  rename("Sample_ID" = `Sample ID`))
}

read_intensities_raw <- function(cnv_path_prefix) {
  return(read_tsv(
    paste(cnv_path_prefix, "intensity_data_frame.csv.gz", sep = "."),
    col_types = cols("Sample ID" = col_character())) %>%
    rename_all(list(~stringr::str_replace_all(., '\\s', '_'))))
}

read_intensities_processed <- function(cnv_path_prefix, values_to) {
  return(read_csv(
    paste(cnv_path_prefix, "csv.gz", sep = "."),
    col_types = cols(`Sample ID` = col_character())) %>%
    rename("Sample_ID" = `Sample ID`) %>%
    pivot_longer(cols = -Sample_ID, values_to = values_to, names_to = "variant"))
}

read_mixture_probabilities <- function(cnv_path_prefix) {
  return(read_csv(
    paste(cnv_path_prefix, "copy_number_assignment.probabilities.csv.gz", sep = "."),
    col_types = cols("Sample ID" = col_character())) %>%
    rename("Sample_ID" = `Sample ID`) %>%
    rowwise() %>%
    mutate(Dosage = rowSums(across(-Sample_ID, ~ as.numeric(cur_column()) * .x)),
           Probability_Max = max(c_across(c(-Sample_ID, -Dosage))),
           Probability_Low = Probability_Max < 0.97,
           Hard_Dosage = round(Dosage)) %>% ungroup()
  )
}

read_naive_clustering <- function(cnv_path_prefix) {
  return(read_csv(
    paste(cnv_path_prefix, "naive_clustering.assignments.csv.gz", sep = "."),
    col_types = cols("Sample ID" = col_character())) %>%
    rename("Sample_ID" = `Sample ID`))
}

load_yaml_config <- function(config_path) {
  return(yaml.load_file(config_path))
}

plot_to_pdf <- function(plot, path, width = 6, height = 6) {
  pdf(path,
      width = width, height = height, useDingbats = F)
  par(xpd = NA)

  print(plot)

  dev.off()
  embed_fonts(path, outfile=path)
}

read_mixture_fit <- function(cnv_path_prefix, names_to = "variant") {
  bind_rows(
    centroids_initial = read_csv(
      paste(cnv_path_prefix, "copy_number_assignment.centroids_initial.csv.gz", sep = ".")) %>%
      pivot_longer(cols = -marker_genotype, names_to = names_to, values_to = "centroids"),
    centroids_fit = read_csv(
      paste(cnv_path_prefix, "copy_number_assignment.centroids_fit.csv.gz", sep = ".")) %>%
      pivot_longer(cols = -marker_genotype, names_to = names_to, values_to = "centroids"), .id = "step")
}

load_cnv_data <- function(cnv_path_prefix, dimensionality_reduction_bed) {

  # Load different data files for a cnv calling dataset:

  # The principal component batch effects,
  batch_effects <- load_batch_effects(cnv_path_prefix)

  # Read mixture probabilities
  mixture_probabilities <- read_mixture_probabilities(cnv_path_prefix)

  # Read mixture centriods
  mixture_fit <- dimensionality_reduction_bed %>%
    rowid_to_column("feature_identifier") %>%
    mutate(feature_identifier = as.character(feature_identifier)) %>%
    inner_join(mixture_fit, by = c("feature_identifier")) %>%
    mutate(variant = paste0(Chrom, ":", Start)) %>%
    mutate(variant = ordered(variant, levels = unique(variant[order(Start)])))

  mixture_fit <- read_mixture_fit(cnv_path_prefix, names_to="variant") %>%
    inner_join(variants, by = c("variant" = "Name")) %>%
    mutate(variant = ordered(variant, levels = levels(variants$Name)))

  # Read naive clustering
  #naive_clustering <- read_naive_clustering(cnv_path_prefix)

  # The corrected intensities and projected intensities for the variants of interest
  intensities_corrected <- read_intensities_processed(
    paste(cnv_path_prefix, "intensity_correction", "corrected", sep = "."), values_to = "R_corrected")
  intensities_projected <- read_intensities_processed(
    paste(cnv_path_prefix, "copy_number_assignment", "projected_intensities", sep = "."), values_to = "R_projected")

  # The corrected intensities dataframe
  corrected_intensities_data_frame <- read_intensities_raw(
    cnv_path_prefix) %>%
    inner_join(intensities_corrected, by = c("Sample_ID", "variant")) %>%
    inner_join(mixture_probabilities, by = "Sample_ID") %>%
    #inner_join(naive_clustering, by = "Sample_ID") %>%
    mutate(theta = atan(Y/X),
           theta_ = atan(sqrt(Y)/sqrt(X)),
           X_corrected = X * (R_corrected / R),
           Y_corrected = Y * (R_corrected / R))

  projected_intensities_data_frame <- dimensionality_reduction_bed %>%
    rowid_to_column("Feature_ID") %>%
    mutate(Feature_ID = as.character(Feature_ID)) %>%
    inner_join(intensities_projected, by = c("Feature_ID" = "variant")) %>%
    inner_join(mixture_probabilities, by = "Sample_ID") %>%
    inner_join(naive_clustering, by = "Sample_ID") %>%
    mutate(variant = paste0(Chrom, ":", Start)) %>%
    mutate(variant = ordered(variant, levels = unique(variant[order(Start)])))

  # The processed intensities for the variants of interest

  # Naive clustering (centroids + assignment)
  # GMM estimatins, means variance
  # GMM probabilities and or dosages per feature and overall.
}

intensity_warped <- function(x, y, labels, method = "ellipsoid") {
  # Get boundaries that separates two clusters
  print(str(labels))
  print(str(x))
  print(str(y))
  r <- x+y
  if (method == "lda") {
    linear_discriminant_model <- MASS::lda(labels ~ x + y, na.action = na.omit)
    linear_discriminant_analysis <- predict(linear_discriminant_model)
    linear_discriminant_components <- linear_discriminant_analysis$x
    print(str(linear_discriminant_components))
    return(linear_discriminant_components)
  } else if (method == "ellipsoid") {
    x_max_centroid <- mean(x[x > quantile(x, prob=0.99)])
    y_max_centroid <- mean(y[y > quantile(y, prob=0.99)])
    x_min_centroid <- mean(x[y > quantile(y, prob=0.99)])
    y_min_centroid <- mean(y[x > quantile(x, prob=0.99)])
    return(list("x_max_centroid" = x_max_centroid,
           "y_max_centroid" = y_max_centroid,
           "x_min_centroid" = x_min_centroid,
           "y_min_centroid" = y_min_centroid))
  }
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
  args <- parser$parse_args(argv)

  # Load variants, settings etc.

  # We can plot two datasets, one with the original fit, and one with a new dataset as an overlay.
  if (!is.null(args$path_ref)) {
    data_reference <- load_cnv_data(args$path_ref)
  }

  config <- load_yaml_config(args$config)
  config <- load_yaml_config(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/tools/asterix/src/main/python/cnvcaller/conf/config.yml")
  dimensionality_reduction_bed <- read_tsv(
    config[["dimensionality reduction"]],
    col_names=c("Chrom", "Start", "End"))
  variants_curated_bed <- read_tsv(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/tools/PGx-passport-pilot/data/cyp2d6/configuration_bed_files/variants_curated_beta.bed",
    col_names=c("Chrom", "Start", "End", "Name"))

  cyp2d6loc <- read_tsv("/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/tools/PGx-passport-pilot/data/cyp2d6/cyp2d6.bed",
                        col_names = c("Chrom", "Start", "End", "Name"))

  variants <- read_tsv(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/results/ugli/analyses/select_corrective_variants/out.locus.bed",
    col_names = c("Chrom", "Start", "End", "Name")) %>%
    mutate(Name = ordered(Name, levels = .$Name[order(.$Start)])) %>%
    rowwise() %>%
    mutate(StartDistance = min(0, Start - cyp2d6loc[1, "Start", drop = T]),
           EndDistance = max(0, End - cyp2d6loc[1, "End", drop = T]),
           Distance = c(StartDistance, EndDistance)[which.max(abs(c(StartDistance, EndDistance)))]) %>%
    ungroup()

  cnv_path_prefix <- "out"

  simple_plot <- ggplot(data = corrected_intensities_data_frame %>% filter(variant %in% c("rs1135840", "DUP-rs1135840")),
         aes(X_corrected, Y_corrected)) +
    rasterize(geom_point(
      shape = 16, alpha = 0.16, size = 0.2,
      aes(colour = R_corrected)), dpi = 300) +
    scale_colour_viridis_c(option = "viridis",
                           name = "R", na.value = "grey50") +
    theme(legend.position = "bottom",
          aspect.ratio = 1, axis.text = element_blank()) +
    facet_wrap(~ variant, ncol = 2, scales = "free")

  plot_to_pdf(
    simple_plot, file.path(".", sprintf(
      "out.simple03.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 4, height = 3)

  simple_plot_r <- ggplot(data = corrected_intensities_data_frame %>% filter(variant %in% c("rs1135840", "DUP-rs1135840")),
                        aes(x=theta, y = R_corrected)) +
    rasterize(geom_jitter(
      shape = 16, alpha = 0.16, size = 0.2,
      aes(colour = R_corrected)), dpi = 300) +
    scale_colour_viridis_c(option = "viridis",
                           name = "R", na.value = "grey50") +
    theme(legend.position = "bottom",
          aspect.ratio = 1, axis.text = element_blank()) +
    facet_wrap(~ variant, ncol = 2, scales = "free")

  plot_to_pdf(
    simple_plot_r, file.path(".", sprintf(
      "out.simple02.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 4, height = 3)

  featurization_plot_2 <- ggplot(data = featurization_data_frame,
                               aes(X_corrected, Y_corrected)) +
    rasterize(geom_point(
      shape = 16, size = 0.2,
      aes(colour = factor(component), alpha = fit)), dpi = 300) +
    scale_colour_caw(name="component") +
    # geom_bin_2d(width=0.8) +
    # scale_colour_viridis_c(option = "viridis",
    #                        name = "RESP", na.value = "grey50") +
    geom_point(data = featurization_summarized, shape = "+", color = "black") +
    theme(legend.position = "null",
          aspect.ratio = 1, axis.text = element_blank()) +
    facet_wrap(~ variant, ncol = 2, scales = "free")

  plot_to_pdf(
    featurization_plot_2, file.path(".", sprintf(
      "out.featurization_2.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 4, height = 3)

  featurization_plot_3 <- ggplot(data = featurization_data_frame,
                                 aes(X_corrected, Y_corrected)) +
    rasterize(geom_point(
      shape = 16, size = 0.2,
      aes(colour = factor(to), alpha = fit)), dpi = 300) +
    scale_colour_caw(name="component") +
    # geom_bin_2d(width=0.8) +
    # scale_colour_viridis_c(option = "viridis",
    #                        name = "RESP", na.value = "grey50") +
    geom_point(data = featurization_summarized, shape = "+", color = "black") +
    theme(legend.position = "null",
          aspect.ratio = 1, axis.text = element_blank()) +
    facet_wrap(~ variant, ncol = 2, scales = "free")

  plot_to_pdf(
    featurization_plot_3, file.path(".", sprintf(
      "out.featurization_3.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 4, height = 3)

  res <- read_csv(
    paste(cnv_path_prefix, "cnv_probabilities.mapped.csv.gz", sep = ".")) %>%
    rename("Sample_ID" = "variant", "variant" = "Sample_ID") %>%
    pivot_longer(cols = -c(Sample_ID, variant), values_to = "probability",
                 names_to = "prediction")

  res_comb <- res %>% pivot_wider(
    names_from = variant,
    values_from = probability) %>%
    mutate(overall = `chr22-42525044:rs1135824` * rs1135840) %>%
    group_by(Sample_ID) %>%
    mutate(ovarall_norm = overall / sum(overall))

  comporison_plot <- res %>% pivot_wider(
    names_from = variant,
    values_from = probability) %>%
    ggplot(aes(`chr22-42525044:rs1135824`, rs1135840)) +
    rasterize(geom_point(
      shape = 16, size = 0.2, alpha = 0.36,
      aes(colour = prediction)), dpi = 300) +
    scale_colour_caw(name="component") +
    theme(legend.position = "null",
          aspect.ratio = 1) +
    facet_wrap(~ prediction, ncol = 2, scales = "free")

  plot_to_pdf(
    comporison_plot, file.path(".", sprintf(
      "out.prediction_comparison.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 6, height = 5)

  res_max <- res %>%
    group_by(Sample_ID, variant) %>%
    filter(probability == max(probability))
  prediction_plot <- ggplot(data = res_max,
                                 aes(prediction, probability)) +
    rasterize(geom_jitter(
      shape = 16, size = 0.2, alpha = 0.36,
      aes(colour = prediction)), dpi = 300) +
    scale_colour_caw(name="component") +
    # geom_bin_2d(width=0.8) +
    # scale_colour_viridis_c(option = "viridis",
    #                        name = "RESP", na.value = "grey50") +
    theme(legend.position = "null",
          aspect.ratio = 1) +
    facet_wrap(~ variant, ncol = 2, scales = "free")

  plot_to_pdf(
    prediction_plot, file.path(".", sprintf(
      "out.prediction.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 4, height = 3)


  featurization <- read_csv(
    paste(cnv_path_prefix, "cnv_probabilities.init.csv.gz", sep = ".")) %>%
    group_by(variant, to) %>% # Group by variant, and the cnv status
    mutate(inner_component = as.numeric(factor(component)))

  featurization_data_frame <- corrected_intensities_data_frame %>%
    inner_join(featurization, by=c("Sample_ID"="Sample_ID", "variant"="variant"))

  featurization_summarized <- featurization_data_frame %>%
    group_by(variant, component) %>%
    filter(!is.na(X_corrected), !is.na(Y_corrected), !is.na(fit)) %>%
    summarise(X_corrected = sum(X_corrected * fit)/sum(fit),
              Y_corrected = sum(Y_corrected * fit)/sum(fit))

  featurization_plot <- ggplot(data = featurization_data_frame %>% filter(variant == "chr22-42525044:rs1135824"),
         aes(X_corrected, Y_corrected)) +
    rasterize(geom_point(
      shape = 16, size = 0.2,
      aes(colour = factor(to), alpha = fit)), dpi = 300) +
    scale_colour_caw(name="component") +
    # geom_bin_2d(width=0.8) +
    # scale_colour_viridis_c(option = "viridis",
    #                        name = "RESP", na.value = "grey50") +
    geom_point(data = featurization_summarized, shape = "+", color = "black") +
    theme(legend.position = "bottom",
          aspect.ratio = 1, axis.text = element_blank()) +
    facet_wrap(~ variant + to + inner_component, ncol = 4, scales = "free")

  plot_to_pdf(
    featurization_plot, file.path(".", sprintf(
      "out.featurization.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 8, height = 11)

  data_new <- load_cnv_data(
    args$path_new,
    dimensionality_reduction_bed)

  # Perform method
  # Generate plots
  # PCA results (no large outliers)
  ggplot(data = data_new$batch_effects, aes(component_x_value, component_y_value, color = dataset)) +
    geom_point(data)

  corrected_R_data_frame <- corrected_intensities_data_frame %>%
    pivot_wider(
      id_cols = c(Sample_ID, marker_genotype),
      names_from = variant,
      values_from = R_corrected
    )

  corrected_R_matrix <- corrected_R_data_frame %>%
    dplyr::select(-Sample_ID) %>%
    select(marker)

  linear_discriminant_model <- MASS::lda(formula = as.formula("marker_genotype~."), data = corrected_R_matrix)
  linear_discriminant_analysis <- predict(linear_discriminant_model, corrected_R_matrix)
  ovarall_projection <- as_tibble(linear_discriminant_analysis$x) %>%
    mutate(
      Sample_ID = corrected_R_data_frame %>% dplyr::pull(Sample_ID),
      lda_class = as.numeric(levels(linear_discriminant_analysis$class))[linear_discriminant_analysis$class]) %>%
    inner_join(mixture_probabilities, by = "Sample_ID") %>%
    inner_join(naive_clustering, by = "Sample_ID")

  qda_outcome <- as_tibble(linear_discriminant_analysis$posterior) %>%
    mutate(
      Sample_ID = corrected_R_data_frame %>% dplyr::pull(Sample_ID),
      qda_class_prediction = linear_discriminant_analysis$class) %>%
    pivot_longer(cols = where(is.numeric), values_to = "qda_class_probability", names_to = "qda_class_label") %>%
    group_by(Sample_ID) %>%
    filter(qda_class_prediction == qda_class_label) %>%
    ungroup()

  # intensities per variable colored by naive clustering (do samples overlay nicely with exisiting clusters)

  overall_projection_plot <- ggplot(data = ovarall_projection, aes(LD1, LD2)) +
    rasterize(geom_point(
      alpha = 0.36, shape = 16, size = 0.2,
      aes(colour = lda_class)), dpi = 300) +
    # geom_bin_2d(width=0.8) +
    scale_colour_viridis_c(option = "viridis",
                           name = "CNV status", na.value = "grey50") +
    theme(legend.position = "bottom",
          aspect.ratio = 1)

  plot_to_pdf(
    overall_projection_plot, file.path(".", sprintf(
    "out.overall_projection.%s.pdf", format(Sys.Date(), "%Y%m%d"))))

  intensities_per_type <- corrected_intensities_data_frame %>%
    pivot_longer(
      cols = c(X, Y, X_corrected, Y_corrected, R, R_corrected),
      names_pattern = "(.)(_corrected)?",
      names_to = c(".value", "corrected")) %>%
    dplyr::select(variant, Sample_ID, corrected, X, Y, Dosage, Hard_Dosage, Probability_Max, marker_genotype, R, theta) %>%
    mutate(corrected = corrected == "_corrected") %>%
    inner_join(variants, by = c("variant" = "Name")) %>%
    mutate(variant = ordered(variant, levels = levels(variants$Name)))

  scales_x <- intensities_per_type %>%
    filter(Distance == 0 & corrected) %>%
    group_by(variant) %>%
    group_map(~ scale_x_continuous(limits = c(
      min(unlist(.x[c("X", "Y")]), na.rm = TRUE),
      max(unlist(.x[c("X", "Y")]), na.rm = TRUE))))

  scales_y <- intensities_per_type %>%
    filter(Distance == 0 & corrected) %>%
    group_by(variant) %>%
    group_map(~ scale_y_continuous(limits = c(
      min(unlist(.x[c("X", "Y")]), na.rm = TRUE),
      max(unlist(.x[c("X", "Y")]), na.rm = TRUE))))

  # intensities per variable colored by naive clustering (do samples overlay nicely with exisiting clusters)
  corrected_intensities_dosage_plot <- ggplot(data = intensities_per_type %>%
    filter(Distance == 0 & corrected) %>%
    arrange(1-Probability_Max), aes(X, Y)) +
    rasterize(geom_point(
      alpha = 0.36, shape = 16, size = 0.2,
      aes(colour = Dosage)), dpi = 300) +
    # geom_bin_2d(width=0.8) +
    scale_colour_viridis_c(option = "viridis",
                           name = "CNV status", na.value = "grey50") +
    theme(legend.position = "bottom",
          aspect.ratio = 1) +
    facet_wrap(~ variant, ncol = 6, scales = "free") +
    facetted_pos_scales(
      y = scales_y, x = scales_x)


  # intensities per variable colored by naive clustering (do samples overlay nicely with exisiting clusters)
  corrected_intensities_dosage_plot_ <- ggplot(data = intensities_per_type %>%
    filter(Distance == 0 & corrected & variant == "DUP-rs1135840") %>%
    arrange(1-Probability_Max), aes(X, Y)) +
    rasterize(geom_point(
      alpha = 0.36, shape = 16, size = 0.5,
      aes(colour = marker_genotype)), dpi = 300) +
    # geom_bin_2d(width=0.8) +
    scale_colour_viridis_c(option = "viridis",
                           name = "CNV status", na.value = "grey50") +
    theme(legend.position = "bottom",
          aspect.ratio = 1)

  plot_to_pdf(corrected_intensities_dosage_plot_, file.path(".", sprintf(
    "out.corrected_intensities_DUP-rs1135840.%s.pdf", format(Sys.Date(), "%Y%m%d"))))

  plot_to_pdf(corrected_intensities_dosage_plot, file.path(".", sprintf(
    "out.corrected_intensities.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
              width = 8, height = 11)

  plot_to_pdf(corrected_intensities_polar_dosage_plot, file.path(".", sprintf(
    "out.corrected_intensities_polar.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
              width = 8, height = 11)

  projected_intensities_per_type <- projected_intensities_data_frame %>%
    pivot_longer(
      cols = c(R_projected),
      names_pattern = "(.)(_projected)?",
      names_to = c(".value", "projected")) %>%
    inner_join(ovarall_projection %>% dplyr::select(Sample_ID, lda_class))

  scales_r <- projected_intensities_per_type %>%
    group_by(variant) %>%
    group_map(~ scale_y_continuous(limits = c(
      min(unlist(.x[c("R")]), na.rm = TRUE),
      max(unlist(.x[c("R")]), na.rm = TRUE))))

  # intensities per variable colored by naive clustering (do samples overlay nicely with exisiting clusters)
  projected_intensities_naive_plot <- ggplot(data = projected_intensities_per_type %>%
    arrange(1 - Probability_Max), aes(marker_genotype, R)) +
    rasterize(geom_point(
      alpha = 0.36, shape = 16, size = 0.2,
      aes(colour = marker_genotype)), dpi = 300) +
    # geom_bin_2d(width=0.8) +
    geom_point(
      data = mixture_fit, colour = "black",
      aes(marker_genotype, centroids, fill = marker_genotype, shape = step), inherit.aes = F) +
    scale_colour_viridis_c(option = "viridis",
                           name = "CNV status", na.value = "grey50") +
    scale_fill_viridis_c(option = "viridis",
                           name = "CNV status", na.value = "grey50") +
    scale_shape_manual(values = c(21, 24)) +
    theme(legend.position = "bottom", axis.text = element_blank(),
          aspect.ratio = 1) +
    facet_wrap(~variant, ncol = 6, scales = "free") +
    facetted_pos_scales(
      y = scales_r)

  # intensities per variable colored by naive clustering (do samples overlay nicely with exisiting clusters)
  projected_intensities_dosage_plot <- ggplot(data = projected_intensities_per_type %>%
    arrange(1 - Probability_Max), aes(Dosage, R)) +
    rasterize(geom_point(
      alpha = 0.36, shape = 16, size = 0.2,
      aes(colour = Dosage)), dpi = 300) +
    # geom_bin_2d(width=0.8) +
    geom_point(
      data = mixture_fit, colour = "black",
      aes(marker_genotype, centroids, fill = marker_genotype, shape = step), inherit.aes = F) +
    scale_colour_viridis_c(option = "viridis",
                           name = "CNV status", na.value = "grey50") +
    scale_fill_viridis_c(option = "viridis",
                         name = "CNV status", na.value = "grey50") +
    scale_shape_manual(values = c(21, 24)) +
    theme(legend.position = "bottom", axis.text = element_blank(),
          aspect.ratio = 1) +
    facet_wrap(~variant, ncol = 6, scales = "free") +
    facetted_pos_scales(
      y = scales_r)


  plot_to_pdf(projected_intensities_naive_plot, file.path(".", sprintf(
    "out.projected_intensities.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 8)

  plot_to_pdf(projected_intensities_dosage_plot, file.path(".", sprintf(
    "out.projected_intensities_dosage.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 8)

  corrected_intensities_per_type <- intensities_per_type %>%
    filter(Distance == 0, corrected, !is.na(X), !is.na(Y)) %>%
    group_by(variant) %>%
    mutate(
      theta_explained_R = predict(lm(R ~ poly(theta, 2, raw=TRUE)))) %>%
    ungroup() %>%
    group_by(Start) %>%
    group_map(~ lda_for_group(.x), .keep = TRUE) %>% bind_rows() %>%
    inner_join(naive_clustering) %>%
    inner_join(mixture_probabilities)

  scales_r <- corrected_intensities_per_type %>%
    group_by(variant) %>%
    group_map(~ scale_y_continuous(limits = c(
      min(unlist(.x[c("R")]), na.rm = TRUE),
      max(unlist(.x[c("R")]), na.rm = TRUE))))

  corrected_intensities_polar_dosage_plot <- ggplot(data = corrected_intensities_per_type %>%
    arrange(1-Probability_Max), aes(marker_genotype, LD1)) +
    rasterize(geom_jitter(
      alpha = 0.36, shape = 16, size = 0.2,
      aes(colour = marker_genotype)), dpi = 300) +
    # geom_bin_2d(width=0.8) +
    # stat_smooth(method = "lm", formula = y ~ poly(x, 2, raw=TRUE), na.rm = T) +
    scale_colour_viridis_c(option = "viridis",
                           name = "CNV status", na.value = "grey50") +
    theme(legend.position = "bottom", axis.text = element_blank(),
          aspect.ratio = 1) +
    facet_wrap(~ Start, ncol = 6, scales = "free")

  plot_to_pdf(corrected_intensities_polar_dosage_plot, file.path(".", sprintf(
    "out.corrected_intensities_polar_dosage.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 8, height = 11)

  # intensities per variable colored by naive clustering (do samples overlay nicely with exisiting clusters)
  corrected_intensities_naive_plot <- ggplot(data = corrected_intensities_per_type %>%
    arrange(1 - Probability_Max), aes(marker_genotype, R)) +
    rasterize(geom_jitter(
      alpha = 0.36, shape = 16, size = 0.2,
      aes(colour = marker_genotype)), dpi = 300) +
    # geom_bin_2d(width=0.8) +
    geom_point(
      data = mixture_fit, colour = "black",
      aes(marker_genotype, centroids, fill = marker_genotype, shape = step), inherit.aes = F) +
    scale_colour_viridis_c(option = "viridis",
                           name = "CNV status", na.value = "grey50") +
    scale_fill_viridis_c(option = "viridis",
                         name = "CNV status", na.value = "grey50") +
    scale_shape_manual(values = c(21, 24)) +
    theme(legend.position = "bottom", axis.text = element_blank(),
          aspect.ratio = 1) +
    facet_wrap(~variant, ncol = 6, scales = "free") +
    facetted_pos_scales(
      y = scales_r)

  # intensities per variable colored by naive clustering (do samples overlay nicely with exisiting clusters)
  corrected_intensities_dosage_plot <- ggplot(data = corrected_intensities_per_type %>%
    arrange(1 - Probability_Max), aes(Dosage, R)) +
    rasterize(geom_point(
      alpha = 0.36, shape = 16, size = 0.2,
      aes(colour = Dosage)), dpi = 300) +
    # geom_bin_2d(width=0.8) +
    scale_colour_viridis_c(option = "viridis",
                           name = "CNV status", na.value = "grey50") +
    scale_shape_manual(values = c(21, 24)) +
    theme(legend.position = "bottom", axis.text = element_blank(),
          aspect.ratio = 1) +
    facet_wrap(~variant, ncol = 6, scales = "free") +
    facetted_pos_scales(
      y = scales_r)

  plot_to_pdf(corrected_intensities_naive_plot, file.path(".", sprintf(
    "out.corrected_intensities_naive.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 8, height = 11)

  plot_to_pdf(corrected_intensities_dosage_plot, file.path(".", sprintf(
    "out.corrected_intensities_dosage.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 8, height = 11)

  # intensities per variable colored by naive clustering (do samples overlay nicely with exisiting clusters)
  concordance_pdf_path <- file.path(".", sprintf(
    "out.concordance.%s.pdf", format(Sys.Date(), "%Y%m%d")))

  pdf(concordance_pdf_path,
      width = 6, height = 6, useDingbats = F)

  par(xpd = NA)

  naive_clustering_expected %>% inner_join(naive_clustering) %>%
    ggplot(aes(marker_genotype, clusters)) +
    rasterize(geom_jitter(
      alpha = 0.36, shape = 16, size = 0.2,
      aes(colour = marker_genotype)), dpi = 300)

  dev.off()
  embed_fonts(concordance_pdf_path, outfile=concordance_pdf_path)

  # intensities per variable colored by naive clustering (do samples overlay nicely with exisiting clusters)

  correction_comparison_plot <- ggplot(data = intensities_per_type %>% filter(Distance == 0), aes(X, Y)) +
    rasterize(geom_point(
    alpha = 0.36, shape = 16, size = 0.25), dpi = 300) +
    # geom_bin_2d(width=0.8) +
  scale_colour_viridis_c(option = "viridis",
                         name = "CNV,\nstatus", na.value = "grey50") +
  coord_equal(ratio = 1) +
  theme(legend.position = "bottom", axis.text = element_blank(),
        panel.background = element_rect(fill = "grey90")) +
  facet_wrap(~ variant + corrected, ncol = 8)

  plot_to_pdf(correction_comparison_plot, file.path(".", sprintf(
    "out.correction_comparison.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 8, height = 11)
  # uncorrected intensites, corrected intensities, corrected intensities with naive clustering

  # intensities per processed variable colored by naive clustering
  # intensities per processed variable, with naive initial and estimated mean, variance, etc...
  # predicted overall dosage against processed intensities for each variable,
  # (failed samples highlighted)
  # predicted dosage per variable against processed intensities
  # (failed samples highlighted)

  intensity_centroids <- corrected_intensities_data_frame %>%
    filter(variant == "DUP-rs1135840", !is.na(X), !is.na(Y)) %>%
    group_by(variant) %>%
    group_modify(~ as_tibble(intensity_warped(.x[["X"]], .x[["Y"]], labels = .x[["marker_genotype"]])), .keep = T) %>%
    pivot_longer(
      cols = c(ends_with("centroid")),
      names_pattern = "(.)_(min|max)_centroid",
      names_to = c(".value", "end"))

  warped_intensities_list


  corrected_intensities_polar_dosage_plot <- ggplot(data = warped_intensities, aes(LD1, LD2)) +
    rasterize(geom_point(
      alpha = 0.36, shape = 16, size = 0.2), dpi = 300) +
    # geom_bin_2d(width=0.8) +
    scale_colour_viridis_c(option = "viridis",
                           name = "CNV status", na.value = "grey50") +
    theme(legend.position = "bottom", axis.text = element_blank(),
          aspect.ratio = 1) +
    facet_wrap(~ variant, ncol = 6, scales = "free")

  plot_to_pdf(corrected_intensities_polar_dosage_plot, file.path(".", sprintf(
    "out.tmptest.%s.pdf", format(Sys.Date(), "%Y%m%d"))), width = 6, height = 6)

  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}