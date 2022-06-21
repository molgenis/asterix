#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(argparse)

# Declare constants
old <- theme_set(theme_minimal(base_size = 12))
theme_update(
  line = element_line(
    colour = "black", size = (1 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  panel.grid.minor = element_blank(),
  text = element_text(family="Lato"),
  title = element_text(colour = "#595A5C", face = "bold")
)

color_vector <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000")

options(ggplot.discrete.colour = color_vector,
        ggplot.discrete.fill = color_vector)

parser <- ArgumentParser(description = 'Visualize cnv calling output.')
parser$add_argument('--path_ref', metavar = 'path', type = 'character',
                    help = 'Path to a directory with output files, used as a reference to compare plot a new dataset to')
parser$add_argument('--path_new', metavar = 'path', type = 'character', required = TRUE,
                    help = 'Path to a directory with output files.')

# Declare function definitions

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

load_cnv_data <- function(cnv_path_prefix) {

  # Load different data files for a cnv calling dataset:

  # The principal component batch effects,
  batch_effects <- load_batch_effects(cnv_path_prefix)

  # The corrected intensities and projected intensities for the variants of interest
  intensities_corrected <- read_intensities_processed(
    paste0(cnv_path_prefix, "intensity_correction", "corrected", sep = "."), values_to = "R_corrected")
  intensities_projected <- read_intensities_processed(
    paste0(cnv_path_prefix, "copy_number_assignment", "projected_intensities", sep = "."), values_to = "R_corrected")

  # The corrected intensities dataframe
  corrected_intensities_data_frame <- read_intensities_raw(
    cnv_path_prefix) %>%
    inner_join(intensities_corrected, by = c("Sample_ID", "variant")) %>%

  # The processed intensities for the variants of interest
  processed_intensities_data_frame <- read_intensities(

  )

  # Naive clustering (centroids + assignment)
  # GMM estimatins, means variance
  # GMM probabilities and or dosages per feature and overall.
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
  data_new <- load_cnv_data(args$path_new)

  # Perform method
  # Generate plots
  # PCA results (no large outliers)
  ggplot(data = data_new$batch_effects, aes(component_x_value, component_y_value, color = dataset)) +
    geom_point(data)

  # intensities per variable colored by naive clustering (do samples overlay nicely with exisiting clusters)


  # intensities per processed variables colored by naive clustering
  # intensities per processed variable, with naive initial and estimated mean, variance, etc...
  # predicted overall dosage against processed intensities for each variable,
  # (failed samples highlighted)
  # predicted dosage per variable against processed intensities
  # (failed samples highlighted)

  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}