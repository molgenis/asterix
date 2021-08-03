#!/usr/bin/env Rscript

# Load libraries
library(Homo.sapiens)
library(biovizBase)
library(ggbio)
library(tidyverse)
library(argparse)
library(cowplot)

# Declare constants
## Argument parser

parser <- ArgumentParser(description='Visualize the principal components of SNP intensities in a cnv locus.')
parser$add_argument('eigenvectors', metavar='file', type = 'character',
                    help='A tab delimited file with principal components of variants in the CNV locus.')

## Plotting defautls
old <- theme_set(theme_classic())
theme_update(line = element_line(
  colour = "black", size = (0.5 / (ggplot2::.pt * 72.27/96)),
  linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  axis.line = element_line(
    colour = "#595A5C", size = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  axis.ticks = element_line(
    colour = "#595A5C", size = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T)
)

# Declare function definitions


# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv=NULL) {
  if (is.null(argv)) {
    # Process arguments
    argv <- commandArgs(trailingOnly = T)
  }
  ## Process input
  args <- parser$parse_args(argv)

  # Check if the file exists and has reading permission
  if (!(file.exists(args$eigenvectors)
        && file.access(args$eigenvectors, 4) == 0)) {
      ##stop(paste0("File '", args$eigenvectors, "', could not be read. Exiting..."))
  }

  # Read the eigenvector matrix
  #eigenvectors <- read_tsv(args$eigenvectors, col_names = F)

  # Check if the eigenvector matrix contains only doubles
  # if (!is.double(eigenvectors)) {
  #     stop("Colnames are not all doubles")
  # }

  ## Plot PCs
  #ggplot(eigenvectors, aes())

  custom_x_scale <- scale_x_continuous(limits = c(42522502, 42522502 + 4382))

  data(genesymbol, package = "biovizBase")
  geneModel <- ggplot() +
      geom_alignment(data = Homo.sapiens, which = genesymbol["CYP2D6"],
                     colour = "grey10", fill = "white", gap.geom = "chevron",
                     stat = "identity", cds.rect.h = 0.1) +
      ylab("CYP2D6 Gene model") +
      custom_x_scale

  intensities <- ggplot() +
      geom_line(aes(x = c(42522502 + 100, 42522502 + 2000, 42522502 + 4000), y = c(2, 3, 50))) +
      custom_x_scale

  plot_grid(intensities, geneModel, ncol = 1, align='v')
}

# Test if main should be run
if (sys.nframe() == 0 && !interactive()) {
  main()
}