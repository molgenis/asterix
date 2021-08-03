#!/usr/bin/env Rscript

# In the CNV Caller class we attempt to call CYP2D6/2D7/2D8 hapolotypes, as described here: by
# Twesigomwe et al., Npj Genomic Medicine (2020), https://doi.org/10.1038/s41525-020-0135-2.
# by making use of raw intensity data.

# The CNV calling principle from raw intensity data is based on a prior paper:
# Franke et al., European Journal of Human Genetics (2016), https://doi.org/10.1038/ejhg.2015.95.
# This principle is based on a couple of steps.

# - First, batch effects are corrected for by randomly selecting a number of SNPs with HWE p < 0.01, with
# a minor allele frequency > 5%, that are equally distributed across all chromosomes, and that are
# approximately in linkage disequilibrium (r2 < 0.4).
# - Using these SNPs we do eigenvector decomposition with PCA
# - Identify principal components are identified that the strongest denominating batch effects.
# - The identified principal components are used to correct the intensities across samples
# - PCA is performed on the variants of the CYP2D6 locus.

# Load libraries
library(XML)
library(illuminaio)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationFilter)
library(ggbio)
library(tidyverse)
library(argparse)
library(tools)

# Declare constants

## Argument parser
parser <- ArgumentParser(description='Visualize the principal components of SNP intensities in a cnv locus.')
parser$add_argument('--config', metavar='file', type = 'character',
                    help='A tab-delimited file with settings.')
parser$add_argument('--input-files', metavar='file', type = 'character',
                    help='A tab-delimited file that lists input idat files.')

## Plotting defaults
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

#' Get the eigenvalues that explain most of the variation
#'
#' @description Gets the eigenvalues that explain most of the variation.
#'
#' @param genotypeData Genotype data as a matrix.
#' @param variants Variants that are to be used in correction.
#'
#' @return Loadings that correct batch effects.
#'
correctionEigenvalues <- function(genotypeData, variants) {
  # Perform eigenvector decomposition
  # Select the principal components explaining the strongest batch effects
  # Return for these correction eigenvalues
}


## Correct intensities for variants in the locus of interest.
correctIntensities <- function(genotypeData, locusOfInterest, config) {
  # Perform eigenvector Decomposition performed on the sample intensity correlation matrix
  eigenvaluesOfTopPrincipalComponents <- correctionEigenvalues(genotypeData, variants)
  # Correct per sample intensities
}


## Perform PCA over variants in locus.
pcaOverVariantsInLocus <- function(correctedGenotypeData) {

}


## Variants in locus
extractLocus <- function(geneId, genomeBuild = "37") {
  # Genome build 38 or 37
  if (genomeBuild == "38") {
    TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.knownGene
  } else if (genomeBuild == "37") {
    TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg19.knownGene
  }

  # Locus of interest from the TxDb object
  locusOfInterest <- genes(Homo.sapiens,
                           filter = list(gene_id=geneId))

  # Return the locus of interest
  return(locusOfInterest)
}


### Main
main <- function(args=NULL) {
  if (is.null(args)) {
    # Process input
    args <- parser$parse_args()
  }

  # Check if the file exists and has reading permission
  if (!(file.exists(args$config)
    && file.access(args$config, 4) == 0)) {
    stop(paste0("Config file '", file_path_as_absolute(args$config),
                "', could not be read. Exiting..."))
  }

  # Check if the file exists and has reading permission
  if (!(file.exists(args$input_files)
    && file.access(args$input_files, 4) == 0)) {
    stop(paste0("Config file '",
                file_path_as_absolute(args$input_files),
                "', could not be read. Exiting..."))
  }

  # Process config file
  config <- xmlParse(file = args$config)

  # Load genotype data

  # Get the variants that are in the locus of interest
  locusOfInterest <- extractLocus(1565)

  # Correct intensities in the locus of interest using a selection of
  # samples outside of the locus of interest
  correctIntensities(genotypeData = NULL,
                     locusOfInterest = locusOfInterest,
                     config = config$variantSelectionCriteria)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}