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

#' Get the eigenvectors that explain most of the variation
#'
#' @description Gets the eigenvectors, or loadings, with which
#' the principal components can be obtained
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

#' Get the variants that are sufficient for identifying batch effects
#'
#' @description Function that gets the variants that are sufficient for
#' batch effects. Variants are selected based on the locus of interest.
#' Variants may not be within this locus. Additionally, SNPs should have
#' a HWE p < 0.01, a minor allele frequency > 5%, they should be equally
#' distributed across all chromosomes and approximately in linkage
#' equilibrium (r2 < 0.4).
#'
#' @param genotype_data Genotype data as a matrix.
#' @param locus_of_interest Locus to exlude variants from
#'
#' @return Loadings that correct batch effects.
#'
get_variants_for_batch_correction <- function(genotype_data, locus_of_interest) {
  # The locus of interest is a GRanges object. In this object,
  # multiple ranges can be included. Therefore, we loop over
  # each of these ranges.
  for (index in seq(elementNROWS(locus_of_interest)))
  locusOfInterest$ranges
}

#' Correct intensities for variants in the locus of interest.
#'
#' @description Corrects the intensities for variants within
#' the locus of interest.
#'
#' @param genotypeData Genotype data as a matrix.
#' @param locusOfInterest The locus of interest.
#' @param config
#'
#' @return Corrected batch effects.
#'
correctIntensities <- function(genotypeData, locusOfInterest, config) {
  variants <- getVariantsForCorrections(gentoypeData, locusOfInterest)
  # Perform eigenvector Decomposition performed on the sample intensity correlation matrix
  eigenvaluesOfTopPrincipalComponents <- correctionEigenvalues(genotypeData, variants)
  # Correct per sample intensities
  return(correctedGenotypeData)
}


#' Identify principal components for locus
#'
#' @description Function that for a specific locus,
#' returns the data projected on principal components.
#'
#' @param correctedGenotypeData genotype data in a locus without batch effects.
#'
#' @return Genotype data projected on principal components.
#'
pcaOverVariantsInLocus <- function(correctedGenotypeData) {
  projectedOnPrincipalComponents <- prcomp(
    correctedGenotypeData, center = TRUE, scale. = TRUE)
  return(projectedOnPrincipalComponents)
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

read_feature_file <- function(haplotype_bed1_file, feature_name) {

  # Load a bed file (1-indexed) that contains the columns,
  # chrom, chromStart, chromEnd and name
  bed_table <- fread(
    haplotype_bed1_file,
    col.names = c("chrom", "chromStart", "chromEnd", "name"))

  # Select the rows that match the requested feature
  bed_table <- bed_table[bed_table$name == feature_name,]

  # Build a GRanges object from the bed file
  locus_of_interest <- GRanges(
    seqnames = bed_table$chrom,
    ranges = IRanges(start = bed_table$chromStart, end = bed_table$chromEnd),
    featureName = bed_table$name)

  return(locus_of_interest)
}

# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv=NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  # Parse arguments
  args <- parser$parse(argv)

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
  if (!is.null(config$geneId)) {
    locusOfInterest <- extractLocus(1565)
  } else {
    locusOfInterest <- read_feature_file(config$bedFilePath, config$featureName)
  }

  # We want to identify all variants that are in the locus of interest.
  variantsInLocus <- extractVariants(locusOfInterest)

  # Correct intensities in the locus of interest using a selection of
  # samples outside of the locus of interest
  correctedGenotypeData <- correctIntensities(
    genotypeData = NULL,
    locusOfInterest = locusOfInterest,
    config = config$variantSelectionCriteria)

  # We expect that the biggest variance is represented by copy number
  # variations. We want te re-identify these copy number variations.
  # To do this, we perform a PCA, projecting data on principal components,
  # each of these represent orthoganol bits of the variation.
  projectedGenotypeData <- pcaOverVariantsInLocus(correctedGenotypeData)


}

if (sys.nframe() == 0 && !interactive()) {
  main()
}