#!/usr/bin/env Rscript

## ----
## Author:  C.A. (Robert) Warmerdam
## Email:   c.a.warmerdam@umcg.nl
##
## Copyright (c) C.A. Warmerdam, 2023
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## A copy of the GNU General Public License can be found in the LICENSE file in the
## root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
## ----

# Load libraries
library(data.table)
library(tidyverse)
library(rtracklayer)

# Declare constants

# Declare function definitions

# Main
read_gen_file <- function(gen_file = NULL) {
  if (is.null(gen_file)) {
    gen_file <- "genotype_probabilities_incl_deletions_2023-03-23.gen"
  }
  return(fread(gen_file))
}

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv=NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  gen <- read_gen_file(gen_file = argv[1])
  ch <- import.chain("/groups/umcg-fg/tmp01/projects/pgx-passport/tools/PGx-pipeline/pgx-imputation-pipeline/data/GRCh37_to_GRCh38_tabs.chain")

  ranges <- makeGRangesFromDataFrame(gen %>% select(c(chrom = V1, start = V4, end = V4)))
  mapped <- as.data.frame(unlist(liftOver(ranges, ch)))

  gen$V4 <- mapped$start

  fwrite(gen, "genotype_probabilities_incl_deletions_liftover_2023-03-23.gen", quote=F, row.names=F, col.names=F, sep=" ")
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}