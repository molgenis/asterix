#!/usr/bin/env Rscript


# Load libraries

# Declare constants

# Declare function definitions

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Process input
  # Load arguments
  argv2 <- list(
    replacements="/groups/umcg-gdio/tmp01/projects/2021001/rawdata/project/2022-10-04_GSAMD-24v3_PGx_QC/output/5_Relatedness/output/sample_swap_mapping_guide_replacement_20221025.txt",
    removals="/groups/umcg-gdio/tmp01/projects/2021001/rawdata/project/2022-10-04_GSAMD-24v3_PGx_QC/output/5_Relatedness/output/sample_swap_mapping_guide_removal_20221025.txt")

  # Read

  message("[INFO] Reading input files")

  mapping.table.replacements <- fread(
    argv2$replacements, data.table=F,
    colClasses=c("character", "character", "character", "character"))

  print(str(mapping.table.replacements))
  mapping.table.removals <- fread(
    argv2$removals, data.table=F,
    colClasses=c("character", "character", "character", "character"))
  print(str(mapping.table.removals))

  if (any(duplicated(mapping.table.removals$IID))) {
    stop(sprintf("Duplicated sample in '%s'!", argv[3]))
  }

  if (any(duplicated(mapping.table.replacements$IID))) {
    stop(sprintf("Duplicated sample in '%s'!", argv[2]))
  }

  # First perform replacements
  sample.table.corrected <- tibble(IID = sample_list) %>%
    left_join(mapping.table.replacements, by=c("IID"="GENOTYPING_ID_INITIAL")) %>%
    mutate(
      IID_original = IID,
      IID_CORRECTED = case_when(
        !is.na(REPLACEMENT) ~ REPLACEMENT,
        TRUE ~ IID),
      EXCLUDE = str_starts(IID_CORRECTED, "exclude"),
      RENAMED = IID != IID_CORRECTED & !EXCLUDE,
      IID = case_when(!EXCLUDE ~ IID_CORRECTED, TRUE ~ IID))

  sample.table.corrected <- sample.table.corrected %>%
    left_join(mapping.table.removals, by=c("IID"="GENOTYPING_ID_INITIAL")) %>%
    mutate(EXCLUDE2 = !is.na(CONCLUSION) & str_starts(CONCLUSION, "exclude"))

  sample.table.toremove <- tibble(fread(sprintf(
    "%s/overall.samples_to_remove_combined_2023-01-27.fam",
    "/groups/umcg-gdio/tmp01/projects/2021001/rawdata/project/2023-01-25_GSAMD-24v3_PGx_QC/output/7_second_it_preparation"
  ), header=F))

  n.removed <- sum(sample.table.corrected$EXCLUDE)
  n.removed2 <- sum(sample.table.corrected$EXCLUDE2)
  n.renamed <- sum(sample.table.corrected$RENAMED)

  message(sprintf("[INFO] Number of samples renamed: %d", n.renamed))
  message(sprintf("[INFO] Number of samples removed in first mapping: %d", n.removed))
  message(sprintf("[INFO] Number of samples removed in second mapping: %d", n.removed2))

  sample.table.corrected2 <- sample.table.corrected  %>%
    left_join(sample.table.toremove %>% distinct(), by=c("IID"="V1")) %>%
    group_by(IID, IID_original) %>%
    summarise(EXCLUDE_SUMMARISED=any(EXCLUDE,EXCLUDE2,!is.na(V2)))

  probabilities_adjusted_corrected <- probabilities_adjusted %>%
    ungroup() %>%
    left_join(sample.table.corrected2, by = c("Sample_ID" = "IID_original")) %>%
    filter(!EXCLUDE_SUMMARISED) %>% select(-Sample_ID) %>% rename(Sample_ID=IID)
  write_delim(samples, "b_dosage.sample", delim=" ", quote=F)

  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}