#!/usr/bin/env Rscript

# Load libraries
library(Homo.sapiens)
library(biomaRt)
library(biovizBase)
library(ggbio)
library(tidyverse)
library(tidyr)
library(dtplyr)
library(broom)
library(argparse)
library(cowplot)
library(ggrepel)
library(viridis)
library(magick)
library(ggforce)
library(ggtext)
library(gridExtra)
library(devtools)
library(ggrastr)
library(GGally)
library(RColorBrewer)
library(LPCM)
library(extrafont)

# Declare constants
## Argument parser

parser <- ArgumentParser(description = 'Visualize the principal components of SNP intensities in a cnv locus.')
parser$add_argument('eigenvectors', metavar = 'file', type = 'character',
                    help = 'A tab delimited file with principal components of variants in the CNV locus.')

## Plotting defautls
old <- theme_set(theme_minimal(base_size = 12))
theme_update(
  line = element_line(
    colour = "black", size = (1 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  # axis.line = element_line(
  #   colour = "#595A5C", size = (1 / (ggplot2::.pt * 72.27/96)),
  #   linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  # axis.ticks = element_line(
  #   colour = "#595A5C", size = (1 / (ggplot2::.pt * 72.27/96)),
  #   linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  panel.grid.minor = element_blank(),
  text = element_text(family="Lato"),
  title = element_text(colour = "#595A5C", face = "bold")
)

color_vector <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000")



variant_annotation_layer <- function(variant_data_frame_missingness, genotype_labels, missingness_variable, x_limits, y_limits) {

  n_missing_samples <- sum(is.na(variant_data_frame_missingness[, missingness_variable, drop = T]))

  genotype_frequencies <- variant_data_frame_missingness %>%
    mutate(ALL_N = n()) %>%
    group_by(genotype, .drop = F) %>%
    summarise(N = n()) %>%
    mutate(percentage = N / sum(N) * 100,
           genotype_labels = genotype_labels[as.character(genotype)],
           label = sprintf("%s: %d (%.2f%%)", genotype_labels, N, percentage))

  variant_data_frame_missingness %>% select(variant, A, B, Distance)

  variant_annotation <- paste0(c(
      sprintf("%s: %s/%s",
              variant_data_frame_missingness$variant[1],
              variant_data_frame_missingness$A[1],
              variant_data_frame_missingness$B[1]),
      sprintf("delta: %d", variant_data_frame_missingness$Distance[1]),
      sprintf("missing samples: %d", n_missing_samples),
      genotype_frequencies$label), collapse = "\n")

  print(variant_annotation)

  return(annotate(
    geom = "text",
    x = x_limits[2], y = y_limits[2], label = variant_annotation,
    hjust = 1, vjust = 1, size = 3,
    # remove label padding, since we have removed the label outline
    parse = F))
}

genotype_layer <- function(color_values, genotype_labels) {
  geom_point(alpha = 0.36, size = 0.75, colour = .data[["genotype"]]) +
    geom_mark_ellipse(aes(label = genotype, filter = !is.na(genotype)),
                      expand = 0.02) +
    scale_colour_manual(
      values = color_values, name = "Genotype",
      breaks = unname(genotype_labels),
      na.value = "grey50")
}

# Declare function definitions
plot_cnv_variant_intensities <- function(
  variant_data_frame, x_variable, y_variable) {

  allele_labels <- get_allele_labels(variant_data_frame)

  xlimits <- get_axis_limits(variant_data_frame, x_variable, y_variable)
  ylimits <- xlimits

  this_plot <- ggplot(
    variant_data_frame,
    aes(x = .data[[x_variable]], y = .data[[y_variable]])) +
    ylab(sprintf("%s intensity", allele_labels["B"])) +
    xlab(sprintf("%s intensity", allele_labels["A"])) +
    coord_equal(ratio = 1, xlim = xlimits, ylim = ylimits)

  return(this_plot)
}

plot_variant_intensities <- function(combined_data_frame_clustered, col_variable = "cluster_aggregated") {

  variant_data_frame <- combined_data_frame_clustered %>% filter(
    variant == var) %>%
    filter(!is.na(!!!x_variable), !is.na(!!!y_variable))

  this_plot <- variant_data_frame %>%
    ggplot(aes(x = .data[[x_variable]], y = .data[[y_variable]],
               colour = .data[[col_variable]])) +
    # geom_point(alpha = 0.36, fill = "black", shape = 16, size = 0.75) +

    geom_text(aes(label = cluster)) +
    scale_colour_manual(
      values = cluster_col,
      na.value = "grey50") +
    # geom_label_repel(
    #   data = cluster_centers %>% filter(variant == var, cluster_size > 2),
    #   aes(label = cluster)) +
    # geom_richtext(
    #   data = variant_data,
    #   aes(x = X, y = Y, label = label),
    #   hjust = 1, vjust = 1,
    #   size = 3,
    #   # remove label background and outline
    #   fill = NA, label.color = NA,
    #   # remove label padding, since we have removed the label outline
    #   label.padding = grid::unit(rep(0, 4), "pt"),
    #   inherit.aes = F
    # ) +
    ylab(sprintf("%s intensity", allele_labels["B"])) +
    xlab(sprintf("%s intensity", allele_labels["A"])) +
    coord_equal(ratio = 1) +
    theme(
      legend.position = "None",
      legend.justification = c(1.0, 1.0),
      legend.box.background = element_rect(
        fill = "white", colour = "black", size = 1))

  pdf("test123.pdf", width = 3, height = 3, useDingbats = F)
  par(xpd = NA)
  plot(this_plot)
  dev.off()
}

get_axis_limits <- function(variant_data_frame, x_variable, y_variable) {
  xlimits <- c(
    floor(min(variant_data_frame[,x_variable, drop=T], variant_data_frame[,y_variable, drop=T], na.rm = T)),
    max(variant_data_frame[,x_variable, drop=T], variant_data_frame[,y_variable, drop=T], na.rm = T))
  return(xlimits)
}

get_allele_labels <- function(variant_data_frame) {
  allele_labels <- c(
    "A" = unique(variant_data_frame$A),
    "B" = unique(variant_data_frame$B))

  return(allele_labels)
}

get_genotype_labels_from_allele_labels <- function(allele_labels) {
  genotype_labels <- c(
    "AA" = paste0(allele_labels[c("A", "A")], collapse = "/"),
    "AB" = paste0(allele_labels[c("A", "B")], collapse = "/"),
    "BB" = paste0(allele_labels[c("B", "B")], collapse = "/"))

  return(genotype_labels)
}

get_genotype_labels <- function(variant_data_frame) {
  allele_labels <- get_allele_labels(variant_data_frame)
  genotype_labels <- c(
    "AA" = paste0(allele_labels[c("A", "A")], collapse = "/"),
    "AB" = paste0(allele_labels[c("A", "B")], collapse = "/"),
    "BB" = paste0(allele_labels[c("B", "B")], collapse = "/"))

  return(genotype_labels)
}


cluster_variants_in_locus <- function(combined_data_frame, variable_to_cluster) {
  intensity_data_frame <- combined_data_frame %>%
    pivot_wider(
      id_cols = "Sample_ID",
      names_from = "variant",
      values_from = all_of(variable_to_cluster))

  intensity_matrix <- intensity_data_frame %>%
    select(-Sample_ID) %>%
    as.matrix()
  rownames(intensity_matrix) <- intensity_data_frame[, "Sample_ID", drop = T]

  variant_R_correlation <- as_tibble(cor(intensity_matrix)^2, rownames = "variant") %>%
    pivot_longer(cols = -variant, names_to = "variant_other", values_to = "correlation") %>%
    inner_join(variants, by = c("variant" = "Name")) %>%
    filter(Distance == 0) %>%
    mutate(
      variant = ordered(variant, levels = variants$Name[order(variants$Start)]),
      variant_other = ordered(variant_other, levels = variants$Name[order(variants$Start)]))

  corr_plot <- ggplot(variant_R_correlation, aes(variant, variant_other, fill = correlation)) +
    geom_tile() +
    scale_fill_viridis(discrete = FALSE) +
    coord_fixed()

  pdf(file_name("", "pdf"),
      useDingbats = F, width = 16, height = 16)

  par(xpd = NA)

  plot(corr_plot)
  #heatmap(cor(intensity_matrix)^2, scale = "none", cexRow=0.5, cexCol = 0.5)

  dev.off()

}

plot_cnv_variant_intensities_with_different_genotype <- function(
  data, x_variable, y_variable, colouring_variant, intensity_variant) {

  intensity_variant_data_frame <- data %>%
    filter(
      variant == intensity_variant, !is.na(!!!x_variable), !is.na(!!!y_variable)) %>%
    select(Sample_ID, variant, A, B, Start, End, Distance, all_of(c(x_variable, y_variable)))

  geno_variant_data_frame <- data %>%
    filter(
      variant == colouring_variant) %>%
    select(Sample_ID, variant, A, B, Start, End, Distance, genotype)

  xlimits <- get_axis_limits(intensity_variant_data_frame, x_variable, y_variable)
  ylimits <- xlimits

  geno_allele_labels <- get_allele_labels(geno_variant_data_frame)
  geno_genotype_labels <- get_genotype_labels_from_allele_labels(geno_allele_labels)
  intensity_allele_labels <- get_allele_labels(intensity_variant_data_frame)

  genotype_frequencies <- geno_variant_data_frame %>%
    mutate(ALL_N = n()) %>%
    group_by(genotype, .drop = F) %>%
    summarise(N = n()) %>%
    mutate(percentage = N / sum(N) * 100,
           genotype_labels = geno_genotype_labels[as.character(genotype)],
           label = sprintf("%s: %d (%.2f%%)", genotype_labels, N, percentage))

  geno_variant_data_frame <- geno_variant_data_frame %>%
    mutate(genotype = recode_factor(genotype, !!!geno_genotype_labels, .ordered = TRUE))

  variant_data_frame <- intensity_variant_data_frame %>%
    inner_join(
      geno_variant_data_frame, by = "Sample_ID",
      suffix = paste0("_", c(
        "int", "geno")
      )
    )

  # geno_variant_data <- geno_variant_data_frame %>%
  #   select(variant, A, B, Distance) %>% distinct() %>%
  #   mutate(label = paste(c(
  #     sprintf("int: %s", !!!intensity_variant),
  #     sprintf("%s: %s/%s", variant, A, B),
  #     sprintf("&Delta;bp: %d", Distance),
  #     genotype_frequencies$label), collapse = "<br>"))
  #
  # intensity_variant_data <- intensity_variant_data_frame %>%
  #   select(variant, A, B, Distance) %>% distinct() %>%
  #   mutate(label = paste(c(
  #     sprintf("%s: %s/%s", variant, A, B),
  #     sprintf("&Delta;bp: %d", Distance)), collapse = "<br>"),
  #          X = xlimits[2], Y = ylimits[2])
  #
  # print(geno_variant_data)
  # print(intensity_variant_data)
  # print(full_join(geno_variant_data, intensity_variant_data,
  #                 by = character(),
  #                 suffix = paste0("_", c("int", "geno"))))
  #
  # variant_data <- full_join(geno_variant_data, intensity_variant_data,
  #                           by = character(),
  #                           suffix = paste0("_", c("int", "geno"))) %>%
  #   mutate(label = sprintf(
  #     "(int) %s<br>(geno) %s<br>",
  #     label_int, label_geno))

  this_plot <- variant_data_frame %>%
    ggplot(aes(x = .data[[x_variable]], y = .data[[y_variable]], colour = genotype)) +
    geom_point(alpha = 0.36, shape = 16, size = 0.75) +
    ylab(sprintf("%s intensity", intensity_allele_labels["B"])) +
    xlab(sprintf("%s intensity", intensity_allele_labels["A"])) +
    # geom_mark_ellipse(aes(label = genotype, filter = !is.na(genotype)),
    #                   expand = 0.02) +
    scale_colour_manual(
      values = plasma(length(geno_genotype_labels), alpha = 1, begin = 0.1, end = 0.9, direction = 1),
      name = "Genotype",
      breaks = unname(geno_genotype_labels),
      na.value = "grey50") +
    # geom_richtext(
    #   data = variant_data,
    #   aes(x=X, y=Y, label = label),
    #   hjust = 1, vjust = 1,
    #   size = 3,
    #   # remove label background and outline
    #   fill = NA, label.color = NA,
    #   # remove label padding, since we have removed the label outline
    #   label.padding = grid::unit(rep(0, 4), "pt"),
    #   inherit.aes = F
    # ) +
    coord_equal(ratio = 1, xlim = xlimits, ylim = ylimits) +
    theme(
      legend.position = "None",
      legend.justification = c(1.0, 1.0),
      legend.box.background = element_rect(
        fill = "white", colour = "black", size = 1))

  return(this_plot)
}

mean_shift_clustering <- function(principal_intensity_matrix, bandwidth = 0.10) {
  result <- ms(
    principal_intensity_matrix,
    h = bandwidth,
    scaled = 1,
    plot = FALSE)

  # assignment
  lpcm_assignment <- tibble(
    Sample_ID = rownames(principal_intensity_matrix),
    cluster = as.character(result$cluster.label))

  min_result <- min(12 + 1, length(unique(result$cluster.label)))

  lpcm_assignment_aggregated <- lpcm_assignment %>%
    group_by(cluster) %>%
    summarise(cluster_size = n()) %>%
    ungroup() %>%
    mutate(
      cluster_size_threshold = sort(cluster_size, decreasing = T)[min_result],
      cluster_size_rank = rank(cluster_size),
      cluster_aggregated = case_when(
        cluster_size < 3 ~ NA_character_,
        cluster_size_threshold <= cluster_size ~ cluster,
        TRUE ~ NA_character_))

  lpcm_assignment_summary <- lpcm_assignment_aggregated %>%
    inner_join(lpcm_assignment, by = "cluster")

  cluster_centers <- result$cluster.center
  colnames(cluster_centers) <- colnames(principal_intensity_matrix)

  clusters <- as_tibble(cluster_centers, rownames = "cluster") %>%
    inner_join(lpcm_assignment_aggregated, by = "cluster")

  return(list("x" = lpcm_assignment_summary,
              "clusters" = clusters))
}

pca_components <- function(
  converted_data_frame, number_of_pcs = NULL) {

  extended_intensity_matrix <- converted_data_frame %>%
    select(-Sample_ID) %>%
    as.matrix()

  rownames(extended_intensity_matrix) <- converted_data_frame[, "Sample_ID", drop = T]

  extended_intensity_matrix_nonmissing <- extended_intensity_matrix[rowSums(is.na(extended_intensity_matrix)) == 0,]

  if (!is.null(number_of_pcs)) {
    princomp_object <- prcomp(
      extended_intensity_matrix_nonmissing,
      center = F, scale. = T, rank. = number_of_pcs)
  } else {
    princomp_object <- prcomp(
      extended_intensity_matrix_nonmissing,
      center = F, scale. = T)
  }

  return(princomp_object)
}

variant_ids <- function(chr, start, snp_mart) {
  variant_data_frame <- tibble(chr, start)
  variant_data_frame_distinct <- variant_data_frame %>% distinct()

  variant_temp <- getBM(
    attributes = 'refsnp_id',
    filters = 'chromosomal_region',
    values = paste0(variant_data_frame_distinct$chr, ":",
                  variant_data_frame_distinct$start, ":",
                  variant_data_frame_distinct$start), mart = snp_mart)

  return(variant_data_frame %>% inner_join(temp) %>% pull(refsnp_id))
}

normalize_component <- function(x, guide) {
  correlation_statistics <- cor.test(x, guide, na.action = na.omit)
  if (correlation_statistics$p.value > 0.5) {
    return(x)
  }
  return(sign(correlation_statistics$estimate) * x)
}

pca_for_group <- function(data_frame) {

  positions <- data_frame %>%
    select(Chrom, Start) %>%
    distinct()

  if (n_distinct(positions) > 1) {
    stop("Positions not unique")
  }

  variable_filtered_intensity_data_frame <- data_frame %>%
    filter(type == "R_corrected") %>%
    pivot_wider(id_cols = "Sample_ID",
                names_from = "var",
                values_from = "value")

  variable_filtered_theta_data_frame <- data_frame %>%
    filter(type == "theta") %>%
    pivot_wider(id_cols = "Sample_ID",
                names_from = "var",
                values_from = "value")

  combined_dataframe_intensities <- as_tibble(
    pca_components(variable_filtered_intensity_data_frame)$x,
    rownames = "Sample_ID") %>%
    pivot_longer(
      cols = -Sample_ID, values_to = "R_intensity", names_to = "component") %>%
    full_join(positions, by = as.character())

  combined_dataframe <- as_tibble(
    pca_components(variable_filtered_theta_data_frame)$x,
    rownames = "Sample_ID") %>%
    pivot_longer(
      cols = -Sample_ID, values_to = "theta", names_to = "component") %>%
    full_join(positions, by = as.character()) %>%
    inner_join(combined_dataframe_intensities, by =
      c("Sample_ID", "component", "Chrom", "Start"))

  return(combined_dataframe)
}


lda_for_group <- function(x) {

  positions <- x %>%
    select(Chrom, Start) %>%
    distinct()

  if (n_distinct(positions) > 1) {
    stop("Positions not unique")
  }

  variable_filtered_intensity_data_frame <- x %>%
    pivot_wider(id_cols = c("Sample_ID", "marker_genotype"),
                names_from = "variant",
                values_from = c("R", "theta_explained_R")) %>%
    filter(across(
      .cols = c(starts_with("R"), starts_with("theta_explained_R")),
      .fns = ~ !is.na(.x)))

  # print(variable_filtered_intensity_data_frame)

  lda_in <- variable_filtered_intensity_data_frame %>% select(-Sample_ID)
  # print(lda_in)

  combined_dataframe_intensities <- as_tibble(
    predict(MASS::lda(marker_genotype ~ ., data = lda_in))$x) %>%
    mutate(Sample_ID = variable_filtered_intensity_data_frame[["Sample_ID"]]) %>%
    full_join(positions, by = as.character())

  return(combined_dataframe_intensities)
}


expected_clusters <- function(x_mat, theoretical_cnv_frequencies) {
  cumulative_cnv_frequencies <- c(0, cumsum(theoretical_cnv_frequencies$freq_g))

  initial_guess_per_feature <- apply(x_mat, 2, function(x) {
    return(as.numeric(as.character(cut(
      x,
      breaks = quantile(x, probs = cumulative_cnv_frequencies),
      labels = theoretical_cnv_frequencies$marker_g,
      include.lowest = T,
      ordered_result = T))))
  })

  mean_initial_guess <- apply(initial_guess_per_feature, 1, function(x) {
    sample(x, size = 1)
  })

  initial_guess_means <- apply(x_mat, 2, function(x) {
    sapply(
      theoretical_cnv_frequencies$marker_g,
      function(level) {
        mean(x[mean_initial_guess == level], na.rm = T)
      })
  })

  return(list(cluster = mean_initial_guess, centroid = initial_guess_means))
}

copy_number_dosage <- function(positions_matrix) {
  expected_cnv_frequencies <- cnv_frequencies()

  dosages <- apply(positions_matrix, 2, function(x) {

    init <- expected_clusters(matrix(x), expected_cnv_frequencies)
    cnv_mod <- normalmixEM(
      x,
      lambda = expected_cnv_frequencies$freq_g,
      mu = init$centroid,
      maxit = 500)

    dosage <- apply(cnv_mod$posterior, 1, function(x) sum(x * c(0, 1, 2, 3, 4)))
  })

  rownames(dosages) <- rownames(positions_matrix)

  return(dosages)
}

clustering <- function(data_frame, iter_n) {
  # Cluster the intensities (possibly for multiple assays).
  # - the number of clusters should be chosen based on the
  # - clustering performance, and the Hardy Weinberg principle.

  # Loop through a number of clusters, 1:6 (or so)
  # - identify to which degree the frequencies are according to the hardy weinberg principle
  # Pick the number of clusters that best adheres to this hardy weinberg principle.

  mapply(function(n) {
    # Perform clustering strategy


  }, iter_n, SIMPLIFY = F, USE.NAMES = F)
}

copy_number_makeup <- function(copy_number) {
  base <- (copy_number %/% 2) - 1
  return(c(base, base + copy_number %% 2))
}

calculate_cnv_hwe_fit <- function(copies, f) {
  # 0 -> -1, -1
  # 1 -> -1, 0 | 0, -1

  # Now calculate based on the following input
  # c(0 = 10, 1 = 20, 2 = 100, 3 = 10
  # What the frequencies would be of the copy number alleles

  allele_frequencies <- table(copy_number_makeup(rep(copies, f))) / (sum(f) * 2)

  expected_frequencies <- outer(allele_frequencies, allele_frequencies)



}

cnv_frequencies <- function(star5 = 0.0295, duplication = 0.0275) {
  no_cnv <- 1 - star5 - duplication
  allele_freq_table <- tibble(freq = c(star5, no_cnv, duplication),
         marker = c(-1, 0, 1))

  geno_freq_table <- full_join(
    allele_freq_table, allele_freq_table, by = as.character()) %>%
    mutate(freq_g = freq.x * freq.y,
           marker_g = marker.x + marker.y) %>%
    group_by(marker_g) %>%
    summarise(freq_g = sum(freq_g))

  return(geno_freq_table)
}


plot_principal_component_correlations <- function(pc_intensities) {

  correlation_matrix <- cor(pc_intensities %>%
                              select(-Sample_ID) %>%
                              as.matrix())^2

  pc_labels <- sapply(
    colnames(correlation_matrix),
    function(column_name) startsWith(column_name, "PC"))
  variants_labels <- sapply(
    colnames(correlation_matrix),
    function(column_name) !startsWith(column_name, "PC"))

  correlation_data_frame <- as_tibble(
    correlation_matrix[variants_labels, pc_labels],
    rownames = "variant") %>%
    pivot_longer(cols = -variant, names_to = "components", values_to = "correlation") %>%
    mutate(
      components = ordered(components, levels = paste0("PC", seq_along(unique(components)))))

  corr_plot <- ggplot(correlation_data_frame, aes(components, variant, fill = correlation)) +
    geom_tile() +
    scale_fill_viridis(discrete = FALSE) +
    coord_fixed()

  pdf(file.path(".", sprintf(
    "out.pc_correlation.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
      useDingbats = F, width = 16, height = 16)

  par(xpd = NA)

  plot(corr_plot)

  dev.off()
}


plot_cnv_correlations <- function(data_frame, pheno) {

  converted_data_frame %>%
    filter(!is.na(CNV_intron_2)) %>%
    group_by(variant, type) %>%
    summarise(
      correlation_CNV_intron_2 = cor.test(value, CNV_intron_2, method = "pearson")$p.value,
      correlation_CNV_intron_6 = cor.test(value, CNV_intron_6, method = "pearson")$p.value,
      correlation_CNV_exon_9 = cor.test(value, CNV_exon_9, method = "pearson")$p.value) %>%
    pivot_longer(
      cols = starts_with("correlation"),
      names_pattern = "(correlation)_(.+)",
      names_to = c(".value", "CNV"),
      values_to = "correlation") -> correlation_data_frame

  corr_plot <- ggplot(correlation_data_frame, aes(CNV, variant, fill = correlation)) +
    geom_tile() +
    scale_fill_viridis(discrete = FALSE) +
    coord_fixed() +
    facet_wrap(~ type, ncol = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  pdf(file.path(".", sprintf(
    "out.cnv_correlation.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
      useDingbats = F, width = 16, height = 16)

  par(xpd = NA)

  plot(corr_plot)

  ggplot(correlation_data_frame, aes(x=-log10(correlation))) +
    geom_histogram(bins = 50)

  dev.off()
}


plot_principal_components <- function(pc_intensities_clustered, pcs, cluster_col, cluster_centers) {

  png(file.path(".", sprintf(
    "out.pc_projections.%s.png", format(Sys.Date(), "%Y%m%d"))),
  width = 3840, height = 3840, res = 300)

  par(xpd = NA)

  this_plot <- pc_intensities_clustered %>%
    ggpairs(
      columns = pcs,
      mapping = aes(colour = cluster_aggregated),
      upper = list(continuous = wrap("density", size = 0.75, alpha = 0.36)),
      diag = list(continuous = wrap("densityDiag", alpha = 0.36, size = 0)),
      lower = list(continuous = wrap("points", alpha = 0.36, shape = 16, size = 0.75))) +
    scale_colour_manual(values = cluster_col, na.value = "grey50") +
    scale_fill_manual(values = cluster_col, na.value = "grey50") +
    theme(
      legend.position = "None",
      legend.justification = c(1.0, 1.0),
      legend.box.background = element_rect(
        fill = "white", colour = "black", size = 1))
  #
  # this_plot[3, 2] <- pc_intensities_clustered %>%
  #   ggplot(aes(x=PC2, y=PC3, colour=cluster_aggregated)) +
  #   geom_point(alpha = 0.36, shape = 16, size = 0.75) +
  #   scale_colour_manual(values = cluster_col, na.value = "grey50") +
  #   geom_label_repel(data = cluster_centers %>% filter(cluster_size > 2), aes(label=cluster),
  #                    box.padding = 0.5, max.overlaps = Inf)

  print(this_plot)

  dev.off()
}


read_intensities <- function(cnv_path_prefix) {

  intensities_corrected <- read_csv(
    paste(cnv_path_prefix, "intensity_correction.corrected.csv.gz", sep = "."),
    col_types = cols(`Sample ID` = col_character())) %>%
    rename("Sample_ID" = `Sample ID`) %>%
    pivot_longer(cols = -Sample_ID, values_to = "R_corrected", names_to = "variant")
  intensities_data_frame <- read_tsv(
    paste(cnv_path_prefix, "intensity_data_frame.csv.gz", sep = "."),
    col_types = cols("Sample ID" = col_character())) %>%
    rename_all(list(~stringr::str_replace_all(., '\\s', '_'))) %>%
    inner_join(intensities_corrected, by = c("Sample_ID", "variant")) %>%
    group_by(Sample_ID) %>%
    mutate(missingness = case_when(
      any(is.na(X)) ~ "incomplete",
      TRUE ~ "complete")) %>%
    ungroup()

  corrected_intensities_data_frame <- intensities_data_frame %>%
    mutate(X_corrected = X * (R_corrected / R),
           Y_corrected = Y * (R_corrected / R))

  return(corrected_intensities_data_frame)
}

process_intensities <- function(corrected_intensities_data_frame, variant_alleles, variants) {
  corrected_intensities_data_frame %>%
    left_join(variant_alleles, by = "SNP") %>%
    left_join(variants, by = c("variant" = "Name")) %>%
    mutate(
      variant = ordered(variant, levels = variants$Name[order(variants$Start)]),
      genotype = ordered(GType, levels = c("AA", "AB", "BB")),
      theta = atan(X_corrected / Y_corrected) * 180 / pi,
      R_alternative = sqrt(X_corrected^2 + Y_corrected^2)) %>%
    select(-GType)
}

select_positions_within_range <- function(variants, start, end) {
  variants %>% filter(Start >= start & End <= end) %>%
    group_by(Chrom, Start, End) %>% summarise()
}


# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
file_name <- function(name, ext) {
  return(file.path(".", sprintf(
    "out.%s.%s.%s", name, format(Sys.Date(), "%Y%m%d"), ext)))
}

main <- function(argv = NULL) {
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

  value_label <- "R"
  # value_label <- "B_Allele_Freq"
  # value_label <- "Log_R_Ratio"

  principal_components_ugli <- read_csv(
    "/groups/umcg-lifelines/tmp01/projects/ov21_0355/pgx-pipeline/analyses/cyp2d6_cnv_calling/out.intensity_correction.batcheffects.csv.gz",
    col_types = cols("Sample ID" = col_character())) %>%
    rename("Sample_ID" = `Sample ID`)

  principal_components_mom <- read_csv(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/results/mom/analyses/cyp2d6_cnv_calling/out.intensity_correction.batcheffects.csv.gz",
    col_types = cols("Sample ID" = col_character())) %>%
    rename("Sample_ID" = `Sample ID`)

  mom_samples_exclude <- read_delim(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/data/medicatie-op-maat/sample_list.exclude.txt",
    col_names = c("Family_ID", "Sample_ID"),
    col_types = cols(.default = col_character()), delim = " ") %>%
    mutate(exclude = T)

  batch_effects <- bind_rows(
    list("mom" = principal_components_mom, "ugli" = principal_components_ugli),
    .id = "dataset") %>%
    left_join(mom_samples_exclude) %>%
    mutate(exclude = !is.na(exclude))

  eigenvalues <- read_csv("/groups/umcg-lifelines/tmp01/projects/ov21_0355/pgx-pipeline/analyses/cyp2d6_cnv_calling/out.intensity_correction.eigenvalues.csv.gz")

  loadings <- read_csv("/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/results/mom/analyses/cyp2d6_cnv_calling/20220602/out.intensity_correction.pca.csv.gz") %>%
    rename(Name = "...1") %>%
    rename_with(.fn = ~ paste0("PC", .x), .cols = -Name) %>%
    inner_join(corrective_variants) %>%
    pivot_longer(cols = starts_with("PC"), names_to = "component",
                 values_to = "loading")

  r_matrix <- read_tsv("/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/results/mom/analyses/cyp2d6_cnv_calling/out.mat.R.csv") %>%
    filter(variant %in% corrective_variants$Name)
  r_as_matrix <- r_matrix %>%
    select(-variant) %>% as.matrix()
  rownames(r_as_matrix) <- r_matrix$variant

  princomp_object <- prcomp(
    t(r_as_matrix),
    center = T, scale. = T, rank. = 100)

  batch_effects_2 <- tibble(princomp)

  intensity_data_frame <- read_tsv(
    paste("/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/results/mom/analyses/cyp2d6_cnv_calling/out", "intensity_data_frame.csv.gz", sep = "."))

  pdf(file.path(".", sprintf(
    "out.principal_component_percent_variance_explained.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
      useDingbats = FALSE, width = 4, height = 4)

  par(xpd = NA)

  ggplot(mytibble, aes(y=percent,x=comp)) +
    geom_line()

  dev.off()


  pdf(file.path(".", sprintf(
    "out.sample_batch_effects.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
      useDingbats = FALSE, width = 8, height = 8)

  par(xpd = NA)

  ggplot(batch_effects %>% filter(dataset == "mom"), aes(`0`, `1`, colour = dataset, shape = exclude)) +
    geom_density2d() +
    geom_point(data = batch_effects %>% filter(dataset == "mom"), alpha = 0.2)

  ggplot(batch_effects %>% filter(dataset == "mom"), aes(`2`, `3`, colour = dataset, shape = exclude)) +
    geom_density2d() +
    geom_point(data = batch_effects %>% filter(dataset == "mom"), alpha = 0.2)

  ggplot(batch_effects %>% filter(dataset == "mom"), aes(`4`, `5`, colour = dataset, shape = exclude)) +
    geom_density2d() +
    geom_point(data = batch_effects %>% filter(dataset == "mom"), alpha = 0.2)

  ggplot(batch_effects %>% filter(dataset == "mom"), aes(`6`, `7`, colour = dataset, shape = exclude)) +
    geom_density2d() +
    geom_point(data = batch_effects %>% filter(dataset == "mom"), alpha = 0.2)

  dev.off()
  #
  # png(file.path(".", sprintf(
  #   "out.sample_correction.%s.png", format(Sys.Date(), "%Y%m%d"))))
  #
  # par(xpd = NA)
  #
  # ggplot(intensities_data_frame, aes(R, R_corrected)) +
  #   geom_point(alpha = 0.2)
  #
  # dev.off()

  variant_alleles <- tibble(
    SNP = c("[T/C]", "[A/G]", "[G/C]", "[I/D]", "[T/G]",
            "[A/C]", "[C/G]", "[A/T]", "[D/I]"),
    A = c("T", "A", "G", "I", "T", "A", "C", "A", "D"),
    B = c("C", "G", "C", "D", "G", "C", "G", "T", "I"))

  cyp2d6loc <- read_tsv("/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/tools/PGx-passport-pilot/data/cyp2d6.bed",
                        col_names = c("Chrom", "Start", "End", "Name"))

  cyp2d6_exons <- read_csv(file.path(
    "~/Documents/projects/pgx-pipeline/tools/PGx_passport_pilot/data/",
    "ExonsSpreadsheet-Homo_sapiens_Transcript_Exons_ENST00000645361.csv")) %>%
    select(c("Exon" = "No."),  "Start", "End")

  cyp2d6_exons <- read_csv(file.path(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/tools/PGx-passport-pilot/data/",
    "ExonsSpreadsheet-Homo_sapiens_Transcript_Exons_ENST00000360608.csv")) %>%
    select(c("Exon" = "No."),  "Start", "End")

  cyp2d6_exons <- read_csv(file.path(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/tools/PGx-passport-pilot/data/",
    "ExonsSpreadsheet-Homo_sapiens_Transcript_Exons_ENST00000645361.csv")) %>%
    select(c("Exon" = "No."),  "Start", "End")

  variants <- read_tsv(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/results/ugli/analyses/select_corrective_variants/out.locus.bed",
    col_names = c("Chrom", "Start", "End", "Name")) %>%
    mutate(Name = ordered(Name, levels = .$Name[order(.$Start)])) %>%
    rowwise() %>%
    mutate(StartDistance = min(0, Start - cyp2d6loc[1, "Start", drop = T]),
           EndDistance = max(0, End - cyp2d6loc[1, "End", drop = T]),
           Distance = c(StartDistance, EndDistance)[which.max(abs(c(StartDistance, EndDistance)))]) %>%
    ungroup()

  corrective_variants <- read_tsv(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/results/ugli/analyses/select_corrective_variants/out.corrective.bed",
    col_names = c("Chrom", "Start", "End", "Name"))

  phenotypes_mom <- read_csv(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/data/medicatie-op-maat/GenotypePhenotypeData_MomCohort_pseudonimized.csv",
    col_types = cols("ID_coded" = col_character())) %>%
    filter(Genotype=="CYP2D6") %>% select(Haplotype, c("Sample_ID" = "ID_coded"))

  cnv_phenotypes_mom <- read_tsv(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/data/medicatie-op-maat/cnv-data.tsv",
    col_types = cols("SampleID_GSA" = "c", .default = "i")) %>%
    select(c(
      "Sample_ID" = "SampleID_GSA",
      "CNV_intron_2" = "CNV intron 2",
      "CNV_intron_6" = "CNV intron 6",
      "CNV_exon_9" = "CNV exon 9"))

  mom_samples_mistakes <- as.character(c(93, 29, 30, 24, 20, 147, 23, 108))

  corrected_intensities_data_frame_ugli <- read_intensities(
    "/groups/umcg-lifelines/tmp01/projects/ov21_0355/pgx-pipeline/analyses/cyp2d6_cnv_calling/out")
  combined_data_frame_ugli <- process_intensities(
    corrected_intensities_data_frame_ugli, variant_alleles, variants)

  corrected_intensities_data_frame_mom <- read_intensities(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/results/mom/analyses/cyp2d6_cnv_calling/out")
  combined_data_frame_mom <- process_intensities(
    corrected_intensities_data_frame_mom, variant_alleles, variants) %>%
    filter(!(Sample_ID %in% mom_samples_exclude$Sample_ID)) %>%
    inner_join(cnv_phenotypes_mom, by = "Sample_ID")

  corrected_intensities_data_frame_mom_2 <- read_intensities(
    "/groups/umcg-wijmenga/tmp01/projects/pgx-pipeline/results/mom/analyses/cyp2d6_cnv_calling/20220602/out")
  combined_data_frame_mom_2 <- process_intensities(
    corrected_intensities_data_frame_mom_2, variant_alleles, variants)

  combined_data_frame <- bind_rows(list(
    "ugli" = combined_data_frame_ugli,
    "mom" = combined_data_frame_mom), .id = "dataset")

  # genotypes <- read_tsv("/groups/umcg-lifelines/tmp01/projects/ov21_0355/pgx-pipeline/analyses/cyp2d6_cnv_calling/out_variants.raw") %>%
  #   pivot_longer(cols = c(-FID, -IID, -PAT, -MAT, -SEX, -PHENOTYPE), names_to = "NameRaw", values_to = "Dosage")
  #
  # genotypes_with_variants <- as_tibble(str_match(unique(genotypes$NameRaw), "^(.+)_([agctAGCT]+)\\(/([agctAGCT]+)\\)$")) %>%
  #   rename(c("NameRaw" = "V1", "Name" = "V2", "Ref" = "V3", "Alt" = "V4")) %>%
  #   rowwise() %>%
  #   mutate(HomRef = paste0(Ref, Ref), Het = paste0(Ref, Alt), HomAlt = paste0(Alt, Alt)) %>%
  #   ungroup() %>%
  #   inner_join(genotypes, by = "NameRaw") %>%
  #   mutate(Dosage = factor(round(Dosage),
  #                            levels = c(0, 1, 2),
  #                            ordered = T))



  # Check if the eigenvector matrix contains only doubles
  # if (!is.double(eigenvectors)) {
  #     stop("Colnames are not all doubles")
  # }

  # intensity_tibble <- intensities_mat %>% inner_join(variants, by=c("variant" = "Name")) %>%
  #   pivot_longer(cols=c(-variant, -Chrom, -Start, -End), names_to="Sample_ID", values_to=value_label) %>%
  #   group_by(variant, Chrom, Start, End) %>%
  #   summarise(variance = var(get(value_label), na.rm = TRUE),
  #             variance_no_missingness = var(get(value_label), na.rm = FALSE)) %>%
  #   pivot_longer(cols=c(variance, variance_no_missingness), names_to="variance_type", values_to="variance_col")

  ## Plot PCs
  #ggplot(eigenvectors, aes())
  #
  # custom_x_scale <- scale_x_continuous(limits = c(42522502 - 100, 42522502 + 4382 + 100))
  #
  # data(genesymbol, package = "biovizBase")
  # geneModel <- ggplot() +
  #   geom_vline(data = variants, aes(xintercept=Start), colour = "#D13B00", inherit.aes = F, alpha = 0.5) +
  #   geom_alignment(data = Homo.sapiens, which = genesymbol["CYP2D6"], inherit.aes = F,
  #                  colour = "grey10", fill = "white", gap.geom = "chevron",
  #                  stat = "identity", cds.rect.h = 0.1) +
  #   ylab("CYP2D6 Gene model") +
  #   custom_x_scale
  #
  # variantLabels <- ggplot(data = variants, aes(x=Start, y=1, label=Name)) +
  #   geom_text_repel(
  #     force_pull = 0, # do not pull toward data points
  #     nudge_y = 0.05,
  #     direction = "x",
  #     angle = 90,
  #     hjust = 0,
  #     segment.size = 0.2,
  #     max.iter = 1e4, max.time = 1
  #   ) + custom_x_scale +
  #   ylim(1, 0.8) +
  #   theme(
  #     axis.line  = element_blank(),
  #     axis.ticks = element_blank(),
  #     axis.text = element_blank(),
  #     axis.title = element_blank(),
  #     panel.spacing.x = unit(.0, "cm")
  #   )
  #
  # pdf(file.path(".", sprintf(
  #   "out.variants_%s.pdf", format(Sys.Date(), "%Y%m%d"))),
  #     useDingbats = FALSE, width = 8, height = 5)
  #
  # par(xpd = NA)
  #
  # plot_grid(variantLabels, geneModel, ncol = 1, align='v', rel_heights=c(2,1))
  #
  # dev.off()
  #
  # intensities_plot <- ggplot(intensity_tibble) +
  #     geom_point(aes(x = Start, y = variance_col, col = variance_type), alpha = 0.5) +
  #     ylab(sprintf("Variance of %s", value_label)) +
  #     custom_x_scale +
  #     theme(legend.position="top")
  #
  # pdf(file.path(".", sprintf(
  #   "out.%s_variance_%s.pdf", value_label, format(Sys.Date(), "%Y%m%d"))),
  #     useDingbats = FALSE, width = 8, height = 6)
  #
  # par(xpd = NA)
  #
  # plot_grid(intensities_plot, geneModel, ncol = 1, align='v')
  #
  # dev.off()

  correlation_data_frame <- combined_data_frame %>%
    select(variant, Sample_ID, R_corrected, theta) %>%
    pivot_longer(c(R_corrected, theta), names_to = "type", values_to = "value", values_drop_na = T) %>%
    filter(!is.na(value)) %>%
    inner_join(x = ., y = ., by = "Sample_ID", suffix = c("_a", "_b")) %>%
    group_by(variant_a, variant_b, type_a, type_b) %>%
    summarise(tidy(cor.test(value_a, value_b))) %>%
    inner_join(variants, by = c("variant_a" = "Name"))

  saveRDS(correlation_data_frame, file_name("variant_correlation", "rds"))
  correlation_data_frame <- readRDS(file_name("variant_correlation", "rds"))

  #
  # corr_plot <- correlation_data_frame_xy %>%
  #   filter(Distance == 0) %>%
  #   ggplot(aes(variant_a, variant_b, fill = abs(estimate))) +
  #   geom_tile() +
  #   scale_fill_viridis(discrete = FALSE) +
  #   coord_fixed() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  #   facet_grid(~ type_a + type_b,
  #              labeller = labeller(type_a = label_both, type_b = label_both))

  pdf(file_name("variant_correlation", "pdf"),
      useDingbats = F, width = 16, height = 16)

  par(xpd = NA)

  plot(corr_plot)
  #heatmap(cor(intensity_matrix)^2, scale = "none", cexRow=0.5, cexCol = 0.5)

  dev.off()

  glued_variant_names <- correlation_data_frame %>%
    ungroup() %>%
    mutate(var_a = str_glue("{variant_a}_type_{type_a}", variant_a = variant_a, type_a = type_a),
           var_b = str_glue("{variant_b}_type_{type_b}", variant_b = variant_b, type_b = type_b))

  var_cor_data_frame <- glued_variant_names %>%
    pivot_wider(id_cols = "var_a",
                names_from = "var_b",
                values_from = "estimate")

  var_cor_matrix <- var_cor_data_frame %>% select(-var_a) %>% as.matrix()
  rownames(var_cor_matrix) <- var_cor_data_frame$var_a

  var_dist_matrix <- as.dist(1 - abs(var_cor_matrix))

  hc1 <- hclust(var_dist_matrix, method = "complete" )

  plot(hc1, cex = 0.6, hang = -1)

  hc2 <- agnes(var_dist_matrix, method = "complete")

  # methods to assess
  m <- c( "average", "single", "complete", "ward")
  names(m) <- c( "average", "single", "complete", "ward")

  # function to compute coefficient
  ac <- function(x) {
    agnes(var_dist_matrix, method = x)$ac
  }

  map_dbl(m, ac)
  ##   average    single  complete      ward
  ## 0.7379371 0.6276128 0.8531583 0.9346210

  hc3 <- agnes(var_dist_matrix, method = "ward")
  pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes")

  hc4 <- diana(var_dist_matrix)
  hc4$dc

  pltree(hc4, cex = 0.6, hang = -1, main = "Dendrogram of diana")

  hc5 <- hclust(var_dist_matrix, method = "ward.D2" )
  #plot(hc5, cex = 0.5, hang = -1, main = "Dendrogram of ward.D2")

  sub_grp <- cutree(hc1, k = 5)

  #rect.hclust(hc1, k = 8, border = 2:5)

  clustered_variant_tibble <- glued_variant_names %>%
    inner_join(tibble(cluster_a = sub_grp, var_a = names(sub_grp))) %>%
    inner_join(tibble(cluster_b = sub_grp, var_b = names(sub_grp)))

  corr_plot <- clustered_variant_tibble %>%
    filter(Distance == 0) %>%
    ggplot(aes(variant_a, variant_b, fill = as.factor(cluster_b))) +
    geom_tile() +
    scale_fill_viridis(discrete = TRUE) +
    coord_fixed() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(~ type_a + type_b,
               labeller = labeller(type_a = label_both, type_b = label_both))

  pdf(file.path(".", sprintf(
    "out.variant_clusters.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
      useDingbats = F, width = 16, height = 16)

  par(xpd = NA)

  plot(corr_plot)
  #heatmap(cor(intensity_matrix)^2, scale = "none", cexRow=0.5, cexCol = 0.5)

  dev.off()
  #
  # vars <- c("seq-rs1135824",
  #           "DUP-seq-rs1135824",
  #           "chr22-42525044:rs1135824")
  # # "seq−rs5030867",
  # # "rs79292917",
  # # "rs72549349")
  #
  # vars <- c(
  #   "seq-rs72549346",
  #   "TOP2-rs1058172",
  #   "TOP-rs1058172"
  #   # "rs1058164",
  #   # "rs28360521",
  #   # "GSA-rs4453786",
  #   # "rs56400210"
  # )
  # intensity <- c("rs199722016", "DUP-rs1135840",
  #           "rs1135840", "seq-rs1135836", "rs72549349",
  #           "rs79292917", "seq-rs5030867", "chr22-42525044:rs1135824",
  #           "DUP-seq-rs1135824", "seq-rs1135824", "rs1135823", "rs1135822", "rs267608305", "rs1081003",
  #           "chr22-42525772:rs28371706", "DUP-seq-rs28371706", "seq-rs28371706",
  #           "DUP-rs1065852", "rs1065852", "rs138100349", "rs72549358",
  #           "rs1080989", "rs1080985", "rs1080983", "rs28360521"
  # )

  whole_gene_set <- tibble(variant = c(
    "rs199722016",  "seq-rs1135836",
    "rs28371725", "rs79292917",
    "seq-rs28371706", "chr22-42525772:rs28371706",
    "DUP-seq-rs28371706", "rs267608305", "rs28360521", "rs1135822",
    "rs1135823", "seq-rs1135824", "DUP-seq-rs1135824",
    "chr22-42525044:rs1135824", "seq-rs5030867", "rs72549349",
    "rs1135740", "rs1135840", "DUP-rs1135840")) %>%
    mutate(type = "R_corrected")

  whole_gene_set <- tibble(
    variant = c("rs199722016", "DUP-rs1135840", "rs1135840",
                "seq-rs1135836", "rs28371725", "rs72549349",
                "rs79292917", "seq-rs5030867", "chr22-42525044:rs1135824",
                "DUP-seq-rs1135824", "seq-rs1135824", "rs1135823", "rs1135822",
                "rs267608305", "rs1081003", "chr22-42525772:rs28371706", "DUP-seq-rs28371706",
                "seq-rs28371706", "DUP-rs1065852", "rs1065852", "rs138100349",
                "rs1080989", "rs1080983", "rs28360521")) %>%
    mutate(type = "R_corrected")

  whole_gene_set <- tibble(
    variant = c("rs199722016", "rs72549349",
                "rs79292917", "seq-rs5030867", "chr22-42525044:rs1135824",
                "DUP-seq-rs1135824", "seq-rs1135824", "rs1135823", "rs1135822",
                "rs267608305", "rs1081003", "chr22-42525772:rs28371706", "DUP-seq-rs28371706",
                "seq-rs28371706", "DUP-rs1065852", "rs1065852", "rs138100349",
                "rs1080989", "rs1080983", "rs28360521")) %>%
    mutate(type = "R_corrected")

  exon9_variant_set <- tibble(
    variant = c("DUP-rs1065852", "rs1065852", "rs1080989", "rs1080983", "rs28360521")) %>%
    full_join(tibble(type = c("R_corrected")), by=character())

  sampled_ugli_samples <- sample(unique(combined_data_frame_ugli$Sample_ID), 10000)
  thinned_data_frame <- combined_data_frame %>%
    filter(dataset == "mom" | Sample_ID %in% sampled_ugli_samples)

  converted_data_frame <- thinned_data_frame %>%
    pivot_longer(c(R_corrected, theta), names_to = "type", values_to = "value", values_drop_na = T) %>%
    mutate(var = str_glue("{variant}_type_{type}", variant = variant, type = type))

  variable_filtered_data_frame <- converted_data_frame %>%
    inner_join(whole_gene_set) %>%
    pivot_wider(id_cols = c("Sample_ID"),
                names_from = "var",
                values_from = "value")

  pca_result <- pca_components(variable_filtered_data_frame, 16)
  principal_intensity_matrix <- pca_result$x

  pc_intensities <- as_tibble(
    principal_intensity_matrix,
    rownames = "Sample_ID") %>%
    inner_join(variable_filtered_data_frame, by = "Sample_ID")

  plot_principal_component_correlations(pc_intensities)

  lpcm_results <- mean_shift_clustering(principal_intensity_matrix[,paste0("PC", 1:3)], bandwidth = 0.12)
  lpcm_assignment <- lpcm_results$x
  lpcm_clusters <- lpcm_results$clusters %>%
    filter(cluster %in% lpcm_assignment$cluster)

  principal_component_means <- colMeans(pca_result$x)
  cluster_centers_principal_component_space <- as.matrix(
    lpcm_clusters %>% select(starts_with("PC")) %>%
      add_column(!!!principal_component_means[!names(principal_component_means) %in% names(.)]))

  cluster_names <- na.omit(lpcm_clusters$cluster_aggregated)
  cluster_col <- colorRampPalette(brewer.pal(n = 8, name = 'Set2'))(length(na.omit(cluster_names)))
  names(cluster_col) <- cluster_names

  pc_intensities_clustered <- pc_intensities %>% inner_join(lpcm_assignment, by = "Sample_ID")

  plot_principal_components(pc_intensities_clustered, paste0("PC", 1:3), cluster_col, lpcm_clusters)

  pdf(file.path(".", sprintf(
    "out.pc_projections.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
      width = 8, height = 8, useDingbats = F)
  par(xpd = NA)

  combined_pc_df <- combined_data_frame %>%
    select(Sample_ID, dataset, CNV_intron_2, CNV_intron_6, CNV_exon_9) %>% distinct() %>%
    inner_join(pc_intensities_clustered, by = "Sample_ID") %>%
    mutate(
      CNV_status = case_when(CNV_intron_2 == CNV_intron_6
                               & CNV_intron_2 == CNV_exon_9 ~ CNV_intron_2,
                             TRUE ~ NA_integer_),
      CNV_label = case_when(
        is.na(CNV_status) & dataset == "mom" ~
          paste(CNV_intron_2, CNV_intron_6, CNV_exon_9, sep = "/")))

  ggplot(combined_pc_df %>% filter(dataset == "ugli"), aes(-PC1, PC2)) +
    geom_point(alpha = 0.36, shape = 16, size = 0.75, aes(colour = -PC1)) +
    scale_colour_viridis_c(
      option = "plasma", name = "PC1", na.value = "grey50") +
    geom_point(data = combined_pc_df %>% filter(dataset == "mom"), shape = 21, colour = "black",
               aes(fill = CNV_status)) +
    # geom_label_repel(data = combined_pc_df %>% filter(dataset == "mom"),
    #                  aes(label = Sample_ID, colour = CNV_status),
    #                  box.padding = 1, max.overlaps = Inf) +
    scale_fill_viridis_c(
      limits = c(0, NA_integer_),
      option = "plasma", name = "CNV,\nstatus", na.value = "grey50")
    scale_colour_viridis_c(
      limits = c(0, NA_integer_),
      option = "plasma", name = "CNV,\nstatus", na.value = "grey50")

  dev.off()
  # clusters_real <- as_tibble(
  #   (cluster_centers_principal_component_space %*% t(pca_result$rotation)) * pca_result$scale,
  #   rownames = "cluster") %>%
  #   pivot_longer(-cluster, names_to = c("variant", ".value"), names_pattern = "(.+)_(R_corrected|theta)") %>%
  #   inner_join(lpcm_clusters, by = "cluster") %>% select(-starts_with("PC")) %>%
  #   mutate(
  #     ratio = tan(theta / (180 / pi)),
  #     X_tmp = R_corrected,
  #     Y_tmp = R_corrected / ratio,
  #     R_tmp = Y_tmp + X_tmp,
  #     X_corrected = X_tmp / R_tmp * R_corrected,
  #     Y_corrected = Y_tmp / R_tmp * R_corrected) %>%
  #   select(-X_tmp, -Y_tmp, -R_tmp, -ratio)

  combined_data_frame_clustered <- combined_data_frame %>%
    inner_join(lpcm_assignment, by = "Sample_ID")

  png_file_list <- base::lapply(levels(droplevels(combined_data_frame %>%
                                         filter(abs(Distance) < 10000) %>%
                                         pull(variant))),
    function(var) {
      png_file_name <- file.path(".", sprintf(
        "out.clustered_intensity_variants_%s_%s.png", var, format(Sys.Date(), "%Y%m%d")))

      this_plot <- plot_cnv_variant_intensities(
        combined_data_frame_clustered, x_variable = "X_corrected", y_variable = "Y_corrected",
        col_variable = "cluster_aggregated", var, cluster_col)
      this_plot_2 <- plot_cnv_variant_intensities(
        combined_data_frame_clustered, x_variable = "X_corrected", y_variable = "Y_corrected",
        col_variable = "genotype", var)

      png(png_file_name, width = 2*540, height = 540, res = 150)

      grid.arrange(this_plot_2, this_plot, ncol = 2)

      par(xpd = NA)

      dev.off()

      return(png_file_name)
    })

  out_file_name <- file.path(".", sprintf(
    "out.clustered_intensity_variants_combined_%s", format(Sys.Date(), "%Y%m%d")))

  images_in <- function(png_file_name) {
    return(png_file_name %>% image_read())
  }

  out <- purrr::map(png_file_list, images_in)
  # image_write(image_append(purrr::lift_dl(c)(out), stack = T),
  #             path = paste(out_file_name, "png", sep = "."), format = "png")

  image_write(image_join(out),
              path = paste(out_file_name, "pdf", sep = "."), format = "pdf")

  # For each bp position with duplicate variants, perform a PCA
  projected_data_frame <- bind_rows(converted_data_frame %>% group_by(Start) %>%
    group_map(~ pca_for_group(.x), .keep = TRUE))

  projected_updated_data_frame <- projected_data_frame %>%
    filter(component == "PC1") %>%
    filter(Start >= cyp2d6loc$Start - 10000 & Start <= cyp2d6loc$End + 10000) %>%
    left_join(cnv_phenotypes_mom, by = "Sample_ID") %>%
    group_by(component, Chrom, Start) %>%
    mutate(
      R_intensity = normalize_component(R_intensity, CNV_intron_6)) %>%
    mutate(
      dataset = case_when(Sample_ID %in% cnv_phenotypes_mom$Sample_ID ~ "mom",
                          TRUE ~ "ugli"),
      CNV_status = case_when(CNV_intron_2 == CNV_intron_6
                               & CNV_intron_2 == CNV_exon_9 ~ CNV_intron_2,
                             TRUE ~ NA_integer_),
      CNV_label = paste0(
        Sample_ID, ": ",
        paste(CNV_intron_2, CNV_intron_6, CNV_exon_9, sep = "/")))

    pdf(file.path(".", sprintf(
      "out.variant_projections2.%s.pdf", format(Sys.Date(), "%Y%m%d"))),
        width = 8, height = 18, useDingbats = F)
    par(xpd = NA)

    # cnv_samples <- projected_updated_data_frame %>%
    #   filter(dataset == "mom", is.na(CNV_status)) %>% pull(Sample_ID)



      ggplot(
        projected_updated_data_frame %>% filter(dataset == "ugli"),
        aes(factor(Start), R_intensity)) +
        rasterize(geom_jitter(
          alpha = 0.18, shape = 16, size = 0.75,
          aes(colour = R_intensity)), dpi = 300) +
        # scale_colour_manual(values = brewer.pal(n = 7, name = 'Set2')) +
        scale_colour_viridis_c(
          option = "plasma", name = "R Intensity", na.value = "grey50") +
        geom_point(
          data = projected_updated_data_frame %>% filter(Sample_ID %in% mom_samples_mistakes),
          shape = 21, colour = "black",
          aes(fill = CNV_status)) +
        geom_label_repel(
          data = projected_updated_data_frame %>% filter(Sample_ID %in% mom_samples_mistakes),
          aes(label = Sample_ID, fill = CNV_status),
          box.padding = 1, max.overlaps = Inf) +
        scale_fill_viridis_c(
          limits = c(0, NA_integer_),
          option = "plasma", name = "CNV,\nstatus", na.value = "grey50") +
        scale_colour_viridis_c(
          limits = c(0, NA_integer_),
          option = "plasma", name = "CNV,\nstatus", na.value = "grey50") +
        theme(legend.position = "bottom", axis.text = element_blank(), axis.line = element_blank()) +
        facet_wrap(~ Start, scales = "free")

    dev.off()

  positions_wide <- projected_data_frame %>%
    filter(component == "PC1") %>%
    filter(Start >= cyp2d6loc$Start - 10000 & Start <= cyp2d6loc$End + 10000) %>%
    left_join(cnv_phenotypes_mom, by = "Sample_ID") %>%
    group_by(component, Chrom, Start) %>%
    mutate(
      R_intensity = normalize_component(R_intensity, CNV_intron_6)) %>%
    pivot_wider(
      id_cols = Sample_ID, values_from = R_intensity, names_from = "Start") %>%
    filter(across(everything(), ~ !is.na(.)))

  positions_matrix <- positions_wide %>%
    select(-Sample_ID) %>%
    as.matrix()

  rownames(positions_matrix) <- positions_wide %>% pull(Sample_ID)


  pca_component <- as.tibble(pca_components(positions_wide)$x, rownames = "Sample_ID") %>%
    select(Sample_ID, PC1, PC2)

  dosages <- copy_number_dosage(positions_matrix)
  copy_number_dosages <- as_tibble(dosages, rownames = "Sample_ID") %>%
    pivot_longer(cols = -Sample_ID,
                 names_to = "Start",
                 values_to = "dosage") %>%
    mutate(Start = as.numeric(Start))







  png_file_list <- base::lapply(combined_data_frame %>%
                                         filter(abs(Distance) < 10000) %>%
                                         distinct(variant) %>%
                                         arrange(variant) %>%
                                         pull(variant), function(variant_label) {

    png_file_name <- file.path(".", sprintf(
      "out.intensity_variants_%s_%s.png", variant_label, format(Sys.Date(), "%Y%m%d")))
    print(png_file_name)

    variant_data_frame_missingness <- combined_data_frame %>% filter(
      variant == variant_label)

    variant_data_frame <- variant_data_frame_missingness %>%
      filter(!is.na(X), !is.na(Y))

    genotype_labels <- get_genotype_labels(variant_data_frame)

    x_limits <- get_axis_limits(variant_data_frame, x_variable = "X", y_variable = "Y")
    y_limits <- x_limits

    annotation_layer <- variant_annotation_layer(
      variant_data_frame_missingness, genotype_labels, missingness_variable = "X",
      x_limits = x_limits, y_limits = y_limits)

    variant_data_frame <- variant_data_frame %>%
      mutate(genotype = recode_factor(genotype, !!!genotype_labels, .ordered = TRUE))

    plot_orig <- plot_cnv_variant_intensities(
      variant_data_frame, x_variable = "X", y_variable = "Y") +
      geom_point(
        aes(colour = CNV_intron_2), alpha = 0.36, size = 0.5) +
      scale_colour_viridis_c(
        option = "plasma", name = "CNV,\nintron 2", na.value = "grey50") +
      theme(legend.position = "None") +
      annotation_layer
    plot_corrected <- plot_cnv_variant_intensities(
      variant_data_frame, x_variable = "X_corrected", y_variable = "Y_corrected") +
      geom_point(
        aes(colour = CNV_intron_2), alpha = 0.36, size = 0.5) +
      scale_colour_viridis_c(
        option = "plasma", name = "CNV,\nintron 2", na.value = "grey50")

    png(png_file_name, width = 1080 + 320, height = 540, res = 150)

    grid.arrange(plot_orig, plot_corrected, ncol = 2)

    par(xpd = NA)

    dev.off()

    return(png_file_name)
  })

  out_file_name <- file.path(".", sprintf(
    "out.intensity_variants_combined_%s", format(Sys.Date(), "%Y%m%d")))

  images_in <- function(png_file_name) {
    return(png_file_name %>% image_read())
  }

  out <- purrr::map(png_file_list, images_in)
  # image_write(image_append(purrr::lift_dl(c)(out), stack = T),
  #             path = paste(out_file_name, "png", sep = "."), format = "png")

  image_write(image_join(out),
              path = paste(out_file_name, "pdf", sep = "."), format = "pdf")

  intensity_variant <- "TOP-rs1058172"
  colouring_variant <- "rs56400210"

  png_file_list <- base::lapply(levels(combined_data_frame$variant), function(colouring_variant) {
    png_file_name <- file.path(".", sprintf(
      "out.intensity_variants_%s_col_%s_%s.png",
      intensity_variant, colouring_variant, format(Sys.Date(), "%Y%m%d")))

    plot_int <- plot_cnv_variant_intensities(
      combined_data_frame, x_variable = "X_corrected", y_variable = "Y_corrected",
      var = colouring_variant)

    plot_geno <- plot_cnv_variant_intensities(
      combined_data_frame, x_variable = "X_corrected", y_variable = "Y_corrected",
      var = intensity_variant)

    plot_combined <- plot_cnv_variant_intensities_with_different_genotype(
      combined_data_frame, x_variable = "X_corrected", y_variable = "Y_corrected",
      intensity_variant = intensity_variant, colouring_variant = colouring_variant)

    png(png_file_name, width = 540 * 3, height = 540, res = 150)

    grid.arrange(plot_int, plot_geno, plot_combined, ncol = 3)

    par(xpd = NA)

    dev.off()

    return(png_file_name)
  })

  out_file_name <- file.path(".", sprintf(
    "out.intensity_variants_combined_%s_%s", intensity_variant, format(Sys.Date(), "%Y%m%d")))

  images_in <- function(png_file_name) {
    return(png_file_name %>% image_read())
  }

  out <- purrr::map(png_file_list, images_in)
  # image_write(image_append(purrr::lift_dl(c)(out), stack = T),
  #             path = paste(out_file_name, "png", sep = "."), format = "png")

  image_write(image_join(out),
              path = paste(out_file_name, "pdf", sep = "."), format = "pdf")

  comparison_data_frame_mom <- combined_data_frame_mom_2 %>%
    inner_join(
      combined_data_frame_mom,
      by = c("variant", "Sample_ID", "SNP", "X",
             "Y", "B_Allele_Freq", "Log_R_Ratio",
             "R", "missingness", "A", "B", "Chrom", "Start", "End", "Distance", "genotype"),
      suffix = c("_insample", "_large")) %>%
    left_join(cnv_phenotypes_mom, by = "Sample_ID")

  png_file_list <- base::lapply(comparison_data_frame_mom %>%
                                  filter(abs(Distance) < 10000) %>%
                                  distinct(variant) %>%
                                  arrange(variant) %>%
                                  pull(variant), function(variant_label) {

    png_file_name <- file.path(".", sprintf(
      "out.intensity_variants_in_or_out_of_sample_%s_%s.png", variant_label, format(Sys.Date(), "%Y%m%d")))
    print(png_file_name)

    variant_data_frame_missingness <- comparison_data_frame_mom %>% filter(
      variant == variant_label)

    variant_data_frame <- variant_data_frame_missingness %>%
      filter(!is.na(X), !is.na(Y))

    genotype_labels <- get_genotype_labels(variant_data_frame)

    x_limits <- get_axis_limits(variant_data_frame, x_variable = "X", y_variable = "Y")
    y_limits <- x_limits

    annotation_layer <- variant_annotation_layer(
      variant_data_frame_missingness, genotype_labels, missingness_variable = "X",
      x_limits = x_limits, y_limits = y_limits)

    variant_data_frame <- variant_data_frame %>%
      mutate(genotype = recode_factor(genotype, !!!genotype_labels, .ordered = TRUE))

    plot_orig <- plot_cnv_variant_intensities(
      variant_data_frame, x_variable = "X", y_variable = "Y") +
      geom_point(
        aes(colour = CNV_intron_2), alpha = 0.36, size = 0.5) +
      scale_colour_viridis_c(
        option = "plasma", name = "CNV,\nintron 2", na.value = "grey50") +
      theme(legend.position = "None") +
      annotation_layer
    plot_corrected_1 <- plot_cnv_variant_intensities(
      variant_data_frame, x_variable = "X_corrected_large", y_variable = "Y_corrected_large") +
      geom_point(
        aes(colour = CNV_intron_2), alpha = 0.36, size = 0.5) +
      scale_colour_viridis_c(
        option = "plasma", name = "CNV,\nintron 2", na.value = "grey50") +
      theme(legend.position = "None") +
      annotation_layer
    plot_corrected_2 <- plot_cnv_variant_intensities(
      variant_data_frame, x_variable = "X_corrected_insample", y_variable = "Y_corrected_insample") +
      geom_point(
        aes(colour = CNV_intron_2), alpha = 0.36, size = 0.5) +
      scale_colour_viridis_c(
        option = "plasma", name = "CNV,\nintron 2", na.value = "grey50")

    png(png_file_name, width = 540*3 + 320, height = 540, res = 150)

    grid.arrange(plot_orig, plot_corrected_1, plot_corrected_2, ncol = 3)

    par(xpd = NA)

    dev.off()

    return(png_file_name)
  })

  out_file_name <- file.path(".", sprintf(
    "out.intensity_variants_combined_in_or_out_of_sample_%s", format(Sys.Date(), "%Y%m%d")))

  images_in <- function(png_file_name) {
    return(png_file_name %>% image_read())
  }

  out <- purrr::map(png_file_list, images_in)
  # image_write(image_append(purrr::lift_dl(c)(out), stack = T),
  #             path = paste(out_file_name, "png", sep = "."), format = "png")

  image_write(image_join(out),
              path = paste(out_file_name, "pdf", sep = "."), format = "pdf")

}

# Test if main should be run
if (sys.nframe() == 0 && !interactive()) {
  main()
}