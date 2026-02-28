#!/usr/bin/env Rscript 
## hmm based haplotype reconstruction 
setwd("/home/fetche-lab/Desktop/gn_remote/HS_rats/HS_data_v4/genotypes/new_processing/tests/chr1/results/rqlt2_inputs/transposed/")

#install.packages("qtl2") 
library(qtl2) 

start_time <- Sys.time() 

cross <- read_cross2("hs_rn8_cross.yaml") 

# Avoid integer overflow in calc_genoprob() on large marker sets by
# processing each chromosome in marker chunks.
marker_chunk_size <- 5000L
max_markers_per_chr <- as.integer(Sys.getenv("QTL2_MAX_MARKERS", "0"))

block_results <- list()

for(chr in names(cross$geno)) {
  geno_chr <- cross$geno[[chr]]
  marker_names <- colnames(geno_chr)
  if (!is.na(max_markers_per_chr) && max_markers_per_chr > 0L) {
    marker_names <- marker_names[seq_len(min(length(marker_names), max_markers_per_chr))]
  }
  individual_names <- rownames(geno_chr)
  n_markers <- length(marker_names)

  founder_call <- matrix(
    NA_character_,
    nrow = n_markers,
    ncol = length(individual_names),
    dimnames = list(marker_names, individual_names)
  )

  chunk_starts <- seq.int(1L, n_markers, by = marker_chunk_size)
  for (chunk_start in chunk_starts) {
    chunk_end <- min(chunk_start + marker_chunk_size - 1L, n_markers)
    idx <- chunk_start:chunk_end

    cross_chunk <- cross
    cross_chunk$geno <- cross$geno[chr]
    cross_chunk$gmap <- cross$gmap[chr]
    if (!is.null(cross$founder_geno)) {
      cross_chunk$founder_geno <- cross$founder_geno[chr]
    }

    cross_chunk$geno[[chr]] <- cross_chunk$geno[[chr]][, idx, drop = FALSE]
    cross_chunk$gmap[[chr]] <- cross_chunk$gmap[[chr]][idx]
    if (!is.null(cross_chunk$founder_geno)) {
      cross_chunk$founder_geno[[chr]] <- cross_chunk$founder_geno[[chr]][, idx, drop = FALSE]
    }

    pr_chunk <- calc_genoprob(
      cross_chunk,
      error_prob = 0.002,
      lowmem = TRUE
    )

    apr_chunk <- genoprob_to_alleleprob(pr_chunk)[[chr]]
    founder_names <- dimnames(apr_chunk)[[3]]

    founder_index <- apply(apr_chunk, c(1,2), which.max)
    max_prob <- apply(apr_chunk, c(1,2), max)

    founder_chunk <- matrix(
      founder_names[founder_index],
      nrow = nrow(founder_index),
      ncol = ncol(founder_index),
      dimnames = dimnames(founder_index)
    )

    threshold <- 0.8
    founder_chunk[max_prob < threshold] <- NA
    founder_call[rownames(founder_chunk), colnames(founder_chunk)] <- founder_chunk
  }

  for(ind in colnames(founder_call)) {

    calls <- founder_call[, ind]
    marker_names <- rownames(founder_call)

    # --- NA-robust block detection ---
    change_points <- c(TRUE, calls[-1] != calls[-length(calls)])
    block_id <- cumsum(change_points)

    blocks <- split(seq_along(calls), block_id)

    for(b in blocks) {

      proximal_idx <- b[1]
      distal_idx   <- b[length(b)]

      block_results[[length(block_results)+1]] <- data.frame(
        chr = chr,
        marker = marker_names[proximal_idx],
        individual = ind,
        founder = calls[proximal_idx],
        position_type = "proximal"
      )

      block_results[[length(block_results)+1]] <- data.frame(
        chr = chr,
        marker = marker_names[distal_idx],
        individual = ind,
        founder = calls[distal_idx],
        position_type = "distal"
      )
    }
  }
}

block_df <- do.call(rbind, block_results)

block_matrix <- tidyr::pivot_wider(
  block_df,
  names_from = individual,
  values_from = founder
)

write.csv(
  block_matrix,
  "smoothed_founder_blocks_proximal_distal_50k.csv",
  row.names = FALSE
)

end_time <- Sys.time() 

runtime = end_time - start_time 

print(paste("Total Runtime:", runtime)) 
