#!/usr/bin/env Rscript

setwd("~/Library/CloudStorage/OneDrive-UWA/Research/Projects/Muscle_Specific_Exons")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

# -------------------------- CONFIG --------------------------
in_tsv  <- "./pext_with_noncanonical_exon_annot.reordered.v49.bed"      # input file
out_dir <- "pext_zscores_out"    # output folder
z_thr   <- 3                     # z-score threshold
mode    <- "high"                # "high" (z >= thr) or "abs" (|z| >= thr)

# How to handle SD=0 within an exon (all tissues same pext):
# "na"  => z = NA (no tissue flagged)
# "zero"=> z = 0
sd0_mode <- "na"

# Optional: ignore non-pext columns that happen to start with "pext_" but aren't tissues
# e.g. pext_exp_prop_mean / pext_pext_union etc. Add patterns you want to DROP here:
drop_pext_cols <- c("^pext_exp_prop_mean$", "^pext_pext_union$", "^pext_pext_union$")

# -------------------------- HELPERS --------------------------
calc_z <- function(x, sd0_mode = c("na","zero")) {
  sd0_mode <- match.arg(sd0_mode)
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    if (sd0_mode == "na") return(rep(NA_real_, length(x)))
    return(rep(0, length(x)))
  }
  (x - m) / s
}

flag_from_z <- function(z, thr = 2, mode = c("high","abs")) {
  mode <- match.arg(mode)
  if (mode == "high") return(z >= thr)
  abs(z) >= thr
}

# -------------------------- LOAD --------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- read_tsv(in_tsv, show_col_types = FALSE)

# Identify pext tissue columns:
pext_cols <- names(df) %>%
  keep(~ str_starts(.x, "pext_"))

# Drop any "pext_" columns you *don't* want treated as tissues (optional):
if (length(drop_pext_cols) > 0) {
  drop_regex <- paste(drop_pext_cols, collapse = "|")
  pext_cols <- pext_cols[!str_detect(pext_cols, drop_regex)]
}

if (length(pext_cols) == 0) {
  stop("No pext_* columns found after filtering. Check your input / drop_pext_cols.")
}

# Define exon identity columns (adjust if yours differ)
exon_id_cols <- c(
  "chrom","start","end",
  "exon_chr","exon_start","exon_end","exon_gene",
  "overlap_type"
)
missing_ids <- setdiff(exon_id_cols, names(df))
if (length(missing_ids) > 0) {
  stop("Missing expected exon id columns: ", paste(missing_ids, collapse = ", "))
}

# -------------------------- LONG + ZSCORES --------------------------
long <- df %>%
  select(all_of(exon_id_cols), all_of(pext_cols)) %>%
  pivot_longer(
    cols = all_of(pext_cols),
    names_to = "tissue",
    values_to = "pext"
  ) %>%
  mutate(
    tissue = str_remove(tissue, "^pext_"),
    pext   = as.numeric(pext)
  ) %>%
  group_by(across(all_of(exon_id_cols))) %>%
  mutate(
    exon_mean = mean(pext, na.rm = TRUE),
    exon_sd   = sd(pext, na.rm = TRUE),
    z         = calc_z(pext, sd0_mode = sd0_mode),
    flagged   = flag_from_z(z, thr = z_thr, mode = mode)
  ) %>%
  ungroup()

# -------------------------- PER-EXON SUMMARY --------------------------
# Tissues meeting the z threshold per exon + a few distribution metrics
exon_summary <- long %>%
  group_by(across(all_of(exon_id_cols))) %>%
  summarise(
    n_tissues = sum(!is.na(pext)),
    mean_pext = first(exon_mean),
    sd_pext   = first(exon_sd),
    min_pext  = min(pext, na.rm = TRUE),
    max_pext  = max(pext, na.rm = TRUE),
    # robust summaries can help when lots of zeros
    median_pext = median(pext, na.rm = TRUE),
    iqr_pext    = IQR(pext, na.rm = TRUE),
    n_flagged   = sum(flagged, na.rm = TRUE),
    flagged_tissues = paste(tissue[which(flagged %in% TRUE)], collapse = ","),
    # for convenience: top tissues by z-score
    top3_tissues_by_z = paste(tissue[order(z, decreasing = TRUE)][1:3], collapse = ","),
    top3_z            = paste(round(sort(z, decreasing = TRUE)[1:3], 3), collapse = ","),
    .groups = "drop"
  )


# Filter for specific tissue outliers


# ---- Filter for specific tissue outliers (Muscle) ----

safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

collapse_same_except_id <- function(x, exon_id_cols) {
  x %>%
    dplyr::distinct(dplyr::across(-dplyr::all_of(exon_id_cols)), .keep_all = TRUE)
}

get_tissue_exon_outliers <- function(
    exon_summary,
    long,
    tissue_pattern,
    exon_id_cols,
    min_pext = 0,
    min_z = NULL,
    match = c("glob", "regex"),
    collapse_dupes = TRUE
) {
  match <- match.arg(match)
  
  # ---- helper: safe max (avoids -Inf warnings) ----
  safe_max <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_real_ else max(x)
  }
  
  
  # ---- optionally collapse duplicates ----
  if (collapse_dupes) {
    long <- collapse_same_except_id(long, exon_id_cols)
  }
  
  # ---- tissue pattern matching ----
  tissue_regex <- if (match == "glob") {
    paste0("^", gsub("\\*", ".*", tissue_pattern), "$")
  } else {
    tissue_pattern
  }
  
  tissue_tbl <- long %>%
    dplyr::filter(grepl(tissue_regex, tissue)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(exon_id_cols))) %>%
    dplyr::summarise(
      pext_tissue     = safe_max(pext),
      tissue_z        = safe_max(z),
      tissue_flagged  = any(flagged %in% TRUE),
      tissues_matched = paste(unique(tissue), collapse = ","),
      .groups = "drop"
    )
  
  exon_summary %>%
    dplyr::inner_join(
      tissue_tbl %>%
        dplyr::filter(
          tissue_flagged,
          pext_tissue > min_pext,
          if (!is.null(min_z)) (!is.na(tissue_z) & tissue_z >= min_z) else TRUE
        ),
      by = exon_id_cols
    )
}

brain_exon_outliers <- get_tissue_exon_outliers(
  exon_summary = exon_summary,
  long         = long,
  tissue       = "Brain_*",
  exon_id_cols = exon_id_cols,
  min_pext     = 0.1,
  min_z        = z_thr,
  collapse_dupes = FALSE
)

muscle_exon_outliers <- get_tissue_exon_outliers(
  exon_summary = exon_summary,
  long         = long,
  tissue       = "Muscle_Skeletal",
  exon_id_cols = exon_id_cols,
  min_pext     = 0.1,
  min_z        = z_thr,
  collapse_dupes = FALSE
)

nerve_exon_outliers <- get_tissue_exon_outliers(
  exon_summary = exon_summary,
  long         = long,
  tissue       = "Nerve_Tibial",
  exon_id_cols = exon_id_cols,
  min_pext     = 0.1,
  min_z        = z_thr,
  collapse_dupes = FALSE
)

get_any_tissue_exon_outliers <- function(
    exon_summary,
    long,
    exon_id_cols,
    min_pext = 0,
    min_z = NULL,
    collapse_dupes = TRUE
) {
  if (collapse_dupes) {
    long <- collapse_same_except_id(long, exon_id_cols)
  }
  
  hits <- long %>%
    dplyr::filter(
      flagged %in% TRUE,
      pext > min_pext,
      if (!is.null(min_z)) z >= min_z else TRUE
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(exon_id_cols))) %>%
    dplyr::summarise(
      n_flagged = dplyr::n(),
      flagged_tissues = paste(tissue, collapse = ","),
      top3_tissues_by_z = paste(tissue[order(z, decreasing = TRUE)][1:min(3, dplyr::n())], collapse = ","),
      top3_z = paste(round(sort(z, decreasing = TRUE)[1:min(3, dplyr::n())], 3), collapse = ","),
      .groups = "drop"
    )
  
  exon_summary %>%
    dplyr::inner_join(hits, by = exon_id_cols)
}

any_outliers <- get_any_tissue_exon_outliers(
  exon_summary = exon_summary,
  long = long,
  exon_id_cols = exon_id_cols,
  min_pext = 0.1,
  min_z = z_thr,
  collapse_dupes = FALSE
)

# -------------------------- SAVE OUTPUTS --------------------------
write_tsv(long, file.path(out_dir, "pext_long_with_zscores_v49.tsv"))
write_tsv(exon_summary, file.path(out_dir, "exon_pext_zscore_summary_v49.tsv"))
write_tsv(muscle_exon_outliers, file.path(out_dir, "muscle_exon_outliers_v49.tsv"))
write_tsv(nerve_exon_outliers, file.path(out_dir, "nerve_exon_outliers_v49.tsv"))
write_tsv(brain_exon_outliers, file.path(out_dir, "brain_exon_outliers_v49.tsv"))
write_tsv(any_outliers, file.path(out_dir, "any_exon_outliers_v49.tsv"))

# Also save only the (exon,tissue) rows that meet the threshold
hits <- long %>% filter(flagged %in% TRUE)
write_tsv(hits, file.path(out_dir, sprintf("exon_tissue_hits_z%s_%s_v49.tsv", z_thr, mode)))
