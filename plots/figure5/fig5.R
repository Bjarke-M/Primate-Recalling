#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)

# ---- paths ----------------------------------------------------------------
BASE_DIR    <- "~alnRegionHet/"
SAMPLE_META <- "~SupTable_Sample_Stats_wGT_QC.tsv"
OUT_PLOT    <- "~fig_5.pdf"

# ---- load metadata --------------------------------------------------------
supptab <- read.table(SAMPLE_META, header = TRUE, sep = "\t")
supptab$remove_as_relative <- as.character(supptab$remove_as_relative)
supptab$remove_as_relative[is.na(supptab$remove_as_relative)] <- "FALSE"
meta <- subset(supptab, finalQC != "fail" & cov_chrA >= 20 & mis_chrA < 0.01 & remove_as_relative == "FALSE")

# ---- discover all species dirs --------------------------------------------
all_species <- list.dirs(BASE_DIR, full.names = FALSE, recursive = FALSE)
all_species <- all_species[all_species != ""]  # drop empty

cat("Found", length(all_species), "species directories\n")

# ---- helper: load one species ---------------------------------------------
load_species <- function(sp) {
  stats_dir    <- file.path(BASE_DIR, sp)
  callable_tsv <- file.path(stats_dir, "callable_lengths.tsv")
  
  if (!file.exists(callable_tsv)) {
    message("Skipping ", sp, ": no callable_lengths.tsv")
    return(NULL)
  }
  
  het_files <- list.files(stats_dir, pattern = "_het_sites\\.tsv$", full.names = TRUE)
  if (length(het_files) == 0) {
    message("Skipping ", sp, ": no het_sites files")
    return(NULL)
  }
  
  # Parse filenames
  parse_filename <- function(f) {
    base <- basename(f)
    parts <- regmatches(base, regexec(
      "(.+)_batch(\\d+)_fploidy(\\d+)_mploidy(\\d+)_(unique|nonunique|noalign)_het_sites\\.tsv",
      base, perl = TRUE
    ))[[1]]
    if (length(parts) == 0) return(NULL)
    data.frame(file = f, SPECIES = parts[2], BATCH = as.integer(parts[3]),
               FPLOIDY = as.integer(parts[4]), MPLOIDY = as.integer(parts[5]),
               BED_TYPE = parts[6], stringsAsFactors = FALSE)
  }
  file_meta <- bind_rows(lapply(het_files, parse_filename))
  if (nrow(file_meta) == 0) return(NULL)
  
  # CHR_TYPE
  file_meta$CHR_TYPE <- "A"
  file_meta$CHR_TYPE[file_meta$FPLOIDY == 2 & file_meta$MPLOIDY == 1] <- "X"
  
  # Read het files
  read_het <- function(row) {
    df <- tryCatch(
      read_tsv(row$file, col_names = c("IND_ID", "CHR", "NO_SNPS"),
               col_types = "cci", skip = 1, show_col_types = FALSE),
      error = function(e) NULL
    )
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$SPECIES  <- row$SPECIES;  df$BATCH    <- row$BATCH
    df$FPLOIDY  <- row$FPLOIDY;  df$MPLOIDY  <- row$MPLOIDY
    df$BED_TYPE <- row$BED_TYPE; df$CHR_TYPE <- row$CHR_TYPE
    df
  }
  het_all <- bind_rows(lapply(seq_len(nrow(file_meta)), function(i) read_het(file_meta[i, ])))
  if (is.null(het_all) || nrow(het_all) == 0) return(NULL)
  
  het_summed <- het_all %>%
    group_by(SPECIES, IND_ID, BED_TYPE, CHR_TYPE) %>%
    summarise(NO_SNPS = sum(NO_SNPS), .groups = "drop")
  
  # Callable
  callable <- read_tsv(callable_tsv, show_col_types = FALSE)
  chr_type_key <- file_meta %>% select(SPECIES, BATCH, FPLOIDY, MPLOIDY, BED_TYPE, CHR_TYPE) %>% distinct()
  callable_total <- callable %>%
    left_join(chr_type_key, by = c("SPECIES", "BATCH", "FPLOIDY", "MPLOIDY", "BED_TYPE")) %>%
    filter(!is.na(CHR_TYPE)) %>%
    group_by(SPECIES, BED_TYPE, CHR_TYPE) %>%
    summarise(TOTAL_CALLABLE = sum(N_CALLABLE), .groups = "drop")
  
  het_merged <- het_summed %>%
    left_join(callable_total, by = c("SPECIES", "BED_TYPE", "CHR_TYPE")) %>%
    mutate(HETEROZYGOSITY = NO_SNPS / TOTAL_CALLABLE) %>%
    left_join(meta, by = c("IND_ID" = "ID")) %>%
    filter(!is.na(finalQC))
  
  het_merged
}

# ---- load all species -----------------------------------------------------
all_data <- bind_rows(lapply(all_species, load_species))

if (nrow(all_data) == 0) stop("No data loaded — check paths and file structure")
cat("Loaded", n_distinct(all_data$SPECIES), "species,", n_distinct(all_data$IND_ID), "individuals\n")

# ---- compute per-genus medians for A and X --------------------------------
family_summary <- all_data %>%
  filter(CHR_TYPE %in% c("A", "X")) %>%
  group_by(family, CHR_TYPE, BED_TYPE) %>%
  summarise(
    MEDIAN_HET = median(HETEROZYGOSITY, na.rm = TRUE),
    MEAN_HET   = mean(HETEROZYGOSITY, na.rm = TRUE),
    Q25        = quantile(HETEROZYGOSITY, 0.25, na.rm = TRUE),
    Q75        = quantile(HETEROZYGOSITY, 0.75, na.rm = TRUE),
    N_IND      = n_distinct(IND_ID),
    N_SPECIES  = n_distinct(SPECIES),
    .groups = "drop"
  )

family_summary$family   <- factor(family_summary$family,   levels = rev(c("Hominidae","Hylobatidae","Cercopithecidae","Callitrichidae","Aotidae",
                                                                    "Cebidae","Atelidae","Pitheciidae","Tarsiidae","Lemuridae","Indriidae",
                                                                    "Cheirogaleidae", "Lepilemuridae", "Daubentoniidae","Lorisidae", "Galagidae")))
family_summary$BED_TYPE <- factor(family_summary$BED_TYPE, levels = c("unique", "nonunique", "noalign"))
family_summary$CHR_TYPE <- factor(family_summary$CHR_TYPE, levels = c("A", "X"))

# ---- colours --------------------------------------------------------------
BED_COLORS <- c("unique"    = "#2196F3",
                "nonunique" = "#FF9800",
                "noalign"   = "#F44336")

# ===== PLOT: all individuals together ======================================
all_data$BED_TYPE <- factor(all_data$BED_TYPE, levels = c("unique", "nonunique", "noalign"))

p1 <- ggplot(filter(all_data, CHR_TYPE %in% c("A", "X")),
             aes(x = BED_TYPE, y = HETEROZYGOSITY, colour = BED_TYPE)) +
  geom_violin(fill = NA, linewidth = 0.6, adjust = 1.2) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.7, colour = "grey20") +
  facet_wrap(~ CHR_TYPE, nrow = 1, 
             labeller = labeller(CHR_TYPE = c("A" = "Autosomes", "X" = "chrX"))) +
  scale_colour_manual(values = BED_COLORS, guide = "none") +
  scale_x_discrete(labels = c("unique"    = "Unique",
                              "nonunique" = "Non-unique",
                              "noalign"   = "No alignment")) +
  scale_y_continuous(labels = scales::scientific) +
  labs(
    tag      = "A",
    x        = "Alignment class",
    y        = "Heterozygosity (het sites / callable bp)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title       = element_text(size = 11, face = "bold"),
    plot.subtitle    = element_text(size = 8, colour = "grey40"),
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_blank(),
    axis.text.x      = element_text(size = 9, angle = 30, hjust = 1)
  )

# ===== PLOT: median het per genus per bed type, faceted by CHR_TYPE ========
p2 <- ggplot(family_summary,
             aes(x = MEAN_HET, y = family, colour = BED_TYPE)) +
  geom_linerange(aes(xmin = Q25, xmax = Q75),
                 position = position_dodge(width = 0.6), linewidth = 0.5, alpha = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 1.8) +
  facet_wrap(~ CHR_TYPE, nrow = 1, 
             labeller = labeller(CHR_TYPE = c("A" = "Autosomes", "X" = "chrX"))) +
  scale_colour_manual(values = BED_COLORS,
                      labels = c("unique" = "Unique", "nonunique" = "Non-unique", "noalign" = "No alignment"),
                      name = "Alignment class") +
  scale_x_continuous(labels = scales::scientific) +
  labs(
    tag      = "B",
    x        = "Heterozygosity (het sites / callable bp)",
    y        = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.y        = element_text(size = 7, face = "italic"),
    legend.position    = "top",
    plot.title         = element_text(size = 11, face = "bold"),
    plot.subtitle      = element_text(size = 8, colour = "grey40"),
    strip.text         = element_text(size = 10, face = "bold"),
    strip.background   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3)
  )

# ===== SAVE ================================================================
pdf(OUT_PLOT, width = 12, height = 5)
p1 + p2 + plot_layout(widths = c(1, 2)) + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 12))
dev.off()
cat("Plots written to:", OUT_PLOT, "\n")