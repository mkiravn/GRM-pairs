#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# -------------------------
# Config / inputs
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript binplot.R <crossprod_file> <grm_file> [title]")
}

crossprod_file <- args[[1]]
grmm_file       <- args[[2]]
plot_title      <- ifelse(length(args) >= 3, args[[3]], "GRM Bin Plot")

y_col <- "yy"

# -------------------------
# Load + join data
# -------------------------
yy <- readr::read_delim(crossprod_file, delim = "\t", show_col_types = FALSE, col_names = TRUE)
colnames(yy) <- c("IID2", "IID1", "yy")

yy <- yy %>%
  filter(as.numeric(IID1) > 0, as.numeric(IID2) > 0) %>%
  mutate(IID1 = as.character(IID1), IID2 = as.character(IID2)) %>%
  filter(!is.na(yy))

x <- readr::read_delim(
  grmm_file,
  delim = "\t",
  col_names = TRUE,
  col_types = readr::cols(.default = readr::col_guess()),
  show_col_types = FALSE
)

if (ncol(x) < 4) stop("GRM file must have at least 4 columns: IID1 IID2 nsnps GRM")

colnames(x)[1:4] <- c("IID1", "IID2", "N_SNPs", "GRM")
x <- x %>%
  mutate(IID1 = as.character(IID1), IID2 = as.character(IID2), GRM = as.numeric(GRM))

xyy <- right_join(x, yy, by = c("IID1", "IID2")) %>%
  filter(!is.na(GRM), !is.na(yy))

# -------------------------
# Binning config
# -------------------------
breaks <- sort(unique(c(seq(-0.02, 0.02, by = 0.001), seq(0.02, 0.2, by = 0.005), seq(0.2, 0.7, by = 0.01))))

# Fortran binning scheme from phenoCovBins
fortran_bins <- sort(unique(c(-0.3, seq(-0.02, 0.02, by = 0.001), 0.021, 0.023, 0.026, 0.03, 0.04, 0.05, 0.065, 0.1, 0.13, 0.2, 0.35, 0.72)))

xyy <- xyy %>%
  mutate(
    GRM_bin = cut(GRM, breaks = breaks, include.lowest = TRUE, right = FALSE, labels = FALSE),
    bin_center = (breaks[GRM_bin] + breaks[GRM_bin + 1]) / 2,
    GRM_bin_fortran = cut(GRM, breaks = fortran_bins, include.lowest = TRUE, right = FALSE, labels = FALSE),
    bin_center_fortran = (fortran_bins[GRM_bin_fortran] + fortran_bins[GRM_bin_fortran + 1]) / 2
  ) 

bin_summary <- xyy %>%
  filter(!is.na(GRM_bin)) %>%
  group_by(GRM_bin, bin_center) %>%
  summarise(
    n_pairs = n(),
    GRM_mean = mean(GRM, na.rm = TRUE),
    YY_mean = mean(yy, na.rm = TRUE),
    YY_median = median(yy, na.rm = TRUE),
    YY_sd = sd(yy, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(n_pairs>100)

# Fortran bin summary
fortran_bin_summary <- xyy %>%
  filter(!is.na(GRM_bin_fortran)) %>%
  group_by(GRM_bin_fortran, bin_center_fortran) %>%
  summarise(
    n_pairs = n(),
    GRM_mean = mean(GRM, na.rm = TRUE),
    YY_mean = mean(yy, na.rm = TRUE),
    YY_median = median(yy, na.rm = TRUE),
    YY_sd = sd(yy, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(n_pairs>100)

# -------------------------
# Plot 1: points per bin
# -------------------------
p1 <- ggplot() +
  geom_point(data = bin_summary, aes(x = bin_center, y = YY_mean), alpha = 1, size = 1, color = "black") +
  geom_point(data = fortran_bin_summary, aes(x = bin_center_fortran, y = YY_mean), alpha = 1, size = 1, color = "violetred", shape = 4) +
  theme_bw() +
  labs(x = "Genetic Covariance", y = "Mean yy", title = paste0(plot_title, " (bin mean points)")) +
  theme(plot.title = element_text(hjust = 0.5))

# -------------------------
# Plot 2: violin distribution by bin
# -------------------------

summary_df <- xyy %>%
  filter(!is.na(GRM_bin)) %>%
  group_by(bin_center) %>%
  summarise(n=n(),
    q25 = quantile(yy, 0.25, na.rm = TRUE),
    q50 = quantile(yy, 0.50, na.rm = TRUE),
    q75 = quantile(yy, 0.75, na.rm = TRUE),
    mean = mean(yy, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(n>100)

fortran_summary_df <- xyy %>%
  filter(!is.na(GRM_bin_fortran)) %>%
  group_by(bin_center_fortran) %>%
  summarise(n=n(),
    q25 = quantile(yy, 0.25, na.rm = TRUE),
    q50 = quantile(yy, 0.50, na.rm = TRUE),
    q75 = quantile(yy, 0.75, na.rm = TRUE),
    mean = mean(yy, na.rm = TRUE),
    .groups = "drop"
  ) %>% filter(n>100)

p2 <- ggplot() +
  geom_linerange(data = summary_df, aes(x = bin_center, ymin = q25, ymax = q75), color = "lightblue", linewidth = 0.6) +
  geom_point(data = summary_df, aes(x = bin_center, y = q50), color = "steelblue", size = 1.5) +   # median
  geom_point(data = summary_df, aes(x = bin_center, y = mean), color = "black", size = 1.2, shape = 4) +  # mean (cross)
  geom_linerange(data = fortran_summary_df, aes(x = bin_center_fortran, ymin = q25, ymax = q75), color = "pink", linewidth = 0.6) +
  geom_point(data = fortran_summary_df, aes(x = bin_center_fortran + 0.0002, y = q50), color = "violetred", size = 1.5, shape = 17) +   # median
  geom_point(data = fortran_summary_df, aes(x = bin_center_fortran+ 0.0002, y = mean), color = "violetred1", size = 1.2, shape = 4) +  # mean (cross)
  theme_bw() +
  facet_wrap(~bin_center > 0.1, scales = "free", ncol=1) +
  labs(
    x = "Genetic Covariance",
    y = "yy",
    title = paste0(plot_title, " (bin distributions)")
  ) +
  theme(plot.title = element_text(hjust = 0.5))
# -------------------------
# Plot 3: zoomed-in point plot
# -------------------------
p3 <- p1 +
  coord_cartesian(xlim = c(-0.02, 0.05)) +
  labs(title = paste0(plot_title, " (zoom -0.02 to 0.05)"))

# -------------------------
# Write output
# -------------------------
out_base <- sub("\\.[^.]*$", "", basename(crossprod_file))

plot1_file <- paste0(out_base, "__bin_mean_points.png")
plot2_file <- paste0(out_base, "__bin_dist.png")
plot3_file <- paste0(out_base, "__bin_zoom.png")

ggsave(plot1_file, p1, width = 6, height = 4)
message(sprintf("[INFO] Wrote: %s", plot1_file))

ggsave(plot2_file, p2, width = 6, height = 6)
message(sprintf("[INFO] Wrote: %s", plot2_file))

ggsave(plot3_file, p3, width = 6, height = 4)
message(sprintf("[INFO] Wrote: %s", plot3_file))