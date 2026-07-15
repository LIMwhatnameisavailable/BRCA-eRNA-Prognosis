# analysis_NET1e_expression.R
# NET1e (ENSR00000023843) expression in High vs Low risk groups
# Coordinates: hg19 chr10:17,440,000-17,460,000 (16 loci)
# Reference: Zhang et al., Nat Commun 2019

library(data.table)
library(ggplot2)
library(ggpubr)

# Step 1: Load eRNA matrix
data_net1e <- fread(
  "Data_Source/TCGA_RPKM_eRNA_300k_peaks_in_Super_enhancer_BRCA.csv",
  data.table = TRUE
)

# Step 2: Standardize column names to TCGA barcode format
orig_cols <- colnames(data_net1e)
new_cols <- ifelse(
  grepl("_tumor",  orig_cols), gsub("_tumor",  "-01A", orig_cols),
  ifelse(
    grepl("_normal", orig_cols), gsub("_normal", "-11A", orig_cols),
    orig_cols
  )
)
colnames(data_net1e) <- new_cols

# Step 3: Retain tumor sample columns (14th character == "0" indicates tumor)
all_sample_cols  <- setdiff(colnames(data_net1e), "V1")
tumor_cols_net1e <- all_sample_cols[substr(all_sample_cols, 14, 14) == "0"]

# Step 4: Extract NET1e core region loci (hg19 chr10:17,440,000-17,460,000)
all_row_ids <- data_net1e$V1

net1e_ids <- all_row_ids[
  grepl("^chr10:", all_row_ids) & {
    pos <- suppressWarnings(as.numeric(sub("chr10:", "", all_row_ids)))
    !is.na(pos) & pos >= 17440000 & pos <= 17460000
  }
]

cat("NET1e loci found:", length(net1e_ids), "(expected 16)\n")
stopifnot("No NET1e loci found, check coordinates or input file" = length(net1e_ids) > 0)

# Step 5: Compute per-sample NET1e score as mean across core loci
net1e_rows <- data_net1e[V1 %in% net1e_ids]
net1e_mat  <- as.matrix(net1e_rows[, ..tumor_cols_net1e])
rownames(net1e_mat) <- net1e_rows$V1

net1e_score <- colMeans(net1e_mat, na.rm = TRUE)

# Step 6: Build sample data frame with 12-character patient IDs
net1e_df <- data.frame(
  sample_id  = substr(tumor_cols_net1e, 1, 12),
  NET1e_RPKM = net1e_score,
  stringsAsFactors = FALSE
)

# Average across multiple tumor samples per patient
net1e_df <- aggregate(NET1e_RPKM ~ sample_id, data = net1e_df, FUN = mean)

# Step 7: Extract risk group labels from combined_data
stopifnot(
  "risk_group column not found in combined_data, run main pipeline first" =
    "risk_group" %in% colnames(combined_data)
)

risk_df <- data.frame(
  sample_id  = rownames(combined_data),
  risk_group = combined_data$risk_group,
  stringsAsFactors = FALSE
)

# Step 8: Merge NET1e expression with risk groups
merged_net1e <- merge(net1e_df, risk_df, by = "sample_id")

if (nrow(merged_net1e) < 100) {
  cat("ID alignment failed.\n")
  cat("net1e_df IDs:", head(net1e_df$sample_id, 3), "\n")
  cat("risk_df  IDs:", head(risk_df$sample_id,  3), "\n")
  stop("Merged sample count too low, check ID format")
}

cat("Merged samples:", nrow(merged_net1e),
    "| High:", sum(merged_net1e$risk_group == "High"),
    "| Low:",  sum(merged_net1e$risk_group == "Low"), "\n")

# Step 9: Descriptive statistics
med_high <- median(merged_net1e$NET1e_RPKM[merged_net1e$risk_group == "High"])
med_low  <- median(merged_net1e$NET1e_RPKM[merged_net1e$risk_group == "Low"])

print(tapply(
  merged_net1e$NET1e_RPKM,
  merged_net1e$risk_group,
  function(x) round(c(median = median(x), mean = mean(x), n = length(x)), 4)
))

# Step 10: Wilcoxon rank-sum test
wilcox_result <- wilcox.test(NET1e_RPKM ~ risk_group,
                             data = merged_net1e, exact = FALSE)
p_val <- wilcox_result$p.value

p_label <- ifelse(p_val < 0.001, "p < 0.001",
           ifelse(p_val < 0.01,  "p < 0.01",
                  paste0("p = ", round(p_val, 3))))

cat("Wilcoxon p =", format.pval(p_val, digits = 3), "\n")

# Step 11: Plot
merged_net1e$risk_group <- factor(merged_net1e$risk_group,
                                  levels = c("Low", "High"))

p_net1e <- ggplot(
  merged_net1e,
  aes(x = risk_group, y = NET1e_RPKM, fill = risk_group)
) +
  geom_boxplot(
    outlier.shape = 21,
    outlier.size  = 1.5,
    outlier.alpha = 0.4,
    width         = 0.5,
    linewidth     = 0.7
  ) +
  geom_jitter(
    aes(color = risk_group),
    width = 0.15, size = 0.7, alpha = 0.25
  ) +
  scale_fill_manual(values  = c("High" = "#D62728", "Low" = "#1F77B4")) +
  scale_color_manual(values = c("High" = "#D62728", "Low" = "#1F77B4")) +
  stat_compare_means(
    method      = "wilcox.test",
    comparisons = list(c("Low", "High")),
    label       = "p.signif",
    size        = 6,
    tip.length  = 0.02
  ) +
  labs(
    title    = "NET1e Expression by Risk Group",
    subtitle = paste0(
      "NET1e (ENSR00000023843) | hg19 chr10:17,440,000-17,460,000 | 16 loci\n",
      "Zhang et al., Nat Commun 2019 | ",
      "n = ", nrow(merged_net1e), " | ", p_label
    ),
    x     = "Risk Group",
    y     = "NET1e Mean RPKM (16 loci, raw)",
    fill  = "Risk Group",
    color = "Risk Group"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9, color = "grey40"),
    panel.grid.major.x = element_blank(),
    axis.text          = element_text(color = "black")
  )

print(p_net1e)

# Step 12: Save outputs
ggsave("Results/NET1e_HighvsLow_RiskGroup.pdf",
       plot = p_net1e, width = 4, height = 5, dpi = 300)
ggsave("Results/NET1e_HighvsLow_RiskGroup.png",
       plot = p_net1e, width = 4, height = 5, dpi = 300)

cat("Analysis complete.\n")
cat("Loci included    :", length(net1e_ids), "\n")
cat("Samples analyzed :", nrow(merged_net1e), "\n")
cat("Median RPKM High :", round(med_high, 4), "\n")
cat("Median RPKM Low  :", round(med_low,  4), "\n")
cat("Wilcoxon p       :", format.pval(p_val, digits = 3), "\n")
