# ==============================================================================
# Script 01: Metabolic Landscape & Differential Analysis

# 1. Setup Environment
packages <- c("pheatmap", "RColorBrewer", "viridis", "MASS", "ggplot2", "eoffice", "scatterplot3d")
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# Set directories
data_dir <- "data"
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ==================== Part A: Global Heatmap ====================
# Load Data
sample_data <- read.table(file.path(data_dir, "discovery_data.txt"), 
                          header = TRUE, sep = "\t", check.names = FALSE)
rownames(sample_data) <- sample_data$SampleID
sample_data_exp <- as.matrix(sample_data[, -(1:3)]) 

# Log2 Transform
sample_data_exp_log <- log2(sample_data_exp + 1)

# Normalization Functions
normalize <- function(x) { (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) }

# Double Z-score Normalization (Row then Column)
sample_data_heatmap <- t(sample_data_exp_log)
data_zscore_row <- t(apply(sample_data_heatmap, 1, normalize))
data_zscore_both <- apply(data_zscore_row, 2, normalize)

# Group Ordering
column_annotation <- sample_data[, "Group2", drop = FALSE]
group_order <- c("HC", "BGD", "EGC", "AGC")
column_annotation$Group2 <- factor(column_annotation$Group2, levels = group_order)
column_annotation <- column_annotation[order(column_annotation$Group2), , drop = FALSE]
data_zscore_ordered <- data_zscore_both[, rownames(column_annotation)]

# Plot Heatmap
annotation_colors <- list(Group2 = c("HC"="#87B60D", "BGD"="#EEC76A", "EGC"="#EF8E1B", "AGC"="#C74543"))
col_fun <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(500)

pdf(file.path(out_dir, "01_Global_Heatmap.pdf"), width = 12, height = 12)
pheatmap(data_zscore_ordered,
         cluster_rows = TRUE, cluster_cols = FALSE,
         color = col_fun, breaks = seq(-2, 2, length.out = 501),
         show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = column_annotation, annotation_colors = annotation_colors,
         main = "Global Metabolic Landscape")
dev.off()

# ==================== Part B: LDA Analysis ====================
# Transpose for LDA 
lda_input <- t(data_zscore_both)
lda_df_input <- data.frame(lda_input, Group2 = sample_data$Group2)

# Run LDA
lda_model <- lda(Group2 ~ ., data = lda_df_input)
lda_coords <- data.frame(predict(lda_model)$x, Group = sample_data$Group2)
lda_coords$Group <- factor(lda_coords$Group, levels = group_order)

# Plot LDA
p_lda <- ggplot(lda_coords, aes(x = LD1, y = LD2, color = Group)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(level = 0.95, size = 1) +
  scale_color_manual(values = annotation_colors$Group2) +
  theme_minimal() +
  labs(title = "LDA Projection")

ggsave(file.path(out_dir, "01_LDA_Plot.pdf"), p_lda, width = 8, height = 6)

# ==================== Part C: Differential Analysis Visualization ====================
# Load 3D Input Data 
data_3d <- read.table(file.path(data_dir, "3D_input.txt"), header = TRUE, sep = "\t")

# Define Colors based on criteria
data_3d$color <- ifelse(abs(data_3d$LogFC) > 0.263 & data_3d$FDR < 0.05 & data_3d$VIP > 1, "#A03434",
                        ifelse(abs(data_3d$LogFC) > 0.263, "#C1D88F",
                               ifelse(data_3d$FDR < 0.05, "#F5E5B4",
                                      ifelse(data_3d$VIP > 1, "#CCEDFC", "#CBCDCC"))))

# Plot 3D Scatter
pdf(file.path(out_dir, "01_3D_Differential_Plot.pdf"), width = 8, height = 8)
scatterplot3d(x = data_3d$LogFC, y = -log10(data_3d$FDR), z = data_3d$VIP,
              color = data_3d$color, pch = 20, cex.symbols = 1.5,
              xlab = "Log2 FC", ylab = "-Log10 FDR", zlab = "VIP Score",
              main = "Differential Metabolites (NGC vs GC)")
dev.off()

message("Script 01 Completed. Check 'output' folder.")