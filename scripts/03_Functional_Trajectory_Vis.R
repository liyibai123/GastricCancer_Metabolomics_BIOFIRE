# ==============================================================================
# Script 03: Visualization of Biomarkers & Trajectories
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

data_dir <- "data"
out_dir <- "output"

# ==================== Trajectory Analysis ====================
# Load Expression and Clusters
# Note: Using model_exp.txt (absolute quant) for plotting trends
raw_exp <- read.table(file.path(data_dir, "model_exp.txt"), header = TRUE)
cluster_info <- read.csv(file.path(out_dir, "02_NMF_Clusters.csv")) # Generated in step 2

# Reshape and Normalize
long_data <- raw_exp %>%
  tidyr::pivot_longer(cols = -c(Index, Group1, Group2), names_to = "Metabolite", values_to = "Value") %>%
  group_by(Metabolite) %>%
  mutate(Z_Score = scale(log2(Value + 1))) %>%
  left_join(cluster_info, by = "Metabolite")

# Calculate Mean Trends per Group per Cluster
trend_data <- long_data %>%
  group_by(Cluster, Group2) %>%
  summarise(Mean_Z = mean(Z_Score, na.rm = TRUE), 
            SE = sd(Z_Score, na.rm=TRUE)/sqrt(n()), .groups = "drop")

# Order Groups
trend_data$Group2 <- factor(trend_data$Group2, levels = c("HC", "BGD", "EGC", "AGC"))

# Plot Trajectories
p_traj <- ggplot(trend_data, aes(x = Group2, y = Mean_Z, group = Cluster, color = as.factor(Cluster))) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = Mean_Z - SE, ymax = Mean_Z + SE, fill = as.factor(Cluster)), alpha = 0.1, color=NA) +
  theme_bw() +
  labs(title = "Metabolic Module Trajectories", y = "Standardized Abundance (Z-score)")

ggsave(file.path(out_dir, "03_Module_Trajectories.pdf"), p_traj, width = 8, height = 6)

# ==================== Individual Boxplots ====================
# Load Boxplot Data
box_input <- read.table(file.path(data_dir, "boxplot_input.txt"), header = TRUE, row.names = 1)
box_group <- read.table(file.path(data_dir, "boxplot_group.txt"), header = TRUE)

# Filter for the 12 selected markers
# (Assuming selected markers are known or loaded from step 2 output)
selected_markers <- read.table(file.path(out_dir, "02_Final_12_Panel.txt"), header = TRUE)$Metabolite
plot_data <- t(box_input[selected_markers, ]) %>% as.data.frame()
plot_data$Group <- box_group$Group[match(rownames(plot_data), box_group$Sample)]
plot_data$Group <- factor(plot_data$Group, levels = c("HC", "BGD", "EGC", "AGC"))

# Loop Plot
pdf(file.path(out_dir, "03_Biomarker_Boxplots.pdf"), width = 5, height = 5)
for(met in selected_markers) {
  p <- ggplot(plot_data, aes(x = Group, y = .data[[met]], fill = Group)) +
    geom_boxplot(alpha = 0.8) +
    geom_jitter(width = 0.2, alpha = 0.3) +
    theme_classic() +
    labs(title = met, y = "Relative Abundance") +
    scale_fill_manual(values = c("#87B60D", "#EEC76A", "#EF8E1B", "#C74543"))
  print(p)
}
dev.off()

message("Visualization Completed.")