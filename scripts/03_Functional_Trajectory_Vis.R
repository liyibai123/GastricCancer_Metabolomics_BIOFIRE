# ==============================================================================
# Script 03: Functional Trajectory Analysis & Biomarker Visualization

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(eoffice)

data_dir <- "data"
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ==================== Part 1: NMF Functional Trajectory Analysis ====================
cat("Starting NMF Trajectory Analysis...\n")

# 1.1 Load Expression Data
raw_exp <- read.table(file.path(data_dir, "model_exp.txt"), header = TRUE, sep="\t", check.names=FALSE)
rownames(raw_exp) <- raw_exp$SampleID
# Extract numeric matrix 
exp_matrix <- as.matrix(raw_exp[, -(1:3)])
model_exp_log <- log2(exp_matrix + 1)
model_nmf_expr <- t(model_exp_log) 
model_nmf_expr <- as.data.frame(model_nmf_expr)

# 1.2 Load Group Info & NMF Clusters
sample_info <- read.csv(file.path(data_dir, "sample_group.csv"))
nmf_groups <- read.csv(file.path(data_dir, "NMF_cluster.csv"), row.names = 1)
# Ensure consistent row names
common_mets <- intersect(rownames(model_nmf_expr), rownames(nmf_groups))
model_nmf_expr <- model_nmf_expr[common_mets, ]
nmf_groups <- nmf_groups[common_mets, , drop=FALSE]

# 1.3 Prepare Long Format Data
model_nmf_expr_long <- model_nmf_expr %>%
  rownames_to_column("Metabolite") %>%
  pivot_longer(cols = -Metabolite, 
               names_to = "Sample", 
               values_to = "Expression") %>%
  left_join(sample_info, by = "Sample") 

# Z-score normalization per metabolite
model_nmf_expr_long <- model_nmf_expr_long %>%
  group_by(Metabolite) %>%
  mutate(Expression_scaled = scale(Expression)) %>%
  ungroup()

# Add Cluster Info
model_nmf_expr_long <- model_nmf_expr_long %>%
  left_join(nmf_groups %>% rownames_to_column("Metabolite"), by = "Metabolite")

# Set Factor Levels
model_nmf_expr_long$Group2 <- factor(model_nmf_expr_long$Group2, 
                                     levels = c("HC", "BGD", "EGC", "AGC"))

# 1.4 Calculate Statistics
# Mean Expression per Metabolite per Group
model_nmf_expr_mean <- model_nmf_expr_long %>%
  group_by(Metabolite, Group2, cluster) %>%  
  summarise(Expression_scaled = mean(Expression_scaled, na.rm = TRUE), .groups = "drop")

# Mean Expression per Cluster per Group (with SE)
cluster_mean_se <- model_nmf_expr_mean %>%
  group_by(cluster, Group2) %>%
  summarise(
    Mean_Expression = mean(Expression_scaled, na.rm = TRUE),
    SE = sd(Expression_scaled, na.rm = TRUE) / sqrt(n()), 
    .groups = "drop"
  )

# ==================== Part 2: Importance-based Plots (Cluster 1-5) ====================
cluster_list <- c(1, 2, 3, 4, 5)
importance_data <- list()

for (cl in cluster_list) {
  file_name <- file.path(data_dir, paste0("NMF_cluster", cl, ".txt"))
  if (file.exists(file_name)) {
    df <- tryCatch({
      read_delim(file_name, delim = "\t", col_names = c("Metabolite", "Importance"), show_col_types = FALSE)
    }, error = function(e) NULL)
    
    if(!is.null(df)) {
      df$Importance <- as.numeric(df$Importance)
      if (all(is.na(df$Importance))) {
        warning(paste("Cluster", cl, "Importance is all NA, skipping normalization"))
      } else {
        # Min-Max Normalization
        min_val <- min(na.omit(df$Importance))
        max_val <- max(na.omit(df$Importance))
        if(max_val > min_val) {
          df$Importance <- (df$Importance - min_val) / (max_val - min_val)
        }
      }
      if(nrow(df) > 0) df <- df[-1, ] 
      
      importance_data[[as.character(cl)]] <- df
    }
  }
}

# Plot Loop
for (cl in cluster_list) {
  cluster_data <- model_nmf_expr_mean %>% filter(cluster == cl)
  cluster_mean_data <- cluster_mean_se %>% filter(cluster == cl)
  
  # Merge Importance info
  if (!is.null(importance_data[[as.character(cl)]])) {
    cluster_data <- cluster_data %>%
      left_join(importance_data[[as.character(cl)]], by = "Metabolite") %>%
      mutate(Importance = ifelse(is.na(Importance), 0, Importance))
  } else {
    cluster_data$Importance <- 0
  }
  
  p <- ggplot() +
    geom_line(data = cluster_data,
              aes(x = Group2, y = Expression_scaled, group = Metabolite, color = Importance),
              alpha = 0.8, linewidth = 0.8) +
    geom_ribbon(data = cluster_mean_data,
                aes(x = Group2, ymin = Mean_Expression - SE, ymax = Mean_Expression + SE, group = 1),
                fill = "gray70", alpha = 0.3) +
    geom_line(data = cluster_mean_data,
              aes(x = Group2, y = Mean_Expression, group = 1),
              color = "black", linewidth = 1.5) +
    geom_text(data = cluster_data %>% filter(Group2 == "AGC"),
              aes(x = Group2, y = Expression_scaled, label = Metabolite, color = Importance),
              hjust = -0.1, vjust = 0, size = 3) +
    scale_color_gradientn(colors = c("#64B5F6", "#42A5F5", "#1E88E5", "#1565C0", "#0D47A1")) +
    theme_bw() +
    labs(x = "Disease Progress", y = "Z-score", title = paste("Cluster", cl, "Trajectory")) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  ggsave(file.path(out_dir, paste0("Cluster_", cl, "_trajectory.pdf")), p, width = 6, height = 4)
  topptx(p, file = file.path(out_dir, paste0("Cluster_", cl, "_with_mean.pptx")), width = 6, height = 4)
}

# ==================== Part 3: Summary Plots (All Clusters) ====================
selected_clusters <- c(1, 2, 3, 4, 5)
filtered_model <- model_nmf_expr_mean %>% filter(cluster %in% selected_clusters)
filtered_se <- cluster_mean_se %>% filter(cluster %in% selected_clusters)

p_summary <- ggplot() +
  geom_line(data = filtered_model, 
            aes(x = Group2, y = Expression_scaled, group = Metabolite, color = factor(cluster)), 
            alpha = 0.3, size = 0.6) + 
  geom_ribbon(data = filtered_se, 
              aes(x = Group2, ymin = Mean_Expression - SE, ymax = Mean_Expression + SE, 
                  group = cluster, fill = factor(cluster)), 
              alpha = 0.3) +  
  geom_line(data = filtered_se, 
            aes(x = Group2, y = Mean_Expression, group = cluster, color = factor(cluster)), 
            size = 1.5) + 
  scale_color_manual(values = c("#C26A4F","#4A9185","#84914E" ,"#7A5299", "#B2903F")) +  
  scale_fill_manual(values = c("#F79E7C","#A7DFD6","#C0D990","#CDBFE2","#FEE69D")) +  
  scale_y_continuous(breaks = seq(-0.8, 0.8, by = 0.2)) +  
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", size = 1.5, fill = NA),  
    panel.grid = element_blank(), 
    axis.text = element_text(face = "bold"),  
    axis.title = element_text(face = "bold"), 
    axis.ticks = element_line(color = "black", size = 1),
    axis.ticks.length = unit(2, "mm")  
  )

ggsave(file.path(out_dir, "NMF_Summary_Trajectory.pdf"), p_summary, width = 6, height = 4)
topptx(p_summary, file = file.path(out_dir, "NMF_tract_summary_1-5.pptx"), width = 6, height = 4)

# ==================== Part 4: Individual Biomarker Boxplots ====================
cat("\nGenerating Individual Boxplots...\n")

# Load selected panel (from previous step or locked file)
if(file.exists(file.path(out_dir, "02_Final_12_Panel.txt"))) {
  selected_panel <- read.table(file.path(out_dir, "02_Final_12_Panel.txt"), header = TRUE)$Metabolite
} else {
  # Fallback if running script independently
  warning("02_Final_12_Panel.txt not found. Please run Script 02 first.")
  selected_panel <- c() 
}

# Load Boxplot Inputs
boxplot_input <- read.table(file.path(data_dir, "boxplot_input.txt"), header = TRUE, sep = "\t", check.names = FALSE)
boxplot_group <- read.table(file.path(data_dir, "boxplot_group.txt"), header = TRUE, sep = "\t", check.names = FALSE)

# Process Data
rownames(boxplot_input) <- boxplot_input$SampleID
boxplot_data <- as.data.frame(t(boxplot_input[, -1]))
# Merge Groups
common_samples <- intersect(rownames(boxplot_data), boxplot_group$Sample)
boxplot_data <- boxplot_data[common_samples, ]
boxplot_data$Group <- boxplot_group$Group[match(common_samples, boxplot_group$Sample)]
boxplot_data$Group <- factor(boxplot_data$Group, levels = c("HC", "BGD", "EGC", "AGC"))

group_colors <- c("HC" = "#87B60d", "BGD" = "#EEC76A", "EGC" = "#EF8E1B", "AGC" = "#C74543")

# Plot Loop
for (met in selected_panel) {
  if(met %in% colnames(boxplot_data)) {
    p <- ggplot(boxplot_data, aes(x = Group, y = .data[[met]], fill = Group)) +
      geom_jitter(aes(color = Group), width = 0.2, size = 2, alpha = 0.4) +
      geom_boxplot(outlier.shape = NA, size = 1.5, alpha = 0.7) +
      labs(title = met, x = "Group", y = "Relative abundance") +
      scale_fill_manual(values = group_colors) +
      scale_color_manual(values = group_colors) +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.line = element_line(colour = "black"), legend.position = "right")
    
    ggsave(file.path(out_dir, paste0("boxplot_", met, ".pdf")), plot = p, width = 7, height = 7)
  }
}

message("Script 03 Completed Successfully.")

