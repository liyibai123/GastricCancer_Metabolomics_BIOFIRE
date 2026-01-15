# ==============================================================================
# Script 02: BIO-FIRE Algorithm Implementation

library(NMF)
library(Boruta)
library(ggplot2)
library(dplyr)

data_dir <- "data"
out_dir <- "output"
set.seed(4028)

# ==================== Step 1: NMF Functional Module Identification ====================
# Load normalized data
nmf_input <- read.table(file.path(data_dir, "model_data_norm.txt"), header = TRUE, sep = "\t")
nmf_matrix <- as.matrix(nmf_input[, -which(names(nmf_input) == "group")]) # Remove group col

# Run NMF (Rank=5 based on prior optimization)
# Note: nrun reduced to 10 for demo; use nrun=100 for full reproduction
nmf_res <- nmf(nmf_matrix, rank = 5, method = "brunet", nrun = 50, seed = 4028)

# Extract features and clusters
features <- extractFeatures(nmf_res, "max")
nmf_clusters <- predict(nmf_res)
cluster_df <- data.frame(Metabolite = names(nmf_clusters), Cluster = nmf_clusters)
write.csv(cluster_df, file.path(out_dir, "02_NMF_Clusters.csv"), row.names = FALSE)

# Plot NMF Map
pdf(file.path(out_dir, "02_NMF_Consensus_Map.pdf"))
consensusmap(nmf_res)
dev.off()

# ==================== Step 2: Boruta Candidate Screening ====================
set.seed(1)
# Ensure group is factor
nmf_input$group <- as.factor(nmf_input$group)

# Run Boruta
boruta_res <- Boruta(group ~ ., data = nmf_input, doTrace = 0, ntree = 500, pValue = 0.001)

# Extract Statistics
boruta_stats <- attStats(boruta_res)
filter_data <- boruta_stats # Keep full stats for boxplot
filter_data_confirm <- boruta_stats[boruta_stats$decision == "Confirmed", ]

# Sort by Importance
filter_data_confirm <- filter_data_confirm[order(filter_data_confirm$meanImp, decreasing = TRUE), ]
filter_data_confirm$Metabolite <- rownames(filter_data_confirm) # Make rowname a column

# Integrate NMF Cluster Info
filter_data_confirm <- merge(filter_data_confirm, cluster_df, by = "Metabolite")
# Re-sort after merge (merge might shuffle order)
filter_data_confirm <- filter_data_confirm[order(filter_data_confirm$meanImp, decreasing = TRUE), ]

# --- Visualization 3.1: Feature Importance Barplot by NMF Cluster ---
p1 <- ggplot(filter_data_confirm, aes(x = reorder(Metabolite, meanImp), 
                                      y = meanImp, fill = as.factor(Cluster))) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  labs(x = "Metabolites", y = "Importance Score", 
       title = "Boruta Feature Selection - Important Metabolites",
       fill = "NMF Cluster") +
  scale_fill_manual(values = c("1" = "#F79E7C", "2" = "#A7DFD6", "3" = "#C0D990", "4" = "#CDBFE2", "5" = "#FEE69D")) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

ggsave(file.path(out_dir, "02_Boruta_Importance_Barplot.pdf"), p1, width = 7, height = 7)

# --- Visualization 3.2: Importance Distribution Boxplot ---
filter_data$Status <- factor(filter_data$decision, levels = c("Confirmed", "Tentative", "Rejected"))

p2 <- ggplot(filter_data, aes(x = Status, y = meanImp, fill = Status)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Confirmed" = "#87B60d", "Tentative" = "#EEC76A", "Rejected" = "#C74543")) +
  labs(title = "Distribution of Metabolite Importance Scores",
       x = "Selection Status", y = "Importance Score (meanImp)") +
  theme_minimal()

ggsave(file.path(out_dir, "02_Boruta_Distribution_Boxplot.pdf"), p2, width = 7, height = 7)

# --- Visualization 3.3: Selection Process Summary ---
filter_summary <- data.frame(
  Step = c("Before Selection", "Confirmed", "Tentative", "Rejected"),
  Metabolites = c(nrow(filter_data),
                  sum(filter_data$decision == "Confirmed"),
                  sum(filter_data$decision == "Tentative"),
                  sum(filter_data$decision == "Rejected"))
)

filter_summary$Step <- factor(filter_summary$Step, levels = c("Before Selection", "Confirmed", "Tentative", "Rejected"))

p4 <- ggplot(filter_summary, aes(x = Step, y = Metabolites, fill = Step)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Before Selection" = "#CDBFE2", "Confirmed" = "#A7DFD6",
                               "Tentative" = "#FEE69D", "Rejected" = "#F79E7C")) +
  geom_text(aes(label = Metabolites), vjust = -0.5, size = 5) +
  labs(title = "Metabolite Selection Process in Boruta", x = "", y = "Number of Metabolites") +
  theme_minimal()

ggsave(file.path(out_dir, "02_Boruta_Selection_Summary.pdf"), p4, width = 7, height = 7)

# ==================== Step 3: Recursive Signature Optimization  ====================
candidates <- read.table(file.path(data_dir, "candidate.txt"), header = TRUE, sep = "\t")

# Sort by Importance (High -> Low)
candidates <- candidates[order(candidates$meanImp, decreasing = TRUE), ]

selected_panel <- c()
covered_clusters <- c()
target_coverage <- 4 

cat(paste("Target Clusters to Cover:", target_coverage, "\n"))

for(i in 1:nrow(candidates)) {
  met <- candidates$Metabolite[i]
  cl <- candidates$group[i]
  
  # Add metabolite to panel
  selected_panel <- c(selected_panel, met)
  
  # Update coverage status
  if(!cl %in% covered_clusters) {
    covered_clusters <- c(covered_clusters, cl)
    cat(sprintf("Step %d: Added %s (New Cluster %s) -> Coverage: %d/%d\n", 
                i, met, cl, length(covered_clusters), target_coverage))
  } else {
    cat(sprintf("Step %d: Added %s (Cluster %s) -> Redundant Coverage\n", i, met, cl))
  }
  
  # STOPPING CRITERION:
  # Stop when the target number of clusters are covered
  if(length(covered_clusters) >= target_coverage) {
    cat(sprintf("\n>>> Optimization Complete. Selected top %d metabolites to cover %d clusters.\n", 
                length(selected_panel), target_coverage))
    break
  }
}
# ==================== Step 4: Save Final Panel ====================
final_df <- data.frame(Metabolite = selected_panel)
write.table(final_df, file.path(out_dir, "02_Final_12_Panel.txt"), 
            row.names = FALSE, quote = FALSE)

cat("BIO-FIRE Algorithm Completed. Final panel saved to output/02_Final_12_Panel.txt\n")
