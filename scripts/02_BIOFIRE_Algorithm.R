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
cat("Running Boruta screening...\n")
nmf_input$group <- as.factor(nmf_input$group)
boruta_res <- Boruta(group ~ ., data = nmf_input, doTrace = 0, ntree = 500)

# ==================== Step 3: Recursive Signature Optimization (Core) ====================
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
