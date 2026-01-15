# ==============================================================================
# Script 02: BIO-FIRE Algorithm Implementation
# Description: NMF Clustering, Boruta Screening, and Recursive Optimization
# ==============================================================================

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
# Ensure group is factor
nmf_input$group <- as.factor(nmf_input$group)

# Run Boruta
boruta_res <- Boruta(group ~ ., data = nmf_input, doTrace = 0, ntree = 500)
boruta_stats <- attStats(boruta_res)
confirmed_vars <- boruta_stats[boruta_stats$decision == "Confirmed", ]
confirmed_vars$Metabolite <- rownames(confirmed_vars)

# Merge with NMF Clusters
candidates <- merge(confirmed_vars, cluster_df, by = "Metabolite")
candidates <- candidates[order(candidates$meanImp, decreasing = TRUE), ]

# Save Candidates
write.csv(candidates, file.path(out_dir, "02_BIOFIRE_Candidates.csv"), row.names = FALSE)

# ==================== Step 3: Recursive Signature Optimization ====================
# Optimization Logic
selected_panel <- c()
covered_clusters <- c()
max_clusters <- length(unique(candidates$Cluster))

cat("Starting Recursive Selection...\n")
for(i in 1:nrow(candidates)) {
  met <- candidates$Metabolite[i]
  cl <- candidates$Cluster[i]
  
  # Logic: Add if it covers a NEW cluster
  if(!cl %in% covered_clusters) {
    selected_panel <- c(selected_panel, met)
    covered_clusters <- c(covered_clusters, cl)
    cat(sprintf("Added: %s (Cluster %s) | Coverage: %d/%d\n", met, cl, length(covered_clusters), max_clusters))
  }
  
  # Stopping Rule
  if(length(covered_clusters) == max_clusters) break
}

# Supplement with top importance if needed (optional logic from original script)
# Or save the core panel
write.table(data.frame(Metabolite=selected_panel), file.path(out_dir, "02_Final_12_Panel.txt"), 
            row.names=FALSE, quote=FALSE)

message("BIO-FIRE Algorithm Completed. Panel saved.")