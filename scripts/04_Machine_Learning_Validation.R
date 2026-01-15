# ==============================================================================
# Script 04: Machine Learning Model Evaluation

# Install and load required packages
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("caret")) install.packages("caret")
if (!require("pROC")) install.packages("pROC")
if (!require("ranger")) install.packages("ranger")
if (!require("ggsci")) install.packages("ggsci")
if (!require("svglite")) install.packages("svglite")
if (!require("eoffice")) install.packages("eoffice")
if (!require("RColorBrewer")) install.packages("RColorBrewer")

library(tidyverse)
library(caret)
library(pROC)
library(ranger)
library(ggsci)
library(svglite)
library(eoffice)
library(RColorBrewer)
library(ggplot2)

# Set Directories
data_dir <- "data"
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir)

# ==================== 1. Data Preparation ====================

# Read original absolute quantitative data
model_data_origin <- read.table(file.path(data_dir, "model_data_origin.txt"), header = TRUE, sep = "\t", check.names = FALSE)
rownames(model_data_origin) <- model_data_origin$Index
model_data_origin <- model_data_origin[, -1]

vali_data_origin <- read.table(file.path(data_dir, "validation_targeted_data.txt"), header = TRUE, sep = "\t", check.names = FALSE)
rownames(vali_data_origin) <- vali_data_origin$Index
vali_data_origin <- vali_data_origin[, -1] 

# ==================== 2. Feature Selection for Machine Learning ====================

# Read selected metabolites for modeling
# 12 metabolites from BIO-FIRE strategy
Metabolite_12M <- read.table(file.path(data_dir, "12M.txt"), header = TRUE, sep = "\t", check.names = FALSE)
Metabolite_12M <- Metabolite_12M[, -1]
Metabolite_12M <- colnames(Metabolite_12M)

# ==================== 3. Prepare 12-Metabolite Dataset ====================
# 3.1 Prepare Modeling Set (Subset 83 -> 12 for Model Training and Testing)
model_data_12M <- model_data_origin %>% select(all_of(Metabolite_12M))
# 3.2 Prepare Validation Set (External Verification)
missing_cols <- setdiff(Metabolite_12M, colnames(vali_data_origin))
if(length(missing_cols) > 0) {
  stop("Error: The 12 targeted metabolites are missing in the validation data file: ", paste(missing_cols, collapse=", "))
}
vali_data_12M <- vali_data_origin %>% select(all_of(Metabolite_12M))

# 3.3 Merge Clinical Group Info
# Load Group Info
model_group <- read.table(file.path(data_dir, "model_group.txt"), header = TRUE, sep = "\t", check.names = FALSE)
vali_group <- read.table(file.path(data_dir, "vali_group.txt"), header = TRUE, sep = "\t", check.names = FALSE)

# Process Group Info
rownames(model_group) <- model_group$Index
rownames(vali_group) <- vali_group$Index

# Keep only necessary Group columns
model_group_target <- model_group %>% select(Group1) 
vali_group_target <- vali_group %>% select(Group1)

# Final Merge for Training/Testing
model_data_final <- merge(model_data_12M, model_group_target, by = "row.names")
rownames(model_data_final) <- model_data_final$Row.names
model_data_final <- model_data_final[, -1] 

# Final Merge for External Validation
vali_data_final <- merge(vali_data_12M, vali_group_target, by = "row.names")
rownames(vali_data_final) <- vali_data_final$Row.names
vali_data_final <- vali_data_final[, -1] 

# Verify Dimensions
cat("\nDataset Dimensions Check:\n")
cat("Modeling Set (83M->12M):", dim(model_data_final)[1], "samples,", dim(model_data_final)[2]-1, "features.\n")
cat("Validation Set (12M):", dim(vali_data_final)[1], "samples,", dim(vali_data_final)[2]-1, "features.\n")

# ==================== 4. Machine Learning Model Evaluation ====================

# Set seed for reproducibility
set.seed(672)

# Split data into training and testing sets (75% training, 25% testing)
train_indices <- createDataPartition(model_data_final$Group1, p = 0.75, list = FALSE)
train_data <- model_data_final[train_indices, ]
test_data <- model_data_final[-train_indices, ]

# Configure cross-validation settings
fitControl <- trainControl(method = "cv", 
                           number = 10, 
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary)

#Define Hyperparameter Search Grids
hyperparameter_grids <- list(
  "RandomForest" = expand.grid(mtry = c(2, 7, 12)),
  "GBM" = expand.grid(interaction.depth = c(1, 2, 3), n.trees = c(50, 100, 150), shrinkage = 0.1, n.minobsinnode = 10),
  "NaiveBayes" = expand.grid(usekernel = c(TRUE, FALSE), fL = 0, adjust = 1),
  "LogisticRegression" = NULL, 
  "SVM" = expand.grid(C = c(0.25, 0.50, 1.00), sigma = 0.221435),
  "XGBoost" = expand.grid(nrounds = c(50, 100, 150), max_depth = c(1, 2, 3), eta = c(0.3, 0.4), gamma = 0, colsample_bytree = c(0.6, 0.8), min_child_weight = 1, subsample = c(0.5, 0.75, 1.0)),
  "KNN" = expand.grid(kmax = c(5, 7, 9), distance = 2, kernel = "optimal"),
  "BaggedCART" = NULL 
)

# Model Training Function 
train_and_evaluate_model <- function(model_name, method, train_data, test_data, vali_data, fitControl, tune_grid = NULL) {
  
  cat("Training", model_name, "model with grid search...\n")
  
  # Train model 
  if (is.null(tune_grid)) {
    model <- train(Group1 ~ .,
                   data = train_data,
                   method = method,
                   trControl = fitControl,
                   metric = "Accuracy") 
  } else {
    model <- train(Group1 ~ .,
                   data = train_data,
                   method = method,
                   trControl = fitControl,
                   tuneGrid = tune_grid, 
                   metric = "Accuracy")
  }
  
  # Calculate AUC for training set
  train_pred <- predict(model, train_data, type = "prob")
  train_roc <- roc(train_data$Group1, train_pred$GC)
  train_auc <- round(train_roc$auc, 3)
  
  # Calculate AUC for testing set
  test_pred <- predict(model, test_data, type = "prob")
  test_roc <- roc(test_data$Group1, test_pred$GC)
  test_auc <- round(test_roc$auc, 3)
  
  # Calculate AUC for validation set
  vali_pred <- predict(model, vali_data, type = "prob")
  vali_roc <- roc(vali_data$Group1, vali_pred$GC)
  vali_auc <- round(vali_roc$auc, 3)
  
  return(list(
    model = model,
    train_roc = train_roc,
    train_auc = train_auc,
    test_roc = test_roc,
    test_auc = test_auc,
    vali_roc = vali_roc,
    vali_auc = vali_auc
  ))
}

# Model Training
models_config <- list(
  "RandomForest" = list(method = "rf", name = "RandomForest"),
  "GBM" = list(method = "gbm", name = "GBM"),
  "NaiveBayes" = list(method = "nb", name = "NaiveBayes"),
  "LogisticRegression" = list(method = "glm", name = "LogisticRegression"),
  "SVM" = list(method = "svmRadial", name = "SVM"),
  "XGBoost" = list(method = "xgbTree", name = "XGBoost"),
  "KNN" = list(method = "kknn", name = "KNN"),
  "BaggedCART" = list(method = "treebag", name = "BaggedCART")
)

cat("Training 8 machine learning models:\n")
cat(paste(names(models_config), collapse = ", "), "\n\n")

# Model Training Loop
model_results <- list()

for (model_name in names(models_config)) {
  config <- models_config[[model_name]]
  current_grid <- hyperparameter_grids[[model_name]]
  
  model_results[[model_name]] <- train_and_evaluate_model(
    model_name = config$name,
    method = config$method,
    train_data = train_data,
    test_data = test_data,
    vali_data = vali_data_final,
    fitControl = fitControl,
    tune_grid = current_grid  
  )
}

# ROC Curve Visualization Function
create_roc_plot <- function(roc_data, auc_values, title, filename) {
  
  custom_palette <- colorRampPalette(brewer.pal(12, "Paired"))(8)
  
  p <- ggplot() +  
    geom_path(data = roc_data, aes(x = fpr, y = tpr, color = Model), size = 0.8) +
    labs(x = "False Positive Rate", y = "True Positive Rate", color = "Model") +
    scale_color_manual(values = custom_palette) + 
    ggtitle(title) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 12, color = "black", face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_line(color = "black")
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 0.8)
  
  y_positions <- seq(0.50, 0.15, by = -0.05)  
  model_names <- names(auc_values)
  
  for (i in 1:length(model_names)) {
    p <- p + annotate("text", x = 0.8, y = y_positions[i], 
                      label = paste0(model_names[i], " (AUC = ", auc_values[[model_names[i]]], ")"), 
                      color = "black", size = 3.5)
  }
  
  topptx(p, file = file.path(out_dir, filename), width = 7, height = 5)
  cat("Saved ROC plot:", filename, "\n")
  return(p)
}

# Prepare ROC Data
roc_data_train <- do.call(rbind, lapply(names(model_results), function(model_name) {
  result <- model_results[[model_name]]
  data.frame(
    fpr = 1 - result$train_roc$specificities,
    tpr = result$train_roc$sensitivities,
    Model = model_name
  )
}))

roc_data_test <- do.call(rbind, lapply(names(model_results), function(model_name) {
  result <- model_results[[model_name]]
  data.frame(
    fpr = 1 - result$test_roc$specificities,
    tpr = result$test_roc$sensitivities,
    Model = model_name
  )
}))

roc_data_vali <- do.call(rbind, lapply(names(model_results), function(model_name) {
  result <- model_results[[model_name]]
  data.frame(
    fpr = 1 - result$vali_roc$specificities,
    tpr = result$vali_roc$sensitivities,
    Model = model_name
  )
}))

train_auc_values <- sapply(model_results, function(x) x$train_auc)
test_auc_values <- sapply(model_results, function(x) x$test_auc)
vali_auc_values <- sapply(model_results, function(x) x$vali_auc)

# Generate ROC Plots 
create_roc_plot(roc_data_train, train_auc_values,
                "ROC Curves - Training Set (8 Machine Learning Models)",
                "train_roc_8ML_models.pptx")

create_roc_plot(roc_data_test, test_auc_values,
                "ROC Curves - Testing Set (8 Machine Learning Models)",
                "test_roc_8ML_models.pptx")

create_roc_plot(roc_data_vali, vali_auc_values,
                "ROC Curves - Validation Cohort (8 Machine Learning Models)",
                "validation_roc_8ML_models.pptx")

# Performance Summary
performance_summary <- data.frame(
  Model = names(model_results),
  Train_AUC = sapply(model_results, function(x) x$train_auc),
  Test_AUC = sapply(model_results, function(x) x$test_auc),
  Validation_AUC = sapply(model_results, function(x) x$vali_auc)
)

performance_summary <- performance_summary[order(-performance_summary$Validation_AUC), ]
write.csv(performance_summary, file.path(out_dir, "04_ML_Performance_Summary.csv"), row.names = FALSE)

best_model_idx <- which.max(performance_summary$Validation_AUC)
best_model <- performance_summary$Model[best_model_idx]
best_auc <- performance_summary$Validation_AUC[best_model_idx]

cat("\nBest Performing Model:\n")
cat("Model:", best_model, "\n")
cat("Validation AUC:", best_auc, "\n")
