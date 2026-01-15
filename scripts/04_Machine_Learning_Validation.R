# ==============================================================================
# Script 04: Machine Learning Model Construction & Validation
# Description: Training 8 ML models and validating on external cohort
# ==============================================================================

library(caret)
library(pROC)
library(ranger)
library(gbm)
library(xgboost)

data_dir <- "data"
out_dir <- "output"
set.seed(672)

# ==================== Data Preparation ====================
# Load Data
train_raw <- read.table(file.path(data_dir, "model_data_origin.txt"), header = TRUE, row.names = 1)
valid_raw <- read.table(file.path(data_dir, "validation_targeted_data.txt"), header = TRUE, row.names = 1)
groups_train <- read.table(file.path(data_dir, "model_group.txt"), header = TRUE, row.names = 1)
groups_valid <- read.table(file.path(data_dir, "vali_group.txt"), header = TRUE, row.names = 1)

# Load Selected Features (12 Markers)
features_12 <- colnames(read.table(file.path(data_dir, "12M.txt"), header = TRUE))

# Subset Data
data_train <- train_raw[, features_12]
data_train$Class <- factor(groups_train$Group1, levels = c("NGC", "GC")) # Ensure factor

data_valid <- valid_raw[, features_12]
data_valid$Class <- factor(groups_valid$Group1, levels = c("NGC", "GC"))

# Train/Test Split (Internal)
inTrain <- createDataPartition(data_train$Class, p = 0.75, list = FALSE)
training_set <- data_train[inTrain, ]
testing_set <- data_train[-inTrain, ]

# ==================== Model Training Loop ====================
ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)

# List of models to train
models <- list(
  RF = list(method = "rf", tuneGrid = expand.grid(mtry = c(2,3,4))),
  GBM = list(method = "gbm", tuneGrid = expand.grid(n.trees = 100, interaction.depth = 3, shrinkage = 0.1, n.minobsinnode = 10)),
  SVM = list(method = "svmRadial", tuneGrid = expand.grid(C = 1, sigma = 0.221)),
  LR = list(method = "glm", tuneGrid = NULL)
  # Add other models as needed...
)

results <- list()

for(alg in names(models)) {
  message(paste("Training", alg, "..."))
  fit <- train(Class ~ ., data = training_set, 
               method = models[[alg]]$method,
               trControl = ctrl,
               tuneGrid = models[[alg]]$tuneGrid,
               metric = "ROC")
  
  # Predictions
  prob_train <- predict(fit, training_set, type = "prob")[, "GC"]
  prob_test <- predict(fit, testing_set, type = "prob")[, "GC"]
  prob_valid <- predict(fit, data_valid, type = "prob")[, "GC"]
  
  # Store AUCs
  results[[alg]] <- data.frame(
    Algorithm = alg,
    Train_AUC = roc(training_set$Class, prob_train)$auc,
    Test_AUC = roc(testing_set$Class, prob_test)$auc,
    Valid_AUC = roc(data_valid$Class, prob_valid)$auc
  )
}

# ==================== Output Results ====================
final_res <- do.call(rbind, results)
write.csv(final_res, file.path(out_dir, "04_ML_Performance_Summary.csv"), row.names = FALSE)

# Plot ROC for Best Model (e.g., RF) in Validation
pdf(file.path(out_dir, "04_Validation_ROC_RF.pdf"), width = 5, height = 5)
rf_model <- train(Class ~ ., data = training_set, method = "rf", trControl = ctrl)
pred_val <- predict(rf_model, data_valid, type = "prob")[, "GC"]
roc_obj <- roc(data_valid$Class, pred_val)
plot(roc_obj, main = paste("External Validation (AUC =", round(roc_obj$auc, 3), ")"), col = "blue")
dev.off()


message("ML Pipeline Completed.")
