# BIO-FIRE: Hybrid Metabolomics Framework for Gastric Cancer Diagnosis

This repository contains the R source code and analysis workflow for the manuscript:
**"Multi-phase hybrid metabolomics framework identifies clinically applicable plasma signatures for early detection of gastric cancer"**.

## Repository Structure

The project is organized into three main directories:

- **data/**: Contains processed input data required to reproduce the analysis.
- **scripts/**: Numbered R scripts corresponding to the analysis steps described in the manuscript.
- **output/**: Directory where generated figures and tables will be saved.

## Analysis Workflow

To reproduce the results, run the scripts in the following order:

### 1. Global Landscape Analysis
- **Script:** `scripts/01_Landscape_and_Differential.R`
- **Function:** Performs global metabolic landscape profiling (Heatmap, LDA) and differential abundance analysis.

### 2. BIO-FIRE Biomarker Discovery
- **Script:** `scripts/02_BIOFIRE_Algorithm.R`
- **Function:** Implements the **BIO-FIRE** algorithm (Biomarker Identification and Optimization via Functional and Importance-based Recursive Enhancement).
- **Key Steps:** 
  1. Functional Module Identification (NMF).
  2. Candidate Screening (Boruta).
  3. Recursive Signature Optimization.

### 3. Functional Trajectory Visualization
- **Script:** `scripts/03_Functional_Trajectory_Vis.R`
- **Function:** Visualizes the longitudinal trends of functional modules and generates boxplots for the identified 12-metabolite panel.

### 4. Machine Learning Validation
- **Script:** `scripts/04_Machine_Learning_Validation.R`
- **Function:** Trains and validates 8 machine learning models (Random Forest, SVM, GBM, etc.).
- **Reproducibility Note:** 
  - To maximize the stability of unsupervised NMF clustering, feature selection (BIO-FIRE) was performed on the entire modeling cohort.
  - Model training and hyperparameter tuning were strictly confined to the training set (75%) via 10-fold cross-validation.
  - Final performance was validated on a completely independent external cohort (blinded to feature selection).

## Requirements

The analysis was performed using **R version 4.4.2**. Please ensure the following packages are installed:

```r
install.packages(c("tidyverse", "caret", "pROC", "NMF", "Boruta", 
                   "randomForest", "gbm", "xgboost", "pheatmap", "eoffice",

                   "scatterplot3d", "MASS"))
