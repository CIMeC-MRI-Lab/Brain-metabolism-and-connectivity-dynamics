---
# README - How to run Behavioral Statistical Analysis

francesca.saviola@epfl.ch
barbara.cassone@unimib.it
asia.ferrari@unige.ch

20252002

This repository contains scripts to perform behavioral statistical analysis on reaction times (RTs), accuracy, and d-prime metrics derived from n-back tasks. Follow the steps below to execute the analysis in the correct order.

---

## Scripts Usage

### Step 1: Extract Behavioral Metrics  
**Script**: `extract_behav_metrics.m`  
This script extracts RTs, accuracy, and d-prime metrics from raw data.

- **Input Files**:  
  - `info_demographics_35s.csv`  
  - `<server_path>\fMRS_fMRI_NBACK_project\sub*\sub-*_task-*nback_*.xlsx`

- **Output Files**:  
  - `EIB_behav_acc_dprime_meanRT_fMRI_35s.csv`  
  - `EIB_behav_acc_dprime_meanRT_fMRS_12s.csv`  
  - `EIB_behav_RT_wide_0back_fMRI_35s.csv`  
  - `EIB_behav_RT_wide_0back_fMRS_12s.csv`  
  - `EIB_behav_RT_wide_1back_fMRI_35s.csv`  
  - `EIB_behav_RT_wide_1back_fMRS_12s.csv`  
  - `EIB_behav_RT_wide_2back_fMRI_35s.csv`  
  - `EIB_behav_RT_wide_2back_fMRS_12s.csv`

---

### Step 2: Analyze Accuracy Differences Between N-Back Tasks  
**Script**: `EIB_behav_acc_fMRI_35s_fMRS_12s.Rmd`

- **Input Files**:  
  - `EIB_behav_acc_dprime_meanRT_fMRI_35s.csv`  
  - `EIB_behav_acc_dprime_meanRT_fMRS_12s.csv`

- **Output File**:  
  - `EIB_behav_acc_fMRI_35s_fMRS_12s.html`

---

### Step 3: Analyze D-Prime Differences Between N-Back Tasks  
**Script**: `EIB_behav_dprime_fMRI_35s_fMRS_12s.Rmd`

- **Input Files**:  
  - `EIB_behav_acc_dprime_meanRT_fMRI_35s.csv`  
  - `EIB_behav_acc_dprime_meanRT_fMRS_12s.csv`

- **Output File**:  
  - `EIB_behav_dprime_fMRI_35s_fMRS_12s.html`

---

### Step 4: Analyze RT Differences Between N-Back Tasks  
**Script**: `EIB_behav_RT_wide_fMRI_35s_fMRS_12s.Rmd`

- **Input Files**:  
  - `EIB_behav_RT_wide_0back_fMRI_35s.csv`  
  - `EIB_behav_RT_wide_0back_fMRS_12s.csv`  
  - `EIB_behav_RT_wide_1back_fMRI_35s.csv`  
  - `EIB_behav_RT_wide_1back_fMRS_12s.csv`  
  - `EIB_behav_RT_wide_2back_fMRI_35s.csv`  
  - `EIB_behav_RT_wide_2back_fMRS_12s.csv`

- **Output File**:  
  - `EIB_behav_RT_wide_fMRI_35s_fMRS_12s.html`

---

### Step 5: Analyze the Effect of EIB Metrics on Accuracy  
**Script**: `EIB_behav_EIB_metrics_acc_fMRS_12s.Rmd`

- **Input Files**:  
  - `EIB_behav_acc_dprime_meanRT_fMRS_12s.csv`  
  - `info_EIB_fMRS_12s.csv`

- **Output File**:  
  - `EIB_behav_EIB_metrics_acc_fMRS_12s.html`

---

### Step 6: Analyze the Effect of EIB Metrics on RTs  
**Script**: `EIB_behav_EIB_metrics_meanRT_fMRS_12s.Rmd`

- **Input Files**:  
  - `EIB_behav_acc_dprime_meanRT_fMRS_12s.csv`  
  - `info_EIB_fMRS_12s.csv`

- **Output File**:  
  - `EIB_behav_EIB_metrics_meanRT_fMRS_12s.html`

---

### Step 7: Analyze the Effect of EIB Metrics on D-Prime
**Script**: `EIB_behav_EIB_metrics_dprime_fMRS_12s.Rmd`

- **Input Files**:  
  - `EIB_behav_acc_dprime_meanRT_fMRS_12s.csv`  
  - `info_EIB_fMRS_12s.csv`

- **Output File**:  
  - `EIB_behav_EIB_metrics_dprime_fMRS_12s.html`

---

## Notes
1. Ensure all input files are correctly placed in their respective directories before running each script.
2. Scripts ending in `.m` are MATLAB files, while `.Rmd` files are R Markdown scripts used for statistical analysis and report generation.
3. Modify file paths as needed to match your system's directory structure.

--- 
