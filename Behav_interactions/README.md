---
# README - How to run Behavioral Statistical Analysis for Interactions with EIB (static and kinetics) and FPN-CAPs

francesca.saviola@epfl.ch
barbara.cassone@unimib.it
asia.ferrari@unige.ch

20252002

This repository contains scripts to analyze the interactions between neurometabolites (GABA, Glx), behavioral metrics (reaction times, accuracy, d-prime), and FPN-CAP network variables (e.g., INDegree, OUTDegree, Resilience). Follow the steps below to execute the analyses in the correct order.

---

## Scripts Usage

### Step 1: Analyze Neurometabolites and FPN-CAP Network Variables  
**Script**: `FPN_CAP_EIB_neurometabolites_12s.Rmd`  

This script examines the relationship between neurometabolite data (GABA, Glx), their ratios (EIB), and FPN-CAP network variables such as INDegree, OUTDegree, and Resilience.

- **Input File**:  
  - `EIB_fMRS_FPN_CAP_12s_norest.csv`

- **Output File**:  
  - `FPN_CAP_EIB_neurometabolites_12s.html`

---

### Step 2: Analyze Behavioral Scores and FPN-CAP Network Variables (12s Data)  
**Script**: `FPN_CAP_behavior_12s.Rmd`  

This script explores the relationship between behavioral scores (reaction times, accuracy, d-prime) and FPN-CAP network variables using 12-second data.

- **Input File**:  
  - `EIB_fMRS_FPN_CAP_12s_norest.csv`

- **Output File**:  
  - `FPN_CAP_behavior_12s.html`

---

### Step 3: Analyze Behavioral Scores and FPN-CAP Network Variables (36s Data)  
**Script**: `FPN_CAP_behavior_36s.Rmd`  

This script investigates the relationship between behavioral scores (reaction times, accuracy, d-prime) and FPN-CAP network variables using 36-second data.

- **Input File**:  
  - `EIB_fMRS_FPN_CAP_36s.csv`

- **Output File**:  
  - `FPN_CAP_behavior_36s.html`

---

## Notes
1. Ensure all input files are correctly placed in their respective directories before running each script.
2. All scripts are written in R Markdown (`.Rmd`) for statistical analysis and report generation.
3. Modify file paths in the scripts as needed to match your system's directory structure.

---
