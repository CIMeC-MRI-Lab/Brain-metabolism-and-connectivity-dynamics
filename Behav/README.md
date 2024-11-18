## Behavioral statistical analysis
## Scripts USAGE

Scripts need to be run in the following order:

1. extract_behav_metrics.m (to extract RTs, accuracy and d prime from raw data)

	input files: 
		- info_demographics_35s.csv
		- <server_path>\fMRS_fMRI_NBACK_project\sub*\sub-*_task-*nback_*.xlsx

	output files:
		- EIB_behav_acc_dprime_meanRT_fMRI_35s.csv
		- EIB_behav_acc_dprime_meanRT_fMRS_12s.csv
		- EIB_behav_RT_wide_0back_fMRI_35s.csv
		- EIB_behav_RT_wide_0back_fMRS_12s.csv
		- EIB_behav_RT_wide_1back_fMRI_35s.csv
		- EIB_behav_RT_wide_1back_fMRS_12s.csv
		- EIB_behav_RT_wide_2back_fMRI_35s.csv
		- EIB_behav_RT_wide_2back_fMRS_12s.csv

2. EIB_behav_acc_fMRI_35s_fMRS_12s.Rmd (run statistical analysis on accuracy differences between n-back tasks)

	input files: 
		- EIB_behav_acc_dprime_meanRT_fMRI_35s.csv
		- EIB_behav_acc_dprime_meanRT_fMRS_12s.csv

	output files:
		- EIB_behav_acc_fMRI_35s_fMRS_12s.html

3. EIB_behav_dprime_fMRI_35s_fMRS_12s.Rmd (run statistical analysis on d prime differences between n-back tasks)

	input files: 
		- EIB_behav_acc_dprime_meanRT_fMRI_35s.csv
		- EIB_behav_acc_dprime_meanRT_fMRS_12s.csv

	output files:
		- EIB_behav_dprime_fMRI_35s_fMRS_12s.html

4. EIB_behav_RT_wide_fMRI_35s_fMRS_12s.Rmd (run statistical analysis on RTs differences between n-back tasks)

	input files: 
		- EIB_behav_RT_wide_0back_fMRI_35s.csv
		- EIB_behav_RT_wide_0back_fMRS_12s.csv
		- EIB_behav_RT_wide_1back_fMRI_35s.csv
		- EIB_behav_RT_wide_1back_fMRS_12s.csv
		- EIB_behav_RT_wide_2back_fMRI_35s.csv
		- EIB_behav_RT_wide_2back_fMRS_12s.csv

	output files:
		- EIB_behav_RT_wide_fMRI_35s_fMRS_12s.html

5. EIB_behav_EIB_metrics_acc_fMRS_12s.Rmd (run statistical analysis on the effect of EIB metrics on accuracy)

	input files: 
		- EIB_behav_acc_dprime_meanRT_fMRS_12s.csv
		- info_EIB_fMRS_12s.csv

	output files:	
		- EIB_behav_EIB_metrics_acc_fMRS_12s.html

6. EIB_behav_EIB_metrics_meanRT_fMRS_12s.Rmd (run statistical analysis on the effect of EIB metrics on RTs)

	input files: 
		- EIB_behav_acc_dprime_meanRT_fMRS_12s.csv
		- info_EIB_fMRS_12s.csv

	output files:	
		- EIB_behav_EIB_metrics_meanRT_fMRS_12s.html

6. EIB_behav_EIB_metrics_dprime_fMRS_12s.Rmd (run statistical analysis on the effect of EIB metrics on d prime)

	input files: 
		- EIB_behav_acc_dprime_meanRT_fMRS_12s.csv
		- info_EIB_fMRS_12s.csv

	output files:	
		- EIB_behav_EIB_metrics_dprime_fMRS_12s.html

