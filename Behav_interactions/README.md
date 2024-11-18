## Behavioral statistical analysis for interactions with EIB static and kinetics and FPN-CAPs

## Scripts USAGE
Scripts need to be run in the following order:

1. FPN_CAP_EIB_neurometabolites_12s.Rmd (run statistical analysis on the relationship between neurometabolites data (GABA, Glx) and its ratio (EIB); and FPN-CAP network variables (INDegree, OUTDegree, Resilience, etc.))

	input files: 
		- EIB_fMRS_FPN_CAP_12s_norest.csv

	output files:
		- FPN_CAP_EIB_neurometabolites_12s.html

2. FPN_CAP_behavior_12s.Rmd (run statistical analysis on the relationship between behavioral scores (Reaction Times, Accuracy, d'); and FPN-CAP network variables (INDegree, OUTDegree, Resilience, etc.))

	input files: 
		- EIB_fMRS_FPN_CAP_12s_norest.csv

	output files:
		- FPN_CAP_behavior_12s.html

3. FPN_CAP_behavior_36s.Rmd (run statistical analysis on the relationship between behavioral scores (Reaction Times, Accuracy, d'); and FPN-CAP network variables (INDegree, OUTDegree, Resilience, etc.))

	input files: 
		- EIB_fMRS_FPN_CAP_36s.csv

	output files:
		- FPN_CAP_behavior_36s.html


