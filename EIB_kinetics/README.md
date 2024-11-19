# README - How to run MRS analysis

stefano.tambalo@unitn.it

20241111

## Organization of raw data folders:

```
<path-to-data-folder>/<keyword>/MRI-MRS/<Load>/*WREF*|*EDIT*
```
where:

path-to-data-folder: root path of your data;

keyword*: arbitrary string included in folder names to select
	the batch of folders to process;

MRI-MRS: fixed folder name (can be edited in the main script)

Load: arbitrary folder names that contain Load-dependent experiments.
	They need to be numerically sorted in order to be processed correctly
	(e.g. load_001/load_002/load_003 and so on...)

*WREF*|*EDIT*: fixed folder names with water reference and
	ON/OFF edited spectra pairs, respectively.


## Organization of Matlab library:

The entire set of scripts can be copied in any folder in the file system
and must be added to the Matlab search path


## Running the scripts:

From Matlab prompt, move to:

```
path-to-data-folder
```

and call the desired function.

**Static MRS**:

```
<matlab prompt>> lnifmri_mrs_static_analysis(fittingmethod, keyword)
```

where

fittingmethod = [1 = OSPREY | 2 = GANNET]

keyword = *keyword* (see "Organization of raw data folders")

**Dynamic (sliding window) MRS**:

```
<matlab prompt>> lnifmri_mrs_dynamic_analysis(fittingmethod, keyword, wsize, wstep)
```

fittingmethod = [1 = OSPREY | 2 = GANNET]

keyword = *keyword* (see "Organization of raw data folders")

wsize = sliding window size (numbers of ON/OFF edited pairs)

wstep = step size of the sliding window (numbers of ON/OFF edited pairs)

## Sample dataset

We included a reduced version of the Matlab workspace to test the code. This set of variables includes all the relevant information to run the processing of MRS data (fitted spectra, subject codes and strings to generate output files). The workspace can be loaded from Matlab prompt:

```
<matlab prompt>> load('dummyworkspace.mat')
```

Then, code cells starting from line XXX can be executed within the Matlab IDE. Before running the code, please create a dummy results folder in the current working directory:

```
mkdir deleteme_deleteme
```