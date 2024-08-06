# First_Step_Muts
Code supplement for "Predicting the First Steps of Evolution in Randomly Assembled Microbial Communities."

All files are either MATLAB executables or .mat files containing variables to be read into MATLAB. When run, the code should read in simulation results and reproduce all figures in the paper, as described below.

System Requirements: MATLAB (tested in R2022a)

FOLDER DIRECTORY:

Simulations/

Functions which can be called to produce raw simulation results.

FILE DIRECTORY:

post_processing.m

Each section will read in raw simulation results (in the Raw_Results folder) and output post-processed data suitable for figure generation (in the Processed_Results folder). This demonstrates how data was analyzed and theory curves generated for each figure. Expected runtime on a standard desktop is several hours for the full Raw_Results dataset.

post_processing_small.m

Same as above, but for truncated raw results (in the Raw_Results_Small folder). Suitable for analysis of the code, but reduced replicates are insufficient to generate accurate figures. Expected runtime is about a minute.

generate_figures.m

Each section will produce a panel of each figure, reading in data from the Processed_Results folder. Runtime should be less than a minute.

ZENODO DATA DIRECTORY:

https://zenodo.org/records/13207468

Raw_Results/

Raw simulation results used to generate figures.

Raw_Results_Small/

Raw simulation results truncated to five replicates, useful for quick read-in when working with data_processing. However, the low number of replicates is insufficient for generating most figure panels.

Processed_Results/

Post-processed simulation results and theory predictions for fastest figure generation.
