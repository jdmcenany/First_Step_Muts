# First_Step_Muts
Code supplement for "Predicting the First Steps of Evolution in Randomly Assembled Microbial Communities"

This repository will be updated with raw simulation results and code to generate figures soon.

FOLDER DIRECTORY:

Simulations/
Functions which can be called to produce raw simulation results. Note that independently running each of these functions to generate the results in the paper will take some time.

Raw_Results/
Raw simulation results used to generate figures.

Processed_Results/
Post-processed simulation results and theory predictions for fastest figure generation.


FILE DIRECTORY:

data_processing.m
Each section will read in raw simulation results (in the Raw_Results folder) and output post-processed data suitable for figure generation (in the Processed_Results folder). This demonstrates how data was analyzed and theory curves generated for each figure.

generate_figures.m
Each section will produce a panel of each figure, reading in data from the Processed_Results folder.
