# Matlab directory

This directory contains all analyses done in Matlab including the computational modeling and the voxel-based morphometry. 
Additionally the directory also contains the scripts for the experiment and the raw anonymized data.


## Directory Structure

The directory is structured in the following way:

```         
Thermal Pain Learning/matlab
├── README.md             # overview of the directory.
│
├── created files                         # Directory containing files created from the computational modeling in MATLAB.
│   └── ... 
│ 
├── csv_files2                            # Directory containing all the raw data stored from the experimental script.
│    └── Data                             
│
├── matlab data                           # Directory containing all the raw data but as .mat files (made in the main matlab scripts).
│    └── Data                    
│
├── model_parameter_recovery.m            # Main MATLAB script for running the model and parameter recovery for the computational modeling.
│
├── pain_main.m                           # Main MATLAB script for running the computational analyses.
│
├── painLearning_ver4                     # Directory containing MATLAB scripts for running the experiment with instructions.
│   └── ... 
│
├── scripts                               # Directory containing written MATLAB scripts used in the Main MATLAB script.
│   └── ... 
│
├── tapas                                 # Directory containing the TAPAS toolbox together with written scripts used in the Main MATLAB script.
│   └── ... 
│
└── VBQ                                   # Directory containing files and script to run the VBQ brain analysis.

```


