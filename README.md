#  Thermosensory predictive coding underpins an illusion of pain

## Table of Contents
1. [Introduction](#introduction)
2. [Directory Structure](#directory-structure)
3. [Access](#access)
4. [Reproducibility](#reproducibility)
5. [Usage](#usage)


## Introduction
This is all the scripts used to run and generate the results, analyses plots and the manuscript for the paper:
"Uncertainty in Thermosensory Expectations shapes an Illusion of Pain".
For complete reproducibility the project is embedded with [OSF](https://osf.io/q5z39/) to store the original results and analyses to big for github. 

The computational environment can be reproduced using the [renv](https://rstudio.github.io/renv/articles/renv.html) package for R environments.


## Directory Structure

The repository is structured in the following way:

```         
Thermal Pain Learning/
├── README.md             # overview of the project.
│
├── Figures/              # Figures generated from code to the from the final manuscript.
│   └── ... 
│
├── Manuscripts/          # Folder containing everything needed to recreate the final manuscript in pdf,docx and html.
│   ├── Manuscript.Rmd                # Rmarkdown for the final manuscript.
│   ├── Supplementary material.Rmd    # Rmarkdown for the supplementary material of the final manuscript.
│   └── Kntting files/                # Files used for knitting the manuscript.
│
├── Markdowns/            # Folder containing markdowns for running the main analyses and plots.
│   ├── Analysis.Rmd                  # Rmarkdown for for running the main analyses in the manuscript.
│   ├── Supplementary Analysis.Rmd    # Rmarkdown for the supplementary analyses in the manuscript.
│   ├── plots.Rmd                     # Rmarkdown for the plots used in the manuscript.
│   └── Control analyses.Rmd          # Rmarkdown for the control analyses done.
│
├── Analysis/            # Folder containing R objects for fitted models.
│   ├── Analysis.Rmd                  # Rmarkdown for for running the main analyses in the manuscript.
│   ├── Supplementary Analysis.Rmd    # Rmarkdown for the supplementary analyses in the manuscript.
│   ├── plots.Rmd                     # Rmarkdown for the plots used in the manuscript.
│   └── Control analyses.Rmd          # Rmarkdown for the control analyses done.
│
│
├── Pictures/             # External pictures used to produce the final figures for the manuscript.
│   └── ... 
│ 
│ 
├── Matlab/               # Folder containing MATLAB scripts plots and data from both the computational modeling, VBQ analyses and the task.
│   ├── Matlab data                  # Folder containing all the raw data stored from the experimental script.
│   │    └── ... 
│   ├── pain_main.m                  # Main MATLAB script for running all the computational analyses.
│   ├── painLearning_ver4            # Folder containing MATLAB scripts for running the experiment with instructions.
│   │   └── ... 
│   ├── scripts                      # Folder containing MATLAB written scripts used in the Main MATLAB script.
│   │   └── ... 
│   └── VBQ                          # Folder containing files and script to run the VBQ brain analysis.
│
├── renv/                 # Directory for the R enviroment
│   └── ... 
│
├── Shiny/              # Directory containing all R scripts in the project
│   ├── utility_shiny.R               # Utility functions for the shiny app.
│   ├── Shiny app.R                   # Shiny app for hightlighting the different learning models.
│   ├── prior_predictive_checks.R     # Prior predictive checks for the learning models.
│   ├── Data/                         # Contingency spaces from the task for the shiny app
│   └── Learning models.R             # Markdown with explanations of the learning models with the shiny app embbeded.
│
│
├── scripts/              # Directory containing all R scripts in the project
│   ├── utils.R                       # Utility functions to gather data, extract summary statistics etc.
│   ├── plots.R                       # Functions used to produced the figures in the final manuscript.
│   └── supplementary_functions.R     # Functions to run the supplementary analysis, produce plots, tables etc.
│
└── renv.lock             # File used to download and install all necessary packages in the versions used to generate the manuscript.


```

## Access

To get access to the data for running the analysis a OSF-token is needed, as the data is stored in the following [OSF-project](https://osf.io/q5z39/). It is recommended to make an osf folder that contains an osf.txt file, in the main directory, which on the first line contains the OSF-token. Make sure that this token is not shared or pushed to github. In the current gitignore an osf folder will be ignored.
To get access to the repository users are recommended to clone the repository with the following command in the terminal

```bash
git clone  https://github.com/Body-Pain-Perception-Lab/Thermal-Pain-Learning.git
```

## Reproducibility

To enhance reproducibility this project is setup with R-package "renv"", ensuring the same packages and versions of these are loaded. This means that users should install the "renv" package and after cloning the repository, run the following code in the console (Not terminal).

```r
renv::restore()
```

This will download and install all the needed packages to run the analysis from scratch with the same version of packages used to generate the manuscript. 


## Usage

You can use this repository in various ways. Choose the one that suits your needs: Note you have to complete the [Access](#access) step to get get the data from OSF. 

1. **Method 1 - Running the main analysis line by line**:
   - Open the Markdown folder and go through the [Analysis.Rmd](./Markdowns/Analysis.Rmd) markdown.

2. **Method 2 - Rerunning the plots line by line**:
   - Open the Markdown folder and go through the [Plots.Rmd](./Markdowns/Plots.Rmd) markdown.
   - To view the plotting functions users are encouraged to check out the [plot.R](./scripts/plots.R) script inside the scripts folder.


