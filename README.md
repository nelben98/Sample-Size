# Sample size calculation PANTHER trial

This is the README file for this modification of the [BATSS](https://batss-dev.github.io/BATSS/) R package, used
in particular to calculate Sample Size calculations of [PANTHER clinical trial](https://panthertrial.org/). 

## This package contains the following folders:

* [Rcode](https://github.com/nelben98/Sample-Size/tree/master/Rcode): All code to be ran and background functions
* [Results](https://github.com/nelben98/Sample-Size/tree/master/Results): Some results (not pushed to repo)
* [log](https://github.com/nelben98/Sample-Size/tree/master/log): Upload of logs when needs to be saved.
* [excel distributions](https://github.com/nelben98/Sample-Size/tree/master/excel_distributions): Underlying probability distribution functions.
* [HPC](https://github.com/nelben98/Sample-Size/tree/master/HPC): Docker definition file - to help run this packages on the HPC.

## How to use:

Running on the HPC is complex and requires some package building and push to the server - this will be uploaded in this repo in the future with a full set of instructions.

For local running - once the INLA and BATSS packages have been installed (ie the complex stage in HPC) - simply run the files on the **Rcode ** folder.

Components:

* bats_glm_breakdown.R: This file contains all the background functions to run the BATSS package using POM - (Not for direct use, back-end only).
* Sample_size_BATSS.R: This will run one set of comparison between two betas (underlying distributions of Days).
* Simulations wrapper: These files are more elaborate versions of Sample_size_BATSS, using the bats_glm_breakdown functions on a group of betas (distributions) - used to compute in parallel.
* Ed_Code_Bayesian_seqdes.R : This file compiles the originial version of the sample size calculations (using Ed's original code).
