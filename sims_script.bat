#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=10:mem=4gb

module load anaconda3/rsims_panth
source activate rsims_panth

cd $HOME/ICTU/Sample_size/Rcode

nchunks=10 #friendly many
Trials=1000 # Big split -> links to m
m=$PBS_ARRAY_INDEX # model index

Rscript Simulations_wrapper_PANTHER.R $nchunks $Trials


