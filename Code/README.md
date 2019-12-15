# Code

This directory contains the MCMC code (mcmc.c) that can be used to perform analyses using the Bayesian WAND model described in the paper. The code is currently configured to analyse the NBA data and can be compiled within a terminal environment using the command

`gcc mcmc.c -lgsl -lgslcblas -lm -O3 -o mcmc_exe`

and executed using the command

`./mcmc_exe`

**WARNING:** Executing the code above will create the directory /outputs in the current working directory if it does not already exist. Alternatively if the directory exists then any previous output files will be overwritten.

## Configuration

The WAND model contains many parameters and the MCMC code must be configured to run on different datasets. Lines 20 to 30 of mcmc.c should be used to specify the number of rankers/entities, the data file location, and the number of desired MCMC iterations. Lines 55 to 72 allow the user to specify different ranking types (complete/partial/top); additional information on these types of rankings are given within the Data directory. The prior distribution is specified though lines 75 to 107. Additional algorithm options, for example, whether it is desired to use the Weighted or standard Plackett-Luce model, or if a fixed seed is to be used, can be defined using lines 110 to 117. The code is configured to initialise at random draw from the prior distribution, however, this can be changed in lines 121 to 222 if desired. Detailed comments within mcmc.c should help guide the user.
