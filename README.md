# BayesianPSHA
Matlab source code to compute dynamically updated seismic hazard. Created by Jack Baker and Abhineet Gupta in support of the following paper.

Baker, J. W., and Gupta, A. (2015). “Bayesian Treatment of Induced Seismicity in Probabilistic Seismic Hazard Analysis.” Bulletin of the Seismological Society of America, (in review).

This repository contains two subfolders. 
-'Change point' contains code to perform change point calculations and reproduce the Example 2 results in the above paper. 
-'Gibbs sampling' contains code to perform Gibbs sampling calculations and reproduce the Example 3 results in the above paper. 

'mainScript.m' in the root directory is a Matlab script that calls the relevant functions in the above two subfolders.