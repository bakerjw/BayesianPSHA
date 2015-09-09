% Master script to call the code associated with change point and Gibbs 
% sampling calculations
%
% Created by Jack Baker
% 9/10/2015


%% Oklahoma City rate updating using change point 

clear; close all; clc;

startDir = [pwd];
cd 'Change point/'
mainCP
cd(startDir)


%% Example 3: fault location updating

clear; close all; clc;

startDir = [pwd];
cd 'Gibbs sampling/'
main_Gibbs_sample % sample fault locations
Gibbs_hazard % compute hazard for sampled fault locations
cd(startDir)

