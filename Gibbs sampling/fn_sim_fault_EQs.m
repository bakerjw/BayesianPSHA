function [faultX, faultY, eqX, eqY] = fn_sim_fault_EQs(fault, numSimsOnFault, numSimsAreal)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% On-fault seismicity
% simulate error terms
simAngle = rand(numSimsOnFault,1) *2*pi;
% simDist = rand(numSimsOnFault,1) * fault.locUnc; % uniform location error
simDist = normrnd(0, fault.locUnc, numSimsOnFault,1); % Gaussian location error


% simulate location on fault
simLoc = rand(numSimsOnFault,1) * fault.Length;

% map results to Cartesian coordinates

    % (corner)    + (length along fault)        +  (detection error)
eqX = fault.XY(1) + simLoc .* cos(fault.Angle)  +  simDist .* cos(simAngle);
eqY = fault.XY(2) + simLoc .* sin(fault.Angle)  +  simDist .* sin(simAngle);

% get Cartesian coordinates of the fault itself
faultX = [fault.XY(1) fault.XY(1) + fault.Length*cos(fault.Angle)];
faultY = [fault.XY(2) fault.XY(2) + fault.Length*sin(fault.Angle)];

%% Areal seismicity
simAngle = rand(numSimsAreal,1)*2*pi;
simDist = rand(numSimsAreal,1)*fault.areaRadius;
arealX = simDist .* cos(simAngle);
arealY = simDist .* sin(simAngle);

eqX = [eqX; arealX];
eqY = [eqY; arealY];


end

