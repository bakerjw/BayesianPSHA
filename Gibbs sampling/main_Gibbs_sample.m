% Generate earthquake catalogs with location uncertainty, and use Gibbs
% Sampling to re-estimate fault paramters
%
% Jack Baker
% Created 5/1/2015
% Updated 5/7/2015 to reparameterize the faults
% Last updated 9/10/2015


% plot specs
label_size = 9; % font size for labels
axis_size = 8;  % font size for axis numbers


maxAxis = 30; % axis range for figures
numSims = 1000; % number of iterations of the Gibbs sampler to run

% fault variables
fault.XY = [3 -8]; % cartesian coordinates relative to the origin
fault.Angle = pi/3; % radians counterclockwise relative to West
fault.Length = 8; % length of fault
fault.areaRadius = 25; % radius of area of interest

% compute coordinates of radius
radiusX = fault.areaRadius * cos(0:0.01:2*pi+0.01);
radiusY = fault.areaRadius .* sin(0:0.01:2*pi+0.01);

% store geometry for later
coordsGibbs = [fault.XY fault.Length fault.Angle]; % put real fault values into a format where we can use the PSHA function later
save('faultData.mat', 'coordsGibbs')


%% Case 1
numSimsOnFault = 10;
fault.locUnc = 5; % standard deviation of Gaussian error (in km)
numSimsAreal = 0;
[faultX, faultY, eqX, eqY] = fn_sim_fault_EQs(fault, numSimsOnFault, numSimsAreal);
coords = [mean(eqX), mean(eqY), 5, pi/5, 0.1]; % make the first guess in the vicinity of the simulated earthquakes [x1 y1 length angle pFault]
coordsGibbs = fn_Gibbs_sample(coords(1:4), eqX, eqY, fault, numSims); % ignore P(fault)
idx = numSims/2:numSims/100:numSims-1; % indices of 50 samples to keep
coordsGibbs = coordsGibbs(idx,:);
save('coordsGibbs1.mat', 'coordsGibbs', 'faultX', 'faultY', 'eqX', 'eqY') % save simulations for later analysis
% load coordsGibbs1

% plot EQ locations
figure
subplot(1,2,1);
h1 = plot(eqX, eqY, '.');
hold on
h2 = plot(faultX, faultY, '-r', 'linewidth', 2);
h3 = plot(0, 0, 'ok', 'linewidth', 2);
plot(radiusX, radiusY, '-k', 'linewidth', 1)
axis(maxAxis*[-1 1 -1 1])
axis square
hx = xlabel('East-west distance (km)', 'FontSize', label_size);
hy = ylabel('North-south distance (km)', 'FontSize', label_size);
legh = legend([h3 h1 h2], 'Location of interest', 'Observed earthquake', '(Unknown) fault location');
set(legh, 'Location', 'northeast', 'FontSize', label_size);
set(gca, 'FontSize', axis_size);



%% Case 2
numSimsOnFault = 30;
fault.locUnc = 2; % standard deviation of Gaussian error (in km)
[faultX, faultY, eqX, eqY] = fn_sim_fault_EQs(fault, numSimsOnFault, numSimsAreal);
coords = [mean(eqX), mean(eqY), 5, pi/5, 0.1]; % make the first guess in the vicinity of the simulated earthquakes [x1 y1 length angle pFault]
coordsGibbs = fn_Gibbs_sample(coords(1:4), eqX, eqY, fault, numSims); % ignore P(fault)
idx = numSims/2:numSims/100:numSims-1; % indices of 50 samples to keep
coordsGibbs = coordsGibbs(idx,:);
save('coordsGibbs2.mat', 'coordsGibbs', 'faultX', 'faultY', 'eqX', 'eqY')
% load coordsGibbs2

% plot EQ locations
subplot(1,2,2);
h1 = plot(eqX, eqY, '.');
hold on
h2 = plot(faultX, faultY, '-r', 'linewidth', 2);
h3 = plot(0, 0, 'ok', 'linewidth', 2);
plot(radiusX, radiusY, '-k', 'linewidth', 1)
axis(maxAxis*[-1 1 -1 1])
axis square
hx = xlabel('East-west distance (km)', 'FontSize', label_size);
hy = ylabel('North-south distance (km)', 'FontSize', label_size);
legh = legend([h3 h1 h2], 'Location of interest', 'Observed earthquake', '(Unknown) fault location');
set(legh, 'Location', 'northeast', 'FontSize', label_size);
set(gca, 'FontSize', axis_size);


set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6.5 2.75]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6.5 2.75]);


print('-dpng', ['Gibbs_observations_w_fault.png']); % save the figure to a file 


%% plot fault estimates

% Case 1
load coordsGibbs1

figure
subplot(1,2,1);
h2 = plot(0, 0, 'ok', 'linewidth', 2);
hold on
for i = 1:size(coordsGibbs,1)
    h1 = plot([coordsGibbs(i,1), coordsGibbs(i,1)+coordsGibbs(i,3)*cos(coordsGibbs(i,4))], ...
              [coordsGibbs(i,2), coordsGibbs(i,2)+coordsGibbs(i,3)*sin(coordsGibbs(i,4))], '-k', 'linewidth', 0.5);
end
h3 = plot(0, 0, 'ok', 'linewidth', 2);
plot(radiusX, radiusY, '-k', 'linewidth', 1)
legh = legend('Location of interest', 'Simulated fault locations');
axis(maxAxis*[-1 1 -1 1])
axis square
hx = xlabel('East-west distance (km)', 'FontSize', label_size);
hy = ylabel('North-south distance (km)', 'FontSize', label_size);
set(gca, 'FontSize', axis_size);


load coordsGibbs2

% Case 2
subplot(1,2,2);
h2 = plot(0, 0, 'ok', 'linewidth', 2);
hold on
for i = 1:size(coordsGibbs,1)
    h1 = plot([coordsGibbs(i,1), coordsGibbs(i,1)+coordsGibbs(i,3)*cos(coordsGibbs(i,4))], ...
              [coordsGibbs(i,2), coordsGibbs(i,2)+coordsGibbs(i,3)*sin(coordsGibbs(i,4))], '-k', 'linewidth', 0.5);
end
h3 = plot(0, 0, 'ok', 'linewidth', 2);
plot(radiusX, radiusY, '-k', 'linewidth', 1)
legh = legend('Location of interest', 'Simulated fault locations');
axis(maxAxis*[-1 1 -1 1])
axis square
hx = xlabel('East-west distance (km)', 'FontSize', label_size);
hy = ylabel('North-south distance (km)', 'FontSize', label_size);
set(gca, 'FontSize', axis_size);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6.5 2.75]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6.5 2.75]);


print('-dpng', ['Gibbs_fault_sims.png']); % save the figure to a file 



