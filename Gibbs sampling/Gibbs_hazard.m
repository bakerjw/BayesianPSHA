% compute hazard associated with Gibbs-sampled faults



% plot specs
label_size = 9; % font size for labels
axis_size = 8;  % font size for axis numbers


T = -1; % Peak Ground Velocity

x = logspace(log(0.2), 2, 100); % PGV values to consider (cm/s)
lamdda_adj = 0.1; % number of earthquakes per year on the source

% Case 1
[ p_PGA_U_M ] = fn_Gibbs_PSHA( 'coordsGibbs1.mat', x, T );
lambda_x1 = lamdda_adj * p_PGA_U_M; % multiply probabilities by the occurence rate of earthquakes
load('coordsGibbs1.mat')

% Case 2
[ p_PGA_U_M ] = fn_Gibbs_PSHA( 'coordsGibbs2.mat', x, T );
lambda_x2 = lamdda_adj * p_PGA_U_M; % multiply probabilities by the occurence rate of earthquakes
load('coordsGibbs2.mat')

% Real fault geometry
[ p_PGA_U_M ] = fn_Gibbs_PSHA( 'faultData.mat', x, T );
lambda_xReal = lamdda_adj * p_PGA_U_M; % multiply probabilities by the occurence rate of earthquakes


%% plot result
axisLimX = [0.3 100];
axisLimY = [0.999e-5 5e-1];
axisLim = [axisLimX axisLimY];

% case 1 
figure
subplot(1,2,1);
h1 = loglog(x, lambda_x1,'-g', 'linewidth', 0.5);
hold on
h3 = loglog(x, mean(lambda_x1),'-k', 'linewidth', 1.25);
hx = xlabel('Peak Ground Velocity (cm/s)', 'FontSize', label_size);
hy = ylabel('Annual rate of exceedance', 'FontSize', label_size);
legh = legend([h1(1), h3], 'Individual curves, Case 1', 'Mean hazard, Case 1');
set(legh, 'Location', 'northeast', 'FontSize', label_size);
axis(axisLim)
set(gca, 'xtick', [0.3 1 10 100])
set(gca, 'xticklabel', [0.3 1 10 100])
set(gca, 'FontSize', axis_size);

% case 2
subplot(1,2,2);
h2 = loglog(x, lambda_x2,'-c', 'linewidth', 0.5);
hold on
h3 = loglog(x, mean(lambda_x1),'-k', 'linewidth', 1.25);
h4 = loglog(x, mean(lambda_x2),'--b', 'linewidth', 1.25);
h5 = loglog(x, lambda_xReal,'-.k', 'linewidth', 1.25);
hx = xlabel('Peak Ground Velocity (cm/s)', 'FontSize', label_size);
hy = ylabel('Annual rate of exceedance', 'FontSize', label_size);
legh = legend([h2(1), h4, h3 h5], 'Individual curves, Case 2', 'Mean hazard, Case 2', 'Mean hazard, Case 1', 'Real fault geometry');
set(legh, 'Location', 'northeast', 'FontSize', label_size);
axis(axisLim)
set(gca, 'xtick', [0.3 1 10 100])
set(gca, 'xticklabel', [0.3 1 10 100])
set(gca, 'FontSize', axis_size);

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6.5 2.75]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6.5 2.75]);


print('-dpng', ['Gibbs_haz.png']); 



