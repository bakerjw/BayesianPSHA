% This script calculates performs change point calculations given a set of
% observed events, computer posterior rate distributions, and evaluates
% implications for seismic hazard

% modified from code by Abhineet Gupta, 
% last modified by Jack Baker, 9/9/2015

% plot specs
label_size = 9; % font size for labels
axis_size = 8;  % font size for axis numbers


%% specify occurence times of earthquakes
startYear = 1974; % start date of catalog
interEventT = [649	1848	10353	173	106	34	29	2	52	510	48	10	12	217	20	197	42	22	87	45	10	87	27	93	12	36	10	45	49	27	30	14	40	68	54	30	44]; % number of days between earthquakes

eventIndices = [3 6 length(interEventT)]; % which dates to highlight in the analysis below


%% Calculate posterior event rates 

params = [0.5, inf; 0.5, inf]; % Prior Gamma distribution parameters for change point model
rateVect = logspace(-5, -1, 100); % occurence rates to compute posterior density at
for i = 1:length(interEventT) % consider the event occurences up to event i
    
    BF(i) = cpBayesFactor(interEventT(1:i), params); % compute Bayes factor
    probVectCP(i,:) = cpPosteriorRate(interEventT(1:i), rateVect, params); % posterior PDF of current rate, if there is a change
    probVectNoCP(i,:) = ncpPosteriorRate(interEventT(1:i), rateVect, params); % posterior PDF, if there is no change
    meanRateCP(i) = trapz(rateVect, rateVect*365 .* probVectCP(i,:)); % mean annual rate (multiplying by 365 to get an annual rate)

    % use the analytical solution for the mean rate in the no-change-point case
    aPost(i) = params(1,1) + i; % posterior parameters
    bPost(i) = 99999 /(sum(interEventT(1:i))*99999+1); % substitute a finite b parameter for numerical convenience
end
meanRateNoCP = aPost .* bPost * 365; % multiply by 365 to get an annual rate
% save('OKCPostRates.mat', 'params', 'rateVect', 'BF', 'probVectCP', 'meanRateCP', 'probVectNoCP', 'meanRateNoCP') % save results

% load OKCPostRates % load previous results (to save time in re-computing)


% use Bayes' Factors to get a composite posterior result (using Change Point rates only when Bayes factor suggests so)
changeFlag = (BF < 0.01);
meanComposite = meanRateCP.*changeFlag + meanRateNoCP.*(1-changeFlag);
probVectComposite = probVectCP .* (changeFlag' * ones(size(rateVect))) + probVectNoCP .* ((1-changeFlag)' * ones(size(rateVect)));

% some supplemental calculations
timesInYrs = cumsum(interEventT)/365 + startYear; % cumulatively sum events, convert to years, and set change time to t=0
numEvents = length(timesInYrs);
freqRate = (1:numEvents) ./ (timesInYrs - startYear); % pure frequentist event rates
exampleT = timesInYrs(eventIndices);


%% plot cumulative counts

xLimVals = [startYear 2017]; % range of years to plot

figure
subplot(2,2,1);
h1 = stairs([0 timesInYrs], [0:numEvents], '-k', 'Linewidth', 2);
hold on
h2 = plot(exampleT(1), eventIndices(1), 'og', 'Linewidth', 2);
h3 = plot(exampleT(2), eventIndices(2), 'ob', 'Linewidth', 2);
h4 = plot(exampleT(3), eventIndices(3), 'or', 'Linewidth', 2);
hx = xlabel('Year', 'FontSize', label_size);
hy = ylabel('Number of M>3 earthquakes', 'FontSize', label_size);
legh = legend([h1 h2 h3 h4], 'Earthquake count', 'Time A', 'Time B', 'Time C');
set(legh, 'Location', 'northwest', 'FontSize', label_size);
axis([xLimVals 0 40])
set(gca, 'xtick', 1975:10:2015)
set(gca, 'FontSize', axis_size);


%% earthquake rates at three points in time

subplot(2,2,2);
h1 = semilogx(rateVect*365, probVectComposite(eventIndices(1),:)./max(probVectComposite(eventIndices(1),:)), '-g', 'Linewidth', 2);
hold on
h2 = semilogx(rateVect*365, probVectComposite(eventIndices(2),:)./max(probVectComposite(eventIndices(2),:)), '-b', 'Linewidth', 2);
h3 = semilogx(rateVect*365, probVectComposite(eventIndices(3),:)./max(probVectComposite(eventIndices(3),:)), '-r', 'Linewidth', 2);
plot([1 1]*meanComposite(eventIndices(1)), [0 0.15], '-g', 'Linewidth', 2);
plot([1 1]*meanComposite(eventIndices(2)), [0 0.15], '-b', 'Linewidth', 2);
plot([1 1]*meanComposite(eventIndices(3)), [0 0.15], '-r', 'Linewidth', 2);
legh = legend([h1 h2 h3], 'Time A', 'Time B', 'Time C');
set(legh, 'Location','northwest', 'FontSize', label_size)
hx = xlabel('Annual rate of M>3 earthquakes', 'FontSize', label_size);
hy = ylabel('Probability density', 'FontSize', label_size);
set(gca, 'xlim', [1e-2 2e1])
set(gca, 'xtick', [ 1e-2 1e-1 1e0 1e1])
set(gca, 'FontSize', axis_size);


%% mean rates over time

yLimVals = [4e-2 2e1]; % range of y axis values to plot

subplot(2,2,3);
h8 = semilogy(timesInYrs, freqRate, '-c', 'Linewidth', 2);
hold on
hCP = semilogy(timesInYrs, meanRateCP, '-c', 'Linewidth', 1);
hNCP = semilogy(timesInYrs, meanRateNoCP, '-k', 'Linewidth', 1);
h1 = semilogy(timesInYrs, meanComposite, '-k', 'Linewidth', 2);
h2 = semilogy(exampleT(1), meanComposite(eventIndices(1)), 'og', 'Linewidth', 2);
h3 = plot(exampleT(2), meanComposite(eventIndices(2)), 'ob', 'Linewidth', 2);
h4 = plot(exampleT(3), meanComposite(eventIndices(3)), 'or', 'Linewidth', 2);
set(gca, 'ylim', yLimVals)
set(gca, 'xlim', xLimVals)
legh = legend([h1 h8 ], 'Mean posterior rate', 'Frequentist rate');
set(legh, 'Location','northwest', 'FontSize', label_size)
hy = ylabel('Annual rate of M>3 earthquakes', 'FontSize', label_size);
hx = xlabel('Year', 'FontSize', label_size);
set(gca, 'FontSize', axis_size);
set(gca, 'xtick', 1975:10:2015)



%% PSHA calculations

x = logspace(log(0.2), 2, 100); % PGV values to consider (cm/s)
axisLimX = [0.3 100]; % figure axis limits
axisLimY = [0.999e-5 5e-1]; % figure axis limits
T = -1; % Peak Ground Velocity
GMPE = 1; % 1=Atkinson (2015)
M_min = 3.0;
M_max = 7;
lambda = 1; % rate of EQ's bigger than M_min
R = 25; % radius of source (km)
Depth = 3; % depth of earthquakes (km)
delta_M = 0.2;
lambda_x = fn_GR_areal_PSHA(lambda, R, Depth, M_min, M_max, x, delta_M, T, GMPE); % compute hazard

% hazard curves
subplot(2,2,4);
loglog(x, lambda_x*meanComposite(eventIndices(1)),'-g', 'linewidth', 2)
hold on
loglog(x, lambda_x*meanComposite(eventIndices(2)),'-b', 'linewidth', 2)
loglog(x, lambda_x*meanComposite(eventIndices(3)),'-r', 'linewidth', 2)
legh = legend('Time A', 'Time B', 'Time C');
set(legh, 'Location','northeast', 'FontSize', label_size)
hx = xlabel('Peak Ground Velocity (cm/s)', 'FontSize', label_size);
hy = ylabel('Annual rate of exceedance', 'FontSize', label_size);
axis([axisLimX axisLimY])
set(gca, 'FontSize', axis_size);
set(gca, 'xtick', [0.3 1 10 100])
set(gca, 'xticklabel', [0.3 1 10 100])

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [6.5 5.5]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6.5 5.5]);

print('-dpng', ['Master_Figure.png']); % save the figure to a file in a format suitable for MS Office files
