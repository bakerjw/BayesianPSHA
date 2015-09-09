function [probVect] = cpPosteriorRate(dataVect, rateVect, params) 
% calculate the probability distribution of the
% event rate after the change point using Raftery and Akman (1986) 
%
%   Raftery, A. E., and Akman, V. E. (1986). "Bayesian analysis of a 
%   Poisson process with a change-point." Biometrika, 73(1), 85?89.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   rateVect is the vector of rates over which probabilities will
%   be caculated.
%
%   changePointRateProbability calculates the probability distribtuion of the
%   pre/post-change point event rate over the time range specified in dataVect. 
%
%   dataVect is the vector with the time between each successive event. 
% 
%   params is a 2x2 matrix where each row defines the shape and scale of 
%   the gamma event rate distribution before and after the change point 
%
%   probVect is the posterior density for lambda_2, evaluated at the rates
%   specified in rateVect
%
% Created by Abhineet Gupta
% Modified by Jack Baker, last modified 9/9/2015

a1 = params(1, 1);
b1 = params(1, 2);
a2 = params(2, 1);
b2 = params(2, 2);

%% Calculate probability density for event rate 
totDays = sum(dataVect)+1;
[numEventsVect, totEvents] = getNumEventsVect(dataVect, totDays);

probVect = zeros(size(rateVect));
probVectUnscaled = zeros(size(rateVect));
logScaleVect = zeros(size(rateVect));   % This vector stores the mean of log values at each rate

% tau = datenum(datesVect) - (datenum(startDate) * ones(size(datesVect)));
tau = [1:totDays-1]';
r1 = numEventsVect + a1;
r2 = totEvents - numEventsVect + a2;
s1 = tau + 1/b1;
s2 = totDays - tau + 1/b2;

for i = 1:length(rateVect)
    probLogVect = -log(totDays) + (r2 - 1).*log(rateVect(i)) - rateVect(i).*s2 + gammaln(r1) - r1.*log(s1);
    [probExpVect, logScaleVect(i)] = scaleExponentials(probLogVect);    % Exponentiate probability vector to account for very large or very small values.
    probVectUnscaled(i) = sum(probExpVect);
end

%% Scale unscaled probability vector such that each value is scaled to the
% same amount, and can be used to calculate probabilities

% The probVect obtained up to this point is a represenatation of the pdf of
% the posterior distribution which is unscaled. To scale it, we make an
% approximation by integrating this pdf at the discrete values provided in
% rateVect, and making the resulting integrand equal to 1 by proper
% scaling.

scale = mean(logScaleVect);
scaleRatio = 10;    %Ratio by which scale is reduced if sum of scaled probabilities is infinity
while true
    probVect = probVectUnscaled.*exp(logScaleVect - scale);
    if ~isinf(trapz(rateVect, probVect))
        break;
    else
        scale = scale + log(scaleRatio);
    end
end

%% Scale probVect to obtain pdf
integrand = trapz(rateVect, probVect);
probVect = probVect./integrand;  %Normalize probability vector
