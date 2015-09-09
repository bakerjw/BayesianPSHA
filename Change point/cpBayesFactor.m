function [bayesFactor] = cpBayesFactor(dataVect, params) 
%  calculate the Bayes factor for no change versus a
%  change point using Raftery and Akman (1986) approach.
%
%   Raftery, A. E., and Akman, V. E. (1986). "Bayesian analysis of a 
%   Poisson process with a change-point." Biometrika, 73(1), 85?89.
%
%   ChangePointBayesFactor calculates the Bayes Factor of the
%   change point over the time range specified in dataVect. dataVect is the
%   vector with the time between each successive event. params is a 2x2 
%   matrix for parameters of the change model prior rates 
%   where each row defines the shape and scale of the event rate
%   before and after the change point in a gamma distribution. The no-change
%   model utilizes only the first rowl of params. The event
%   rates are assumed with gamma distribution priors and the change point
%   is assumed with a continuous uniform distribtuion over the total
%   duration of data. A Bayes Factor less than one implies a preference for
%   the change point model.
%
% Created by Abhineet Gupta
% Modified by Jack Baker, last modified 9/9/2015

a1 = params(1, 1);
b1 = params(1, 2);
a2 = params(2, 1);
b2 = params(2,2);

a0 = params(1,1);
b0 = params(1,2);


% Use Eq 3.2 in Raftery Akman when certian conditions are met
if isinf(b0) && isinf(b1) && isinf(b2) && a0 == 0.5 && a1 == 0.5 && a2 == 0.5
    % initialize date and event vectors 
    totDays = sum(dataVect) + 1;    % 1 is added to account for startDate
    [numEventsVect, totEvents] = getNumEventsVect(dataVect, totDays);
    
    % numerater of Bayes Factor
    numLog = gammaln(totEvents + a0) - gammaln(a0) - (totEvents + a0)*log(1/b0 + totDays);
    
    % denominator of Bayes Factor
    tau = [1:totDays-1]';
    numEvents = numEventsVect;
    r1 = numEvents + a1;
    r2 = totEvents - numEvents + a2;
    s1 = tau + 1/b1;
    s2 = totDays - tau + 1/b2;
    denVect = -log(totDays) + gammaln(r1) + gammaln(r2) - r1.*log(s1) - r2.*log(s2);
    [denVect, scale] = scaleExponentials(denVect);   % Exponentiate probability vector to account for very large or very small values.
    denLog = log(sum(denVect)) + scale - gammaln(a1) - gammaln(a2);
    
    bayesFactor = exp(numLog - denLog);
    bayesFactor = bayesFactor*((4*pi)^0.5)*(totDays^-0.5);
else
    warning(['The constant in Bayes factor calculation is evaluated only for', char(10),...
        'the case of gamma priors of 0.5 and infinity. For other prior values,' char (10),...
        'Bayes factor may only be used for comparison between datasets.']);
end

