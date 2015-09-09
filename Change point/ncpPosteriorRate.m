function [probVect] = ncpPosteriorRate(interEventT, rateVect, params) 
% calculate the posterior probability distribution 
% of the, event rate, assuming a Gamma distribution prior and a Poisson
% likelihood function. 
%
% Created by Jack Baker
% last modified 9/9/2015


a = params(1, 1);
b = params(1, 2);

%% Calculate probability of event rate before the change point occured
totDays = sum(interEventT) + 1;    % 1 is added to account for startDate
totEvents = length(interEventT) + 1;   % 1 is added to account for first event on day 0

probVect = zeros(size(rateVect));
probVectUnscaled = zeros(size(rateVect));
logScaleVect = zeros(size(rateVect));   % This vector stores the mean of log values at each rate
for i = 1:length(rateVect)
    probLogVect = (totEvents+a-1)*log(rateVect(i)) - rateVect(i)*(totDays+1/b);
    probVectUnscaled(i) = exp(probLogVect);
end


%% Scale probVect to obtain pdf
integrand = trapz(rateVect, probVectUnscaled);
probVect = probVectUnscaled./integrand;  %Normalize probability vector
