function [expVect, scale] = scaleExponentials(logVect)
%scaleExponentials scales the vector logVect such that when each term in
%the vector is exponentiated and summed, the sum is less than infinity i.e.
%the sum is within the range of floating point calculations in MATLAB.
%   The unscaled value can be obtained by expVect*exp(scale)

maxScale = 1e250; % The number to which the maximum of the exp(logVect - mean) is scaled.
scaleRatio = 1e10;  % Ratio by which the maxScale is reduced if overflow is not eliminated

% Subtract mean of logs to get first scale
scaleMean = mean(logVect);
logVect = logVect - scaleMean;

% Scale the largest value such that the sum is not infinity.
maxLog = max(logVect);
scale = 0;
while true
    expVect = exp(logVect - scale);
    if ~isinf(sum(expVect))
        break;
    else
        scale = maxLog - log(maxScale);
        maxScale = maxScale/scaleRatio;
    end
end

scale = scale + scaleMean;