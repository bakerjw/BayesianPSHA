function [ p_PGA_U_M ] = fn_Gibbs_PSHA( input_filename, x, T )
% function to compute hazard for a Gibbs-sampled set of faults


load(input_filename)

% hard-coded parameters
Depth = 3; % depth of earthquakes (km)
GMPE = 1; % 1=Atkinson (2015)
M_min = 3.0;
M_max = 7; % maximum possible range--actual faults will be less than this
delta_M = 0.1;
GR_b = 1;
M = M_min:delta_M:M_max; % potential magnitude values
% magnitude probabilities
for i = 1:(length(M)-1)
    p_M_initial(i) = ( (1-10^(-GR_b*(M(i+1) - M_min))) - (1-10^(-GR_b*(M(i) - M_min))) ) / (1-10^(-GR_b*(M_max - M_min)));
end
p_M_initial(length(M)) =  0;

% initialize output matrix
numSims = size(coordsGibbs,1);
p_PGA_U_M = zeros(numSims,length(x));

% discretized fault geometries
numDisc = 10;
faultXdisc = linspaceNDim(coordsGibbs(:,1), coordsGibbs(:,1)+coordsGibbs(:,3).*cos(coordsGibbs(:,4)), numDisc);
faultYdisc = linspaceNDim(coordsGibbs(:,2), coordsGibbs(:,2)+coordsGibbs(:,3).*sin(coordsGibbs(:,4)), numDisc);

    
% PSHA calc
for i = 1:numSims
    
    % adjust magnitude probabilities based on fault length
    p_M = p_M_initial;
    p_M(M > 4.33+1.49*log10(coordsGibbs(i,3))) = 0; % don't allow magnitudes bigger than implied by fault length (from Wells and Coppersmith 1994, table 2a for RLD, SS)
    p_M = p_M./sum(p_M); % renormalize probabilities after truncation
    
    % compute closest distance
    Dist = sqrt(min((faultXdisc(i,:)).^2 + (faultYdisc(i,:)).^2));
    R(i) = sqrt(Dist.^2 + Depth.^2);

    for j = 1:length(x)
        % evaluate GMPE
        switch GMPE
            case 1
                [median, sigma] = A_2015_small_M(M, T, R(i));
            case 2
                for k = 1:length(M)
                    [median(k), sigma(k)] = AB_2006_ENA(M(k), T, R(i), 140, 760);
                end
            case 3
                for k = 1:length(M)
                    [median(k), sigma(k)] = BSSA_2014_nga(M(k), T, R(i), 0, 0, 1, 760);
                end
        end
        p_given_M = 1 - normcdf(log(x(j)),log(median),sigma);
        p_PGA_U_M(i,j) = p_PGA_U_M(i,j) + sum(p_M .* p_given_M); % add contribution of this R value to running total
    end
end

fprintf(['Flt length = ' num2str(mean(coordsGibbs(:,3)),3) ', dist = ' num2str(mean(R),2) ' \n'])

end

