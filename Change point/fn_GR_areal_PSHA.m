function [lambda_x, M_deagg, M] = fn_GR_areal_PSHA(lambda, R_max, Depth, M_min, M_max, x, delta_M, T, GMPE)

% Created by Jack Baker
% Modified 7/5/2015 to make lambda referenced to M_min
%
% Perform PSHA calculation for a site with an areal source of radius R_max 
% and a Gutenberg-Richter magnitude distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters
%
% lambda        rate of M>M_min earthquakes
%
% R_max         radius of areal source (km)
%
% Depth         depth of areal source (km)
%
% M_min         minimum magnitude
%
% M_min         maximum magnitude
%
% x             Vector of IM values of interest
%
% delta_M       interval for magnitude discretization 
%
% T             period of interest for the spectral acceleration
%               calculation
% 
% GMPE          = 1 for A_2015, =2 for AB_2006, =3 for BSSA NGAW2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output parameters
%
% lambda_x      vector of annual rates of exceedance associated with each
%               value of x
%
% M_deagg       matrix of deaggregation values associate with each
%               magnitude and amplitude considered
%
% M             vector of considered magnitude values (corresponding to
%               deaggregation results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create vector of distances and associated probabilities
numDistances = 30; % number of distance values to consider
R = linspace(0, R_max, numDistances); % equally spaced distance values
if GMPE == 1 % convert R_epi to R_hypo
    R = sqrt(R.^2 + Depth.^2);
elseif GMPE == 2 % convert R_jb to R_clst
    h_eff = 9; % effective depth (km)
    R = sqrt(R.^2 + h_eff.^2);
end
F_R = R.^2 ./ R_max.^2; % CDF
p_R = [0 diff(F_R)]; % use discrete differences in CDF to approximate the PDF

% create vector of magnitudes and associated probabilities
GR_b = 1;
GR_a = log10(lambda) + GR_b * M_min; % find a value implied by this lambda
M = M_min:delta_M:M_max;
for i = 1:(length(M)-1)
    p_M(i) = ( (1-10^(-GR_b*(M(i+1) - M_min))) - (1-10^(-GR_b*(M(i) - M_min))) ) / (1-10^(-GR_b*(M_max - M_min)));
end
p_M(length(M)) =  0;

% p(exceeding each x threshold value)
p_PGA_U_M = zeros(1,length(x));
for j = 1:length(x)
    for i = 1:length(R)
        
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
        p_given_MR = 1 - normcdf(log(x(j)),log(median),sigma);
        p_PGA_U_M(j) = p_PGA_U_M(j) + sum(p_M .* p_R(i) .* p_given_MR); % add contribution of this R value to running total
    end
end

lambda_x = lambda * p_PGA_U_M; % multiply probabilities by the occurence rate of earthquakes

% % deaggregation (TODO)
% for j = 1:length(x)
%     for i = 1:length(R)
% 
%         [median, sigma] = A_2015_small_M(M, T, R);
%         p_given_M = 1 - normcdf(log(x(j)),log(median),sigma);
%         
%         M_deagg(j,:) = (p_given_M .* p_M) / p_PGA_U_M(j);
%     end
% end




