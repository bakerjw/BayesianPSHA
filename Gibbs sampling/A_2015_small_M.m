function [median, sigma] = A_2015_small_M(M, T, R_hyp, h_eff)
%
% coded by Jack Baker, 2/26/2015
%
% updated 3/3/2015 to correct an error in the original electronically
% posted version of the model
%
% updated 4/14/2015 to change default h_eff from Yenier and Atkinson (2014)
% to alternate h_eff proposed in Atkinson (2015) based on personal 
% communication with Atkinson.
%
% Atkinson 2015 GMPE, as defined in the following document. Note that 
% output standard deviations are converted to natural log values from
% the base-10 values reported in the paper.
%
% Atkinson, G. M. (2015). “Ground-Motion Prediction Equation for 
% Small-to-Moderate Events at Short Hypocentral Distances, with Application 
% to Induced-Seismicity Hazards.” Bulletin of the Seismological Society of 
% America.
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Variables
%
% M                     = Moment Magnitude
%
% T                     = Period (sec); Use Period = -1 for PGV 
%
% R_hyp                 = hypocentral distance (km)
% 
% h_eff                 = effective depth (optional). Will be estimated if
%                         not provided by the user
%
% Output Variables
%
% median                = Median amplitude prediction (units of g)
%
% sigma                 = NATURAL LOG standard deviation 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Coefficients
% period = [          -1        0       0.03	0.05	0.1     0.2     0.3     0.5     1       2       3       5       ];
% c0	 = [            -4.198	-2.427	-2.313	-2.337	-2.839	-3.918	-2.076	-4.128	-2.009	-4.503	-3.869	-4.374	];
% c1	 = [            1.818    1.877	 1.84	1.902	1.905	2.112	1.889	1.792	1.89	1.532	1.11	1.134	];
% c2	 = [            -0.1009	-0.1214	-0.1119	-0.1252	-0.1134	-0.1266	-0.1257	-0.0791	-0.1248	-0.043	0.0039	0.0038	];
% c3	 = [            -1.721	-1.806	-1.708	-1.838	-1.658	-1.591	-1.886	-1.526	-1.828	-1.404	-1.447	-1.426	];
% c4   = [           -0.0006  -0.002  -0.002  -0.002  -0.002  -0.0014	-0.00105 -0.0006 0       0       0       0      ];
% sigmaIntra	 = [	0.28	0.29	0.29	0.29	0.3     0.31	0.31	0.3     0.27	0.25	0.25	0.26	];
% sigmaInter	 = [	0.18	0.24	0.26	0.29	0.25    0.2     0.18	0.19	0.21	0.22	0.21	0.17	];
% sigmaTotal	 = [	0.33	0.37	0.39	0.41	0.39    0.37	0.36	0.35	0.34	0.33	0.33	0.31	];

period	= [ 	-1	0	0.03	0.05	0.1	0.2	0.3	0.5	1	2	3	5	];
c0	= [ 	-4.151	-2.376	-2.283	-2.018	-1.954	-2.266	-2.794	-3.873	-4.081	-4.462	-3.827	-4.321	];
c1	= [ 	1.762	1.818	1.842	1.826	1.83	1.785	1.852	2.06	1.742	1.485	1.06	1.08	];
c2	= [ 	-0.09509	-0.1153	-0.1189	-0.1192	-0.1185	-0.1061	-0.1078	-0.1212	-0.07381	-0.03815	0.009086	0.009376	];
c3	= [ 	-1.669	-1.752	-1.785	-1.831	-1.774	-1.657	-1.608	-1.544	-1.481	-1.361	-1.398	-1.378	];
c4	= [ 	-0.0006	-0.002	-0.002	-0.002	-0.002	-0.0014	-0.001	-0.0006	0	0	0	0	];
sigmaIntra	= [ 	0.27	0.28	0.28	0.28	0.29	0.3	0.3	0.29	0.26	0.24	0.24	0.25	];
sigmaInter	= [ 	0.19	0.24	0.27	0.3	0.25	0.21	0.19	0.2	0.22	0.23	0.22	0.18	];
sigmaTotal	= [ 	0.33	0.37	0.39	0.41	0.39	0.37	0.36	0.35	0.34	0.33	0.32	0.31	];

if nargin < 4 % use estimated effective depth
    h_eff = max(1, 10.^(-0.28+0.19.*M)); % effective depth based on alternative method proposed in Atkinson (2015)
    
    % h_eff = max(1, 10.^(-1.72+0.43.*M)); % alternate effective depth from Yenier and Atkinson (2014) -- not preferred by Atkinson
end


% interpolate between periods if neccesary    
if (isempty(find(period == T, 1))) % if the period of interest doesn't match one of the provided periods
    index_low = sum(period<T);
    T_low = period(index_low);
    T_high = period(index_low+1);
    
    % get predictions at closest periods to the period of interest
    [median_low, sigma_low] = A_2015_small_M(M, T_low, R_hyp, h_eff);
    [median_high, sigma_high] = A_2015_small_M(M, T_high, R_hyp, h_eff);
    
    % interpolate results
    x = [log(T_low) log(T_high)];
    Y_sa = [log(median_low) log(median_high)];
    Y_sigma = [sigma_low sigma_high];
    sa = exp(interp1(x,Y_sa,log(T)));
    sigma = interp1(x,Y_sigma,log(T));
    
else % compute median and sigma for the given index value
    i = find(period == T);
    
    % effective distance
    R = sqrt(R_hyp.^2 + h_eff.^2);
    
    % Compute median and sigma
    logY = c0(i) + c1(i)*M + c2(i)*M.^2 + c3(i).*log10(R) + c4(i).*R;
    if T>=0
        median = 10.^logY/981; % convert to units of g
    else
        median = 10.^logY;
    end
    sigma = sigmaTotal(i) * log(10) * ones(size(median)); % convert to natural log, and make a vector if needed

end





