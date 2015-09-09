function [ coords ] = fn_Gibbs_sample(coords, eqX, eqY, fault, numSims)

% sample fault realizations from posterior, using Gibbs sampler

% coords = [x1 y1 length angle p(fault)]

numDisc = 30; % number of points to discretize the fault into


coordValues{1} = [-fault.areaRadius:1:fault.areaRadius]'; % column vector of candidate parameter values 
coordValues{2} = [-fault.areaRadius:1:fault.areaRadius]'; % column vector of candidate parameter values 
coordValues{3} = [3:0.5:13]'; % column vector of candidate parameter values 
coordValues{4} = [pi/6:pi/24:pi/2]'; % column vector of candidate parameter values 

prior{1} = ones(size(coordValues{1})); % uniform prior distribution
prior{2} = ones(size(coordValues{2})); % uniform prior distribution
 prior{3} = ones(size(coordValues{3}));     % uniform distribution
% prior{3} = coordValues{3}.^-1;     % power law distribution
% prior{3} = exp(-coordValues{3}/2);     % exponential prior distribution
prior{4} = ones(size(coordValues{4})); % uniform prior distribution



numValues = [length(coordValues{1}) length(coordValues{2}) length(coordValues{3}) length(coordValues{4})];

for i = 2:numSims
    
    coords(i,:) = coords(i-1,:); % copy over previous parameter values as initial values
    for j = 1:4 % for each component of coords
        distances = zeros(numValues(j), length(eqX)); % initialize
        clstDist = zeros(numValues(j), length(eqX)); % initialize
        log_fx = zeros(numValues(j), length(eqX)); % initialize
        
        
        tempCoords = repmat(coords(i,:),numValues(j),1); 
        tempCoords(:,j) = coordValues{j}; % fill in candidates for new parameter
        
        faultXdisc = linspaceNDim(tempCoords(:,1), tempCoords(:,1)+tempCoords(:,3).*cos(tempCoords(:,4)), numDisc);
        faultYdisc = linspaceNDim(tempCoords(:,2), tempCoords(:,2)+tempCoords(:,3).*sin(tempCoords(:,4)), numDisc);
        
        % find conditional posterior distribution for component j
        for k=1:length(eqX)
            
            % check distance to each discretized point on teh fault
            distances = sqrt((faultXdisc-eqX(k)).^2 + (faultYdisc-eqY(k)).^2);
             
%             log_fx(:,k) = log(sum(normpdf(distances, 0, fault.locUnc),2)); % Base likelihood on average distance to candidate fault
            log_fx(:,k) = log(max(normpdf(distances', 0, fault.locUnc))'); % Base likelihood on closest distance to candidate fault (reasonable, but leads to long faults)
%             log_fx(:,k) = log(max(normpdf(distances', 0, fault.locUnc))'./ (tempCoords(:,3)).^0.3 ); % base likelihood on closest distance, but penalize long faults slightly
            
            
        end
        fx = exp(sum(log_fx,2)) .* prior{j}; % posterior distribution
        %fxFinal = fxFinal ./ sum(fxFinal); % normalize
        Fx = cumsum(fx);
        FxNorm = Fx ./ Fx(end);
        
        % sample a new value from the posterior distribution
        U = rand; % generate a uniform random number
        [~,idx] = histc(U,[0; FxNorm] ); % sample from the CDF using the inverse method 
        coords(i,j) = coordValues{j}(idx);  % store the new sampled value
    end
    
    
end


