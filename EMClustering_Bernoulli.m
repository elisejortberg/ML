function [clusterParameters, dataProbabilities, estimatedLabels, logLikelihood] = EMClustering_Bernoulli(inputData, params)
% Runs EM
% Inputs:
%       inputData (Nxd)             - matrix with data
%       params 
%           .numberOfClusters       - number of clusters
%           .maxNumberOfClusters    - if numberOfClusters is not specified,
%                                     use this to run EM with the number of clusters from 1:maxNumberOfClusters
%           .numberOfRuns           - number of random initializations
%
% Outputs
%       clusterParameters           - struct with Bernoulli Distribution
%                                     parameters
%       dataProbabilities (NxK)     - array with the probabilities of each
%                                     data point for each cluster
%       estimatedLabels (Nx1)       - label based on maximum probability
%       logLikelihood (1xiterN)     - log-likelihood as a function of
%                                     iteration


    if isfield(params, 'numberOfClusters')
        % If number of clusters is specified, run once
        [clusterParameters, dataProbabilities, estimatedLabels, logLikelihood] = runEM_Bernoulli(inputData, params);  
    end
end

function [clusterParameters, dataProbabilities, estimatedLabels, outputLog] = runEM_Bernoulli(inputData, params)

% Initialize cells for random runs
clusterParameters = cell(params.numberOfRuns, 1);
dataProbabilities = cell(params.numberOfRuns, 1);
estimatedLabels = cell(params.numberOfRuns, 1);
outputLog = cell(params.numberOfRuns, 1);
logLikelihood = zeros(params.numberOfRuns, 1);
initClusterParameters = cell(params.numberOfRuns, 1);

% Run EM with random initializations in parallel with parfor
%parfor idxR = 1:params.numberOfRuns
parfor idxR = 1:params.numberOfRuns
    
    fprintf('\tOn run %d\n', idxR);
    
    % Random means from data
    %indexRandomInit = randperm(size(inputData,1));
    b=0.75; %pick initMeans between 0.25 and 0.75- images binary
    a=0.25;

    indexRandomInit = (b-a).*rand(params.numberOfClusters, size(inputData,2))+ a;
    initMeans = indexRandomInit; 
    
    %GMM
    %initClusterParameters{idxR} = struct('prior', num2cell(ones(params.numberOfClusters,1)/params.numberOfClusters), 'mu', mat2cell(initMeans, ones(params.numberOfClusters, 1), size(inputData, 2)), 'covar', eye(size(inputData, 2)));
    %Bernoulli
    initClusterParameters{idxR} = struct('prior', num2cell(ones(params.numberOfClusters,1)/params.numberOfClusters), 'mu', mat2cell(initMeans, ones(params.numberOfClusters, 1), size(inputData, 2))); 
    
    [clusterParameters{idxR}, dataProbabilities{idxR}, estimatedLabels{idxR}, outputLog{idxR}] = EM_Bernoulli(inputData, initClusterParameters{idxR});
    
    logLikelihood(idxR) = outputLog{idxR}(end);
end

% Choose run with best log-likelihood
[~, bestRunIndex] = max( logLikelihood );

clusterParameters = clusterParameters{bestRunIndex};
dataProbabilities = dataProbabilities{bestRunIndex};
estimatedLabels = estimatedLabels{bestRunIndex};
outputLog = outputLog{bestRunIndex};


end

function [clusterParameters, dataProbabilities, estimatedLabels, logLikelihood] = EM_Bernoulli(inputData, initClusterParameters)

% Harcoded tolerance 
tolerance = 0.001;

% Initialize
stopCondition = false;
iter = 0;

Ezij = zeros(size(inputData, 1), numel(initClusterParameters));

clusterParameters = initClusterParameters; 

logLikelihood = [];

while ~stopCondition
    iter = iter + 1;
    
    % E step
    for idxJ = 1:numel(clusterParameters) %mu = 1xd (784), data= nxd
        mu = clusterParameters(idxJ).mu;
        prior = clusterParameters(idxJ).prior;
        Ezij(:, idxJ) = prior*exp(sum(inputData*log((mu+realmin).') + (1-inputData)*log(1-(mu+realmin).'), 2));
    end
    Ezij = bsxfun(@rdivide, Ezij, sum(Ezij,2));
    
    % Get rid of NaNs due to likelihood being too far off. In such case,
    % denominator is 0. With epsilon trick, we realize Ezij has to be 1/K
    Ezij(isnan(Ezij)) = 1/numel(clusterParameters);
    
    % M step
    for idxJ = 1:numel(clusterParameters)
        % Prior (pi)
        clusterParameters(idxJ).prior = mean(Ezij(:,idxJ));
        
        % Mean (mu)
        clusterParameters(idxJ).mu = mean(bsxfun(@times,Ezij(:,idxJ),inputData)) / clusterParameters(idxJ).prior;
        
    end
    
    % Stopping criteria
    sumProbabilities = 0;
    for idxJ = 1:numel(clusterParameters)
        mu = clusterParameters(idxJ).mu;
        prior = clusterParameters(idxJ).prior;
        sumProbabilities = sumProbabilities + prior*exp(sum(inputData*log((mu+realmin).') + (1-inputData)*log(1-(mu+realmin).'), 2));
    end
    
    logLikelihood(iter) = sum(log(sumProbabilities)); 
    
    % Check likelihood increase
    if iter >= 2
        percentChange = 100*(logLikelihood(iter)-logLikelihood(iter-1))/abs(logLikelihood(iter));
        stopCondition = percentChange <= tolerance;
        
        % Discard run when likelihood decreases (happens for small, badly initialized datasets)
        if sign(percentChange) == -1
            clusterParameters = [];
            dataProbabilities = [];
            estimatedLabels = [];
            logLikelihood = -Inf;
            return
        end
        
        % Printing
        if mod(iter,5) == 2
            fprintf('\t\t On iteration %d with percent change = %f\n', iter, percentChange);
        end
    end
    
    if iter > 500 
        break
    end
    
end

% Compute probabilities
dataProbabilities = zeros(size(inputData, 1), numel(initClusterParameters));
for idxJ = 1:numel(clusterParameters)
    mu = clusterParameters(idxJ).mu;
    prior = clusterParameters(idxJ).prior;
    dataProbabilities(:, idxJ) = prior*exp(sum(inputData*log((mu+realmin).') + (1-inputData)*log(1-(mu+realmin).'), 2));
end

% Estimate labels through hard assignment
[~, estimatedLabels] = max(dataProbabilities,[],2);

end