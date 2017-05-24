function [mu_out, cov_out, prior_out, estimatedLabels, logLikelihood, costvComplex] = ExpectationMaximization(inputData, numClust, stopTolerance, numRuns)
%Inputs: 
% data: nSamples by d array with the data distributed along the rows
% numClust: k number of clusters to be identified
% stopTol: stop tolerance for convergence
% numRuns: number of runs in which algorithm is reinitialized
%
%Outputs:
% mu: means of clusters
% cov: covariances of clusters
% prior: prior probabilities

%Note: save outputs as cells within cells per run, so at end can select
%best logLikelihood option as final output

%n points in d dimensional space
    [n, d] = size(inputData);

    for i=1:numRuns
        log_prior_list = [];
        maxIter=50;
        mu = {};
        cov = {};
        prior = zeros(1, numClust);
        prior(:) = 1/numClust;
        likelihood = prior;
        %Step 1: Randomly initialize mu, cov, prior
        centroids = inputData(randsample(n, numClust), :);
        for j=1:numClust
            mu{j} = centroids(j,:);
            cov{j} = eye(d);
        end
        %Step 2: Gaussian Mixture Model
        %Step 2a. Compute Initial Log Likelihood
        log_prior = loglikelihood(inputData, mu, cov, prior, numClust);
        likelihoodNew = log_prior;
        %create list of log_priors across iterations, pick best at end
        log_prior_list=[log_prior_list log_prior];
        
        cntr=0;
        mu_new = mu;
        cov_new = cov;
        while cntr < maxIter
            cntr= cntr+1;
            likelihoodOld = likelihoodNew;
            %Step 2b. Expectation (E) Step
            Q = ExpectationStep(inputData, likelihood, mu_new, cov_new, numClust);
            
            mu_old = mu_new;
            %Step 2c. Maximization (M) Step
            [mu_new,cov_new,likelihood,log_prior]=MaximizationStep(Q,inputData,numClust, mu_old);
            likelihoodNew = log_prior;
            log_prior_list=[log_prior_list log_prior];
            
            %Loop Exit Mechanism - if reach stopTolerance
            thresh = abs((likelihoodNew - likelihoodOld)/likelihoodOld);
            if thresh < stopTolerance
                break
            end
            
        end %end while loop
        
    %Here need to track mu, covariance, Q, etc per run and pick BEST for final output
    mu_per_run{i} = mu_new;
    cov_per_run{i} = cov_new;
    log_prior_per_run{i} = log_prior_list;
    %figure; plot(log_prior_list);
    best_log_priors(i) = log_prior_list(end);
    Q_per_run{i} = Q;
    likelihood_list{i} = likelihood;
    
    end %end for number of runs
    
    %pick best option out of all runs for final return
    [maxlogs, idx] = max(best_log_priors); 
    best_Q_estimate = Q_per_run{idx};
    cov_out = cov_per_run{idx};
    mu_out = mu_per_run{idx};
    prior_out = likelihood_list{idx};
    
    logLikelihood = log_prior_per_run{idx}; %1 x iterations long (cntr)
    
    estimatedLabels = zeros(1, numClust);
    for i = 1:n
        label = max(best_Q_estimate(i,:)); 
        estimatedLabels(i) = find(best_Q_estimate(i,:) == label);
    end
    
    m= numClust*(length(mu)+ d + numClust); %number of estimated parameters 
    costvComplex = logLikelihood(end) - (m/2)*log(n); 
end