function [mu_new, cov_new, likelihood, log_prior] = MaximizationStep(Q, inputData, numClust, mu)
    %maximize log p(y,z|phi)
    
    mu_new = {};
    [n, d] = size(inputData);
    %mu_new will be a cell k long with d values
    norm_Q = zeros(1, numClust);
    for i =1:n
        norm_Q = norm_Q + Q(i, :);
    end
    %1. Estimate new means
    % norm constant * sum over E[z]*y
    for j = 1:numClust
        sum_Q = zeros(1, d);
        for i = 1:n
            sum_Q = sum_Q + Q(i, j)*inputData(i,:);
        end
        mu_new{j} = sum_Q/norm_Q(j);
    end
    
    %2. Estimate new covariances
    %norm constant * sum over E[z]*(y-mu_new)*(y-mu_new)'
    for j = 1:numClust
        covsum = zeros(d, d);
        norm_c = 0;
        for i = 1:n
            norm_c = norm_c + Q(i, j);
            covsum = covsum + Q(i, j)*((inputData(i,:)-mu{j})'*(inputData(i,:)-mu{j}));
        end 
        cov_new{j} = covsum/norm_c;
    end
    
    %3. Update Likelihood Estimates
    likelihood = zeros(1, numClust);
    for i =1:n
        likelihood = likelihood + Q(i,:);
    end
    likelihood = likelihood/n;
    log_prior = loglikelihood(inputData, mu_new, cov_new, likelihood, numClust);
    
end