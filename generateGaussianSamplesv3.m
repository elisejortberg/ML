function [data, classIndex_out] = generateGaussianSamplesv3(mu, sigma, nsamples, prior)
% Function to simulate data from k Gaussian densities (1 for each class) in d dimensions.
%
% INPUTS:
% mu : k by 1 cell with the class dependent d?dimensional mean vector
% sigma: k by 1 cell with the class dependent d?by?d covariance matrix
% nSamples: scalar indicating number of samples to be generated.
% prior: k by 1 vector with class dependent mean, sum of priors must equal
% 1
%
% OUTPUTS:
% data:  nSamples by d array with the simulated data distributed along the rows
% classIndex:  vector of length nSamples with the class index for each datapoint

    classIndex = [];
    
    d = length(mu{1}); %d = len(mu)
    I = eye(d); 

    classIndex = mnrnd(nsamples, prior); %generate how many samples should be in each class
    nclasses = length(prior); %nclasses = k
    
    data= {};
    classIndex_out = {};
    for j=(1:nclasses)
        data{j} = mvnrnd(mu{j}, sigma{j}, classIndex(j)); %generate randomly distributed data
        classIndex_out{j} = j*ones([classIndex(j), 1]); %assign as class 1, 2, 3..,j
    end
    
    %plot(data{1}(:,1), data{1}(:, 2), '.'); hold on, plot(data{2}(:,1), data{2}(:, 2), 'r.');
    
end