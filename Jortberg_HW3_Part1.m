%Homework 3 Draft EECE 5644
%load datasetP2Train1.mat
load datasetP2Train2.mat
load datasetP2Test.mat
%2.1 Parameter Estimation
%a) Estimate mu and cov for datasetP2Train1 (ML Gaussian)

c = length(unique(classIndex));
d = size(data,2);
class = {};
mu = {};
sigma = {};
n={};
for h= 1:c
    idx = find(classIndex==h);
    class{h} = data(idx, :); 
    mu{h} = [mean(class{h}(:,1)); mean(class{h}(:,2)); mean(class{h}(:,3))]; %d long
    sigma{h} = cov(class{h}); 
    n{h} = length(class{h});
end
 
%part b, use shrinkage equation for Sigma.  Plot training error as function
%of alpha
ncom = length(data); %should n be number of features, 3 (instead of 20/10000)?
cov_com = cov(data); 
prior = [(1/3); (1/3); (1/3)]; %should I set priors to be roughly equal?
alpha = [0.01:0.01:0.99];
accuracy_train = zeros(1, length(alpha));
accuracy_test = zeros(1, length(alpha));
for i=1:length(alpha)
    for h=1:c;
        sigma_shrink{h} = ((1-alpha(i))*n{h}*sigma{h} + alpha(i)*ncom*cov_com)/((1-alpha(i))*n{h} + alpha(i)*ncom);
    end
    [score] = gaussianDiscriminantAnalysis(data, mu, sigma_shrink, prior);
    [~,classIndexEstimate] = max(score, [], 2);
    accuracy_train(i) = 100*nnz(classIndexEstimate == classIndex)/length(classIndexEstimate);
    
    [score] = gaussianDiscriminantAnalysis(dataTest, mu, sigma_shrink, prior);
    [~,classIndexEstimate] = max(score, [], 2);
    accuracy_test(i) = 100*nnz(classIndexEstimate == classIndexTest)/length(classIndexEstimate);
end

figure; plot(alpha, 1-(accuracy_train/100)); xlabel('alpha'); ylabel('training error'); title('Training Set Classification Error')
figure; plot(alpha, 1-(accuracy_test/100)); xlabel('alpha'); ylabel('test error'); title('Test Set Classification Error')

idx_op = find(min(1-(accuracy_test/100)) == 1-(accuracy_test/100));
alpha_op = alpha(idx_op);
%Part c)
    %Plot using test data is much smoother (not step like) since there are more
    %data points.
    %I'd estimate that the optimal alpha is 0.13
%Part f)
    %Now, with nsamples in training and test set equal, the classification
    %error plots look similar
    %Optimal alpha is 0.03
    %From c to f, the training and testing error plot became more similar
    %as nsamples in training = nsamples in testing
    %Shrinkage is useful because when n is small, ML for sigma can be
    %unstable. Shrinkage helps minimize estimation error.
    
 

