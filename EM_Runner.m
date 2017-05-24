%% Tester Data Set: 4 Separated Classes
% nsamples = 400;
% mu{1} = [-1,-1]; 
% mu{2} = [3,3];
% mu{3} = [8, 8];
% mu{4} = [2, 10];
% I = eye(2); 
% sigma{1} = I; 
% sigma{2} = I;
% sigma{3} = I;
% sigma{4} = I;
% prior = [(1/4), (1/4), (1/4), (1/4)];
% [data, classIndex] = generateGaussianSamplesv3(mu, sigma, nsamples, prior);
% figure(1); scatter(data{1}(:,1), data{1}(:,2)); hold on, scatter(data{2}(:,1), data{2}(:,2));
% hold on, scatter(data{3}(:,1), data{3}(:,2)); scatter(data{4}(:,1), data{4}(:,2))

% inputData = vertcat(data{1}, data{2}, data{3}, data{4});

%% dataset1.mat
load dataset1.mat

classes = unique(labels);
class = {};
for i=1:length(classes)
    class_idx = find(labels==classes(i));
    class{i} = data(class_idx, :);
end
%figure(1); scatter(class{1}(:,1), class{1}(:,2)); hold on, scatter(class{2}(:,1), class{2}(:,2));

numClust = 2;
stopTolerance = 0.00001;
numRuns =10;
BICflag =0;
BIC = [];
if BICflag
    numClustmin = 2;
    numClustmax = 10;
    for i =numClustmin:numClustmax
        disp(i);
        [mu_out, cov_out, prior_out, estimatedLabels, logLikelihood, costvComplex] = ExpectationMaximization(data, i, stopTolerance, numRuns);
        BIC = [BIC costvComplex];
        %for each run, make a cell, then pick max costvComplex 
    end 
    figure; plot(numClustmin:numClustmax, BIC)
else
    [mu_out, cov_out, prior_out, estimatedLabels, logLikelihood, costvComplex] = ExpectationMaximization(data, numClust, stopTolerance, numRuns);
end
    
classes = unique(estimatedLabels);
class_estimated = {};
for i=1:length(classes)
    class_idx = find(estimatedLabels==classes(i));
    class_estimated{i} = data(class_idx, :);
end
%figure(2); scatter(class_estimated{1}(:,1), class_estimated{1}(:,2)); hold on, scatter(class_estimated{2}(:,1), class_estimated{2}(:,2));
%figure(3); plot(logLikelihood);

G = gmdistribution(mu_out{1},cov_out{1})
F = @(x,y) pdf(G,[x y])
G1 = gmdistribution(mu_out{2},cov_out{2})
F1 = @(x,y) pdf(G1,[x y])

figure; 
subplot(3,1,1);
scatter(class{1}(:,1), class{1}(:,2)); hold on, scatter(class{2}(:,1), class{2}(:,2));
title('True Class');
subplot(3,1,2);
scatter(class_estimated{1}(:,1), class_estimated{1}(:,2)); hold on, scatter(class_estimated{2}(:,1), class_estimated{2}(:,2));
ezcontour(F); hold on, ezcontour(F1)
title('Estimated Classes');
subplot(3,1,3);
plot(logLikelihood);
title('Log Likelihood over iterations');

%% dataset2.mat

% load dataset2.mat
% 
% classes = unique(labels);
% class = {};
% for i=1:length(classes)
%     class_idx = find(labels==classes(i));
%     class{i} = data(class_idx, :);
% end
% 
% numClust = 3;
% stopTolerance = 0.00001;
% numRuns =10;
% BIC_flag=1;
% BIC = [];
% if BIC_flag
%     numClustmin = 2;
%     numClustmax = 10;
%     for i =numClustmin:numClustmax
%         disp('cluster'); disp(i);
%         [mu_out, cov_out, prior_out, estimatedLabels, logLikelihood, costvComplex] = ExpectationMaximization(data, i, stopTolerance, numRuns);
%         BIC = [BIC costvComplex];
%     end 
%     figure; plot(numClustmin:numClustmax, BIC); title('dataset2: numClust v BIC');
% else
%     [mu_out, cov_out, prior_out, estimatedLabels, logLikelihood, costvComplex] = ExpectationMaximization(data, numClust, stopTolerance, numRuns)
% end
% 
% classes = unique(estimatedLabels);
% class_estimated = {};
% for i=1:length(classes)
%     class_idx = find(estimatedLabels==classes(i));
%     class_estimated{i} = data(class_idx, :);
% end
% 
% G = gmdistribution(mu_out{1},cov_out{1})
% F = @(x,y) pdf(G,[x y])
% G1 = gmdistribution(mu_out{2},cov_out{2})
% F1 = @(x,y) pdf(G1,[x y])
% G2 = gmdistribution(mu_out{3},cov_out{3})
% F2 = @(x,y) pdf(G2,[x y])
% 
% figure; 
% subplot(3,1,1);
% scatter(class{1}(:,1), class{1}(:,2)); hold on, scatter(class{2}(:,1), class{2}(:,2));
% hold on, scatter(class{3}(:,1), class{3}(:,2));
% title('True Class');
% subplot(3,1,2);
% scatter(class_estimated{1}(:,1), class_estimated{1}(:,2)); hold on, scatter(class_estimated{2}(:,1), class_estimated{2}(:,2));
% scatter(class_estimated{3}(:,1), class_estimated{3}(:,2)); hold on,
% ezcontour(F); hold on, ezcontour(F1); hold on, ezcontour(F2)
% title('Estimated Classes');
% subplot(3,1,3);
% plot(logLikelihood);
% title('Log Likelihood over iterations');

%% dataset3.mat

% load dataset3.mat
% classes = unique(labels);
% class = {};
% for i=1:length(classes)
%     class_idx = find(labels==classes(i));
%     class{i} = data(class_idx, :);
% end
% %figure(1); scatter(class{1}(:,1), class{1}(:,2)); hold on, scatter(class{2}(:,1), class{2}(:,2));
% 
% numClust = 2;
% stopTolerance = 0.0001;
% numRuns =10;
% BIC_flag =1;
% BIC = [];
% if BIC_flag
%     numClustmin = 2;
%     numClustmax = 10;
%     for i =numClustmin:numClustmax
%         disp('cluster'); disp(i);
%         [mu_out, cov_out, prior_out, estimatedLabels, logLikelihood, costvComplex] = ExpectationMaximization(data, i, stopTolerance, numRuns);
%         BIC = [BIC costvComplex];
%     end 
%     figure; plot(numClustmin:numClustmax, BIC); title('dataset3: numClust v BIC');
% else
%     [mu_out, cov_out, prior_out, estimatedLabels, logLikelihood, costvComplex] = ExpectationMaximization(data, numClust, stopTolerance, numRuns)
% end
% 
% classes = unique(estimatedLabels);
% class_estimated = {};
% for i=1:length(classes)
%     class_idx = find(estimatedLabels==classes(i));
%     class_estimated{i} = data(class_idx, :);
% end
% 

% 
% G = gmdistribution(mu_out{1},cov_out{1})
% F = @(x,y) pdf(G,[x y])
% G1 = gmdistribution(mu_out{2},cov_out{2})
% F2 = @(x,y) pdf(G1,[x y])
% 
% figure; 
% subplot(3,1,1);
% scatter(class{1}(:,1), class{1}(:,2)); hold on, scatter(class{2}(:,1), class{2}(:,2));
% title('dataset3, Prob 1dTrue Class');
% subplot(3,1,2);
% scatter(class_estimated{1}(:,1), class_estimated{1}(:,2)); hold on, scatter(class_estimated{2}(:,1), class_estimated{2}(:,2));
% hold on, ezcontour(F); hold on, ezcontour(F2); xlim([-5 5])
% title('Estimated Classes');
% subplot(3,1,3);
% plot(logLikelihood);
% title('Log Likelihood over iterations');
% 
% 


%% dataset4.mat

% load dataset4.mat
% classes = unique(labels);
% class = {};
% for i=1:length(classes)
%     class_idx = find(labels==classes(i));
%     class{i} = data(class_idx, :);
% end
% 
% numClust = 3;
% stopTolerance = 0.00001;
% numRuns =10;
% BIC_flag=1;
% BIC = [];
% if BIC_flag
%     numClustmin = 2;
%     numClustmax = 10;
%     for i =numClustmin:numClustmax
%         disp('cluster'); disp(i);
%         [mu_out, cov_out, prior_out, estimatedLabels, logLikelihood, costvComplex] = ExpectationMaximization(data, i, stopTolerance, numRuns);
%         BIC = [BIC costvComplex];
%     end 
%     figure; plot(numClustmin:numClustmax, BIC); title('dataset4: numClust v BIC');
% else
%     [mu_out, cov_out, prior_out, estimatedLabels, logLikelihood, costvComplex] = ExpectationMaximization(data, numClust, stopTolerance, numRuns);
% end

% 
% 
% classes = unique(estimatedLabels);
% class_estimated = {};
% for i=1:length(classes)
%     class_idx = find(estimatedLabels==classes(i));
%     class_estimated{i} = data(class_idx, :);
% end
% 
% G = gmdistribution(mu_out{1}(1:2),cov_out{1}(1:2,1:2))
% F = @(x,y) pdf(G,[x y])
% G1 = gmdistribution(mu_out{2}(1:2),cov_out{2}(1:2, 1:2))
% F1 = @(x,y) pdf(G1,[x y])
% G2 = gmdistribution(mu_out{3}(1:2),cov_out{3}(1:2, 1:2))
% F2 = @(x,y) pdf(G2,[x y])
% 
% figure(1); 
% subplot(3,1,1);
% scatter(class{1}(:,1), class{1}(:,2)); hold on, scatter(class{2}(:,1), class{2}(:,2)); hold on, scatter(class{3}(:,1), class{3}(:,2));
% title('dataset4, Prob 1e: True Class');
% subplot(3,1,2);
% scatter(class_estimated{1}(:,1), class_estimated{1}(:,2)); hold on, scatter(class_estimated{2}(:,1), class_estimated{2}(:,2));
% scatter(class_estimated{3}(:,1), class_estimated{3}(:,2));
% hold on, ezcontour(F); hold on, ezcontour(F1); hold on, ezcontour(F2);xlim([4 8])
% title('Estimated Classes');
% subplot(3,1,3);
% plot(logLikelihood);
% title('Log Likelihood over iterations');

%% dataset5.mat

% load dataset5.mat
% numClust = 5;
% stopTolerance = 0.00001;
% numRuns =10;
% BIC_flag=1;
% BIC = [];
% if BIC_flag
%     numClustmin = 2;
%     numClustmax = 10;
%     for i =numClustmin:numClustmax
%         disp('cluster'); disp(i);
%         [mu_out, cov_out, prior_out, estimatedLabels, logLikelihood, costvComplex] = ExpectationMaximization(data, i, stopTolerance, numRuns);
%         BIC = [BIC costvComplex];
%     end 
%     figure; plot(numClustmin:numClustmax, BIC); title('dataset4: numClust v BIC');
% else
%     [mu_out, cov_out, prior_out, estimatedLabels, logLikelihood, costvComplex] = ExpectationMaximization(data, numClust, stopTolerance, numRuns);
% end
% 
% classes = unique(estimatedLabels);
% class_estimated = {};
% for i=1:length(classes)
%     class_idx = find(estimatedLabels==classes(i));
%     class_estimated{i} = data(class_idx, :);
% end
% 
% G = gmdistribution(mu_out{1},cov_out{1})
% F = @(x,y) pdf(G,[x y])
% G1 = gmdistribution(mu_out{2},cov_out{2})
% F1 = @(x,y) pdf(G1,[x y])
% G2 = gmdistribution(mu_out{3},cov_out{3})
% F2 = @(x,y) pdf(G2,[x y])
% G3 = gmdistribution(mu_out{4},cov_out{4})
% F3 = @(x,y) pdf(G3,[x y])
% G4 = gmdistribution(mu_out{5},cov_out{5})
% F4 = @(x,y) pdf(G4,[x y])
% 
% 
% figure(1); 
% subplot(3,1,1);
% scatter(data(:,1), data(:,2));
% title('dataset5, Prob 2e');
% subplot(3,1,2);
% scatter(class_estimated{1}(:,1), class_estimated{1}(:,2)); hold on, scatter(class_estimated{2}(:,1), class_estimated{2}(:,2));
% scatter(class_estimated{3}(:,1), class_estimated{3}(:,2)); scatter(class_estimated{4}(:,1), class_estimated{4}(:,2)); hold on,
% scatter(class_estimated{5}(:,1), class_estimated{5}(:,2)); 
% hold on, ezcontour(F); hold on, ezcontour(F1); hold on, ezcontour(F2); hold on,
% ezcontour(F3); hold on, ezcontour(F4); xlim([0 6]); ylim([0 120])
% title('Estimated Classes');
% subplot(3,1,3);
% plot(logLikelihood);
% title('Log Likelihood over iterations');

%% Plotting for test data set 
% c1 = zeros(1, 2);
% c2 = zeros(1,2);
% c3 = zeros(1, 2);
% c4 = zeros(1,2);
% for i = 1:length(estimatedLabels)
%     if estimatedLabels(i) == 1
%         c1_tmp = inputData(i, :);
%         c1 = vertcat(c1, c1_tmp);
%     end
%     if estimatedLabels(i) == 2
%         c2_tmp = inputData(i, :);
%         c2 = vertcat(c2, c2_tmp);
%     end
%     if estimatedLabels(i) == 3
%         c3_tmp = inputData(i, :);
%         c3 = vertcat(c3, c3_tmp);
%     end
%     if estimatedLabels(i) == 4
%         c4_tmp = inputData(i, :);
%         c4 = vertcat(c4, c4_tmp);
%     end
% end
% figure(2); scatter(c1(:,1), c1(:,2)); hold on, scatter(c2(:,1), c2(:,2));
% hold on, scatter(c3(:,1), c3(:,2)); scatter(c4(:,1), c4(:,2));

%figure(3); plot(logLikelihood);