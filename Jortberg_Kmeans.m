function [estimatedLabels, mu_estimated, MSE_array] = Jortberg_Kmeans(inputData, numClust, stopTolerance, numberOfRuns, true_class)
% Inputs:
    % inputData: nSamples x dDimensions array
    % numberOfClusters: number of clusters to be estimated
    % stopTolerance: convergence criteria (when RMSE < x?)
    % numberOfRuns: number times algorithm runs with random initializations
% Outputs:
    % estimatedLabels: nSamples x 1 vector with estimated cluster assg't
    % estimatedMeans: numberOfClusters x dDimensions with est'd cluster mus
    % MSE: nIterations x 1 vector with MSE as function of iteration #
   
    %for now initialize once, repeat random initialization after
    %while err > 0.2
    
    %1. initialize n, c, mu's 1:c
    
    %n points in d dimensional space
    %MSE = 1;
    [n, d] = size(inputData);
    
    %initialize centroids numberOfRuns times
    for t =1:numberOfRuns
        MSE_array = [];
        MSE_prior = 0;
        delta_MSE = 10;
        centroids = inputData(randsample(n, numClust), :);
        mu_estimated= {};
        for i = 1:numClust
            mu_estimated{i} = centroids(i,:);
        end
        
       cntr = 0;
       while abs(delta_MSE) >= stopTolerance 
           
        %2. classify, for each x_i, assign x_i to closest mean (i=1:n)
            % 2a) calculate closest distance use euclidean distance
            euc_dist = zeros(n, numClust);
            for i = 1:numClust
              edist = inputData - repmat(mu_estimated{i}, [n 1]); %shouldnt be subtracting inputData but updated classes???
              edist = sqrt(sum((edist).^2, 2));
              euc_dist(:,i) = edist;
            end
            % 2b) reclassify based on smallet distance
            estimatedLabels = zeros(n, 1);
            for i=1:n
                estimatedLabels(i) = find(euc_dist(i,:) == min(euc_dist(i,:)), 1);
            end

        %3. recalculate means
        %separate out by predicted class then calculate the means(numClust x d)
            class_updated = {};
            mu_estimated = {};
            mu_t = zeros(1, d);
            for i=1:numClust
                idx = find(estimatedLabels==i); 
                class_updated{i} = inputData(idx, :);
                for j =1:d
                    mu_t(j) = mean(class_updated{i}(:, j)); %find mean per dimension
                end
                mu_estimated{i} = mu_t; %update means
            end

        %4. check convergence. repeat 2, 3 if not converged
            %Sum Squared Error SSE 
            %SSE = sum over all c, sum over all x in c, norm(x-mu)^2
           MSE = MSEcalculator(class_updated, mu_estimated);

           %MSE_array = [MSE_array delta_MSE];
           MSE_array = [MSE_array MSE];
           delta_MSE = MSE_prior - MSE;
           MSE_prior = MSE;
        cntr = cntr+1;
       end %end while loop
       
       figure; subplot(3,1,1);%input real classification for plots
       for i=1:numClust
            title('K Means Clustered Data')
            scatter(class_updated{i}(:,1), class_updated{i}(:,2)); hold on,
            scatter(mu_estimated{i}(1), mu_estimated{i}(2), 50, 'kx', 'linewidth', 3.0);
       end
       subplot(3,1,2); plot(MSE_array); title('MSE as function of iteration'); hold on,
       subplot(3,1,3); 
       k_true = length(true_class);
       colors = ['r', 'g', 'm'];
       for m=1:k_true
            title('True Classes')
            meantrue1 = mean(true_class{m}(:,1));
            meantrue2 = mean(true_class{m}(:,2));
            scatter(true_class{m}(:,1), true_class{m}(:,2), colors(m)); hold on,
            scatter(meantrue1, meantrue2, 50, 'kx', 'linewidth', 3.0);
       end
       
    end %end for loop, number of runs
    
    %Run = user specified number of re initializations
    %Iteration = number of times while loop runs
    
    %Plot - pick plot with smallest MSE
    %if MSE_array(end) < MSE_array(end) run before then plot, else don't
%     figure; 
%     for i=1:numClust
%         title('Clustered Data')
%         scatter(class_updated{i}(:,1), class_updated{i}(:,2)); hold on,
%         scatter(mu_estimated{i}(1), mu_estimated{i}(2), 50, 'kx', 'linewidth', 3.0);
%     end
%     figure; plot(MSE_array); title('MSE as function of iteration');
end

%% NOTES
%         euc_dist_c1 = inputData-repmat(centroids(1,:), [length(inputData) 1]);
%         euc_dist_c1 = sqrt(sum((euc_dist_c1).^2, 2));
%         euc_dist_c2 = inputData-repmat(centroids(2,:), [length(inputData) 1]);
%         euc_dist_c2 = sqrt(euc_dist_c2(:,1).^2 + euc_dist_c2(:,2).^2); 
%         class1 = inputData[euc_dist_c1>euc_dist_c2];
%         class2 = inputData[euc_dist_c2>euc_dist_c1];

%            for i=1:numClust
%                for j = 1:d
%                    MSE_t = (norm(inputData(:,d) - repmat(mu_estimated{i}(:,d), [n, 1])).^2); %is this what I want, max(abs(sum))^2?
%                    %MSE_t = abs(inputData- repmat(mu_estimated{i},[n,1])).^2;
%                    MSE= MSE + MSE_t;
%                    if i==numClust
%                      delta_MSE = MSE_prior-MSE;
%                    end
%                end
%            end 


