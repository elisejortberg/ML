load digitsSmall.mat
trained=1;
if trained
    SVMModels = train(inputData, labels); %train models, one v all
    classes = unique(labels);
end
load digitsLarge.mat

score_total = zeros(length(inputData), length(classes));

for i =1:length(classes)
    score = test(SVMModels{i}, inputData);
    score_total(:,i) = score;
end

label_classified= zeros(length(inputData), 1);
for i =1:length(inputData)
    [~, idx] = max(score_total(i,:)); % pick max score as final classification
    label_classified(i) = idx-1; 
end

C = confusionmat(labels, label_classified);
true_cntr=0;
for i=1:length(labels)
    if labels(i) == label_classified(i)
        true_cntr= true_cntr+1;
    end
end
accuracy= true_cntr/length(labels);
figure; image(C); title(['Accuracy = ', num2str(accuracy)]); colorbar; hold on, 
xlabel('Predicted Labels'); ylabel('Actual Labels')

