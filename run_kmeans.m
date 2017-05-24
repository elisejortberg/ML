%load dataset1
classes = unique(labels);
true_class = {};
for i=1:length(classes)
    class_idx = find(labels==classes(i));
    true_class{i} = data(class_idx, :);
end

k = 2;
stopTolerance = 0.000001;
runs = 10;
[estimatedLabels, mu_estimated, MSE_array] = Jortberg_Kmeans(data, k, stopTolerance, runs, true_class);