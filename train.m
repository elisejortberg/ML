function SVMModel = train(data, labels)
    %train 10 times against each digit
    numbers = unique(labels);
    SVMModel = cell(1, length(numbers));
    %create 10 models, one for each value
    for i=1:length(numbers) %matlab doesn't do 0 indices
        true_class= zeros(1, length(labels));
        disp(i-1);

        true_class(labels==numbers(i)) = 1;
        %part 1) Regular
        %SVMModel{i} = fitcsvm(data, true_class); 
        %part 2) RBF Kernel %better if cross validate
        %SVMModel{i} = fitcsvm(data, true_class,'Standardize',true,'KernelFunction','RBF','KernelScale','auto');
        %part 3) polynomial, 3rd order (default) %bad performance without auto
        SVMModel{i} = fitcsvm(data, true_class,'Standardize',true,'KernelFunction','polynomial','KernelScale','auto');
        
    end
end