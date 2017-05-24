function Q = ExpectationStep(inputData, likelihood, mu, cov, numClust)
    %Estimate missing data
    [n, d] = size(inputData);
    Q = zeros(n, numClust);
    for i =1:n
        norm_c = 0; %normalization (pi factor)
        datapoint = inputData(i, :);
        for j = 1:numClust
            %Q(i,j) = likelihood(j)*mvnpdf(datapoint, mu{j}, cov{j});
            my_gaussian = (((((2*pi)^(d))*det(cov{j}))^(-(1/2)))   *  exp((-1/2)*(datapoint-mu{j}) * inv(cov{j}) * (datapoint-mu{j})'));
            Q(i,j) = likelihood(j)*my_gaussian;
            norm_c = norm_c + Q(i,j);
        end
        %Normalize each datapoints calculated Q by total
        for j=1:numClust %normalize per cluster
            Q(i,j) = Q(i,j)/norm_c; 
        end
        
    end
end