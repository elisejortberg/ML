function log_prior = loglikelihood(data, mu, cov, prior, numClust)
    sum_per_point = 0;
    [n,d] = size(data);
    
    for i =1:n %for each data point, find likelihood to belong to each cluster
        sum_total = 0;
        datapoint = data(i, :);
        for j =1:numClust
            sigma = abs(cov{j});
            [~, err] = cholcov(sigma,0); %check that sigma is valid
            if err~=0
                %sigma = eye(d)*0.001;
                sigma = sigma + .001 * eye(d);
                %sigma = eye(d);
            end
            %try
                %sum_total = sum_total + prior*mvnpdf(datapoint, mu{j}, sigma); 
            %catch
                my_gaussian = (((((2*pi)^(d))*det(sigma))^(-(1/2)))   *  exp((-1/2)*(datapoint-mu{j}) * inv(sigma) * (datapoint-mu{j})'));
                sum_total = sum_total + prior*my_gaussian;
            %end
        end
        sum_per_point = sum_per_point + log(sum_total);
        
    end
    log_prior = sum_per_point;
end
