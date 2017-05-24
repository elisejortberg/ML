function [MSE_fin] = MSEcalculator(classData, mu)
    k = length(mu);
    MSE_fin = 0;
    MSE_per_dimension=0;
    num_samp = 0;
    for i=1:k
        [n, d] = size(classData{i});
        num_samp = num_samp+n;
        for j=1:d
            %MSE_t = norm(classData{i}(:,j) - mu{i}(j))^2; 
            MSE_t = sqrt(sum((classData{i}(:,j) - mu{i}(j)).^2));
            MSE_per_dimension = MSE_per_dimension+MSE_t;
        end
        MSE_fin = (MSE_fin + MSE_per_dimension);
    end 
    MSE_fin = MSE_fin/num_samp;
end
