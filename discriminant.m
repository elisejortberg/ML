function [g, class1, class2, xbound, ybound1, ybound2, accuracy] = discriminant(data, classIndex, mu, sigma, nsamples, prior, part)
%inputs: similar to generateGaussianSamples
%outputs: 
%g = nsamples x k array, with score of class in columns, score of data
    %point in rows
%accuracy = how many times output matches classIndex
%class1 = class1 data points
%class2 = class 2 data points
%xbound = xboundary
%ybound1 = yboundary for class 1
%ybound2 = yboundary for class 2

data_all = vertcat(data{1}, data{2});
class_true = vertcat(classIndex{1}, classIndex{2});
d = length(mu{1});

%g = (-1/2)*(x-mu)'*inv(sigma)*(x-mu)  - (1/2)*ln(det(sigma)) + ln(prior);
%g = (-1/2)*(x)'*inv(sigma)*(x)  - (1/2)*log(det(sigma)) + log(prior(1));

sigma1 = sigma{1};
sigma2 = sigma{2};
mu1 = mu{1}';
mu2 = mu{2}';
prior1 = prior(1);
prior2 = prior(2);

W1 = (-1/2)*inv(sigma1);
W2 = (-1/2)*inv(sigma2);
w1 = inv(sigma1)*mu1;
w2 = inv(sigma2)*mu2;
wo1 = (-1/2)*mu1'*inv(sigma1)*mu1 - (1/2)*log(det(sigma1)) + log(prior1);
wo2 = (-1/2)*mu2'*inv(sigma2)*mu2 - (1/2)*log(det(sigma2)) + log(prior2);

g1 = zeros(nsamples,1);
g2 = zeros(nsamples,1);
class = zeros(nsamples,1);
class1 = [];
class2 = [];
TP=0;

%part a
xbound = [-10:0.1:10];
if part ==1
    ybound1 = 3-xbound;
    ybound2 = [];
%part b
% case I: same as part a)
% case II: 
elseif part ==2
    ybound1 = 3-xbound;
    ybound2 = (99/60)-(xbound/10);
% part c
elseif part==3
    ybound1 = 2-xbound;
    ybound2 = (64/27) - (37/27)*xbound;
else
    ybound1 = [];
    ybound2 = [];
end

%figure; plot(xbound, ybound)
for i=1:nsamples
    x = data_all(i, :)';
    g1(i) = x'*W1*x + w1'*x + wo1;
    g2(i) = x'*W2*x + w2'*x + wo2;
    
    if g1(i)>g2(i)
        class(i) = 1;
        class1 = vertcat(class1, x');
    else
        class(i) = 2;
        class2 = vertcat(class2, x');
    end
    if class(i)==class_true(i)
       TP = TP+1;
    end
end

accuracy = TP/nsamples;

% disp('size1, 2'); disp(size(class1)); disp(size(class2));
% figure; plot(class1(:,1), class1(:,2), 'b.'); hold on, plot(class2(:,1), class2(:,2), 'r.'); title(sprintf('accuracy = %0.5g', accuracy));
% hold on, xlabel('x1'); ylabel('x2'); plot(xbound, ybound1, 'g-');
% hold on, plot(xbound, ybound2, 'r'); legend('class1', 'class2', 'case1 boundary', 'case2 boundary');
g = horzcat(g1, g2);

% W = (-1/2)*inv(sigma)
% w = inv(sigma)*mu
% wo = (-1/2)*mu'*inv(sigma)*mu - (1/2)*log(det(sigma)) + log(prior)
% g = x'*W*x + w'x + wo
%boundary = when priors equal, g1(x) = g2(x)



end
