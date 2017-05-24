%1.3 
%a)
nsamples = 400;
mu{1} = [0,0]; 
mu{2} = [3,3];
I = eye(2); 
sigma{1} = I; 
sigma{2} = I;
prior = [0.5, 0.5];
[data, classIndex] = generateGaussianSamplesv3(mu, sigma, nsamples, prior);
figure(1); plot(data{1}(:,1), data{1}(:, 2), '.'); hold on, plot(data{2}(:,1), data{2}(:, 2), 'r.');
hold on, xlabel('x1'), ylabel('x2'), title('Distribution a'), legend('class 1', 'class 2')
[g, class1, class2, xbound, ybound1, ~, accuracy] = discriminant(data, classIndex, mu, sigma, nsamples, prior, 1);

figure(2); plot(class1(:,1), class1(:,2), 'b.'); hold on, plot(class2(:,1), class2(:,2), 'r.'); title(sprintf('accuracy = %0.5g', accuracy));
hold on, xlabel('x1'); ylabel('x2'); plot(xbound, ybound1, 'g-');
legend('class1', 'class2', 'case1 boundary');


%b)
sigma{1} = [3, 1; 1, 0.8]; 
sigma{2} = sigma{1};
[data, classIndex] = generateGaussianSamplesv3(mu, sigma, nsamples, prior);
figure(3); plot(data{1}(:,1), data{1}(:, 2), '.'); hold on, plot(data{2}(:,1), data{2}(:, 2), 'r.');
hold on, xlabel('x1'), ylabel('x2'), title('Distribution b'), legend('class 1', 'class 2')
[g, class1, class2, xbound, ybound1, ybound2, accuracy] = discriminant(data, classIndex, mu, sigma, nsamples, prior, 2);

figure(4); plot(class1(:,1), class1(:,2), 'b.'); hold on, plot(class2(:,1), class2(:,2), 'r.'); title(sprintf('accuracy = %0.5g', accuracy));
hold on, xlabel('x1'); ylabel('x2'); plot(xbound, ybound1, 'g-');
hold on, plot(xbound, ybound2, 'r'); legend('class1', 'class2', 'case1 boundary', 'case2 boundary');

% %c) 
sigma{1} = [2, 0.5; 0.5, 1];
sigma{2} = [2 -1.9; -1.9 5];
mu{2} = [2, 2];
[data, classIndex] = generateGaussianSamplesv3(mu, sigma, nsamples, prior);
figure(5); plot(data{1}(:,1), data{1}(:, 2), '.'); hold on, plot(data{2}(:,1), data{2}(:, 2), 'r.');
hold on, xlabel('x1'), ylabel('x2'), title('Distribution c'), legend('class 1', 'class 2')
[g, class1, class2, xbound, ybound1, ybound2, accuracy] = discriminant(data, classIndex, mu, sigma, nsamples, prior, 3);
figure(6); plot(class1(:,1), class1(:,2), 'b.'); hold on, plot(class2(:,1), class2(:,2), 'r.'); title(sprintf('accuracy = %0.5g', accuracy));
hold on, xlabel('x1'); ylabel('x2'); plot(xbound, ybound1, 'g-');
hold on, plot(xbound, ybound2, 'r'); legend('class1', 'class2', 'case1 boundary', 'case2 boundary');


%d)
mu{1} = [0,0]; 
mu{2} = [3,3];
I = eye(2); 
sigma{1} = I; 
sigma{2} = I;
priors = [0.05, 0.95];
[data, classIndex] = generateGaussianSamplesv3(mu, sigma, nsamples, priors);
figure(7); plot(data{1}(:,1), data{1}(:, 2), '.'); hold on, plot(data{2}(:,1), data{2}(:, 2), 'r.');
hold on, xlabel('x1'), ylabel('x2'), title('Distribution d'), legend('class 1', 'class 2');
acc = classify_c3(nsamp, mu, sigma, priors);

%e)
sigma{1} = [3, 1; 1, 0.8]; 
sigma{2} = sigma{1};
[data, classIndex] = generateGaussianSamplesv3(mu, sigma, nsamples, priors);
figure(8); plot(data{1}(:,1), data{1}(:, 2), '.'); hold on, plot(data{2}(:,1), data{2}(:, 2), 'r.');
hold on, xlabel('x1'), ylabel('x2'), title('Distribution e'), legend('class 1', 'class 2');
acc = classify_c3(nsamp, mu, sigma, priors);

%f)
sigma{1} = [2, 0.5; 0.5, 1];
sigma{2} = [2 -1.9; -1.9 5];
mu{2} = [2, 2];
[data, classIndex] = generateGaussianSamplesv3(mu, sigma, nsamples, prior);
figure(9); plot(data{1}(:,1), data{1}(:, 2), '.'); hold on, plot(data{2}(:,1), data{2}(:, 2), 'r.');
hold on, xlabel('x1'), ylabel('x2'), title('Distribution f'); legend('class 1', 'class 2');
acc = classify_c3(nsamp, mu, sigma, priors);


