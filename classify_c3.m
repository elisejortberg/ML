function acc = classify_c3(nsamp, mu, sigma, priors)
% Create space of x = (x1, x2)
x1 = linspace(-5,5,nsamp);
x2 = linspace(-5,5,nsamp);

[X1, X2] = meshgrid(x1,x2);

%priors for class
priorClass1 = priors(1);
priorClass2 = priors(2);

mu1 = mu{1};
mu2 = mu{2};
sig1 = sigma{1};
sig2 = sigma{2};

% Evaluate class likelihoods in grid
probClass1 = reshape(mvnpdf([X1(:), X2(:)], mu1, sig1), [], nsamp );
probClass2 = reshape(mvnpdf([X1(:), X2(:)], mu2, sig2), [], nsamp );

% Compute posterior: posterior = likelihood*prior / sum(likelihood*prior for all classes)
posteriorClass1 = probClass1*priorClass1./(probClass1*priorClass1 + probClass2*priorClass2);
posteriorClass2 = probClass2*priorClass2./(probClass1*priorClass1 + probClass2*priorClass2);

% Plot likelihoods and posteriors
% figure(1)
% imagesc(x1, x2, probClass1);
% title('p(x | \omega_1)')
% xlabel('x1')
% xlabel('x2')
% 
% figure(2)
% imagesc(x1, x2, probClass2);
% title('p(x | \omega_2)')
% xlabel('x1')
% xlabel('x2')

% Create decision boundary map by choosing the class index that corresponds
% to the maximum posterior probability
[~, decionBoundaryMap] = max([posteriorClass1(:), posteriorClass2(:)], [], 2);
[~, accuracyMap] = max([probClass1(:), probClass2(:)], [], 2);
decionBoundaryMap = reshape(decionBoundaryMap, [], nsamp);
accuracyMap = reshape(accuracyMap, [], nsamp);

% Plot
acc =(accuracyMap-decionBoundaryMap);
acc = 1+mean(acc(:));

figure;
imagesc(x1, x2, decionBoundaryMap);
title(sprintf('Boundary. Acc = %g', acc))
xlabel('x1')
xlabel('x2')
end