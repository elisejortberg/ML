%hw5 bonus EM for Bernoulli Distribution

load digitsSmall.mat;
%run for digits 0, 1, 4
%index = digit
%1 = 5
%2 = 0
%3 = 4
%4 = 1

% data = reshape(inputData(2,:),28,[]);

%data = inputData(labels==0 | labels==1 |labels==4, :) >0.5;
data = inputData;


params = [];
params.numberOfRuns = 10;
params.numberOfClusters = 20;

[clusterParameters, dataProbabilities, estimatedLabels, outputLog] = EMClustering_Bernoulli(data, params);

%% Part D: k=3
% delta_k1 = clusterParameters(1).mu - mean(inputData(labels==4,:)>0.5);
% delta_k2 = clusterParameters(2).mu - mean(inputData(labels==0,:)>0.5);
% delta_k3 = clusterParameters(3).mu - mean(inputData(labels==1,:)>0.5);
% 
% figure; title('Top: EM Results, Bottom: Estimated cluster means - True means')
% subplot(2,3,1); imagesc(reshape(clusterParameters(1).mu, 28, [])); colorbar;
% subplot(2,3,2); imagesc(reshape(clusterParameters(2).mu, 28, [])); colorbar;
% subplot(2,3,3); imagesc(reshape(clusterParameters(3).mu, 28, [])); colorbar;
% subplot(2,3,4); imagesc(reshape(delta_k1, 28, [])); colorbar;
% subplot(2,3,5); imagesc(reshape(delta_k2, 28, [])); colorbar;
% subplot(2,3,6); imagesc(reshape(delta_k3, 28, [])); colorbar;

%% Part e: k=5
% figure; 
% subplot(1, 5, 1); imagesc(reshape(clusterParameters(1).mu, 28, [])); colorbar;
% subplot(1, 5, 2); imagesc(reshape(clusterParameters(2).mu, 28, [])); colorbar;
% subplot(1, 5, 3); imagesc(reshape(clusterParameters(3).mu, 28, [])); colorbar;
% subplot(1, 5, 4); imagesc(reshape(clusterParameters(4).mu, 28, [])); colorbar;
% subplot(1, 5, 5); imagesc(reshape(clusterParameters(5).mu, 28, [])); colorbar;

%% Part e: k=10
% figure; 
% subplot(2, 5, 1); imagesc(reshape(clusterParameters(1).mu, 28, [])); colorbar;
% subplot(2, 5, 2); imagesc(reshape(clusterParameters(2).mu, 28, [])); colorbar;
% subplot(2, 5, 3); imagesc(reshape(clusterParameters(3).mu, 28, [])); colorbar;
% subplot(2, 5, 4); imagesc(reshape(clusterParameters(4).mu, 28, [])); colorbar;
% subplot(2, 5, 5); imagesc(reshape(clusterParameters(5).mu, 28, [])); colorbar;
% subplot(2, 5, 6); imagesc(reshape(clusterParameters(6).mu, 28, [])); colorbar;
% subplot(2, 5, 7); imagesc(reshape(clusterParameters(7).mu, 28, [])); colorbar;
% subplot(2, 5, 8); imagesc(reshape(clusterParameters(8).mu, 28, [])); colorbar;
% subplot(2, 5, 9); imagesc(reshape(clusterParameters(9).mu, 28, [])); colorbar;
% subplot(2, 5, 10); imagesc(reshape(clusterParameters(10).mu, 28, [])); colorbar;
%% Part e: k=20
figure; 
subplot(4, 5, 1); imagesc(reshape(clusterParameters(1).mu, 28, [])); colorbar;
subplot(4, 5, 2); imagesc(reshape(clusterParameters(2).mu, 28, [])); colorbar;
subplot(4, 5, 3); imagesc(reshape(clusterParameters(3).mu, 28, [])); colorbar;
subplot(4, 5, 4); imagesc(reshape(clusterParameters(4).mu, 28, [])); colorbar;
subplot(4, 5, 5); imagesc(reshape(clusterParameters(5).mu, 28, [])); colorbar;
subplot(4, 5, 6); imagesc(reshape(clusterParameters(6).mu, 28, [])); colorbar;
subplot(4, 5, 7); imagesc(reshape(clusterParameters(7).mu, 28, [])); colorbar;
subplot(4, 5, 8); imagesc(reshape(clusterParameters(8).mu, 28, [])); colorbar;
subplot(4, 5, 9); imagesc(reshape(clusterParameters(9).mu, 28, [])); colorbar;
subplot(4, 5, 10); imagesc(reshape(clusterParameters(10).mu, 28, [])); colorbar;
subplot(4, 5, 11); imagesc(reshape(clusterParameters(11).mu, 28, [])); colorbar;
subplot(4, 5, 12); imagesc(reshape(clusterParameters(12).mu, 28, [])); colorbar;
subplot(4, 5, 13); imagesc(reshape(clusterParameters(13).mu, 28, [])); colorbar;
subplot(4, 5, 14); imagesc(reshape(clusterParameters(14).mu, 28, [])); colorbar;
subplot(4, 5, 15); imagesc(reshape(clusterParameters(15).mu, 28, [])); colorbar;
subplot(4, 5, 16); imagesc(reshape(clusterParameters(16).mu, 28, [])); colorbar;
subplot(4, 5, 17); imagesc(reshape(clusterParameters(17).mu, 28, [])); colorbar;
subplot(4, 5, 18); imagesc(reshape(clusterParameters(18).mu, 28, [])); colorbar;
subplot(4, 5, 19); imagesc(reshape(clusterParameters(19).mu, 28, [])); colorbar;
subplot(4, 5, 20); imagesc(reshape(clusterParameters(20).mu, 28, [])); colorbar;