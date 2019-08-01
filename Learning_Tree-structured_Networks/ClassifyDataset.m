function accuracy = ClassifyDataset(dataset, labels, P, G)
% returns the accuracy of the model P and graph G on the dataset 
%
% Inputs:
% dataset: N x 10 x 3, N test instances represented by 10 parts
% labels:  N x 2 true class labels for the instances.
%          labels(i,j)=1 if the ith instance belongs to class j 
% P: struct array model parameters (explained in PA description)
% G: graph structure and parameterization (explained in PA description) 
%
% Outputs:
% accuracy: fraction of correctly classified instances (scalar)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

N = size(dataset, 1);
accuracy = 0.0;
K = size(labels,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loglikelihood = zeros(N,K);
dimlen = size(G);
mu = zeros(dimlen(1), 3);
sigma = zeros(dimlen(1), 3);

for n = 1:N
    for k = 1:K
        mu = zeros(dimlen(1), 3);
        sigma = zeros(dimlen(1), 3);
        
        for i = 1:dimlen(1)
            if G(i,1) == 0
                mu_y = P.clg(i).mu_y(k);
                mu_x = P.clg(i).mu_x(k);             
                mu_angle = P.clg(i).mu_angle(k);
            else
                mu_y = P.clg(i).theta(k,1:4) * [1,squeeze(dataset(n,G(i,2),:))']';
                mu_x = P.clg(i).theta(k,5:8) * [1,squeeze(dataset(n,G(i,2),:))']';
                mu_angle = P.clg(i).theta(k,9:12) * [1,squeeze(dataset(n,G(i,2),:))']';                
            end   
             sigma_y = P.clg(i).sigma_y(k);
             sigma_x = P.clg(i).sigma_x(k);
             sigma_angle = P.clg(i).sigma_angle(k);
             
             mu(i,:) = [mu_y, mu_x, mu_angle];
             sigma(i,:) = [sigma_y, sigma_x, sigma_angle];
        end
        logll_n_k = lognormpdf(squeeze(dataset(n,:,:)), mu, sigma);
        loglikelihood(n,k) = log(P.c(k)) + sum(sum(logll_n_k));
    end
end
[maxnum, predict] = max(loglikelihood,[],2); 
[maxnum, groudtruth] = max(labels,[],2);
accuracy = sum(predict == groudtruth) / length(groudtruth);
fprintf('Accuracy: %.2f\n', accuracy);
end