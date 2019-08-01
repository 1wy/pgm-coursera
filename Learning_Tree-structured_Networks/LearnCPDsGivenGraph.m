function [P loglikelihood] = LearnCPDsGivenGraph(dataset, AllG, labels)
%
% Inputs:
% dataset: N x 10 x 3, N poses represented by 10 parts in (y, x, alpha)
% G: graph parameterization as explained in PA description
% labels: N x 2 true class labels for the examples. labels(i,j)=1 if the 
%         the ith example belongs to class j and 0 elsewhere        
%
% Outputs:
% P: struct array parameters (explained in PA description)
% loglikelihood: log-likelihood of the data (scalar)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

N = size(dataset, 1);
K = size(labels,2);

loglikelihood = 0;
P.c = zeros(1,K);
dimG = size(AllG);
dimFea = size(dataset, 3);


if length(dimG) == 2
    AllG = repmat(AllG,1,1,2);
end
% estimate parameters
% fill in P.c, MLE for class probabilities
% fill in P.clg for each body part and each class
% choose the right parameterization based on G(i,1)
% compute the likelihood - you may want to use ComputeLogLikelihood.m
% you just implemented.
%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE

% These are dummy lines added so that submit.m will run even if you 
% have not started coding. Please delete them.
P.c = mean(labels);
init_zeros = zeros(1,K);
P.clg = repmat(struct('mu_y',init_zeros,'sigma_y',init_zeros,'mu_x',init_zeros,'sigma_x',init_zeros,...
                      'mu_angle',init_zeros,'sigma_angle',init_zeros,'theta',zeros(K,(dimFea+1)*dimFea)),dimG(1),1);

for k = 1:K
    idx = logical(labels(:,k));
    G = squeeze(AllG(:,:,k));
    for i = 1:dimG(1)
        if G(i,1) == 0
            traindata_y = squeeze(dataset(idx, i, 1)); 
            [P.clg(i).mu_y(k), P.clg(i).sigma_y(k)] = FitGaussianParameters(traindata_y);
            
            traindata_x = squeeze(dataset(idx, i, 2)); 
            [P.clg(i).mu_x(k), P.clg(i).sigma_x(k)] = FitGaussianParameters(traindata_x);
            
            traindata_angle = squeeze(dataset(idx, i, 3)); 
            [P.clg(i).mu_angle(k), P.clg(i).sigma_angle(k)] = FitGaussianParameters(traindata_angle);
        else
            parent = G(i,2);
            traindata = squeeze(dataset(idx, parent,:));
            
            trainlabel_y = squeeze(dataset(idx, i, 1));
            [P.clg(i).theta(k,1:4), P.clg(i).sigma_y(k)] = FitLinearGaussianParameters(trainlabel_y, traindata);
            
            trainlabel_x = squeeze(dataset(idx, i, 2));
            [P.clg(i).theta(k,5:8), P.clg(i).sigma_x(k)] = FitLinearGaussianParameters(trainlabel_x, traindata);
            
            trainlabel_angle = squeeze(dataset(idx, i, 3));
            [P.clg(i).theta(k,9:12), P.clg(i).sigma_angle(k)] = FitLinearGaussianParameters(trainlabel_angle, traindata);
            
            P.clg(i).theta(k,1:4) = P.clg(i).theta(k,[4,1,2,3]);
            P.clg(i).theta(k,5:8) = P.clg(i).theta(k,[8,5,6,7]);
            P.clg(i).theta(k,9:12) = P.clg(i).theta(k,[12,9,10,11]);
        end
    end
end
for k = 1:K
    idx = logical(labels(:,k));
    G = squeeze(AllG(:,:,k));
    loglikelihood = loglikelihood + ComputeLogLikelihood(P, G, dataset(idx,:,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('log likelihood: %f\n', loglikelihood);
end

