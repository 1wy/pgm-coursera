function [P G loglikelihood] = LearnGraphAndCPDs(dataset, labels)

% dataset: N x 10 x 3, N poses represented by 10 parts in (y, x, alpha) 
% labels: N x 2 true class labels for the examples. labels(i,j)=1 if the 
%         the ith example belongs to class j
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

N = size(dataset, 1);
K = size(labels,2);

G = zeros(10,2,K); % graph structures to learn
% initialization
for k=1:K
    G(2:end,:,k) = ones(9,2);
end

% estimate graph structure for each class
for k=1:K
    % fill in G(:,:,k)
    % use ConvertAtoG to convert a maximum spanning tree to a graph G
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%
    idx = logical(labels(:,k));
    [A W] = LearnGraphStructure(dataset(idx,:,:));
    G(:,:,k) = ConvertAtoG(A);
end

% estimate parameters

% P.c = zeros(1,K);
% compute P.c
% the following code can be copied from LearnCPDsGivenGraph.m
% with little or no modification
%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE

[P loglikelihood] = LearnCPDsGivenGraph(dataset, G, labels);

% These are dummy lines added so that submit.m will run even if you 
% have not started coding. Please delete them.
% P.clg.sigma_x = 0;
% P.clg.sigma_y = 0;
% P.clg.sigma_angle = 0;
% loglikelihood = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('log likelihood: %f\n', loglikelihood);
end