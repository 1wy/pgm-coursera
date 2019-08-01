function loglikelihood = ComputeLogLikelihood(P, G, dataset)
% returns the (natural) log-likelihood of data given the model and graph structure
%
% Inputs:
% P: struct array parameters (explained in PA description)
% G: graph structure and parameterization (explained in PA description)
%
%    NOTICE that G could be either 10x2 (same graph shared by all classes)
%    or 10x2x2 (each class has its own graph). your code should compute
%    the log-likelihood using the right graph.
%
% dataset: N x 10 x 3, N poses represented by 10 parts in (y, x, alpha)
%
% Output:
% loglikelihood: log-likelihood of the data (scalar)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

N = size(dataset,1); % number of examples
K = length(P.c); % number of classes

loglikelihood = 0;
% You should compute the log likelihood of data as in eq. (12) and (13)
% in the PA description
% Hint: Use lognormpdf instead of log(normpdf) to prevent underflow.
%       You may use log(sum(exp(logProb))) to do addition in the original
%       space, sum(Prob).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimlen = size(G);
dim3 = [];
if dimlen ==  3
    dim3 = 1:dimlen(3);
else
    dim3 = 1;
end


for n = 1:N
    logll_n = zeros(K,1);
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
        logll_n(k,1) = log(P.c(k)) + sum(sum(logll_n_k));
    end
    loglikelihood = loglikelihood + logsumexp(logll_n); 
end
end