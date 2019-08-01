% File: EM_cluster.m
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

function [P loglikelihood ClassProb] = EM_cluster(poseData, AllG, InitialClassProb, maxIter)

% INPUTS
% poseData: N x 10 x 3 matrix, where N is number of poses;
%   poseData(i,:,:) yields the 10x3 matrix for pose i.
% G: graph parameterization as explained in PA8
% InitialClassProb: N x K, initial allocation of the N poses to the K
%   classes. InitialClassProb(i,j) is the probability that example i belongs
%   to class j
% maxIter: max number of iterations to run EM

% OUTPUTS
% P: structure holding the learned parameters as described in the PA
% loglikelihood: #(iterations run) x 1 vector of loglikelihoods stored for
%   each iteration
% ClassProb: N x K, conditional class probability of the N examples to the
%   K classes in the final iteration. ClassProb(i,j) is the probability that
%   example i belongs to class j

% Initialize variables
N = size(poseData, 1);
K = size(InitialClassProb, 2);

ClassProb = InitialClassProb;

loglikelihood = zeros(maxIter,1);

dimG = size(AllG);
dimFea = size(poseData,3);
init_zeros = zeros(1,K);
P.c = zeros(1,K);
P.clg = repmat(struct('mu_y',init_zeros,'sigma_y',init_zeros,'mu_x',init_zeros,'sigma_x',init_zeros,...
    'mu_angle',init_zeros,'sigma_angle',init_zeros,'theta',zeros(K,(dimFea+1)*dimFea)),dimG(1),1);
if length(dimG) == 2
    AllG = repmat(AllG,1,1,K);
end
% EM algorithm
for iter=1:maxIter
    
    % M-STEP to estimate parameters for Gaussians
    %
    % Fill in P.c with the estimates for prior class probabilities
    % Fill in P.clg for each body part and each class
    % Make sure to choose the right parameterization based on G(i,1)
    %
    % Hint: This part should be similar to your work from PA8
    
    P.c = zeros(1,K);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P.c = mean(ClassProb);
    for k = 1:K
        W = ClassProb(:,k);
        G = squeeze(AllG(:,:,k));
        for i = 1:dimG(1)
            if G(i,1) == 0
                traindata_y = squeeze(poseData(:, i, 1));
                [P.clg(i).mu_y(k), P.clg(i).sigma_y(k)] = FitG(traindata_y, W);
                
                traindata_x = squeeze(poseData(:, i, 2));
                [P.clg(i).mu_x(k), P.clg(i).sigma_x(k)] = FitG(traindata_x, W);
                
                traindata_angle = squeeze(poseData(:, i, 3));
                [P.clg(i).mu_angle(k), P.clg(i).sigma_angle(k)] = FitG(traindata_angle, W);
            else
                parent = G(i,2);
                traindata = squeeze(poseData(:, parent,:));
                
                trainlabel_y = squeeze(poseData(:, i, 1));
                [P.clg(i).theta(k,1:4), P.clg(i).sigma_y(k)] = FitLG(trainlabel_y, traindata, W);
                
                trainlabel_x = squeeze(poseData(:, i, 2));
                [P.clg(i).theta(k,5:8), P.clg(i).sigma_x(k)] = FitLG(trainlabel_x, traindata, W);
                
                trainlabel_angle = squeeze(poseData(:, i, 3));
                [P.clg(i).theta(k,9:12), P.clg(i).sigma_angle(k)] = FitLG(trainlabel_angle, traindata, W);
                
                P.clg(i).theta(k,1:4) = P.clg(i).theta(k,[4,1,2,3]);
                P.clg(i).theta(k,5:8) = P.clg(i).theta(k,[8,5,6,7]);
                P.clg(i).theta(k,9:12) = P.clg(i).theta(k,[12,9,10,11]);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % E-STEP to re-estimate ClassProb using the new parameters
    %
    % Update ClassProb with the new conditional class probabilities.
    % Recall that ClassProb(i,j) is the probability that example i belongs to
    % class j.
    %
    % You should compute everything in log space, and only convert to
    % probability space at the end.
    %
    % Tip: To make things faster, try to reduce the number of calls to
    % lognormpdf, and inline the function (i.e., copy the lognormpdf code
    % into this file)
    %
    % Hint: You should use the logsumexp() function here to do
    % probability normalization in log space to avoid numerical issues
    
    ClassProb = zeros(N,K);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    logll = zeros(N, K);
    for n = 1:N
        for k = 1:K
            mu = zeros(dimG(1), 3);
            sigma = zeros(dimG(1), 3);

            for i = 1:dimG(1)
                if G(i,1) == 0
                    mu_y = P.clg(i).mu_y(k);
                    mu_x = P.clg(i).mu_x(k);
                    mu_angle = P.clg(i).mu_angle(k);
                else
                    mu_y = P.clg(i).theta(k,1:4) * [1,squeeze(poseData(n,G(i,2),:))']';
                    mu_x = P.clg(i).theta(k,5:8) * [1,squeeze(poseData(n,G(i,2),:))']';
                    mu_angle = P.clg(i).theta(k,9:12) * [1,squeeze(poseData(n,G(i,2),:))']';
                end
                sigma_y = P.clg(i).sigma_y(k);
                sigma_x = P.clg(i).sigma_x(k);
                sigma_angle = P.clg(i).sigma_angle(k);

                mu(i,:) = [mu_y, mu_x, mu_angle];
                sigma(i,:) = [sigma_y, sigma_x, sigma_angle];
            end
            logll_n_k = lognormpdf(squeeze(poseData(n,:,:)), mu, sigma);
            logll(n,k) = log(P.c(k)) + sum(sum(logll_n_k));
        end
        tmp_logll = logll(n,:) - max(logll(n,:));
        ClassProb(n,:) = exp(tmp_logll) / sum(exp(tmp_logll));
        loglikelihood(iter) = loglikelihood(iter) + logsumexp(logll(n,:));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compute log likelihood of dataset for this iteration
    % Hint: You should use the logsumexp() function here
    %     loglikelihood(iter) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Print out loglikelihood
    disp(sprintf('EM iteration %d: log likelihood: %f', ...
        iter, loglikelihood(iter)));
    if exist('OCTAVE_VERSION')
        fflush(stdout);
    end
    
    % Check for overfitting: when loglikelihood decreases
    if iter > 1
        if loglikelihood(iter) < loglikelihood(iter-1)
            break;
        end
    end
    
end

% Remove iterations if we exited early
loglikelihood = loglikelihood(1:iter);
end
