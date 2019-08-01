% File: EM_HMM.m
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

function [P loglikelihood ClassProb PairProb] = EM_HMM(actionData, poseData, G, InitialClassProb, InitialPairProb, maxIter)

% INPUTS
% actionData: structure holding the actions as described in the PA
% poseData: N x 10 x 3 matrix, where N is number of poses in all actions
% G: graph parameterization as explained in PA description
% InitialClassProb: N x K matrix, initial allocation of the N poses to the K
%   states. InitialClassProb(i,j) is the probability that example i belongs
%   to state j.
%   This is described in more detail in the PA.
% InitialPairProb: V x K^2 matrix, where V is the total number of pose
%   transitions in all HMM action models, and K is the number of states.
%   This is described in more detail in the PA.
% maxIter: max number of iterations to run EM

% OUTPUTS
% P: structure holding the learned parameters as described in the PA
% loglikelihood: #(iterations run) x 1 vector of loglikelihoods stored for
%   each iteration
% ClassProb: N x K matrix of the conditional class probability of the N examples to the
%   K states in the final iteration. ClassProb(i,j) is the probability that
%   example i belongs to state j. This is described in more detail in the PA.
% PairProb: V x K^2 matrix, where V is the total number of pose transitions
%   in all HMM action models, and K is the number of states. This is
%   described in more detail in the PA.

% Initialize variables
N = size(poseData, 1);
K = size(InitialClassProb, 2);
L = size(actionData, 2); % number of actions
V = size(InitialPairProb, 1);
NPoseClass = sqrt(size(InitialPairProb, 2));

ClassProb = InitialClassProb;
PairProb = InitialPairProb;

loglikelihood = zeros(maxIter,1);

P.c = [];
P.clg.sigma_x = [];
P.clg.sigma_y = [];
P.clg.sigma_angle = [];


init_zeros = zeros(1,K);
P.c = zeros(1,K);
P.clg = repmat(struct('mu_y',init_zeros,'sigma_y',init_zeros,'mu_x',init_zeros,'sigma_x',init_zeros,...
    'mu_angle',init_zeros,'sigma_angle',init_zeros,'theta',zeros(K,(size(poseData,3)+1)*size(poseData,3))),size(G,1),1);

% EM algorithm
for iter=1:maxIter
    
    % M-STEP to estimate parameters for Gaussians
    % Fill in P.c, the initial state prior probability (NOT the class probability as in PA8 and EM_cluster.m)
    % Fill in P.clg for each body part and each class
    % Make sure to choose the right parameterization based on G(i,1)
    % Hint: This part should be similar to your work from PA8 and EM_cluster.m
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    P.c = zeros(1,K);
    for i =1:L
        P.c = P.c + ClassProb(actionData(i).marg_ind(1),:);
    end
    P.c = P.c / L;
    
    for k = 1:K
        W = ClassProb(:,k);
        for i = 1:size(G,1)
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
    
    % M-STEP to estimate parameters for transition matrix
    % Fill in P.transMatrix, the transition matrix for states
    % P.transMatrix(i,j) is the probability of transitioning from state i to state j
    P.transMatrix = zeros(K,K);
    
    % Add Dirichlet prior based on size of poseData to avoid 0 probabilities
    P.transMatrix = P.transMatrix + size(PairProb,1) * .05;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    P.transMatrix = P.transMatrix + reshape(sum(PairProb),NPoseClass,NPoseClass);
    Denominator = repmat(sum(P.transMatrix,2),1,NPoseClass);
    P.transMatrix = P.transMatrix ./ Denominator;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % E-STEP preparation: compute the emission model factors (emission probabilities) in log space for each
    % of the poses in all actions = log( P(Pose | State) )
    % Hint: This part should be similar to (but NOT the same as) your code in EM_cluster.m
    
    logEmissionProb = zeros(N,K);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n = 1:N
        for k = 1:K
            mu = zeros(size(G,1), 3);
            sigma = zeros(size(G,1), 3);

            for i = 1:size(G,1)
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
            logEmissionProb(n,k) = sum(sum(logll_n_k)); % the difference with EM_cluster
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % E-STEP to compute expected sufficient statistics
    % ClassProb contains the conditional class probabilities for each pose in all actions
    % PairProb contains the expected sufficient statistics for the transition CPDs (pairwise transition probabilities)
    % Also compute log likelihood of dataset for this iteration
    % You should do inference and compute everything in log space, only converting to probability space at the end
    % Hint: You should use the logsumexp() function here to do probability normalization in log space to avoid numerical issues
    
    ClassProb = zeros(N,K);
    PairProb = zeros(V,K^2);
    loglikelihood(iter) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:L
        NPose = length(actionData(i).marg_ind);
        F = repmat(struct('var',[],'card',[],'val',[]), 1, 2*NPose);
        F(1).var = 1;    F(1).card = [K];    F(1).val = log(P.c);
        logTranMax = log(reshape(P.transMatrix',1,9));
        
        for j = 2:NPose
            F(j).var = [j, j-1];    F(j).card = [K,K];    F(j).val = logTranMax;
        end
        
        for j = 1:NPose
            F(j+NPose).var = [j];    F(j+NPose).card = [K];    F(j+NPose).val = logEmissionProb(actionData(i).marg_ind(j),:);
        end
        
        [M, PCalibrated] = ComputeExactMarginalsHMM(F);
        
        for j = 1:NPose
            ClassProb(actionData(i).marg_ind(j),:) = M(j).val;
        end
        
        for j = 1:NPose-1
            temp = PCalibrated.cliqueList(j).val;
            PairProb(actionData(i).pair_ind(j),:) = temp - logsumexp(temp);
        end
        
        loglikelihood(iter) = loglikelihood(iter) + logsumexp(PCalibrated.cliqueList(1).val);
    end

    ClassProb = exp(ClassProb); PairProb = exp(PairProb);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Print out loglikelihood
    disp(sprintf('EM iteration %d: log likelihood: %f', ...
        iter, loglikelihood(iter)));
    if exist('OCTAVE_VERSION')
        fflush(stdout);
    end
    
    % Check for overfitting by decreasing loglikelihood
    if iter > 1
        if loglikelihood(iter) < loglikelihood(iter-1)
            break;
        end
    end
    
end

% Remove iterations if we exited early
loglikelihood = loglikelihood(1:iter);
end
