%CLIQUETREECALIBRATE Performs sum-product or max-product algorithm for 
%clique tree calibration.

%   P = CLIQUETREECALIBRATE(P, isMax) calibrates a given clique tree, P 
%   according to the value of isMax flag. If isMax is 1, it uses max-sum
%   message passing, otherwise uses sum-product. This function 
%   returns the clique tree where the .val for each clique in .cliqueList
%   is set to the final calibrated potentials.
%
% Copyright (C) Daphne Koller, Stanford University, 2012

function P = CliqueTreeCalibrate(P, isMax)

% Number of cliques in the tree.
N = length(P.cliqueList);

if isMax
    for i = 1:N
        P.cliqueList(i).val = log(P.cliqueList(i).val);
    end
end
% Setting up the messages that will be passed.
% MESSAGES(i,j) represents the message going from clique i to clique j. 
MESSAGES = repmat(struct('var', [], 'card', [], 'val', []), N, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We have split the coding part for this function in two chunks with
% specific comments. This will make implementation much easier.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% YOUR CODE HERE
% While there are ready cliques to pass messages between, keep passing
% messages. Use GetNextCliques to find cliques to pass messages between.
% Once you have clique i that is ready to send message to clique
% j, compute the message and put it in MESSAGES(i,j).
% Remember that you only need an upward pass and a downward pass.
%
Nsend = sum(sum(P.edges));
cnt = 0;
while cnt < Nsend
    [i, j] = GetNextCliques(P, MESSAGES);
    phi = P.cliqueList(i);
    neighbor = setdiff(find(P.edges(:,i)), j)';
    if ~isempty(neighbor)
        for ne = neighbor
            if isMax
                phi = FactorSum(phi, MESSAGES(ne,i));
            else
                phi = FactorProduct(phi, MESSAGES(ne,i));
            end
        end
    end    
    V = setdiff(P.cliqueList(i).var, P.cliqueList(j).var);
    if isMax
        msg = FactorMaxMarginalization(phi, V);
    else
%     normalize the message, my understanding is that clique represents
%     joint probabilities, so sum(msg.val) equals to 1.
        msg = FactorMarginalization(phi, V);
        msg.val = msg.val / sum(msg.val);
    end
    MESSAGES(i,j) = msg;
    cnt = cnt+1;
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%
% Now the clique tree has been calibrated. 
% Compute the final potentials for the cliques and place them in P.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:N
    neighbor = find(P.edges(:,i))';
    beta = P.cliqueList(i);
    if ~isempty(neighbor)
       for ne = neighbor
            if isMax
                beta = FactorSum(beta, MESSAGES(ne,i));
            else
                beta = FactorProduct(beta, MESSAGES(ne,i));
            end
       end
    end
    P.cliqueList(i) = beta;
end
return
