% function [nll, grad] = InstanceNegLogLikelihood(X, y, theta, modelParams)
% returns the negative log-likelihood and its gradient, given a CRF with parameters theta,
% on data (X, y).
%
% Inputs:
% X            Data.                           (numCharacters x numImageFeatures matrix)
%              X(:,1) is all ones, i.e., it encodes the intercept/bias term.
% y            Data labels.                    (numCharacters x 1 vector)
% theta        CRF weights/parameters.         (numParams x 1 vector)
%              These are shared among the various singleton / pairwise features.
% modelParams  Struct with three fields:
%   .numHiddenStates     in our case, set to 26 (26 possible characters)
%   .numObservedStates   in our case, set to 2  (each pixel is either on or off)
%   .lambda              the regularization parameter lambda
%
% Outputs:
% nll          Negative log-likelihood of the data.    (scalar)
% grad         Gradient of nll with respect to theta   (numParams x 1 vector)
%
% Copyright (C) Daphne Koller, Stanford Univerity, 2012

function [nll, grad] = InstanceNegLogLikelihood(X, y, theta, modelParams)

% featureSet is a struct with two fields:
%    .numParams - the number of parameters in the CRF (this is not numImageFeatures
%                 nor numFeatures, because of parameter sharing)
%    .features  - an array comprising the features in the CRF.
%
% Each feature is a binary indicator variable, represented by a struct
% with three fields:
%    .var          - a vector containing the variables in the scope of this feature
%    .assignment   - the assignment that this indicator variable corresponds to
%    .paramIdx     - the index in theta that this feature corresponds to
%
% For example, if we have:
%
%   feature = struct('var', [2 3], 'assignment', [5 6], 'paramIdx', 8);
%
% then feature is an indicator function over X_2 and X_3, which takes on a value of 1
% if X_2 = 5 and X_3 = 6 (which would be 'e' and 'f'), and 0 otherwise.
% Its contribution to the log-likelihood would be theta(8) if it's 1, and 0 otherwise.
%
% If you're interested in the implementation details of CRFs,
% feel free to read through GenerateAllFeatures.m and the functions it calls!
% For the purposes of this assignment, though, you don't
% have to understand how this code works. (It's complicated.)

featureSet = GenerateAllFeatures(X, modelParams);

% Use the featureSet to calculate nll and grad.
% This is the main part of the assignment, and it is very tricky - be careful!
% You might want to code up your own numerical gradient checker to make sure
% your answers are correct.
%
% Hint: you can use CliqueTreeCalibrate to calculate logZ effectively.
%       We have halfway-modified CliqueTreeCalibrate; complete our implementation
%       if you want to use it to compute logZ.

nll = 0;
grad = zeros(size(theta));
%%%
% Your code here:
load Part2Sample.mat;
%     % convert features to Factors
varDict = {repmat('', length(featureSet.features), 1)};
for i = 1:length(featureSet.features)
    varDict(i) = {num2str(featureSet.features(i).var)};
end
varDict = unique(varDict);
F = repmat(struct('var',[],'card', [], 'val',[]), length(varDict), 1);
for i = 1:length(varDict)
    F(i).var = str2num(varDict{i});
    F(i).card = repmat(modelParams.numHiddenStates, 1, length(F(i).var));
    F(i).val = zeros(1, prod(F(i).card));
end
for fea = featureSet.features
    idx = find(strcmp(varDict, num2str(fea.var)), 1);
    assignIdx = AssignmentToIndex(fea.assignment, F(idx).card);
    F(idx).val(assignIdx) = F(idx).val(assignIdx) + theta(fea.paramIdx);
end
for i = 1:length(F)
    %         F(i).val = F(i).val - max(F(i).val);
    F(i).val = exp(F(i).val);
end
%
P = CreateCliqueTree(F);
[P, logZ] = CliqueTreeCalibrate(P, 0);

fContrib = 0;
for fea = featureSet.features
    if all(y(fea.var) == fea.assignment)
        fContrib = fContrib + theta(fea.paramIdx);
    end
end

regContrib = modelParams.lambda / 2 * (theta * theta');
nll = logZ - fContrib + regContrib;

%
for f=1:length(F)
    for i=1:length(P.cliqueList)
        if all(ismember(F(f).var, P.cliqueList(i).var))
            V = setdiff(P.cliqueList(i).var, F(f).var);
            F(f) = FactorMarginalization(P.cliqueList(i), V);
            F(f).val = F(f).val./sum(F(f).val);
            break;
        end
    end;
end;
% compute the gradient
featureCount = zeros(size(theta));
modelFeatureCount = zeros(size(theta));

for fea = featureSet.features
    for f = 1:length(F)
        if strcmp(num2str(F(f).var), num2str(fea.var))
            idx = AssignmentToIndex(fea.assignment, F(f).card);
            modelFeatureCount(fea.paramIdx) = modelFeatureCount(fea.paramIdx) + F(f).val(idx);
            break;
        end
    end
    if all(y(fea.var) == fea.assignment)
        featureCount(fea.paramIdx) = featureCount(fea.paramIdx) + 1;
    end
end

regContribGrad = modelParams.lambda * theta;
grad = modelFeatureCount - featureCount + regContribGrad;

end
