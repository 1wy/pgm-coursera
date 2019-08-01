% Copyright (C) Daphne Koller, Stanford University, 2012

function [MEU OptimalDecisionRule] = OptimizeMEU( I )

  % Inputs: An influence diagram I with a single decision node and a single utility node.
  %         I.RandomFactors = list of factors for each random variable.  These are CPDs, with
  %              the child variable = D.var(1)
  %         I.DecisionFactors = factor for the decision node.
  %         I.UtilityFactors = list of factors representing conditional utilities.
  % Return value: the maximum expected utility of I and an optimal decision rule 
  % (represented again as a factor) that yields that expected utility.
  
  % We assume I has a single decision node.
  % You may assume that there is a unique optimal decision.
  D = I.DecisionFactors(1);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % YOUR CODE HERE...
  % 
  % Some other information that might be useful for some implementations
  % (note that there are multiple ways to implement this):
  % 1.  It is probably easiest to think of two cases - D has parents and D 
  %     has no parents.
  % 2.  You may find the Matlab/Octave function setdiff useful.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  EUFactor = CalculateExpectedUtilityFactor(I);
  OptimalDecisionRule = D;
  OptimalDecisionRule.val = zeros(size(D.val));
  
  if length(D.var) == 1 % without parents
    [MEU, max_index] = max(EUFactor.val);
    OptimalDecisionRule.val(max_index) = 1;
  else % with parents
    MEU = 0;
    card = EUFactor.card(1);
    numCondition = length(EUFactor.val) / card;
    for i = 1:numCondition
        [MEU_i, max_index] = max(EUFactor.val(1+(i-1)*card:i*card));
        MEU = MEU + MEU_i;
        OptimalDecisionRule.val(max_index + (i-1)*card) = 1;
    end
    
  end


end
