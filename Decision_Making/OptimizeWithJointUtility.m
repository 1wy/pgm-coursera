% Copyright (C) Daphne Koller, Stanford University, 2012

function [MEU OptimalDecisionRule] = OptimizeWithJointUtility( I )
  % Inputs: An influence diagram I with a single decision node and one or more utility nodes.
  %         I.RandomFactors = list of factors for each random variable.  These are CPDs, with
  %              the child variable = D.var(1)
  %         I.DecisionFactors = factor for the decision node.
  %         I.UtilityFactors = list of factors representing conditional utilities.
  % Return value: the maximum expected utility of I and an optimal decision rule 
  % (represented again as a factor) that yields that expected utility.
  % You may assume that there is a unique optimal decision.
    
  % This is similar to OptimizeMEU except that we must find a way to 
  % combine the multiple utility factors.  Note: This can be done with very
  % little code.
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % YOUR CODE HERE
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  D = I.DecisionFactors(1);
  U = I.UtilityFactors(1);
  for U_i = I.UtilityFactors(2:end)
      U = FactorAdd(U, U_i);
  end
  I.UtilityFactors = U;
 
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
