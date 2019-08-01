% Copyright (C) Daphne Koller, Stanford University, 2012

function [MEU OptimalDecisionRule] = OptimizeLinearExpectations( I )
  % Inputs: An influence diagram I with a single decision node and one or more utility nodes.
  %         I.RandomFactors = list of factors for each random variable.  These are CPDs, with
  %              the child variable = D.var(1)
  %         I.DecisionFactors = factor for the decision node.
  %         I.UtilityFactors = list of factors representing conditional utilities.
  % Return value: the maximum expected utility of I and an optimal decision rule 
  % (represented again as a factor) that yields that expected utility.
  % You may assume that there is a unique optimal decision.
  %
  % This is similar to OptimizeMEU except that we will have to account for
  % multiple utility factors.  We will do this by calculating the expected
  % utility factors and combining them, then optimizing with respect to that
  % combined expected utility factor.  
  MEU = [];
  OptimalDecisionRule = [];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % YOUR CODE HERE
  %
  % A decision rule for D assigns, for each joint assignment to D's parents, 
  % probability 1 to the best option from the EUF for that joint assignment 
  % to D's parents, and 0 otherwise.  Note that when D has no parents, it is
  % a degenerate case we can handle separately for convenience.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  D = I.DecisionFactors(1);
  AllU = I.UtilityFactors;
  nU = length(AllU);
  EUFs = repmat(struct('var', [], 'card', [], 'val', []), 1,nU);
  for i = 1:nU
      
      I.UtilityFactors = AllU(i);
      
      EUFactor = CalculateExpectedUtilityFactor(I);
      EUFs(1, i) = EUFactor;
  end

  EUFactor = EUFs(1);
  for EUF_i = EUFs(2:end)
      EUFactor = FactorAdd(EUFactor, EUF_i);
  end
  
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
