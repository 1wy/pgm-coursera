% Copyright (C) Daphne Koller, Stanford University, 2012

function EUF = CalculateExpectedUtilityFactor( I )

  % Inputs: An influence diagram I with a single decision node and a single utility node.
  %         I.RandomFactors = list of factors for each random variable.  These are CPDs, with
  %              the child variable = D.var(1)
  %         I.DecisionFactors = factor for the decision node.
  %         I.UtilityFactors = list of factors representing conditional utilities.
  % Return value: A factor over the scope of the decision rule D from I that
  % gives the conditional utility given each assignment for D.var
  %
  % Note - We assume I has a single decision node and utility node.
  EUF = [];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % YOUR CODE HERE...
  %
%   for i = 1:I.DecisionFactors(1).card
%       decisianFac = I.DecisionFactors(1);
% 
%       F = [I.RandomFactors ];
%       U = I.UtilityFactors(1);
% 
%       VarF = unique([F(:).var]);
%       VarU = unique(U.var);
%       Z = setdiff(VarF, VarU);
%       Fnew = VariableElimination(F, Z);
% 
%       U.val = reshape(U.val, U.card);
%       [lia,lib]=ismember(U.var, Fnew.var);
%       U.val = permute(U.val,lib);
%       EUF(i) = sum(U.val(:) .* Fnew.val');
%   end
  F = [I.RandomFactors I.UtilityFactors];
  
  VarF = unique([F(:).var]);
  VarD = unique(I.DecisionFactors(1).var);
  Z = setdiff(VarF, VarD);
  Fnew = VariableElimination(F, Z);
  
  f = Fnew(1);
  for f_i = Fnew(2:end)
      f = FactorProduct(f, f_i);
  end

  EUF = f;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
end  
