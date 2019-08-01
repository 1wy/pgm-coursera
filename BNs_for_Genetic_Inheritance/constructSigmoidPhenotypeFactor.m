function phenotypeFactor = constructSigmoidPhenotypeFactor(alleleWeights, geneCopyVarOneList, geneCopyVarTwoList, phenotypeVar)
% This function takes a cell array of alleles' weights and constructs a 
% factor expressing a sigmoid CPD.
%
% You can assume that there are only 2 genes involved in the CPD.
%
% In the factor, for each gene, each allele assignment maps to the allele
% whose weight is at the corresponding location.  For example, for gene 1,
% allele assignment 1 maps to the allele whose weight is at
% alleleWeights{1}(1) (same as w_1^1), allele assignment 2 maps to the
% allele whose weight is at alleleWeights{1}(2) (same as w_2^1),....  
% 
% You may assume that there are 2 possible phenotypes.
% For the phenotypes, assignment 1 maps to having the physical trait, and
% assignment 2 maps to not having the physical trait.
%
% THE VARIABLE TO THE LEFT OF THE CONDITIONING BAR MUST BE THE FIRST
% VARIABLE IN THE .var FIELD FOR GRADING PURPOSES
%
% Input:
%   alleleWeights: Cell array of weights, where each entry is an 1 x n 
%   of weights for the alleles for a gene (n is the number of alleles for
%   the gene)
%   geneCopyVarOneList: m x 1 vector (m is the number of genes) of variable 
%   numbers that are the variable numbers for each of the first parent's 
%   copy of each gene (numbers in this list go in the .var part of the
%   factor)
%   geneCopyVarTwoList: m x 1 vector (m is the number of genes) of variable 
%   numbers that are the variable numbers for each of the second parent's 
%   copy of each gene (numbers in this list go in the .var part of the
%   factor) -- Note that both copies of each gene are from the same person,
%   but each copy originally came from a different parent
%   phenotypeVar: Variable number corresponding to the variable for the 
%   phenotype (goes in the .var part of the factor)
%
% Output:
%   phenotypeFactor: Factor in which the values are the probabilities of 
%   having each phenotype for each allele combination (note that this is 
%   the FULL CPD with no evidence observed)

phenotypeFactor = struct('var', [], 'card', [], 'val', []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSERT YOUR CODE HERE
% Note that computeSigmoid.m will be useful for this function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% Fill in phenotypeFactor.var.  This should be a 1-D row vector.
% Fill in phenotypeFactor.card.  This should be a 1-D row vector.
phenotypeFactor.var = [phenotypeVar, geneCopyVarOneList', geneCopyVarTwoList'];
numAlleles = zeros(1, length(alleleWeights));
idx = cell(1, length(alleleWeights));
for i = 1:length(alleleWeights)
    numAlleles(i) = length(alleleWeights{1});
    idx{i} = 1:numAlleles(i);
end
phenotypeFactor.card = [2, numAlleles, numAlleles];
idx = [idx, idx];
comb = cartesian(idx);

table = zeros(1, length(comb));
% for i = 1:length(comb)
%     accumlateInfluence = 0;
%     for j = 1:numAlleles
%         wgt = alleleWeights{j};
%         accumlateInfluence = accumlateInfluence + wgt(comb(i,j)) + wgt(comb(i,j+length(alleleWeights)));
%     end
%     table(i) = computeSigmoid(accumlateInfluence);
% end
contrib = zeros(size(comb));
for g = 1:length(alleleWeights)
    wgt = alleleWeights{g};
    contrib(:,g) = wgt(comb(:,g));
    contrib(:,g+length(alleleWeights)) = wgt(comb(:,g+length(alleleWeights)));
end
table = computeSigmoid(sum(contrib, 2));
% table(2,:) = 1 - table(1,:);
% table = table';
% table = table(:);
phenotypeFactor.val = zeros(1, prod(phenotypeFactor.card));
phenotypeFactor.val(1:2:end) = table;
phenotypeFactor.val(2:2:end) = 1-table;
% Replace the zeros in phentoypeFactor.val with the correct values.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function G = cartesian(varargin)
    args = varargin{1};
    n = length(args);

    [F{1:n}] = ndgrid(args{:});

    for i=n:-1:1
        G(:,i) = F{i}(:);
    end

%     C = unique(G , 'rows');
end