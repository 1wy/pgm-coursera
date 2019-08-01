%COMPUTEINITIALPOTENTIALS Sets up the cliques in the clique tree that is
%passed in as a parameter.
%
%   P = COMPUTEINITIALPOTENTIALS(C) Takes the clique tree skeleton C which is a
%   struct with three fields:
%   - nodes: cell array representing the cliques in the tree.
%   - edges: represents the adjacency matrix of the tree.
%   - factorList: represents the list of factors that were used to build
%   the tree. 
%   
%   It returns the standard form of a clique tree P that we will use through 
%   the rest of the assigment. P is struct with two fields:
%   - cliqueList: represents an array of cliques with appropriate factors 
%   from factorList assigned to each clique. Where the .val of each clique
%   is initialized to the initial potential of that clique.
%   - edges: represents the adjacency matrix of the tree. 
%
% Copyright (C) Daphne Koller, Stanford University, 2012


function P = ComputeInitialPotentials(C)

% number of cliques
N = length(C.nodes);

% initialize cluster potentials 
P.cliqueList = repmat(struct('var', [], 'card', [], 'val', []), N, 1);
P.edges = C.edges;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YOUR CODE HERE
%
% First, compute an assignment of factors from factorList to cliques. 
% Then use that assignment to initialize the cliques in cliqueList to 
% their initial potentials. 

% C.nodes is a list of cliques.
% So in your code, you should start with: P.cliqueList(i).var = C.nodes{i};
% Print out C to get a better understanding of its structure.
%

Nfac = length(C.factorList);
Vars = [];
for i = 1:Nfac
    Vars = union(Vars, C.factorList(i).var);
end
clique_assign = zeros(Nfac, N);
choice = cell(Nfac,1);
hasVar = cell(N, 1);
card = zeros(1, Vars(end));
% assign the factor that has only one choice
% for i = 1:Nfac
%     cho = [];
%     card(C.factorList(i).var) = C.factorList(i).card; 
%     for j = 1:N
%         is_mem = ismember(C.factorList(i).var, C.nodes{j});
%         if all(is_mem) 
%             cho = [cho j]; 
%         end
%     end
%     choice{i} = cho;
%     if length(cho) == 1
%        clique_assign(i, cho(1)) = 1;
%        hasVar{cho(1)} = union(hasVar{cho(1)}, C.factorList(i).var);
%     end
% end

% for i = 1:Nfac
%     if any(clique_assign(i,:)) == 0
%         assign = choice{i}(1);
%         for j = choice{i}(2:end)
%             if sum(clique_assign(:,j)) < sum(clique_assign(:,assign))
%                 assign = j;
%             end
%         end
%         clique_assign(i,assign) = 1;
%     end
% end

% assign the factor that has only one choice
for i = 1:Nfac
    card(C.factorList(i).var) = C.factorList(i).card; 
    for j = 1:N
        is_mem = ismember(C.factorList(i).var, C.nodes{j});
        if all(is_mem)  
            clique_assign(i,j) = 1;
            break;
        end
    end
end

for i = 1:N
    idx = find(clique_assign(:,i))';
    if ~isempty(idx)
    %     P.cliqueList(i) = struct('var', C.nodes{i}, 'card', );
        factor = C.factorList(idx(1));
        for id = idx(2:end)
            factor = FactorProduct(factor, C.factorList(id));
        end
        remainAddVar = setdiff(C.nodes{i},factor.var);
        if ~isempty(remainAddVar)
            Nrep = prod(card(remainAddVar));
            factor.var = [factor.var, remainAddVar];
            factor.card = [factor.card, card(remainAddVar)];
            factor.val = repmat(factor.val, 1, Nrep);
        end
        
        % adjust the order to be consistent with that of node
        if ~isempty(setdiff(factor.var, C.nodes{i})) || ~isempty(setdiff(C.nodes{i},factor.var))
            a = 1;
        end
        factor.val = reshape(factor.val, factor.card);
        [lia,lib]=ismember(C.nodes{i}, factor.var);
%         fprintf('%d    =================================\n', i);
%         disp(factor.var);
%         disp(C.nodes{i});
        factor.var = C.nodes{i};
        factor.card = factor.card(lib);
        factor.val = permute(factor.val,lib);
        factor.val = factor.val(:)';
    else
        factor.var = C.nodes{i};
        factor.card = card(C.nodes{i});
        factor.val = ones(1,prod(factor.card));
    end
    P.cliqueList(i) = factor;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

