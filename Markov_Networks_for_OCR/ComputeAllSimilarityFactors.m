function factors = ComputeAllSimilarityFactors (images, K)
% This function computes all of the similarity factors for the images in
% one word.
%
% Input:
%   images: An array of structs containing the 'img' value for each
%     character in the word.
%   K: The alphabet size (accessible in imageModel.K for the provided
%     imageModel).
%
% Output:
%   factors: Every similarity factor in the word. You should use
%     ComputeSimilarityFactor to compute these.
%
% Copyright (C) Daphne Koller, Stanford University, 2012

n = length(images);
nFactors = nchoosek (n, 2);

factors = repmat(struct('var', [], 'card', [], 'val', []), nFactors, 1);

[x, y] = meshgrid(1:n, 1:n);
comb = [x(:), y(:)];
comb = comb(comb(:,1) < comb(:,2),:);
% Your code here:
for i = 1:length(comb)
    factors(i) = ComputeSimilarityFactor(images, K, comb(i,1), comb(i,2));
end
end

