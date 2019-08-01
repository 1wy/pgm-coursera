load PA3SampleCases.mat;
load PA3Models.mat;

factors = ComputeSingletonFactors(Part1SampleImagesInput, imageModel);

factors = ComputeEqualPairwiseFactors(Part2SampleImagesInput, imageModel.K);

factors = ComputePairwiseFactors(Part2SampleImagesInput, pairwiseModel, imageModel.K);

factors = ComputeTripletFactors(Part3SampleImagesInput, tripletList, imageModel.K);

factor = ComputeSimilarityFactor(Part4SampleImagesInput, 26, 1, 2);

factors = ComputeAllSimilarityFactors(Part5SampleImagesInput, 26);

factors = ChooseTopSimilarityFactors(factors, 2);

a = 1;