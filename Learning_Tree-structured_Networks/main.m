
load PA8SampleCases.mat;

% [mu sigma] = FitGaussianParameters(exampleINPUT.t1a1);
% [beta, sigma] = FitLinearGaussianParameters(exampleINPUT.t2a1, exampleINPUT.t2a2);

% load PA8Data.mat;
load PA8Data.mat;
% VisualizeDataset(trainData.data);

% logll = ComputeLogLikelihood(exampleINPUT.t3a1, exampleINPUT.t3a2, exampleINPUT.t3a3); 

[P1 L1] = LearnCPDsGivenGraph(exampleINPUT.t4a1, exampleINPUT.t4a2, exampleINPUT.t4a3);
% [P1 L1] = LearnCPDsGivenGraph(trainData.data, G1, trainData.labels);
% accuracy1 = ClassifyDataset(testData.data, testData.labels, P1, G1)
% 
% [P2 L2] = LearnCPDsGivenGraph(trainData.data, G2, trainData.labels);
% accuracy2 = ClassifyDataset(testData.data, testData.labels, P2, G2)
% [A W] = LearnGraphStructure(exampleINPUT.t6a1);
[P G L] = LearnGraphAndCPDs(exampleINPUT.t7a1, exampleINPUT.t7a2);
% acc = ClassifyDataset(exampleINPUT.t5a1, exampleINPUT.t5a2, exampleINPUT.t5a3, exampleINPUT.t5a4);
a = 1;