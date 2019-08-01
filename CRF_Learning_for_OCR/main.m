
% load Train1X.mat;
% load Train1Y.mat;
% load Validation1X.mat;
% load Validation1Y.mat;
% load Part1Lambdas.mat;
% load ValidationAccuracy;
% 
% theta = LRTrainSGD(Train1X, Train1Y, 0);
% prediction = sigmoid (Train1X * theta);
% idx = prediction>=0.5;
% prediction(idx) = 1;
% prediction(~idx) = 0;
% LRAccuracy(Train1Y, prediction)
% 
% allAcc = LRSearchLambdaSGD(Train1X, Train1Y, Validation1X, Validation1Y, Part1Lambdas);

load Part2Sample.mat;

[nll, grad] = InstanceNegLogLikelihood(sampleX, sampleY, sampleTheta, sampleModelParams);

a = 1;