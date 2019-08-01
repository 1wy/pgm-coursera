load('PA4Sample.mat');
% 
% P = ComputeInitialPotentials(InitPotential.INPUT);
% 
% [OUTPUT1, OUTPUT2] = GetNextCliques(GetNextC.INPUT1, GetNextC.INPUT2);
% 
% P = CliqueTreeCalibrate(SumProdCalibrate.INPUT, 0);
% 
% M = ComputeExactMarginalsBP(ExactMarginal.INPUT, [], 0);
% 
% % check answer
% load('PA4Sample.mat', 'SixPersonPedigree');
% i = 8;
% M1 = ComputeMarginal(i, SixPersonPedigree, []);
% disp(M1);
% M2 = ComputeExactMarginalsBP(SixPersonPedigree, [], 0);
% disp(M2(i));
% 
% B = FactorMaxMarginalization(FactorMax.INPUT1, FactorMax.INPUT2);
% P = CliqueTreeCalibrate(MaxSumCalibrate.INPUT, 1);
% M = ComputeExactMarginalsBP(MaxMarginals.INPUT, [], 1);
% A = MaxDecoding(MaxDecoded.INPUT);
% 
load('PA4Sample.mat', 'OCRNetworkToRun');
maxMarginals = ComputeExactMarginalsBP(OCRNetworkToRun, [], 1);
MAPAssignment = MaxDecoding(maxMarginals);
% DecodedMarginalsToChars(MAPAssignment)

load InitPotential;
P = ComputeInitialPotentials(InitPotential.INPUT);
a = 1;