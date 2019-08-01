
load PA9SampleCases.mat;
load PA9Data.mat;
% [P ll CP] = EM_cluster2(exampleINPUT.t1a1, exampleINPUT.t1a2, exampleINPUT.t1a3, exampleINPUT.t1a4);
% [P loglikelihood ClassProb] = EM_cluster(poseData1, G, InitialClassProb2, 20);
% [maxprob, assignments] = max(ClassProb, [], 2);
% VisualizeDataset(poseData1(assignments == 1, :, :));
[P ll CP PP] = EM_HMM2(exampleINPUT.t2a1, exampleINPUT.t2a2, exampleINPUT.t2a3, exampleINPUT.t2a4, exampleINPUT.t2a5, exampleINPUT.t2a6);

a = 1;