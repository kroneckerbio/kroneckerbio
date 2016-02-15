%% Load various SBML and SimBiology models
% Important Note: Remember to call FinalizeModel after loading and before
% running kroneckerbio analysis functions.

opts = [];
opts.Verbose = 2;

%% Simple model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1 = LoadModelSbmlAnalytic('../Testing/enzyme-catalysis-basic.xml', opts);
m1 = FinalizeModel(m1);

%% More complicated model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.UseNames = true;
m2 = LoadModelSbmlAnalytic('Brown_EGFNGF.xml', opts);

% Add test outputs
%   Note that outputs always use Names, not IDs since they are always added by
%   the user independent of SBML
m2 = AddOutput(m2, 'Out1', '("RasGapActive" + kSos*RapGapActive)/2 + sqrt(AktActive)^(kRap1ToBRaf)');
m2 = AddOutput(m2, 'Out2', 'EGF + 2*NGF');
m2 = FinalizeModel(m2, opts);

%% Simple mass action model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = [];
opts.Verbose = 2;
m3 = LoadModelSbmlMassAction('../Testing/simple_massaction.xml', opts);
m3 = FinalizeModel(m3);
% m3 = LoadModelSbmlAnalytic('simple_massaction.xml', opts);
% m3 = FinalizeModel(m3, opts);

%% Bigger mass action model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warning: loading as a massaction model is slow
m4 = LoadModelSbmlAnalytic('Chen2009_ErbB_A431.xml', opts);
m4 = FinalizeModel(m4, opts);
% m4 = LoadModelSbmlMassAction('../Models/Chen2009_ErbB_A431.xml', opts);
% m4 = FinalizeModel(m4);

%% Load SimBiology model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../Testing/simple_massaction_simbio_model.mat')

% Test SimBio -> analytic model
m5 = LoadModelSimBioAnalytic(simbiomodel);
m5 = FinalizeModel(m5);

% Test SimBio -> massaction model
m6 = LoadModelSimBioMassAction(simbiomodel);
m6 = FinalizeModel(m5);

