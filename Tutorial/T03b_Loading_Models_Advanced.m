%% T03b Advanced Loading Models
% Load various SBML and SimBiology models
%
% Important Note: Remember to call FinalizeModel after loading and before
%   running kroneckerbio analysis functions.

%% Simple model
current_path = fileparts(mfilename('fullpath'));
opts = [];
opts.Verbose = 0;
m1 = LoadModelSbmlAnalytic(fullfile(current_path, '../Testing/enzyme-catalysis-basic.xml'), opts);
m1 = FinalizeModel(m1);

%% More complicated model
opts.ICsAsSeeds = true;
m2 = LoadModelSbmlAnalytic(fullfile(current_path, 'Brown_EGFNGF.xml'), opts);

% Add test outputs
%   Note that outputs always use Names, not IDs since they are always added by
%   the user independent of SBML
m2 = AddOutput(m2, 'Out1', '("RasGapActive" + kSos*RapGapActive)/2 + sqrt(AktActive)^(kRap1ToBRaf)');
m2 = AddOutput(m2, 'Out2', 'EGF + 2*NGF');
m2 = FinalizeModel(m2, opts);

%% Simple mass action model
opts = [];
opts.Verbose = 0;
m3 = LoadModelSbmlMassAction(fullfile(current_path, '../Testing/simple_massaction.xml'), opts);
m3 = FinalizeModel(m3);
% m3 = LoadModelSbmlAnalytic(fullfile(current_path, '../Testing/simple_massaction.xml'), opts);
% m3 = FinalizeModel(m3, opts);

%% Bigger mass action model
% Warning: loading as a massaction model is slow
%
% Note: When using libSBML's Matlab bindings TranslateSBML, an interactive
% prompt will appear asking whether to discard warnings. These prompts will
% be suppressed for opts.Verbose <= 1, which is useful for scripts.
opts = [];
opts.Verbose = 0;
m4 = LoadModelSbmlAnalytic(fullfile(current_path, 'Chen2009_ErbB_A431.xml'), opts);
m4 = FinalizeModel(m4, opts);
% m4 = LoadModelSbmlMassAction(fullfile(current_path, 'Chen2009_ErbB_A431.xml'), opts);
% m4 = FinalizeModel(m4);

%% Load SimBiology model
load(fullfile(current_path, '../Testing/simple_massaction_simbio_model.mat'))

% Test SimBio -> analytic model
m5 = LoadModelSimBioAnalytic(simbiomodel);
m5 = FinalizeModel(m5);

% Test SimBio -> massaction model
m6 = LoadModelSimBioMassAction(simbiomodel);
m6 = FinalizeModel(m5);

