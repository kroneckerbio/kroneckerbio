%% Load various SBML models of increasing size
opts = [];
opts.Verbose = 2;

%% Simple model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m1 = LoadModelSbmlAnalytic('enzyme-catalysis-basic.xml', opts);
m1 = FinalizeModel(m1);

%% More complicated model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.UseNames = true;
m2 = LoadModelSbmlAnalytic('../Models/Brown_EGFNGF.xml', opts);

% Add test outputs
%   Note that outputs always use Names, not IDs since they are always added by
%   the user independent of SBML
m2 = AddOutput(m2, 'Out1', '("RasGapActive" + kSos*RapGapActive)/2 + sqrt(AktActive)^(kRap1ToBRaf)');
m2 = AddOutput(m2, 'Out2', 'EGF + 2*NGF');
m2 = FinalizeModel(m2);

%% Supposed to be mass action %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m3 = LoadModelSbmlAnalytic('../Models/Chen2009_ErbB_A431.xml');
m3 = FinalizeModel(m3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
