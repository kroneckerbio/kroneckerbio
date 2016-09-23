% Test analytic -> SBML using libSBML
clear; close all; clc

% load model - expensive for testing so cache it
loadModel = true;
if loadModel
    % Loaded SBML model
%     m = LoadModelSbmlAnalytic('test.xml');
%     m = FinalizeModel(m);
    
    % MM model
    m = michaelis_menten_model();

    save('test_model_2_.mat', 'm');
end

loaded = load('test_model_2_.mat');
m = loaded.m;

%%
ExportModelAnalyticSBML2(m, 'test_model_2_exported_.xml');