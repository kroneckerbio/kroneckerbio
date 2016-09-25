% Test analytic -> SimBio and SBML code
clear; close all; clc

% load model - expensive for testing so cache it
loadModel = true;
if loadModel
    % Loaded SBML model
    %     m = LoadModelSbmlAnalytic('test.xml');
    %     m = FinalizeModel(m);
    
    % MM model
    %     m = michaelis_menten_model();
    
    % Quoted identifiers model
%     m = InitializeModelAnalytic('QuotedIdentifiersModel');
%     m = AddCompartment(m, 'v', 3, 1);
%     m = AddCompartment(m, 'v!', 3, 1);
%     m = AddInput(m, 'E#', 'v', 1);
%     m = AddState(m, 'S', 'v', 'S0^2');
%     m = AddState(m, 'P:P', 'v', 'P0 + 10');
%     m = AddState(m, 'D', 'v!', 1.1);
%     m = AddParameter(m, 'K_m', 10);
%     m = AddParameter(m, 'kc@', 2);
%     m = AddSeed(m, 'S0', 5);
%     m = AddSeed(m, 'P0', 1.5);
%     m = AddReaction(m, 'r1', 'S', 'P:P', '"kc@"*"E#"*S/(K_m+S)');
%     m = AddReaction(m, 'r2', 'D', [], '"kc@"*D');
%     m = FinalizeModel(m);

    % Simple massaction model
%     m = LoadModelMassAction('Equilibrium.txt');
    
    % 
    % Note that the 'pump' species has a cell matrix for its InitialValue,
    %   but it's just set to a constant
    % TODO: Gives different results when simulating!!!
    m = LoadModelMassAction('Simple.txt');
    
    save('test_model_1_.mat', 'm');
end

loaded = load('test_model_1_.mat');
m = loaded.m;

%% Convert to SimBio
% simbio = ExportModelAnalyticSimBio(m);
simbio = ExportModelMassActionSimBio(m);

%% Simulate with SimBio tools
% Also open in the SimBio GUI to validate
cs = getconfigset(simbio);
cs.RuntimeOptions.StatesToLog = 'all';

result = sbiosimulate(simbio);

sbioplot(result);

%% Convert to SBML
filename = 'test_model_1_exported_.xml';
% ExportModelAnalyticSBML(m, filename);
ExportModelMassActionSBML(m, filename);

%% Load and see if it's close to what we started with
% Of course, the outputs are now rules and seeds are now params (which is fine in analytic models)
m2 = LoadModelSbmlAnalytic(filename);
% m2 = LoadModelSbmlMassAction(filename);
m2 = FinalizeModel(m2);

% Simulate
con2 = experimentInitialValue(m2, [], [], [], 'ExptLoaded');
tF = 10;
t = linspace(0, tF, 100);
sim2 = SimulateSystem(m2, con2, tF);

figure
plot(t, sim2.x(t))
legend({m2.States.Name})
title('Loaded Model States')
% TODO: How do we want to handle the outputs?

%% Compare to original model
con = experimentInitialValue(m, [], [], [], 'ExptOrig');
sim = SimulateSystem(m, con, tF);

figure
plot(t, sim.x(t))
legend({m.States.Name})
title('Original Model States')