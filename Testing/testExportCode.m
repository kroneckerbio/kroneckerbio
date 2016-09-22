% Test analytic -> SimBio and SBML code
clear; close all; clc

% load model - expensive for testing so cache it
loadModel = true;
if loadModel
    % Loaded SBML model
%     m = LoadModelSbmlAnalytic('test.xml');
%     m = FinalizeModel(m);
    
    % MM model
    m = michaelis_menten_model();

    save('test_model_1_.mat', 'm');
end

loaded = load('test_model_1_.mat');
m = loaded.m;

%% Convert to SimBio
simbio = ExportModelAnalyticSimBio(m);

%% Simulate with SimBio tools
% Also open in the SimBio GUI to validate
cs = getconfigset(simbio);
cs.RuntimeOptions.StatesToLog = 'all';

result = sbiosimulate(simbio);

sbioplot(result);

%% Convert to SBML
filename = 'test_model_1_exported_.xml';
ExportModelAnalyticSBML(m, filename);

%% Load and see if it's close to what we started with
% Of course, the outputs are now rules and seeds are now params (which is fine in analytic models)
m2 = LoadModelSbmlAnalytic(filename);
m2 = FinalizeModel(m2);

% Simulate
con = experimentInitialValue(m2, [], [], [], 'InitialValueExperiment');
tF = 10;
t = linspace(0, tF, 100);
sim1 = SimulateSystem(m2, con, tF);

figure
plot(t, sim1.x(t))
ylim([0,25])
% TODO: How do we want to handle the outputs?
