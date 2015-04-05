% Global optimization example
%   Could use a better test system with known multiple local minima
clear; close all; clc

%% Construct equilibrium experiment A + B <-> C
m = InitializeModel('Equilibrium');

m = AddCompartment(m, 'Solution', 3, 1);

m = AddState(m, 'A', 'Solution', 1);
m = AddState(m, 'B', 'Solution', 2);
m = AddState(m, 'C', 'Solution', 0);

m = AddOutput(m, 'A', 'A');
m = AddOutput(m, 'B', 'B');
m = AddOutput(m, 'C', 'C');

m = AddParameter(m, 'kf', 5);
m = AddParameter(m, 'kr', 3);

m = AddReaction(m, '', '', 'A', 'B', 'C', '', 'kf', 'kr');

m = FinalizeModel(m);

%% Construct experiment
con = InitialValueExperiment(m, 1, [], [], [], 'InitialValueExperiment');

%% Simulate
sim = Simulate(m, con);

figure
plot(sim.sol.x', sim.sol.y')
legend('A','B','C','Location','best')
title('Simulated System')

%% Generate sample data to fit to
nTimes = 10;
times = linspace(0, 1, nTimes);
data = zeros(3, nTimes);
for i = 1:3
    data(i,:) = max(sim.x(times, i) + normrnd(0, 0.1, [1,nTimes]), 0);
end

figure
set(0,'DefaultAxesColorOrder',[0    0.4470    0.7410; 0.8500    0.3250    0.0980; 0.9290    0.6940    0.1250]);
plot(sim.sol.x', sim.sol.y', times', data', '+')
legend('A','B','C','AGen','BGen','CGen','Location','best')
title('Simulated System with Randomly Generated Datapoints')

%% Construct objective
outputlist = [ones(nTimes,1); ones(nTimes,1)*2; ones(nTimes,1)*3];
timeslist = repmat(times, 1, 3)';
dataCols = data';
measurements = dataCols(:);
sd = sdLinear(0.05, 0.1);

obj = objectiveWeightedSumOfSquaresNonNeg(outputlist, timeslist, sd, measurements);

%% Local Fit
mFit = FitObjective(m, con, obj);

fprintf('Local optimization results:\n')
for i = 1:m.nk
    fprintf('%s\tActual: %5.3f\tFit: %5.3f\n', m.Parameters(i).Name, m.Parameters(i).Value, mFit.Parameters(i).Value)
end

simLocalFit = Simulate(mFit, con);
figure
plot(sim.sol.x', sim.sol.y', times', data', '+', simLocalFit.sol.x', simLocalFit.sol.y', '--')
legend('A','B','C', ...
    'AGen','BGen','CGen', ...
    'AFit','BFit','CFit','Location','best')
title('Local Fit Comparison')

%% Linearized parameter uncertainty
% F = ObjectiveInformation(mfit, con, obj);

%% Global Fit
opts = [];
opts.Verbose = 1;
opts.LowerBounds = [1e-1, 1e-1];
opts.UpperBounds = [1e2,  1e2 ];
opts.GlobalOptimization = true;
% opts.GlobalOpts.Algorithm = 'multistart';
% opts.GlobalOpts.Algorithm = 'patternsearch';
% opts.GlobalOpts.UseParallel = true; % Can't currently run multistart in parallel due to global vars in obj funs

mfitGlobal = FitObjective(m, con, obj, opts);

fprintf('Global optimization results:\n')
for i = 1:m.nk
    fprintf('%s\tActual: %5.3f\tFit: %5.3f\n', m.Parameters(i).Name, m.Parameters(i).Value, mfitGlobal.Parameters(i).Value)
end

simGlobal = Simulate(mfitGlobal, con);
figure
plot(sim.sol.x', sim.sol.y', times', data', '+', simLocalFit.sol.x', simLocalFit.sol.y', '--')
legend('A','B','C', ...
    'AGen','BGen','CGen', ...
    'AFit','BFit','CFit','Location','best')
title('Global Fit Comparison')

