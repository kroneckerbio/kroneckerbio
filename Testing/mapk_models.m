function [m, con, obs, obj, opts, objPrior] = mapk_models()
%% Load topologies
% Four models of the One-Step MAPK pathway
mDKDP = LoadModel('../Tutorial/Ferrell_MAPK_DKDP.txt');
mDKPP = LoadModel('../Tutorial/Ferrell_MAPK_DKPP.txt');
mPKDP = LoadModel('../Tutorial/Ferrell_MAPK_PKDP.txt');
mPKPP = LoadModel('../Tutorial/Ferrell_MAPK_PKPP.txt');

m = [mDKDP; mDKPP; mPKDP; mPKPP];

n_top = numel(m);

%% Options
% Global
opts = [];
opts.UseModelSeeds = false;
opts.UseModelInputs = false;
opts.RelTol = 1e-4;
opts.AbsTol = 1e-8;
opts.Normalized = true;
opts.Verbose = 0;

% Fitting
opts.LowerBound = 1e-10;
opts.UpperBound = 1e8;
opts.TolOptim = 1e-2;
opts.UseAdjoint = false;
opts.MaxIter = 200;
opts.Restart = 0;

% Topology probability
opts.NeedFit = true;

optsTop = opts;
optsTop.UseParams = {1:m(1).nk; 1:m(2).nk; 1:m(3).nk; 1:m(4).nk};

%% Do one simple experiment for the models to fit to
opts.populate = true;
outputs = [1;2;3];
tF = 100;
n = 8;
lintimes = linspace(tF/n,tF,n);
sd = sdLinear(0.01,0.1);

us = [1;2];
n_con = numel(us);

% Define measurement structure
outputlist = vec(repmat(outputs', [numel(lintimes),1]));
timelist = repmat(lintimes, [numel(outputs),1]);

con = experimentZero(2);
obs = observationZero([1,2]);
obj = objectiveZero([1,2]);
rand_state = rng(1);
for i_con = 1:n_con
    con(i_con) = experimentInitialValue(m(1), [], us(1), [], 'InitialValueExperiment');

    obs = observationLinearWeightedSumOfSquares(outputlist, timelist, sd, 'Fitting Data');
    sim = SimulateSystem(m(1), con(i_con), obs, opts);
    obj(i_con) = obs.Objective(sim.measurements);
end
rng(rand_state);

%% Define priors
% Topologies prior is uniform
opts.PriorTopology = zeros(n_top,1) + 1/n_top;

% Parameter prior is a spherical Gaussian
CI = 100; % Fold change to 95% CI
variance = (log(CI)/1.96).^2; % Log space variance

mux = cell(n_top,1);
Vx = cell(n_top,1);
for i = 1:n_top
    % Mean
    mux{i} = zeros(m(i).nk,1) + 1e-1;
    
    % Variance (normalized with no covariance)
    Vx{i} = diag(mux{i} .* (zeros(m(i).nk,1) + variance) .* mux{i});
end

%% Create prior objectives
objPrior = objectiveZero(n_top);
for i = 1:n_top
    objPrior(i) = objectiveLogNormalPriorOnKineticParameters(mux{i}, Vx{i});
end
