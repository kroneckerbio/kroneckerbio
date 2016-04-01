% Load the One-Step MAPK suite
%% Options
% Global
opts.Verbose = 1;

% Integration
opts.RelTol = 1e-4;
opts.AbsTol = 1e-7;
opts.UseModelSeeds = true;
opts.UseModelInputs = false;

% Parameters
opts.UseAdjoint = false;
opts.Normalized = true;

% Fitting
opts.LowerBound = 1e-5;
opts.UpperBound = 1e3;
opts.TolOptim = 1e-4;
opts.MaxIter = 200;
opts.Restart = 0;

boundmean = 1e-1;
seed = 1337;

%% Load topologies
% Four models of the One-Step MAPK pathway
mDKDP = LoadModel('Ferrell_MAPK_DKDP.txt');
mDKPP = LoadModel('Ferrell_MAPK_DKPP.txt');
mPKDP = LoadModel('Ferrell_MAPK_PKDP.txt');
mPKPP = LoadModel('Ferrell_MAPK_PKPP.txt');

m = [mDKDP; mDKPP; mPKDP; mPKPP];
nTop = numel(m);

%% Do one simple experiment for the models to fit to
outputs = [1;2;3];
tF = 10000;
n = 10;
lintimes = linspace(tF/n,tF,n);
sd = sdLinear(0.01, 0.1);

con = Experiment(m(1), tF);

% Define measurement structure
outputlist = vec(repmat(outputs', [numel(lintimes),1]));
timelist = repmat(lintimes, [numel(outputs),1]);

obj = objectiveWeightedSumOfSquaresNonNeg(outputlist, timelist, sd, [], 'Fitting Data');

% Create fake data
sim = SimulateSelect(m(1), con, lintimes, opts);
rand_state = rng(1);
obj = obj.AddData(sim.sol);
rng(rand_state);

clear n

%% Fit all models to experiment with tight bounds to get true models
trueOpts = opts;
trueOpts.LowerBound = boundmean / 10;
trueOpts.UpperBound = boundmean * 10;

for i = 1:nTop
    % Scramble parameters
    m(i) = m(i).Update(opts.LowerBound + zeros(m(i).nk,1), m(i).s, m(i).q);
    % FitObjective
	m(i) = FitObjective(m(i), con, obj, trueOpts);
end

clear i trueOpts
