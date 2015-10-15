% TODO: currently outdated - fix

% Global
opts = [];
opts.UseModelSeeds = false;
opts.UseModelInputs = false;
opts.RelTol = 1e-4;
opts.AbsTol = 1e-8;
opts.Normalized = true;
opts.Verbose = 1;

% Fitting
opts.LowerBound = 1e-10;
opts.UpperBound = 1e8;
opts.TolOptim = 1e-6;
opts.UseAdjoint = false;
opts.MaxIter = 200;
opts.Restart = 0;

% Topology probability
opts.NeedFit = false;

%% Load topologies
% Four models of the One-Step MAPK pathway
mDKDP = LoadModel('../Models/Ferrell_MAPK_DKDP.txt');
mDKPP = LoadModel('../Models/Ferrell_MAPK_DKPP.txt');
mPKDP = LoadModel('../Models/Ferrell_MAPK_PKDP.txt');
mPKPP = LoadModel('../Models/Ferrell_MAPK_PKPP.txt');

m = [mDKDP; mDKPP; mPKDP; mPKPP];

nTop = numel(m);

%% Do one simple experiment for the models to fit to
opts.populate = true;
outputs = [1;2;3];
tF = 100;
n = 8;
lintimes = linspace(tF/n,tF,n);
sd = sdLinear(0.01,0.1);

con = experimentInitialValue(mDKDP, [], [], [], 'InitialValueExperiment');

% Define measurement structure
outputlist = vec(repmat(outputs', [numel(lintimes),1]));
timelist = repmat(lintimes, [numel(outputs),1]);

% obj = objectiveWeightedSumOfSquaresNonNeg(outputlist, timelist, sd, [], 'Fitting Data');
% obs = observationLinearWeightedSumOfSquares(outputlist, timelist, sd, 'Fitting Data');
% obj = obs.Objective([]);

% Create test data
sims = cell(1,nTop);
for i = 1:nTop
    sims{i} = SimulateSystem(m(i), con, tF);
end

sim = SimulateSelect(m(1), con, lintimes, opts);
rand_state = rng(1);
obj = obj.AddData(sim.sol);
rng(rand_state);

clear n

%% Define priors
% Topologies prior is uniform
opts.PriorTopology = zeros(nTop,1) + 1/nTop;

% Parameter prior is a spherical Gaussian
CI = 100; % Fold change to 95% CI
variance = (log(CI)/1.96).^2; % Log space variance

mux = cell(nTop,1);
Vx = cell(nTop,1);
for i = 1:nTop
    % Mean
    mux{i} = zeros(m(i).nk,1) + 1e-1;
    
    % Variance (normalized with no covariance)
    Vx{i} = diag(mux{i} .* (zeros(m(i).nk,1) + variance) .* mux{i});
end

%% Create prior objectives
objPrior = Gzero(nTop);
for i = 1:nTop
    objPrior(i) = objectiveLogNormalPriorOnKineticParameters(mux{i}, Vx{i});
end

%% Fit all models
for i = 1:nTop
    mstart = Update(m(i), zeros(m(i).nk,1) + 1e-1);
    optsFit = opts;
    optsFit.UseParams = 1:m(i).nk;
    m(i) = FitObjective(mstart, con, [obj; objPrior(i)], optsFit);
end

clear i optsFit

%% Topological probabilities
optsTop = opts;
optsTop.UseParams = {1:m(1).nk; 1:m(2).nk; 1:m(3).nk; 1:m(4).nk};
pmy = TopologyProbability(m, con, obj, objPrior, [], [], optsTop);

%% Optimal experimental design for topology uncertainty
target = @entropy;
[best, data] = BestTopologyExperiment(m, con, obj, objPrior, [], [], [con;con], [obj,obj], target, optsTop);
