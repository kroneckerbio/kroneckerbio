clear; close all; clc

%% Options
% Global
opts.Verbose = 1;
opts.Order = 2;
opts.UseParams = 1:48;
opts.UseSeeds = [];
opts.UseModelSeeds = false;
nTk = numel(opts.UseParams);
nTs = numel(opts.UseSeeds);
nV = nTk + nTs;
opts.RelTol = 1e-6;
opts.AbsTol = 1e-9;
opts.NewJacobianMethod = true;
opts.UseMEX = false;
opts.MEXDirectory = 'mexfuns';

% Generate data
opts.useSeed = true;
opts.Seed = 1337;
opts.populate = true;
opts.NonNegMeasurements = true;

% Fitting
opts.TolOptim = 1e-6;
opts.Normalized = true;
opts.UseAdjoint = false;
opts.AdaptAbsTol = true;
opts.MaxIter = 300;
opts.MaxFunEvals = 5000;
opts.Restart = 0;
opts.RestartNoise = 0.001;
opts.GlobalOptimization = false;

% Loading SBML
opts.UseNames = true;

%% Load model
yNames = {
    'EGF'
    'NGF'
    'freeEGFReceptor'
    'boundEGFReceptor'
    'freeNGFReceptor'
    'boundNGFReceptor'
    'SosInactive'
    'SosActive'
    'P90RskInactive'
    'P90RskActive'
    'RasInactive'
    'RasActive'
    'Raf1Inactive'
    'Raf1Active'
    'BRafInactive'
    'BRafActive'
    'MekInactive'
    'MekActive'
    'ErkInactive'
    'ErkActive'
    'PI3KInactive'
    'PI3KActive'
    'AktInactive'
    'AktActive'
    'C3GInactive'
    'C3GActive'
    'Rap1Inactive'
    'Rap1Active'
    'RasGapActive'
    'RapGapActive'
    'PP2AActive'
    'Raf1PPtase'
    };

yNames2 = {'Out1'; 'Out2'};
% yExprs2 = {'EGF + 2*NGF', '(RasGapActive + kSos*RapGapActive)/2 + sqrt(AktActive)^(kRap1ToBRaf)'};
yExprs2 = {'("RasGapActive" + kSos*RapGapActive)/2 + sqrt(AktActive)^(kRap1ToBRaf)', 'EGF + 2*NGF'};

symModel = sbml2Symbolic('Brown_EGFNGF.xml', opts);

% symModel = AddOutputsToSymbolic(symModel, yNames, [], opts);
symModel = AddOutputsToSymbolic(symModel, yNames2, yExprs2, opts);

% m = LoadModelSbmlAnalytic('Brown_EGFNGF.xml', yNames, [], [], opts);
m = symbolic2PseudoKronecker(symModel, opts);

if opts.UseMEX
    compileMEXFunctions(opts.MEXDirectory);
end

%% Create nominal experiment

% Set the final simulation time
tF = 120;

% Set up the input functions
u0 = m.u; % Default input values. The input is a constant value for all time.
nq = 0; % Number of input parameters
nu = m.nu; % Number of inputs
ufun = @(t,q) repmat(u0,1,numel(t)); % Sets up a function that returns a matrix consisting of u0 repeated numel(t) times
dudqfun = @(t,q) zeros(nu, nq); % The derivatives of u wrt q are zero, since the values are constant. dudq is then a nu-by-nq matrix of zeros
d2udq2fun = @(t,q) zeros(nu*nq,nq); % The second derivative of u wrt q is also zero.
discontinuities = [];
q = zeros(nq,1);
inp = Input(m,ufun,discontinuities,q,dudqfun,d2udq2fun); % Creates a struct containing the input functions

% Change the initial concentrations of EGF and NGF
s = m.s;
s(1) = 1000; % EGF
s(2) = 4560; % NGF

% Create the experiment
d = []; % Dosing
con = experimentInitialValue(m, s, inp, d, 'nominal');

%% Simulate nominal experiment and plot states

% Simulate the model under the experiment
tic; sim = SimulateSystem(m,con,tF); toc

%%%%% Plot all the states %%%%%%
t = 0:0.1:tF; % Choose which times to plot
xall = sim.x(t); % Get all the state concentrations for the desired times. x is an nx-by-nt matrix of concentrations, with each column listing the concentrations of all the states in the model at a particular time.
figure; plot(t,xall)
title('All state values')

%%%%% Plot only select states %%%%%
% Create a list of the state names
listofstates = {m.States.Name};

% Set states to plot
statestoplot = {'ErkActive'}; 

% Get the indices associated with each state
nstatestoplot = length(statestoplot);
stateindices = zeros(nstatestoplot,1);
for xi = 1:nstatestoplot
    isthisstate = strcmp(listofstates,statestoplot{xi});
    stateindices(xi) = find(isthisstate);
end
xselect = sim.x(t,stateindices);
figure; plot(t,xselect);
title('Select state values')
legend(statestoplot,'Location','Best')
    
%% Generate data to fit model to

% Choose to generate data for these two outputs in the model
outputnames = {'Out1';'Out2'};
noutputstofit = length(outputnames);
outputlist = {m.Outputs.Name};
outputs = zeros(noutputstofit,1);
for yi = 1:noutputstofit
    outputs(yi) = find(strcmp(outputlist,outputnames{yi}));
end

% Generate 8 time points
ntimesearly = 6;
ntimeslate = 2;
ntimes = ntimesearly + ntimeslate;
tlate = linspace(20,tF,ntimeslate+1).';
tlate = tlate(2:end);
lintimes = [linspace(20/ntimesearly,20,ntimesearly).'; tlate]; % 6 time points early, where interesting dynamics exist, and 2 time points later

% Create the function that returns the standard error of the measurements
sd = @samplesd; % Standard deviation of measurement error is 5% floored at a value of 1

% Simulate the model to generate data to fit to
nonNegMeasurements = true; % Concentrations being fit will be non-negative
populate = true; % I don't have the data points yet, so I want the function to perform the simulations for me
seed = []; % Let the function generate the RNG seed for me
[obj, dataFit, obs] = generateFakeData(m, con, outputs, lintimes, sd, nonNegMeasurements, populate, seed, opts);

%% Scramble parameters and set lower and upper bounds

% Starting (scrambled) parameters are the geometric mean of each class of
% parameters. Parameter classes include binding and unbinding parameters,
% kcat values, and Km values.
% Bounds are 0.001X below the smallest and 1000X the largest value in each 
% class
opts.LowerBound = zeros(nV,1) + min(m.k(opts.UseParams));
opts.UpperBound = zeros(nV,1) + max(m.k(opts.UseParams));

% Set up indices for classes of parameters
konIndex = [1;3];
koffIndex = [2;4];
kcatIndex  = 5:2:47;
KmIndex    = 6:2:48;

roundtonearestoom = @(x) 10.^round(log10(x));
kstart = m.k;
kstart = roundtonearestoom(kstart);

%kstart(konindex) = geomean(kstart(konindex));
opts.LowerBound(konIndex) = zeros(size(konIndex)) + min(m.k(konIndex)) / 1000;
opts.UpperBound(konIndex) = zeros(size(konIndex)) + max(m.k(konIndex)) * 1000;

%kstart(koffindex) = geomean(kstart(koffindex));
opts.LowerBound(koffIndex) = zeros(size(koffIndex)) + min(m.k(koffIndex)) / 1000;
opts.UpperBound(koffIndex) = zeros(size(koffIndex)) + max(m.k(koffIndex)) * 1000;

%kstart(kcatindex) = geomean(kstart(kcatindex));
opts.LowerBound(kcatIndex) = zeros(size(kcatIndex)) + min(m.k(kcatIndex)) / 1000;
opts.UpperBound(kcatIndex) = zeros(size(kcatIndex)) + max(m.k(kcatIndex)) * 1000;

%kstart(kmindex) = geomean(kstart(kmindex));
opts.LowerBound(KmIndex) = zeros(size(KmIndex)) + min(m.k(KmIndex)) / 1000;
opts.UpperBound(KmIndex) = zeros(size(KmIndex)) + max(m.k(KmIndex)) * 1000;

opts.MaxIter = 100;

mstart = m.Update(kstart);

clear konindex koffindex kcatindex kmindex pstart

%% Fit to the objective function to recover the original parameter values

[mfit, confit, G, D] = FitObjective(mstart,con,obj,opts);

%% Plot the resulting fit

% Get data being fit to
ydata = reshape(dataFit.Measurements, noutputstofit, ntimes);
t = dataFit.Times';

% Get fitted model data
simfit = SimulateSystem(mfit, confit, obs);
ymodel = reshape(simfit.true_measurements, noutputstofit, ntimes);

% Get standard error values
sigma = zeros(size(ymodel));
for yi = 1:numel(ydata)
    [yind,tind] = ind2sub([noutputstofit,ntimes],yi);
    yval = ymodel(yi);
    sigma(yi) = sd(t(tind),yind,yval);
end

% Plot the quality of fit
for yi = 1:noutputstofit
    figure; hold all
    plot(t,ymodel(yi,:));
    errorbar(t,ydata(yi,:),sigma(yi,:),'ro');
    hold off
end
