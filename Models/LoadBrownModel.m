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
opts.Restart = 4;
opts.RestartNoise = 0.001;

%% Load model
m = LoadModelSbmlAnalytic('Brown_EGFNGF.xml', [], 1:32, [], opts);

ic = m.ic;
ic(1) = 1000; %EGF
ic(2) = 4560; %NGF
m = m.update(m.p, ic);

clear ic

%% Create nominal experiment
tF = 120;
u = [];
con = constructExperiment(m, tF, u, [], [], [], 'nominal');

%% Generate data
outputs = (1:m.nY).';
ntimes = 100;
lintimes = linspace(tF/ntimes, tF, ntimes).'; % 100 evenly spaced time points
sd = @(t,yInd,yVal)(max(0.1*yVal,1)); % Standard deviation of measurement error is 10% floored at 1
[obj dataFit] = generateFakeData(m, con, outputs, lintimes, sd, opts);

clear ntimes

%% Prepare for fitting
% Starting (scrambled) parameters are the geometric mean of each class
% Bounds are 0.001X below the smallest and 1000X the largest value in each 
% class
opts.LowerBound = zeros(nV,1) + min(m.p(opts.UseParams));
opts.UpperBound = zeros(nV,1) + max(m.p(opts.UseParams));

kbindindex = 1:4;
kcatindex  = 5:2:47;
kmindex    = 6:2:48;

pstart = m.p;
pstart(kbindindex) = geomean(pstart(kbindindex));
opts.LowerBound(kbindindex) = zeros(size(kbindindex)) + min(m.p(kbindindex)) / 1000;
opts.UpperBound(kbindindex) = zeros(size(kbindindex)) + max(m.p(kbindindex)) * 1000;

pstart(kcatindex) = geomean(pstart(kcatindex));
opts.LowerBound(kcatindex) = zeros(size(kcatindex)) + min(m.p(kcatindex)) / 1000;
opts.UpperBound(kcatindex) = zeros(size(kcatindex)) + max(m.p(kcatindex)) * 1000;

pstart(kmindex) = geomean(pstart(kmindex));
opts.LowerBound(kmindex) = zeros(size(kmindex)) + min(m.p(kmindex)) / 1000;
opts.UpperBound(kmindex) = zeros(size(kmindex)) + max(m.p(kmindex)) * 1000;

mstart = m.update(pstart, m.ic);

clear kbindindex kcatindex kmindex pstart
