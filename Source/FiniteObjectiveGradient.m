function D = FiniteObjectiveGradient(m, con, obj, opts)
%FiniteObjectiveGradient Evaluate the gradient of a set of objective
%   functions using the finite difference approximation
%
%   D = FiniteObjectiveGradient(m, con, obj, opts)
%
%   This function is identical to ObjectiveGradient except that the
%   gradient is calculated differently. This is useful for testing that
%   custom-made objective functions have been coded correctly and that
%   integration tolerances have been set appropriately.
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix ]
%       The objective structures defining the objective functions to be
%       evaluated.
%       .UseModelSeeds [ logical scalar {false} ]
%           Indicates that the model's seed parameters should be used
%           instead of those of the experimental conditions
%       .UseModelInputs [ logical scalar {false} ]
%           Indicates that the model's inputs should be used instead of
%           those of the experimental conditions
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Which kinetic parameters the gradient will be calculated on
%       .UseSeeds [ logical matrix nx by nCon | logical vector nx |
%                   positive integer vector {[]} ]
%           Which seed parameters the gradient will be calculated on
%       .UseControls [ cell vector nCon of logical vectors or positive 
%                      integer vectors | logical vector nq | positive 
%                      integer vector {[]} ]
%           Which input control parameters the gradient will be calculated
%           on
%     	.ObjWeights [ real matrix nObj by nCon {ones(nObj,nCon)} ]
%           Applies a post evaluation weight on each objective function
%           in terms of how much it will contribute to the final objective
%           function value
%       .Normalized [ logical scalar {true} ]
%           Indicates if the gradient should be computed in log parameters
%           space
%       .RelTol [ nonnegative scalar {1e-6} ]
%           Relative tolerance of the integration
%       .AbsTol [ cell vector of nonnegative vectors | nonnegative vector |
%                 nonegative scalar {1e-9} ]
%           Absolute tolerance of the integration. If a cell vector is
%           provided, a different AbsTol will be used for each experiment.
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%
%   Outputs
%       D: [ real vector nT ]
%           The sum of all objective function gradients

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inVuts
assert(nargin >= 3, 'KroneckerBio:FiniteObjectiveGradient:AtLeastThreeInputs', 'FiniteObjectiveGradient requires at least 3 input arguments.')
if nargin < 4
    opts = [];
end

assert(isscalar(m), 'KroneckerBio:FiniteObjectiveGradient:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelSeeds  = false;
defaultOpts.UseModelInputs = false;

defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseSeeds       = [];
defaultOpts.UseControls    = [];

defaultOpts.ObjWeights     = ones(size(obj));

defaultOpts.Normalized     = true;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
ns = m.ns;
nk = m.nk;
nCon = numel(con);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseICs is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, opts.UseModelSeeds, ns, nCon);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseControls, nTq] = fixUseControls(opts.UseControls, opts.UseModelInputs, nCon, m.nq, cat(1,con.nq));

nT = nTk + nTs + nTq;

% Store starting parameter sets
T0 = collectActiveParameters(m, con, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams, opts.UseSeeds, opts.UseControls);

% Refresh conditions and objectives
con = refreshCon(m, con);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 1, opts.continuous, nx, nCon);

%% Loop through conditions
D = zeros(nT,1);

% Initial value
if verbose; fprintf('Initial round\n'); end
G = computeObj(m, con, obj, opts);

for iT = 1:nT
    if verbose; fprintf('Step %d of %d\n', iT, nT); end
    
    % Set baseline parameters
    Ti = T0(iT);
    Tup = T0;
    Tdown = T0;
    
    % Change current parameter by finite amount
    if opts.Normalized
        diff = Ti * 1e-8;
    else
        diff = 1e-8;
    end
    
    % Compute objective values
    Tup(iT) = Tup(iT) + diff;
    [m, con] = updateAll(m, con, Tup, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams, opts.UseSeeds, opts.UseControls);
    Gup = computeObj(m, con, obj, opts);

    Tdown(iT) = Tdown(iT) - diff;
    [m, con] = updateAll(m, con, Tdown, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams, opts.UseSeeds, opts.UseControls);
    Gdown = computeObj(m, con, obj, opts);

    % Compute D
    if opts.Normalized
        D(iT) = Ti * ( (Gup - G) / diff + (G - Gdown) / diff ) / 2;
    else
        D(iT) = ( (Gup - G) / diff + (G - Gdown) / diff ) / 2;
    end
end
