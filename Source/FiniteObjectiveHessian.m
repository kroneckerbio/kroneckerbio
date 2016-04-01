function H = FiniteObjectiveHessian(m, con, obj, opts)
%FiniteObjectiveHessian Evaluate the hessian of a set of objective
%   functions using the finite difference approximation
%
%   H = FiniteObjectiveHessian(m, con, obj, opts)
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix n_obj by n_con ]
%       The objective structures defining the objective functions to be
%       evaluated. Note that this matrix must have a number of columns
%       equal to numel(con) (e.g. one objective for each experimental
%       condition is a row vector and multiple objective structures for
%       a single experimental conditions is a column vector).
%   opts: [ options struct scalar {} ]
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
%    	.UseAdjoint [ logical scalar {false} ]
%           Indicates whether the gradient should be calculated via the
%           adjoint method or the forward method.
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
%       H: [ real vector nT ]
%           The sum of all objective function hessians

% (c) 2015 David R Hagen and Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 4
    opts = [];
end

assert(nargin >= 3, 'KroneckerBio:ObjectiveGradient:TooFewInputs', 'ObjectiveGradient requires at least 3 input arguments')
assert(isscalar(m), 'KroneckerBio:ObjectiveGradient:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = [];
defaultOpts.AbsTol         = [];

defaultOpts.ComplexStep  = false;

defaultOpts.UseParams        = 1:m.nk;
defaultOpts.UseSeeds         = [];
defaultOpts.UseInputControls = [];
defaultOpts.UseDoseControls  = [];

defaultOpts.ObjWeights     = ones(size(obj));

defaultOpts.Normalized     = true;
defaultOpts.UseAdjoint     = false;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
ns = m.ns;
nk = m.nk;

% Ensure structures are proper sizes
[con, n_con] = fixCondition(con);
[obj, n_obj] = fixObjective(obj, n_con);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, ns, n_con);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, n_con, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, n_con, cat(1,con.nh));

nT = nTk + nTs + nTq + nTh;

% Store starting parameter sets
T0 = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nx, n_con, opts.UseAdjoint, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

%% Loop through conditions
H = zeros(nT,nT);

% Initial value
if verbose; fprintf('Initial round\n'); end
[~, D] = computeObjGrad(m, con, obj, opts);

for iT = 1:nT
    if verbose; fprintf('Step %d of %d\n', iT, nT); end
    
    % Set baseline parameters
    T_i = T0(iT);
    T_up = T0;
    
    % Change current parameter by finite amount
    step_size = 1e-8;
    if opts.Normalized
        norm_factor = T_i;
    else
        norm_factor = 1;
    end
    if opts.ComplexStep
        imag_factor = 1i;
    else
        imag_factor = 1;
    end
    diff = step_size * norm_factor * imag_factor;
    
    % Compute objective values
    T_up(iT) = T_up(iT) + diff;
    [m, con] = updateAll(m, con, T_up, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
    [~, D_up] = computeObjGrad(m, con, obj, opts);

    % Compute D
    if opts.ComplexStep
        H(:,iT) = imag(D_up) ./ step_size;
    else
        H(:,iT) = (D_up - D) ./ step_size;
    end
end
