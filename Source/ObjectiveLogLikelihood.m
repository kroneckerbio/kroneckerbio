function logp = ObjectiveLogLikelihood(m, con, obj, opts)
%ObjectiveLogLikelihood Evaluate the log likelihood of a set of 
%   information-theory-based objective functions
%
%   p = ObjectiveLogLikelihood(m, con, obj, opts)
%   
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix n_obj by n_con ]
%       The information-theory-based objective structures defining the
%       objective functions to be evaluated. Note that this matrix must
%       have a number of columns equal to numel(con) (e.g. one objective
%       for each experimental condition is a row vector and multiple
%       objective structures for a single experimental conditions is a
%       column vector).
%   opts: [ options struct scalar ]
%       Optional
%       .UseModelSeeds [ logical scalar {false} ]
%           Indicates that the model's seeds should be used instead of
%           those of the experimental conditions
%       .UseModelInputs [ logical scalar {false} ]
%           Indicates that the model's inputs should be used instead of
%           those of the experimental conditions
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Indicates the kinetic parameters that will be considered by the
%           prior objective functions
%       .UseSeeds [ logical matrix ns by nCon | logical vector ns |
%                   positive integer vector {[]} ]
%           Indicates the seeds that will be considered by the prior
%           objective function. If UseModelSeeds is true then UseSeeds can
%           be a vector of linear indexes or a vector of logicals length of
%           ns. If UseModelSeeds is false then UseSeeds can be a matrix of
%           logicals size ns by nCon. It can also be a vector of length ns,
%           and every experiment will be considered to have the same active
%           seed parameters. It can also be a vector of linear indexes into
%           the ns vector and assumed the same for all conditions.
%       .UseControls [ cell vector nCon of logical vectors or positive 
%                      integer vectors | logical vector nq | positive 
%                      integer vector {[]} ]
%           Indicates the input control parameters that will be considered
%           by the prior objective functions
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
%       p: [ real nonnegative scalar ]
%           The product of all objective function probabilities

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
assert(nargin >= 3, 'KroneckerBio:ObjectiveLogLikelihood:TooFewInputs', 'ObjectiveValue requires at least 3 input arguments')
if nargin < 4
    opts = [];
end

assert(isscalar(m), 'KroneckerBio:ObjectiveLogLikelihood:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = [];
defaultOpts.AbsTol         = [];

defaultOpts.UseParams        = 1:m.nk;
defaultOpts.UseSeeds         = [];
defaultOpts.UseInputControls = [];
defaultOpts.UseDoseControls  = [];

defaultOpts.ObjWeights     = ones(size(obj));

opts = mergestruct(defaultOpts, opts);

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

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 1, opts.continuous, nx, n_con);

%% Probability
% Initialize variables
logp = 0;

for i_con = 1:n_con
    opts_i = opts;
    opts_i.AbsTol = opts.AbsTol{i_con};
    opts_i.UseSeeds = opts.UseSeeds(:,i_con);
    opts_i.UseInputControls = opts.UseInputControls{i_con};
    opts_i.UseDoseControls = opts.UseDoseControls{i_con};
    opts_i.ObjWeights = opts.ObjWeights(:,i_con);

    % Integrate
    ints = integrateAllSys(m, con(i_con), obj(:,i_con), opts_i);
    
    % Add fields for prior objectives
    [ints.UseParams] = deal(opts.UseParams);
    [ints.UseSeeds] = deal(opts_i.UseSeeds);
    [ints.UseInputControls] = deal(opts_i.UseInputControls);
    [ints.UseDoseControls] = deal(opts_i.UseDoseControls);
    
    % Probability
    for i_obj = 1:n_obj
        logp = logp + obj(i_obj,i_con).logp(ints(i_obj));
    end
end
