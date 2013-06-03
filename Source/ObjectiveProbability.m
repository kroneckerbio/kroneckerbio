function p = ObjectiveProbability(m, con, obj, opts)
%ObjectiveProbability Evaluate the likelihood of a set of 
%   information-theory-based objective functions
%
%   p = ObjectiveProbability(m, con, obj, opts)
%   
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix ]
%       The information-theory-based objective structures defining the
%       probability distribution to be evaluated.
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
assert(nargin >= 3, 'KroneckerBio:ObjectiveProbability:TooFewInputs', 'ObjectiveProbability requires at least 3 input arguments')
if nargin < 4
    opts = [];
end

assert(isscalar(m), 'KroneckerBio:ObjectiveProbability:MoreThanOneModel', 'The model structure must be scalar')

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
defaultOpts.UseAdjoint     = false;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
ns = m.ns;
nk = m.nk;
nCon = numel(con);
nObj = size(obj,1);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, opts.UseModelSeeds, ns, nCon);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseControls, nTq] = fixUseControls(opts.UseControls, opts.UseModelInputs, nCon, m.nq, cat(1,con.nq));

nT = nTk + nTs + nTq;

% Refresh conditions and objectives
con = refreshCon(m, con);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 1, opts.continuous, nx, nCon);

%% Probability
% Initialize variables
p = 1;

for iCon = 1:nCon

    intOpts = opts;
    
    % If opts.UseModelSeeds is false, the number of variables can change
    if opts.UseModelSeeds
        UseSeeds_i = opts.UseSeeds;
    else
        UseSeeds_i = opts.UseSeeds(:,iCon);
    end
    intOpts.UseSeeds = UseSeeds_i;
    inTs = nnz(UseSeeds_i);
    
    % If opts.UseModelInputs is false, the number of variables can change
    if opts.UseModelInputs
        UseControls_i = opts.UseControls{1};
    else
        UseControls_i = opts.UseControls{iCon};
    end
    intOpts.UseControls = UseControls_i;
    inTq = nnz(UseControls_i);
    
    inT = nTk + inTs + inTq;

    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    intOpts.ObjWeights = opts.ObjWeights(:,iCon);
    tGet = opts.tGet{iCon};
        
    % Integrate
    if opts.continuous(iCon) && opts.complex(iCon)
        sol = integrateObj(m, con(iCon), obj(:,iCon), intOpts);
    elseif opts.complex(iCon)
        sol = integrateSys(m, con(iCon), intOpts);
    elseif opts.continuous(iCon)
        sol = integrateObjSelect(m, con(iCon), obj(:,iCon), tGet, intOpts);
    else
        sol = integrateSysSelect(m, con(iCon), tGet, intOpts);
    end
    
    % Add fields for prior objectives
    sol.UseParams = opts.UseParams;
    sol.UseSeeds = UseSeeds_i;
    sol.UseControls = UseControls_i;
    
    % Probability
    for iObj = 1:nObj
        p = p * obj(iObj,iCon).p(sol);
    end
end
