function sim = SimulateLna(m, con, obs, opts)
%SimulateLna Integrate the concentration of every species and its variance
%   according to the linear noise approximation
%
%   Mathematically: x = Integral([f; vec(Vdot)], t=0:tF)
%   
%   [...] = Simulate(m, con, opts)
%   
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   opts: [ options struct scalar {} ]
%       .V0 [ matrix nx by nx {zeros(nx,nx)} ]
%           The initial variance matrix of the state species
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
%   SimulateLna(m, con, opts)
%   	Plots the concentrations under each condition
%
%   sim = SimulateLna(m, con, opts)
%   	A vector of structures with each entry being the simulation
%       under one of the conditions.
%       .t [ sorted nonnegative row vector ]
%           Timepoints chosen by the ode solver
%       .y [ handle @(t,y) returns matrix numel(y) by numel(t) ]
%           This function handle evaluates some outputs y of the system at
%           some particular time points t. The user may exclude y, in which
%           case all outputs are returned.
%       .Vy [ handle @(t,y) returns matrix numel(y) by numel(t) ]
%           This function handle evaluates the variance of some outputs y
%           at some particular time points t. The user may exclude y, in
%           which case the variance of all states is returned.
%       .x [ handle @(t,x) returns matrix numel(x) by numel(t) ]
%           This function handle evaluates some states x of the system at
%           some particular time points t. The user may exclude x, in which
%           case all states are returned.
%       .Vx [ handle @(t,x) returns matrix numel(x) by numel(t) ]
%           This function handle evaluates the variance of some states x at
%           some particular time points t. The user may exclude x, in which
%           case the variance of all states is returned.
%       .sol [ odesolver struct scalar ]
%           The integrator solution to the system

% (c) 2015 David R Hagen
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 4
    opts = [];
end

assert(nargin >= 3, 'KroneckerBio:SimulateLna:TooFewInputs', 'SimulateLna requires at least 3 input arguments')
assert(isscalar(m), 'KroneckerBio:SimulateLna:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose = 1;

defaultOpts.RelTol  = [];
defaultOpts.AbsTol  = [];

defaultOpts.V0      = zeros(m.nx,m.nx);

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;

% Ensure structures are proper sizes
[con, n_con] = fixCondition(con);
[obs, n_obs] = fixObservation(obs, n_con);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTolLna(opts.AbsTol, 1, false(n_con,1), nx, n_con);

%% Run integration for each experiment
sim = emptystruct([n_obs,n_con]);

for i_con = 1:n_con
    opts_i = opts;
    opts_i.AbsTol = opts.AbsTol{i_con};
    
    if verbose; fprintf(['Integrating linear noise approximation for ' con(i_con).Name '...']); end
    ints = integrateAllLna(m, con(i_con), obs(:,i_con), opts_i);
    if verbose; fprintf('done.\n'); end
    
    for i_obs = 1:n_obs
        sim = insertstruct(sim, obs(i_obs,i_con).Lna(ints(i_obs)), i_obs,i_con);
    end
end
