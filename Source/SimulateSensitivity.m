function sim = SimulateSensitivity(m, con, obs, opts)
%SimulateSensitivity Integrate the sensitivities of every species with
%   respect to every parameter over all time
%
%   Mathematically: dx/dT = Integral(df/dx * dx/dT + df/dT, t=0:tF)
%   
%   sim = SimulateSensitivity(m, con, opts)
%   
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   opts: [ options struct scalar {} ]
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Indicates the kinetic parameters whose sensitivities are
%           desired
%       .UseSeeds [ logical matrix ns by nCon | logical vector ns |
%                   positive integer vector {[]} ]
%           Indicates the seed parameters whose sensitivities are desired
%       .UseInputControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the input control parameters whose sensitivites are
%           desired
%       .UseDoseControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the dose control parameters whose sensitivites are
%           desired
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
%   SimulateSensitivity(m, con, opts)
%   	Plots the sensitivities under each condition
%
%   sim = SimulateSensitivity(m, con, opts)
%   	A vector of structures with each entry being the simulation
%       under one of the conditions.
%       .t [ sorted nonnegative row vector ]
%           Timepoints chosen by the ode solver
%       .y [ handle @(t,y) returns matrix numel(y) by numel(t) ]
%           This function handle evaluates some outputs y of the system at
%           some particular time points t. The user may exclude y, in which
%           case all outputs are returned.
%       .x [ handle @(t,x) returns matrix numel(x) by numel(t) ]
%           This function handle evaluates some states x of the system at
%           some particular time points t. The user may exclude x, in which
%           case all states are returned.
%       .dydT [ handle @(t,y) returns matrix numel(y)*nT by numel(t) ]
%           This function handle evaluates the sensitivity of some outputs
%           y to the active parameters of the system at some particular
%           time points t. The user may exclude y, in which case all
%           outputs are returned.
%       .dxdT [ handle @(t,x) returns matrix numel(x) by numel(t) ]
%           This function handle evaluates the sensitivity of some states
%           x to the active parameters of the system at some particular
%           time points t. The user may exclude x, in which case all
%           states are returned.
%       .sol [ odesolver struct scalar ]
%           The integrator solution to the system

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 4
    opts = [];
end

assert(nargin >= 3, 'KroneckerBio:SimulateSensitivity:TooFewInputs', 'SimulateSensitivity requires at least 2 input arguments')
assert(isscalar(m), 'KroneckerBio:SimulateSensitivity:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose          = 1;

defaultOpts.RelTol           = [];
defaultOpts.AbsTol           = [];

defaultOpts.Normalized       = true;
defaultOpts.UseParams        = 1:m.nk;
defaultOpts.UseSeeds         = [];
defaultOpts.UseInputControls = [];
defaultOpts.UseDoseControls  = [];

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nk = m.nk;

% Ensure structures are proper sizes
[con, n_con] = fixCondition(con);
[obs, n_obs] = fixObservation(obs, n_con);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTx] = fixUseSeeds(opts.UseSeeds, nx, n_con);

% Ensure UseControls are cell vectors of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, n_con, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, n_con, cat(1,con.nh));

nT = nTk + nTx + nTq + nTh;

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, false(n_con,1), nx, n_con, false, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

%% Run integration for each experiment
sim = emptystruct([n_obs,n_con]);

for i_con = 1:n_con
    % Modify opts structure
    opts_i = opts;
    opts_i.AbsTol = opts.AbsTol{i_con};
    opts_i.UseSeeds = opts.UseSeeds(:,i_con);
    opts_i.UseInputControls = opts.UseInputControls{i_con};
    opts_i.UseDoseControls = opts.UseDoseControls{i_con};
    
    % Integrate [x; dx/dT] over time
    if verbose; fprintf(['Integrating sensitivities for ' con(i_con).Name '...']); end
    ints = integrateAllSens(m, con(i_con), obs(:,i_con), opts_i);
    if verbose; fprintf('done.\n'); end
    
    for i_obs = 1:n_obs
        sim = insertstruct(sim, obs(i_obs,i_con).Sensitivity(ints(i_obs)), i_obs,i_con);
    end
end
