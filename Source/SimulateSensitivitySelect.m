function [varargout] = SimulateSensitivitySelect(m, con, tGet, opts)
%SimulateSensitivitySelect Integrate the sensitivities of every species with
%   respect to every parameter over all time returns the values at select
%   time points
%
%   Mathematically: dx/dT = Integral(df/dx * dx/dT + df/dT, t=0:tF)
%   
%   sim = SimulateSensitivitySelect(m, con, tGet, opts)
%   
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   tGet: [ nonegative vector ]
%       Indicates which time points will be returned. This does not need
%       be sorted. Times larger than con.tF will return NaN for all values.
%   opts: [ options struct scalar {} ]
%       .UseModelSeeds [ logical scalar {false} ]
%           Indicates that the model's seed parameters should be used
%           instead of those of the experimental conditions
%       .UseModelInputs [ logical scalar {false} ]
%           Indicates that the model's inputs should be used instead of
%           those of the experimental conditions
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Indicates the kinetic parameters whose sensitivities are
%           desired
%       .UseSeeds [ logical matrix ns by nCon | logical vector ns |
%                   positive integer vector {[]} ]
%           Indicates the seed parameters whose sensitivities are desired
%       .UseControls [ cell vector nCon of logical vectors or positive 
%                      integer vectors | logical vector nq | positive 
%                      integer vector {[]} ]
%           Indicates the input control parameters whose sensitivites are
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
%   SimulateSensitivitySelect(m, con, opts)
%   	Plots the sensitivities under each condition
%
%   sim = SimulateSensitivitySelect(m, con, opts)
%   	A vector of structures with each entry being the simulation
%       under one of the conditions.
%       .t tGet
%       .y [ matrix ny by numel(tGet) ]
%           The value of the outputs at each selected time point
%       .x [ matrix nx by numel(tGet) ]
%           The value of the states at each selected time point
%       .dydT [ matrix ny*nT by numel(tGet) ]
%           The value of the sensitivites of the outputs at each selected
%           time point
%       .dxdT [ matrix nx by numel(tGet) ]
%           The value of the sensitivities of the states at each selected
%           time point
%       .sol [ struct scalar ]
%           The discrete integrator solution to the system

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 4
    opts = [];
end

assert(nargin >= 3, 'KroneckerBio:SimulateSensitivitySelect:TooFewInputs', 'SimulateSensitivitySelect requires at least 3 input arguments')
assert(isscalar(m), 'KroneckerBio:SimulateSensitivitySelect:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelSeeds  = false;
defaultOpts.UseModelInputs = false;

defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseSeeds       = [];
defaultOpts.UseControls    = [];

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
ny = m.ny;
nk = m.nk;
nCon = numel(con);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTx] = fixUseSeeds(opts.UseSeeds, opts.UseModelSeeds, nx, nCon);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseControls nTq] = fixUseControls(opts.UseControls, opts.UseModelInputs, nCon, m.nq, cat(1,con.nq));

nT = nTk + nTx + nTq;

% Refresh conditions
con = refreshCon(m, con);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, false(nCon,1), nx, nCon, false, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams, opts.UseSeeds, opts.UseControls);

%% Run integration for each experiment
sim = emptystruct(nCon, 'Type', 'Name', 't', 'y', 'x', 'dydT', 'dxdT', 'sol');

for iCon = 1:nCon
    % Modify opts structure
    intOpts = opts;
    intOpts.AbsTol = opts.AbsTol{iCon};
    
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
    
    % Integrate [x; dx/dT] over time
    if verbose; fprintf(['Integrating sensitivities for ' con(iCon).Name '...']); end
    sol = integrateSensSelect(m, con(iCon), tGet, intOpts);
    if verbose; fprintf('done.\n'); end
    
    % Store results
    sim(iCon).Type  = 'Simulation.Sensitivity.SelectPoints';
    sim(iCon).Name  = [m.Name ' in ' con(iCon).Name];
    sim(iCon).t     = sol.x;
    sim(iCon).y     = bsxfun(@plus, sol.C1*sol.y(1:nx,:) + sol.C2*sol.u, sol.c);
    sim(iCon).x     = sol.y(1:nx,:);
    sim(iCon).dydT  = reshape(sol.C1*reshape(sol.y(nx+1:end,:), nx,inT*numel(sol.x)), ny*inT,numel(sol.x));
    sim(iCon).dxdT  = sol.y(nx+1:end,:);
    sim(iCon).sol   = sol;
end

%% Work-down
if nargout == 0
    % Draw each result
    for iCon = 1:nCon
        subplot(nCon,1,iCon)
        plotSensitivityExperiment(m, sim(iCon), 'o-', 'Linewidth', 2);
    end
else
    varargout{1} = sim;
end
