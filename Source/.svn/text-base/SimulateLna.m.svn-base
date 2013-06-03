function [varargout] = SimulateLna(m, con, opts)
%SimulateLna Integrate the concentration of every species and its variance
%   according to the lienar noise approximation
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
%   opts: [ options struct scalar ]
%       Optional
%       .V0 [ matrix nx by nx {zeros(nx,nx)} ]
%           The initial variance matrix of the state species
%       .UseModelICs [ logical scalar {false} ]
%           Indicates that the model's initial conditions should be used
%           instead of those of the experimental conditions
%       .UseModelInputs [ logical scalar {false} ]
%           Indicates that the model's inputs should be used instead of
%           those of the experimental conditions
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
%   Simulate(m, con, opts)
%   	Plots the concentrations under each condition
%
%   sim = Simulate(m, con, opts)
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
%       
%   The following are permitted if nCon = 1
%   [t, y, Vy] = Simulate(m, con, opts)
%       sim.t, sim.y, sim.Vy
%
%   [t, y, Vy, x, Vx] = Simulate(m, con, opts)
%     	sim.t, sim.y, sim.Vy, sim.x, sim.Vx
%
%   [t, y, Vy, x, Vx, sol] = Simulate(m, con, opts)
%     	sim.t, sim.y, sim.Vy, sim.x, sim.Vx, sim.sol
%       
%	Special
%   This function can also be requested to return an empty array of
%   structures with the same fields as simulation. This may be necessary to
%   initialize an array that the user will later fill.
%
% 	empty = Simulate(m)
%   	For integer m, create an empty simulation array m by 1
%
%   empty = Simulate([m,n,p...])
%       For integers array, create an empty simulation array size
%       [m,n,p...]

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 3
    opts = [];
    if nargin < 2
        con = [];
        if nargin < 1
            m = [];
        end
    end
end

% Special case: return empty structure array if inputs are numeric
if isnumeric(m)
    varargout{1} = emptystruct(m, 'Type', 'Name', 't', 'y', 'Vy', 'x', 'Vx', 'sol');
    return
end

assert(nargin >= 2, 'KroneckerBio:SimulateLna:TooFewInputs', 'SimulateLna requires at least 2 input arguments')
assert(nargout <= 4, 'KroneckerBio:SimulateLna:FourOrFewerOutputs', 'SimulateLna must have between 0 and 4 outputs')
assert(isscalar(m), 'KroneckerBio:SimulateLna:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;

defaultOpts.V0             = zeros(m.nx,m.nx);

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nCon = numel(con);

% Refresh conditions
con = refreshCon(m, con);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTolLna(opts.AbsTol, 1, false(nCon,1), nx, nCon);

%% Run integration for each experiment
sim = emptystruct(nCon, 'Type', 'Name', 't', 'y', 'Vy', 'x', 'Vx', 'sol');
intOpts = opts;

for iCon = 1:nCon
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};

    % Integrate system
    if verbose; fprintf(['Integrating linear noise approximation for ' con(iCon).Name '...']); end
    sol = integrateLna(m, con(iCon), intOpts);
    if verbose; fprintf('done.\n'); end
    
    % Store results
    sim(iCon).Type  = 'Simulation.LinearNoiseApproximation';
    sim(iCon).Name  = [m.Name ' in ' con.Name];
    sim(iCon).t     = sol.x;
    sim(iCon).y     = @(t, varargin)evaluateOutputs(sol, t, varargin{:});
    sim(iCon).Vy    = @(t, varargin)evaluateOutputVariance(sol, t, varargin{:});
    sim(iCon).x     = @(t, varargin)evaluateStates(sol, t, varargin{:});
    sim(iCon).Vx    = @(t, varargin)evaluateStateVariance(sol, t, varargin{:});
    sim(iCon).sol   = sol;
end

%% Work-down
switch (nargout)
    case 0
        % Save the hold state of the figure
        holdState = ishold;
        % Draw each result
        for iCon = 1:nCon
            % TODO: write a function for plotting LNA
            plotExperiment(m, sim(iCon));
            hold on;
        end
        % Reset the hold state
        if ~holdState
            hold off;
        end
    case 1
        varargout{1} = sim;
    case 2
        varargout{1} = sim.t;
        varargout{2} = sim.y;
    case 3
        varargout{1} = sim.t;
        varargout{2} = sim.y;
        varargout{3} = sim.Vy;
    case 4
        varargout{1} = sim.t;
        varargout{2} = sim.y;
        varargout{3} = sim.Vy;
        varargout{4} = sim.x;
    case 5
        varargout{1} = sim.t;
        varargout{2} = sim.y;
        varargout{3} = sim.Vy;
        varargout{4} = sim.x;
        varargout{5} = sim.Vx;
    case 6
        varargout{1} = sim.t;
        varargout{2} = sim.y;
        varargout{3} = sim.Vy;
        varargout{4} = sim.x;
        varargout{5} = sim.Vx;
        varargout{5} = sim.sol;
    otherwise
        error(nargoutchk(0, 6, nargout, 'struct'));
end

end
% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Evaluation functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = evaluateOutputs(sol, t, ind)
[ny, nx] = size(sol.C1);

if nargin < 3
    val = sol.C1 * deval(sol, t,1:nx) + sol.C2 * sol.u(t) + repmat(sol.c, 1,numel(t));
else
    assert(all(ind <= ny), 'KroneckerBio:SimulateLna:InvalidOutputIndex', 'The requested output is larger than the total number of outputs')
    val = sol.C1(ind,:) * deval(sol, t,1:nx) + sol.C2(ind,:) * sol.u(t) + repmat(sol.c(ind,:), 1,numel(t));
end
end

function val = evaluateOutputVariance(sol, t, ind)
error('The formula for the output variance needs to be implemented')
nx = size(sol.C1, 2);
ind1 = upperInd(nx);
ind2 = lowerInd(nx);

val = zeros(nx*nx, numel(t));
if nargin < 3
    extract = deval(sol, t,nx+1:nx+nVx);
    val(ind1,:) = extract;
    val(ind2,:) = extract;
    val = reshape(val, nx,nx,numel(t));
    % TODO: conversion formula
else
    assert(all(ind <= ny), 'KroneckerBio:SimulateLna:InvalidOutputIndex', 'The requested output is larger than the total number of outputs')
    extract = deval(sol, t,nx+1:nx+nVx);
    val(ind1,:) = extract;
    val(ind2,:) = extract;
    val = reshape(val, nx,nx,numel(t));
    % TODO: conversion formula
    val = val(ind,ind);
end
end

function val = evaluateStates(sol, t, ind)
nx = size(sol.C1, 2);

if nargin < 3
    val = deval(sol, t,1:nx);
else
    assert(all(ind <= nx), 'KroneckerBio:SimulateLna:InvalidStateIndex', 'The requested state is larger than the total number of states')
    val = deval(sol, t, ind);
end
end

function val = evaluateStateVariance(sol, t, ind)
nx = size(sol.C1, 2);
ind1 = upperInd(nx);
ind2 = lowerInd(nx);

val = zeros(nx*nx, numel(t));
if nargin < 3
    extract = deval(sol, t,nx+1:nx+nVx);
    val(ind1,:) = extract;
    val(ind2,:) = extract;
    val = reshape(val, nx,nx,numel(t));
else
    assert(all(ind <= nx), 'KroneckerBio:SimulateLna:InvalidStateIndex', 'The requested state is larger than the total number of states')
    extract = deval(sol, t,nx+1:nx+nVx);
    val(ind1,:) = extract;
    val(ind2,:) = extract;
    val = reshape(val, nx,nx,numel(t));
    val = val(ind,ind);
end
end
