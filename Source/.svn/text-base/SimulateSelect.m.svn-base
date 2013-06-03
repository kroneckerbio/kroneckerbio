function [varargout] = SimulateSelect(m, con, tGet, opts)
%SimulateSelect Integrate the concentration of every species over time
%   using mass action kinetics and return the values at select time points
%
%   Mathematically: x = Integral(f, t=0:tF)
%   
%   [...] = SimulateSelect(m, con, tGet, opts)
%   
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   tGet: [ nonegative vector ]
%       Indicates which time points will be returned. This does not need to
%       be sorted. Times larger than con.tF will return NaN for all values.
%   opts: [ options struct scalar ]
%       Optional
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
%   SimulateSelect(m, con, opts)
%   	Plots the concentrations under each condition
%
%   sim = SimulateSelect(m, con, opts)
%   	A vector of structures with each entry being the simulation
%       under one of the conditions.
%       .t tGet
%       .y [ matrix ny by numel(tGet) ]
%           The value of the outputs at each selected time point
%       .x [ matrix nx by numel(tGet) ]
%           The value of the states at each selected time point
%       .sol [ struct scalar ]
%           The discrete integrator solution to the system
%       
%   The following are permitted if nCon = 1
%   [t, y] = SimulateSelect(m, con, opts)
%       sim.t, sim.y
%
%   [t, y, x] = SimulateSelect(m, con, opts)
%     	sim.t, sim.y, sim.x
%
%   [t, y, x, sol] = SimulateSelect(m, con, opts)
%     	sim.t, sim.y, sim.x, sim.sol
%       
%	Special
%   This function can also be requested to return an empty array of
%   structures with the same fields as simulation. This may be necessary to
%   initialize an array that the user will later fill.
%
% 	empty = SimulateSelect(m)
%   	For integer m, create an empty simulation array m by 1
%
%   empty = SimulateSelect([m,n,p...])
%       For integers array, create an empty simulation array size
%       [m,n,p...]

% (c) 2010 David R Hagen, Joshua F Apgar, Jared E Toettcher, & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
if nargin < 4
    opts = [];
    if nargin < 3
        tGet = [];
        if nargin < 2
            con = [];
            if nargin < 1
                m = [];
            end
        end
    end
end

% Special case: return empty structure array if inputs are numeric
if isnumeric(m)
    varargout{1} = emptystruct(m, 'Type', 'Name', 't', 'y', 'x', 'sol');
    return
end

assert(nargin >= 3, 'KroneckerBio:SimulateSelect:TooFewInputs', 'SimulateSelect requires at least 2 input arguments')
assert(nargout <= 4, 'KroneckerBio:SimulateSelect:FourOrFewerOutputs', 'SimulateSelect must have between 0 and 4 outputs')
assert(isscalar(m), 'KroneckerBio:SimulateSelect:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;

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
opts.AbsTol = fixAbsTol(opts.AbsTol, 1, false(nCon,1), nx, nCon);

%% Run integration for each experiment
sim = emptystruct(nCon, 'Type', 'Name', 't', 'y', 'x', 'sol');
intOpts = opts;

for iCon = 1:nCon
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};

    % Integrate system
    if verbose; fprintf(['Integrating system for ' con(iCon).Name '...']); end
    sol = integrateSysSelect(m, con(iCon), tGet, intOpts);
    if verbose; fprintf('done.\n'); end
    
    % Store results
    sim(iCon).Type   = 'Simulation.MassActionKinetics.SelectPoints';
    sim(iCon).Name   = [m.Name ' in ' con.Name];
    sim(iCon).t      = sol.x;
    sim(iCon).y      = m.C1*sol.y + m.C2*sol.u + repmat(sol.c, 1,numel(tGet));
    sim(iCon).x      = sol.y;
    sim(iCon).sol    = sol;
end

%% Work-down
switch (nargout)
    case 0
        % Save the hold state of the figure
        holdState = ishold;
        % Draw each result
        for iCon = 1:nCon
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
        varargout{3} = sim.x;
    case 4
        varargout{1} = sim.t;
        varargout{2} = sim.y;
        varargout{3} = sim.x;
        varargout{4} = sim.sol;
    otherwise
        error(nargoutchk(0, 4, nargout, 'struct'));
end
end
