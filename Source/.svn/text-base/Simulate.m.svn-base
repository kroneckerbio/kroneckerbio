function [varargout] = Simulate(m, con, opts)
%Simulate Integrate the concentration of every species over time using
%   mass action kinetics
%
%   Mathematically: x = Integral(f, t=0:tF)
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
%       .x [ handle @(t,x) returns matrix numel(x) by numel(t) ]
%           This function handle evaluates some states x of the system at
%           some particular time points t. The user may exclude x, in which
%           case all states are returned.
%       .sol [ odesolver struct scalar ]
%           The integrator solution to the system
%       
%   The following are permitted if nCon = 1
%   [t, y] = Simulate(m, con, opts)
%       sim.t, sim.y
%
%   [t, y, x] = Simulate(m, con, opts)
%     	sim.t, sim.y, sim.x
%
%   [t, y, x, sol] = Simulate(m, con, opts)
%     	sim.t, sim.y, sim.x, sim.sol
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

% (c) 2010 David R Hagen, Joshua F Apgar, Jared E Toettcher, & Bruce Tidor
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
    varargout{1} = emptystruct(m, 'Type', 'Name', 't', 'y', 'x', 'sol');
    return
end

assert(nargin >= 2, 'KroneckerBio:Simulate:TooFewInputs', 'Simulate requires at least 2 input arguments')
assert(nargout <= 4, 'KroneckerBio:Simulate:FiveOrFewerOutputs', 'Simulate must have between 0 and 4 outputs')
assert(isscalar(m), 'KroneckerBio:Simulate:MoreThanOneModel', 'The model structure must be scalar')

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
    sol = integrateSys(m, con(iCon), intOpts);
    if verbose; fprintf('done.\n'); end
    
    % Store results
    sim(iCon).Type   = 'Simulation.MassActionKinectics';
    sim(iCon).Name   = [m.Name ' in ' con.Name];
    sim(iCon).t      = sol.x;
    sim(iCon).y      = @(t, varargin)evaluateOutputs(sol, t, varargin{:});
    sim(iCon).x      = @(t, varargin)evaluateStates(sol, t, varargin{:});
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
% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Evaluation functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = evaluateOutputs(sol, t, ind)
if nargin < 3
    val = sol.C1 * deval(sol, t) + sol.C2 * sol.u(t) + repmat(sol.c, 1,numel(t));
else
    val = sol.C1(ind,:) * deval(sol,t) + sol.C2(ind,:) * sol.u(t) + repmat(sol.c(ind,:), 1,numel(t));
end
end

function val = evaluateStates(sol, t, ind)
if nargin < 3
    val = deval(sol, t);
else
    val = deval(sol, t, ind);
end
end
