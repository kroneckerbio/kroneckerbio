function [varargout] = SimulateSensitivity(m, con, opts)
%SimulateSensitivity Integrate the sensitivities of every species with
%   respect to every parameter over all time in the mass action kinetics
%   framework
%
%   Mathematically: dx/dT = Integral(df/dx * dx/dT + df/dT, t=0:tF)
%   
%   [...] = SimulateSensitivity(m, con, opts)
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
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Indicates the kinetic parameters whose sensitivities are
%           desired
%       .UseICs [ logical matrix nx by nCon | logical vector nx |
%                 positive integer vector {[]} ]
%           Indicates the initial conditions of the state species whose
%           sensitivites are desired
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
%   SimulateSensitivity(m, con, opts)
%   	Plots the sensitivities under each condition
%
%   sim = SimulateSensitivity(m, con, opts)
%   	A vector of structures with each entry being the simulation
%       under one of the conditions.
%       .t [ sorted nonnegative row vector ]
%           Timepoints chosen by the ode solver
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
%       
%   The following are permitted if nCon = 1
%   [t, dydT] = SimulateSensitivity(m, con, tGet, opts)
%     	sim.t, sim.dydT
%
%   [t, dydT, dxdT] = SimulateSensitivity(m, con, tGet, opts)
%     	sim.t, sim.dydT, sim.dxdT
%
%   [t, dydT, dxdT, sol] = SimulateSensitivity(m, con, tGet, opts)
%     	sim.t, sim.dydT, sim.dxdT, sim.sol
%       
%	Special
%   This function can also be requested to return an empty array of
%   structures with the same fields as simulation. This may be necessary to
%   initialize an array that the user will later fill.
%
% 	empty = SimulateSensitivity(m)
%   	For integer m, create an empty simulation array m by 1
%
%   empty = SimulateSensitivity([m,n,p...])
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
    varargout{1} = emptystruct(m, 'Type', 'Name', 't', 'dydT', 'dxdT', 'sol');
    return
end

assert(nargin >= 2, 'KroneckerBio:SimulateSensitivity:TooFewInputs', 'SimulateSensitivity requires at least 2 input arguments')
assert(nargout <= 4, 'KroneckerBio:SimulateSensitivity:FourOrFewerOutputs', 'SimulateSensitivity must have between 0 and 4 outputs')
assert(isscalar(m), 'KroneckerBio:SimulateSensitivity:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;

defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseICs         = [];
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

% Ensure UseICs is a logical matrix
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseControls nTq] = fixUseControls(opts.UseControls, opts.UseModelInputs, nCon, m.nq, cat(1,con.nq));

nT = nTk + nTx + nTq;

% Refresh conditions
con = refreshCon(m, con);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, false(nCon,1), nx, nCon, false, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);

%% Run integration for each experiment
sim = emptystruct(nCon, 'Type', 'Name', 't', 'dydT', 'dxdT', 'sol');
intOpts = opts;

for iCon = 1:nCon
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    
    % Integrate dx/dp over time
    if verbose; fprintf(['Integrating sensitivities for ' con(iCon).Name '...']); end
    sol = integrateSens(m, con(iCon), intOpts);
    if verbose; fprintf('done.\n'); end
    
    % Store results
    sim(iCon).Type = 'Simulation.MassActionKineticsSensitivity';
    sim(iCon).Name = [m.Name ' in ' con.Name];
    sim(iCon).t    = sol.x;
    sim(iCon).dydT = @(varargin)evaluateOutputs(sol, varargin{:});
    sim(iCon).dxdT = @(varargin)evaluateStates(sol, varargin{:});
    sim(iCon).sol  = sol;
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
        varargout{2} = sim.dydT;
    case 3
        varargout{1} = sim.t;
        varargout{2} = sim.dydT;
        varargout{3} = sim.dxdT;
    case 4
        varargout{1} = sim.t;
        varargout{2} = sim.dydT;
        varargout{3} = sim.dxdT;
        varargout{4} = sim.sol;
    otherwise
        error(nargoutchk(0, 4, nargout, 'struct'));
end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Evaluation functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function val = evaluateOutputs(sol, t, ind)
        nT = (size(sol.y, 1) - nx) / nx;
        nt = length(sol.x);
        
        if nargin < 3
            val = deval(sol, t, nx+1:nx+nx*nT); % xT_t
            val = reshape(val, nx,nT*nt); % x_Tt
            val = sol.C1*val; % y_Tt
            val = reshape(val, ny*nT,nt); % yT_t
        else
            val = deval(sol, t, nx+1:nx+nx*nT); % xT_t
            val = reshape(val, nx,nT*nt); % x_Tt
            val = sol.C1(ind,:)*val; % y_Tt
            val = reshape(val, nnz(ind)*nT,nt); % yT_t
        end
    end

    function val = evaluateStates(sol, t, ind)
        nT = (size(sol.y, 1) - nx) / nx;
        nt = numel(sol.x);
        
        if nargin < 3
            val = deval(sol, t, nx+1:nx+nx*nT); % xT_t
        else
            val = deval(sol, t, nx+1:nx+nx*nT); % xT_t
            val = reshape(val, nx,nT*nt); % x_Tt
            val = val(ind,:); % x_Tt chopped out rows
            val = reshape(val, nnz(ind)*nT,nt); % xT_t
        end
    end

end