function [varargout] = SimulateDoublesensitivitySelect(m, con, tGet, opts)
%SimulateDoublesensitivitySelect Integrate the second-order sensitivities 
%   of every species with respect to every parameter over all time in the
%   mass action kinetics framework and return the values at select time
%   points
%   
%   Mathematically: d2x/dT2 = Integral(df/dx * d2x/dT2 +
%                                      2 * d2f/dTx * dx/dT +
%                                      (d2f/dx2 * dx/dT) * dx/dT +
%                                      d2f/dT2, t=0:tF)
%   
%   [...] = SimulateDoublesensitivitySelect(m, con, tGet, opts)
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
%   SimulateDoublesensitivitySelect(m, con, tGet, opts)
%   	Plots the second-rder sensitivities under each condition
%
%   sim = SimulateDoublesensitivitySelect(m, con, tGet, opts)
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
%       .d2ydT2 [ matrix ny*nT by numel(tGet) ]
%           The value of the second order sensitivites of the outputs at
%           each selected time point
%       .d2xdT2 [ matrix nx by numel(tGet) ]
%           The value of the second-order sensitivities of the states at
%           each selected time point
%       .sol [ struct scalar ]
%           The discrete integrator solution to the system
%       
%   The following are permitted if nCon = 1
%   [t, d2ydT2] = SimulateDoublesensitivitySelect(m, con, tGet, opts)
%     	sim.t, sim.d2ydT2
%
%   [t, d2ydT2, d2xdT2] = SimulateDoublesensitivitySelect(m, con, tGet, opts)
%     	sim.t, sim.d2ydT2, sim.d2xdT2
%
%   [t, d2ydT2, d2xdT2, sol] = SimulateDoublesensitivitySelect(m, con, tGet, opts)
%     	sim.t, sim.d2ydT2, sim.d2xdT2, sim.sol
%       
%	Special
%   This function can also be requested to return an empty array of
%   structures with the same fields as simulation. This may be necessary to
%   initialize an array that the user will later fill.
%
% 	empty = SimulateDoublesensitivitySelect(m)
%   	For integer m, create an empty simulation array m by 1
%
%   empty = SimulateDoublesensitivitySelect([m,n,p...])
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
    varargout{1} = emptystruct(m, 'Type', 'Name', 't', 'y', 'x', 'dydT', 'dxdT', 'd2ydT2', 'd2xdT2', 'sol');
    return
end

assert(nargin >= 3, 'KroneckerBio:SimulateSensitivitySelect:TooFewInputs', 'SimulateSensitivitySelect requires at least 3 input arguments')
assert(nargout <= 4, 'KroneckerBio:SimulateSensitivitySelect:FourOrFewerOutputs', 'SimulateSensitivitySelect must have between 0 and 4 outputs')
assert(isscalar(m), 'KroneckerBio:SimulateSensitivitySelect:MoreThanOneModel', 'The model structure must be scalar')

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
opts.AbsTol = fixAbsTol(opts.AbsTol, 3, false(nCon,1), nx, nCon, false, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);

%% Run integration for each experiment
sim = emptystruct(nCon, 'Type', 'Name', 't', 'y', 'x', 'dydT', 'dxdT', 'sol');
intOpts = opts;

for iCon = 1:nCon
    % If opts.UseModelICs is false, the number of variables can change
    if opts.UseModelICs
        inTx = nTx;
    else
        inTx = sum(opts.UseICs(:,iCon));
    end
    inT = nTk + inTx;
    
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    
    % Integrate [x; dx/dT; d2x/dT2] over time
    if verbose; fprintf(['Integrating double sensitivities for ' con(iCon).Name '...']); end
    sol = integrateDblsensSelect(m, con(iCon), tGet, intOpts);
    if verbose; fprintf('done.\n'); end
    
    % Store results
    sim(iCon).Type = 'Simulation.MassActionKineticsDoublesensitivity.SelectPoints';
    sim(iCon).Name = [m.Name ' in ' con.Name];
    sim(iCon).t    = sol.x;
    sim(iCon).y    = bsxfun(@plus, sol.C1*sol.y(1:nx,:) + sol.C2*sol.u, sol.c);
    sim(iCon).x    = sol.y(1:nx,:);
    sim(iCon).dydT = reshape(sol.C1*reshape(sol.y(nx+1:nx+nx*inT,:), nx,inT*numel(sol.x)), ny*inT,numel(sol.x));
    sim(iCon).dxdT = sol.y(nx+1:nx+nx*inT);
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
end