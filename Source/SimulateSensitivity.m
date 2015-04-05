function [varargout] = SimulateSensitivity(m, con, opts)
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
if nargin < 3
    opts = [];
end

assert(nargin >= 2, 'KroneckerBio:SimulateSensitivity:TooFewInputs', 'SimulateSensitivity requires at least 2 input arguments')
assert(isscalar(m), 'KroneckerBio:SimulateSensitivity:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose          = 1;

defaultOpts.RelTol           = [];
defaultOpts.AbsTol           = [];

defaultOpts.UseParams        = 1:m.nk;
defaultOpts.UseSeeds         = [];
defaultOpts.UseInputControls = [];
defaultOpts.UseDoseControls  = [];

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
[opts.UseSeeds, nTx] = fixUseSeeds(opts.UseSeeds, nx, nCon);

% Ensure UseControls are cell vectors of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, nCon, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, nCon, cat(1,con.nh));

nT = nTk + nTx + nTq + nTh;

% Refresh conditions
con = refreshCon(m, con);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, false(nCon,1), nx, nCon, false, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

%% Run integration for each experiment
sim = emptystruct(nCon, 'Type', 'Name', 't', 'y', 'x', 'dydT', 'dxdT', 'sol');

for iCon = 1:nCon
    % Modify opts structure
    intOpts = opts;
    intOpts.AbsTol = opts.AbsTol{iCon};
    
    UseSeeds_i = opts.UseSeeds(:,iCon);
    intOpts.UseSeeds = UseSeeds_i;
    inTs = nnz(UseSeeds_i);
    
    intOpts.UseInputControls = opts.UseInputControls{iCon};
    inTq = nnz(intOpts.UseInputControls);
    
    intOpts.UseDoseControls = opts.UseDoseControls{iCon};
    inTh = nnz(intOpts.UseDoseControls);
    
    inT = nTk + inTs + inTq + inTh;
    
    % Integrate [x; dx/dT] over time
    if verbose; fprintf(['Integrating sensitivities for ' con(iCon).Name '...']); end
    sol = integrateSens(m, con(iCon), intOpts);
    if verbose; fprintf('done.\n'); end
    
    % Store results
    sim(iCon).Type  = 'Simulation.Sensitivity';
    sim(iCon).Name  = [m.Name ' in ' con(iCon).Name];
    sim(iCon).t     = sol.x;
    sim(iCon).y     = @(t, varargin)evaluateOutputs(sol, t, varargin{:});
    sim(iCon).x     = @(t, varargin)evaluateStates(sol, t, varargin{:});
    sim(iCon).dydT  = @(t, varargin)evaluateOutputSensitivities(sol, t, varargin{:});
    sim(iCon).dxdT  = @(t, varargin)evaluateStateSensitivities(sol, t, varargin{:});
    sim(iCon).sol   = sol;
end

%% Work-down
if nargout == 0
    % Draw each result
    for iCon = 1:nCon
        subplot(nCon,1,iCon)
        plotSensitivityExperiment(m, sim(iCon), 'Linewidth', 2);
    end
else
    varargout{1} = sim;
end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Evaluation functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function val = evaluateOutputs(sol, t, ind)
        nt = numel(t);
        
        if nargin < 3
            val = sol.C1 * deval(sol, t) + sol.C2 * sol.u(t) + repmat(sol.c, 1,nt);
        else
            val = sol.C1(ind,:) * deval(sol,t) + sol.C2(ind,:) * sol.u(t) + repmat(sol.c(ind,:), 1,nt);
        end
    end

    function val = evaluateStates(sol, t, ind)
        if nargin < 3
            val = deval(sol, t);
        else
            val = deval(sol, t, ind);
        end
    end

    function val = evaluateOutputSensitivities(sol, t, ind)
        nt = numel(t);
        
        val = zeros(ny*inT,nt);
        for it = 1:nt
            x = deval(sol, t(it), 1:nx);
            u = sol.u(t(it));
            dydx = sol.dydx(t,x,u);
            if nargin >= 3;
                dydx = dydx(ind,:);
            end
            
            dxdT = reshape(deval(sol, t(it), nx+1:nx+nx*inT), nx,inT);
            
            val(:,it) = reshape(dydx*dxdT, nnz(ind)*inT,1);
        end
    end

    function val = evaluateStateSensitivities(sol, t, ind)
        nx = size(sol.C1,2);
        nT = (size(sol.y, 1) - nx) / nx;
        nt = numel(t);
        
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
