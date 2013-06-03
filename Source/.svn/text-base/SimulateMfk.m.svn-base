function [varargout] = SimulateMfk(m, con, opts)
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

assert(nargin >= 2, 'KroneckerBio:SimulateMfk:TooFewInputs', 'SimulateMfk requires at least 2 input arguments')
assert(nargout <= 4, 'KroneckerBio:SimulateMfk:FourOrFewerOutputs', 'SimulateMfk must have between 0 and 4 outputs')
assert(isscalar(m), 'KroneckerBio:SimulateMfk:MoreThanOneModel', 'The model structure must be scalar')

% Constants
nx = m.nx;
nCon = numel(con);

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;

defaultOpts.V0             = zeros(nx);

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Refresh conditions
con = refreshCon(m, con);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 1, false(nCon,1), nx + nx*nx, nCon);

%% Run integration for each experiment
sim(nCon) = struct('Type', [], 'Name', [], 't', [], 'y', [], 'x', [], 'sol', []);
intOpts = opts;

for iCon = 1:nCon
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};

    % Integrate system
    if verbose; fprintf(['Integrating system for ' con(iCon).Name '...']); end
    sol = integrateMfk(m, con(iCon), intOpts);
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
