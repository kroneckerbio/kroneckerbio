function [varargout] = SimulateHessian(m, con, opts, xSol, dxdpSol)
%SimulateHessian Integrate the sensitivities of the sensitivities over all
%   time
%
%   Mathematically: d2x/dp2 = Integral(df/dx * d2x/dp2 + (2*d2f/dpdx +
%                               d2p/dx2 * dx/dp) * dx/dp + d2f/dp2)
%   
%   [...] = SimulateHessian(m, con, opts, xSol, dxdpSol)
%   
%   Inputs:
%       m    - The KroneckerBio model that will be simulated
%       con  - An array of experimental conditions under which the model
%              will be simulated
%       opts - Optional function options
%           .useParams - Vector of indexes indicating the rate constants
%                        whose sensitivities will be considered
%           .useICs    - Vector of indexes indicating the initial
%                        concentrations whose sensitivities will be
%                        considered
%           .AbsTol    - Absolute tolerance of the integration
%           .RelTol    - Relative tolerance of the integration
%           .verbose   - Print progress to command window
%       xSol    - The solution to the model simulation. Optional, but
%                 speeds up the calculation. Consider providing this if the
%                 model simulation was needed in addition to these
%                 sensitivities.
%       dxdpSol - The solution to the model sensitivies. Optional, but
%                 speeds up the calculation. It would be useful to provide
%                 this if the sensitivites were needed in addition to these
%                 double sensitivities.
%
%   Outputs:
%       SimulateHessian(m, con, opts, xSol, dxdpSol)
%        - Plots the double sensitivities under each condition
%
%       simulation = SimulateHessian(m, con, opts, xSol, dxdpSol)
%        - An array of structures with each entry being the sensitivities
%          under one of the conditions.
%           simulation.sol - The integrator solution to the sensitivies
%           simulation.t   - A vector of timepoints sampled
%           simulation.y   - A function handle @(t, y) that evaluates the
%                            sensitivities at some particular timepoints
%                            (t) with respect to some particular output
%                            indexes (y)
%           simulation.x   - A function handle @(t, x) that evaluates the
%                            sensitivities at some particular timspoints
%                            (t) with respect to some particular output
%                            indexes (y)
%           simulation.te  - Empty
%           simulation.ye  - Empty
%           simulation.xe  - Empty
%           simulation.ie  - Empty
%       
%       [t, y] = SimulateHessian(m, con, opts, xSol, dxdpSol)
%        - simulation.t, simulation.y
%
%       [t, y, x] = SimulateHessian(m, con, opts, xSol, dxdpSol)
%        - simulation.t, simulation.y, simulation.x
%
%       [t, y, x, sol] = SimulateHessian(m, con, opts, xSol, dxdpSol)
%        - simulation.t, simulation.y, simulation.x, simulation.sol
%       
%	Special:
%       This function can also be requested to return an empty array of
%       structures with the same fields as simulation. This may be
%       necessary to initialize an array that the user will later fill.
%
%       empty = SimulateHessian(m)
%           For integer m, create an empty simulation array m by 1
%
%       empty = SimulateHessian(m,n)
%           For integers m and n, create an empty simulation array m by n
%
%       empty = SimulateHessian([m,n,p...])
%           For integers array, create an empty simulation array size
%           [m,n,p...]

% (c) 2009 David R Hagen and Bruce Tidor
% This software is released under the GNU GPLv2 or later.

%% Work-up
% Clean up inputs
if nargin < 5
    dxdpSol = [];
    if nargin < 4
        xSol = [];
        if nargin < 3
            opts = [];
            if nargin < 2
                con = [];
                if nargin < 1
                    m = [];
                end
            end
        end
    end
end

% Special case: return empty structure array if inputs are numeric
if isnumeric(m)
    assert(nargout <= 1, 'KroneckerBio:SimulateHessian:OneOutputOnEmpty', 'SimulationHessian must have only one output when an empty simulation structure is requested.')
    temp = struct('sol',[],'t',[],'y',[],'x',[],'te',[],'ye',[],'xe',[],'ie',[]);
    if isempty(m)
        m = 1;
    end
    if ~isscalar(m)
        varargout{1}(prod(m),1) = temp;
        varargout{1} = reshape(varargout{1}, m);
        return
    end
    if isempty(con)
        con = 1;
    end
    varargout{1}(m, con) = temp;
    return
end

assert(nargin >= 2, 'KroneckerBio:SimulateHessian:TooFewInputs', 'SimulateHessian requires at least 2 input arguments')
assert(nargout <= 4, 'KroneckerBio:SimulateHessian:FiveOrFewerOutputs', 'SimulationHessian must have between 0 and 4 outputs')
assert(isscalar(m), 'KroneckerBio:SimulateHessian:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;

opts = mergestruct(defaultOpts, opts);

% Constants
nX = m.nX;
nPK = length(opts.useParams);
nPX = length(opts.useICs);
nP = nPK + nPX;

%% Run integration for each experiment
nCon = length(con);
simulation(nCon) = struct('sol',[],'t',[],'y',[],'x',[],'te',[],'ye',[],'xe',[],'ie',[]);

for iCon = 1:nCon
    % *Integrate System*
    % If xSol is not provided, simulate the system. Otherwise, use it.
    if isempty(xSol) || isempty(xSol(iCon).y)
        sim = simulate(m, con(iCon), opts);
        ixSol = sim.sol;
    else
        ixSol = xSol(iCon);
    end
    
    % *Integrate First Derivative*
    % If dxdpSol is not provided, simulate the system. Otherwise, use it.
    if isempty(dxdpSol) || isempty(dxdpSol(iCon).y)
        sim = SimulateSensitivity(m, con(iCon), opts, ixSol);
        idxdpSol = sim.sol;
    else
        idxdpSol = dxdpSol(iCon);
    end
    
    % *Integrate Second Dervative*
    % Initial values are zero
    dx2dp20 = zeros(nX*nP*nP,1);
    
    % Construct the odes and jacobian for integration of the 2nd derivative
    [d2xdp2Der, d2xdp2Jac] = constructSystem(m, con(iCon), ixSol, idxdpSol);
    
    % Get the times that the solver will need to stop
    stopTimes = unique([0, con(iCon).tF, con(iCon).tStop(con(iCon).tStop < con(iCon).tF)]);

    % Integrate d2x/dp2 over time
    if opts.verbose; fprintf(['Integrating hessian for ' con(iCon).name '...']); end
    d2xdp2Sol = accumulateSol(d2xdp2Der, d2xdp2Jac, [], [], dx2dp20, con(iCon), opts.AbsTol, opts.RelTol, stopTimes);
    if opts.verbose; fprintf('done.\n'); end
    
    % Store results
    simulation(iCon).sol	= d2xdp2Sol;
    simulation(iCon).t      = d2xdp2Sol.x;
    simulation(iCon).y      = @(varargin)evaluateOutputs(d2xdp2Sol, varargin{:});
    simulation(iCon).x      = @(varargin)evaluateSpecies(d2xdp2Sol, varargin{:});
    simulation(iCon).te     = d2xdp2Sol.xe;
    if ~isempty(d2xdp2Sol.ye)
        simulation(iCon).ye = m.c*d2xdp2Sol.ye;
    else
        simulation(iCon).ye = [];
    end
    simulation(iCon).xe     = d2xdp2Sol.ye;
    simulation(iCon).ie     = d2xdp2Sol.ie;
end

%% Work-down

switch (nargout)
    case 0
        % Save the hold state of the figure
        holdState = ishold;
        % Draw each result
        for iCon = 1:nCon
            plotExperiment(m, simulation(iCon));
            hold on;
        end
        % Reset the hold state
        if ~holdState
            hold off;
        end
    case 1
        varargout{1} = simulation;
    case 2 
        varargout{1} = simulation.t;
        varargout{2} = simulation.y;
    case 3
        varargout{1} = simulation.t;
        varargout{2} = simulation.y;
        varargout{3} = simulation.x;
    case 4
        varargout{1} = simulation.t;
        varargout{2} = simulation.y;
        varargout{3} = simulation.x;
        varargout{4} = simulation.sol;
end

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Evaluation functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function val = evaluateOutputs(sol, t, ind)
        if nargin < 3
            ind = [];
        end
        
        nP2 = size(sol.y, 1) / m.nX;
        nT = length(sol.x);

        if isempty(ind)
            val = deval(sol, t); % xpp_t
            val = reshape(val, m.nX, nP2*nT); % x_ppt
            val = m.c*val; % y_ppt
            val = reshape(val, m.nY*nP2,nT); % ypp_t
        else
            val = deval(sol, t); % xpp_t
            val = reshape(val, m.nX, nP2*nT); % x_ppt
            val = m.c(ind,:)*val; % y_ppt
            val = reshape(val, length(ind)*nP2,nT); % ypp_t
        end
    end

    function val = evaluateSpecies(sol, t, ind)
        if nargin < 3
            ind = [];
        end
        
        if isempty(ind)
            val = deval(sol, t);
        else
            val = deval(sol, t, ind);
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% The system for integrating d2xdp2 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [d2xdp2Der, d2xdp2Jac] = constructSystem(m, con, xSol, dxdpSol)
        
        Ip2 = speye(nP*nP);
        
        fpUseParams = zeros(1,nX*nPK); % indexes of useParams when f and p are merged
        for i = 1:nPK
            p = opts.useParams(i);
            fpUseParams(1,((i-1)*nX+1):(i*nX)) = ((p-1)*nX+1):(p*nX);
        end
        
        dfdx = m.dfdx;
        d2fdx2 = m.d2fdx2;
        d2fdp2 = @d2fdp2_sub; %non-kronecker
%         d2fdpdx = @d2fdpdx_sub; %unused
        d2fdxdp = @d2fdxdp_sub;
        
        d2xdp2Der = @der;
        d2xdp2Jac = @jac;
        
        function val = der(t, d2xdp2, u)
            u = u(t);
            x = deval(xSol, t, 1:nX);
            dxdp = reshape( deval(dxdpSol, t, 1:nX*nP), nX,nP); % x_p
            val = reshape( d2fdx2(t,x,u), nX*nX,nX) * dxdp; % fx_p : d2f/dx2 * dx/dp
            val = full(val); % fx_p, to guarentee fullness
            val = reshape(val, nX,nX,nP); % f_x_p
            val = permute(val, [1,3,2]); % f_p_x
            val = reshape(val, nX*nP,nX); % fp_x
            val = val + 2 * d2fdxdp(t,x,u); % fp_x :2*d2f/dpdx + d2f/dx2 * dx/dp
            val = val * dxdp; % fp_p : (2*d2f/dpdx + d2f/dx2 * dx/dp) * dx/dp
%            val = reshape( dfdx(t,x,u) * reshape( d2xdp2, nX,nP*nP), nX*nP,nP) + val; % fp_p : df/dx * d2x/dp2 + (2*d2f/dpdx + d2f/dx2 * dx/dp) * dx/dp %kronecker only
            val = reshape( dfdx(t,x,u) * reshape( d2xdp2, nX,nP*nP), nX*nP,nP) + val + d2fdp2(t,x,u); % fp_p : df/dx * d2x/dp2 + (2*d2f/dpdx + d2f/dx2 * dx/dp) * dx/dp + d2f/dp2 %non-kronecker
            val = reshape(val, nX*nP*nP,1); % fpp
        end
        
        function val = jac(t, d2xdp2, u)
            u = u(t);
            x = deval(xSol, t, 1:nX);
            val = kron(Ip2,dfdx(t,x,u));
        end
        
        function val = d2fdxdp_sub(t, x, u)
            val = m.d2fdxdp(t,x,u);
            val = [val(fpUseParams,:); zeros(nX*nPX, nX)];
        end
                
%         %This function is only needed for non-kronecker models, as
%         d2f/dp2 is zero for all kronecker models. Not fully tested.
        function val = d2fdp2_sub(t, x, u)
            val = m.d2fdp2(t,x,u); % fk_k
            % Keep only the columns and rows corresponding to parameters
            val = [val(fpUseParams, opts.useParams) zeros(nX*nPK,nPX);
                   zeros(nX*nPX, nP)                                 ];% fp_p
        end
          
%         % This function can be used instead of d2fdxdp by changing the
%         order of the reshape calls in fwd. Tests showed this way to be
%         slower because m.d2fdpdx must make the matrix full in order to
%         invert the axes.
%         function val = d2fdpdx_sub(t, x, u)
%             val = m.d2fdpdx(t,x,u);
%             val = [val(:, opts.useParams) zeros(nX*nX, nPX)];
%         end
    end
end