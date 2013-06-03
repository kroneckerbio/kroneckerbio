function [H, All] = ObjectiveHessian(m, con, obj, opts, xSol, dxdpSol, d2xdp2Sol)
%OBJECTIVEHESSIAN Compute the hessian of an objective function.
%
%   Given a model and experimental conditions, this function computes the
%   double derivative of an objective function with respect to the
%   parameters.
%
%   Mathematically: H = dG/dp2 or H = dG/dlnp2
%
%   [...] = ObjectiveHessian(m, con, obj, opts, xSol, dxdpSol, d2xdp2Sol)
%
%   Inputs:
%       m    - The KroneckerBio model that will be simulated
%       con  - A structure vector of the experimental conditions under
%              which the hessian will be evaluated
%       obj  - A structure array of the objective functions under which the
%              hessian will be evaluated. Each row of obj is matched to the
%              corresponding entry in con.
%       opts - Optional function options
%           .useParams  - Vector of indexes indicating the rate constants
%                         whose sensitivities will be considered
%           .useICs     - Vector of indexes indicating the initial
%                         concentrations whose sensitivities will be
%                         considered
%           .Normalized - Logical that determines if the simple hessian
%                         or the normalized hessian will be computed. The
%                         normalized hessian is normalized with respect to
%                         the values of the parameters. Default = true
%           .verbose    - Print progress to command window
%       xSol      - A structure vector containing the solution to the
%                   system under each condition. The solution to the model
%                   simulation. Optional, but speeds up the calculation.
%                   Consider providing this and the following solutions if
%                   the model simulations are needed in addition to the
%                   hessian.
%       dxdpSol   - A structure vector containing the solution to the model
%                   sensitivities under each condition. Optional, but
%                   speeds up the calculation.
%       d2xdp2Sol - A stucture vector containing the solutions to the model
%                   double sensitivities under each condition. Optional, 
%                   but speeds up the calculation.
%   Outputs:
%       H = ObjectiveHessian(m, con, obj, ...)
%           H - An array
%
%       [H, All] = ObjectiveHessian(m, con, obj, ...)
%           All - The hessian is the sum of all hessians, so the cell array
%                 of size(obj) is returned containing each individual
%                 hessian for the entries in obj.
%
%   Additional info:
%   - The experimental condition vector can also be a cell vector
%   - The objective function array can also be a cell array. Empty entries 
%   in the cell array and entries in the structure array with empty
%   values are ignored. This way, conditions can have different numbers of
%   objective functions associated with them.
%   - The optional solutions, xSol, dxdpSol, and d2xdp2, can also be cell
%   vectors. Empty entries in the cell vector and entries in the structure
%   vector with empty values are ignored.

% (c) 2009 David R Hagen and Bruce Tidor
% This software is released under the GNU GPLv2 or later.

%% Work-up
% Clean up inputs
assert(nargin >= 3, 'KroneckerBio:ObjectiveHessian:AtLeastThreeInputs', 'ObjectiveHessian requires at least 3 input arguments.')
if nargin < 7
    d2xdp2Sol = [];
    if nargin < 6
        dxdpSol = [];
        if nargin < 5
            xSol = [];
            if nargin < 4
                opts = [];
            end
        end
    end
end

% Options
defaultOpts.useParams   = 1:m.nP;
defaultOpts.useICs      = [];
defaultOpts.UseModelICs = true;
defaultOpts.Normalized  = true;
defaultOpts.verbose     = false;

opts = mergestruct(defaultOpts, opts);

% Constants
nPK = length(opts.useParams);
nPX = length(opts.useICs);
nP = nPK + nPX;
nCon = length(con);
nObj = size(obj, 2);

% Make conditions consistent as cell vectors
if isstruct(con)
    con = num2cell(con);
end

% Make objectives consistent as cell arrays
if isstruct(obj)
    obj = num2cell(obj);
end

% Make solutions consistent as cell vectors
if isstruct(xSol)
    xSol = num2cell(xSol);
end
if isstruct(dxdpSol)
    dxdpSol = num2cell(dxdpSol);
end
if isstruct(d2xdp2Sol)
    d2xdp2Sol = num2cell(d2xdp2Sol);
end

%% Loop through conditions
H = zeros(nP,nP);

% Initialize All array if requested
if nargout >= 2
    All = cell(nCon,nObj);
end

for iCon = 1:nCon
    % First integration if not provided
    if isempty(xSol) || isempty(xSol{iCon})
        sim = simulate(m, con{iCon}, opts);
        xSoliCon = sim.sol;
    else
        xSoliCon = xSol{iCon};
    end
    
    % Second integration if not provided
    if isempty(dxdpSol) || isempty(dxdpSol{iCon})
        sim = SimulateSensitivity(m, con{iCon}, opts, xSoliCon);
        dxdpSoliCon = sim.sol;
    else
        dxdpSoliCon = dxdpSol{iCon};
    end
    
    % Third Integration if not provided
    if isempty(d2xdp2Sol) || isempty(d2xdp2Sol{iCon})
        sim = SimulateHessian(m, con{iCon}, opts, xSoliCon, dxdpSoliCon);
        d2xdp2SoliCon = sim.sol;
    else
        d2xdp2SoliCon = d2xdp2Sol{iCon};
    end
    
    % Sum all hessians as computed by each objective function
    if opts.Normalized
        if opts.verbose; fprintf(['Computing hessian (normalized) for' con{iCon}.name '...']); end
        % Loop through each objective function in the current row
        for iObj = 1:nObj
            % Ignore empty structures
            if ~isempty(obj{iCon,iObj}) && ~isempty(obj{iCon,iObj}.Hn)
                p = m.p(opts.useParams); % rate constants
                if opts.UseModelICs
                    p = [p; m.ic(opts.useICs)]; % model initial conditions
                else
                    p = [p; con{iCon}.ic(opts.useICs)]; % con initial conditions
                end
                curH = obj{iCon,iObj}.Hn(xSoliCon, dxdpSoliCon, d2xdp2SoliCon, p);
                H = H + curH;
                
                % Store hessian if requested
                if nargout >= 2
                    All{iCon, iObj} = curH;
                end
            end
        end
    else
        if opts.verbose; fprintf(['Computing hessian (non-normalized) for' con{iCon}.name '...']); end
        % Loop through each objective function in the current row
        for iObj = 1:nObj
            % Ignore empty structures
            if ~isempty(obj{iCon,iObj}) && ~isempty(obj{iCon,iObj}.H)
                curH = obj{iCon,iObj}.H(xSoliCon, dxdpSoliCon, d2xdp2SoliCon);
                H = H + curH;
                
                % Store hessian if requested
                if nargout >= 2
                    All{iCon, iObj} = curH;
                end
            end
        end
    end
    if opts.verbose; fprintf('done.\n'); end
end
