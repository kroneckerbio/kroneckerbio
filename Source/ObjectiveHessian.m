function H = ObjectiveHessian(m, con, obj, opts, xSol, dxdpSol, d2xdp2Sol)
warning('KroneckerBio:OutOfDate', 'ObjectiveHessian is radically out of date and needs to be fixed')
%ObjectiveHessian Evaluate the hessian of a set of objective functions
%
%   D = ObjectiveHessian(m, con, obj, opts)
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix ]
%       The objective structures defining the objective functions to be
%       evaluated.
%       .UseModelSeeds [ logical scalar {false} ]
%           Indicates that the model's seed parameters should be used
%           instead of those of the experimental conditions
%       .UseModelInputs [ logical scalar {false} ]
%           Indicates that the model's inputs should be used instead of
%           those of the experimental conditions
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Which kinetic parameters the gradient will be calculated on
%       .UseSeeds [ logical matrix nx by nCon | logical vector nx |
%                   positive integer vector {[]} ]
%           Which seed parameters the gradient will be calculated on
%       .UseControls [ cell vector nCon of logical vectors or positive 
%                      integer vectors | logical vector nq | positive 
%                      integer vector {[]} ]
%           Which input control parameters the gradient will be calculated
%           on
%     	.ObjWeights [ real matrix nObj by nCon {ones(nObj,nCon)} ]
%           Applies a post evaluation weight on each objective function
%           in terms of how much it will contribute to the final objective
%           function value
%       .Normalized [ logical scalar {true} ]
%           Indicates if the gradient should be computed in log parameters
%           space
%    	.UseAdjoint [ logical scalar {false} ]
%           Indicates whether the gradient should be calculated via the
%           adjoint method or the forward method.
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
%       H: [ real vector nT ]
%           The sum of all objective function hessians

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

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
