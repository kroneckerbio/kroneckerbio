function [GTolRatio DTolRatio] = OptimalAbsTol(m, con, obj, opts)
% Clean up inputs
assert(nargin >= 3, 'KroneckerBio:OptimalAbsTol:TooFewInputs', 'OptimalAbsTol requires at least 3 input arguments')
if nargin < 4
    opts = [];
end

assert(isscalar(m), 'KroneckerBio:OptimalAbsTol:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;

defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseICs         = [];
defaultOpts.UseControls    = [];

defaultOpts.ObjWeights     = ones(size(obj));

defaultOpts.Normalized     = true;
defaultOpts.UseAdjoint     = true;

defaultOpts.Coverage       = 1;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nk = m.nk;
nCon = numel(con);
nObj = size(obj,1);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseICs is a logical matrix
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseControls nTq] = fixUseControls(opts.UseControls, opts.UseModelInputs, nCon, m.nq, cat(1,con.nq));

nT = nTk + nTx + nTq;

%% Integration type: simple, continuous, complex, or both
% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(obj);

%% Tolerances
% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nx, nCon, opts.UseAdjoint, opts.UseParams, opts.UseICs, opts.UseModelICs, true);

%% Determine the order of the integration
if nargout <= 1
    order = 0;
elseif nargout <= 2
    order = 1;
end

intOpts = opts;

GTolRatio = cell(nCon,1);
DTolRatio = cell(nCon,1);

for iCon = 1:nCon
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    intOpts.tGet = opts.tGet{iCon};
    intOpts.ObjWeights = opts.ObjWeights(:,iCon);
    
    if verbose; fprintf(['Computing optimal tolerance ratio for ' con(iCon).Name '...\n']); end
    % Integrate the basic system
    if order == 0
        if opts.continuous(iCon)
            %sol = integrateObjFwd(m, con, obj, intOpts);
        else
            %sol = integrateObjFwdComplex(m, con, obj, intOpts);
        end
    elseif order == 1
        if opts.continuous(iCon)
            %sol = integrateObjSensFwd(m, con, obj, intOpts);
        else
            [GTolRatio{iCon} DTolRatio{iCon}] = integrateOptimalAbsTolSensSimple(m, con, obj, intOpts);
        end
    end
    
end
