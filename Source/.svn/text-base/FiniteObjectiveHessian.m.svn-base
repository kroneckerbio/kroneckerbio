function [H, All] = FiniteObjectiveHessian(m, con, obj, opts)
%FiniteObjectiveHessian Compute the hessian of an objective function by
%   the finite difference approximation
%
%   Mathematically: H = dG/dp2 or H = dG/dlnp2
%
%   [H, All] = FiniteObjectiveHessian(m, con, obj, opts)
%
%   Inputs:
%       m    - The KroneckerBio model that will be simulated
%       con  - A structure vector of the experimental conditions under
%              which the hessian will be evaluated
%       obj  - A structure array of the objective functions under which the
%              hessian will be evaluated. Each row of obj is matched to the
%              corresponding entry in con
%       opts - Optional function options
%           .useParams   - Vector of indexes indicating the rate constants
%                          whose sensitivities will be considered
%           .useICs      - Vector of indexes indicating the initial
%                          concentrations whose sensitivities will be
%                          considered
%           .UseModelICs - Boolean that determines whether to use the
%                          initial concentrations of the model or the
%                          conditions. Default = true
%           .Normalized  - Boolean that determines if the simple hessian or
%                          the normalized hessian will be computed. The
%                          normalized hessian is normalized with respect to
%                          the values of the parameters. Default = true
%           .Verbose     - Print progress to command window
%   Outputs:
%       H = FiniteObjectiveHessian(m, con, obj, ...)
%           H - The hessian as an array
%
%       [H, All] = FiniteObjectiveHessian(m, con, obj, ...)
%           All - The hessian is the sum of hessians from all experiements,
%           but All returns the hessians in a size(obj) cell array
%
%   Additional info:
%   - The experimental condition vector can also be a cell vector
%   - The objective function array can also be a cell array. Empty entries 
%   in the cell array and entries in the structure array with empty
%   values are ignored. This way, conditions can have different numbers of
%   objective functions associated with them.

% (c) 2009 David R Hagen and Bruce Tidor
% This software is released under the GNU GPLv2 or later.

% TODO: get 0s instead of NaNs when parameters are zero when using
% Normalized.
% TODO: Updating of con and obj when nargin=1 appears to be incomplete

%% Work-up
% Clean up inputs
assert(nargin >= 3, 'KroneckerBio:FiniteObjectiveHessian:AtLeastThreeInputs', 'FiniteObjectiveHessian requires at least 3 input arguments')
if nargin < 4
	opts = [];
end

assert(isscalar(m), 'KroneckerBio:FiniteObjectiveHessian:MoreThanOneModel', 'The model structure must be scalar')

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

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

opts = mergestruct(defaultOpts, opts);

% Constants
nx = m.nx;
nk = m.nk;
nCon = numel(con);
nObj = size(obj, 2);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);
UseParamsInd = find(opts.UseParams);

% Ensure UseICs is a logical matrix
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);
UseICsInd = find(opts.UseParams);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseControls nTq] = fixUseControls(opts.UseControls, opts.UseModelInputs, nCon, m.nq, cat(1,con.nq));

nT = nTk + nTx + nTq;

assert(nCon == 1 || opts.UseModelICs || ~nTx, 'KroneckerBio:FiniteObjectiveGradient:VectorCon', 'A vector for con is not yet fully supported variable initial conditions. Please contribute!')

% Add missing fields to structure
con = pastestruct(Uzero(m), con);
obj = pastestruct(Gzero(m), obj);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(obj);

%% Tolerances
% RelTol
if isempty(opts.RelTol) || isnan(opts.RelTol)
    opts.RelTol = 1e-6;
end

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 1, opts.continuous, nx, nCon);

%% Loop through conditions
H = zeros(nT,nT);

if nargout <= 1
    % Initial value
    if verbose; fprintf('Initial round\n'); end
    [unused D] = computeObjSens(m, con, obj, opts);
    
    % Finite stepping
    for i = 1:nT
        if verbose; fprintf('Step %d of %d\n', i, nk); end
        
        % Set baseline parameters
        kup = m.p;
        kdown = m.p;
        
        if opts.UseModelICs
            xup = m.ic;
            xdown = m.ic;
        else
            xup = con.ic;
            xdown = con.ic;
        end
        
        % Change current parameter by finite amount
        if i <= nTk
            pi = kup(UseParamsInd(i));
            diff = kup(UseParamsInd(i)) * 1e-8;
            kup(UseParamsInd(i)) = pi + diff;
            kdown(UseParamsInd(i)) = pi - diff;
        else
            pi = xup(UseICsInd(i-nTk));
            diff = xup(UseICsInd(i-nTk)) * 1e-8;
            xup(UseICsInd(i-nTk)) = pi + diff;
            xdown(UseICsInd(i-nTk)) = pi - diff;
        end
        
        % Run models to get goal function
        if opts.UseModelICs
            mtemp = m.update(kup, xup);
            contemp = con.update(m);
            objtemp = obj.update(mtemp);
            Dup = computeObjSens(mtemp, contemp, objtemp, opts);
            mtemp = m.update(kdown, xdown);
            contemp = con.update(m);
            objtemp = obj.update(mtemp);
            Ddown = computeObjSens(mtemp, contemp, objtemp, opts);
        else
            mtemp = m.update(kup, m.ic);
            contemp = con.update(m, xup);
            objtemp = obj.update(mtemp);
            Dup = computeObjSens(mtemp, contemp, objtemp, opts);
            mtemp = m.update(kdown, m.ic);
            contemp = con.update(m, xdown);
            objtemp = obj.update(mtemp);
            Ddown = computeObjSens(mtemp, contemp, objtemp, opts);
        end
        
        % Compute H
        if opts.Normalized
            H(i,:) = pi * p .* ( (Dup - D) / diff + (D - Ddown) / diff ) / 2;
        else
            H(i,:) = ( (Dup - D) / diff + (D - Ddown) / diff ) / 2;
        end
    end
    
else
    % Initialize All array
    All = cell(nCon,nObj);
    
    % Initial value
    for iCon = 1:nCon
        for iObj = 1:nObj
            if verbose; fprintf('Initial round\n'); end
            [unused D] = computeObjSens(m, con(iCon), obj(iCon,iObj), opts);
            
            % Finite stepping
            curH = zeros(nT,nT);
            
            for i = 1:nT
                if verbose; fprintf('Step %d of %d\n', i, nk); end
                
                % Set baseline parameters
                kup = m.p;
                kdown = m.p;
                
                if opts.UseModelICs
                    xup = m.ic;
                    xdown = m.ic;
                else
                    xup = con{iCon}.ic;
                    xdown = con{iCon}.ic;
                end
                
                % Change current parameter by finite amount
                if i <= nTk
                    pi = kup(UseParamsInd(i));
                    diff = kup(UseParamsInd(i)) * 1e-8;
                    kup(UseParamsInd(i)) = pi + diff;
                    kdown(UseParamsInd(i)) = pi - diff;
                else
                    pi = xup(UseICsInd(i-nTk));
                    diff = xup(UseICsInd(i-nTk)) * 1e-8;
                    xup(UseICsInd(i-nTk)) = pi + diff;
                    xdown(UseICsInd(i-nTk)) = pi - diff;
                end
                
                % Run models to get goal function
                if opts.UseModelICs
                    mtemp = m.update(kup, xup);
                    contemp = con(iCon).update(m);
                    objtemp = obj(iCon,iObj).update(mtemp);
                    Dup = computeObjSens(mtemp, contemp, objtemp, opts);
                    mtemp = m.update(kdown, xdown);
                    contemp = con(iCon).update(m);
                    objtemp = obj(iCon,iObj).update(mtemp);
                    Ddown = computeObjSens(mtemp, contemp, objtemp, opts);
                else
                    mtemp = m.update(kup, m.ic);
                    contemp = con(iCon).update(m, xup);
                    objtemp = obj(iCon,iObj).update(mtemp);
                    Dup = computeObjSens(mtemp, contemp, objtemp, opts);
                    mtemp = m.update(kdown, m.ic);
                    contemp = con(iCon).update(m, xdown);
                    objtemp = obj(iCon,iObj).update(mtemp);
                    Ddown = computeObjSens(mtemp, contemp, objtemp, opts);
                end
                
                % Compute H
                if opts.Normalized
                    curH(i,:) = pi * p .* ( (Dup - D) / diff + (D - Ddown) / diff ) / 2;
                else
                    curH(i,:) = ( (Dup - D) / diff + (D - Ddown) / diff ) / 2;
                end
            end
            
            H = H + curH;
            
            All{iCon, iObj} = curH;
        end
    end
end
