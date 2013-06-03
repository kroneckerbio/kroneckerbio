function [F All] = ObjectiveInformation(m, con, obj, opts, dxdTSol)
%ObjectiveInformation Compute the Fisher information matrix of a set of
%   objective functions
%
%   Mathematically: F = dy/dp' * V^-1 * dy/dp
%
%   [...] = ObjectiveInformation(m, con, obj, opts, dxdTSol)
%
%   Inputs:
%       m    - The KroneckerBio model that will be simulated
%       con  - A structure vector of the experimental conditions under
%              which the hessian will be evaluated
%       obj  - A structure array of the objective functions under which the
%              hessian will be evaluated. Each row of obj is matched to the
%              corresponding entry in con.
%       opts - Optional function options
%           .UseParams   - Vector of indexes indicating the rate constants
%                          whose sensitivities will be considered
%           .UseICs      - Vector of indexes indicating the initial
%                          concentrations whose sensitivities will be
%                          considered
%           .UseModelICs - Boolean that determines whether to use the
%                          initial concentrations of the model or the
%                          conditions. Default = true
%           .Normalized  - Logical that determines if the simple
%                          information or the normalized information will
%                          be computed. The normalized information is
%                          normalized with respect to the values of the
%                          parameters. Default = true
%           .Verbose     - Print progress to command window
%       dxdTSol   - A structure vector containing the solution to the model
%                   sensitivities under each condition. Optional, but
%                   speeds up the calculation.
%
%   Outputs:
%       F = ObjectiveInformation(m, con, obj, ...)
%           F - A matrix
%
%       [F, All] = ObjectiveInformation(m, con, obj, ...)
%           All - The fisher information matrix (FIM) is the sum of all
%                 fisher information matrices assuming there is no
%                 covariance between errors in seperate experiments. All
%                 provides the individual FIMs for each experiment.
%
%   Additional info:
%   - The experimental condition vector can also be a cell vector
%   - The objective function array can also be a cell array. Empty entries 
%   in the cell array and entries in the structure array with empty
%   values are ignored. This way, conditions can have different numbers of
%   objective functions associated with them.
%   - The optional solutions, dxdTSol, can also be cell vectors. Empty
%   entries in the cell vector and entries in the structure vector with
%   empty values are ignored.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
assert(nargin >= 3, 'KroneckerBio:ObjectiveInformation:AtLeastThreeInputs', 'ObjectiveInformation requires at least 3 input arguments')
if nargin < 5
    dxdTSol = [];
    if nargin < 4
        opts = [];
    end
end

assert(isscalar(m), 'KroneckerBio:ObjectiveInformation:MoreThanOneModel', 'The model structure must be scalar')

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

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nk = m.nk;
nCon = length(con);
nObj = size(obj, 1);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseICs is a logical matrix
[opts.UseICs, nTx] = fixUseICs(opts.UseICs, opts.UseModelICs, nx, nCon);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseControls nTq] = fixUseControls(opts.UseControls, opts.UseModelInputs, nCon, m.nq, cat(1,con.nq));

nT = nTk + nTx + nTq;

% Standardize structures
con = refreshCon(m, con);
obj = refreshObj(m, con, obj, opts.UseParams, opts.UseICs, opts.UseControls);

% Construct starting variable parameter set
T = collectActiveParameters(m, con, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nx, nCon, false, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);

%% Loop through conditions
F = zeros(nT,nT);
Txind = nTk; % Stores the position in F where the first x0 parameter goes for each iCon
Tqind = nTk+nTx; % Stores the position in F where the first q parameter goes for each iCon
intOpts = opts;

% Initialize All array if requested
if nargout >= 2
    All = cell(nObj,nCon);
end

for iCon = 1:nCon
    % If opts.UseModelICs is false, the number of variables can change
    if opts.UseModelICs
        inTx = nTx;
    else
        intOpts.UseICs = opts.UseICs(:,iCon);
        inTx = sum(intOpts.UseICs);
    end
    
    % If opts.UseModelInputs is false, the number of variables can change
    if opts.UseModelInputs
        inTq = nTq;
    else
        intOpts.UseControls = opts.UseControls(iCon);
        inTq = sum(intOpts.UseControls{1});
    end
    
    inT = nTk + inTx + inTq;
    
    % Sensitivity integration if not provided
    if isempty(dxdTSol) || isempty(dxdTSol{iCon})
        if ~opts.UseModelICs
            intOpts.UseICs = opts.UseICs(:,iCon);
        end
        
        % Modify opts structure
        intOpts.AbsTol = opts.AbsTol{iCon};
        tGet = opts.tGet{iCon};
        
        % Integrate
        if verbose; fprintf(['Integrating sensitivities for ' m.Name ' in ' con(iCon).Name '...']); end
        if opts.continuous(iCon) && opts.complex(iCon)
            sol = integrateObjSens(m, con(iCon), obj(:,iCon), intOpts);
        elseif opts.complex(iCon)
            sol = integrateSens(m, con(iCon), intOpts);
        elseif opts.continuous(iCon)
            sol = integrateObjSensSelect(m, con(iCon), obj(:,iCon), tGet, intOpts);
        else
            sol = integrateSensSelect(m, con(iCon), tGet, intOpts);
        end
        if verbose; fprintf('done.\n'); end
    else
        sol = dxdTSol{iCon};
    end
    
    % Sum all FIMs as computed by each objective function
    if opts.Normalized
        % Loop through each objective function in the current row
        for iObj = 1:nObj
            Fi = obj(iObj,iCon).Fn(sol, T);
            % Add to correct place in F
            linInd = [1:nTk, Txind+1:Txind+inTx, Tqind+1:Tqind+inTq]; % linear indexes of parameters
            F(linInd,linInd) = F(linInd,linInd) + Fi;
            
            % Store FIM if requested
            if nargout >= 2
                All{iObj,iCon} = Fi;
            end
        end
    else%~opts.Normalized
        % Loop through each objective function in the current row
        for iObj = 1:nObj
            Fi = obj(iObj,iCon).F(sol);
            % Add to correct place in F
            linInd = [1:nTk, Txind+1:Txind+inTx, Tqind+1:Tqind+inTq]; % linear indexes of parameters
            F(linInd,linInd) = F(linInd,linInd) + Fi;
            
            % Store FIM if requested
            if nargout >= 2
                All{iObj,iCon} = Fi;
            end
        end
    end
    
    % Update condition x0 position
    if ~opts.UseModelICs
        Txind = Txind + inTx;
    end
    % Update condition q position
    if ~opts.UseModelInputs
        Tqind = Tqind + inTq;
    end
end