function [F All] = ObjectiveInformation(m, con, obj, opts, dxdTSol)
%ObjectiveInformation Compute the Fisher information matrix of a set of
%   objective functions
%
%   Mathematically: F = E(d2logpdT2(y,Y), y, p_y|T(y,T))
%
%   [...] = ObjectiveInformation(m, con, obj, opts, dxdTSol)
%
%   Inputs:
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix ]
%       The objective structures defining the objective functions to be
%       evaluated.
%   opts: [ options struct scalar ]
%       Optional
%       .UseModelSeeds [ logical scalar {false} ]
%           Indicates that the model's seed parameters should be used
%           instead of those of the experimental conditions. This will
%           determine both which parameters are used for simulation as well
%           as what parameters will be varied in the optimization.
%       .UseModelInputs [ logical scalar {false} ]
%           Indicates that the model's inputs should be used instead of
%           those of the experimental conditions. This will determine both
%           which parameters are used for simulation as well as what
%           parameters will be varied in the optimization.
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Indicates the kinetic parameters that will be allowed to vary
%           during the optimization
%       .UseSeeds [ logical matrix ns by nCon | logical vector ns |
%                   positive integer vector {[]} ]
%           Indicates the seeds that will be allowed to vary during the
%           optimzation. If UseModelSeeds is true then UseSeeds can be a
%           vector of linear indexes or a vector of logicals length of ns.
%           If UseModelSeeds is false then UseSeeds can be a matrix of
%           logicals size ns by nCon. It can also be a vector of length ns,
%           and every experiment will be considered to have the same active
%           seed parameters. It can also be a vector of linear indexes into
%           the ns vector and assumed the same for all conditions.
%       .UseControls [ cell vector nCon of logical vectors or positive 
%                      integer vectors | logical vector nq | positive 
%                      integer vector {[]} ]
%           Indicates the input control parameters that will be allowed to
%           vary during the optimization
%       .RelTol [ nonnegative scalar {1e-6} ]
%           Relative tolerance of the integration
%       .AbsTol [ cell vector of nonnegative vectors | nonnegative vector |
%                 nonegative scalar {1e-9} ]
%           Absolute tolerance of the integration. If a cell vector is
%           provided, a different AbsTol will be used for each experiment.
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%   dxdTSol: [ sensitivity solution struct vector nCon {} ]
%       A structure vector containing the solution to the model
%       sensitivities under each condition can be provided to prevent this
%       method from recalculating it
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

% (c) 2013 David R Hagen & Bruce Tidor
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
defaultOpts.UseModelSeeds  = false;
defaultOpts.UseModelInputs = false;

defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseSeeds       = [];
defaultOpts.UseControls    = [];

defaultOpts.Normalized     = true;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nk = m.nk;
ns = m.ns;
nCon = numel(con);
nObj = size(obj,1);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, opts.UseModelSeeds, ns, nCon);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseControls, nTq] = fixUseControls(opts.UseControls, opts.UseModelInputs, nCon, m.nq, cat(1,con.nq));

nT = nTk + nTs + nTq;

% Refresh conditions and objectives
con = refreshCon(m, con);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nx, nCon, false, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams, opts.UseSeeds, opts.UseControls);

%% Loop through conditions
F = zeros(nT,nT);
Tsind = nTk; % Stores the position in F where the first x0 parameter goes for each iCon
Tqind = nTk+nTs; % Stores the position in F where the first q parameter goes for each iCon

% Initialize All array if requested
if nargout >= 2
    All = cell(nObj,nCon);
end

for iCon = 1:nCon

    intOpts = opts;
    
    % If opts.UseModelSeeds is false, the number of variables can change
    if opts.UseModelSeeds
        UseSeeds_i = opts.UseSeeds;
    else
        UseSeeds_i = opts.UseSeeds(:,iCon);
    end
    intOpts.UseSeeds = UseSeeds_i;
    inTs = nnz(UseSeeds_i);
    
    % If opts.UseModelInputs is false, the number of variables can change
    if opts.UseModelInputs
        UseControls_i = opts.UseControls{1};
    else
        UseControls_i = opts.UseControls{iCon};
    end
    intOpts.UseControls = UseControls_i;
    inTq = nnz(UseControls_i);
    
    inT = nTk + inTs + inTq;
    
    % Modify opts structure
    intOpts.AbsTol = opts.AbsTol{iCon};
    tGet = opts.tGet{iCon};
    
    % Sensitivity integration if not provided
    if isempty(dxdTSol) || isempty(dxdTSol{iCon})
        % Integrate
        if opts.continuous(iCon) && opts.complex(iCon)
            sol = integrateObjSens(m, con(iCon), obj(:,iCon), intOpts);
        elseif opts.complex(iCon)
            sol = integrateSens(m, con(iCon), intOpts);
        elseif opts.continuous(iCon)
            sol = integrateObjSensSelect(m, con(iCon), obj(:,iCon), tGet, intOpts);
        else
            sol = integrateSensSelect(m, con(iCon), tGet, intOpts);
        end
    else
        sol = dxdTSol{iCon};
    end
    
    % Add fields for prior objectives
    sol.UseParams = opts.UseParams;
    sol.UseSeeds = UseSeeds_i;
    sol.UseControls = UseControls_i;
    
    % Sum all FIMs as computed by each objective function
    if opts.Normalized
        % Loop through each objective function in the current row
        for iObj = 1:nObj
            Fi = obj(iObj,iCon).Fn(sol);
            
            % Add to correct place in F
            linInd = [1:nTk, Tsind+1:Tsind+inTs, Tqind+1:Tqind+inTq]; % linear indexes of parameters
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
            linInd = [1:nTk, Tsind+1:Tsind+inTs, Tqind+1:Tqind+inTq]; % linear indexes of parameters
            F(linInd,linInd) = F(linInd,linInd) + Fi;
            
            % Store FIM if requested
            if nargout >= 2
                All{iObj,iCon} = Fi;
            end
        end
    end
    
    % Update condition s position
    if ~opts.UseModelSeeds
        Tsind = Tsind + inTs;
    end
    % Update condition q position
    if ~opts.UseModelInputs
        Tqind = Tqind + inTq;
    end
end
