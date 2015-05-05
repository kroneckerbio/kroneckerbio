function [F, All] = ObjectiveInformation(m, con, obj, opts, dxdTSol)
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
defaultOpts.Verbose          = 1;

defaultOpts.RelTol           = [];
defaultOpts.AbsTol           = [];

defaultOpts.UseParams        = 1:m.nk;
defaultOpts.UseSeeds         = [];
defaultOpts.UseInputControls = [];
defaultOpts.UseDoseControls  = [];

defaultOpts.Normalized     = true;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = m.nx;
nk = m.nk;
ns = m.ns;
n_con = numel(con);
n_obj = size(obj,1);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, ns, n_con);

% Ensure UseControls are cell vectors of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, n_con, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, n_con, cat(1,con.nh));

nT = nTk + nTs + nTq + nTh;

% Refresh conditions and objectives
con = refreshCon(m, con);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nx, n_con, false, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

%% Loop through conditions
F = zeros(nT,nT);

% Initialize All array if requested
if nargout >= 2
    provide_all = true;
    All = cell(n_obj,n_con);
else
    provide_all = false;
end

for i_con = 1:n_con
    opts_i = opts;
    opts_i.AbsTol = opts.AbsTol{i_con};
    opts_i.UseSeeds = opts.UseSeeds(:,i_con);
    opts_i.UseInputControls = opts.UseInputControls{i_con};
    opts_i.UseDoseControls = opts.UseDoseControls{i_con};

    inTs = nnz(opts_i.UseSeeds);
    inTq = nnz(opts_i.UseInputControls);
    inTh = nnz(opts_i.UseDoseControls);
    
    inT = nTk + inTs + inTq + inTh;
    
    % Sensitivity integration if not provided
    if isempty(dxdTSol) || isempty(dxdTSol{iCon})
        ints = integrateAllSens(m, con(i_con), obj(:,i_con), opts_i);
    else
        ints = dxdTSol{iCon};
    end
    
    % Sum all FIMs as computed by each objective function
    Fi = zeros(inT,inT);
    Allout = cell(n_obj,1);
    if opts.Normalized
        % Loop through each objective function in the current row
        for i_obj = 1:n_obj
            Fij = obj(i_obj,i_con).Fn(ints(i_obj));
            Fi = Fi + Fij;
            
            % Store FIM if requested
            if provide_all
                Allout{i_obj} = Fij;
            end
        end
    else%~opts.Normalized
        % Loop through each objective function in the current row
        for i_obj = 1:n_obj
            Fij = obj(i_obj,experimentZero([nCon,nTop])).F(ints(i_obj));
            Fi = Fi + Fij;
            
            % Store FIM if requested
            if provide_all
                Allout{i_obj} = Fij;
            end
        end
    end
    
    if provide_all
        All(:,iCon) = Allout;
    end
    
    % Determine where this experiment's parameters go among all experiments
    linInd = zeros(inT,1);
    
    % Rate parameters are always the same
    linInd(1:nTk) = 1:nTk;
    
    % s parameters are different
    Tsind = nnz(opts.UseSeeds(:,1:i_con-1));
    linInd(nTk+1:nTk+inTs) = nTk+Tsind+1:nTk+Tsind+inTs;
    
    % q parameters are different
    Tqind = nnz(cat(1,opts.UseInputControls{1:i_con-1}));
    linInd(nTk+inTs+1:nTk+inTs+inTq) = nTk+nTs+Tqind+1:nTk+nTs+Tqind+inTq;
    
    % h parameters are different
    Thind = nnz(cat(1,opts.UseDoseControls{1:i_con-1}));
    linInd(nTk+inTs+inTq+1:nTk+inTs+inTq+inTh) = nTk+nTs+nTq+Thind+1:nTk+nTs+nTq+Thind+inTq;

    % Reshape Fi
    Fout = zeros(nT,nT);
    Fout(linInd,linInd) = Fi;
    
    % Add this information to the total
    F = F + Fout;
end
