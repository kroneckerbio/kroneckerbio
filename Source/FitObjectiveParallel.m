function [m, con, G, D] = FitObjectiveParallel(m, con, obj, opts)
%FitObjective Optimize the parameters of a model to minimize a set of
%   objective functions
%
%   Mathematically: T = argmin(G(T))
%
%   [m, con, G, D] = FitObjective(m, con, obj, opts)
%
%   FitObjective uses the derivatives in the Kronecker model and in the
%   objective functions to build a function that can not only evaluate the
%   objective functions at particular parameter sets, but also evaluate the
%   gradient at those parameter sets, thereby pointing in the direction of
%   a more optimum parameter set. This function is built around Matlab's
%   fmincon, which is a gradient descent minimizer. It varies the
%   parameters attempting to find the parameter set that will minimize the
%   objective function while keeping with the bounds.
%
%   Inputs
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector n_con ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix n_obj by n_con ]
%       The objective structures defining the objective functions to be
%       evaluated. Note that this matrix must have a number of columns
%       equal to numel(con) (e.g. one objective for each experimental
%       condition is a row vector and multiple objective structures for
%       a single experimental conditions is a column vector).
%   opts: [ options struct scalar {} ]
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
%       .UseInputControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the input control parameters that will be allowed to
%           vary during the optimization
%       .UseDoseControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the input control parameters that will be allowed to
%           vary during the optimization
%       .LowerBound [ nonegative vector {0} ]
%           The lower bound on the fitted parameters. It can be length
%           nk+n_con*nx, nk+nx, nT, just nk if nTx = 0, or a scalar. The
%           bounds will be interpreted in that order if the length matches
%           multiple orders.
%       .UpperBound [ nonegative vector {0} ]
%           The upper bound for the fitted parameters. It must be the same
%           length as LowerBound.
%     	.ObjWeights [ real matrix nObj by nCon {ones(nObj,nCon)} ]
%           Applies a post evaluation weight on each objective function
%           in terms of how much it will contribute to the final objective
%           function value.
%       .Normalized [ logical scalar {true} ]
%           Indicates if the optimization should be done in log parameters
%           space
%    	.UseAdjoint [ logical scalar {false} ]
%           Indicates whether the gradient should be calculated via the
%           adjoint method or the forward method
%     	.TolOptim [ positive scalar {1e-5} ]
%           The objective tolerance. The optimization stops when it is
%           predicted that the objective function cannot be improved more
%           than this in the next iteration.
%     	.Restart [ nonnegative integer scalar {0} ]
%           A scalar integer determining how many times the optimzation
%           should restart once optimization has stopped.
%     	.RestartJump [ handle @(iter,G) returns nonnegative vector nT or
%                      scalar | nonnegative vector nT or scalar {0.001} ]
%           This function handle controls the schedule for the noise that
%           will be added to the parameters before each restart. The
%           parameters for the next iteration will be normally distributed
%           in log space with a mean equal to the previous iteration and a
%           standard deviation equal to the value returned by this
%           function. The value returned should usually be a scalar, but it
%           can also be a vector with length equal to the number of active
%           parameters. It can also be numeric, and the noise will be
%           treated as this constant value.
%      	.TerminalObj [ real scalar {-inf} ]
%           Optimization is halted when this objective function value is
%           reached
%       .MaxStepSize [ nonegative scalar {1} ]
%           Scalar fraction indicator of the maximum relative step size
%           that any parameter can take in a single interation
%     	.Algorithm [ string {active-set} ]
%           Option for fmincon. Which optimization algorithm to use
%     	.MaxIter [ postive scalar integer {1000} ]
%           Option for fmincon. Maximum number of iterations allowed before
%           optimization will be terminated.
%     	.MaxFunEvals [ postive scalar integer {5000} ]
%           Option for fmincon. Maximum number of objective function
%           evaluations allowed before optimization will be terminated.
%       .RelTol [ nonnegative scalar {1e-6} ]
%           Relative tolerance of the integration
%       .AbsTol [ cell vector of nonnegative vectors | nonnegative vector |
%                 nonegative scalar {1e-9} ]
%           Absolute tolerance of the integration. If a cell vector is
%           provided, a different AbsTol will be used for each experiment.
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%       .GlobalOptimization [ logical scalar {false} ]
%           Use global optimization in addition to fmincon
%       .GlobalOpts [ options struct scalar {} ]
%           TODO: API in progress
%       .NumWorkers [ positive integer scalar {[]} ]
%           Number of parallel workers to use. If [], the value is set to
%           the maximum number of workers specified in the 'local' cluster
%           profile, which is typically the number of cores on the machine.
%           If a parallel pool is not initialized, a pool will be
%           initialized on the local cluster with the provided number of
%           workers. If a pool is initialized, but has the wrong number of
%           workers, the existing pool will be closed, and a new pool will
%           be opened with the provided number of workers. If a pool is
%           already open with the provided number of workers, no new pool
%           will be opened.
%
%   Outputs
%       m: [ model scalar ]
%           The model with all the optimum kinetic parameters applied, as
%           well as the IC parameters if UseModelSeeds = true.
%       con: [ experiment vector ]
%           The experimental conditions which will have the optimum IC
%           parameters applied if UseModelSeeds = false.
%       G: [ real scalar ]
%           The optimum objective function value
%       D: [ real vector nT ]
%           The objective gradient at the optimum parameter set

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Check for optimization toolbox
% Note: this isn't necessary for simple local optimzation methods like
%   fminsearch but we don't use that
assert(logical(license('test','optimization_toolbox')), 'KroneckerBio:FitObjectiveParallel:OptimizationToolboxMissing', 'FitObjectiveParallel requires the optimization toolbox for fmincon.')

% Check for parallel toolbox
assert(~isempty(ver('distcomp')), 'KroneckerBio:FitObjectiveParallel:ParallelToolboxMissing', 'FitObjectiveParallel requires the Parallel Computing toolbox.')

% Clean up inputs
assert(nargin >= 3, 'KroneckerBio:FitObjectiveParallel:TooFewInputs', 'FitObjective requires at least 3 input arguments')
if nargin < 4
    opts = [];
end

assert(isscalar(m), 'KroneckerBio:FitObjectiveParallel:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose          = 1;

defaultOpts.RelTol           = [];
defaultOpts.AbsTol           = [];

defaultOpts.Normalized       = true;
defaultOpts.UseParams        = 1:m.nk;
defaultOpts.UseSeeds         = [];
defaultOpts.UseInputControls = [];
defaultOpts.UseDoseControls  = [];

defaultOpts.ObjWeights       = ones(size(obj));

defaultOpts.UseAdjoint       = true;

defaultOpts.LowerBound       = 0;
defaultOpts.UpperBound       = inf;
defaultOpts.Aeq              = [];
defaultOpts.beq              = [];
defaultOpts.TolOptim         = 1e-5;
defaultOpts.Restart          = 0;
defaultOpts.RestartJump      = 0.001;
defaultOpts.TerminalObj      = -inf;

defaultOpts.MaxStepSize      = 1;
defaultOpts.Algorithm        = 'active-set';
defaultOpts.MaxIter          = 1000;
defaultOpts.MaxFunEvals      = 5000;

defaultOpts.GlobalOptimization = false;
defaultOpts.GlobalOpts         = [];

defaultOpts.NumWorkers         = [];

% Global fit default options
defaultGlobalOpts.Algorithm = 'globalsearch';
defaultGlobalOpts.StartPointsToRun = 'bounds-ineqs';
defaultGlobalOpts.nStartPoints = 10;
defaultGlobalOpts.UseParallel = false;
defaultGlobalOpts.MaxIter = 1000; % for pattern search; fix to better default

% Assign default values and make final options structs
opts = mergestruct(defaultOpts, opts);
opts.GlobalOpts = mergestruct(defaultGlobalOpts, opts.GlobalOpts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Check for global optimization toolbox only if global optimization is specified
%   Note: this isn't necesssary for all global optimization methods, but we depend
%   on this functionality for the current implementation
if opts.GlobalOptimization
    error('KroneckerBio:FitObjectiveParallel:GlobalOptimizationNotImplemented', 'Parallel global optimization is not yet implemented.')
    assert(logical(license('test','gads_toolbox')), 'KroneckerBio:FitObjectiveParallel:GlobalOptimizationToolboxMissing', 'Global optimization requires the global optimization (gads) toolbox.')
end

% Constants
nx = m.nx;
nk = m.nk;
ns = m.ns;

% Ensure structures are proper sizes
[con, n_con] = fixCondition(con);
[obj, n_obj] = fixObjective(obj, n_con);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, ns, n_con);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, n_con, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, n_con, cat(1,con.nh));

nT = nTk + nTs + nTq + nTh;

% Ensure Restart is a positive integer
if ~(opts.Restart >= 0)
    opts.Restart = 0;
    warning('KroneckerBio:FitObjectiveParallel:NegativeRestart', 'opts.Restart was not nonegative. It has been set to 0.')
end

if ~(opts.Restart == floor(opts.Restart))
    opts.Restart = floor(opts.Restart);
    warning('KroneckerBio:FitObjectiveParallel:NonintegerRestart', 'opts.Restart was not a whole number. It has been floored.')
end

% Ensure RestartJump is a function handle
if isnumeric(opts.RestartJump)
    opts.RestartJump = @(iter,G)(opts.RestartJump);
end

% Construct starting variable parameter set
T0 = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nx, n_con, opts.UseAdjoint, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

% Bounds
opts.LowerBound = fixBounds(opts.LowerBound, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
opts.UpperBound = fixBounds(opts.UpperBound, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

%% Options structure for integration
intOpts = opts;

%% Local optimization options
localOpts = optimoptions('fmincon');
localOpts.Algorithm               = opts.Algorithm;
localOpts.TolFun                  = opts.TolOptim;
localOpts.TolX                    = 0;
localOpts.OutputFcn               = @isTerminalObj;
localOpts.GradObj                 = 'on';
localOpts.Hessian                 = 'off'; % unused
localOpts.MaxFunEvals             = opts.MaxFunEvals;
localOpts.MaxIter                 = opts.MaxIter;
localOpts.RelLineSrchBnd          = opts.MaxStepSize;
localOpts.RelLineSrchBndDuration  = Inf;
localOpts.TolCon                  = 1e-6;

if verbose
    localOpts.Display = 'iter';
else
    localOpts.Display = 'off';
end

%% Global optimization options
% TODO: make sure options are relevant for solver
globalOpts = opts.GlobalOpts;

% Sanity checking
if strcmpi(globalOpts.Algorithm, 'multistart') && globalOpts.UseParallel
    warning('KroneckerBio:FitObjectiveParallel:InvalidMultistartOpts', 'Using multistart with UseParallel is not supported at this time (due to global variable in obj fun usage).')
end

%% Normalize parameters
if opts.Normalized
    % Normalize starting parameters and bounds
    T0 = log(T0);
    opts.LowerBound = log(opts.LowerBound);
    opts.UpperBound = log(opts.UpperBound);
    
    % Change relative line search bound to an absolute scale in log space
    % Because fmincon lacks an absolute option, this hack circumvents that
    localOpts.TypicalX = zeros(nT,1) + log(1 + opts.MaxStepSize)*log(realmax);
    localOpts.RelLineSrchBnd = 1 / log(realmax);
end

%% Apply bounds to starting parameters before optimizing
% fmincon will choose a wierd value if a starting parameter is outside the bounds
T0(T0 < opts.LowerBound) = opts.LowerBound(T0 < opts.LowerBound);
T0(T0 > opts.UpperBound) = opts.UpperBound(T0 > opts.UpperBound);

%% Abort in rare case of no optimization
if numel(T0) == 0
    [G, D] = objective(T0);
    return
end

%% Run optimization
% Initialize loop
That = T0;
Gbest = inf;
Tbest = T0;

for iRestart = 1:opts.Restart+1
    % Init abort parameters
    aborted = false;
    Tabort = That;
    
    % Run specified optimization
    if opts.GlobalOptimization
        
        % Create local optimization problem
        %   Used as a subset/refinement of global optimization
        localProblem = createOptimProblem('fmincon', 'objective', @objective, ...
            'x0', That, 'Aeq', opts.Aeq, 'beq', opts.beq, ...
            'lb', opts.LowerBound, 'ub', opts.UpperBound, ...
            'options', localOpts);
        
        if opts.Verbose
            fprintf('Beginning global optimization with %s...\n', globalOpts.Algorithm)
        end
        
        switch globalOpts.Algorithm
            case 'globalsearch'
                gs = GlobalSearch('StartPointsToRun', globalOpts.StartPointsToRun);
                [That, G, exitflag] = run(gs, localProblem);
            case 'multistart'
                ms = MultiStart('StartPointsToRun', globalOpts.StartPointsToRun, ...
                    'UseParallel', globalOpts.UseParallel);
                [That, G, exitflag] = run(ms, localProblem, globalOpts.nStartPoints);
            case 'patternsearch'
                psOpts = psoptimset('MaxIter', globalOpts.MaxIter, 'UseParallel', globalOpts.UseParallel);
                [That, G, exitflag] = patternsearch(@objective, That, [], [], ...
                    opts.Aeq, opts.beq, opts.LowerBound, opts.UpperBound, [], psOpts);
            otherwise
                error('KroneckerBio:FitObjectiveParallel:InvalidGlobalOptAlgorithm', 'Global optimization algorithm %s not recognized.', globalOpts.Algorithm)
        end
        
        [~, D] = objective(That); % since global solvers don't return gradient at endpoint
        
    else
        
        % Start parallel pool with the specified number of workers, if necessary
        NumWorkers = opts.NumWorkers;
        NumWorkers = setParallelPoolSize(NumWorkers);
        
        % Split experiments into worker groups
        nCon = numel(con);
        baseNumPerWorker = floor(nCon/NumWorkers);
        remainder = nCon - baseNumPerWorker*NumWorkers;
        numPerWorker = repmat(baseNumPerWorker, NumWorkers, 1);
        numPerWorker(1:remainder) = numPerWorker(1:remainder) + 1;
        icon_worker = cell(NumWorkers, 1);
        endi = 0;
        for ii = 1:NumWorkers
            starti = endi+1;
            endi = sum(numPerWorker(1:ii));
            icon_worker{ii} = (starti:endi)';
        end
            
        % Determine which parameters are fit by which experiments
        T_experiment = zeros(nT, 1);
        T_experiment(1:nTk) = 0;
        nTs_con = sum(opts.UseSeeds,1);
        nTq_con = cellfun(@sum, opts.UseInputControls(:)');
        nTh_con = cellfun(@sum, opts.UseDoseControls(:)');
        nT_con = {nTs_con, nTq_con, nTh_con};
        endi = nTk;
        for i_type = 1:3
            for j_con = 1:nCon
                starti = endi + 1;
                endi = endi + nT_con{i_type}(j_con);
                T_experiment(starti:endi) = j_con;
            end
        end
        
        % Generate objective function parts for each worker
        slave_objectives = cell(1,NumWorkers);
        for ii = 1:NumWorkers
            slave_objectives{ii} = generateSlaveObjective(m, con, obj, opts, intOpts, T_experiment, icon_worker{ii});
        end
        
        % Distribute slave objective functions to workers
        spmd
            slave_objectives = codistributed(slave_objectives);
        end
        
        if opts.Verbose; fprintf('Beginning gradient descent...\n'); end
        [That, G, exitflag, ~, ~, D] = fmincon(@masterObjective, That, [], [], opts.Aeq, opts.beq, opts.LowerBound, opts.UpperBound, [], localOpts);
        
    end
    
    % Check abortion status
    % Abortion values are not returned by fmincon and must be retrieved
    if aborted
        That = Tabort;
        [G, D] = objective(That);
    end
    
    % Re-apply stiff bounds
    That(That < opts.LowerBound) = opts.LowerBound(That < opts.LowerBound);
    That(That > opts.UpperBound) = opts.UpperBound(That > opts.UpperBound);
    
    % See if these parameters are better than previous iterations
    if G < Gbest
        Gbest = G;
        Tbest = That;
    end

    % Terminate if goal has been met
    if G <= opts.TerminalObj
        break
    end
    
    % Jump parameters before restarting
    % Retain the deterministic nature of fitting by fixing the stream
    rng_state = rng;
    prodThat = prod(That); % model fingerprint
    rng(mod(prodThat/eps(prodThat)*iRestart,2^32));
    if opts.Normalized
        That = Tbest + randn(nT,1) .* vec(opts.RestartJump(iRestart,G));
    else
        That = exp(log(Tbest) + randn(nT,1) .* vec(opts.RestartJump(iRestart,G)));
    end
    rng(rng_state);
    
    % Prevent jumps from leaving bounds
    while any(That < opts.LowerBound) || any(That > opts.UpperBound)
        That(That < opts.LowerBound) = 2*(opts.LowerBound(That < opts.LowerBound)) - That(That < opts.LowerBound);
        That(That > opts.UpperBound) = 2*(opts.UpperBound(That > opts.UpperBound)) - That(That > opts.UpperBound);
    end
end

% Unnormalize
if opts.Normalized
    Tbest = exp(Tbest);
end

% Update parameter sets
[m, con] = updateAll(m, con, Tbest, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Halt optimization on terminal goal %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function stop = isTerminalObj(x, optimValues, state)
        if optimValues.fval  < opts.TerminalObj
            aborted = true;
            Tabort = x;
            stop = true;
        else
            stop = false;
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Master objective function, calls and aggregates slave objective functions %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [G,D] = masterObjective(T)
        
        if nargout == 1
            
            spmd
                this_slave_objective = getLocalPart(slave_objectives);
                G_d = this_slave_objective{1}(T);
            end
            
        else
            
            spmd
                this_slave_objective = getLocalPart(slave_objectives);
                [G_d, D_d] = this_slave_objective{1}(T);
            end
           
            % Sum slave objective function gradients to get total gradient
            D = zeros(numel(T),1);
            for wi = 1:NumWorkers
                D = D + D_d{wi};
            end
            
        end
        
        % Sum slave objective function values to get total value
        G = 0;
        for wi = 1:NumWorkers
            G = G + G_d{wi};
        end
        
    end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Slave objective function generator %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function objectivefun = generateSlaveObjective(m, con, obj, opts, intOpts, T_experiment, i_cons)
% T_experiment:
%   nT-by-1 vector of experiment indices indicating which experiment fits
%   each parameter. 0 indicates model parameter fit by all experiments.
% i_cons:
%   Vector of experimental indices to be simulated by the slave function.

    nT = numel(T_experiment);
    
    % If worker has no experiments assigned, assign a 0-valued objective
    % function to it
    if isempty(i_cons)
        % Clear variables to avoid wasting memory in closure
        clear m con obj opts intOpts T_experiment
        
        objectivefun = @emptyobjective;
        
        return
    end

    % Get information on which parameters need to be used for the worker's
    % experiments
    TisWorker = T_experiment == 0 | any(bsxfun(@eq, T_experiment, i_cons(:)'), 2);

    % Filter down input arguments to those relevant to the worker
    con = con(i_cons);
    obj = obj(:, i_cons);
    opts.UseSeeds = opts.UseSeeds(:,i_cons);
    opts.UseInputControls = opts.UseInputControls(i_cons);
    opts.UseDoseControls = opts.UseDoseControls(i_cons);
    opts.ObjWeights = opts.ObjWeights(:,i_cons);
    if iscell(opts.AbsTol)
        opts.AbsTol = opts.AbsTol(i_cons);
    end
    intOpts.UseSeeds = intOpts.UseSeeds(:,i_cons);
    intOpts.UseInputControls = intOpts.UseInputControls(i_cons);
    intOpts.UseDoseControls = intOpts.UseDoseControls(i_cons);
    intOpts.ObjWeights = intOpts.ObjWeights(:,i_cons);
    if iscell(intOpts.AbsTol)
        intOpts.AbsTol = intOpts.AbsTol(i_cons);
    end
    
    % Return worker-specific objective function
    objectivefun = @objective;

    function [G, D] = objective(T)
        
        % Get portion of T that is relevant to the worker
        T_worker = T(TisWorker);
        
        % Reset answers
        G = 0;
        D_worker = zeros(numel(T_worker),1);
        
        % Unnormalize
        if opts.Normalized
            T_worker = exp(T_worker);
        else
            % If fmincon chooses values that violate the lower bounds, force them to be equal to the lower bounds
            T_worker(T_worker < opts.LowerBound) = opts.LowerBound(T_worker < opts.LowerBound);
        end
        
        % Update parameter sets
        [m, con] = updateAll(m, con, T_worker, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
        
        % Integrate system to get objective function value
        if nargout == 1
            G = computeObj(m, con, obj, intOpts);
        end
        
        % Integrate sensitivities or use adjoint to get objective gradient
        if nargout == 2
            [G, D_worker] = computeObjGrad(m, con, obj, intOpts);
        end
        
        % Assign back to vector sized for all workers
        D = zeros(nT, 1);
        D(TisWorker) = D_worker;
    end

    function [G, D] = emptyobjective(T)
        % Function to use if worker has no experiments assigned to it
        G = 0;
        D = zeros(nT, 1);
    end

end
