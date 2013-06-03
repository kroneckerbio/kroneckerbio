function [m, con, G, D] = FitObjective(m, con, obj, opts)
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
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix ]
%       The objective structures defining the objective functions to be
%       evaluated.
%   opts: [ options struct scalar ]
%       Optional
%       .UseModelICs [ logical scalar {false} ]
%           Indicates that the model's initial conditions should be used
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
%       .UseICs [ logical matrix nx by nCon | logical vector nx |
%                 positive integer vector {[]} ]
%           Indicates the initial conditions of the state species that will
%           be allowed to vary during the optimzation. If UseModelICs is
%           true then UseICs can be a vector of linear indexes or a vector
%           of logicals length of nx. If UseModelICs is false then UseICs
%           can be a matrix of logicals size nx by nCon. It can also be a
%           vector of length nx, and every experiment will be considered to
%           have the same active IC parameters. It can also be a vector of
%           linear indexes into the nx vector and assumed the same for all
%           conditions.
%       .UseControls [ cell vector nCon of logical vectors or positive 
%                      integer vectors | logical vector nq | positive 
%                      integer vector {[]} ]
%           Indicates the input control parameters that will be allowed to
%           vary during the optimization
%       .LowerBound [ nonegative vector {0} ]
%           The lower bound on the fitted parameters. It can be length
%           nk+nCon*nx, nk+nx, nT, just nk if nTx = 0, or a scalar. The
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
%
%   Outputs
%       m: [ model scalar ]
%           The model with all the optimum kinetic parameters applied, as
%           well as the IC parameters if UseModelICs = true.
%       con: [ experiment vector ]
%           The experimental conditions which will have the optimum IC
%           parameters applied if UseModelICs = false.
%       G: [ real scalar ]
%           The optimum objective function value
%       D: [ real vector nT ]
%           The objective gradient at the optimum parameter set

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
assert(nargin >= 3, 'KroneckerBio:FitObjective:TooFewInputs', 'FitObjective requires at least 3 input arguments')
if nargin < 4
    opts = [];
end

assert(isscalar(m), 'KroneckerBio:FitObjective:MoreThanOneModel', 'The model structure must be scalar')

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

defaultOpts.LowerBound     = zeros(m.nk+m.nx, 1);
defaultOpts.UpperBound     = inf(m.nk+m.nx, 1);
defaultOpts.TolOptim       = 1e-5;
defaultOpts.AdaptAbsTol    = true;
defaultOpts.Restart        = 0;
defaultOpts.RestartJump    = 0.001;
defaultOpts.TerminalObj    = -inf;

defaultOpts.MaxStepSize    = 1;
defaultOpts.Algorithm      = 'active-set';
defaultOpts.MaxIter        = 1000;
defaultOpts.MaxFunEvals    = 5000;

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

% Ensure Restart is a positive integer
if ~(opts.Restart >= 0)
    opts.Restart = 0;
    warning('KroneckerBio:FitObjective', 'opts.Restart was not nonegative. It has been set to 0.')
end

if ~(opts.Restart == floor(opts.Restart))
    opts.Restart = floor(opts.Restart);
    warning('KroneckerBio:FitObjective', 'opts.Restart was not a whole number. It has been floored.')
end

% Ensure RestartJump is a function handle
if isnumeric(opts.RestartJump)
    opts.RestartJump = @(iter,G)(opts.RestartJump);
end

% Refresh conditions and objectives
con = refreshCon(m, con);
obj = refreshObj(m, con, obj, opts.UseParams, opts.UseICs, opts.UseControls);

% Construct starting variable parameter set
T0 = collectActiveParameters(m, con, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);

% Fix integration type
[opts.continuous, opts.complex, opts.tGet] = fixIntegrationType(con, obj);

% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
opts.AbsTol = fixAbsTol(opts.AbsTol, 2, opts.continuous, nx, nCon, opts.UseAdjoint, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);

% Bounds
opts.LowerBound = fixBounds(opts.LowerBound, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);
opts.UpperBound = fixBounds(opts.UpperBound, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);

%% Options structure for integration
intOpts = opts;

%% Minimization options
minOpts.Algorithm               = opts.Algorithm;
minOpts.TolFun                  = opts.TolOptim;
minOpts.TolX                    = 0;
minOpts.OutputFcn               = @isTerminalObj;
minOpts.GradObj                 = 'on';
minOpts.Hessian                 = 'off'; % unused
minOpts.HessMult                = []; % unused
minOpts.MaxFunEvals             = opts.MaxFunEvals;
minOpts.MaxIter                 = opts.MaxIter;
minOpts.RelLineSrchBnd          = opts.MaxStepSize;
minOpts.RelLineSrchBndDuration  = inf;
minOpts.TolCon                  = 1e-6;

if verbose
    minOpts.Display = 'iter';
else
    minOpts.Display = 'off';
end

%% Normalize parameters
if opts.Normalized
    % Normalize starting parameters and bounds
    T0 = log(T0);
    opts.LowerBound = log(opts.LowerBound);
    opts.UpperBound = log(opts.UpperBound);
    
    % Change relative line search bound to an absolute scale in log space
    % Because fmincon lacks an absolute option, this hack circumvents that
    minOpts.TypicalX = zeros(nT,1) + log(1 + opts.MaxStepSize)*log(realmax);
    minOpts.RelLineSrchBnd = 1 / log(realmax);
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
    
    if opts.Verbose; fprintf('Beginning gradient descent...\n'); end
    [That, G, unused, unused, unused, D] = fmincon(@objective, That, [], [], [], [], opts.LowerBound, opts.UpperBound, [], minOpts);
    
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
    defaultStream = RandStream.getDefaultStream;
    prodThat = prod(That); % model fingerprint
    RandStream.setDefaultStream(RandStream.create('mt19937ar', 'seed', mod(prodThat/eps(prodThat)*iRestart,2^32)));
    if opts.Normalized
        That = Tbest + randn(nT,1) .* vec(opts.RestartJump(iRestart,G));
    else
        That = exp(log(Tbest) + randn(nT,1) .* vec(opts.RestartJump(iRestart,G)));
    end
    RandStream.setDefaultStream(defaultStream);
    
    % Prevent jumps from leaving bounds
    while any(That < opts.LowerBound) || any(That > opts.UpperBound)
        That(That < opts.LowerBound) = 2*(opts.LowerBound(That < opts.LowerBound)) - That(That < opts.LowerBound);
        That(That > opts.UpperBound) = 2*(opts.UpperBound(That > opts.UpperBound)) - That(That > opts.UpperBound);
    end
end

% Unnormalize
if opts.Normalized
    Tbest = exp(Tbest);
    D = D ./ Tbest;
end

% Update parameter sets
[m, con] = updateAll(m, con, [], Tbest, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Objective function for fmincon %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [G, D] = objective(T)
        % Reset answers
        G = 0;
        D = zeros(nT,1);
        
        % Unnormalize
        if opts.Normalized
            T = exp(T);
        else
            % If fmincon chooses negative values, force them to be zero
            T(T < 0) = 0;
        end
        
        % Update parameter sets
        [m, con, obj] = updateAll(m, con, obj, T, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);
        
        % Integrate system to get objective function value
        if nargout == 1
            G = computeObj(m, con, obj, intOpts);
        end
        
        % Integrate sensitivities or use adjoint to get objective gradient
        if nargout == 2
            [G, D] = computeObjSens(m, con, obj, intOpts);
            
            % Normalize gradient
            if opts.Normalized
                D = vec(D) .* T;
            end
        end
    end

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

end
