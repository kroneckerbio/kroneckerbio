function Ts = SampleParameterSpace(m, con, obj, n, opts)
%SampleParameterSpace Sample the parameter space according to the
%   probability
%
%   Ts = SampleParameterSpace(m, con, obj, n, opts)
%
%   Inputs:
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be simulated
%   con: [ experiment struct vector ]
%       The experimental conditions under which the model will be simulated
%   obj: [ objective struct matrix ]
%       The objective structures defining the objective functions to be
%       evaluated. The objective structures must be an information
%       theory-based type. Specifically, they must have obj.p defined.
%   n: [ nonnegative integer scalar ]
%       The number of samples to be returned
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
%       .Normalized [ logical scalar {true} ]
%           Indicates if the optimization should be done in log parameters
%           space
%       .RecomputeFIM [ logical scalar {true} ]
%           Indicates if the Fisher information matrix should be recomputed
%           as the sampler wanders through parameter space
%       .AdaptAcceptance [ logical scalar {true} ]
%           If the acceptance is ratio is outside of LowerAcceptance and
%           UpperAcceptance, then MaxStepSize will be adjusted after each
%           StepsPerCheck so that the acceptance ratio remains optimal even
%           as the shape of the region being explored changes
%       .MaxStepSize [ nonegative scalar {1} ]
%           This is the maximum relative change that is allowed in any
%           parameter in a single step. Because all steps are Guassian,
%           this will be violated with probability 0.05.
%       .StepsPerCheck [ nonnegative scalar  {100} ]
%           The Fisher information matrix will be recalulated and the
%           MaxStepSize adjusted (if applicable) after this many steps
%       .RelTol [ nonnegative scalar {1e-6} ]
%           Relative tolerance of the integration
%       .AbsTol [ cell vector of nonnegative vectors | nonnegative vector |
%                 nonegative scalar {1e-9} ]
%           Absolute tolerance of the integration. If a cell vector is
%           provided, a different AbsTol will be used for each experiment.
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%
%   This function uses the Metropolis-Hastings algorithm to sample the
%   parameter space according to the probability given by the objective
%   functions that the parameter sets are true. The space is assumed to be
%   approximately Guassian, the shape of which is given by the Fisher
%   information matrix. The step direction and magnitude are chosen to be
%   optimal in this Gaussian framework (acceptance ratio of 0.2 for many
%   parameters), but the nonlinear nature of biological systems can prevent
%   efficient movement in the space. The most common source of problems is
%   that large steps relative to the size of the parameter will almost
%   always be rejected. The MaxStepSize option limits the distance that can
%   be traveled in highly uncertain directions. Of course, this means that
%   highly uncertain directions will not be sampled very efficiently.
%
%   The shape of the parameter space can change as the sampler wanders
%   around. To allow for continuous efficient sampling, the Fisher
%   information matrix can be recomputed after a set number of runs
%   (StepsPerCheck). This feature is set via the RecomputeFIM option. This
%   is a minor violation of the Metropolis criterion, so it should be
%   turned off if numerically perfect results are required.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
assert(nargin >= 4, 'KroneckerBio:SampleParameterSpace:TooFewInputs', 'SampleParameterSpace requires at least 4 input arguments')
if nargin < 5
    opts = [];
end

assert(isscalar(m), 'KroneckerBio:SampleParameterSpace:MoreThanOneModel', 'The model structure m must be scalar')
assert(isscalar(n) && n >= 0 && round(n) == n, 'KroneckerBio:SampleParameterSpace:n', 'The number of samples n must be a nonnegative integer')

%% Options
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

defaultOpts.AdaptAcceptance = true; % Adapt the acceptance ratio
defaultOpts.LowerAcceptance = 0.15; % The lowest acceptance ratio that is good
defaultOpts.UpperAcceptance = 0.5;  % The highest acceptance ratio that is good
defaultOpts.MaxStepSize     = 1;    % Uncertainty vectors with 95% CI stetching more than this fold change are truncated
defaultOpts.StepsPerCheck   = 100;  % Number of steps between checks on the acceptance rate

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

opts.LowerBound = fixBounds(opts.LowerBound, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);
opts.UpperBound = fixBounds(opts.UpperBound, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);

%% Starting parameter set
T = collectActiveParameters(m, con, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);
T = T'; % Sampler must have a row vector

%% Variance of sample space
F = ObjectiveInformation(m, con, obj, opts);
%truncation = 1 ./ (abs(log(opts.InformationTrun)) ./ 2.576).^2; % V_log for 99% CI at bounds
truncation = 1 ./ (abs(log(1 + opts.MaxStepSize)) ./ 1.96).^2; % V_log for 95% CI at bounds
F = posdef(F, truncation);
V = (2.38)^2*infoinv(F)/nT;

%% Run Metropolis-Hastings sampler
Ts = zeros(n, nT);
TsStartIndex = 1; % Starting index of where we are putting in the next draws
TsEndIndex = 0; % Ending index for next draws

isOnCurrent = true;     % True if the next pdf to be calculated is the current position
currentT = T;           % The last parameter set used for the current position
triedT   = currentT;    % The last parameter set that was tried for a new location
currentp = computep(T); % Return this value for the current probability when the move is rejected
triedp   = currentp;    % Return this value for the current probability when the move is accepted
nextT    = T;           % A place to store the next tried parameter set until the current move is over
nextp    = currentp;    % For the probability

while true % dowhile
    % Reset last step
    isOnCurrent = true;
    currentT = T;
    currentp = pdf(T);
    isOnCurrent = false; % mhsample starts on the new parameter set
    
    % Metropolis-Hastings sampler
    [drawn accept] = mhsample(T, opts.StepsPerCheck, 'pdf', @pdf, 'symmetric', true, 'proprnd', @mhrnd);
    
    % Update parameter set
    T = drawn(end,:);

    % Decide how many of this run to keep
    if accept ~= 0 && accept ~= 1 % Don't keep any if it is 0 or 1
        keepstep = ceil(1 / min(accept, 0.305-0.305*accept)); % Assume 0.234 is optimal acceptance ratio
        drawn = drawn(1:keepstep:opts.StepsPerCheck,:);
        
        % Store the draws
        TsEndIndex = min(n, TsStartIndex + size(drawn,1) - 1);
        Ts(TsStartIndex:TsEndIndex,:) = drawn(1:TsEndIndex-TsStartIndex+1,:);
        TsStartIndex = TsEndIndex + 1;
    end
    
    % Terminate when all are drawn
    if TsEndIndex == n
        if verbose; fprintf('%d/%d complete\n', TsEndIndex, n); end
        break
    end
    
    % Update information matrix in this region
    [m, con, obj] = updateAll(m, con, obj, T', opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);
    try
        % If this fails, just use the last information matrix
        F = ObjectiveInformation(m, con, obj, opts);
    end
    if opts.AdaptAcceptance && accept < opts.LowerAcceptance
        % Acceptance is too low, enhance truncation
        if verbose; fprintf('Acceptance ratio %4.2f is too low, truncating uncertainty directions\n', accept); end
        opts.MaxStepSize = opts.MaxStepSize / 2;
    elseif opts.AdaptAcceptance && accept > opts.UpperAcceptance
        % Acceptance is too high, enhance wandering
        if verbose; fprintf('Acceptance ratio %4.2f is too high, reducing truncation\n', accept); end
        opts.MaxStepSize = opts.MaxStepSize * 2;
    else
        if verbose; fprintf('Acceptance ratio %4.2f\n', accept); end
    end
    truncation = 1 ./ (abs(log(1 + opts.MaxStepSize)) ./ 1.96).^2; % V_log for 95% CI at bounds
    F = posdef(F, truncation);
    V = (2.38)^2*infoinv(F)/nT;
    
    if verbose; fprintf('%d/%d complete\n', TsEndIndex, n); end
end

Ts = Ts'; % Kronecker needs column vectors

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p = pdf(T)
        if isOnCurrent
            if all(currentT == T)
                % The last step was rejected, return previous probability
                p = currentp;
            elseif all(triedT == T)
                % The last move was accepted, tried parameter set is now current
                currentT = triedT;
                currentp = triedp;
                
                p = currentp;
            else
                % It shouldn't be possible to get here
                p = computep(T);
                currentT = T;
                currentp = p;
            end
            
            % Move next tried parameters into place
            triedT = nextT;
            triedp = nextp;
        else
            % This is a potential location, store and compute it
            p = computep(T);
            nextT = T;
            nextp = p;
        end
        
        % Alternate step
        isOnCurrent = ~isOnCurrent;
    end

    function p = computep(T)
        T = T'; % Kronecker needs column vectors
        
        % Update everything
        [m, con, obj] = updateAll(m, con, obj, T, opts.UseModelICs, opts.UseModelInputs, opts.UseParams, opts.UseICs, opts.UseControls);
        
        % Compute probability
        try
            p = ObjectiveProbability(m, con, obj, opts);
        catch
            % That parameter set crashed the integrator, return 0
            p = 0;
        end
    end

    function T = mhrnd(T)
        T = T'; % Kronecker needs column vectors
        if opts.Normalized
            T = exp(mvnbndrndgibbs(log(T), V, log(opts.LowerBound), log(opts.UpperBound), log(T), 1, 0, 0));
        else
            T = mvnbndrndgibbs(T, V, opts.LowerBound, opts.UpperBound, T, 1, 0, 0);
        end
        T = T';
    end

end