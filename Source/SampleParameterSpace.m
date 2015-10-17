function [Ts, data] = SampleParameterSpace(m, con, obj, n, opts)
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
%           nk+nCon*nx, nk+nx, nT, just nk if nTx = 0, or a scalar. The
%           bounds will be interpreted in that order if the length matches
%           multiple orders.
%       .UpperBound [ nonegative vector {0} ]
%           The upper bound for the fitted parameters. It must be the same
%           length as LowerBound.
%       .Normalized [ logical scalar {true} ]
%           Indicates if the optimization should be done in log parameters
%           space
%       .StepsPerCheck [ nonnegative scalar  {100} ]
%           The Fisher information matrix will be recalulated and the
%           MaxSampleStep adjusted (if applicable) after this many steps
%       .MaxSampleStep [ nonegative scalar {1} ]
%           This is the maximum relative change that is allowed in any
%           parameter in a single step. Because all steps are Guassian,
%           this will be violated with probability 0.05.
%       .AdaptMaxSampleStep [ logical scalar {true} ]
%           If the acceptance is ratio is outside of LowerAcceptance and
%           UpperAcceptance, then MaxSampleStep will be adjusted after each
%           StepsPerCheck so that the acceptance ratio remains optimal even
%           as the shape of the region being explored changes
%       .AdaptProposal [ logical scalar {true} ]
%           After every StepsPerCheck the proposal distribution will be
%           recomputed by linearization at the current location.
%       .SampleThinning [ numeric scalar {0.234} ]
%           A number between 0 and 1. This fraction of parameter samples is
%           retained
%       .AdaptThinning [ logical scalar {false} ]
%           After every StepsPerCheck, the thinning is recomputed based on
%           the acceptance ratio. The formula is 
%           SampleThinning = min(accept, 0.305-0.305*accept)
%       .LowerAcceptance [ numeric scalar {0.15} ]
%           A number between 0 and 1. The lowest acceptance ratio that is
%           acceptable when AdaptMaxSampleStep is true
%       .UpperAcceptance [ numeric scalar {0.5} ]
%           A number between 0 and 1. The highest acceptance ratio that is
%           acceptable when AdaptMaxSampleStep is true
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
%   always be rejected. The MaxSampleStep option limits the distance that can
%   be traveled in highly uncertain directions. Of course, this means that
%   highly uncertain directions will not be sampled very efficiently.
%
%   The shape of the parameter space can change as the sampler wanders
%   around. To allow for continuous efficient sampling, the Fisher
%   information matrix can be recomputed after a set number of runs
%   (StepsPerCheck). This feature is set via the RecomputeFIM option. This
%   is a minor violation of the Metropolis criterion, so it should be
%   turned off if numerically perfect results are required.

% (c) 2015 David R Hagen & Bruce Tidor
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
defaultOpts.Verbose            = 1;

defaultOpts.RelTol             = [];
defaultOpts.AbsTol             = [];

defaultOpts.Normalized         = true;
defaultOpts.UseParams          = 1:m.nk;
defaultOpts.UseSeeds           = [];
defaultOpts.UseInputControls   = [];
defaultOpts.UseDoseControls    = [];

defaultOpts.ObjWeights         = ones(size(obj));

defaultOpts.LowerBound         = 0;
defaultOpts.UpperBound         = inf;

defaultOpts.StepsPerCheck      = 100;   % Number of steps between checks on the acceptance rate
defaultOpts.MaxSampleStep      = 1;     % Uncertainty vectors with 95% CI stetching more than this fold change are truncated
defaultOpts.AdaptMaxSampleStep = true;  % Change the max sample step based on the acceptance ratio
defaultOpts.AdaptProposal      = true;  % Re-linearize the system every StepsPerCheck to get a new proposal distribution
defaultOpts.SampleThinning     = 0.234; % The fraction of draws that will be retained after thinning
defaultOpts.AdaptThinning      = false; % Change the thinning based on the acceptance ratio
defaultOpts.LowerAcceptance    = 0.15;  % The lowest acceptance ratio that is good
defaultOpts.UpperAcceptance    = 0.5;   % The highest acceptance ratio that is good

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
ns = m.ns;
nk = m.nk;
n_con = numel(con);

% Ensure UseParams is logical vector
[opts.UseParams, nTk] = fixUseParams(opts.UseParams, nk);

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, ns, n_con);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, n_con, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, n_con, cat(1,con.nh));

nT = nTk + nTs + nTq + nTh;

% Bounds
opts.LowerBound = fixBounds(opts.LowerBound, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
opts.UpperBound = fixBounds(opts.UpperBound, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

%% Starting parameter set
T = collectActiveParameters(m, con, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);

%% Variance of sample space
F = ObjectiveInformation(m, con, obj, opts);
%truncation = 1 ./ (abs(log(1 + opts.MaxSampleStep)) ./ 2.576).^2; % V_log for 99% CI at bounds
truncation = 1 ./ (abs(log(1 + opts.MaxSampleStep)) ./ 1.96).^2; % V_log for 95% CI at bounds
F = posdef(F, truncation);
V = (2.38)^2/nT*infoinv(F); % Assume optimal scaling is 2.38^2/nT

%% Run Metropolis-Hastings sampler
Ts = zeros(nT,n);
TsStartIndex = 1; % Starting index of where we are putting in the next draws
TsEndIndex = 0; % Ending index for next draws

recent_T    = [T,T,T];
recent_logp = repmat(compute_logp(T), [1,3]);

while true % dowhile
    % Metropolis-Hastings sampler
    [drawn, accept] = mhsample(row(T), opts.StepsPerCheck, 'logpdf', @logpdf, 'symmetric', true, 'proprnd', @mhrnd);
    drawn = drawn.';
    
    % Update parameter set
    T = drawn(:,end);
    
    % Decide how much to thin
    if opts.AdaptThinning
        thin_fraction = min(accept, 0.305-0.305*accept);
    else
        thin_fraction = opts.SampleThinning;
    end
    
    keepstep = ceil(1 / thin_fraction);
    
    if isinf(keepstep)
        keepstep = 0;
    end
    
    % Thin the draws
    drawn = drawn(:,1:keepstep:opts.StepsPerCheck);
    
    % Store the draws
    TsEndIndex = min(n, TsStartIndex + size(drawn,2) - 1);
    Ts(:,TsStartIndex:TsEndIndex) = drawn(:,1:TsEndIndex-TsStartIndex+1);
    TsStartIndex = TsEndIndex + 1;
    
    % Adapt max sample step
    if opts.AdaptMaxSampleStep && accept < opts.LowerAcceptance
        % Acceptance is too low, enhance truncation
        if verbose; fprintf('Acceptance ratio %4.2f is too low, truncating uncertainty directions\n', accept); end
        opts.MaxSampleStep = opts.MaxSampleStep / 2;
    elseif opts.AdaptMaxSampleStep && accept > opts.UpperAcceptance
        % Acceptance is too high, enhance wandering
        if verbose; fprintf('Acceptance ratio %4.2f is too high, reducing truncation\n', accept); end
        opts.MaxSampleStep = opts.MaxSampleStep * 2;
    else
        if verbose; fprintf('Acceptance ratio %4.2f\n', accept); end
    end
    
    % Terminate when all are drawn
    if TsEndIndex == n
        if verbose; fprintf('%d/%d complete\n', TsEndIndex, n); end
        break
    end
    
    % Adapt proposal distribution
    if opts.AdaptProposal
        [m, con] = updateAll(m, con, T, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
        
        try
            % If this fails, just use the last information matrix
            F = ObjectiveInformation(m, con, obj, opts);
        catch ME
            if strcmp(ME.identifier, 'KroneckerBio:accumulateOde:IntegrationFailure')
                % That parameter set crashed the integrator
                % Leave F unchanged
            else
                rethrow(ME)
            end
        end
        truncation = 1 ./ (abs(log(1 + opts.MaxSampleStep)) ./ 1.96).^2; % V_log for 95% CI at bounds
        F = posdef(F, truncation);
        V = (2.38)^2/nT*infoinv(F);
    end
    
    if verbose; fprintf('%d/%d complete\n', TsEndIndex, n); end
end

data.FinalAcceptance = accept;
data.FinalMaxSampleStep = opts.MaxSampleStep;
data.FinalProposal = V;

% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function logp = logpdf(T)
        T = vec(T); % Kronecker needs column vectors

        if all(recent_T(:,1) == T)
            % Matches first set of parameters
            logp = recent_logp(1);
        elseif all(recent_T(:,2) == T)
            % Matches second set of parameters
            logp = recent_logp(2);
        elseif all(recent_T(:,3) == T)
            % Matches third set of parameters
            logp = recent_logp(3);
        else
            % New parameter set
            logp = compute_logp(T);
        end
        
        % Shift recent parameter sets
        recent_T(:,1:2) = recent_T(:,2:3);
        recent_T(:,3) = T;
        recent_logp(1:2) = recent_logp(2:3);
        recent_logp(3) = logp;
    end

    function logp = compute_logp(T)
        T = vec(T); % Kronecker needs column vectors
        
        % Update everything
        [m, con] = updateAll(m, con, T, opts.UseParams, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
        
        % Compute probability
        try
            logp = ObjectiveLogLikelihood(m, con, obj, opts);
        catch ME
            if strcmp(ME.identifier, 'KroneckerBio:accumulateOde:IntegrationFailure')
                % That parameter set crashed the integrator, return 0
                logp = -inf;
            else
                rethrow(ME)
            end
        end
    end

    function T = mhrnd(T)
        T = vec(T); % Kronecker needs column vectors
        if opts.Normalized
            T = exp(mvnbndrndgibbs(log(T), V, log(opts.LowerBound), log(opts.UpperBound), log(T), 1, 0, 0));
        else
            T = mvnbndrndgibbs(T, V, opts.LowerBound, opts.UpperBound, T, 1, 0, 0);
        end
        T = row(T); % mhsample needs row vectors
    end
end
