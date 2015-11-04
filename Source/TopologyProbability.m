function [pmy, data] = TopologyProbability(m, con, obj, objPriorParams, objPriorSeeds, objPriorInputControls, objPriorDoseControls, opts, log_pyTm, F)
% TopologyProbability Compute the relative probability that each member of
%   a set of topologies is true according to a set of
%   information-theory-based objective functions
%
%   Mathematically: 
%       p_m|y(m,y) = p_y|m(y,m) * p_m(m) / sum(p_y|m(y,i) * p_m(i), i)
%       p_y|m(y,m) = Integral(p_y|T,m(y,T,m), T, -inf, inf)
%
%   pmy = TopologyProbability(m, con, obj, objPrior, opts, F)
%
%   Inputs
%   m: [ model struct vector ]
%       The KroneckerBio topologies that will be evaluated. All topologies
%       must have the same number of inputs, controls, seeds, and outputs.
%   con: [ experiment struct vector ]
%       The experimental conditions under which the topologies will be
%       evaluated
%   obj: [ objective struct matrix ]
%       The objective structures defining the data likelihood functions to
%       be evaluated.
%   objPriorParams: [ objective struct vector nTop | {[]} ]
%       The priors of the kinetic parameters. Optional and can be empty if
%       UseParams specifies no parameters.
%   objPriorSeeds: [ objective struct vector nCon |
%                    objective struct scalar | {[]} ]
%       The priors of the seeds. Optional and can be empty if UseSeeds
%       specifies no seeds.
%   objPriorInputControls: [ objective struct vector nCon |
%                            objective struct scalar | {[]} ]
%       The priors of the input control parameters. Optional and can be
%       empty if UseInputControls specifies no controls.
%   objPriorDoseControls: [ objective struct vector nCon |
%                           objective struct scalar | {[]} ]
%       The priors of the dose control parameters. Optional and can be
%       empty if UseDoseControls specifies no controls.
%   opts: [ options struct scalar | {[]} ]
%       .PriorTopology [ nonnegative vector nTop {zeros(nTop,1) + 1/nTop} ]
%           The discrete prior on the topologies
%       .TopologyMethod [ {'linear'} | 'path' | 'harmonic' ]
%           Linear makes a linear approximation around the maximum a
%           posterior parameters
%           Path is a Monte Carlo method that will converge in a few weeks
%           for a simple model
%           Harmonic is a Monte Carlo method that will converge sometime
%           next millenium for a simple model
%       .NeedFit [ {true} | false ]
%           Do the parameters need to be fit before linearized? This should
%           only be set to false if the parameters were fit beforehand
%       .UseParams [ logical vector nk | positive integer vector {1:nk} ]
%           Indicates the kinetic parameters that will be allowed to vary
%           during the optimization and will be considered free parameters
%           in the topology evaluation
%       .UseSeeds [ logical matrix ns by nCon | logical vector ns |
%                   positive integer vector {[]} ]
%           Indicates the seeds that will be allowed to vary during the
%           optimzation and will be considered free parameters in the
%           topology evaluation. If UseModelSeeds is true then UseSeeds can
%           be a vector of linear indexes or a vector of logicals length of
%           ns. If UseModelSeeds is false then UseSeeds can be a matrix of
%           logicals size ns by nCon. It can also be a vector of length ns,
%           and every experiment will be considered to have the same active
%           seed parameters. It can also be a vector of linear indexes into
%           the ns vector and assumed the same for all conditions.
%       .UseInputControls [ cell vector nCon of logical vectors or positive 
%                           integer vectors | logical vector nq | positive 
%                           integer vector {[]} ]
%           Indicates the input control parameters that will be allowed to
%           vary during the optimization and will be considered free
%           parameters in the topology evaluation
%       .UseDoseControls [ cell vector nCon of logical vectors or positive 
%                          integer vectors | logical vector nq | positive 
%                          integer vector {[]} ]
%           Indicates the dose control parameters that will be allowed to
%           vary during the optimization and will be considered free
%           parameters in the topology evaluation
%       All other options available for FitObjective to fit the model to
%       the data before processing
%   log_pyTm: [ vector nTop | {[]} ]
%       The log likelihood of each topology. This should be supplied only
%       if these values were calculated beforehand.
%   F: [ cell vector nTop of matrix nT by nT | {[]} ]
%       The Fisher information matrices at the linearization point. This
%       should be supplied only if these matrices were calculated
%       beforehand.
%
%   Outputs
%   pmy: [ nonnegative vector nTop ]
%       The posterior probability of each topology given the data

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
assert(nargin >= 3, 'KroneckerBio:TopologyProbability:TooFewInputs', 'TopologyProbability requires at least 3 input arguments')
if nargin < 10
    F = [];
    if nargin < 9
        log_pyTm = [];
        if nargin < 8
            opts = [];
            if nargin < 7
                objPriorDoseControls = [];
                if nargin < 6
                    objPriorInputControls = [];
                    if nargin < 5
                        objPriorSeeds = [];
                        if nargin < 4
                            objPriorParams = [];
                        end
                    end
                end
            end
        end
    end
end

% Constants
nTop = numel(m);
nCon = numel(con);
nObj = size(obj,1);

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = [];
defaultOpts.AbsTol         = [];

defaultOpts.UseParams      = cell(nTop,1); 
for i=1:nTop; defaultOpts.UseParams{i} = 1:m(i).nk;end
defaultOpts.UseSeeds         = [];
defaultOpts.UseInputControls = [];
defaultOpts.UseDoseControls  = [];

defaultOpts.Normalized     = true;
defaultOpts.UseAdjoint     = true;
defaultOpts.LowerBound     = 0;
defaultOpts.UpperBound     = inf;

% Sampling
defaultOpts.MaxSampleStep   = 1;     % Uncertainty vectors with 95% CI stetching more than this fold change are truncated
defaultOpts.AdaptMaxSampleStep = true; % Change the max sample step based on the acceptance ratio
defaultOpts.AdaptProposal   = true;  % Re-linearize the system every StepsPerCheck to get a new proposal distribution
defaultOpts.AdaptThinning   = false; % Change the thinning based on the acceptance ratio
defaultOpts.SampleThinning  = 0.234; % The fraction of draws that will be retained after thinning
defaultOpts.LowerAcceptance = 0.15;  % The lowest acceptance ratio that is good
defaultOpts.UpperAcceptance = 0.5;   % The highest acceptance ratio that is good
defaultOpts.StepsPerCheck   = 100;   % Number of steps between checks on the acceptance rate

% Topology probability
defaultOpts.PriorTopology  = zeros(nTop,1) + 1/nTop; % Uniform prior
defaultOpts.TopologyMethod = 'Linear'; % Harmonic, Path
defaultOpts.NeedFit        = true; % Fit the parameters
defaultOpts.TopologyTol    = 0.05;

% Path sampling
defaultOpts.PathProposalDist = 'Linear'; % The initial proposal distribution for path sampling
defaultOpts.PathSchedule     = 'Variable'; % Linear, Harmonic. How should the change from proposal to the posterior be determined
defaultOpts.VariableScheduleTarget = 0.4; % For PathSchedule = Variable, what should the target variance per step be
defaultOpts.PathEvalMethod   = 'Importance'; % How should the samples be converted into the answer
defaultOpts.PathMinimumSamples = 0;
defaultOpts.PathMaximumSamples = inf;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nk = zeros(nTop,1);
ns = m(1).ns;
nx = zeros(nTop,1);
for iTop = 1:nTop
    nk(iTop) = m(iTop).nk;
    nx(iTop) = m(iTop).nx;
end

%% Standardize structures
if isempty(con)
    % Create a con and obj to fit with the prior
    con = Uzero(1);
    obj = objectiveZero;
    nCon = 1;
    nObj = 1;
end

% Experimental conditions
con = vec(con);

% Objective functions
assert(size(obj,2) == nCon, 'KroneckerBio:TopologyProbability:ObjSize', 'Second dimension of "obj" must be equal to nCon')

%% Active Parameters
% UseParams is same for all topologies if only one is provided
if isnumeric(opts.UseParams)
    opts.UseParams = repmat({opts.UseParams}, [nTop,1]);
end

% Ensure UseParams is vector of logical indexes within a cell array
nTk = zeros(nTop,1);
for iTop = 1:nTop
    [opts.UseParams{iTop}, nTk(iTop)] = fixUseParams(opts.UseParams{iTop}, nk(iTop));
end

% Ensure UseSeeds is a logical matrix
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, ns, nCon);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, nCon, cat(1,con.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, nCon, cat(1,con.nh));

nT = nTk + nTs + nTq + nTh;

%% Priors
if nTk > 0
    assert(numel(objPriorParams) == nTop, ...
        'KroneckerBio:TopologyProbability:ObjPriorParamsSize', 'Input "objPriorParams" must be a vector of length numel(m)')
    objPriorParams = [reshape(objPriorParams, [1,1,nTop]), objectiveZero([1, nCon-1, nTop])];
else
    objPriorParams = objectiveZero([0,nCon,nTop]);
end

if nTs > 0
    assert(numel(objPriorSeeds) == nCon, ...
        'KroneckerBio:TopologyProbability:ObjPriorSeedsSize', 'Input "objPriorSeeds" must be a vector of length numel(con)')
    objPriorSeeds = repmat(reshape(objPriorSeeds, [1,nCon,1]), [1,1,nTop]);
else
    objPriorSeeds = objectiveZero([0,nCon,nTop]);
end

if nTq > 0
    assert(numel(objPriorInputControls) == nCon, ...
        'KroneckerBio:TopologyProbability:ObjPriorControlsSize', 'Input "objPriorInputControls" must be a vector of length numel(con)')
    objPriorInputControls = repmat(reshape(objPriorInputControls, [1,nCon,1]), [1,1,nTop]);
else
    objPriorInputControls = objectiveZero([0,nCon,nTop]);
end

if nTh > 0
    assert(numel(objPriorDoseControls) == nCon, ...
        'KroneckerBio:TopologyProbability:ObjPriorControlsSize', 'Input "objPriorControls" must be a vector of length numel(con)')
    objPriorDoseControls = repmat(reshape(objPriorDoseControls, [1,nCon,1]), [1,1,nTop]);
else
    objPriorDoseControls = objectiveZero([0,nCon,nTop]);
end

% Match dimensions of objPrior to obj
objPrior = [objPriorParams; objPriorSeeds; objPriorInputControls; objPriorDoseControls];

%% Tolerances
% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% Fix AbsTol to be a cell array of vectors appropriate to the problem
tempAbsTol = cell(nTop,1);
for iTop = 1:nTop
    if isnumeric(opts.AbsTol)
        tempAbsTol{iTop} = opts.AbsTol;
    end
    
    if isstruct(opts.AbsTol)
        if size(opts.AbsTol,2) == nTop
            tempAbsTol{iTop} = opts.AbsTol(:,iTop);
        else
            assert(size(opts.AbsTol,2) == 1)
            tempAbsTol{iTop} = opts.AbsTol;
        end
    end
end
opts.AbsTol = tempAbsTol;

%% Fix bounds
if isnumeric(opts.LowerBound)
    opts.LowerBound = repmat({opts.LowerBound}, nTop,1);
end
if isnumeric(opts.UpperBound)
    opts.UpperBound = repmat({opts.UpperBound}, nTop,1);
end

%% Distribute options for topology specific functions
assert(~isfield(opts, 'ObjWeights'), 'ObjWeights is not supported for Topology Probability.')
optsTop = repmat(opts, nTop,1);
for iTop = 1:nTop
    optsTop(iTop).UseParams   = opts.UseParams{iTop};
    optsTop(iTop).AbsTol      = opts.AbsTol{iTop};
    optsTop(iTop).LowerBound  = opts.LowerBound{iTop};
    optsTop(iTop).UpperBound  = opts.UpperBound{iTop};
end

%% Fit
if opts.NeedFit
    for iTop = 1:nTop
        if verbose; fprintf(['Fitting ' m(iTop).Name ' to objectives...\n']); end
        m(iTop) = FitObjective(m(iTop), con, [obj; objPrior(:,:,iTop)], optsTop(iTop));
    end
end

% Store parameter sets that were fit
T = cell(nTop,1);
for iTop = 1:nTop
    T{iTop} = collectActiveParameters(m(iTop), con, optsTop(iTop).UseParams, optsTop(iTop).UseSeeds, optsTop(iTop).UseInputControls, optsTop(iTop).UseDoseControls);
end

%% Choose topology method
if strcmpi(opts.TopologyMethod, 'Linear')
    %% Information
    % Compute the log likelihood if nor provided
    if isempty(log_pyTm)
        log_pyTm = zeros(nTop,1);
        for iTop = 1:nTop
            if verbose; fprintf(['Computing likelihood of ' m(iTop).Name '...\n']); end
            log_pyTm(iTop) = ObjectiveLogLikelihood(m(iTop), con, [obj; objPrior(:,:,iTop)], optsTop(iTop));
        end
    end
    
    % Compute the information if not provided
    if isempty(F)
        F = cell(nTop,1);
        for iTop = 1:nTop
            if verbose; fprintf(['Computing information of ' m(iTop).Name '...\n']); end
            F{iTop} = ObjectiveInformation(m(iTop), con, [obj; objPrior(:,:,iTop)], optsTop(iTop));
        end
    end
    
    % Extract eigenvalues
    lambda = cell(nTop,1);
    for iTop = 1:nTop
        % Eigendecompose
        lambda{iTop} = infoeig(F{iTop});
        
        % Bend flat directions
        lambda{iTop}(lambda{iTop} < 1e-16) = 1e-16;
    end
    
    %% Compute p_y|m for each topology
    % p_y|m will be computed in log space in order to give more digits to the
    % exponent for computing this number, which is often too small to be
    % represented by float64.
    log_pym = zeros(nTop,1); % In log space
    for iTop = 1:nTop
        % Equation for linearization about maximum a posteriori
        log_pym(iTop) = log_pyTm(iTop) + nT(iTop)/2 * log(2*pi) + -1/2 * sum(log(lambda{iTop}));
    end
    
    %% Compute p_m|y for the set
    % Apply topology prior
    log_pmy = log_pym + log(opts.PriorTopology);
    
    % Rescale p_y|m
    % Since the entire distribution is normalized at the end, this operation
    % has no mathematical effect on the answer. However, it ensures that all
    % the probabilities are representable by float64 when numbers are
    % exponentiated back into regular space.
    log_pmy = log_pmy - max(log_pmy);
    
    % Return to regular space
    pmy = exp(log_pmy);
    
    % Normalize
    pmy = pmy ./ sum(pmy);
    
elseif strcmpi(opts.TopologyMethod, 'Harmonic')
    %% Sample
    % Start from the best fit parameter set
    nSam = opts.StepsPerCheck;
    Ts = cell(nTop,1);
    pym = zeros(nTop,1);
    pyTm = cell(nTop,1);
    pyTminv = cell(nTop,1);
    needsSampled = true(nTop,1); % Can this topology be sampled
    
    step = ones(nTop,1);
    unc_log_pym = zeros(nTop,1);
    
    for iTop = 1:nTop
        Ts{iTop} = collectActiveParameters(m(iTop), con, optsTop(iTop).UseModelSeeds, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseSeeds, optsTop(iTop).UseControls);
        pyTm{iTop} = ObjectiveProbability(m(iTop), con, obj(:,:,iTop), optsTop(iTop)); % p_y|mT
        pyTminv{iTop} = pyTm{iTop}.^-1;
        needsSampled(iTop) = logical(pyTm{iTop} * ObjectiveProbability(m(iTop), Uzero(m(iTop)), objPrior(iTop), optsTop(iTop))) & (nT(iTop) >= 1); % p_T|y
    end
    
    while true % dowhile
        % Draw parameter sets for each topology
        for iTop = find(needsSampled)'
            % Sample the posterior parameter space
            mtemp = updateAll(m(iTop), con, Ts{iTop}(:,end), optsTop(iTop).UseModelSeeds, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseSeeds, optsTop(iTop).UseControls);
            Tsi = SampleParameterSpace(mtemp, con(:,iTop), [obj(:,:,iTop); objPrior(:,:,iTop)], nSam, optsTop(iTop));
            Ts{iTop} = [Ts{iTop}, Tsi];
            
            % Compute the likelihood
            pyTmi = zeros(nSam,1);
            for i_sample = 1:nSam
                mtemp = updateAll(m(iTop), con(:,iTop), Tsi(:,i_sample), optsTop(iTop).UseModelSeeds, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseSeeds, optsTop(iTop).UseControls);
                pyTmi(i_sample) = ObjectiveProbability(mtemp, con, obj, optsTop(iTop));
            end
            pyTm{iTop} = [pyTm{iTop}; pyTmi];
            pyTminv{iTop} = [pyTminv{iTop}; pyTmi.^-1];
            
            % Compute autocorrelation in data
            first = zeros(nT(iTop),1);
            for iT = 1:nT(iTop)
                acf = autocorrelation(Ts{iTop}(iT,:));
                
                % Find first non-autocorrelation point
                %first(iT) = find(acf > 0.05, 1, 'last') + 1; % Tends to pick up noise near the end
                found = find(acf < 0.05, 1, 'first');
                if isempty(first)
                    first(iT) = size(Ts{iTop},2);
                else
                    first(iT) = found;
                end
            end
            step(iTop) = max(first);
        end
        
        % Uncertainty in estimates
        for iTop = 1:nTop
            % Marginal likelihood
            pym(iTop) = mean(pyTminv{iTop}).^-1;
            
            % Compute absolute standard error
            if pyTm{iTop}(1) > 0
                % s^2_H = p_y|m^4 * var(1/p_y|Tm) / (n-1)
                unc_log_pym(iTop) = sqrt(pym(iTop) .^ 4 * var(pyTminv{iTop}(1:step:end),1) / (numel(1:step:numel(pyTminv{iTop})) - 1));
            else
                unc_log_pym(iTop) = 0;
            end
        end
        
        % Error relative to largest probability
        relerr = unc_log_pym ./ max(pym.* opts.PriorTopology);

        % Stop sampling those with low error
        needsSampled(~isnan(relerr) & relerr < opts.TopologyTol) = false;
        
        % Stop if uncertainty is low enough
        if all(~needsSampled)
            break
        end
    end
    
    % Include topology prior
    pmy = pym .* opts.PriorTopology;
    
    % Normalize
    pmy = pmy ./ sum(pmy);
    
elseif strcmpi(opts.TopologyMethod, 'Path')
    %% Proposal distribution
    if strcmpi(opts.PathProposalDist, 'Prior')
        % Initial proposal distribution is the prior
        objApprox = objPrior;
    elseif strcmpi(opts.PathProposalDist, 'Linear')
        % Initial proposal distribution is the linearization at the maximum a posteriori parameters
        objApprox = objectiveZero([1,nCon,nTop]);
        for iTop = 1:nTop
            % Fisher information for posterior
            F_approx = ObjectiveInformation(m(iTop), con, [obj; objPrior(:,:,iTop)], optsTop(iTop));
            
            if opts.Normalized
                V_approx = zeros(m(iTop).nk,m(iTop).nk);
                V_approx(optsTop(iTop).UseParams,optsTop(iTop).UseParams) = infoinv(F_approx);
                V_approx = spdiags(m(iTop).k, 0, m(iTop).nk, m(iTop).nk) * V_approx * spdiags(m(iTop).k, 0, m(iTop).nk, m(iTop).nk);
                
                objApprox(1,1,iTop) = objectiveLogNormalPriorOnKineticParameters(m(iTop).k, V_approx);
            else
                V_approx = zeros(m(iTop).nk,m(iTop).nk);
                V_approx(optsTop(iTop).UseParams,optsTop(iTop).UseParams) = infoinv(F_approx);
                
                objApprox(1,1,iTop) = objectiveNormalPriorOnKineticParameters(m(iTop).k, V_approx);
            end
        end
    else
        error('KroneckerBio:TopologyProbability:PathProposalDist', [opts.PathProposalDist ' is not valid for opts.PathProposalDist; must be "Prior" or "Linear"'])
    end
    
    %% Path Schedule
    if strcmpi(opts.PathSchedule, 'Linear')
        % Number of bridges
        % The span of each bridge is equal to the tolerance to cause the
        % step size to go down as the tolerance does
        n_bridges = ceil(opts.TopologyTol.^-1);
        bridge_points = repmat({linspace(0, 1, n_bridges)}, nTop,1);
        
        % Initialize MCMC sampler
        samples = repmat({cell(n_bridges,1)}, nTop,1);
        sample_step_size = repmat({zeros(nTop,n_bridges)}, nTop,1);
        for iTop = 1:nTop
            % Stores the last parameter set of the previous sampling as a
            % starting point for the next run of sampling
            % Initialize with the best fit parameters
            T_previous = collectActiveParameters(m(iTop), con, optsTop(iTop).UseModelSeeds, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseSeeds, optsTop(iTop).UseControls);

            for i_bridge = n_bridges:-1:1
                if verbose; fprintf('Sampling for autocorrelation for topology %d, bridge %d...\n', iTop, i_bridge); end
                % Scale factor for this bridge point
                current_bridge_point = bridge_points{iTop}(i_bridge);
                
                [samples_i_bridge, step_i_bridge] = determine_deautocorrelation_step_size_on_geometric_mix(m(iTop), con, objApprox(:,:,iTop), [obj; objPrior(:,:,iTop)], optsTop(iTop), current_bridge_point, T_previous, opts.TopologyTol, opts.TopologyTol);
                if verbose; fprintf('Autocorrelation dies after %d steps...\n', step_i_bridge); end

                % Store sampling results
                samples{iTop}{i_bridge} = samples_i_bridge;
                sample_step_size{iTop}(i_bridge) = step_i_bridge;

                % Set loop values for next iteration
                T_previous = samples_i_bridge(:,end);
            end
        end
    elseif strcmpi(opts.PathSchedule, 'Harmonic')
        error('Harmonic path schedule is not implemented')
    elseif strcmpi(opts.PathSchedule, 'Variable')
        % Initialize smart sampler
        bridge_points = cell(nTop,1);
        samples = cell(nTop,1);
        sample_step_size = cell(nTop,1);
        max_sample_step = cell(nTop,1);
        for iTop = 1:nTop
            samples{iTop} = cell(0,1);
            
            % Updated with the latest bridge point that needs sampled and will
            % be used to find the next bridge point
            % Initialize at beta = 0
            current_bridge_point = 0;
            
            % Stores the last parameter set of the previous sampling as a
            % starting point for the next run of sampling
            % Initialize with the best fit parameters
            T_previous = collectActiveParameters(m(iTop), con, optsTop(iTop).UseModelSeeds, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseSeeds, optsTop(iTop).UseControls);
            
            % Step up until beta = 1 is reached
            while true
                if verbose; fprintf('Current bridge beta = %d\n', current_bridge_point); end
                
                % Determine optimal max step size
                max_sample_step_i = determine_optimal_max_step_size(m(iTop), con, objApprox(:,:,iTop), [obj; objPrior(:,:,iTop)], optsTop(iTop), current_bridge_point, T_previous);
                if verbose; fprintf('Optimal max sample step is %d...\n', max_sample_step_i); end
                
                % Determine autocorrelation at this bridge point
                sampleOpts = optsTop(iTop);
                sampleOpts.MaxSampleStep = max_sample_step_i;
                [samples_i_bridge, step_i_bridge, max_sample_step_i] = determine_deautocorrelation_step_size_on_geometric_mix(m(iTop), con, objApprox(:,:,iTop), [obj; objPrior(:,:,iTop)], sampleOpts, current_bridge_point, T_previous, opts.TopologyTol, opts.TopologyTol);
                if verbose; fprintf('Autocorrelation dies after %d steps...\n', step_i_bridge); end
                
                % Store sampling results
                bridge_points{iTop} = [bridge_points{iTop}; current_bridge_point];
                samples{iTop} = [samples{iTop}; {samples_i_bridge}];
                sample_step_size{iTop} = [sample_step_size{iTop}; step_i_bridge];
                max_sample_step{iTop} = [max_sample_step{iTop}; max_sample_step_i];
                
                % Terminate if there is no next bridge point to find
                if current_bridge_point == 1
                    break
                end
                
                % Thin samples by autocorrelation
                samples_i_bridge = samples_i_bridge(:,1:step_i_bridge:end);
                
                % Find the next bridge point that has the appropriate distance
                % from the current bridge point
                n_samples = size(samples_i_bridge, 2);
                
                % Evaluate the approximation and the posterior at each parameter set
                logl0_evaluated = zeros(n_samples,1);
                logl1_evaluated = zeros(n_samples,1);
                for i_sample = 1:n_samples
                    % Update to parameter set of this sample
                    [mtemp, contemp] = updateAll(m(iTop), con, samples_i_bridge(:,i_sample), optsTop(iTop).UseModelSeeds, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseSeeds, optsTop(iTop).UseControls);
                    logl0_evaluated(i_sample) = ObjectiveLogLikelihood(mtemp, contemp, objApprox(:,:,iTop), optsTop(iTop));
                    
                    [mtemp, contemp] = updateAll(m(iTop), con, samples_i_bridge(:,i_sample), optsTop(iTop).UseModelSeeds, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseSeeds, optsTop(iTop).UseControls);
                    logl1_evaluated(i_sample) = ObjectiveLogLikelihood(mtemp, contemp, [obj; objPrior(:,:,iTop)], optsTop(iTop));
                end
                
                % Find proposed beta that will give the weight a variance equal
                % to a desired variance
                optimization_function = @(proposed_beta)(variance_in_weight(current_bridge_point, proposed_beta, logl0_evaluated, logl1_evaluated) - opts.VariableScheduleTarget.^2);
                
                if optimization_function(1) < 0
                    % Final window
                    % Running the optimizer here will crash
                    optimal_beta = 1;
                else
                    optimal_beta = fzero(optimization_function, [current_bridge_point, 1]);
                end
                
                % Set loop values for next iteration
                current_bridge_point = optimal_beta;
                T_previous = samples_i_bridge(:,end);
            end
        end
    else
        error('KroneckerBio:TopologyProbability:PathSchedule', [opts.PathSchedule ' is not valid for opts.PathSchedule; must be "Linear" or "Harmonic" or "Variable"'])
    end
    
    %% Computation method
    if strcmpi(opts.PathEvalMethod, 'Importance')
        % Sample until error is less than tolerance
        n_bridges = zeros(nTop,1);
        i_sample = cell(nTop,1);
        weights = cell(nTop,1);
        pmy = zeros(nTop,1);
        unc_pmy = nan(nTop,1);
        mean_weights = cell(nTop,1);
        unc_mean_weights = cell(nTop,1);
        reached_minimum_sampling = false(nTop,1);
        reached_maximum_sampling = false(nTop,1);
        
        for iTop = 1:nTop
            n_bridges(iTop) = numel(bridge_points{iTop});
            i_sample{iTop} = zeros(n_bridges(iTop)-1,1);
            weights{iTop} = cell(n_bridges(iTop)-1,1);
            mean_weights{iTop} = zeros(n_bridges(iTop)-1,1);
            unc_mean_weights{iTop}  = nan(n_bridges(iTop)-1,1);
        end
        
        while true % dowhile
            % Take additional samples
            % Look to see if any slots have not reached their minimum
            iTop = NaN;
            i_bridge = NaN;
            for jTop = 1:nTop
                for j_bridge = 1:n_bridges(jTop)-1
                    if numel(weights{jTop}{j_bridge}) < max(opts.PathMinimumSamples,2) % Must have at least two in order to do error estimate
                        iTop = jTop;
                        i_bridge = j_bridge;
                    end
                    if ~isnan(iTop); break; end
                end
                if ~isnan(iTop); break; end
            end
            
            % If the minimums are all met, then chose a topology an
            % bridge weighted by the size of the variance
            if isnan(iTop)
                % Find the topology with the largest uncertainty and the
                % bridge with the largest uncertainty and try to reduce it
                [~, iTop] = max(unc_pmy);
                [~, i_bridge] = max(unc_mean_weights{iTop});
            end
            
            % Fetch actual bridge values
            current_bridge_point = bridge_points{iTop}(i_bridge);
            next_bridge_point = bridge_points{iTop}(i_bridge + 1);
            
            % Index into sample by autocorrelation step
            i_sample{iTop}(i_bridge) = i_sample{iTop}(i_bridge) + sample_step_size{iTop}(i_bridge);
            
            % Replenish parameter sample supply when empty
            while i_sample{iTop}(i_bridge) > size(samples{iTop}{i_bridge}, 2)
                sampleOpts = optsTop(iTop);
                sampleOpts.MaxSampleStep = max_sample_step{iTop}(i_bridge);
                [more_samples, data] = sample_geometric_mix(m(iTop), con, objApprox(:,:,iTop), [obj; objPrior(:,:,iTop)], sampleOpts, current_bridge_point, samples{iTop}{i_bridge}(:,end));
                samples{iTop}{i_bridge} = [samples{iTop}{i_bridge}, more_samples];
                
                % Update sample step
                max_sample_step{iTop}(i_bridge) = data.FinalMaxSampleStep;
            end
            
            % Evaluate weight at sample and store
            [mtemp, contemp] = updateAll(m(iTop), Uzero(m(iTop)), samples{iTop}{i_bridge}(:,i_sample{iTop}(i_bridge)), optsTop(iTop).UseModelSeeds, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseSeeds, optsTop(iTop).UseControls);
            log_likelihood_0 = ObjectiveLogLikelihood(mtemp, contemp, objApprox(:,:,iTop), optsTop(iTop));
            
            [mtemp, contemp] = updateAll(m(iTop), con, samples{iTop}{i_bridge}(:,i_sample{iTop}(i_bridge)), optsTop(iTop).UseModelSeeds, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseSeeds, optsTop(iTop).UseControls);
            log_likelihood_1 = ObjectiveLogLikelihood(mtemp, contemp, [obj; objPrior(:,:,iTop)], optsTop(iTop));
            
            % Likelihood at current bridge and next bridge
            log_likelihood_i = (1 - current_bridge_point) * log_likelihood_0 + current_bridge_point * log_likelihood_1;
            log_likelihood_ip1 = (1 - next_bridge_point) * log_likelihood_0 + next_bridge_point * log_likelihood_1;
            
            % Recenter log likelihoods before exponentiating
            mean_log_likelihood = mean([log_likelihood_i; log_likelihood_ip1]);
            log_likelihood_i = log_likelihood_i - mean_log_likelihood;
            log_likelihood_ip1 = log_likelihood_ip1 - mean_log_likelihood;
            
            % Change in weight
            weights{iTop}{i_bridge} = [weights{iTop}{i_bridge}; exp(log_likelihood_ip1) ./ exp(log_likelihood_i)];
            
            % Check for having reached maximum or minimum on all bridges
            reached_minimum_sampling = true(nTop,1);
            reached_maximum_sampling = true(nTop,1);
            for iTop = 1:nTop
                for i_bridge = 1:n_bridges(iTop)-1
                    reached_minimum_sampling(iTop) = reached_minimum_sampling(iTop) & (numel(weights{iTop}{i_bridge}) >= opts.PathMinimumSamples);
                    reached_maximum_sampling(iTop) = reached_maximum_sampling(iTop) & (numel(weights{iTop}{i_bridge}) >= opts.PathMaximumSamples);
                end
            end
            
            % Compute mean and uncertainty weight at each bridge point
            for iTop = 1:nTop
                for i_bridge = 1:numel(mean_weights{iTop})
                    mean_weights{iTop}(i_bridge) = mean(weights{iTop}{i_bridge});
                    unc_mean_weights{iTop}(i_bridge) = sqrt(sum((weights{iTop}{i_bridge} - mean_weights{iTop}(i_bridge)).^2) ./ (numel(weights{iTop}{i_bridge}) - 1));
                end
            end
            
            % Compute total weight
            log_total_weights = zeros(nTop,1);
            unc_log_total_weights = zeros(nTop,1);
            for iTop = 1:nTop
                log_total_weights(iTop) = sum(log(mean_weights{iTop}));
                unc_log_total_weights(iTop) = sqrt(sum((unc_mean_weights{iTop} ./ mean_weights{iTop}).^2));
            end
            
            % Convert to topology probability
            % Include topology prior
            log_pmy = log_total_weights + log(opts.PriorTopology);
            
            % Rescale log(pmy) before returning to real space
            log_pmy = log_pmy - max(log_pmy);
            
            % Return to regular space
            pmy = exp(log_pmy);
            
            % Normalize
            pmy = pmy ./ sum(pmy);
            
            % Final error
            unc_pmy = pmy .* sqrt(unc_log_total_weights.^2 .* (1 - 2 .* pmy) + sum(pmy.^2 .* unc_log_total_weights.^2));
            
            if verbose; fprintf('Current estimate:\n'); fprintf('%d +/- %d : %d +/- %d\n', [row(pmy); row(unc_pmy); row(log_total_weights + log(opts.PriorTopology)); row(unc_log_total_weights)]); end
            
            % Terminate is minimum is met and either maximum or tolerance is met
            if all(reached_minimum_sampling) && all(reached_maximum_sampling | unc_pmy < opts.TopologyTol)
                break
            end
        end
        
        % Store relevent data
        data.Bridges = bridge_points;
        data.AllSamples = samples;
        data.SampleStepSize = sample_step_size;
        data.MaxSampleStep = max_sample_step;
        data.AllWeights = weights;
        data.MeanWeights = mean_weights;
        data.MeanWeightsUnc = unc_mean_weights;
        data.LogTotalWeights = log_total_weights;
        data.LogTotalWeightsUnc = unc_log_total_weights;
        data.Posterior = pmy;
        data.PosteriorUnc = unc_pmy;
        
    elseif strcmpi(opts.PathEvalMethod, 'Integration')
        % Sample until error is less than tolerance
        n_bridges = zeros(nTop,1);
        i_sample = cell(nTop,1);
        log_likelihoods = cell(nTop,1);
        pmy = zeros(nTop,1);
        unc_pmy = nan(nTop,1);
        
        for iTop = 1:nTop
            n_bridges(iTop) = numel(bridge_points{iTop});
            i_sample{iTop} = zeros(n_bridges(iTop),1);
            log_likelihoods{iTop} = cell(n_bridges(iTop),1);
        end
        
        while numel(log_likelihoods{1}{1}) < opts.MinimumSamples || (numel(log_likelihoods{1}{1}) < opts.MinimumSamples && ~all(err_pmy < opts.TopologyTol))
            % Sample each uncertainty topology
            for iTop = 1:nTop
                % Sample each bridge point
                for i_bridge = 1:n_bridges(iTop)
                    % Index into sample by autocorrelation step
                    i_sample{iTop}(i_bridge) = i_sample{iTop}(i_bridge) + sample_step_size{iTop}(i_bridge);
                    
                    % Replenish parameter sample supply when empty
                    while i_sample{iTop}(i_bridge) > size(samples{iTop}{i_bridge}, 2)
                        current_bridge_point = bridge_points{iTop}(i_bridge);

                        more_samples = sample_geometric_mix(m(iTop), con, objApprox(:,:,iTop), [obj; objPrior(:,:,iTop)], optsTop(iTop), current_bridge_point, samples{iTop}{i_bridge}(:,end));
                        
                        samples{iTop}{i_bridge} = [samples{iTop}{i_bridge}, more_samples];
                    end
                    
                    % Evaluate log-likelihood at sample and store
                    [mtemp, contemp] = updateAll(m(iTop), Uzero(m(iTop)), samples{iTop}{i_bridge}(:,i_sample{iTop}(i_bridge)), optsTop(iTop).UseModelSeeds, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseSeeds, optsTop(iTop).UseControls);
                    log_likelihood_0 = ObjectiveLogLikelihood(mtemp, contemp, objApprox(:,:,iTop), optsTop(iTop));
                    
                    [mtemp, contemp] = updateAll(m(iTop), con, samples{iTop}{i_bridge}(:,i_sample{iTop}(i_bridge)), optsTop(iTop).UseModelSeeds, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseSeeds, optsTop(iTop).UseControls);
                    log_likelihood_1 = ObjectiveLogLikelihood(mtemp, contemp, [obj; objPrior(:,:,iTop)], optsTop(iTop));
                    
                    log_likelihoods{iTop}{i_bridge} = [log_likelihoods{iTop}{i_bridge}; log_likelihood_1 - log_likelihood_0];
                end
            end
            
            % Mean and uncertainty at each bridge
            mean_likelihood = cell(nTop,1);
            unc_likelihood  = cell(nTop,1);
            for iTop = 1:nTop
                mean_likelihood{iTop} = zeros(n_bridges(iTop),1);
                unc_likelihood{iTop} = zeros(n_bridges(iTop),1);
                
                for i_bridge = 1:n_bridges(iTop)
                    mean_likelihood{iTop}(i_bridge) = mean(log_likelihoods{iTop}{i_bridge});
                    unc_likelihood{iTop}(i_bridge) = sqrt(sum((log_likelihoods{iTop}{i_bridge} - mean_likelihood{iTop}(i_bridge)).^2) / (numel(log_likelihoods{iTop}{i_bridge}) - 1));
                end
            end
            
            % Log marginal likelihood is integral over bridges
            % Use trapezoidal rule
            log_pym = zeros(nTop,1);
            unc_log_pym = nan(nTop,1);
            for iTop = 1:nTop
                variance_sum = 0;
                for i_bridge = 1:n_bridges(iTop)-1
                    log_pym(iTop) = log_pym(iTop) + (bridge_points{iTop}(i_bridge+1) - bridge_points{iTop}(i_bridge)) * (mean_likelihood{iTop}(i_bridge+1) + mean_likelihood{iTop}(i_bridge)) / 2;
                    variance_sum = variance_sum + (bridge_points{iTop}(i_bridge+1) - bridge_points{iTop}(i_bridge)) .^ 2 * (unc_likelihood{iTop}(i_bridge+1).^2 + unc_likelihood{iTop}(i_bridge).^2) / 4;
                end
                
                unc_log_pym(iTop) = sqrt(variance_sum);
            end
            
            % Include topology prior
            log_pmy = log_pym + log(opts.PriorTopology);
            
            % Rescale log(pmy) before returning to real space
            log_pmy = log_pmy - max(log_pmy);
            
            % Return to regular space
            pmy = exp(log_pmy);
            
            % Normalize
            pmy = pmy ./ sum(pmy);
            
            % Final error
            unc_pmy = sqrt(pmy.^2 .* (unc_log_pym.^2 .* (1 - 2 .* pmy) + sum(pmy.^2 .* unc_log_pym.^2)));
            
            if verbose; fprintf('Current estimate:\n'); fprintf('%d +/- %d : %d +/- %d\n', [row(pmy); row(unc_pmy); row(log_pym + log(opts.PriorTopology)); row(unc_log_pym)]); end
        end
    else
        error()
    end
else
    error('KroneckerBio:TopologyProbability:InvalidMethod', 'That is not a valid method for computing topology probability.')
end
end

function obj = objectiveRescaled(objOld, scale)
    obj.Type = objOld.Type;
    obj.Name = [objOld.Name ' rescaled'];
    obj.Continuous = objOld.Continuous;
    obj.Complex = objOld.Complex;
    obj.DiscreteTimes = objOld.DiscreteTimes;
    obj.g = @(t,x,u)(objOld.g(t,x,u) * scale);
    obj.dgdx = @(t,x,u)(objOld.dgdx(t,x,u) * scale);
    obj.dgdk = @(t,x,u)(objOld.dgdk(t,x,u) * scale);
    obj.d2gdx2 = @(t,x,u)(objOld.d2gdx2(t,x,u) * scale);
    obj.d2gdk2 = @(t,x,u)(objOld.d2gdk2(t,x,u) * scale);
    obj.d2gdxdk = @(t,x,u)(objOld.d2gdxdk(t,x,u) * scale);
    obj.d2gdkdx = @(t,x,u)(objOld.d2gdkdx(t,x,u) * scale);
    obj.G = @G;
    obj.dGdx = @(t,sol)(objOld.dGdx(t,sol) * scale);
    obj.dGdk = @(t,sol)(objOld.dGdk(t,sol) * scale);
    obj.dGdx2 = @(t,sol)(objOld.dGdx2(t,sol) * scale);
    obj.dGdk2 = @(t,sol)(objOld.dGdk2(t,sol) * scale);
    obj.dGdxdk = @(t,sol)(objOld.dGdxdk(t,sol) * scale);
    obj.dGdkdx = @(t,sol)(objOld.dGdkdx(t,sol) * scale);
    obj.p = @(sol)(objOld.p(sol) ^ scale);
    obj.logp = @(sol)(objOld.logp(sol) * scale);
    obj.F = @(sol)(objOld.F(sol) * scale);
    obj.Fn = @(sol)(objOld.Fn(sol) * scale);
    obj.AddData = @(sol)(objectiveRescaled(objOld.AddData(sol), scale));
    obj.Update = @update;
    
    obj = pastestruct(objectiveZero, obj);
    
    function [val, discrete_times] = G(sol)
        [val, discrete_times] = objOld.G(sol);
        val = val * scale;
    end
    
    function objNew = update(m, con, UseParams, UseSeeds, UseControls)
        objNew = objectiveRescaled(objOld.Update(m, con, UseParams, UseSeeds, UseControls), scale);
    end
end

% Function that computes the variance in the weight
function var_log_weight = variance_in_weight(current_beta, proposed_beta, log_l0_evaluated, log_l1_evaluated)
    n = numel(log_l0_evaluated);

    log_l_current_beta = (1 - current_beta) .* log_l0_evaluated - current_beta .* log_l1_evaluated;

    log_l_proposed_beta = (1 - proposed_beta) .* log_l0_evaluated - proposed_beta .* log_l1_evaluated;
    
    log_weight = log_l_current_beta - log_l_proposed_beta;
    
    mean_log_weight = mean(log_weight);
    
    var_log_weight = 1/(n-1) * sum((log_weight - mean_log_weight).^2);
end

function [samples, deautocorrelation_step_size, max_sample_step] = determine_deautocorrelation_step_size_on_geometric_mix(m, con, obj0, obj1, opts, bridge_point, T_start, autocorrelation_tolerance, minimum_autocorrelation_sample_size)
nT = numel(T_start);

% Rescale objective functions
for iObj = 1:numel(obj0)
    obj0(iObj) = objectiveRescaled(obj0(iObj), 1 - bridge_point);
end

for iObj = 1:numel(obj1)
    obj1(iObj) = objectiveRescaled(obj1(iObj), bridge_point);
end

% Initialize sampler
samples = zeros(nT,0);
deautocorrelation_step_size = NaN;

% Sample until autocorrelation is computed
while true %dowhile
    % Continue sampling from last MCMC parameter set
    [mtemp, contemp] = updateAll(m, con, T_start, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams, opts.UseSeeds, opts.UseControls);
    [more_samples, data] = SampleParameterSpace(mtemp, contemp, [obj0; obj1], opts.StepsPerCheck, opts);
    samples = [samples, more_samples];
    
    % Update starting parameter set
    T_start = samples(:,end);
    
    % Update max sample step
    max_sample_step = data.FinalMaxSampleStep;
    opts.MaxSampleStep = max_sample_step;
    
    % Compute autocorrelation for each parameter
    first_deautocorrelation = nan(nT,1);
    for iT = 1:nT
        % Find autocorrelation in this parameter
        % Use tolerance to control how deep to compute the autocorrelation
        if opts.Normalized
            acf = autocorrelation(log(samples(iT,:)));
        else
            acf = autocorrelation(samples(iT,:));
        end
        
        % Find first non-autocorrelation point
        found = find(acf < autocorrelation_tolerance, 1, 'first');
        if isempty(found)
            % No non-autocorrelation found
            first_deautocorrelation(iT) = NaN;
        else
            first_deautocorrelation(iT) = found;
        end
    end
    
    % Autocorrelation is largest autocorrelation
    % Correct for Matlab max([]) flaw by putting a 1 on the end
    % Correct for Matlab max(nan) flaw by setting to NaN if any are NaN
    deautocorrelation_step_size = max([first_deautocorrelation; 1]);
    if opts.Verbose
        fprintf('ACF estimate: [');
        for iT = 1:nT-1
            fprintf('%d,', first_deautocorrelation(iT))
        end
        fprintf('%d]...\n', first_deautocorrelation(nT))
    end
    
    % If all autocorrelations are reliable
    if all(first_deautocorrelation < size(samples,2) * minimum_autocorrelation_sample_size)
        break
    end
end
end

function optimal_max_sample_step = determine_optimal_max_step_size(m, con, obj0, obj1, opts, bridge_point, T_start)
opts.AdaptMaxSampleStep = true;
opts.AdaptThinning = false;
opts.SampleThinning = 1;
opts.Verbose = 0;

% Rescale objective functions
for iObj = 1:numel(obj0)
    obj0(iObj) = objectiveRescaled(obj0(iObj), 1 - bridge_point);
end

for iObj = 1:numel(obj1)
    obj1(iObj) = objectiveRescaled(obj1(iObj), bridge_point);
end

% Sample until max sample size stabilizes
while true %dowhile
    % Perform one batch of sampling according to opts.StepsPerCheck
    [mtemp, contemp] = updateAll(m, con, T_start, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams, opts.UseSeeds, opts.UseControls);
    [samples, data] = SampleParameterSpace(mtemp, contemp, [obj0; obj1], opts.StepsPerCheck, opts);
    
    % Start next round at stop point
    T_start = samples(:,end);
    
    % Extract acceptance ratio
    max_sample_step = data.FinalMaxSampleStep;
    
    if max_sample_step == opts.MaxSampleStep
        % This sample step size was appropriate, so keep it
        optimal_max_sample_step = max_sample_step;
        break
    else
        % Update MaxSampleStep and keep going
        opts.MaxSampleStep = max_sample_step;
    end
end
end

function [samples, data] = sample_geometric_mix(m, con, obj0, obj1, opts, bridge_point, T_start)
% Rescale objective functions
for iObj = 1:numel(obj0)
    obj0(iObj) = objectiveRescaled(obj0(iObj), 1 - bridge_point);
end

for iObj = 1:numel(obj1)
    obj1(iObj) = objectiveRescaled(obj1(iObj), bridge_point);
end

[mtemp, contemp] = updateAll(m, con, T_start, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams, opts.UseSeeds, opts.UseControls);
[samples, data] = SampleParameterSpace(mtemp, contemp, [obj0; obj1], opts.StepsPerCheck, opts);
end

function random_index = randomindexweighted(weights)
cumulative_weights = cumsum(weights);

weight_sum = cumulative_weights(end);

random_weight = rand() * weight_sum;

random_index = find(random_weight <= cumulative_weights, 1, 'last');
end