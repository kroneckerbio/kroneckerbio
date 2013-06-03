function [pmy pyTm T F] = TopologyProbability(m, con, obj, objPrior, opts, F)
% TopologyProbability Compute the relative probability that each member of
%   a set of topologies is true according to a set of
%   information-theory-based objective functions
%
%   [pmy pyTm T F] = TopologyProbability(m, con, obj, objPrior, opts, F)

% Clean up inputs
assert(nargin >= 4, 'KroneckerBio:TopologyProbability:TooFewInputs', 'TopologyProbability requires at least 4 input arguments')
if nargin < 6
    F = [];
    if nargin < 5
        opts = [];
    end
end

% Constants
nTop = numel(m);
nCon = size(con,1);
nObj = size(obj,1);

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;

defaultOpts.UseParams      = NaN;
defaultOpts.UseICs         = NaN;
defaultOpts.UseControls    = NaN;

defaultOpts.Normalized     = true;
defaultOpts.UseAdjoint     = true;
defaultOpts.LowerBound     = 0;
defaultOpts.UpperBound     = inf;

defaultOpts.TopologyMethod = 'linear';
defaultOpts.TopologyTol    = 0.05;
defaultOpts.PriorTopology  = zeros(nTop,1) + 1/nTop; % Uniform prior
defaultOpts.NeedFit        = true; % Fit the parameters
defaultOpts.StepsPerCheck  = 100;  % Number of steps between checks on the acceptance rate

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nx = zeros(nTop,1);
nk = zeros(nTop,1);
for iTop = 1:nTop
    nx(iTop) = m(iTop).nx;
    nk(iTop) = m(iTop).nk;
end

%% Standardize structures
if isempty(con)
    % Create a con and obj to fit with the prior
    con = Uzero([0 nTop]);
    for iTop = 1:nTop
        con(iTop) = Uzero(m(iTop));
    end
    nCon = 1;
    
    obj = Gzero([0 nCon nTop]);
end

% Experimental conditions
assert((size(con,2) == 1 && (opts.UseModelICs || all(nx == nx(1)))) || size(con,2) == nTop, 'KroneckerBio:TopologyProbability:ConSize', 'Second dimension of input "con" must be equal to numel(m) or 1 if opts.UseModelICs is false or every model has the same number of species')
if size(con,2) == 1
    con = repmat(con, 1,nTop);
end

% Objective functions
assert(size(obj,3) == 1 || size(obj,3) == nTop, 'KroneckerBio:TopologyProbability:ObjSize', 'Third dimension of input "obj" must be equal to numel(m) or 1')
if size(obj,3) == 1
    obj = repmat(obj, [1,1,nTop]);
end

% Objective priors
assert(size(objPrior,3) == nTop, 'KroneckerBio:TopologyProbability:ObjPriorSize', 'Third dimension of input "objPrior" must be equal to numel(m)')

% Match dimensions of objPrior to obj
filler = Gzero([size(objPrior,1), nCon-size(objPrior,2), nTop]);
for iTop = 1:nTop
    filler(:,:,iTop) = Gzero(m(iTop));
end
objPrior = [objPrior, filler];

%% Active Parameters
% Default UseParams is all kinetic parameters
if isnumeric(opts.UseParams) && (isempty(opts.UseParams) || any(isnan(opts.UseParams)))
    opts.UseParams = cell(nTop,1);
    for iTop = 1:nTop
        opts.UseParams{iTop} = 1:m(iTop).nk;
    end
end

% Default UseICs is no initial condition parameters
if isnumeric(opts.UseICs) && (isempty(opts.UseICs) || any(isnan(opts.UseICs)))
    opts.UseICs = cell(nTop,1);
    for iTop = 1:nTop
        opts.UseICs{iTop} = [];
    end
end

% Default UseControls is no controls
if isnumeric(opts.UseControls) && (isempty(opts.UseControls) || any(isnan(opts.UseControls)))
    opts.UseControls = cell(nCon,nTop);
    for iTop = 1:nTop
        opts.UseControls{iTop} = [];
    end
end

% Ensure UseParams is vector of logical indexes within a cell array
nTk = zeros(nTop,1);
for iTop = 1:nTop
    [opts.UseParams{iTop} nTk(iTop)] = fixUseParams(opts.UseParams{iTop}, nk(iTop));
end

% Ensure UseICs is a matrix of logical indexes within a cell array
nTx = zeros(nTop,1);
for iTop = 1:nTop
    [opts.UseICs{iTop} nTx(iTop)] = fixUseICs(opts.UseICs{iTop}, opts.UseModelICs, nx(iTop), nCon);
end

% Ensure UseControls is a cell array of logical vectors
nTq = zeros(nTop,1);
for iTop = 1:nTop
    [opts.UseControls{iTop} nTq(iTop)] = fixUseControls(opts.UseControls{iTop}, opts.UseModelInputs, nCon, m(iTop).nq, cat(1,con(:,iTop).nq));
end

nT = nTk + nTx + nTq;

%% Refresh structures
con = refreshCon(m, con);
obj = refreshObj(m, con, obj, opts.UseParams, opts.UseICs, opts.UseControls);

objPrior = refreshObj(m, con, objPrior, opts.UseParams, opts.UseICs, opts.UseControls);

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
for iTop = 1:nTop
    opts.LowerBound{iTop} = fixBounds(opts.LowerBound{iTop}, opts.UseModelICs, opts.UseModelInputs, opts.UseParams{iTop}, opts.UseICs{iTop}, opts.UseControls{iTop});
    opts.UpperBound{iTop} = fixBounds(opts.UpperBound{iTop}, opts.UseModelICs, opts.UseModelInputs, opts.UseParams{iTop}, opts.UseICs{iTop}, opts.UseControls{iTop});
end

%% Distribute options for topology specific functions
assert(~isfield(opts, 'ObjWeights'), 'ObjWeights is not supported for Topology Probability.')
optsTop = repmat(opts, nTop,1);
for iTop = 1:nTop
    optsTop(iTop).UseParams   = opts.UseParams{iTop};
    optsTop(iTop).UseICs      = opts.UseICs{iTop};
    if ~opts.UseModelICs
        optsTop(iTop).UseICs  = optsTop(iTop).UseICs(:,1:nCon);
    end
    optsTop(iTop).UseControls = opts.UseControls{iTop};
    optsTop(iTop).AbsTol      = opts.AbsTol{iTop};
    optsTop(iTop).LowerBound  = opts.LowerBound{iTop};
    optsTop(iTop).UpperBound  = opts.UpperBound{iTop};
end

%% Fit
if opts.NeedFit
    for iTop = 1:nTop
        if verbose; fprintf(['Fitting ' m(iTop).Name ' to objectives...\n']); end
        m(iTop) = FitObjective(m(iTop), con(:,iTop), [obj(:,:,iTop); objPrior(:,:,iTop)], optsTop(iTop));
    end
end

% Store parameter sets that were fit
T = cell(nTop,1);
for iTop = 1:nTop
    T{iTop} = collectActiveParameters(m(iTop), con(:,iTop), optsTop(iTop).UseModelICs, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseICs, optsTop(iTop).UseControls);
end

if strcmp(opts.TopologyMethod, 'linear')
    %% Information
    % Compute the information if not provided
    if isempty(F)
        F = cell(nTop,1);
        for iTop = 1:nTop
            F{iTop} = ObjectiveInformation(m(iTop), con(:,iTop), [obj(:,:,iTop); objPrior(:,:,iTop)], optsTop(iTop));
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
    pym = zeros(nTop,1); % In log space
    pyTm = ones(nTop,1);
    for iTop = 1:nTop
        if verbose; fprintf(['Computing probability of ' m(iTop).Name '...\n']); end
        % Likelihood
        pyTm(iTop) = ObjectiveProbability(m(iTop), con(:,iTop), [obj(:,:,iTop); objPrior(:,:,iTop)], optsTop(iTop));
        
        % Equation for linearization about maximum a posteriori
        %pym(iTop) = pyTm(iTop) * (2*pi).^(nT(iTop)/2) * det(F{iTop}).^(-1/2); % Without log space
        pym(iTop) = log(pyTm(iTop)) + (nT(iTop)/2) * log(2*pi) + -1/2 * sum(log(lambda{iTop}));
    end
    
    %% Compute p_m|y for the set
    % Apply topology prior
    %pmy = pym .* opts.PriorTopology; % Without log space
    pmy = pym + log(opts.PriorTopology); % In log space
    
    % Rescale p_y|m
    % Since the entire distribution is normalized at the end, this operation
    % has no mathematical effect on the answer. However, it ensures that all
    % the probabilities are representable by float64 when numbers are
    % exponentiated back into regular space.
    pmy = pmy - max(pmy);
    
    % Return to regular space
    pmy = exp(pmy);
    
    % Normalize
    pmy = pmy ./ sum(pmy);
elseif strcmp(opts.TopologyMethod, 'harmonic')
    %% Sample
    % Start from the best fit parameter set
    nSam = opts.StepsPerCheck;
    Ts = cell(nTop,1);
    pym = zeros(nTop,1);
    pyTm = cell(nTop,1);
    pyTminv = cell(nTop,1);
    needsSampled = true(nTop,1); % Can this topology be sampled
    
    step = ones(nTop,1);
    error = zeros(nTop,1);
    
    for iTop = 1:nTop
        Ts{iTop} = collectActiveParameters(m(iTop), con(:,iTop), optsTop(iTop).UseModelICs, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseICs, optsTop(iTop).UseControls);
        pyTm{iTop} = ObjectiveProbability(m(iTop), con(:,iTop), obj(:,:,iTop), optsTop(iTop)); % p_y|mT
        pyTminv{iTop} = pyTm{iTop}.^-1;
        needsSampled(iTop) = logical(pyTm{iTop} * ObjectiveProbability(m(iTop), Uzero(m(iTop)), objPrior(iTop), optsTop(iTop))) & (nT(iTop) >= 1); % p_T|y
    end
    
    while true % dowhile
        % Draw parameter sets for each topology
        for iTop = find(needsSampled)'
            % Sample the posterior parameter space
            mtemp = updateAll(m(iTop), con(:,iTop), [], Ts{iTop}(:,end), optsTop(iTop).UseModelICs, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseICs, optsTop(iTop).UseControls);
            Tsi = SampleParameterSpace(mtemp, con(:,iTop), [obj(:,:,iTop); objPrior(:,:,iTop)], nSam, optsTop(iTop));
            Ts{iTop} = [Ts{iTop}, Tsi];
            
            % Compute the likelihood
            pyTmi = zeros(nSam,1);
            for iSam = 1:nSam
                mtemp = updateAll(m(iTop), con(:,iTop), [], Tsi(:,iSam), optsTop(iTop).UseModelICs, optsTop(iTop).UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseICs, optsTop(iTop).UseControls);
                pyTmi(iSam) = ObjectiveProbability(mtemp, con(:,iTop), obj(:,:,iTop), optsTop(iTop));
            end
            pyTm{iTop} = [pyTm{iTop}; pyTmi];
            pyTminv{iTop} = [pyTminv{iTop}; pyTmi.^-1];
            
            % Compute autocorrelation in data
            first = zeros(nT(iTop),1);
            for iT = 1:nT(iTop)
                acf = autocorr(Ts{iTop}(iT,:), size(Ts{iTop},2)-1);
                
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
                error(iTop) = sqrt(pym(iTop) .^ 4 * var(pyTminv{iTop}(1:step:end),1) / (numel(1:step:numel(pyTminv{iTop})) - 1));
            else
                error(iTop) = 0;
            end
        end
        
        % Error relative to largest probability
        relerr = error ./ max(pym.* opts.PriorTopology);

        % Stop sampling those with low error
        needsSampled(~isnan(relerr) & relerr < opts.TopologyTol) = false;
        
        % Stop if uncertainty is low enough
        if all(~needsSampled)
            break
        end
    end
    
    % Include topology prior
    pmy = pym .* opts.TopologyProbability;
    
    % Normalize
    pmy = pmy ./ sum(pmy);
elseif strcmp(opts.TopologyMethod, 'annealed')
    
end