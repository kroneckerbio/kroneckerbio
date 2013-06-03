function [best data] = BestTopologyExperiment(m, con, obj, objPrior, conPos, objPos, target, opts, pmyStart)
%[best data] = BestTopologyExperiment(m, con, obj, objPrior, conPos,
%objPos, target, opts, pmyStart)

% Clean-up inputs
if nargin < 9
    pmyStart = [];
    if nargin < 8
        opts = [];
    end
end

% Constants
nTop = numel(m);
nCon = size(con,1);
nObj = size(obj,1);
nObjPrior = size(objPrior,1);
nConPos = size(conPos,1);
nObjPos = size(objPos,1);

% Options
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

defaultOpts.PriorTopology  = zeros(nTop,1) + 1/nTop; % Uniform prior
defaultOpts.NeedFit        = true; % Fit the parameters

defaultOpts.StepsPerCheck  = 100;

defaultOpts.BestTopologyMethod = 'hastings';
defaultOpts.TargetTol = 0.05;
defaultOpts.MinTargetIter = 0;
defaultOpts.MaxTargetIter = 200;

opts = mergestruct(defaultOpts, opts);
verbose = logical(max(opts.Verbose, 0));
opts.Verbose = max(opts.Verbose-1, 0);

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
    
    obj = Gzero([1 nCon nTop]);
    nObj = 1;
end

% Experimental conditions
assert((size(con,2) == 1 && (opts.UseModelICs || all(nx == nx(1)))) || size(con,2) == nTop, 'KroneckerBio:BestTopologyExperiment:ConSize', 'Second dimension of input "con" must be equal to numel(m) or 1 if opts.UseModelICs is false or every model has the same number of species')
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

% Possible experimental conditions
assert((size(conPos,2) == 1 && (opts.UseModelICs || all(nx == nx(1)))) || size(conPos,2) == nTop, 'KroneckerBio:BestTopologyExperiment:ConPosSize', 'Second dimension of input "conPos" must be equal to numel(m) or 1 if opts.UseModelICs is false or every model has the same number of species')
if size(conPos,2) == 1
    conPos = repmat(conPos, 1,nTop);
end

% Possible objective functions
assert(size(objPos,3) == 1 || size(objPos,3) == nTop, 'KroneckerBio:TopologyProbability:PosObjSize', 'Third dimension of input "objPos" must be equal to numel(m) or 1')
if size(objPos,3) == 1
    objPos = repmat(objPos, [1,1,nTop]);
end

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
    [opts.UseICs{iTop} nTx(iTop)] = fixUseICs(opts.UseICs{iTop}, opts.UseModelICs, nx(iTop), nCon+nConPos);
end

% Ensure UseControls is a cell array of logical vectors
nTq = zeros(nTop,1);
for iTop = 1:nTop
    [opts.UseControls{iTop} nTq(iTop)] = fixUseControls(opts.UseControls{iTop}, opts.UseModelInputs, nCon+nConPos, m(iTop).nq, cat(1, con(:,iTop).nq, conPos(:,iTop).nq));
end

nT = nTk + nTx + nTq;

%% Refresh structures
con = refreshCon(m, con);
obj = refreshObj(m, con, obj, opts.UseParams, opts.UseICs, opts.UseControls);
objPrior = refreshObj(m, con, objPrior, opts.UseParams, opts.UseICs, opts.UseControls);

conPos = refreshCon(m, conPos);
objPos = refreshObj(m, conPos, objPos, opts.UseParams, opts.UseICs, opts.UseControls);

%% Integration type: simple, continuous, complex, or both
[continuous, complex, tGet] = fixIntegrationType(con, obj);
[continuousPos, complexPos, tGetPos] = fixIntegrationType(conPos, objPos);

%% Tolerances
% RelTol
opts.RelTol = fixRelTol(opts.RelTol);

% AbsTol should be a cell vector nTop
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

% Create a universal AbsTol structure
simAbsTol = cell(nCon+nConPos,nTop);
opts.AbsTol = cell(nTop,1);
for iTop = 1:nTop
    opts.AbsTol{iTop} = emptystruct(nCon+nConPos, 'System', 'Adjoint', 'Sensitivity');
    temp = fixAbsTol(tempAbsTol{iTop}, 1, false(nCon+nConPos,1), nx(iTop), nCon+nConPos);
    [opts.AbsTol{iTop}.System] = deal(temp{:});
    temp = fixAbsTol(tempAbsTol{iTop}, 2, false(nCon+nConPos,1), nx(iTop), nCon+nConPos, false, opts.UseModelICs, opts.UseModelInputs, opts.UseParams{iTop}, opts.UseICs{iTop}, opts.UseControls{iTop});
    [opts.AbsTol{iTop}.Sensitivity] = deal(temp{:});
    temp = fixAbsTol(tempAbsTol{iTop}, 2, false(nCon+nConPos,1), nx(iTop), nCon+nConPos, true, opts.UseModelICs, opts.UseModelInputs, opts.UseParams{iTop}, opts.UseICs{iTop}, opts.UseControls{iTop});
    [opts.AbsTol{iTop}.Adjoint] = deal(temp{:});
    simAbsTol(:,iTop) = fixAbsTol(tempAbsTol{iTop}, 1, [continuous; continuousPos], nx(iTop), nCon+nConPos, opts.UseAdjoint, opts.UseModelICs, opts.UseModelInputs, opts.UseParams{iTop}, opts.UseICs{iTop}, opts.UseControls{iTop});
end

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
assert(~isfield(opts, 'ObjWeights') || all(vec(opts.ObjWeights) == 1), 'ObjWeights is not supported for BestTopologyExperiments.')
optsTop = repmat(opts, nTop,1);
for iTop = 1:nTop
    optsTop(iTop).UseParams   = opts.UseParams{iTop};
    optsTop(iTop).UseICs      = opts.UseICs{iTop};
    if ~opts.UseModelICs
        optsTop(iTop).UseICs  = optsTop(iTop).UseICs(:,1:nCon);
    end
    optsTop(iTop).UseControls = opts.UseControls{iTop}(1:nCon);
    optsTop(iTop).AbsTol      = opts.AbsTol{iTop}(1:nCon);
    optsTop(iTop).LowerBound  = opts.LowerBound{iTop};
    optsTop(iTop).UpperBound  = opts.UpperBound{iTop};
end

%% Fit
if opts.NeedFit
    for iTop = 1:nTop
        if verbose; fprintf(['Fitting ' m(iTop).Name ' to existing data...\n']); end
        m(iTop) = FitObjective(m(iTop), con(:,iTop), [obj(:,:,iTop); objPrior(:,:,iTop)], optsTop(iTop));
    end
    opts.NeedFit = false;
end

%% Starting p_m|y
if isempty(pmyStart)
    if verbose; fprintf('Computing starting probability...\n'); end
    optsSub = opts;
    for iTop = 1:nTop
        optsSub.AbsTol{iTop} = optsSub.AbsTol{iTop}(1:nCon);
    end
    if ~opts.UseModelICs
        for iTop = 1:nTop
            optsSub.UseICs{iTop} = optsSub.UseICs{iTop}(:,1:nCon);
        end
    end
    if ~opts.UseModelInputs
        for iTop = 1:nTop
            optsSub.UseControls{iTop} = optsSub.UseControls{iTop}(1:nCon);
        end
    end
    
    pmyStart = TopologyProbability(m, con, obj, objPrior, optsSub);
end

%% Expected target function value for each possible experiment
if strcmp(opts.BestTopologyMethod, 'hastings')
    % Initialize containers
    Etarget    = zeros(nObjPos,nConPos);
    trueTopAll = cell(nObjPos,nConPos);
    trueTAll   = cell(nObjPos,nConPos);
    pmyAll     = cell(nObjPos,nConPos);
    targetAll  = cell(nObjPos,nConPos);
    errAll     = cell(nObjPos,nConPos);
    
    % Initialize Metropolis-Hastings sampler
    iSam = zeros(nTop,1);
    sam = cell(nTop,1);
    step = nan(nTop,1);
    for iTop = 1:nTop
        sam{iTop} = collectActiveParameters(m(iTop), con(:,iTop), opts.UseModelICs, opts.UseModelInputs, optsTop(iTop).UseParams, optsTop(iTop).UseICs, optsTop(iTop).UseControls);
    end
    
    for iPosCon = 1:nConPos
        for iPosObj = 1:nObjPos
            if verbose; fprintf(['Monte Carlo sampling for ' objPos(iPosObj,iPosCon,1).Name '...\n']); end
            
            % Reset growing variables
            index = 0;
            trueTops = zeros(0,1);
            trueTs   = cell(0,1);
            pmys     = zeros(nTop,0);
            targets  = zeros(0,1);
            errs     = zeros(0,1);
            
            while true %dowhile
                % Draw random topology
                trueTop = randp(pmyStart);
                if verbose; fprintf('Model %d (%s) was randomly chosen\n', trueTop, m(trueTop).Name); end
                
                % Sample until autocorrelation is computed
                while isnan(step(trueTop))
                    if verbose; fprintf('Sampling for autocorrelation...\n'); end
                    [mDraw, conDraw, objDraw] = updateAll(m(trueTop), con(:,trueTop), [obj(:,:,trueTop); objPrior(:,:,trueTop)], sam{trueTop}(:,end), opts.UseModelICs, opts.UseModelInputs, optsTop(trueTop).UseParams, optsTop(trueTop).UseICs, optsTop(trueTop).UseControls);
                    sam{trueTop} = [sam{trueTop}, SampleParameterSpace(mDraw, conDraw, objDraw, opts.StepsPerCheck, optsTop(trueTop))];
                    
                    % Autocorrelation for each parameter
                    first = zeros(nT(trueTop),1);
                    for iT = 1:nT(trueTop)
                        acf = autocorr(sam{trueTop}(iT,:), size(sam{trueTop},2)-1);
                        
                        % Find first non-autocorrelation point
                        found = find(acf < opts.TargetTol, 1, 'first');
                        if isempty(first)
                            % No non-autocorrelation found
                            first(iT) = nan;
                        else
                            first(iT) = found;
                        end
                    end
                    step(trueTop) = max(first);
                    
                    % Correct for Matlab max(nan) flaw
                    if any(isnan(first))
                        step(trueTop) = nan;
                    end
                    
                    if verbose; fprintf('Autocorrelation dies after %d steps according to a sample of %d...\n', step(trueTop), size(sam{trueTop},2)); end
                end
                
                % Index ahead into sample by autocorrelation step
                iSam(trueTop) = iSam(trueTop) + step(trueTop);
                
                % Replenish parameter sample supply when empty
                if iSam(trueTop) > size(sam{trueTop}, 2)
                    if verbose; fprintf('Generating additional parameter sampling...\n'); end
                    [mDraw, conDraw, objDraw] = updateAll(m(trueTop), con(:,trueTop), [obj(:,:,trueTop); objPrior(:,:,trueTop)], sam{trueTop}(:,end), opts.UseModelICs, opts.UseModelInputs, optsTop(trueTop).UseParams, optsTop(trueTop).UseICs, optsTop(trueTop).UseControls);
                    sam{trueTop} = [sam{trueTop}, SampleParameterSpace(mDraw, conDraw, objDraw, opts.StepsPerCheck, optsTop(trueTop))];
                end
                
                % Construct true parameterized model
                Ttrue = sam{trueTop}(:,iSam(trueTop));
                mRand = updateAll(m(trueTop), con(:,trueTop), [obj(:,:,trueTop); objPrior(:,:,trueTop)], Ttrue, opts.UseModelICs, opts.UseModelInputs, optsTop(trueTop).UseParams, optsTop(trueTop).UseICs, optsTop(trueTop).UseControls);
                
                % Generate data according to possible experiment
                if verbose; fprintf('Generating data from Monte Carlo model...\n'); end
                optsSub = optsTop(iTop);
                optsSub.AbsTol = simAbsTol{nCon+iPosCon,trueTop};
                if continuousPos(iPosCon) && complexPos(iPosCon)
                    sol = integrateObj(mRand, conPos(iPosCon), objPos(iObjPos,iConPos,trueTop), optsSub);
                elseif complexPos(iPosCon)
                    sol = integrateSys(mRand, conPos(iPosCon,trueTop), optsSub);
                elseif continuousPos(iPosCon)
                    sol = integrateObjSelect(mRand, conPos(iPosCon,trueTop), objPos(iObjPos,iConPos,trueTop), tGetPos{iPosCon,trueTop}, optsSub);
                else
                    sol = integrateSysSelect(mRand, conPos(iPosCon,trueTop), tGetPos{iPosCon,trueTop}, optsSub);
                end
                
                % Create objective function appropriate to each topology
                newObj = pastestruct(Gzero(mRand), objPos(nObjPos,nConPos).AddData(sol));
                
                % Fit the new data
                mfit = m;
                for iTop = 1:nTop
                    if verbose; fprintf(['Fitting ' m(iTop).Name ' to existing data plus Monte Carlo data...\n']); end
                    optsSub = optsTop(iTop);
                    if ~opts.UseModelICs
                        optsSub.UseICs = [optsSub.UseICs, opts.UseICs{iTop}(:,nCon+iPosCon)];
                    end
                    if ~opts.UseModelInputs
                        optsSub.UseControls = [optsSub.UseControls; opts.UseControls{iTop}(nCon+iPosCon)];
                    end
                    optsSub.AbsTol = [optsSub.AbsTol; opts.AbsTol{iTop}(nCon+iPosCon)];
                    mfit(iTop) = FitObjective(mfit(iTop), [con(:,iTop); conPos(iPosCon,iTop)], [[obj(:,:,iTop); objPrior(:,:,iTop)], [newObj; repmat(Gzero(m(iTop)), nObj+nObjPrior-1,1)]], optsSub);
                end
                
                % Predicted p_m|y
                if verbose; fprintf('Computing topology probability after adding Monte Carlo data...\n'); end
                optsSub = opts;
                for iTop = 1:nTop
                    optsSub.AbsTol{iTop} = optsSub.AbsTol{iTop}([1:nCon,iPosCon]);
                end
                if ~opts.UseModelICs
                    for iTop = 1:nTop
                        optsSub.UseICs{iTop} = optsSub.UseICs{iTop}(:,[1:nCon,iPosCon]);
                    end
                end
                if ~opts.UseModelInputs
                    for iTop = 1:nTop
                        optsSub.UseControls{iTop} = optsSub.UseControls{iTop}([1:nCon,iPosCon]);
                    end
                end
                pmy = TopologyProbability(mfit, [con; conPos(iPosCon,:)], [obj, [repmat(newObj, [1,1,nTop]); repmat(Gzero(m(trueTop)), [nObj-1,1,nTop])]], objPrior, optsSub);
                
                % Grow Monte Carlo variables if necessary
                index = index + 1;
                if index > numel(trueTops)
                    trueTops = [trueTops; zeros(numel(trueTops),1)];
                    trueTs   = [trueTs;   zeros(numel(trueTs),1)];
                    pmys     = [pmys,     zeros(nTop, size(pmys,2))];
                    targets  = [targets;  zeros(numel(targets),1)];
                    errs     = [errs;     zeros(numel(errs),1)];
                end
                
                % Store iteration values
                trueTops(index) = trueTop;
                trueTs{index}   = Ttrue;
                pmys(:,index)   = pmy;
                
                % Evaluate target function
                targets(index) = target(pmy);
                if verbose; fprintf('Target Value = %-6.4g\n', targets(index)); end
                
                % Compute error in expected target function value
                err = sqrt(var(targets(1:index)) ./ (index-1));
                errs(index) = err;
                if verbose; fprintf('Uncertainty = %-6.4g\n', err); end                
                
                % Termination criteria
                if ((err <= opts.TargetTol) && (index >= opts.MinTargetIter)) || (index >= opts.MaxTargetIter)
                    if verbose; fprintf('Monte Carlo terminated\n'); end
                    
                    % Store results for this experiment
                    Etarget(iPosObj,iPosCon) = mean(targets(1:index));
                    
                    trueTopAll{iPosObj,iPosCon} = trueTops(1:index);
                    trueTAll{iPosObj,iPosCon}   = trueTs(1:index);
                    pmyAll{iPosObj,iPosCon}     = pmys(:,1:index);
                    targetAll{iPosObj,iPosCon}  = targets(1:index);
                    errAll{iPosObj,iPosCon}     = errs(1:index);
                    break
                end
            end % Monte Carlo
        end % iPosObj
    end % iPosCon
elseif strcmp(opts.BestTopologyMethod, 'importance')
    % Linearized parameter variance
    V = cell(nTop,1);
    for iTop = 1:nTop
        if verbose; fprintf(['Computing posterior variance of ' m(iTop).Name '...\n']); end
        V{iTop} = ObjectiveInformation(m(iTop), con(:,iTop), [obj(:,:,iTop); objPrior(:,:,iTop)], optsTop(iTop));
        V{iTop} = posdef(V{iTop}); % Posterior information
        V{iTop} = inv(V{iTop});    % Posterior variance
    end
    
    % Monte Carlo parameters
    nCarlo = opts.MaxTargetIter; % Maximum number of Monte Carlo iterations
    pdrawmax = zeros(nTop,1); % Maximum value coming from mvnpdf on Monte Carlo parameters
    prealmax = zeros(nTop,1); % Maximum value coming from likelihood for fit to data
    
    for iTop = 1:nTop
        pdrawmax(iTop) = mvnpdf(zeros(nT(iTop),1), zeros(nT(iTop),1), V{iTop});
        prealmax(iTop) = ObjectiveProbability(m(iTop), con(:,iTop), obj(:,:,iTop), optsTop(iTop));
    end
    
    % Initialize vectors
    Etarget = zeros(nObjPos,nConPos);
    pmyAll    = cell(nObjPos,nConPos,nTop);
    targetAll = cell(nObjPos,nConPos,nTop);
    weightAll = cell(nObjPos,nConPos,nTop);
    wmeanAll  = cell(nObjPos,nConPos);
    wvarAll   = cell(nObjPos,nConPos);
    errAll    = cell(nObjPos,nConPos);
    for iPosCon = 1:nConPos
        for iPosObj = 1:nObjPos
            wmeanAll{nObjPos,nConPos}  = zeros(nCarlo,1);
            wvarAll{nObjPos,nConPos}   = zeros(nCarlo,1);
            errAll{nObjPos,nConPos}    = zeros(nCarlo,1);
            for iTop = 1:nTop
                pmyAll{nObjPos,nConPos,iTop} = zeros(nTop,nCarlo);
            end
        end
    end
    
    % Initialize Gibbs sampler
    iSam = inf(nTop,1); % Start with replenishment
    sam = cell(nTop,1);
    mu  = cell(nTop,1);
    lb  = cell(nTop,1);
    ub  = cell(nTop,1);
    for iTop = 1:nTop
        sam{iTop} = collectActiveParameters(m(iTop), con(:,iTop), optsTop(iTop).UseParams, optsTop(iTop).UseICs, opts.UseModelICs);
        mu{iTop}  = m(iTop).k(optsTop(iTop).UseParams);
        lb{iTop}  = optsTop(iTop).LowerBound;
        ub{iTop}  = optsTop(iTop).UpperBound;
        if opts.Normalized
            sam{iTop} = log(sam{iTop});
            mu{iTop}  = log(mu{iTop});
            lb{iTop}  = log(lb{iTop});
            ub{iTop}  = log(ub{iTop});
        end
    end
    
    for iPosCon = 1:nConPos
        for iPosObj = 1:nObjPos
            if verbose; fprintf(['Monte Carlo sampling for ' objPos(iPosObj,iPosCon,1).Name '...\n']); end
            
            index = 0;
            indexes = zeros(nTop,1);
            
            % Reset growing variables
            targets = cell(iTop,1);
            weights = cell(iTop,1);
            for iTop = 1:nTop
                targets{iTop} = zeros(nCarlo,1);
                weights{iTop} = zeros(nCarlo,1);
            end
            
            while true %dowhile
                % Draw random topology
                if verbose; fprintf('Drawing random model...\n'); end
                trueTop = randp(pmyStart);
                index = index + 1; % Total count for this experiment
                indexes(trueTop) = indexes(trueTop) + 1; % Count for a particular topology
                if verbose2; fprintf('Model %d (%s) was chosen\n', trueTop, m(trueTop).Name); end
                
                % Draw parameter set
                while true %dowhile
                    % Replenish sample supply when empty
                    if iSam(trueTop) > size(sam{trueTop}, 2)
                        sam{trueTop} = mvnbndrndgibbs(mu{trueTop}, V{trueTop}, lb{trueTop}, ub{trueTop}, sam{trueTop}(:,end), nCarlo, 0, 0);
                        sam{trueTop} = sam{trueTop}(:,randperm(nCarlo)); % Scramble to reduce serial correlation
                        iSam(trueTop) = 1;
                    end
                    
                    % Fetch random parameter set to try
                    Trand = sam{trueTop}(:,iSam(trueTop));
                    chi2Trand = (mu{trueTop} - Trand).' * (V{trueTop} \ (mu{trueTop} - Trand)); % Needed to test if this draw is safe
                    if opts.Normalized
                        Trand = exp(Trand);
                    end
                    iSam(trueTop) = iSam(trueTop) + 1;
                    
                    % Probability of drawing this set
                    if opts.Normalized
                        pdraw = mvnpdf(log(Trand), log(m(trueTop).k(optsTop(trueTop).UseParams)), V{trueTop}) / pdrawmax(trueTop); % Chance of drawing in monte carlo
                    else
                        pdraw = mvnpdf(Trand, m(trueTop).k(optsTop(trueTop).UseParams), V{trueTop}) / pdrawmax(trueTop); % Chance of drawing in monte carlo
                    end
                    if verbose2; fprintf('pdraw = %-6.4g\n', pdraw); end
                    
                    if chi2pvalue(chi2Trand, nT(trueTop)) > opts.TargetTol % pdraw should be safe
                        % Update with new parameter set
                        [mRand conRand objRand] = updateAll(m(trueTop), con(:,trueTop), obj(:,:,trueTop), Trand, opts.UseModelICs, opts.UseModelInputs, optsTop(trueTop).UseParams, optsTop(trueTop).UseICs, optsTop(trueTop).UseControls);
                        
                        % Probability that this set is real
                        preal = ObjectiveProbability(mRand, conRand, objRand, optsTop(trueTop)) / prealmax(trueTop);
                        if verbose2; fprintf('preal = %-6.4g\n', preal); end
                        
                        % Weight of this draw
                        weight = preal / pdraw;
                        
                        if weight > opts.TargetTol % preal should be meaningful
                            % Contributable parameter set found
                            break
                        end
                    end
                end % Gibbs
                
                % Generate data
                if verbose; fprintf('Generating data from Monte Carlo model...\n'); end
                optsSub = opts;
                optsSub.AbsTol = opts.AbsTol{nCon+iPosCon,trueTop};
                optsSub.tGet   = postGet{iPosCon};
                if continuousPos(iPosCon) && complexPos(iPosCon)
                    sol = integrateObj(mRand, conPos(iPosCon), objPos(nObjPos,nConPos), optsSub);
                elseif complexPos(iPosCon)
                    sol = integrateSys(mRand, conPos(iPosCon,trueTop), optsSub);
                elseif continuousPos(iPosCon)
                    sol = integrateObjDisc(mRand, conPos(iPosCon,trueTop), objPos(nObjPos,nConPos), optsSub);
                else
                    sol = integrateSysDisc(mRand, conPos(iPosCon,trueTop), optsSub);
                end
                
                % Create objective function appropriate to each topology
                sol.c = mRand.c;
                newObj = pastestruct(Gzero(mRand), objPos(nObjPos,nConPos).AddData(sol, conPos(iPosCon,iTop).u));
                
                % Fit the new data
                mfit = m;
                for iTop = 1:nTop
                    if verbose; fprintf(['Fitting ' m(iTop).Name ' to existing data plus Monte Carlo data...\n']); end
                    optsSub = optsTop(iTop);
                    if ~opts.UseModelICs
                        optsSub.UseICs = [optsSub.UseICs, opts.UseICs{iTop}(:,nCon+iPosCon)];
                    end
                    optsSub.AbsTol = [optsSub.AbsTol; opts.AbsTol(nCon+iPosCon,iTop)]; % Final experiment has no tolerance needs
                    mfit(iTop) = FitObjective(mfit(iTop), [con(:,iTop); conPos(iPosCon,iTop)], [[obj(:,:,iTop); objPrior(:,:,iTop)], [newObj; repmat(Gzero(m(iTop)), nObj-1,1)]], optsSub);
                end
                
                % Predicted p_m|y
                if verbose; fprintf('Computing topology probability after adding Monte Carlo data...\n'); end
                optsSub = opts;
                optsSub.AbsTol = optsSub.AbsTol([1:nCon,iPosCon],:);
                if ~opts.UseModelICs
                    for iTop = 1:nTop
                        optsSub.UseICs{iTop} = optsSub.UseICs{iTop}(:,[1:nCon,iPosCon]);
                    end
                end
                pmy = TopologyProbability(mfit, [con; conPos(iPosCon,:)], [[obj; objPrior], [repmat(newObj, [1,1,nTop]), repmat(Gzero(m(1)), [1,nObj-1,nTop])]], optsSub);
                
                % Evaluate target function
                targets{trueTop}(indexes(trueTop)) = target(pmy);
                if verbose; fprintf('Target Value = %-6.4g\n', targets{trueTop}(indexes(trueTop))); end
                
                % Importance sampling weight
                weights{trueTop}(indexes(trueTop)) = weight; % Keep all the weights
                
                % Reweight the importance weights to balance the topologies
                targetStacked     = zeros(index,1); % Stack targets into a vector
                reweightedWeights = zeros(index,1);
                rwwInd = 0;
                for iTop = 1:nTop
                    % Compute mean topology weight
                    topoWeight = mean(weights{iTop}(1:indexes(iTop)));
                    if isnan(topoWeight)
                        % This topology has no draws yet, ignore it
                        continue
                    end
                    
                    % Reweight accordingly
                    nextrwwInd = rwwInd+indexes(iTop);
                    targetStacked(rwwInd+1:nextrwwInd) = targets{iTop}(1:indexes(iTop));
                    reweightedWeights(rwwInd+1:nextrwwInd) = weights{iTop}(1:indexes(iTop)) ./ topoWeight;
                    rwwInd = nextrwwInd;
                end
                
                wmean  = weightedmean(targetStacked, reweightedWeights, 1);
                wvar   = weightedvar(targetStacked, reweightedWeights, 1);
                err    = sqrt(wvar * sum((reweightedWeights./sum(reweightedWeights)).^2));
                if verbose; fprintf('Mean = %-6.4g Error = %-6.4g\n', wmean, err); end
                
                % Store iteration results
                pmyAll{iPosObj,iPosCon,trueTop}(:,indexes(trueTop)) = pmy;
                wmeanAll{iPosObj,iPosCon}(index) = wmean;
                wvarAll{iPosObj,iPosCon}(index)  = wvar;
                errAll{iPosObj,iPosCon}(index)   = err;
                
                if ((err <= opts.TargetTol) && (index >= opts.MinTargetIter)) || (index >= opts.MaxTargetIter)
                    if verbose; fprintf('Monte Carlo terminated\n'); end
                    Etarget(iPosCon,iPosObj) = wmeanAll{iPosCon,iPosObj}(index);
                    for iTop = 1:nTop
                        pmyAll{iPosObj,iPosCon,iTop} = pmyAll{iPosObj,iPosCon,iTop}(:,1:indexes(iTop));
                        targetAll{iPosObj,iPosCon,iTop} = targets{iTop}(1:indexes(iTop));
                        weightAll{iPosObj,iPosCon,iTop} = weights{iTop}(1:indexes(iTop));
                    end
                    wmeanAll{iPosObj,iPosCon}  = wmeanAll{iPosObj,iPosCon}(1:index);
                    wvarAll{iPosObj,iPosCon}   = wvarAll{iPosObj,iPosCon}(1:index);
                    errAll{iPosObj,iPosCon}    = errAll{iPosObj,iPosCon}(1:index);
                    break
                end
            end % Monte Carlo
        end % iPosObj
    end % iPosCon
end

[unused best] = min(Etarget);

if nargout >= 2
    data.Targets               = Etarget;
    data.TopologiesChosen      = trueTopAll;
    data.ParametersChosen      = trueTAll;
    data.Startingpmy           = pmyStart;
    data.Allpmy                = pmyAll;
    data.AllTargets            = targetAll;
    data.AllWeightedMeanErrors = errAll;
end
