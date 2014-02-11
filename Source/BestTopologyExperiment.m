function [best data] = BestTopologyExperiment(m, con, obj, objPriorParams, objPriorSeeds, objPriorControls, conPos, objPos, target, opts)
%[best data] = BestTopologyExperiment(m, con, obj, objPrior, conPos, objPos, target, opts, pmyStart)

% Clean-up inputs
if nargin < 8
    opts = [];
end

% Constants
nTop = numel(m);
nCon = numel(con);
nObj = size(obj,1);
nConPos = numel(conPos);
nObjPos = size(objPos,1);

% Deafult options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelSeeds  = false;
defaultOpts.UseModelInputs = false;

defaultOpts.UseParams      = cell(nTop,1); 
for i=1:nTop; defaultOpts.UseParams{i} = 1:m(i).nk;end
defaultOpts.UseSeeds       = [];
defaultOpts.UseControls    = [];

defaultOpts.Normalized     = true;
defaultOpts.UseAdjoint     = true;
defaultOpts.LowerBound     = 0;
defaultOpts.UpperBound     = inf;

% Topology probability
defaultOpts.PriorTopology  = zeros(nTop,1) + 1/nTop; % Uniform prior
defaultOpts.TopologyMethod = 'Linear'; % Harmonic, Path
defaultOpts.NeedFit        = true; % Fit the parameters, only applies to starting pmy
defaultOpts.TopologyTol    = 0.05;

% Sampling
defaultOpts.StepsPerCheck  = 100;

% Best topology experiment
defaultOpts.BestTopologyMethod = 'hastings';
defaultOpts.TargetTol = 0.05;
defaultOpts.MinTargetIter = 0;
defaultOpts.MaxTargetIter = 200;

opts = mergestruct(defaultOpts, opts);
verbose = logical(max(opts.Verbose, 0));
opts.Verbose = max(opts.Verbose-1, 0);

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
    obj = Gzero;
    nCon = 1;
    nObj = 1;
end

% Experimental conditions
con = vec(con);

% Objective functions
assert(size(obj,2) == nCon, 'KroneckerBio:TopologyProbability:ObjSize', 'Second dimension of "obj" must be equal to nCon')

% Possible experimental conditions
conPos = vec(conPos);

% Possible objective functions
assert(size(objPos,2) == nConPos, 'KroneckerBio:TopologyProbability:ObjPosSize', 'Second dimension of "objPos" must be equal to nCon')

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
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, opts.UseModelSeeds, ns, nCon+nConPos);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseControls, nTq] = fixUseControls(opts.UseControls, opts.UseModelInputs, nCon+nConPos, m(1).nq, cat(1,con.nq,conPos.nq));

nT = nTk + nTs + nTq;

%% Priors
if nTk > 0
    assert(numel(objPriorParams) == nTop, ...
        'KroneckerBio:TopologyProbability:ObjPriorParamsSize', 'Input "objPriorParams" must be a vector of length numel(m)')
    objPriorParams = [reshape(objPriorParams, [1,1,nTop]), Gzero([1, nCon-1, nTop])];
else
    objPriorParams = Gzero([0,nCon,nTop]);
end

if nTs > 0
    if opts.UseModelSeeds
        assert(numel(objPriorSeeds) == 1, ...
            'KroneckerBio:TopologyProbability:ObjPriorSeedsSize', 'Input "objPriorSeeds" must be a scalar if UseModelSeeds is true')
        objPriorSeeds = [repmat(objPriorSeeds, [1,1,nTop]), Gzero([1, nCon-1, nTop])];
    else
        assert(numel(objPriorSeeds) == nCon, ...
            'KroneckerBio:TopologyProbability:ObjPriorSeedsSize', 'Input "objPriorSeeds" must be a vector of length numel(con) if UseModelSeeds is false')
        objPriorSeeds = repmat(reshape(objPriorSeeds, [1,nCon,1]), [1,1,nTop]);
    end
else
    objPriorSeeds = Gzero([0,nCon,nTop]);
end

if nTq > 0
    if opts.UseModelControls
        assert(numel(objPriorSeeds) == 1, ...
            'KroneckerBio:TopologyProbability:ObjPriorControlsSize', 'Input "objPriorControls" must be a scalar if UseModelControls is true')
        objPriorSeeds = [repmat(objPriorSeeds, [1,1,nTop]), Gzero([1, nCon-1, nTop])];
    else
        assert(numel(objPriorSeeds) == nCon, ...
            'KroneckerBio:TopologyProbability:ObjPriorControlsSize', 'Input "objPriorControls" must be a vector of length numel(con) if UseModelControls is false')
        objPriorSeeds = repmat(reshape(objPriorSeeds, [1,nCon,1]), [1,1,nTop]);
    end
else
    objPriorControls = Gzero([0,nCon,nTop]);
end

% Match dimensions of objPrior to obj
objPrior = [objPriorParams; objPriorSeeds; objPriorControls];
nObjPrior = size(objPrior,1);

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
opts.AbsTol = cell(nTop,1);
for iTop = 1:nTop
    opts.AbsTol{iTop} = emptystruct(nCon+nConPos, 'System', 'Objective', 'Sensitivity', 'Gradient', 'Adjoint');
    temp = fixAbsTol(tempAbsTol{iTop}, 1, false(nCon+nConPos,1), nx(iTop), nCon+nConPos);
    [opts.AbsTol{iTop}.System] = deal(temp{:});
    temp = fixAbsTol(tempAbsTol{iTop}, 1, true(nCon+nConPos,1), nx(iTop), nCon+nConPos);
    [opts.AbsTol{iTop}.Objective] = deal(temp{:});
    temp = fixAbsTol(tempAbsTol{iTop}, 2, false(nCon+nConPos,1), nx(iTop), nCon+nConPos, false, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams{iTop}, opts.UseSeeds, opts.UseControls);
    [opts.AbsTol{iTop}.Sensitivity] = deal(temp{:});
    temp = fixAbsTol(tempAbsTol{iTop}, 2, true(nCon+nConPos,1), nx(iTop), nCon+nConPos, false, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams{iTop}, opts.UseSeeds, opts.UseControls);
    [opts.AbsTol{iTop}.Gradient] = deal(temp{:});
    temp = fixAbsTol(tempAbsTol{iTop}, 2, false(nCon+nConPos,1), nx(iTop), nCon+nConPos, true, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams{iTop}, opts.UseSeeds, opts.UseControls);
    [opts.AbsTol{iTop}.Adjoint] = deal(temp{:});
end

%% Fix bounds
if isnumeric(opts.LowerBound)
    opts.LowerBound = repmat({opts.LowerBound}, nTop,1);
end
if isnumeric(opts.UpperBound)
    opts.UpperBound = repmat({opts.UpperBound}, nTop,1);
end
for iTop = 1:nTop
    opts.LowerBound{iTop} = fixBounds(opts.LowerBound{iTop}, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams{iTop}, opts.UseSeeds, opts.UseControls);
    opts.UpperBound{iTop} = fixBounds(opts.UpperBound{iTop}, opts.UseModelSeeds, opts.UseModelInputs, opts.UseParams{iTop}, opts.UseSeeds, opts.UseControls);
end

%% Distribute options for topology specific functions
assert(~isfield(opts, 'ObjWeights') || all(vec(opts.ObjWeights) == 1), 'ObjWeights is not supported for BestTopologyExperiments.')
optsAll = repmat(opts, nTop,1);
for iTop = 1:nTop
    optsAll(iTop).UseParams       = opts.UseParams{iTop};
    if ~opts.UseModelSeeds
        optsAll(iTop).UseSeeds    = opts.UseSeeds(:,1:nCon);
    end
    if ~opts.UseModelInputs
        optsAll(iTop).UseControls = opts.UseControls(1:nCon);
    end
    optsAll(iTop).AbsTol          = opts.AbsTol{iTop}(1:nCon);
    optsAll(iTop).LowerBound      = opts.LowerBound{iTop};
    optsAll(iTop).UpperBound      = opts.UpperBound{iTop};
end

%% Fit
if opts.NeedFit
    for iTop = 1:nTop
        if verbose; fprintf(['Fitting ' m(iTop).Name ' to existing data...\n']); end
        m(iTop) = FitObjective(m(iTop), con, [obj; objPrior(:,:,iTop)], optsAll(iTop));
    end
    opts.NeedFit = false;
end

%% Starting p_m|y
if verbose; fprintf('Computing starting probability...\n'); end
optsTemp = opts;
optsTemp.NeedFit = false;
for iTop = 1:nTop
    optsTemp.AbsTol{iTop} = optsTemp.AbsTol{iTop}(1:nCon);
end
if ~opts.UseModelSeeds
    for iTop = 1:nTop
        optsTemp.UseSeeds = optsTemp.UseSeeds(:,1:nCon);
    end
end
if ~opts.UseModelInputs
    for iTop = 1:nTop
        optsTemp.UseControls = optsTemp.UseControls(1:nCon);
    end
end

pmyStart = TopologyProbability(m, con, obj, objPriorParams, objPriorSeeds, objPriorControls, optsTemp);

%% Expected target function value for each possible experiment
% Initialize containers
Etarget    = zeros(nObjPos,nConPos);
trueTopAll = cell(nObjPos,nConPos);
trueTAll   = cell(nObjPos,nConPos);
trueyAll   = cell(nObjPos,nConPos);
fitsTAll   = cell(nObjPos,nConPos);
pmyAll     = cell(nObjPos,nConPos);
targetAll  = cell(nObjPos,nConPos);
errAll     = cell(nObjPos,nConPos);

% Initialize Metropolis-Hastings sampler
iSam = zeros(nTop,1);
sam = cell(nTop,1);
step = nan(nTop,1);
for iTop = 1:nTop
    sam{iTop} = collectActiveParameters(m(iTop), con, opts.UseModelSeeds, opts.UseModelInputs, optsAll(iTop).UseParams, optsAll(iTop).UseSeeds, optsAll(iTop).UseControls);
end

for iPosCon = 1:nConPos
    for iPosObj = 1:nObjPos
        if verbose; fprintf(['Monte Carlo sampling for ' objPos(iPosObj,iPosCon,1).Name '...\n']); end
        
        % Reset growing variables
        index = 0;
        trueTops = zeros(0,1);
        trueTs   = cell(0,1);
        trueys   = cell(0,1);
        fitsTs   = cell(nTop,0);
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
                [mDraw, conDraw] = updateAll(m(trueTop), con, sam{trueTop}(:,end), opts.UseModelSeeds, opts.UseModelInputs, optsAll(trueTop).UseParams, optsAll(trueTop).UseSeeds, optsAll(trueTop).UseControls);
                sam{trueTop} = [sam{trueTop}, SampleParameterSpace(mDraw, conDraw, [obj; objPrior(:,:,trueTop)], opts.StepsPerCheck, optsAll(trueTop))];
                
                % Autocorrelation for each parameter
                first = zeros(nT(trueTop),1);
                for iT = 1:nT(trueTop)
                    acf = autocorrelation(sam{trueTop}(iT,:));
                    
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
                [mDraw, conDraw] = updateAll(m(trueTop), con, sam{trueTop}(:,end), opts.UseModelSeeds, opts.UseModelInputs, optsAll(trueTop).UseParams, optsAll(trueTop).UseSeeds, optsAll(trueTop).UseControls);
                sam{trueTop} = [sam{trueTop}, SampleParameterSpace(mDraw, conDraw, [obj; objPrior(:,:,trueTop)], opts.StepsPerCheck, optsAll(trueTop))];
            end
            
            % Construct true parameterized model
            Ttrue = sam{trueTop}(:,iSam(trueTop));
            mRand = updateAll(m(trueTop), con, Ttrue, opts.UseModelSeeds, opts.UseModelInputs, optsAll(trueTop).UseParams, optsAll(trueTop).UseSeeds, optsAll(trueTop).UseControls);
            
            try
                % Generate data according to possible experiment
                if verbose; fprintf('Generating data from Monte Carlo model...\n'); end
                optsTemp = optsAll(iTop);
                if continuousPos(iPosCon) && complexPos(iPosCon)
                    optsTemp.AbsTol = opts.AbsTol{trueTop}(nCon+iPosCon).Objective;
                    sol = integrateObj(mRand, conPos(iPosCon), objPos(iObjPos,iConPos), optsTemp);
                elseif complexPos(iPosCon)
                    optsTemp.AbsTol = opts.AbsTol{trueTop}(nCon+iPosCon).System;
                    sol = integrateSys(mRand, conPos(iPosCon), optsTemp);
                elseif continuousPos(iPosCon)
                    optsTemp.AbsTol = opts.AbsTol{trueTop}(nCon+iPosCon).Objective;
                    sol = integrateObjSelect(mRand, conPos(iPosCon), objPos(iObjPos,iConPos), tGetPos{iPosCon}, optsTemp);
                else
                    optsTemp.AbsTol = opts.AbsTol{trueTop}(nCon+iPosCon).System;
                    sol = integrateSysSelect(mRand, conPos(iPosCon), tGetPos{iPosCon}, optsTemp);
                end
                
                % Create objective function appropriate to each topology
                [new_obj, new_data] = objPos(nObjPos,nConPos).AddData(sol);
                
                % Fit the new data
                mfit = m;
                for iTop = 1:nTop
                    if verbose; fprintf(['Fitting ' m(iTop).Name ' to existing data plus Monte Carlo data...\n']); end
                    optsTemp = optsAll(iTop);
                    optsTemp.AbsTol = [optsTemp.AbsTol; opts.AbsTol{iTop}(nCon+iPosCon)];
                    if ~opts.UseModelSeeds
                        optsTemp.UseSeeds = [optsTemp.UseSeeds, opts.UseSeeds(:,nCon+iPosCon)];
                    end
                    if ~opts.UseModelInputs
                        optsTemp.UseControls = [optsTemp.UseControls; opts.UseControls(nCon+iPosCon)];
                    end
                    mfit(iTop) = FitObjective(mfit(iTop), [con; conPos(iPosCon)], [[obj; objPrior(:,:,iTop)], [new_obj; Gzero([nObj+nObjPrior-1,1])]], optsTemp);
                end
                
                % Compute the topology probability
                if verbose; fprintf('Computing topology probability after adding Monte Carlo data...\n'); end
                optsTemp = opts;
                opts.NeedFit = false;
                for iTop = 1:nTop
                    optsTemp.AbsTol{iTop} = optsTemp.AbsTol{iTop}([1:nCon,iPosCon]);
                end
                if ~opts.UseModelSeeds
                    for iTop = 1:nTop
                        optsTemp.UseSeeds = optsTemp.UseSeeds(:,[1:nCon,iPosCon]);
                    end
                end
                if ~opts.UseModelInputs
                    for iTop = 1:nTop
                        optsTemp.UseControls = optsTemp.UseControls([1:nCon,iPosCon]);
                    end
                end
                pmy = TopologyProbability(mfit, [con; conPos(iPosCon)], [obj, [new_obj; Gzero([nObj-1,1])]], objPriorParams, objPriorSeeds, objPriorControls, optsTemp);
            catch me
                if strcmp(me.identifier, 'KroneckerBio:accumulateSol:IntegrationFailure')
                    % The integrator crashed unexpectedly
                    % Just ignore this draw and start over
                    if verbose; fprintf('Integrator crashed--discarding draw and resuming...\n'); end
                    continue
                else
                    me.rethrow()
                end
            end
            
            % Grow Monte Carlo variables if necessary
            index = index + 1;
            if index > numel(trueTops)
                trueTops = [trueTops; zeros(numel(trueTops),1)];
                trueTs   = [trueTs;   cell(numel(trueTs),1)];
                trueys   = [trueys;   cell(numel(trueys),1)];
                fitsTs   = [fitsTs,   cell(nTop, size(fitsTs,2))];
                pmys     = [pmys,     zeros(nTop, size(pmys,2))];
                targets  = [targets;  zeros(numel(targets),1)];
                errs     = [errs;     zeros(numel(errs),1)];
            end
            
            % Store iteration values
            trueTops(index) = trueTop;
            trueTs{index}   = Ttrue;
            trueys{index}   = new_data;
            for iTop = 1:nTop
                fitsTs{iTop,index} = collectActiveParameters(mfit(iTop), con, opts.UseModelSeeds, opts.UseModelInputs, optsAll(iTop).UseParams, optsAll(iTop).UseSeeds, optsAll(iTop).UseControls);
            end
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
                trueyAll{iPosObj,iPosCon}   = trueys(1:index);
                fitsTAll{iPosObj,iPosCon}   = fitsTs(:,1:index);
                pmyAll{iPosObj,iPosCon}     = pmys(:,1:index);
                targetAll{iPosObj,iPosCon}  = targets(1:index);
                errAll{iPosObj,iPosCon}     = errs(1:index);
                break
            end
        end % Monte Carlo
    end % iPosObj
end % iPosCon

[unused, best] = min(Etarget);

data.Targets               = Etarget;
data.TopologiesChosen      = trueTopAll;
data.ParametersChosen      = trueTAll;
data.MeasurementsChosen    = trueyAll;
data.FittedParameters      = fitsTs;
data.Startingpmy           = pmyStart;
data.Allpmy                = pmyAll;
data.AllTargets            = targetAll;
data.AllWeightedMeanErrors = errAll;
