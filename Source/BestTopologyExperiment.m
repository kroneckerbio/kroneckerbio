function [best, data] = BestTopologyExperiment(m, con, obj, objPriorParams, objPriorSeeds, objPriorInputControls, objPriorDoseControls, con_pos, obs_pos, target, opts)
%BestTopologyExperiment Determine which experiments will most efficiently
%   minimize a target function of the posterior probability distribution of
%   the topologies. Traditionally, this algorithm is used to minimize some
%   uncertainty in the topology.
%
%   best = BestTopologyExperiment(m, con, obj, obj_prior_params, 
%               obj_prior_seeds, obj_prior_input_controls, 
%               obj_prior_dose_controls, con_pos, obs_pos, target, opts)
%
%   This algorithm generates random models according to the current
%   probability distribution according to con and obj and priors. Each
%   random model is simulated according to each con_pos and observed
%   according to each obs_pos. Then the resulting objective function is
%   used to recompute the posterior probability. The target function is
%   used to evalute the quality of that posterior probability. (A sane
%   target function returns a smaller scalar for distributions that
%   represent smaller topological uncertainty. The entropy function is
%   recommended.) The index to the best observation scheme is returned,
%   where best is defined as having the smallest average target function
%   value over the sample of parameters.
%
%   Inputs:
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be used
%   con: [ experiment struct vector ]
%       The experimental conditions for the known data
%   obj: [ objective struct matrix ]
%       The objective structures defining the information that is already
%       known
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
%   con_pos: [ experiment struct vector ]
%       The candidate experimental conditions
%   obs_pos: [ objective struct matrix ]
%       The candidate observation structures corresponding to the
%       measurement technique that will be applied under the candidate
%       experimental conditions
%   target: [ handle @(pmy) returns scalar ]
%       The goal function quantifies the value of a particular information
%       matrix
%   opts: [ options struct scalar {} ]
%       .TargetTol [ positive scalar {0.05} ]
%           The acceptable uncertainty in the estimate of the target
%           function
%       .MinTargetIter [ nonnegative integer {0} ]
%           The minimum number of samples used to consider each observation
%           scheme
%       .MaxTargetIter [ nonnegative integer {200} ]
%           The maximum number of sample used to consider each observation
%           scheme
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
%           vary during the optimization and will be considered free
%           parameters whose uncertainty will be optimized
%       .UseDoseControls [ cell vector nCon of logical vectors or positive 
%                          integer vectors | logical vector nq | positive 
%                          integer vector {[]} ]
%           Indicates the dose control parameters that will be allowed to
%           vary during the optimization and will be considered free
%           parameters whose uncertainty will be optimized
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
% 
%   Outputs
%   best: [ positive integer ]
%       Linear index into obs_pos
%   data: [ struct ]
%       Additional information

% (c) 2015 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean-up inputs
if nargin < 11
    opts = [];
end

if isempty(con)
    % Create a con to fit with the prior
    con = experimentZero(m(1));
end
if isempty(obj)
    obj = objectiveZero([1,numel(con)]);
end

% Constants
nTop = numel(m);
nCon = numel(con);
nObj = size(obj,1);
nConPos = numel(con_pos);
nObjPos = size(obs_pos,1);

% Deafult options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = [];
defaultOpts.AbsTol         = [];

defaultOpts.UseParams      = cell(nTop,1); 
for i=1:nTop; defaultOpts.UseParams{i} = 1:m(i).nk;end
defaultOpts.UseSeeds       = [];
defaultOpts.UseInputControls = [];
defaultOpts.UseDoseControls = [];

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
% Experimental conditions
con = vec(con);

% Objective functions
assert(size(obj,2) == nCon, 'KroneckerBio:TopologyProbability:ObjSize', 'Second dimension of "obj" must be equal to nCon')

% Possible experimental conditions
con_pos = vec(con_pos);

% Possible objective functions
assert(size(obs_pos,2) == nConPos, 'KroneckerBio:TopologyProbability:ObsPosSize', 'Second dimension of "obs_pos" must be equal to n_con')

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
[opts.UseSeeds, nTs] = fixUseSeeds(opts.UseSeeds, ns, nCon+nConPos);

% Ensure UseControls is a cell vector of logical vectors
[opts.UseInputControls, nTq] = fixUseControls(opts.UseInputControls, nCon+nConPos, cat(1,con.nq,con_pos.nq));
[opts.UseDoseControls, nTh] = fixUseControls(opts.UseDoseControls, nCon+nConPos, cat(1,con.nh,con_pos.nh));

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
nObjPrior = size(objPrior,1);

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
    temp = fixAbsTol(tempAbsTol{iTop}, 2, false(nCon+nConPos,1), nx(iTop), nCon+nConPos, false, opts.UseParams{iTop}, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
    [opts.AbsTol{iTop}.Sensitivity] = deal(temp{:});
    temp = fixAbsTol(tempAbsTol{iTop}, 2, true(nCon+nConPos,1), nx(iTop), nCon+nConPos, false, opts.UseParams{iTop}, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
    [opts.AbsTol{iTop}.Gradient] = deal(temp{:});
    temp = fixAbsTol(tempAbsTol{iTop}, 2, false(nCon+nConPos,1), nx(iTop), nCon+nConPos, true, opts.UseParams{iTop}, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
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
    opts.LowerBound{iTop} = fixBounds(opts.LowerBound{iTop}, opts.UseParams{iTop}, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
    opts.UpperBound{iTop} = fixBounds(opts.UpperBound{iTop}, opts.UseParams{iTop}, opts.UseSeeds, opts.UseInputControls, opts.UseDoseControls);
end

%% Distribute options for topology specific functions
assert(~isfield(opts, 'ObjWeights') || all(vec(opts.ObjWeights) == 1), 'ObjWeights is not supported for BestTopologyExperiments.')
optsAll = repmat(opts, nTop,1);
for iTop = 1:nTop
    optsAll(iTop).UseParams       = opts.UseParams{iTop};
    optsAll(iTop).UseSeeds        = opts.UseSeeds(:,1:nCon);
    optsAll(iTop).UseInputControls = opts.UseInputControls(1:nCon);
    optsAll(iTop).UseDoseControls = opts.UseDoseControls(1:nCon);
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
    optsTemp.UseSeeds = optsTemp.UseSeeds(:,1:nCon);
    optsTemp.UseInputControls = optsTemp.UseInputControls(1:nCon);
    optsTemp.UseDoseControls = optsTemp.UseDoseControls(1:nCon);
end

pmyStart = TopologyProbability(m, con, obj, objPriorParams, objPriorSeeds, objPriorInputControls, objPriorDoseControls, optsTemp);

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
    sam{iTop} = collectActiveParameters(m(iTop), con, optsAll(iTop).UseParams, optsAll(iTop).UseSeeds, optsAll(iTop).UseInputControls, optsAll(iTop).UseDoseControls);
end

for iPosCon = 1:nConPos
    for iPosObj = 1:nObjPos
        if verbose; fprintf(['Monte Carlo sampling for ' obs_pos(iPosObj,iPosCon,1).Name '...\n']); end
        
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
                [mDraw, conDraw] = updateAll(m(trueTop), con, sam{trueTop}(:,end), optsAll(trueTop).UseParams, optsAll(trueTop).UseSeeds, optsAll(trueTop).UseInputControls, optsAll(trueTop).UseDoseControls);
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
                [mDraw, conDraw] = updateAll(m(trueTop), con, sam{trueTop}(:,end), optsAll(trueTop).UseParams, optsAll(trueTop).UseSeeds, optsAll(trueTop).UseInputControls, optsAll(trueTop).UseDoseControls);
                sam{trueTop} = [sam{trueTop}, SampleParameterSpace(mDraw, conDraw, [obj; objPrior(:,:,trueTop)], opts.StepsPerCheck, optsAll(trueTop))];
            end
            
            % Construct true parameterized model
            Ttrue = sam{trueTop}(:,iSam(trueTop));
            mRand = updateAll(m(trueTop), con, Ttrue, optsAll(trueTop).UseParams, optsAll(trueTop).UseSeeds, optsAll(trueTop).UseInputControls, optsAll(trueTop).UseDoseControls);
            
            try
                % Generate data according to possible experiment
                if verbose; fprintf('Generating data from Monte Carlo model...\n'); end
                optsTemp = optsAll(iTop);
                sim = SimulateSystem(mRand, con_pos(iPosCon), obs_pos(iPosObj,iPosCon), optsTemp);
                new_data = sim.measurements;
                new_obj = obs_pos(nObjPos,nConPos).Objective(new_data);
                
                % Fit the new data
                mfit = m;
                for iTop = 1:nTop
                    if verbose; fprintf(['Fitting ' m(iTop).Name ' to existing data plus Monte Carlo data...\n']); end
                    optsTemp = optsAll(iTop);
                    optsTemp.AbsTol = [optsTemp.AbsTol; opts.AbsTol{iTop}(nCon+iPosCon)];
                    optsTemp.UseSeeds = [optsTemp.UseSeeds, opts.UseSeeds(:,nCon+iPosCon)];
                    optsTemp.UseInputControls = [optsTemp.UseInputControls; opts.UseInputControls(nCon+iPosCon)];
                    optsTemp.UseDoseControls = [optsTemp.UseDoseControls; opts.UseDoseControls(nCon+iPosCon)];
                    mfit(iTop) = FitObjective(mfit(iTop), [con; con_pos(iPosCon)], [[obj; objPrior(:,:,iTop)], [new_obj; objectiveZero([nObj+nObjPrior-1,1])]], optsTemp);
                end
                
                % Compute the topology probability
                if verbose; fprintf('Computing topology probability after adding Monte Carlo data...\n'); end
                optsTemp = opts;
                opts.NeedFit = false;
                for iTop = 1:nTop
                    optsTemp.AbsTol{iTop} = optsTemp.AbsTol{iTop}([1:nCon,iPosCon]);
                    optsTemp.UseSeeds = optsTemp.UseSeeds(:,[1:nCon,iPosCon]);
                    optsTemp.UseInputControls = optsTemp.UseInputControls([1:nCon,iPosCon]);
                    optsTemp.UseDoseControls = optsTemp.UseDoseControls([1:nCon,iPosCon]);
                end
                pmy = TopologyProbability(mfit, [con; con_pos(iPosCon)], [obj, [new_obj; objectiveZero([nObj-1,1])]], objPriorParams, objPriorSeeds, objPriorInputControls, objPriorDoseControls, optsTemp);
            catch me
                if strcmp(me.identifier, 'KroneckerBio:accumulateOde:IntegrationFailure')
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
                fitsTs{iTop,index} = collectActiveParameters(mfit(iTop), con, optsAll(iTop).UseParams, optsAll(iTop).UseSeeds, optsAll(iTop).UseInputControls, optsAll(iTop).UseDoseControls);
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

[~, best] = min(Etarget);

data.Targets               = Etarget;
data.TopologiesChosen      = trueTopAll;
data.ParametersChosen      = trueTAll;
data.MeasurementsChosen    = trueyAll;
data.FittedParameters      = fitsTs;
data.Startingpmy           = pmyStart;
data.Allpmy                = pmyAll;
data.AllTargets            = targetAll;
data.AllWeightedMeanErrors = errAll;
