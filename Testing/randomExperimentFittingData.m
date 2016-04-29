function [obj, opts, mVariant, conVariant] = randomExperimentFittingData(m, con, opts, tF, nExperiments, nTotalDataPoints)
% Generate random experiments for testing objective-based Kronecker
% functions with multiple experiments
% 
% Takes a model-experiment combination (m and con), chooses random k, s, q,
% and h from those indicated in opts.Use*, and varies them randomly up and
% down up to twofold to generate one variant model and a set of
% nExperiments variant experiments. The function then simulates the variant
% model and experiments up to time tF, gathers data from nTotalDataPoints random outputs and
% time points, and creates an objective function array obj fitting to the
% data with added error. Also returns (1) opts struct directing
% FitObjective to fit the altered parameters, and (2) the variant model and
% experiments used to generate the data in mVariant and conVariant,
% respectively.
%
% Warning: AbsTol is set to default value of 1e-9 no matter what to avoid
% dealing with standardization, which is complicated.
%
% (c) 2016 David Flowers
% This work is released under the MIT license.

% Input checks
assert(isscalar(con), 'con must be a scalar.')

conOrig = con;

usefields = {'UseParams','UseSeeds','UseInputControls','UseDoseControls'};
hasfields = isfield(opts, usefields);
assert(all(hasfields), 'getRandomModelFittingData:MissingUseFields',...
    'Fields %s are required in opts but were not found',...
    strjoin(usefields(~hasfields), ', '))

% Get standardization functions from private directory (bad practice here)
try
    currentdir = pwd;
    cd(fullfile('..','Source','private'))
    fixUseParams_ = @fixUseParams;
    fixUseSeeds_ = @fixUseSeeds;
    fixUseControls_ = @fixUseControls;
    cd(currentdir)
catch ME % Avoids getting stuck in a different directory if an error occurs
    cd(currentdir)
    rethrow(ME)
end

% Standardize Use* as ns-by-nCon logical array for UseSeeds or nCon-by-1
% cell array of n(q or h)-by-1 logical vectors otherwise
opts.UseParams = fixUseParams_(opts.UseParams, m.nk);
opts.UseSeeds = fixUseSeeds_(opts.UseSeeds, m.ns, numel(con));
opts.UseInputControls = fixUseControls_(opts.UseInputControls, numel(con), [con.nq]);
opts.UseDoseControls = fixUseControls_(opts.UseDoseControls, numel(con), [con.nh]);

% Set scalar AbsTol to avoid the complexity (obviously this isn't a great
% solution)
opts.AbsTol = 1e-9;

% Function that selects Use* parameters at random and generates a new Use*
% logical vector
selectUses = @(use) use & (randi([0,1],size(use,1),1) == 1);

% Assign random model parameters to be fit
opts.UseParams = selectUses(opts.UseParams);

% Choose parameters to randomize and fit
UseNames = {'UseSeeds','UseInputControls','UseDoseControls'};
nTypes = numel(UseNames);
for ii = 1:nTypes
    
    % Replicate Use* arrays to match number of experiments
    isSeed = strcmp(UseNames{ii}, 'UseSeeds');
    if isSeed
        opts.(UseNames{ii}) = repmat(opts.(UseNames{ii})(:,1), 1, nExperiments);
    else
        opts.(UseNames{ii}) = repmat(opts.(UseNames{ii})(1), nExperiments, 1);
    end
    
    % Choose subsets of parameters to fit for each experiment
    for i_con = 1:nExperiments
        if isSeed
            opts.(UseNames{ii})(:,i_con) = selectUses(opts.(UseNames{ii})(:,i_con));
        else
            opts.(UseNames{ii}){i_con} = selectUses(opts.(UseNames{ii}){i_con});
        end
    end
    
end

% Generate random values for parameters by varying provided parameters up
% and down two-fold
origValues = m.k(opts.UseParams);
randomMults = 2.^(2*rand(sum(opts.UseParams),1)-1);
mVariant = m;
mVariant.k(opts.UseParams) = randomMults.*origValues;
mVariant = mVariant.Update(mVariant.k);

UseTypes = {'s','q','h'};
conVariant = repmat(conOrig, nExperiments, 1);
for ii = 1:nTypes
    
    UseType = UseTypes{ii};
    UseName = UseNames{ii};
    isSeed = strcmp(UseName, 'UseSeeds');
    
    for i_con = nExperiments:-1:1
        
        % Pull out use*
        if isSeed
            Use_di = opts.(UseName)(:,i_con);
        else
            Use_di = opts.(UseName){i_con};
        end
        
        % Get original values for parameters of this type
        origValues = conOrig.(UseType)(Use_di);
        
        randomMults = 2.^(2*rand(sum(Use_di),1)-1);
        
        if any(Use_di) % Avoids problems with empties when Use_di is all false
            conVariant(i_con).(UseType)(Use_di) = randomMults.*origValues;
        end
        
        % Update experiments on last parameter type
        if ii == nTypes
            conVariant(i_con) = conVariant(i_con).Update(...
                conVariant(i_con).s,...
                conVariant(i_con).q,...
                conVariant(i_con).h);
        end
        
    end
end

% Determine random time points, experiments, and outputs for the objective
% function
timelist = tF*rand(nTotalDataPoints, 1);
outputlist = randi(m.ny, nTotalDataPoints, 1);
experimentlist = randi(nExperiments, nTotalDataPoints, 1);

sd = sdLinear(0.1, 0.01);

% Get observations
for ei = nExperiments:-1:1
    isExperiment = experimentlist == ei;
    obs(ei) = observationLinearWeightedSumOfSquares(...
        outputlist(isExperiment),...
        timelist(isExperiment),...
        sd,...
        sprintf('Experiment %d', ei));
end

% Simulate new experiments to get measurement values
% Random noise is added according to sd when calling sim.measurements
sim = SimulateSystem(mVariant, conVariant, obs, opts);
measurements = {sim.measurements};

% Create objective functions
obj = objectiveZero([nExperiments nExperiments]);
for ei = nExperiments:-1:1
    obj(ei,ei) = obs(ei).Objective(measurements{ei});
end

end
