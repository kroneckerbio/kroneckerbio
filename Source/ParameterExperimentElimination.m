function UseExperiments = ParameterExperimentElimination(m, con, obj, con_pos, obj_pos, opts, F, EFs)
%ParameterExperimentElimination removes experiments from consideration in
%   optimal experimental design based on the ability of some experiments to
%   provide greater information in all directions.
%
%   UseExperiments = ParameterExperimentElimination(m, con, obj, con_pos, obj_pos, opts, F, EFs)
%
%   Inputs:
%   m: [ model struct scalar ]
%       The KroneckerBio model that will be used
%   con: [ experiment struct vector ]
%       The experimental conditions for the known data
%   obj: [ objective struct matrix ]
%       The objective structures defining the information that is already
%       known
%   con_pos: [ experiment struct vector ]
%       The candidate experimental conditions
%   obj_pos: [ objective struct matrix ]
%       The candidate objective structures corresponding to the measurement
%       technique that will be applied under the candidate experimental
%       conditions
%   opts: [ options struct scalar {} ]
%       .UseModelSeeds [ logical scalar {false} ]
%           Indicates that the model's seed parameters should be used
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
%       .UseControls [ cell vector nCon of logical vectors or positive 
%                      integer vectors | logical vector nq | positive 
%                      integer vector {[]} ]
%           Indicates the input control parameters that will be allowed to
%           vary during the optimization
%       .UseExperiments [ logical matrix size of obj_pos |
%                         linear index vector into obj_pos
%                         {true(size(obj_pos))} ]
%           Indicates which candidate experiments are actually available
%       .Cost [ matrix size of obj_pos ]
%           The cost of each candidate experiment
%       .Verbose [ nonnegative integer scalar {1} ]
%           Bigger number displays more progress information
%   F: [ matrix nT by nT {} ]
%       The Fisher information for the known data can be supplied so that
%       this function will not recalculate it
%   EFs: [ cell matrix size of objPos {} ]
%       The exprected information matrices for the candidate experiments
%       can be supplied so that this function will not recalculate them
% 
%   Outputs:
%   UseExperiments = ParameterExperimentElimination(m, con, obj, con_pos, obj_pos, opts, F, EFs)
%       A vector of the indexes to the best experiments

% (c) 2013 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
assert(nargin >= 5, 'KroneckerBio:ParameterExperimentElimination:TooFewInputs', 'ParameterExperimentElimination requires at least 5 input arguments')
if nargin < 8
    EFs = [];
    if nargin < 7
        F = [];
        if nargin < 6
            opts = [];
        end
    end
end

% Default options
defaultOpts.Verbose        = 1;

defaultOpts.RelTol         = NaN;
defaultOpts.AbsTol         = NaN;
defaultOpts.UseModelSeeds  = false;
defaultOpts.UseModelInputs = false;

defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseSeeds       = [];
defaultOpts.UseControls    = [];

defaultOpts.ObjWeights     = ones(size(obj));

defaultOpts.Normalized     = true;

defaultOpts.UseExperiments      = true(size(obj));
defaultOpts.Cost                = zeros(size(obj));

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nPosCon = numel(con_pos);
nPosObj = size(obj_pos,1);

%% Fix UseExperiments
UseExperiments = fixUseExperiments(opts.UseExperiments, nPosObj, nPosCon);
remainingCons = find(opts.UseExperiments); % The indexes of the conditions allowed to be chosen
nPos = numel(remainingCons);

%% Compute existing information
if isempty(F)
    F = ObjectiveInformation(m, con, obj, opts);
end
nT = size(F,1);

%% Fetch expected hessians
if isempty(EFs)
    [~, EFs] = ObjectiveInformation(m, con_pos, obj_pos, opts);
end

%% Single elimination
conList = remainingCons(combinator(nPos, 2, 'c', 'n')); % Generate list of pairs of experiments, no repeats
for iPos = 1:size(conList,1)
    firstCon = conList(iPos,1);
    secondCon = conList(iPos,2);
    [iObj1, iCon1] = ind2sub([nPosObj nPosCon], firstCon);
    [iObj2, iCon2] = ind2sub([nPosObj nPosCon], secondCon);
    
    firstZeros  = all(EFs{firstCon} == 0);
    secondZeros = all(EFs{secondCon} == 0);
    firstInf    = any(EFs{firstCon} == inf);
    secondInf   = any(EFs{secondCon} == inf);
    
    % Domination: 0 tied, 1 first dominates, 2 second dominates, -1 complimentary
    if all(firstZeros == secondZeros)
        zeroDom = 0;
    elseif all(firstZeros >= secondZeros)
        zeroDom = 2;
    elseif all(firstZeros <= secondZeros)
        zeroDom = 1;
    else
        zeroDom = -1;
    end
    if all(firstInf == secondInf)
        infDom = 0;
    elseif all(firstInf >= secondInf)
        infDom = 1;
    elseif all(firstInf <= secondInf)
        infDom = 2;
    else
        infDom = -1;
    end
    
    % Only continue if they are not already mutually exclusive
    if ((zeroDom == 1 || zeroDom == 0) && (infDom == 1 || infDom == 0) && (opts.Cost(firstCon) <= opts.Cost(secondCon))) || ...
            (zeroDom == 2 || zeroDom == 0) && (infDom == 2 || infDom == 0) && (opts.Cost(firstCon) >= opts.Cost(secondCon))
        % Find difference in information between experiments
        Fdiff = EFs{firstCon} - EFs{secondCon};
        
        % Eigendecompose the difference
        lambdadiff = infoeig(Fdiff);
        
        % Compare for dominace
        if all(lambdadiff >= -sqrt(opts.RelTol*32)) && (zeroDom == 1 || zeroDom == 0) && (infDom == 1 || infDom == 0) && (opts.Cost(firstCon) <= opts.Cost(secondCon))
            if verbose; fprintf('%s %s (#%i) dominates %s %s (#%i)\n', con_pos(iCon1).Name, obj_pos(firstCon).Name, firstCon, con_pos(iCon2).Name, obj_pos(secondCon).Name, secondCon); end
            UseExperiments(secondCon) = false;
        elseif all(lambdadiff <= sqrt(opts.RelTol*32)) && (zeroDom == 2 || zeroDom == 0) && (infDom == 2 || infDom == 0) && (opts.Cost(firstCon) >= opts.Cost(secondCon))
            if verbose; fprintf('%s %s (#%i) dominates %s %s (#%i)\n', con_pos(iCon2).Name, obj_pos(secondCon).Name, secondCon, con_pos(iCon1).Name, obj_pos(firstCon).Name, firstCon); end
            UseExperiments(firstCon) = false;
        end
    end
end
