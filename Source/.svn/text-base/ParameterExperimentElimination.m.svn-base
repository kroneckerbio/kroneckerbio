function useExperiments = ParameterExperimentElimination(m, con, obj, conPos, objPos, opts, F, EFs)
%ParameterExperimentElimination Remove experiments from consideration in
%   optimal experimental design based on the ability of some experiments to
%   provide greater information in all directions.

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
defaultOpts.UseModelICs    = false;
defaultOpts.UseModelInputs = false;

defaultOpts.UseParams      = 1:m.nk;
defaultOpts.UseICs         = [];
defaultOpts.UseControls    = [];

defaultOpts.ObjWeights     = ones(size(obj));

defaultOpts.Normalized     = true;
defaultOpts.UseAdjoint     = true;

defaultOpts.UseExperiments      = true(size(obj));
defaultOpts.Cost                = zeros(size(obj));
defaultOpts.Budget              = inf;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nPosCon = numel(conPos);
nPosObj = size(objPos,1);

%% Fix UseExperiments
useExperiments = fixUseExperiments(opts.UseExperiments, nPosObj, nPosCon);
remainingCons = find(opts.UseExperiments); % The indexes of the conditions allowed to be chosen
nPos = numel(remainingCons);

%% Compute existing information
if isempty(F)
    F = ObjectiveInformation(m, con, obj, opts);
end
nT = size(F,1);

%% Fetch expected hessians
if isempty(EFs)
    [unused EFs] = ObjectiveInformation(m, conPos, objPos, opts);
end

%% Single elimination
conList = remainingCons(combinator(nPos, 2, 'c', 'n')); % Generate list of pairs of experiments, no repeats
for iPos = 1:size(conList,1)
    firstCon = conList(iPos,1);
    secondCon = conList(iPos,2);
    [iObj1 iCon1] = ind2sub([nPosObj nPosCon], firstCon);
    [iObj2 iCon2] = ind2sub([nPosObj nPosCon], secondCon);
    
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
            if verbose; fprintf('%s %s (#%i) dominates %s %s (#%i)\n', conPos(iCon1).Name, objPos(firstCon).Name, firstCon, conPos(iCon2).Name, objPos(secondCon).Name, secondCon); end
            useExperiments(secondCon) = false;
        elseif all(lambdadiff <= sqrt(opts.RelTol*32)) && (zeroDom == 2 || zeroDom == 0) && (infDom == 2 || infDom == 0) && (opts.Cost(firstCon) >= opts.Cost(secondCon))
            if verbose; fprintf('%s %s (#%i) dominates %s %s (#%i)\n', conPos(iCon2).Name, objPos(secondCon).Name, secondCon, conPos(iCon1).Name, objPos(firstCon).Name, firstCon); end
            useExperiments(firstCon) = false;
        end
    end

end