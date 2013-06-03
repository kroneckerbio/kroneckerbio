function [bestCons data] = BestParameterExperiment(m, con, obj, conPos, objPos, goal, opts, F, EFs)
%BESTPARAMETEREXPERIMENT Determine which experiments will most efficiently
%   minimize a goal function of the fisher information matrix (FIM) of the
%   objective function. Traditionally, this algorithm is used to minimize
%   some uncertainty in the parameters directions. (The variance-covariance
%   matrix of the parameters is approximated by the inverse of the
%   FIM if the objective function is some form of least squares.) 
%
%   This algorithm computes the current FIM according to the current
%   experiments defined by con and obj. It then adds the expected FIM from
%   each of the possible experiments defined by conPos and objPos and asks
%   the goal function to evaluate the resulting FIM. Finally, the function
%   returns the index of the experiment that has the lowest expected goal
%   function value.
%
%   [...] = BestParameterExperiment(m, con, obj, conPos, objPos, goal, opts)
%
%   Inputs:
%       m   - The KroneckerBio model of interest
%       con - A vector of structures, with each item corresponding to the
%             experimental conditions of the experiments that have already
%             been performed
%       obj - A vector of structures, with each item corresponding to the
%             objective function containing the measurements for the
%             experiment
%       conPos - A vector of experimental conditions that will be examined
%       objPos - A vector of objective functions corresponding to the
%                measurement techique that will be applied under the
%                possible experimental conditions
%       goal   - A function handle @(F) that returns a scalar evaluation of
%                a hessian. This function chooses the best experiments
%                based on which set has the lowest goal function.
%       opts - Optional function options
%           UseParams - Vector of indexes indicating the rate constants
%                       that will be considered in the hessian.
%           UseICs    - Vector of indexes indicating the initial
%                       concentrations whose sensitivities will be
%                       considered
%           UseModelICs   - Boolean that determines whether to use the
%                           initial conditions of the model or the
%                           conditions. This is used both for simulating
%                           the model and calculating the sensitivities of
%                           the IC parameters.
%           ReturnCount   - Integer indicating the number of experiments
%                           to return
%           MaxGreedySize - Integer indicating how many experiments should
%                           be considered simultaneously. If it set to 1,
%                           the algorithm is perfectly greedy, choosing
%                           the best first experiment and then choosing
%                           the second best, given that the first was
%                           already chosen. If it is greater than or equal
%                           to ReturnCount, then the algorithm is
%                           perfectly optimal, choosing the combination
%                           that all together minimizes the goal function.
%           TerminalGoal  - Scalar that terminates the search if the
%                           expected goal is less than or equal to this
%                           value. Default = -inf
%           AllowRepeats  - Boolean that determines whether or to allow
%                           experiments to be repeated
%           Verbose       - Print progress to command window
%
%   Outputs:
%       bestCons = FindBestExperiment(m, con, obj, conPos, objPos, goal, opts)
%       	bestCons - A vector of the indexes to the best experiments
%
%       [bestCons data] = FindBestExperiment(m, con, obj, conPos, objPos, goal, opts)
%        	data - Includes a structure that holds the hessians and
%                  goal function values at each iteration. This data
%                  may be useful for analysis. It also can require a lot of
%                  memory.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean up inputs
assert(nargin >= 6, 'KroneckerBio:BestParameterExperiment:TooFewInputs', 'BestParameterExperiment requires at least 6 input arguments')
if nargin < 9
    EFs = [];
    if nargin < 8
        F = [];
        if nargin < 7
            opts = [];
        end
    end
end

assert(isscalar(m), 'KroneckerBio:BestParameterExperiment:MoreThanOneModel', 'The model structure must be scalar')

% Default options
defaultOpts.Verbose         = 1;

defaultOpts.RelTol          = NaN;
defaultOpts.AbsTol          = NaN;
defaultOpts.UseModelICs     = false;
defaultOpts.UseModelInputs  = false;

defaultOpts.UseParams       = 1:m.nk;
defaultOpts.UseICs          = [];
defaultOpts.UseControls     = [];

defaultOpts.ObjWeights      = ones(size(obj));

defaultOpts.Normalized      = true;
defaultOpts.UseAdjoint      = true;

defaultOpts.UseExperiments  = true(size(objPos));
defaultOpts.Cost            = zeros(size(objPos));
defaultOpts.ReturnCount     = 1;        % Number of experiments to return
defaultOpts.MaxGreedySize   = 1;        % Number of experiments to consider at once for the greedy search. Inf = not greedy
defaultOpts.Budget          = inf;
defaultOpts.MaxGreedyBudget = inf;      % Amount of budget to consider in a single step. Inf = not greedy
defaultOpts.TerminalGoal    = -inf;
defaultOpts.AllowRepeats    = true;

opts = mergestruct(defaultOpts, opts);

verbose = logical(opts.Verbose);
opts.Verbose = max(opts.Verbose-1,0);

% Constants
nPosCon = numel(conPos);
nPosObj = size(objPos,1);
nPos = numel(objPos);

%% Fix UseExperiments
opts.UseExperiments = fixUseExperiments(opts.UseExperiments, nPosObj, nPosCon);
remainingCons = vec(find(opts.UseExperiments)); % The indexes of the conditions allowed to be chosen
nPos = numel(remainingCons);
% The block size is bounded by MaxGreedySize, ReturnCount, and the total
% number of experiments that can be completed in the budget
if opts.AllowRepeats
    maxInBudget = floor(opts.Budget / min(opts.Cost(remainingCons)));
    maxReturnCount = min(opts.ReturnCount, maxInBudget);
    maxInBudgetStep = floor(opts.MaxGreedyBudget / min(opts.Cost(remainingCons)));
    blockSize = min([opts.MaxGreedySize, maxReturnCount maxInBudgetStep]);
    maxIterations = ceil(maxReturnCount / blockSize);
    maxSearchSize = nPos.^blockSize;
else
    maxInBudget = find(cumsum(sort(opts.Cost(remainingCons))) <= opts.Budget, 1, 'last'); % Bounded by nPos
    if isempty(maxInBudget); maxInBudget = 0; end % Budget may not allow for any experiments
    maxReturnCount = min(opts.ReturnCount, maxInBudget);
    maxInBudgetStep = floor(opts.MaxGreedyBudget / min(opts.Cost(remainingCons)));
    blockSize = min([opts.MaxGreedySize, maxReturnCount maxInBudgetStep]);
    maxIterations = ceil(maxReturnCount / blockSize);
    maxSearchSize = nchoosek(nPos, blockSize);
end

%% Compute old information
if isempty(F)
    F = ObjectiveInformation(m, con, obj, opts);
end

%% Fetch expected information
if isempty(EFs)
    [unused EFs] = ObjectiveInformation(m, conPos, objPos, opts);
end

%% Find best experiment
if verbose; fprintf('Combining information to find best experiment...\n'); end
    
bestCons = zeros(maxReturnCount, 1);
bestGoal = goal(F);

if verbose
    fprintf('Starting goal is %d\n', bestGoal)
end

remainingCount  = maxReturnCount; % The number of best experiments remaining to be found
remainingBudget = opts.Budget;    % The total budget for the experiments
currentF = F; % Updates after each block to hold the current expected information
currentIteration = 0; % Keeps track of which iteration the loop is on

% Allocate algorithm storage if it is requested
if nargout >= 2
    allInds = cell(maxIterations, 1);
    allFIMs = cell(maxIterations, 1);
    allGoals = cell(maxIterations, 1);
    bestFIMs = cell(maxIterations, 1);
    bestGoals = zeros(maxIterations, 1);
end

while remainingCount > 0 && bestGoal > opts.TerminalGoal && any(remainingBudget >= opts.Cost(remainingCons))
    if verbose; fprintf('%d experiments and %d budget remaining\n', remainingCount, remainingBudget); end
    currentIteration = currentIteration + 1;
    
    % CurrentBlockSize is the set size of the experiments to optimize over.
    % It is equal to blockSize until the number of experiments remaining to
    % be found is less than blockSize.
    currentBudget    = min(opts.MaxGreedyBudget, remainingBudget);
    if opts.AllowRepeats
        currentMaxInBudget = floor(currentBudget / min(opts.Cost(remainingCons)));
    else
        currentMaxInBudget = find(cumsum(sort(opts.Cost(remainingCons))) <= currentBudget, 1, 'last');
    end
    currentBlockSize = min([blockSize, remainingCount, currentMaxInBudget]);
    
    % Combinatorics
    clear conList goalValues
    conList = generateBlock(nPos, remainingCons, currentBlockSize, opts.AllowRepeats, opts.Cost, currentBudget);
    nSearch = size(conList, 1); % Total number of possible sets
    
    if nargout >= 2
        % Initialize search variables
        goalValues = zeros(nSearch, 1);
%         newEFs = cell(nSearch,1);
    end
    
    bestGoal = inf;
    bestSetInd = 1;
    for iSearch = 1:nSearch
        % Compute expected hessian
        newEF = currentF;
        for iBlock = 1:nnz(conList(iSearch,:))
            newEF = newEF + EFs{conList(iSearch,iBlock)};
        end
        
        % Compute goal and keep it if it is better
        goalValue = goal(newEF);
        if goalValue < bestGoal
            bestGoal = goalValue;
            bestSetInd = iSearch;
        end
        
        if nargout >= 2
            % Store all results
            goalValues(iSearch) = goalValue;
%             newEFs{iSearch} = newEF;
        end
    end
    
    % Find corresponding best set
    bestSet = nonzeros(conList(bestSetInd,:));                    % Experiments that are part of the best set (indexes of remainingCons)
    blockInd1 = maxReturnCount - remainingCount + 1;              % Where in bestCons to starting storing these experiment indexes
    blockInd2 = maxReturnCount - remainingCount + numel(bestSet); % Where in bestCons is the last index to store
    bestCons(blockInd1:blockInd2) = bestSet;                      % Store values of the best set
    
    if verbose
        for i = 1:numel(bestSet)
            [iObj iCon] = ind2sub([nPosObj,nPosCon], bestSet(i));
            fprintf([conPos(iCon).Name ' ' objPos(bestSet(i)).Name ' was chosen\n']);
        end
        fprintf('Expected goal is %d\n', bestGoal);
    end
    
    % Apply information for best set
    for i = 1:numel(bestSet)
        currentF = currentF + EFs{bestSet(i)};
    end
    
    % Remove conditions if repeats are not allowed
    if ~opts.AllowRepeats
        remainingCons(logical(lookup(remainingCons,bestSet))) = [];
    end
    
    % Store iteration data only if it is requested
    if nargout >= 2
        allInds{currentIteration}   = conList;
%        allFIMs{currentIteration}   = newEFs;
        allGoals{currentIteration}  = goalValues;
        bestFIMs{currentIteration}  = currentF;
        bestGoals(currentIteration) = bestGoal;
    end
    
    % Decrement
    remainingCount  = remainingCount - currentBlockSize;
    remainingBudget = remainingBudget - sum(opts.Cost(bestSet));
end

%% Work-down
% Clean up if terminated on TerminalGoal rather than ReturnCount
actual = (bestCons ~= 0);
bestCons = bestCons(actual);

if nargout >= 2
    data.StartingFIM = F;                             % FIM of system before possible experiments are applied
    data.AllFIMs     = EFs;                           % FIMs for each possible experiment
    data.SearchSets  = allInds(1:currentIteration);   % Each block of experiments searched
%    data.SearchFIMs  = allFIMs(1:currentIteration);   % Sum of FIMs in each search
    data.SearchGoals = allGoals(1:currentIteration);  % All the goal values of the possible iterations
    data.BestFIMs    = bestFIMs(1:currentIteration);  % FIM after each iteration
    data.BestGoals   = bestGoals(1:currentIteration); % Goal after each iteration
end

end

function list = generateBlock(nPos, conPos, blockSize, allowRepeats, cost, budget)
% Minimize memory footprint of list
if nPos < 2^8
    conPos = uint8(conPos);
elseif nPos < 2^16
    conPos = uint16(conPos);
elseif nPos < 2^32
    conPos = uint32(conPos);
elseif nPos < 2^64
    conPos = uint64(conPos);
end

minCost = min(min(cost(conPos)));

% Create set of single experiments
list = cell(blockSize,1);
list{1} = conPos; % The list of all blocks that will be grown
setCost = cost(conPos); % The cost of each block in the current set
overpriced = (setCost > budget); % Blocks that will be removed and discarded
list{1}(overpriced) = []; % The list of all blocks that will be grown
setCost(overpriced) = []; % The cost of each block in the current set

% Loop over blockSize
for iDim = 2:blockSize
    growable = (setCost + minCost <= budget); % Blocks that will be removed and carried to the next stage
    nGrow = sum(growable);
    addable = (cost(conPos) <= max(budget - setCost));
    nAdd = sum(addable);
    list{iDim} = reshape(repmat(reshape(list{iDim-1}(growable,:), 1,(iDim-1)*nGrow), nAdd,1), nAdd*nGrow,iDim-1); % Repeat each line nPos to prepare to receive conPos
    list{iDim} = [list{iDim}, repmat(conPos(addable), nGrow,1)]; % Concatenate new indexes
    
    % Remove grown blocks from previous set
    list{iDim-1}(growable,:) = [];
    
    % Remove permutations by forcing an ascending block
    repeats = list{iDim}(:,iDim-1) >= list{iDim}(:,iDim);
    list{iDim}(repeats,:) = [];
    
    % Remove repeats
    if ~allowRepeats
        repeats = any(bsxfun(@eq, list{iDim}(:,1:iDim-2), list{iDim}(:,iDim)),2); % Find all rows that have repeats with new indexes
        list{iDim}(repeats,:) = []; % Remove repeats
    end
    
    % Remove overpriced blocks
    setCost = sum(cost(list{iDim}),2);
    overpriced = (setCost > budget);
    list{iDim}(overpriced,:) = [];
    setCost(overpriced) = [];
end

% Concatenate blocks into single list, padding with zeros
list = cellfun(@(A)cat(2, A, zeros(size(A,1), blockSize-size(A,2))), list, 'UniformOutput', false);
list = cat(1, list{:});
end

% function list = generateBlock(conPos, blockSize, allowRepeats, cost, budget)
% 
% nPos = numel(conPos);
% listSize = repmat(nPos, blockSize,1);
% 
% if allowRepeats
%     list = false(listSize);
%     inds = combinator(nPos, blockSize, 'c', 'r'); % Generate all possible sets that fit in blockSize
%     inds = mat2cell(inds, size(inds,1), ones(size(inds,1),1)); % Distribute each column over a cell vector
%     list(sub2ind(listSize, inds{:})) = true; % Mark each valid set of experiments in logical hypercube
% else
%     list = false(listSize);
%     inds = combinator(nPos, blockSize, 'c', 'n'); % Generate all possible sets that fit in blockSize
%     inds = mat2cell(inds, size(inds,1), ones(size(inds,1),1)); % Distribute each column over a cell vector
%     list(sub2ind(listSize, inds{:})) = true; % Mark each valid set of experiments in logical hypercube
% end
% 
% for iDim = blockSize:-1:1
% end
% 
% end

% function list = generateBlock(conPos, blockSize, allowRepeats, costPos, budget)
% nElements = 0;
% list = cell(0,1);
% 
% % Run the recusive list builder
% recursiveBuild(zeros(0,1), conPos, budget);
% 
% % Remove empty elements
% list = list(1:nElements);
% 
%     function recursiveBuild(recursiveList, remainingCons, remainingBudget)
%         % Recusively generate the tree that adds each experiment until no more can be added
%         canBeExtended = false;
%         if numel(recursiveList) < blockSize
%             % There is room for more experiments, see if there is enough budget
%             for iPos = 1:numel(remainingCons)
%                 if remainingBudget >= costPos(iPos)
%                     % Budget allows for this experiment, add it and try to add more recursively
%                     if allowRepeats
%                         % Send the full list of experiments to the next iteration
%                         recursiveBuild([recursiveList; remainingCons(iPos)], remainingCons, remainingBudget - costPos(iPos));
%                     else
%                         % Remove repeats from possible experiments
%                         nextRemainingCons = remainingCons;
%                         nextRemainingCons(iPos) = [];
%                         recursiveBuild([recursiveList; remainingCons(iPos)], nextRemainingCons, remainingBudget - costPos(iPos));
%                     end
%                     canBeExtended = true;
%                 end
%             end
%         end
%         
%         if ~canBeExtended
%             % This is a full length set. Add it to the list of sets.
%             % Increment counter
%             nElements = nElements + 1;
%             if numel(list) < nElements
%                 % Double size of list
%                 list = [list; cell(max(nElements,1),1)];
%             end
%             
%             % Add element
%             list{nElements} = recursiveList;
%         end
%     end
% 
% end