function [best, data] = BestParameterExperiment(m, con, obj, con_pos, obs_pos, goal, opts, F, EFs)
%BestParameterExperiment Determine which experiments will most efficiently
%   minimize a goal function of the fisher information matrix (FIM) of the
%   objective function. Traditionally, this algorithm is used to minimize
%   some uncertainty in the parameter. (The variance-covariance matrix of
%   the parameters is approximated by the inverse of the FIM if the
%   objective function is some form of least squares.)
%
%   This algorithm computes the current FIM according to the current
%   experiments defined by con and obj. It then adds the expected FIM from
%   each of the possible experiments defined by conPos and objPos and asks
%   the goal function to evaluate the resulting FIM. Finally, the function
%   returns the index of the experiment that has the lowest expected goal
%   function value.
%
%   best = BestParameterExperiment(m, con, obj, conPos, objPos, goal, opts)
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
%   obs_pos: [ objective struct matrix ]
%       The candidate observation structures corresponding to the
%       measurement technique that will be applied under the candidate
%       experimental conditions
%   goal: [ handle @(F) returns scalar ]
%       The goal function quantifies the value of a particular information
%       matrix
%   opts: [ options struct scalar {} ]
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
%       .ReturnCount [ whole scalar {1} ]
%           The number of a best experiments to return
%       .MaxGreedySize [ natural scalar {1} ]
%           The number of experiments that should be considered in a single
%           iteration. If it set to 1, the algorithm is perfectly greedy,
%           choosing the best first experiment and then choosing the second
%           best, given that the first was already chosen. If it is greater
%           than or equal to ReturnCount, then the algorithm is perfectly
%           optimal, choosing the combination that all together minimizes
%           the goal function.
%       .UseExperiments [ logical matrix size of obj_pos |
%                         linear index vector into obj_pos
%                         {true(size(obj_pos))} ]
%           Indicates which candidate experiments are actually available
%       .TerminalGoal [ scalar {-inf} ]
%           The search is terminated if the expected goal is less than or
%           equal to this value
%       .Cost [ matrix size of obj_pos ]
%           The cost of each candidate experiment
%       .Budget [ nonegative scalar {inf} ]
%           The total budget that can be spent on the experiments
%       .MaxGreedyBudget [ nonnegative scalar {inf} ]
%           The maximum budget that should be considered in a single
%           iteration
%       .AllowRepeats [ logical scalar {true} ]
%           Determines if experiments are allowed to be repeated
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
%   best = FindBestExperiment(m, con, obj, conPos, objPos, goal, opts)
%       A vector of the indexes to the best experiments
%
%   [best, data] = FindBestExperiment(m, con, obj, conPos, objPos, goal, opts)
%       Returns intermediate results from the algorithm

% (c) 2015 David R Hagen & Bruce Tidor
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

if isempty(con)
    con = experimentZero(0);
end
if isempty(obj)
    obj = objectiveZero([0,numel(con)]);
end

% Default options
defaultOpts.Verbose         = 1;

defaultOpts.RelTol         = [];
defaultOpts.AbsTol         = [];

defaultOpts.UseParams        = 1:m.nk;
defaultOpts.UseSeeds         = [];
defaultOpts.UseInputControls = [];
defaultOpts.UseDoseControls  = [];

defaultOpts.ObjWeights      = ones(size(obj));

defaultOpts.Normalized      = true;
defaultOpts.UseAdjoint      = true;

defaultOpts.UseExperiments  = true(size(obs_pos));
defaultOpts.Cost            = zeros(size(obs_pos));
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
n_pos_con = numel(con_pos);
n_pos_obs = size(obs_pos,1);

%% Fix UseExperiments
opts.UseExperiments = fixUseExperiments(opts.UseExperiments, n_pos_obs, n_pos_con);
remaining_cons = vec(find(opts.UseExperiments)); % The indexes of the conditions allowed to be chosen
n_pos = numel(remaining_cons);
% The block size is bounded by MaxGreedySize, ReturnCount, and the total
% number of experiments that can be completed in the budget
if opts.AllowRepeats
    max_in_budget = floor(opts.Budget / min(opts.Cost(remaining_cons)));
    max_return_count = min(opts.ReturnCount, max_in_budget);
    max_in_budget_step = floor(opts.MaxGreedyBudget / min(opts.Cost(remaining_cons)));
    block_size = min([opts.MaxGreedySize, max_return_count max_in_budget_step]);
    max_iterations = ceil(max_return_count / block_size);
    max_search_size = n_pos.^block_size;
else
    max_in_budget = find(cumsum(sort(opts.Cost(remaining_cons))) <= opts.Budget, 1, 'last'); % Bounded by nPos
    if isempty(max_in_budget); max_in_budget = 0; end % Budget may not allow for any experiments
    max_return_count = min(opts.ReturnCount, max_in_budget);
    max_in_budget_step = floor(opts.MaxGreedyBudget / min(opts.Cost(remaining_cons)));
    block_size = min([opts.MaxGreedySize, max_return_count max_in_budget_step]);
    max_iterations = ceil(max_return_count / block_size);
    max_search_size = nchoosek(n_pos, block_size);
end

%% Compute old information
if isempty(F)
    F = ObjectiveInformation(m, con, obj, opts);
end

%% Fetch expected information
if isempty(EFs)
    [~, EFs] = ObjectiveInformation(m, con_pos, obs_pos, opts);
end

%% Find best experiment
best = zeros(max_return_count, 1);
best_goal = goal(F);

if verbose; fprintf('Starting goal is %d\n', best_goal); end

remaining_count   = max_return_count; % The number of best experiments remaining to be found
remaining_budget  = opts.Budget;    % The total budget for the experiments
current_F         = F; % Updates after each block to hold the current expected information
current_iteration = 0; % Keeps track of which iteration the loop is on

% Allocate algorithm storage if it is requested
if nargout >= 2
    all_inds   = cell(max_iterations, 1);
%     allFIMs    = cell(maxIterations, 1);
    all_goals  = cell(max_iterations, 1);
    best_Fs    = cell(max_iterations, 1);
    best_goals = zeros(max_iterations, 1);
end

while remaining_count > 0 && best_goal > opts.TerminalGoal && any(remaining_budget >= opts.Cost(remaining_cons))
    if verbose; fprintf('%d experiments and %d budget remaining\n', remaining_count, remaining_budget); end
    current_iteration = current_iteration + 1;
    
    % CurrentBlockSize is the set size of the experiments to optimize over.
    % It is equal to blockSize until the number of experiments remaining to
    % be found is less than blockSize.
    current_budget    = min(opts.MaxGreedyBudget, remaining_budget);
    if opts.AllowRepeats
        current_max_in_budget = floor(current_budget / min(opts.Cost(remaining_cons)));
    else
        current_max_in_budget = find(cumsum(sort(opts.Cost(remaining_cons))) <= current_budget, 1, 'last');
    end
    current_block_size = min([block_size, remaining_count, current_max_in_budget]);
    
    % Combinatorics
    clear conList goalValues
    con_list = generateBlock(n_pos, remaining_cons, current_block_size, opts.AllowRepeats, opts.Cost, current_budget);
    n_search = size(con_list, 1); % Total number of possible sets
    
    if nargout >= 2
        % Initialize search variables
        goal_values = zeros(n_search, 1);
%         newEFs      = cell(nSearch,1);
    end
    
    best_goal = inf;
    best_set_ind = 1;
    for iSearch = 1:n_search
        % Compute expected information
        newEF = current_F;
        for iBlock = 1:nnz(con_list(iSearch,:))
            newEF = newEF + EFs{con_list(iSearch,iBlock)};
        end
        
        % Compute goal and keep it if it is better
        goalValue = goal(newEF);
        if goalValue < best_goal
            best_goal = goalValue;
            best_set_ind = iSearch;
        end
        
        if nargout >= 2
            % Store all results
            goal_values(iSearch) = goalValue;
%             newEFs{iSearch}      = newEF;
        end
    end
    
    % Find corresponding best set
    best_set = nonzeros(con_list(best_set_ind,:));                    % Experiments that are part of the best set (indexes of remainingCons)
    block_ind1 = max_return_count - remaining_count + 1;              % Where in bestCons to starting storing these experiment indexes
    block_ind2 = max_return_count - remaining_count + numel(best_set); % Where in bestCons is the last index to store
    best(block_ind1:block_ind2) = best_set;                      % Store values of the best set
    
    if verbose
        for i = 1:numel(best_set)
            [iObj, iCon] = ind2sub([n_pos_obs,n_pos_con], best_set(i));
            fprintf([con_pos(iCon).Name ' ' obs_pos(best_set(i)).Name ' was chosen\n']);
        end
        fprintf('Expected goal is %d\n', best_goal);
    end
    
    % Apply information for best set
    for i = 1:numel(best_set)
        current_F = current_F + EFs{best_set(i)};
    end
    
    % Remove conditions if repeats are not allowed
    if ~opts.AllowRepeats
        remaining_cons(logical(lookupmember(remaining_cons,best_set))) = [];
    end
    
    % Store iteration data only if it is requested
    if nargout >= 2
        all_inds{current_iteration}   = con_list;
%        allFIMs{currentIteration}     = newEFs;
        all_goals{current_iteration}  = goal_values;
        best_Fs{current_iteration}    = current_F;
        best_goals(current_iteration) = best_goal;
    end
    
    % Decrement
    remaining_count  = remaining_count - current_block_size;
    remaining_budget = remaining_budget - sum(opts.Cost(best_set));
end

%% Work-down
% Clean up if terminated on TerminalGoal rather than ReturnCount
actual = (best ~= 0);
best = best(actual);

if nargout >= 2
    data.StartingFIM = F;                             % FIM of system before possible experiments are applied
    data.AllFIMs     = EFs;                           % FIMs for each possible experiment
    data.SearchSets  = all_inds(1:current_iteration);   % Each block of experiments searched
%    data.SearchFIMs  = allFIMs(1:currentIteration);   % Sum of FIMs in each search
    data.SearchGoals = all_goals(1:current_iteration);  % All the goal values of the possible iterations
    data.BestFIMs    = best_Fs(1:current_iteration);  % FIM after each iteration
    data.BestGoals   = best_goals(1:current_iteration); % Goal after each iteration
end

end

function list = generateBlock(nPos, con_pos, block_size, allow_repeats, cost, budget)
% Minimize memory footprint of list
if nPos < 2^8
    con_pos = uint8(con_pos);
elseif nPos < 2^16
    con_pos = uint16(con_pos);
elseif nPos < 2^32
    con_pos = uint32(con_pos);
elseif nPos < 2^64
    con_pos = uint64(con_pos);
end

minCost = min(min(cost(con_pos)));

% Create set of single experiments
list = cell(block_size,1);
list{1} = con_pos; % The list of all blocks that will be grown
set_cost = cost(con_pos); % The cost of each block in the current set
overpriced = (set_cost > budget); % Blocks that will be removed and discarded
list{1}(overpriced) = []; % The list of all blocks that will be grown
set_cost(overpriced) = []; % The cost of each block in the current set

% Loop over blockSize
for i_dim = 2:block_size
    growable = (set_cost + minCost <= budget); % Blocks that will be removed and carried to the next stage
    n_grow = sum(growable);
    addable = (cost(con_pos) <= max(budget - set_cost));
    n_add = sum(addable);
    list{i_dim} = reshape(repmat(reshape(list{i_dim-1}(growable,:), 1,(i_dim-1)*n_grow), n_add,1), n_add*n_grow,i_dim-1); % Repeat each line nPos to prepare to receive conPos
    list{i_dim} = [list{i_dim}, repmat(con_pos(addable), n_grow,1)]; % Concatenate new indexes
    
    % Remove grown blocks from previous set
    list{i_dim-1}(growable,:) = [];
    
    % Remove permutations by forcing an ascending block
    repeats = list{i_dim}(:,i_dim-1) >= list{i_dim}(:,i_dim);
    list{i_dim}(repeats,:) = [];
    
    % Remove repeats
    if ~allow_repeats
        repeats = any(bsxfun(@eq, list{i_dim}(:,1:i_dim-2), list{i_dim}(:,i_dim)),2); % Find all rows that have repeats with new indexes
        list{i_dim}(repeats,:) = []; % Remove repeats
    end
    
    % Remove overpriced blocks
    set_cost = sum(cost(list{i_dim}),2);
    overpriced = (set_cost > budget);
    list{i_dim}(overpriced,:) = [];
    set_cost(overpriced) = [];
end

% Concatenate blocks into single list, padding with zeros
list = cellfun(@(A)cat(2, A, zeros(size(A,1), block_size-size(A,2))), list, 'UniformOutput', false);
list = cat(1, list{:});
end
