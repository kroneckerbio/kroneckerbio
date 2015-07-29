function seed = fixStateSeed(seed)
% Allowed MassAction model seed inputs:
% - empty cell(0,2) for no seed, default initial condition = 0
% - number for no seed, specified initial condition
% - string for seed, coefficient 1
% - cell array of seed name strings, all coefficients 1
% - nSeeds x 2 cell matrix of rows of [seed, coefficient pairs]

if isempty(seed)
    seed = cell(0,2);
elseif isnumeric(seed)
    assert(numel(seed) == 1, 'KroneckerBio:State:RepeatedStateSeed', 'States cannot have the same seed appear twice')
    seed = {'', seed};
elseif ischar(seed)
    seed = {row(seed), 1};
elseif iscellstr(seed)
    assert(numel(seed) == numel(unique(seed)), 'KroneckerBio:State:RepeatedStateSeed', 'States cannot have the same seed appears twice')
    seed = [vec(seed), num2cell(ones(numel(seed),1))];
else
    assert(size(seed,2) == 2, 'KroneckerBio:State:InvalidSeed', 'State has an invalid seed')
    assert(iscellstr(seed(:,1)), 'KroneckerBio:State:RepeatedStateSeed', 'States cannot have the same seed appear twice')
    assert(isnumeric([seed{:,2}]) && all([seed{:,2}] >= 0), 'KroneckerBio:State:NegativeSeed', 'States must have nonegative contributions from their seeds')
end
