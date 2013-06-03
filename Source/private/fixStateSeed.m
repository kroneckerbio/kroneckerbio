function seed = fixStateSeed(seed)

if isempty(seed)
    seed = cell(0,2);
elseif isnumeric(seed)
    assert(numel(seed) == 1, 'KroneckerBio:State:RepeatedStateSeed', 'States cannot have the same seed appears twice')
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
