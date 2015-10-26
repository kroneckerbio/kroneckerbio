function ic = StateInitialConditionMassAction(ic)
% Allowed MassAction model state initial condition inputs:
% - empty cell(0,2) for no seed, default initial condition = 0
% - number for no seed, specified initial condition
% - string for seed, coefficient 1
% - cell array of seed name strings, all coefficients 1
% - nSeeds x 2 cell matrix of rows of [seed, coefficient pairs]

if isempty(ic)
    ic = cell(0,2);
elseif isnumeric(ic)
    assert(numel(ic) == 1, 'KroneckerBio:State:RepeatedStateSeed', 'States cannot have the same seed appear twice')
    ic = {'', ic};
elseif ischar(ic)
    ic = {row(ic), 1};
elseif iscellstr(ic)
    assert(numel(ic) == numel(unique(ic)), 'KroneckerBio:State:RepeatedStateSeed', 'States cannot have the same seed appears twice')
    ic = [vec(ic), num2cell(ones(numel(ic),1))];
else
    assert(size(ic,2) == 2, 'KroneckerBio:State:InvalidSeed', 'State has an invalid seed')
    assert(iscellstr(ic(:,1)), 'KroneckerBio:State:RepeatedStateSeed', 'States cannot have the same seed appear twice')
    assert(isnumeric([ic{:,2}]) && all([ic{:,2}] >= 0), 'KroneckerBio:State:NegativeSeed', 'States must have nonegative contributions from their seeds')
end
