function m = RemoveInput(m, name)
% Find last species by this name
indMain = strcmp(name, {m.Species.Name});
indAdd = strcmp(name, {m.add.Species.Name});
indLast = find(indAdd, 'last');

% Assert that the last species by this name is an input
assert((~isempty(indLast) && m.add.Species(indLast).IsInput) || (isempty(indLast) && ~isempty(indMain) && m.Species(indMain).IsInput) ||(isempty(indLast) && ~isempty(indMain)), 'KroneckerBio:RemoveInput:NotAnInput', 'Species %s is not an input', name)

% Remove all mention of this species
m.Species(indMain) = [];
m.add.Inputs(indAdd) = [];
m.add.nxu = m.add.nxu - nnz(indAdd);

m.Ready = false;