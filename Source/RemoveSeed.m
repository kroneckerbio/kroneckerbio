function m = RemoveSeed(m, name)
% RemoveSeed Remove a seed by name

% Find added instances of this name
ind = strcmp(name, {m.Seeds.Name});
assert(sum(ind)==1, 'KroneckerBio:RemoveSeed:SeedNotFound', 'Seed with name %s not found in model', name)

% Remove all mention of this seed
m.Seeds(ind,:) = [];
m.ns = m.ns - 1;

m.Ready = false;