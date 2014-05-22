function m = RemoveSeed(m, name)
% RemoveSeed Remove a seed by name

% Find added instances of this name
ind_main = strcmp(name, {m.Seeds.Name});
ind_add = strcmp(name, {m.add.Seeds.Name});

% Remove all mention of this seed
m.Seeds(ind_main,:) = [];
m.add.Seeds(ind_add,:) = [];
m.add.nk = m.add.nk - nnz(ind_add);

m.Ready = false;
