function m = RemoveInput(m, name)

% Find inputs by this name
ind_main = strcmp(name, {m.Inputs.Name});
ind_add = strcmp(name, {m.add.Inputs.Name});

% Remove all mention of this input
m.Inputs(ind_main) = [];
m.add.Inputs(ind_add) = [];
m.add.nu = m.add.nu - nnz(ind_add);

m.Ready = false;
