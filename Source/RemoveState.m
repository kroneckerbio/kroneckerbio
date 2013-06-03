function m = RemoveState(m, name)

% Find states by this name
ind_main = strcmp(name, {m.States.Name});
ind_add = strcmp(name, {m.add.States.Name});

% Remove all mention of this state
m.States(ind_main) = [];
m.add.States(ind_add) = [];
m.add.nx = m.add.ns - nnz(ind_add);

m.Ready = false;
