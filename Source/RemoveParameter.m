function m = RemoveParameter(m, name)
% RemoveParameter Remove a parameter by name

% Find added instances of this name
ind_main = strcmp(name, {m.Parameters.Name});
ind_add = strcmp(name, {m.add.Parameters.Name});

% Remove all mention of this parameter
m.Parameters(ind_main,:) = [];
m.add.Parameters(ind_add,:) = [];
m.add.nk = m.add.nk - nnz(ind_add);

m.Ready = false;
