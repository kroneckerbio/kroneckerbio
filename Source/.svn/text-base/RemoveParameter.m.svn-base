function m = RemoveParameter(m, name)

% Find added instances of this name
indAdd = strcmp(name, {m.add.Parameters.Name});

% Remove all mention of this parameter
m.Parameters(strcmp(name, {m.Parameters.Name})) = [];
m.add.Parameters(indAdd) = [];
m.add.nk = m.add.nk - nnz(indAdd);

m.Ready = false;