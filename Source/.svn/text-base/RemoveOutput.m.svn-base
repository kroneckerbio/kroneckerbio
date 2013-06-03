function m = RemoveOutput(m, name)

% Find added instances of this name
indAdd = strcmp(name, {m.add.Outputs.Name});

% Remove all mention of this output
m.Outputs(strcmp(name, {m.Outputs.Name})) = [];
m.add.Outputs(indAdd) = [];
m.add.ny = m.add.ny - nnz(indAdd);

m.Ready = false;