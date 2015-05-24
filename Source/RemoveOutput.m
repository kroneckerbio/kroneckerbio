function m = RemoveOutput(m, name)
% RemoveOutput Remove an output by name

% Find added instances of this name
ind_main = strcmp(name, {m.Outputs.Name});
ind_add = strcmp(name, {m.add.Outputs.Name});

% Remove all mention of this output
m.Outputs(ind_main,:) = [];
m.add.Outputs(ind_add,:) = [];
m.add.ny = m.add.ny - nnz(ind_add);

m.Ready = false;
