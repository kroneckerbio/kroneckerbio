function m = RemoveCompartment(m, name)
% RemoveCompartment Remove a compartment by name

% Find added instances of this name
ind_main = strcmp(name, {m.Compartments.Name});
ind_add = strcmp(name, {m.add.Compartments.Name});

% Remove all mention of this compartment
m.Compartments(ind_main,:) = [];
m.add.Compartments(ind_add,:) = [];
m.add.nv = m.add.nv - nnz(ind_add);

m.Ready = false;
