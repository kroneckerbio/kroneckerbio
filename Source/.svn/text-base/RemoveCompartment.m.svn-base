function m = RemoveCompartment(m, name)

% Find added instances of this name
indAdd = strcmp(name, {m.add.Compartments.Name});

% Remove all mention of this compartment
m.Compartments(strcmp(name, {m.Compartments.Name})) = [];
m.add.Compartments(indAdd) = [];
m.add.nv = m.add.nv - nnz(indAdd);

m.Ready = false;