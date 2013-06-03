function m = RemoveSpecies(m, name)

% Find added instances of this name
indAdd = strcmp(name, {m.add.Species.Name});

% Remove all mention of this species
m.Species(strcmp(name, {m.Species.Name})) = [];
m.add.Species(indAdd) = [];
m.add.nxu = m.add.nxu - nnz(indAdd);

m.Ready = false;