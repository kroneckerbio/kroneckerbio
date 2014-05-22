function m = RemoveState(m, name)

if isFullSpeciesName(name)
    existing_full_names = strcat(vec({m.States.Compartment}), '.', vec({m.States.Name}));
    new_full_names = strcat(vec({m.add.States.Compartment}), '.', vec({m.add.States.Name}));

    ind_main = strcmp(name, existing_full_names);
    ind_add = strcmp(name, new_full_names);

    assert(nnz(ind_main) + nnz(ind_add) ~= 0, 'KroneckerBio:RemoveState:NoSuchState', 'No state was found with the name "%s"', name)

    m.States(ind_main,:) = [];
    m.add.States(ind_add,:) = [];
    m.add.nx = m.add.nx - nnz(ind_add);
else
    % Find inputs by this name
    ind_main = strcmp(name, {m.States.Name});
    ind_add = strcmp(name, {m.add.States.Name});
    
    assert(nnz(ind_main) + nnz(ind_add) ~= 0, 'KroneckerBio:RemoveState:NoSuchState', 'No state was found with the name "%s"', name)
    assert(nnz(ind_main) + nnz(ind_add) == 1, 'KroneckerBio:RemoveState:AmbiguousRemoval', 'Multiple states were found with the name "%s", use the full name ("compartment.%s") to remove it', name, name)
    
    % Remove all mention of this input
    m.States(ind_main,:) = [];
    m.add.States(ind_add,:) = [];
    m.add.nx = m.add.nx - nnz(ind_add);
end

m.Ready = false;
