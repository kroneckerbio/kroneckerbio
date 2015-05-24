function m = RemoveInput(m, name)
% RemoveInput Remove an input by name

if isFullSpeciesName(name)
    existing_full_names = strcat(vec({m.Inputs.Compartment}), '.', vec({m.Inputs.Name}));
    new_full_names = strcat(vec({m.add.Inputs.Compartment}), '.', vec({m.add.Inputs.Name}));

    ind_main = strcmp(name, existing_full_names);
    ind_add = strcmp(name, new_full_names);

    assert(nnz(ind_main) + nnz(ind_add) ~= 0, 'KroneckerBio:RemoveInput:NoSuchInput', 'No input was found with the name "%s"', name)

    m.Inputs(ind_main,:) = [];
    m.add.Inputs(ind_add,:) = [];
    m.add.nu = m.add.nu - nnz(ind_add);
else
    % Find inputs by this name
    ind_main = strcmp(name, {m.Inputs.Name});
    ind_add = strcmp(name, {m.add.Inputs.Name});
    
    assert(nnz(ind_main) + nnz(ind_add) ~= 0, 'KroneckerBio:RemoveInput:NoSuchInput', 'No input was found with the name "%s"', name)
    assert(nnz(ind_main) + nnz(ind_add) == 1, 'KroneckerBio:RemoveInput:AmbiguousRemoval', 'Multiple inputs were found with the name "%s", use the full name ("compartment.%s") to remove it', name, name)
    
    % Remove all mention of this input
    m.Inputs(ind_main,:) = [];
    m.add.Inputs(ind_add,:) = [];
    m.add.nu = m.add.nu - nnz(ind_add);
end

m.Ready = false;
