function m = RemoveState(m, name)
% RemoveState Remove a state by name. Name can be an unqualified species name or
% a qualified compartment.species. An unqualified species will be removed if it
% is unique in the model.

if isFullSpeciesName(name)
    existing_full_names = strcat(vec({m.States.Compartment}), '.', vec({m.States.Name}));
    ind = strcmp(name, existing_full_names);
    
    assert(nnz(ind) ~= 0, 'KroneckerBio:RemoveState:StateNotFound', 'State with name %s not found in model', name)
else
    ind = strcmp(name, {m.States.Name});
    
    assert(nnz(ind) ~= 0, 'KroneckerBio:RemoveState:StateNotFound', 'State with name %s not found in model', name)
    assert(nnz(ind) == 1, 'KroneckerBio:RemoveState:AmbiguousStateRemoval', 'Multiple states were found with the name "%s", use the full name ("compartment.%s") to remove it', name, name)
end

m.States(ind,:) = [];
m.nx = m.nx - 1;

m.Ready = false;