function m = RemoveInput(m, name)
% RemoveInput Remove an input by name. Name can be an unqualified species name or
% a qualified compartment.species. An unqualified species will be removed if it
% is unique in the model.

if isFullSpeciesName(name)
    existing_full_names = strcat(vec({m.Inputs.Compartment}), '.', vec({m.Inputs.Name}));
    ind = strcmp(name, existing_full_names);

    assert(nnz(ind) ~= 0, 'KroneckerBio:RemoveInput:InputNotFound', 'Input with name %s not found in model', name)
else
    ind = strcmp(name, {m.Inputs.Name});
    
    assert(nnz(ind) ~= 0, 'KroneckerBio:RemoveInput:InputNotFound', 'Input with name %s not found in model', name)
    assert(nnz(ind) == 1, 'KroneckerBio:RemoveInput:AmbiguousInputRemoval', 'Multiple inputs were found with the name "%s", use the full name ("compartment.%s") to remove it', name, name)
end

m.Inputs(ind,:) = [];
m.nu = m.nu - 1;

m.Ready = false;