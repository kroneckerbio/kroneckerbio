function m = RemoveCompartment(m, name)
% RemoveCompartment Remove a compartment by name

% Find added instances of this name
ind = strcmp(name, {m.Compartments.Name});
assert(sum(ind)==1, 'KroneckerBio:RemoveCompartment:CompartmentNotFound', 'Compartment with name %s not found in model', name)

% Remove all mention of this compartment
m.Compartments(ind,:) = [];
m.nv = m.nv - 1;

m.Ready = false;