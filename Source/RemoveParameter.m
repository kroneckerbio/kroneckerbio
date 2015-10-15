function m = RemoveParameter(m, name)
% RemoveParameter Remove a parameter by name

% Find added instances of this name
ind = strcmp(name, {m.Parameters(1:m.nk).Name});
assert(sum(ind)==1, 'KroneckerBio:RemoveParameter:ParameterNotFound', 'Parameter with name %s not found in model', name)

% Remove all mention of this parameter
m.Parameters(ind,:) = [];
m.nk = m.nk - 1;

m.Ready = false;
