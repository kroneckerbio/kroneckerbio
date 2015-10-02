function m = RemoveOutput(m, name)
% RemoveOutput Remove an output by name

% Find added instances of this name
ind = strcmp(name, {m.Outputs.Name});
assert(sum(ind)==1, 'KroneckerBio:RemoveOutput:OutputNotFound', 'Output with name %s not found in model', name)

% Remove all mention of this output
m.Outputs(ind,:) = [];
m.ny = m.ny - 1;

m.Ready = false;
