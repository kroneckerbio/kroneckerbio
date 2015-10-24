function m = RemoveRule(m, name)
% RemoveRule Remove a rule by name

% Find added instances of this name
ind = strcmp(name, {m.Rules(1:m.nz).Name});
assert(sum(ind)==1, 'KroneckerBio:RemoveRule:RuleNotFound', 'Rule with name %s not found in model', name)

% Remove all mention of this rule
m.Rules(ind,:) = [];
m.nz = m.nz - 1;

m.Ready = false;
