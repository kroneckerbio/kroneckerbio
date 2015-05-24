function m = RemoveReaction(m, name)
% RemoveReaction Remove a reaction by name

% Find added instances of this name
ind_main = strcmp(name, {m.Reactions.Name});
ind_add = strcmp(name, {m.add.Reactions.Name});

% Remove all mention of this seed
m.Reactions(ind_main,:) = [];
m.add.Reactions(ind_add,:) = [];
m.add.nr = m.add.nr - nnz(ind_add);

m.Ready = false;
