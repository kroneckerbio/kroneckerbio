function m = RemoveReaction(m, name)
% RemoveReaction Remove a reaction by name. Since reactions don't have to have
% names, this obviously only works when you give them useful names.

% Find added instances of this name
ind = strcmp(name, {m.Reactions.Name});
assert(sum(ind)==1, 'KroneckerBio:RemoveReaction:ReactionNotFound', 'Reaction with name %s not found in model', name)

% Remove all mention of this seed
m.Reactions(ind,:) = [];
m.nr = m.nr - 1;

m.Ready = false;
