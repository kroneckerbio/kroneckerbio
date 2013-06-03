function m = RemoveReaction(m, name)

m.Reactions(strcmp(name, {m.Reactions.Name})) = [];
for ir = 1:m.add.nr
    m.add.Reactions(ir).Parameters(strcmp(name, m.add.Reactions(ir).Names)) = [];
end

m.Ready = false;