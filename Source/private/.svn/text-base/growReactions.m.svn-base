function reactions = growReactions(reactions, nr)

% Add more room in vector if necessary
current = numel(reactions);
add = nr - current;
if add > 0 || isempty(reactions)
    % Double length
    add = max(current,add);
    reactions = [reactions; struct('Names', cell(add,1), 'Compartment', cell(add,1), 'Reactants', cell(add,1), 'Products', cell(add,1), 'Parameters', cell(add,1))];
end