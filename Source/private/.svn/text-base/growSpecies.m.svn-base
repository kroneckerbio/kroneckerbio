function states = growSpecies(states, nx)

% Add more room in vector if necessary
current = numel(states);
add = nx - current;
if add > 0 || isempty(states)
    % Double length
    add = max(current,add);
    states = [states; struct('Name', cell(add,1), 'Compartment', cell(add,1), 'IsInput', cell(add,1), 'Value', cell(add,1), 'Units', cell(add,1))];
end