function inputs = growInputs(inputs, nu)

% Add more room in vector if necessary
current = numel(inputs);
add = nu - current;
if add > 0 || isempty(inputs)
    % Double length
    add = max(current,add);
    inputs = [inputs; struct('Name', cell(add,1), 'Compartment', cell(add,1), 'ValueFunction', cell(add,1), 'Parameters', cell(add,1), 'Units', cell(add,1))];
end