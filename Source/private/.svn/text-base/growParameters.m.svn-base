function parameters = growParameters(parameters, nk)

% Add more room in vector if necessary
current = numel(parameters);
add = nk - current;
if add > 0 || isempty(parameters)
    % Double length
    add = max(current,add);
    parameters = [parameters; struct('Name', cell(add,1), 'Value', cell(add,1), 'Units', cell(add,1))];
end