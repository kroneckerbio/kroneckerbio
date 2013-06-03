function outputs = growOutputs(outputs, ny)

% Add more room in vector if necessary
current = numel(outputs);
add = ny - current;
if add > 0 || isempty(outputs)
    % Double length
    add = max(current,add);
    outputs = [outputs; struct('Name', cell(add,1), 'Expressions', cell(add,1), 'Values', cell(add,1))];
end