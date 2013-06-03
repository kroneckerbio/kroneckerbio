function compartments = growCompartments(compartments, nc)

% Add more room in vector if necessary
current = numel(compartments);
add = nc - current;
if add > 0 || isempty(compartments)
    % Double length
    add = max(current,add);
    compartments = [compartments; struct('Name', cell(add,1), 'Dimension', cell(add,1), 'Expressions', cell(add,1), 'Values', cell(add,1))];
end