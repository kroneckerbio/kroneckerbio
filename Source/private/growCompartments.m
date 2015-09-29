function compartments = growCompartments(compartments, nv)

if nargin < 2
    nv = 0;
    if nargin < 1
        compartments = [];
    end
end

% Add more room in vector if necessary
current = numel(compartments);
add = nv - current;
if add > 0 || isempty(compartments)
    % Double length
    add = max(current,add);
    compartments = [compartments; struct('Name', cell(add,1), 'Dimension', cell(add,1), 'Size', cell(add,1))];
end
