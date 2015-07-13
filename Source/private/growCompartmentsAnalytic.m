function compartments = growCompartmentsAnalytic(compartments, nc)

if nargin < 2
    nc = 0;
    if nargin < 1
        compartments = [];
    end
end

current = numel(compartments);
add = nc - current;
if add > 0 || isempty(compartments)
    % Double length
    add = max(current,add);
    compartments = [compartments; struct('Name', cell(add,1), 'ID', cell(add,1), 'Dimension', zeros(1), 'Size', zeros(1))];
end