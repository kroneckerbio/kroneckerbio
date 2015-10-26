function growCompartments(this)

% Add more room in vector if necessary
current = numel(this.Compartments);
add = this.nv - current;
if add > 0 || isempty(this.Compartments)
    % Double length
    add = max(current, add);
    this.Compartments = [this.Compartments; struct('Name', cell(add,1), 'Dimension', cell(add,1), 'Size', cell(add,1))];
end
