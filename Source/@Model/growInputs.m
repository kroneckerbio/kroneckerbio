function growInputs(this)

% Add more room in vector if necessary
current = numel(this.Inputs);
add = this.nu - current;
if add > 0 || isempty(this.Inputs)
    % Double length
    add = max(current,add);
    this.Inputs = [this.Inputs; struct('Name', cell(add,1), 'Compartment', cell(add,1), 'DefaultValue', cell(add,1))];
end
