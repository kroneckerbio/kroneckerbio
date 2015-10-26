function growStates(this)

% Add more room in vector if necessary
current = numel(this.States);
add = this.nx - current;
if add > 0 || isempty(this.States)
    % Double length
    add = max(current,add);
    this.States = [this.States; struct('Name', cell(add,1), 'Compartment', cell(add,1), 'InitialValue', cell(add,1))];
end
