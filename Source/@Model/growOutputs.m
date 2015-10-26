function growOutputs(this)

% Add more room in vector if necessary
current = numel(this.Outputs);
add = this.ny - current;
if add > 0 || isempty(this.Outputs)
    % Double length
    add = max(current,add);
    this.Outputs = [this.Outputs; struct('Name', cell(add,1), 'Expression', cell(add,1))];
end
