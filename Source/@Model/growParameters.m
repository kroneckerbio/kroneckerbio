function growParameters(this)

% Add more room in vector if necessary
current = numel(this.Parameters);
add = this.nk - current;
if add > 0 || isempty(this.Parameters)
    % Double length
    add = max(current,add);
    this.Parameters = [this.Parameters; struct('Name', cell(add,1), 'Value', cell(add,1))];
end
