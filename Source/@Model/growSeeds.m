function growSeeds(this)

% Add more room in vector if necessary
current = numel(this.Seeds);
add = this.ns - current;
if add > 0 || isempty(this.Seeds)
    % Double length
    add = max(current,add);
    this.Seeds = [this.Seeds; struct('Name', cell(add,1), 'Value', cell(add,1))];
end
