function growReactions(this)

% Add more room in vector if necessary
current = numel(this.Reactions);
add = this.nr - current;
if add > 0 || isempty(this.Reactions)
    % Double length
    add = max(current,add);
    this.Reactions = [this.Reactions; struct('Name', cell(add,1), 'Reactants', cell(add,1), 'Products', cell(add,1), 'Parameter', cell(add,1), 'Compartment', cell(add,1))];
end
