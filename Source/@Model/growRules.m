function growRules(this)

current = numel(this.Rules);
add = this.nz - current;
if add > 0 || isempty(this.Rules)
    % Double length
    add = max(current,add);
    this.Rules = [this.Rules; struct('Name', cell(add,1), 'Expression', cell(add,1), 'Type', cell(add,1))];
end
