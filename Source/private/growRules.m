function rules = growRules(rules, nz)

if nargin < 2
    nz = 0;
    if nargin < 1
        rules = [];
    end
end

current = numel(rules);
add = nz - current;
if add > 0 || isempty(rules)
    % Double length
    add = max(current,add);
    rules = [rules; struct('Name', cell(add,1), 'ID', cell(add,1), 'Expression', cell(add,1))];
end
