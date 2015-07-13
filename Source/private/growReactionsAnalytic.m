function reactions = growReactionsAnalytic(reactions, nr)

if nargin < 2
    nr = 0;
    if nargin < 1
        reactions = [];
    end
end

current = numel(reactions);
add = nr - current;
if add > 0 || isempty(reactions)
    % Double length
    add = max(current,add);
    reactions = [reactions; struct('Name', cell(add,1), 'ID', cell(add,1), 'Reactants', cell(add,1), 'Products', cell(add,1), 'Rate', cell(add,1))];
end