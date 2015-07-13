function inputs = growInputsAnalytic(inputs, nu)

if nargin < 2
    nu = 0;
    if nargin < 1
        inputs = [];
    end
end

current = numel(inputs);
add = nu - current;
if add > 0 || isempty(inputs)
    % Double length
    add = max(current,add);
    inputs = [inputs; struct('Name', cell(add,1), 'ID', cell(add,1), 'Compartment', cell(add,1), 'DefaultValue', zeros(1))];
end