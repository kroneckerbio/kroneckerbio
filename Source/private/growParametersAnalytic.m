function parameters = growParametersAnalytic(parameters, nk)

if nargin < 2
    nk = 0;
    if nargin < 1
        parameters = [];
    end
end

current = numel(parameters);
add = nk - current;
if add > 0 || isempty(parameters)
    % Double length
    add = max(current,add);
    parameters = [parameters; struct('Name', cell(add,1), 'ID', cell(add,1), 'Value', zeros(1))];
end