function outputs = growOutputsAnalytic(outputs, ny)

if nargin < 2
    ny = 0;
    if nargin < 1
        outputs = [];
    end
end

current = numel(outputs);
add = ny - current;
if add > 0 || isempty(outputs)
    % Double length
    add = max(current,add);
    outputs = [outputs; struct('Name', cell(add,1), 'ID', cell(add,1), 'Expression', cell(add,1))];
end