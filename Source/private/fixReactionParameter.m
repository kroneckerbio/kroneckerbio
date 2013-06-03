function name = fixReactionParameter(name)

% Empty parameter is empty string
if ischar(name) || isempty(name)
    if isempty(name) || strcmp(name, '0')
        name = {'', 1};
    else
        name = {row(name), 1};
    end
else
    % Cell {name, modifier}
    assert(iscell(name) && numel(name) == 2 && ischar(name{1}) && isnumeric(name{2}), ...
        'KroneckerBio:Reaction:InvalidParameter', 'Invalid parameter provided to a reaction')
    
    if isempty(name{1}) || strcmp(name{1}, '0')
        name{1} = '';
    end
    
    name = row(name);
end
