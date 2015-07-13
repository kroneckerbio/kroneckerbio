function name = fixReactionSpecies(name, compartment)

% Empty species is empty string
if isempty(name) || strcmp(name, '0')
    name = '';
    return
end

if nargin == 1
    if isempty(name)
        name = '';
    end
else
    if any(name == '.')
        % Name is complete, do not modify
    else
        % Add the compartment to the name
        name = [compartment '.' name];
    end
end
