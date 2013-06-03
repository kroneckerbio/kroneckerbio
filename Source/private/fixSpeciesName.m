function name = fixSpeciesName(name)
% Standardize the compartment name

if isempty(name)
    name = '';
else
    name = row(name);
end

assert(ischar(name) && isempty(regexp(name, '\.|,', 'once')), ...
    'KroneckerBio:Species:Name', 'Species name must be a string without "." or ","')
