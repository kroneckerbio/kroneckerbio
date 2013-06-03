function name = fixCompartmentName(name)
% Standardize the compartment name

if isempty(name)
    name = '';
else
    name = row(name);
end

assert(ischar(name) && isempty(regexp(name, '\.|,', 'once')), ...
    'KroneckerBio:Compartment:Name', 'Compartment name must be a string without "." or ","')
