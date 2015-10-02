function name = fixParameterName(name)
% Standardize the parameter name

if isempty(name)
    name = '';
else
    name = row(name);
end

assert(ischar(name) && isempty(regexp(name, '\.|"', 'once')), ...
    'KroneckerBio:Parameter:Name', 'Parameter name must be a string without a dot or double quote')
