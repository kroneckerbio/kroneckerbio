function name = ParameterName(name)
% Standardize the parameter name

if isempty(name)
    name = '';
else
    name = row(name);
end

assert(ischar(name) && isempty(regexp(name, '\.|"', 'once')), ...
    'KroneckerBio:FieldValidator:ParameterName', 'Parameter name must be a string without a dot or double quote')
