function name = OutputName(name)

name = row(name);
if isempty(name)
    name = '';
end

assert(ischar(name), 'KroneckerBio:FieldValidator:OutputName', 'The output name must be a string')

name = strtrim(name);
