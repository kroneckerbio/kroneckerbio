function name = ModelName(name)

name = row(name);
if isempty(name)
    name = '';
end

assert(ischar(name), 'KroneckerBio:FieldValidator:ModelName', 'The model name must be a string')

name = strtrim(name);
