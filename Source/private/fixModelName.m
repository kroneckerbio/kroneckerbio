function name = fixModelName(name)

name = row(name);
if isempty(name)
    name = '';
end

assert(ischar(name), 'KroneckerBio:Model:Name', 'The model name must be a string')

name = strtrim(name);
