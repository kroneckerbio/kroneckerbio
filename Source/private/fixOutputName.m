function name = fixOutputName(name)

name = row(name);
if isempty(name)
    name = '';
end

assert(ischar(name), 'KroneckerBio:Output:Name', 'The output name must be a string')

name = strtrim(name);
