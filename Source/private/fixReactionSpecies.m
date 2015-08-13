function names = fixReactionSpecies(names)
% Standardize a vector of species as a cell row of strings

if isempty(names) || ischar(names) && strcmp(names, '0')
    names = cell(1,0);
    return
end

if ischar(names)
    names = {names};
end

assert(iscellstr(names), 'KroneckerBio:fixReactionSpecies:names', 'names must be a cell array of strings or a string or empty')
names = row(names);
