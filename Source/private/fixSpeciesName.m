function name = fixSpeciesName(name)
% Standardize the species name

assert(ischar(name) && isempty(regexp(name, '\.|"', 'once')), ...
    'KroneckerBio:Species:Name', 'Species name must be a string without a dot or double quote')
