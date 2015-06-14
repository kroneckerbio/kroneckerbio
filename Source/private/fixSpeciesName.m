function name = fixSpeciesName(name)
% Standardize the species name

assert(ischar(name), ...
    'KroneckerBio:Species:Name', 'Species name must be a string')
