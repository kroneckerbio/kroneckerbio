function name = SpeciesName(name)
% Standardize the species name

assert(ischar(name) && isempty(regexp(name, '\.|"', 'once')), ...
    'KroneckerBio:FieldValidator:SpeciesName', 'Species name must be a string without a dot or double quote')
