function xu_names = getUnqualifiedSpeciesNames(xu_full_names)
% Extract unqualified species from qualified compartment.species
%
% Inputs:
%   xu_full_names [ cell vector of strings ]
%       Names of all compartment.species in model
%
% Outputs:
%   xu_names [ cell vector of strings ]
%       Names of unqualified species in model. May contain duplicates

xu_full_names = vec(xu_full_names);
nxu = length(xu_full_names);
xu_names = cell(nxu,1);
for i = 1:nxu
    parts = strsplit(xu_full_names{i}, '.');
    assert(length(parts) == 2, 'KroneckerBio:Parser:getUnqualifiedSpeciesNames:InvalidSpecies', 'xu_full_names must contain all compartment.species')
    xu_names{i} = parts{2};
end