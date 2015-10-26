function species = getSpeciesFromExpr(expr, xu_full_names)
% Extract species from expression. Needed to get species that appear in 
%   rate expression but aren't reactants or products Looks for unqualified 
%   species and qualified compartment.species.
%
% Inputs:
%   expr [ string ]
%       Arbitrary expression string
%   xu_full_names [ cell vector of strings ]
%       Names of all compartment.species in model
%
% Outputs:
%   species [ 1 x n cell vector of strings ]
%       List of species found in expr. Can be unqualified species or qualified
%       compartment.species
%
% Note: Names in expr must be double-quoted if they contain invalid
%   characters. Note: this implies qualified compartment.species must
%   always be quoted because the dot is an invalid character.

import Parser.*

% Make species lookup list, quoting things that contain invalid characters
xu_names = getUnqualifiedSpeciesNames(xu_full_names);
all_names = [xu_names; xu_full_names];
for i = 1:length(all_names)
    names_i = all_names{i};
    if regexp(names_i, '[^\w.]')
        all_names{i} = ['"' names_i '"'];
    end
end

% Tokenize expression into potentially substitutable parts and see if
%   they're valid species names
species = {};
parts = regexp(expr, '[\w\.]+|"[^"]+"', 'match');
for iPart = 1:length(parts)
    part = parts{iPart};
    if any(ismember(all_names, part))
        part = strrep(part, '"', ''); % strip double-quotes because parts will go into species lists, which don't require them
        species = [species, part];
    end
end
end