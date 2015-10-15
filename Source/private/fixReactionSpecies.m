function names = fixReactionSpecies(names)
% Standardize a vector of species as a cell row of strings
%
% Input:
%   names [ string | cell array of strings ]
%       Reactant or product species names in various forms. Empty species can be
%       represented by anything with isempty == true. A single species can
%       be represented by a string or a cell array with 1 string element. A cell
%       array of species must contain all valid species.
% Output:
%   names: [ 1 x nSpecies cell vector of strings ]
%       Standardized list of species

% % Make sure cell array contains only species, no blanks.
% if iscell(names)
%     
% end

% Set default empty species form
if isempty(names)
    names = cell(1,0);
    return
end

% Wrap single species string in cell array
if ischar(names)
    names = {names};
end

% Check for valid names in cell array
assert(~any(cellfun(@isempty, names)), 'KroneckerBio:fixReactionSpecies:InvalidBlankName', 'names may not mix strings and blanks');
assert(iscellstr(names), 'KroneckerBio:fixReactionSpecies:InvalidName', 'names must be a cell array of strings or a string or empty')
names = row(names);

% Basic checks on species names:
%   Have 0 (unqualified) or 1 (qualified compartment.species) dot
%   Can't contain double quote "
for i = 1:length(names)
    name_i = names{i};
    assert(sum(ismember(name_i, '.')) < 2, 'KroneckerBio:fixReactionSpecies:InvalidNameDots', 'species name %s must have at most 1 dot, indicating a qualified compartment.species', name_i)
    assert(~any(ismember(name_i, '"')), 'KroneckerBio:fixReactionSpecies:InvalidNameDoubleQuotes', 'species name %s must not contain double quotes', name_i)
end