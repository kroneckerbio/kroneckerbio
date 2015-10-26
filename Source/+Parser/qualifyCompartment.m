function species = qualifyCompartment(species, xu_full_names, compartment)
% Qualify input species string with compartment according to compartment
% qualification rules
%
% Inputs:
%   species [ string ]
%       Single species string, with or without compartment, which will be
%       converted to qualified form
%   xu_full_names [ cell vector of strings ]
%       Names of all compartment.species in model
%   compartment [ string {''} ]
%       Default compartment (reaction compartment)
%
% Outputs:
%   species [ string ]
%       Validated species string qualified with compartment

import Parser.*

if nargin < 3
    compartment = [];
end

xu_names = getUnqualifiedSpeciesNames(xu_full_names);

if ismember('.', species) % qualified - make sure this species exists
    
    assert(ismember(species, xu_full_names), 'KroneckerBio:Parser:qualifyCompartments:MissingQualifiedSpeciesName', 'The qualified name %s does not exist in the model', species)
    
else % unqualified - apply default compartment rules
    
    assert(ismember(species, xu_names), 'KroneckerBio:Parser:qualifyCompartments:MissingUnqualifiedSpeciesName', 'The unqualified name %s does not exist in the model', species)
    
    if isempty(compartment) % no reaction compartment
        % Make sure species is unique - if not, it's ambiguous and throws an error
        speciesPos = ismember(xu_names, species);
        assert(sum(speciesPos) == 1, 'KroneckerBio:Parser:qualifyCompartments:AmbiguousSpeciesName', 'The species name %s appears in multiple compartments and no default compartment is specified', species)
        species = xu_full_names{speciesPos};
    else % reaction compartment present
        species = [compartment '.' species];
        assert(ismember(species, xu_full_names), 'KroneckerBio:Parser:qualifyCompartments:MissingSpeciesInReactionCompartment', 'The species %s not found in compartment %s', species, compartment)
    end
    
end
