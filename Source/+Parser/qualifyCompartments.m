function [speciesList, qualified] = qualifyCompartments(speciesList, xu_full_names, compartment)
% Qualify input species string with compartment according to compartment
% qualification rules
%
% Inputs:
%   speciesList [ cell vector of strings ]
%       Species strings, with or without compartment, which will be
%       converted to qualified form
%   xu_full_names [ cell vector of strings ]
%       Names of all compartment.species in model
%   compartment [ string {''} ]
%       Default compartment (reaction compartment)
%
% Outputs:
%   speciesList [ cell vector of strings ]
%       Validated species strings qualified with compartments
%   qualified [ logical vector ]
%       Whether each member of input speciesList came with a qualified
%       compartment

import Parser.*

if nargin < 3
    compartment = [];
end

speciesList = vec(speciesList);
nSpecies = length(speciesList);
qualified = false(nSpecies,1);
for i = 1:nSpecies
    speciesOld = speciesList{i};
    speciesNew = qualifyCompartment(speciesOld, xu_full_names, compartment);
    speciesList{i} = speciesNew;
    
    if strcmp(speciesOld, speciesNew)
        qualified(i) = true;
    end
end