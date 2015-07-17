function fullName = fixSpeciesFullName(species, compartment, m)
% Get full compartment.species name for input species. Errors if species not
% found in compartment. Returns empty [] if species is missing or has
% an empty string for a name.
% Inputs:
%   species [ string | [] ]
%       Species name as 'species' short name or 'compartment.species' full name
%       or empty scalar double []
%   compartment [ string {1st component in model} ]
%       Compartment to search in if species doesn't contain one
%   m [ Model.Analytic struct ]
%       Model to search for species and compartments in
% Outputs:
%   fullName [ string | [] ]
%       Full species name in 'compartment.species' form if a species string is
%       entered or empty scalar double [] if no species is entered
% Note: make sure m.* and m.add.* are consistent - i.e., components only appear
% in 1 place when this is called.

%% Ignore if empty
if isempty(species) || (ischar(species) && strcmp(species, ''))
    fullName = [];
    return
end

%% Get all species and compartments
xuNames = [{m.States.Name}, {m.add.States.Name}, {m.Inputs.Name}, {m.add.Inputs.Name}];
xuNames = xuNames(~cellfun('isempty',xuNames));
xuCompartments = [{m.States.Compartment}, {m.add.States.Compartment}, {m.Inputs.Compartment}, {m.add.Inputs.Compartment}];
xuCompartments = xuCompartments(~cellfun('isempty',xuCompartments));

%% If name is already a full name, ignore compartment argument, check validity and return as is
if any(species == '.')
    nameParts = strsplit(species, '.');
    compartment = nameParts{1};
    name = nameParts{2};
    if any(ismember(xuCompartments, compartment) == ismember(xuNames, name))
        fullName = species;
        return
    else
        error('fixReactionSpeciesAnalytic: species %s not found in model', species)
    end
end

%% Species is a short name only - get compartment
% Check if species is unique - ignore default compartment and return with
% correct compartment
nameMask = ismember(xuNames, species);
if sum(nameMask) == 1 % unique
    fullName = [xuCompartments{nameMask} '.' species];
    return
end

% Get default (1st) compartment if not supplied
if isempty(compartment)
    vNames = [{m.Compartments.Name}, {m.add.Compartments.Name}];
    vNames = unique(vNames(~cellfun('isempty',vNames)));
    if isempty(vNames)
        error('fixReactionSpeciesAnalytic: model has no compartments for state to reside in')
    end
    compartment = vNames{1};
end

% Check if species exists in compartment
xuInds = ismember(xuCompartments, compartment);
if ~any(xuInds)
    error('fixReactionSpeciesAnalytic: model does not have compartment %s', compartment)
end
xuNames = xuNames(xuInds);
xuInd = ismember(xuNames, species);
switch sum(xuInd)
    case 1 % found
        fullName = [compartment '.' species];
        return
    case 0
        error('fixReactionSpeciesAnalytic: species %s not found in compartment %s', species, compartment)
    otherwise
        error('fixReactionSpeciesAnalytic: species %s found multiple times in compartment %s', species, compartment)
end
    

