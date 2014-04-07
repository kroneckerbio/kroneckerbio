function name = fixReactionCompartment(name)
% Standardize compartment names as a cell array of strings

if isempty(name)
    name = cell(0,1);
elseif ischar(name)
    name = {name};
elseif iscellstr(name)
    name = vec(name);
else
    error('KroneckerBio:Reaction:InvalidPossibleCompartments', 'The possible compartments for a reaction must be a string or cell array of strings')
end
