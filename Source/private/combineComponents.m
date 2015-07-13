function combinedCompnents = combineComponents(oldComponents, newComponents, verbosity)
% Combine existing oldComponents with new newComponents. Elements in
% newComponents that share an ID with elements in oldComponents are discarded.
% Inputs:
%   oldComponents [ Components struct vector ]
%   newComponents [ Components struct vector ]
%   verbosity [ scalar integer ]
% Outputs:
%   combinedComponents [ Components struct vector ]

% Clean up empty spaces in components first
oldComponents(cellfun('isempty',{oldComponents.ID})) = [];
newComponents(cellfun('isempty',{newComponents.ID})) = [];

oldNames = {oldComponents.Name};
newNames = {newComponents.Name};

newMask = ismember(newNames, oldNames); % already exists in oldComponents
newInds = find(newMask);

if any(newMask) && verbosity > 0
    for i = 1:sum(newMask)
        warning('Component %s already exists in model, ignoring', newNames{newInds(i)})
    end
end

newComponents = newComponents(~newMask);
combinedCompnents = [oldComponents; newComponents];