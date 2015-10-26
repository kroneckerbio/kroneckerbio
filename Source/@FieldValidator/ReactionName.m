function [name1, name2] = ReactionName(name)
% Standardize reaction name. Reaction names don't have to be unique so the
% forward name will be identical to the reverse if just 1 is provided. The
% default name if nothing is provided is an empty string.
%
% Inputs:
%   name [ string | 1 x 2 cell array of strings ]
%       Reaction name string or different forward and reverse reaction names in cell
%       array
% Outputs:
%   name1 [ string {''} ]
%       Standardized forward reaction name.
%   name2 [ string {''} ]
%       Standardized reaction reaction name

if isempty(name)
    name = '';
end

if ischar(name)
    name1 = name;
    name2 = name;
end

if iscell(name)
    if numel(name)
        name1 = name{1};
        name2 = name{1};
    elseif numel(name)
        name1 = name{1};
        name2 = name{2};
    else
        error('KroneckerBio:FieldValidator:ReactionName:TooManyNames', 'The name of reaction %s was given as a cell array, but there are more than two entries, the maximum', name)
    end
end