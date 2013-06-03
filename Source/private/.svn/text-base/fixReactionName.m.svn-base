function [name1 name2] = fixReactionName(name, kForward, kReverse)

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
        error('KroneckerBio:fixReactionName:TooManyNames', 'The name of reaction %s was given as a cell array, but there are more than two entries, the maximum', name)
    end
end