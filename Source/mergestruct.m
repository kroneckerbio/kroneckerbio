function combostruct = mergestruct(oldstruct, newstruct)
%MERGESRUCT Merge one structure into another.
%   combostruct = mergestruct(oldstruct, newstruct) adds the fields and
%   values of newstruct into oldstruct and returns the merged structure.
%   For field names that are the same between oldstruct and newstruct, the
%   values in newstruct are returned.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
if nargin < 2
    newstruct = [];
end

if isempty(newstruct)
    newstruct = struct();
end

% Initialize combostruct
combostruct = oldstruct;

% Fetch fields
newfields = fieldnames(newstruct);
oldfields = fieldnames(oldstruct);

% Replace fields in combostruct with newstruct
for i = 1:length(newfields)
    % If fields are the same regardless of case, keep the old case and new
    % value.
    index = find(strcmpi(oldfields, newfields{i}), 1, 'first');
    if ~isempty(index)
        % Match was found
        combostruct.(oldfields{index}) = newstruct.(newfields{i});
    else
        % Is new field
        combostruct.(newfields{i}) = newstruct.(newfields{i});
    end
end