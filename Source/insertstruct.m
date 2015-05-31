function array = insertstruct(array, elements, varargin)
%insertstruct Put a struct array into a struct array, adding fields as
%   needed
%
%   new = insertstruct(old, elements, ...)
%
%   The old array is updated with elements at position ... to create a new
%   array. The size of elements must be assignable to old(...). If the
%   fields of old and elements do not match, then the missing fields are
%   added to old and elements (with a value of []) so that the elements can
%   be added. The order of the fields in new begin with the fields of old
%   and continue with the new fields of elements.

% (c) 2015 David R Hagen
% This work is released under the MIT license.

all_old_names = fieldnames(array);
all_new_names = fieldnames(elements);

% Determine names that need to be added
i_new_names = ~ismember(all_new_names, all_old_names);
new_names = all_new_names(i_new_names);

i_old_names = ~ismember(all_old_names, all_new_names);
old_names = all_old_names(i_old_names);

% Add each new name to array
for name = row(new_names)
    [array.(name{1})] = deal([]);
end

% Add each old name to elements
for name = row(old_names)
    [elements.(name{1})] = deal([]);
end
elements = orderfields(elements, array);

% Insert elements
array(varargin{:}) = elements;
