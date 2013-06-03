function out = pastestruct(old, new)
%PASTESTRUCT Paste an new structure over an old, tossing new fields
%
%   out = pastestruct(old, new)
%
%   This returns a structure that has the dimensions of new, the fields of
%   old, and the values of new. When a field is unique to old, it is
%   retained. When a field is unique to new, it is discarded.
%
%   This function is useful for standardizing a structure to a given
%   template so that it can be concatenated with similar structures.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Template is old
out = repmat(old, size(new));

% Fetch the fields
fields = fieldnames(old);

% Loop over all fields, copy if they exist in new
for i = 1:length(fields)
    if isfield(new, fields{i})
        [out.(fields{i})] = deal(new.(fields{i}));
    end
end

