function varargout = dealfield(obj, field)
%DEALFIELD From a vector of objects or structs, deal the contents of a
%   field to the output
%
%   [varargout{1:n}] = dealfield(obj, field)
%
%   Normally, one can do [contents{1:n}] = deal(mystruct(1:n).field) to
%   dump the specified field of each member of the struct vector. Matlab
%   failed to make the obvious extension to Java objects. (Accessing the
%   field of vector of objects returns an opaque error.) 
%
%   This function extends this capability to objects by looping over the
%   vector for you. It will also work for accessing the field of a
%   structure by string name at runtime in case you didn't know that
%   mystruct.('myfield') would do that, too.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

n = numel(obj);

assert(nargout == n, 'dealfield:nargout', 'The number of outputs should match the number of inputs.')

for i = 1:n
    varargout{i} = obj(i).(field);
end