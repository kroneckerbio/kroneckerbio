function ind = lookup(x, d)
%LOOKUP Map a vector onto a dictionary
%
%   I = LOOKUP(X, D)
%
%   This function compares each element in matrix X to vector D and returns
%   the linear index in D corresponding to the location of the matching
%   value. The indexes are returned in matrix I, which has the same
%   dimensions as X. If D has multiple values that match, the index of the
%   first is returned. If D does not have a matching value, then an index
%   of 0 is returned.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

ind = zeros(size(x));

if iscellstr(x) || iscellstr(d)
    % Use strcmp
    for i = 1:numel(x)
        % Find the index of the first matching value
        tmp = find(strcmp(x(i), d), 1);
        
        % If the value does not exist, return a 0 for the index
        if isempty(tmp)
            tmp = 0;
        end
        
        % Store the final result
        ind(i) = tmp;
    end
else
    % Use eq
    for i = 1:numel(x)
        % Find the index of the first matching value
        tmp = find(x(i) == d, 1);
        
        % If the value does not exist, return a 0 for the index
        if isempty(tmp)
            tmp = 0;
        end
        
        % Store the final result
        ind(i) = tmp;
    end
end