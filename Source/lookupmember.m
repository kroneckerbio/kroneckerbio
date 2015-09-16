function ind = lookupmember(x, d)
%LOOKUPMEMBER Map a vector onto a dictionary
%
%   I = LOOKUPMEMBER(X, D)
%
%   This function compares each element in matrix X to vector D and returns
%   the linear index in D corresponding to the location of the matching
%   value. The indexes are returned in matrix I, which has the same
%   dimensions as X. If D has multiple values that match, the index of the
%   first is returned. If D does not have a matching value, then an index
%   of 0 is returned.

if iscellstr(x) || iscellstr(d)
    if ischar(x)
        % Standardize x as cell array of strings
        x = {x};
    end
    if ischar(d)
        % Standardize d as cell array of strings
        d = {d};
    end
    
    % Use strcmp
    ind = zeros(size(x));
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
elseif isstruct(x) && isstruct(d)
    % Use isequal
    ind = zeros(size(x));
    for i = 1:numel(x)
        for j = 1:numel(d)
            % Find the index of the first matching value
            tmp = isequal(x(i), d(j));
            
            % If the value does not exist, return a 0 for the index
            if tmp
                ind(i) = j;
                break
            end
        end
    end
else
    % Use eq
    ind = zeros(size(x));
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