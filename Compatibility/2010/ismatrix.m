function tf = ismatrix(A)
% ISMATRIX Returns true if size(A) returns [m, n] with nonzero m and n,
%   otherwise returns false

dims = size(A);

if numel(dims) == 2 && dims(1) > 0 && dims(2) > 0
    tf = true;
else
    tf = false;
end
