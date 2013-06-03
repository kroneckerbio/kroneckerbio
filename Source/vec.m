function a = vec(a)
%VEC Vectorize a matrix
%   This function converts a matrix into a column vector with the elements
%   taken columnwise from the matrix.
%
%   vector = vec(matrix)
%
%   is short for the matlab operation
%
%   vector = matrix(:)

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

a = a(:);