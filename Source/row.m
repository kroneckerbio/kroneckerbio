function a = row(a)
%ROW Row vectorize a matrix
%   This function converts a matrix into a row vector with the elements
%   taken columnwise from the matrix.
%
%   vector = row(matrix)
%
%   is short for the matlab operation
%
%   vector = matrix(:).'

% (c) 2012 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

a = a(:).';