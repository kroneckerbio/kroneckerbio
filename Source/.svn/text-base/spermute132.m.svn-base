function A = spermute132(A, dim, final)
%SPERMUTE132 Permute the second and third axis of a three-dimensional array
%   stored as a two-dimenionsional sparse matrix
%
%   A = spermute132(A, dim, final)
%
%   Input dim is a 3 element vector indicating the size of the stored
%   array. Input final is a 2 element vector indicating the size of A as it
%   will be returned.
%
%   Matlab cannot store sparse arrays greater than 2 dimensions. Therefore,
%   clever methods must be used to store the n-dimensional arrays in
%   2-dimensional ones. Because n-dimensional arrays cannot be made, it is
%   not possible to use permute() to permute axes. This function works
%   around that by making clever use of linear indexing.

% (c) 2011 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Extract dimensions
d1 = dim(1);
d2 = dim(2);
d3 = dim(3);

% Put the swapping columns together
A = reshape(A, d1,d2*d3); % 1_23

% Permute
A = A(:, reshape(1:(d2*d3), d2,d3)'); % 1_32

% Reshape to desired dimensions
A = reshape(A, final);