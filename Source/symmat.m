function A = symmat(A)
%SYMMAT Symmetrize a matrix.
%
%   A = symmat(A) symmetrizes square matrix A. It is meant to replace code
%   like this:
%
%   A = (A + A.') / 2;

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

A = (A + A.') ./ 2;