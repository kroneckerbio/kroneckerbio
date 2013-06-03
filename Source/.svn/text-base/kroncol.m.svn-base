function val = kroncol(v1, v2)
% KRONCOL Fast kronecker product of two column vectors
%
%   val = kroncol(v1, v2)
%
%   If v1 and v2 are column vectors, val is exactly the same as if 
%   val = kron(v1, v2) had been called. However, kron is not optimized in
%   Matlab. This function uses a trick with the matrix product (which is
%   heavily optimized) to arrive at the exact same answer.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

val = v2 * v1.';
val = val(:);