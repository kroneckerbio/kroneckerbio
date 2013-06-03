function S = entropy(p, dim)
%ENTROPY Entropy of a discrete probability density.
%
%   S = entropy(p)
%   This function computes the entropy S associated with vector p. Vector p
%   must be normalized; that is, its entries must sum to 1.
%
%   S = entropy(p, dim)
%   Compute the entropy along dimensional dim. If dim is excluded or empty,
%   the first non-singular dimension is used
%   
%   Mathematically: S = -sum(p .* log(p))

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if nargin < 2
    dim = [];
end

% Compute entropy
S = -nansum(p .* log(p), dim);