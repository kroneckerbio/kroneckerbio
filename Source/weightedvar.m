function xvar = weightedvar(x, w, dim)
%WEIGHTEDVAR Weighted variance
%
%   var = weightedvar(x, w, dim)
%
%   The variance of the sample x is computed according the weight w on each
%   point. The size of x and w must be the same. If dimension dim is
%   missing or empty, the variance is computed along the first
%   non-singleton dimension.
%
%   This function computes the weighted variance differently from MATLAB's
%   var(), which, in the opinion of David R Hagen, computes it incorrectly.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if nargin < 3
    dim = find(size(x) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
end

e = x - weightedmean(x, w, dim);
w = w ./ sum(w, dim);

xvar = 1 ./ (1 - sum(w.^2, dim)) .* sum(w .* e.^2, dim);