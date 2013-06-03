function xbar = weightedmean(x, w, dim)
%WEIGHTEDMEAN Weighted mean
%
%   mean = weightedmean(x, w, dim)
%
%   The mean of vector x is computed according to the weight w on each
%   point. The size of x and w must be the same. If dimension dim is
%   missing or empty, the mean is taken along the first non-singleton
%   dimension.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

if nargin < 3
    dim = find(size(x) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
end

xbar = sum(x .* w, dim) ./ sum(w, dim);