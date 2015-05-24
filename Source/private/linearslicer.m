function sliced = linearslicer(dim, varargin)
%LINEARSLICER Convert several slice indexes along different dimensions into
%   a slice index assuming all the dimensions are now stacked
%
%   sliced = linearslicer(dim, ...)
%
%   Inputs
%   dim: [ nonnegative integer vector ]
%       The size of each dimension
%   ...: [ index vector ]
%       The number of these must be equal to the number of dimensions
%       according to dim and each must be a legal index into a dimension of
%       the corresponding size
%
%   Outputs
%   sliced: [ linear index vector ]
%       The linear index of every element that would be selected by the
%       dimensional indexes supplied

% (c) 2015 David R Hagen
% This work is released under the MIT license.

box = reshape(1:prod(dim), dim);
sliced = vec(box(varargin{:}));
