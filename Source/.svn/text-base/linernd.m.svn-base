function r = linernd(p1, p2, n)
%LINERND choses a random point along the line segment connecting two points
%
%   r = linernd(p1, p2, n)
%
%   Each column in matrix r is a random point along the the line segment
%   connecting the two points p1 and p2. The first dimension of r is the
%   same length as p1 and p2. There are n columns in r, each of which is an
%   independent draw.
%
%   r = linernd(p1, p2) is the same as r = linernd(p1, p2, 1)
%   r = linernd(p1) is the same as r = linernd(p1, zeros(size(p1)))

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

%% Work-up
% Clean-up inputs
if nargin < 3
    n = [];
    if nargin < 2
        p2 = [];
    end
end

% Get length of vector
dim = numel(p1);

% If p2 is missing, then it is the origin
if isempty(p2)
    p2 = zeros(dim,1);
end

% If n is missing then it is 1
if isempty(n)
    n = 1;
end

% Vectorize points
p1 = p1(:);
p2 = p2(:);

% Get length of vector
dim = numel(p1);

%% Random number generator
r = zeros(dim,n);
for i = 1:n
    r(:,i) = (p2 - p1) * rand + p1;
end
