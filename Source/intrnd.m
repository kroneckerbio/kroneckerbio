function r = intrnd(a,b,varargin)
%INTRND Draw random integers
%
%   r = intrnd(a, b)
%   Draw a random integer over the range of a:b, inclusive.
%   
%   r = intrnd(a, b, m)
%   Draw a square matrix m by m independent random integers
%
%   r = intrnd(a, b, m, n)
%   r = intrnd(a, b, [m,n])
%   Draw m by n independent random integers
%
%   r = intrnd(a, b, m,...)
%   r = intrnd(a, b, [m,...])
%   Draw an arbitrary dimensional array of independent random integers
%
%   This function is deprecated and will be removed in future versions. Use
%   Matlab's randi function instead for versions 2008b and newer.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Concatenate to get the dimension vector
dims = [varargin{:}];
dims = num2cell(dims);

% Generate random integers
r = ceil((a-1) + rand(dims{:})*(b-a+1));