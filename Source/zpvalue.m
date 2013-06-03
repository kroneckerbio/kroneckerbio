function p = zpvalue(x, mu, sigma)
%ZPVALUE Upper-tail univariate normal p-value
%
%   p = zpvalue(x)
%   p = zpvalue(x, mu)
%   p = zpvalue(x, mu, sigma)
%
%   Given the univariant normal distribution described by mu and sigma,
%   this function computes the upper-tail p-value of x.
%
%   Mathematically: Z = (x-mu)/sigma
%                   p = 1 - 1/2*(1+erf(Z/sqrt(2)))

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
if nargin < 3
    sigma = [];
    if nargin < 2
        mu = [];
    end
end

n = size(x);

if isempty(mu)
    mu = zeros(n);
end
if isempty(sigma)
    sigma = ones(n);
end

% Convert to standard normal statistic
Z = (x - mu) ./ sigma;

% Compute the p-value
p = 1 - 0.5 * (1 + erf(Z/sqrt(2)));