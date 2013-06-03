function [p, chi2] = mvnpvalue(x, mu, V)
%MVNPVALUE Multivariant normal p-value
%   p = mvnpvalue(x)
%   p = mvnpvalue(x, mu)
%   p = mvnpvalue(x, mu, V)
%   
%	Given the multivariant normal distribution described by mu and V,
%   this function computes the probability that vector x would end up so
%   many standard deviations away in the distribution. If mu is empty or
%   not provided, it is assumed to be a zero vector length of x. If V is
%   empty or missing, it is assumed to be the identity maxtrix of
%   appropriate size.
%
%   Mathematically: chi2 = (x-mu)' * V^-1 * (x-mu)
%                   p    = chi2pvalue(chi2, length(x))
%
%   [p, chi2] = pvalue(...) also returns the corresponding chi-square value

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

% Clean up inputs
if nargin < 3
    V = [];
    if nargin < 2
        mu = [];
    end
end

k = length(x);

if isempty(mu)
    mu = zeros(k,1);
end
if isempty(V)
    V = eye(k);
end

x = vec(x);
mu = vec(mu);

% Compute chi-square value
deltax = x - mu;
chi2 = deltax.' * (V \ deltax);

% Compute p value of chi2
p = chi2pvalue(chi2, k);