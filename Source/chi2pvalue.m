function p = chi2pvalue(chi2, dof)
%CHI2PVALUE Chi-square p-value
%   p = chi2pvalue(chi2, dof)
%
%   Given a chi-square distribution with dof degrees of freedom, this
%   function computes the upper-tail p-value of chi2.
%
%   Mathematically: p = 1 - gammainc(chi2/2, dof/2)
%
%   Note: gammainc(x,a) in matlab is the regularized gamma function, which
%   appears in math textbooks as P(s,x) with s being a (the exponent) and x
%   being the upper integration limit.

% (c) 2010 David R Hagen & Bruce Tidor
% This work is released under the MIT license.

p = 1 - gammainc(chi2/2, dof/2);